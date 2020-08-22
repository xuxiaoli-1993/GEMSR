module gems_bound
   use gems_disc
   use mpi
contains

   !********************************************************************
   !  set boundary condition
   !********************************************************************
   subroutine set_boundary(faces, nodes, cells, interf)
      implicit none
      type(face), pointer :: faces(:), cf
      type(itf) :: interf
      type(node), pointer :: nodes(:)
      type(cell), pointer :: cells(:), cc
      !  type(cell),pointer,save::cg(:)
      type(vector), pointer :: qvr, qvl
      real(rfp) :: fac, vecn(ndim), unl, dl2, rho0, c0, c2, tr, unr, omega(3)
      integer :: i, j, ist, n
      !
      real(rfp), allocatable :: rb(:)
      integer :: id, ib, nd, ierr, st(mpi_status_size)
      !
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      !
      ib = interf % nitf + 1
      do i = ib, size(faces)
         cf => faces(i)
         if (cf % itype == 0) exit ! no boundary faces, skip
         !
         !     allocate(cf%right_cell)
         cc => cf % right_cell
         allocate(cc % weight(size(cf % f2n)), cc % gradient(neq, s_nrec), stat = ist)
         if (s_idt > 0) allocate(cc % qvn(s_idt), stat = ist)
         if (neqm > 0) allocate(cc % jv(3))
         if (ist /= 0) print *, 'error in allocate cell_right', i, ist
         fac = dot_product(cf % left_cell % centp - cf % centp, cf % vecn)
         !
         cc % centp = cf % left_cell % centp - two * fac * cf % vecn
         cc % vol = cf % left_cell % vol
         cc % gradient = zero
         cc % rp = zero
         cc % dqv % v = zero
         cc % itype = cf % left_cell % itype
         !     
      end do
      !
      call setup_special_faces(faces, interf) !pressure equlibrium condition
      !
   end subroutine set_boundary

   subroutine setup_special_faces(faces, interf)
      implicit none
      type(face), pointer :: faces(:)
      type(itf) :: interf
      integer :: it, i
      !
      do i = 1, size(bc)
         it = mod(bc(i) % itype, 1000)
         select case(bc(i) % igrp)
         case(outlet_type)
            if (it == 2) &
               call collect_face(bc(i), i, faces, interf % nitf)
         case(geom_topology)
         end select
      end do
   end subroutine setup_special_faces

   subroutine collect_face(mybc, ib, faces, nb)
      type(bc_type) :: mybc
      real(rfp), pointer :: r(:)
      type(face), pointer :: faces(:), cf
      integer :: nf, i, j, ib, nb
      if (mybc % n == 0) return ! no swirll outlet in this partition
      nf = mybc % n
      ispecialface = 1 !nf
      !
      allocate(mybc % facelist(nf), r(nf))
      !
      nf = 0
      do i = nb + 1, size(faces)
         cf => faces(i)
         if (cf % itype == 0) exit ! no boundary faces, skip
         if (cf % itype /= ib) cycle
         nf = nf + 1
         mybc % facelist(nf) = i
         if (ndim == 2) then
            r(nf) = cf % centp(s_iaxis)
         else
            r(nf) = sqrt(sum(cf % centp**2) - cf % centp(s_iaxis)**2)
            r(nf) = r(nf) + atan2(cf % centp(ndim), cf % centp(ndim - 1)) * 1.e-8_rfp ! theta
         end if
      end do
      !
      call sort(r, mybc % facelist)

      !   if(ndim == 3)  then   ! spefical treatment for 3d case
      !    iswirl =     bc(ib)%vars(2)
      !    ib = 1
      !     cf => faces(facelist(1))
      !     r(1) = sqrt(sum(cf%centp*2) - cf%centp(s_iaxis)**2)
      !     do i = 2 , nf 
      !     cf => faces(facelist(i))
      !     r(2) = sqrt(sum(cf%centp*2) - cf%centp(s_iaxis)**2)
      !     if(abs(r(2)/r(1) - one) < 1.e-6_rfp) then
      !      ib = ib + 1
      !     else
      !      exit
      !     end if
      !    end do
      !    if(mod(nf,ib) /= 0) print *,'the mesh is not cylindrical distribution',nf,ib
      !    iswirl = ib
      !   end if 
      deallocate(r)
   end subroutine collect_face
   !

   !********************************************************************
   !  boundary condition
   !********************************************************************
   subroutine boundary_condition(faces, nodes, cells, interf)
      implicit none
      type(face), pointer :: faces(:), cf
      type(node), pointer :: nodes(:)
      type(cell), pointer :: cells(:), lc
      type(itf) :: interf
      type(bc_type), pointer :: mybc
      integer :: i, it, ib, im(1), imass = 0, dbdummy
      real(rfp) :: dd(ndim)

      ib = interf % nitf + 1
      call setup_initial_bc(bc)
      !
      do i = ib, size(faces)
         cf => faces(i)
         lc => cf%left_cell

         it = cf % itype
         if (it == 0) exit ! no boundary faces, exit
         mybc => bc(it)
         select case(mybc % igrp) ! boundary type
         case(inlet_type)
            call inlet_boundary(cf, mybc)
         case(outlet_type)
            call outlet_boundary(cf, mybc)
            !     if(mybc%vars(3) > zero) imass = 1   ! need to balance mass flow
         case(farfield_type)
            call farfield_boundary(cf, mybc)
         case(wall_type)
            call wall_boundary(cf, mybc, faces, interf)
         case(geom_topology)
            call geom_boundary(cf, mybc, cells)
         case(mhd_type)
            call mhd_boundary(cf, mybc, mybc % itype)
         case default
            print *, ' Donot know the boundary type', bc(it) % itype, it, bc(i) % igrp
         end select
         call radiation_boundary(cf, mybc)
         ! calculate the gradiant at ghost cell
         !    dd = cf%left_cell%centp-cf%right_cell%centp
         !    im = maxloc(abs(dd))
         !    cf%right_cell%gradient(:,im(1)) = (cf%left_cell%qv%v - cf%right_cell%qv%v) / dd(im(1))
         !    if(s_ivis == 2) cf%right_cell%gradient(ike:iom,:) = zero

      end do
      !
      call special_bc_conditions(faces)
      if (imass > 0 .and.mod(s_nstep, 100) == 0) call balance_mass(faces)

   end subroutine boundary_condition

   subroutine special_bc_conditions(faces)
      type(face), pointer :: faces(:), cf
      integer :: i, it
      if (ispecialface == 0) return
      do i = 1, size(bc)
         it = mod(bc(i) % itype, 1000)
         select case(bc(i) % igrp)
         case(outlet_type)
            if (it == 2) &
               call swirl_outlet1(faces, bc(i) % facelist, int(bc(i) % vars(2)))
         case(geom_topology)
            if (it == 4) bc(i) % n = 0 ! wall temp distribution
         end select
      end do
   end subroutine special_bc_conditions

   subroutine setup_initial_bc(bc)
      type(bc_type) :: bc(:)
      integer :: i, it
      !  bc%n = 0
      do i = 1, size(bc)
         it = mod(bc(i) % itype, 1000)
         select case(bc(i) % igrp)
         case(inlet_type)
            if (it == 6) bc(i) % n = 0 ! inlet fan
            if (it == 40) bc(i) % n = 0 ! inlet fan
         case(wall_type)
            if (it == 20) bc(i) % n = 0 ! wall temp distribution
         case(outlet_type)
            if (it == 5) bc(i) % n = 0 ! wall temp distribution
         end select
      end do
   end subroutine setup_initial_bc

   subroutine balance_mass(faces)
      implicit none
      type(face), pointer :: faces(:), cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qv
      real(rfp) :: mass, mass0
      integer :: i, j, il = 4, ih = 6
      do i = 1, size(bc)
         if (bc(i) % igrp /= outlet_type) cycle
         if (bc(i) % n == 0) cycle
         mass = zero
         mass0 = bc(i) % vars(3)
         do j = 1, size(faces)
            cf => faces(j)
            if (cf % itype == 0) exit
            if (cf % itype /= i) cycle
            qv = (cf % left_cell % qv + cf % right_cell % qv) * 0.5_rfp
            mass = mass + rhof(qv) * sum(qv % v(imb:ime) * cf % vecn) * cf % area
         end do
         ! convergen
         if (abs(mass - mass0)/mass0 < 1.e-4_rfp) exit

         if (mass > mass0) then
            bc(i) % vars(ih) = mass
            bc(i) % vars(ih + 1) = bc(i) % vars(1) ! back pressure
         else
            bc(i) % vars(il) = mass
            bc(i) % vars(il + 1) = bc(i) % vars(1) ! back pressure
         end if
         ! initialize
         if (bc(i) % vars(il) <= -19621011._rfp) then
            bc(i) % vars(il) = mass
            bc(i) % vars(il + 1) = bc(i) % vars(1)
         end if
         if (bc(i) % vars(ih) <= -19621011._rfp) then
            bc(i) % vars(ih) = mass
            bc(i) % vars(ih + 1) = bc(i) % vars(1)
         end if

         if (bc(i) % vars(il) > mass0) then
            bc(i) % vars(1) = two * bc(i) % vars(il + 1) - g_pbase
         else if (bc(i) % vars(ih) < mass0) then
            bc(i) % vars(1) = half * (bc(i) % vars(ih + 1) - g_pbase)
         else
            bc(i) % vars(1) = half * (bc(i) % vars(ih + 1) - bc(i) % vars(il + 1))
         end if

      end do

   end subroutine balance_mass
   !
   subroutine inlet_boundary(cf, bc)
      implicit none
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      real(rfp) :: h0, s, vl2, h, p, t, dh, ds, dp, dt, rho, hp, ht, u
      real(rfp) :: vl(ndim), vr(ndim), omega(3)
      real(rfp) :: rhol, rhor, pl, pr, c0, rho0, pa, rhoa
      integer :: istop, it
      !
      ! variables for computation
      !
      cl => cf % left_cell
      cr => cf % right_cell
      qvl => cl % qv
      qvr => cr % qv
      ! EB boundary condition
      if (neqm > 0) call mhd_boundary(cf, bc, bc % itype / 1000)
      ! combination boundary
      it = mod(bc % itype, 1000) ! > 1000 means EB boundary
      !
      vl = qvl % v(imb:ime)
      vr = qvr % v(imb:ime)
      select case(it)
      case(0, 1, 8)
         vl2 = sum(vl**2)
         if (vl2 > sound_speed(qvl)) then
            qvr % v(:neqf) = bc % vars(:neqf)
            return
         end if
         u = sqrt(vl2)
         h0 = bc % vars(neq + 1)
         s = bc % vars(neq + 2)
         if (it == 8) then
            omega = zero
            omega(s_iaxis) = s_omega
            omega = acrossb(omega, avto3(cf % centp))
            h0 = h0 + half * sum(omega**2)
         end if
         !    call hs_condition(h0,s,vl,vr,qvr)
         call qvfromhs(h0 - half * vl2, s, qvr)
         !
         if (it == 0) then
            qvr % v(imb:ime) = u * bc % vars(neq4:neq3d)
         else if (it == 1) then
            qvr % v(imb:ime) = -u * cf % vecn
         else if (it == 8) then
            u = sqrt(abs(u**2 - sum(omega**2)))
            qvr % v(imb:ime) = u * bc % vars(neq4:neq3d) - omega(:ndim)
         end if

         if (neqf > ien) qvr % v(ien + 1:neqf) = bc % vars(ien + 1:neqf)
         !
      case(2)
         qvr % v(:neqf) = bc % vars(:neqf)
         !dumping the boundary conditions
         !    qvr%v(:neqf) = half * ( qvr%v(:neqf) + bc%vars(:neqf))
         qvr % v(ico) = qvl % v(ico)
      case(10)
         qvr % v(:neqf) = bc % vars(:neqf)
      case(20)
         qvr % v = bc % vars(:neq)
         qvr % v(imb:ime) = -bc % vars(neq + 3) * cf % vecn
         qvr % v(ico) = qvl % v(ico)
      case(21) ! 1/7 law full developped pipe or duct flow
         qvr % v(:neqf) = bc % vars(:neqf)
         qvr % v(ico) = qvl % v(ico)
         p = bc % vars(neq + 1) ! thickness
         t = bc % vars(neq + 2) ! center location
         u = bc % vars(neq + 3) ! max velocity
         h = bc % vars(neq + 4 + ndim) ! max velocity
         qvr % v(imb) = u * min(one, abs((p - abs(cf % centp(2) - t))/h))**(one/7.0_rfp)
      case(22)
         qvr % v(:neqf) = bc % vars(:neqf)
         !dumping the boundary conditions
         !    qvr%v(:neqf) = half * ( qvr%v(:neqf) + bc%vars(:neqf))
         qvr % v(ien) = qvl % v(ien)
         qvr % v(ico) = qvl % v(ico)
      case(3) ! mass flow
         qvr % v(ico) = qvl % v(ico)
         !    qvr%v(imb:ime) = qvr%v(imb:ime) / rhof(qvr)
         qvr % v(ien:) = bc % vars(ien:neq)
         qvr % v(imb:ime) = -bc % vars(neq + 3) / rhof(qvr) * cf % vecn
      case(30) ! pulse
         qvr % v(ico) = qvl % v(ico)
         dt = int(s_elapsed_time / bc % vars(neq + 2)) * bc % vars(neq + 2)
         dt = s_elapsed_time - dt
         if (dt <= bc % vars(neq + 1)) then
            qvr % v(imb:ime) = -bc % vars(neq + 3) / rhof(qvr) * cf % vecn
            qvr % v(ien:) = bc % vars(ien:neq)
         else
            qvr % v(imb:ime) = -qvl % v(imb:ime)
            qvr % v(ien) = bc % vars(ien)
            if (ien < neq) qvr % v(ien + 1:) = qvl % v(ien + 1:)
         end if
      case(5)
         ! user define function
         call inlet_user_define(cf, bc)
      case(6)
         call inlet_fan(cf, bc)
      case(7)
         qvr % v(:neqf) = bc % vars(:neqf)
         qvr % v(ico) = qvl % v(ico)
         omega = zero
         omega(s_iaxis) = s_omega
         omega = acrossb(omega, avto3(cf % centp))
         qvr % v(imb:ime) = qvr % v(imb:ime) - omega(:ndim)
      case(40)
         call inlet_profile(cf, bc)
         return
      end select

      if (s_ivis == 2.and.s_visb == 1) then ! k-omega inlet condition
         u = max(sqrt(sum(qvr % v(imb:ime)**2)), 1.e-10_rfp)
         qvr % v(iom) = t_omegainf * u / b_lenref
         qvr % v(ike) = t_keinf * mumix(qvr) * qvr % v(iom) / rhof(qvr)
      end if

      if (s_iaux > 0) then ! swirl
         qvr % v(iaub) = bc % vars(iaub) * cf % centp(s_iaxis)
      end if
      !    if(nrad > 0) qvr%v(irad) = qvl%v(irad)

   end subroutine inlet_boundary

   subroutine hs_condition(h0, s0, vl, vr, qvr)
      type(vector) :: qvr
      integer :: istop
      real(rfp) :: h0, s0
      real(rfp) :: u, h, dh, s, ds, dp, dt, vl(:), vr(:)
      !    dt = (h0 - half * sum(vl**2))/htf(qvr)
      !    dp = exp((htf(qvr) * log(dt) - s0)/(g_runi / g_mwi(1)))
      !    qvr%v(ien) = dt
      !    qvr%v(ico) = dp - g_pbase
      !    return

      istop = -2
      h = one - h0f(qvr) / h0
      if (abs(h) <= 1.e-12_rfp) then ! care about round-off error
         dh = zero
      else
         dh = h0 * h
      end if
      dh = dh - half * sum((vl + vr) * (vl - vr))
      if (abs(dh) < 1.e-10_rfp) istop = istop + 1
      s = one - entropyf(qvr) / s0
      if (abs(s) <= 1.e-12_rfp) then
         ds = zero
         istop = istop + 1
      else
         ds = s0 * s
      end if
      if (istop == 0) return
      !
      dp = rhof(qvr) * (dh - ds * qvr % v(ien))
      dt = (dh - hpf(qvr) * dp) / htf(qvr)
      !
      qvr % v(ico) = qvr % v(ico) + dp
      qvr % v(ien) = qvr % v(ien) + dt
      !
   end subroutine hs_condition
   !
   subroutine inlet_profile(cf, bc)
      implicit none
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr
      integer :: i
      !
      cr => cf % right_cell
      qvr => cr % qv
      i = bc % n
      qvr % v(:neqf) = bc % vars(i + ico:i + neqf)
      bc % n = bc % n + neq
      !
   end subroutine inlet_profile

   subroutine inlet_fan(cf, bc)
      implicit none
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector) :: qvr, qvl
      type(vector) :: qv
      real(rfp) :: h0, s0, u, vl2, h, p, t, r0, alpha, beta, vd(3)
      integer :: i
      !
      cl => cf % left_cell
      cr => cf % right_cell
      qvl = cl % qv
      qvr = cr % qv
      i = bc % n
      qv % v(ico) = bc % vars(i + ico)
      qv % v(ien) = bc % vars(i + ien)
      qv % v(imb:ime) = zero
      vd(:ndim) = bc % vars(i + imb:i + ime)
      h0 = h0f(qv)
      s0 = entropyf(qv)
      u = sqrt(sum(qvl % v(imb:ime)**2))
      h = h0 - half * u * u
      call qvfromhs(h, s0, qvr)
      qvr % v(imb:ime) = u * vd(:ndim)
      if (s_ivis == 2) then ! k-omega inlet condition
         qvr % v(iom) = t_omegainf * u / b_lenref
         qvr % v(ike) = t_keinf * mumix(qvr) * qvr % v(iom) / rhof(qvr)
      end if
      cr % qv = qvr
      bc % n = bc % n + neq
      !
   end subroutine inlet_fan

   subroutine inlet_user_define(cf, bc)
      implicit none
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      real(rfp) :: h0, s, vl2, h, p, t, dt = 1.e-5_rfp
      real(rfp) :: r, theta0, U0, Uc, R0, T0, sos, x0
      real(rfp) :: V0, D0, K, rho, nu, mu, V, y0, eta, Vr(2)
      !
      cl => cf % left_cell
      cr => cf % right_cell
      qvl => cl % qv
      qvr => cr % qv
      !
      qvr % v(:neqf) = bc % vars(:neqf)
      qvr % v(ico) = qvl % v(ico)
      sos = sound_speed(qvl)

      ! nozzle-like inlet
      Uc = 500.
      T0 = 3500.
      R0 = 75e-6
      r = sqrt(cf%centp(1)**2 + cf%centp(3)**2)
      U0 = Uc

      if(r <= R0) then
         qvr%v(imb) = zero
         qvr%v(imb+2) = zero

         qvr%v(imb+1) = U0  
         qvr%v(ien) = T0  

         if(U0**2 > sos) then
            qvr%v(ico) = bc%vars(1)
         end if 
      end if
   end subroutine inlet_user_define

   subroutine swirl_outlet(faces, facelist, ni)
      implicit none
      type(face), pointer :: faces(:)
      type(vector), pointer :: qvr
      integer :: facelist(:)
      real(rfp) :: r1, r2, vt1, vt2, pb, v2
      integer :: i, ib, nf, j, ni
      !
      nf = size(facelist)
      do j = 1, ni
         ib = facelist(j)
         qvr => faces(ib) % right_cell % qv
         r1 = radius(faces(ib) % centp)
         vt1 = rhof(qvr) * vtheta(qvr, faces(ib) % centp, r1)**2 / r1
         pb = bc(faces(ib) % itype) % vars(ico)
         qvr % v(ico) = pb
         do i = ni + 1, nf, ni
            ib = facelist(i + j - 1)
            qvr => faces(ib) % right_cell % qv
            r2 = radius(faces(ib) % centp)
            vt2 = rhof(qvr) * vtheta(qvr, faces(ib) % centp, r2)**2 / r2
            pb = pb + half * (vt1 + vt2) * (r2 - r1)
            !  v2 = sum(qvr%v(imb:ime)**2)
            !  if(v2 < sound_speed(qvr)) qvr%v(ico) = pb
            v2 = sum(qvr % v(imb:ime) * faces(ib) % vecn)
            !  if(v2 < zero) qvr%v(imb:ime) = - qvr%v(imb:ime)
            if (v2 * v2 < sound_speed(qvr)) qvr % v(ico) = pb
            !  qvr%v(ico) = pb
            r1 = r2; vt1 = vt2
         end do
      end do
   end subroutine swirl_outlet

   subroutine swirl_outlet1(faces, facelist, ni)
      ! try average value
      implicit none
      type(face), pointer :: faces(:)
      type(vector), pointer :: qvr
      integer :: facelist(:), ni
      real(rfp) :: r1, r2, vt1, vt2, pb, v2
      integer :: i, ib, nf, j, k
      !
      nf = size(facelist)
      k = 0
      do i = 1, nf, ni
         vt2 = zero
         do j = 1, ni
            k = j + i - 1
            ib = facelist(k)
            qvr => faces(ib) % right_cell % qv
            r2 = radius(faces(ib) % centp)
            vt2 = vt2 + rhof(qvr) * vtheta(qvr, faces(ib) % centp, r2)**2 / r2
         end do
         !
         vt2 = vt2 / real(ni, rfp)
         if (i == 1) then
            pb = bc(faces(ib) % itype) % vars(ico)
         else
            pb = pb + half * (vt1 + vt2) * (r2 - r1)
         end if
         !
         do j = 1, ni
            k = j + i - 1
            ib = facelist(k)
            qvr => faces(ib) % right_cell % qv
            qvr % v(ico) = pb
         end do
         vt1 = vt2; r1 = r2
      end do
   end subroutine swirl_outlet1

   function radius(p)result(r)
      real(rfp) :: p(:), r
      if (ndim == 2) then
         r = p(s_iaxis)
      else
         r = sqrt(sum(p**2) - p(s_iaxis)**2)
      end if
   end function radius

   function vtheta(qv, p, r)result(v)
      type(vector) :: qv
      real(rfp) :: v, p(:), pc(3), r
      if (ndim == 2) then
         v = qv % v(iaub)
      else
         pc(:ndim) = p / r
         v = -qv % v(imb + 1) * pc(3) + qv % v(ime) * pc(2) ! u sin(q)-v cos(q)  
         !    v = v + s_omega * r
      end if
   end function vtheta

   subroutine outlet_boundary(cf, bc)
      implicit none
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      real(rfp) :: v2, r2, r1, r
      real(rfp) :: rhol, rhor, pl, pr, c0, rho0
      integer :: it
      !
      ! variables for computation
      !
      cl => cf % left_cell
      cr => cf % right_cell
      qvl => cl % qv
      qvr => cr % qv
      !
      if (neqf > 0) then
         it = mod(bc % itype, 1000)
         select case(it)
         case(0) ! back pressure
            qvr % v(:neqf) = qvl % v(:neqf)
            !v2 = sum(qvl % v(imb:ime)**2)
            !   if(v2 > sound_speed(qvl).and.qvl%v(ico) > bc%vars(1)) return  ! supersonic
            !if (v2 <= sound_speed(qvl)) qvr % v(ico) = bc % vars(1) ! subsonic flow
            qvr%v(ico) = bc%vars(1)  ! subsonic flow
         case(1)
            !   full developed flow
            qvr % v(:neqf) = qvl % v(:neqf)
         case(2) ! for farfield
            qvr % v(:neqf) = qvl % v(:neqf)
            qvr%v(ico) = bc%vars(1)
         case(3) ! makeup turbo-outlet
            qvr % v(:neqf) = qvl % v(:neqf)
            r = sum(cf % centp**2) - cf % centp(s_iaxis)**2
            r1 = bc % vars(4)**2; r2 = bc % vars(5)**2
            !    qvr%v(ico) = bc%vars(1) + half * rhof(qvl) * ( r2 - r0**2) * s_omega**2
            r = 3._rfp / 4._rfp * (r - r1) / (r2 - r1)
            qvr % v(ico) = bc % vars(1) *(half + r) + g_pbase * (r - half)
         case(5) ! given back pressure distribution
            qvr % v(:neqf) = qvl % v(:neqf)
            bc % n = bc % n + 1
            qvr % v(ico) = bc % vars(bc % n)
         end select
      end if
      !
      !    if(nrad > 0) qvr%v(irad) = qvl%v(irad)
      ! EB boundary condition
      it = bc % itype / 1000
      if (neqm > 0) call mhd_boundary(cf, bc, it)
   end subroutine outlet_boundary

   subroutine farfield_boundary(cf, bc)
      implicit none
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      real(rfp) :: un
      !
      ! variables for computation
      !
      cl => cf % left_cell
      qvl => cl % qv
      !
      un = sum(cf % vecn * qvl % v(imb:ime)) ! normal velocity
      if (un < zero) then ! inlet
         call inlet_boundary(cf, bc)
      else
         call outlet_boundary(cf, bc)
      end if
   end subroutine farfield_boundary

   subroutine wall_boundary(cf, bc, faces, interf)
      implicit none
      type(face), pointer :: faces(:)
      type(itf) :: interf
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      type(vector) :: qvf ! the values on th e face
      real(rfp) :: dd(ndim), dn, tw, alpha, omega(3)
      integer :: i, im(1), it, gravity_direction

      cl => cf % left_cell
      cr => cf % right_cell
      qvl => cl % qv
      qvr => cr % qv
      qvr = qvl
      qvf = qvl
      i = s_iaxis
      it = mod(bc % itype, 1000)
      if(ndim==1) then
         gravity_direction=1
      else
         gravity_direction=2
      end if

      select case(it)
      case(0)
         call symmetry(cf, bc)
      case(10)
         call symmetry(cf, bc)
         qvr % v(ien) = two * bc % vars(ndim + 1) - qvl % v(ien)
         !    qvr%v(ien) = bc%vars(ndim+1)
      case(1)
         if (nmeq > 0) qvr % v(imb:ime) = two * bc % vars(:ndim) - qvl % v(imb:ime)
         if (s_ivis == 2) call no_wall_function(cf, bc % vars(ndim + 2))
         qvr % v(ien) = qvl % v(ien) + &
            sum((cr % centp - cl % centp) * cf % vecn) * bc % vars(ndim + 3) / lamda(qvl, cl)
         if (s_iaux == 1) then ! swirl
            qvr % v(iaub) = two * bc % vars(ndim + 4) - qvl % v(iaub)
            !   dp/dr=vq^2/r
            qvr % v(ico) = qvl % v(ico) + rhof(qvl) * (cr % centp(i) - cl % centp(i)) * &
               bc % vars(ndim + 4)**2 /cf % centp(i)
         end if
         !    qvr%v(ico) = qvr%v(ico) + viscous_correction(cf,nodes)
         !    qvr%v(ico) = qvl%v(ico) + sum(cl%gradient(1,:) * (cr%centp - cl%centp))
         !     qvr%v(ico) = viscous_correction(cf,nodes) * two - qvl%v(ico)
         !    qvf%v(imb:ime) = zero
         !    dd = cf%centp-cr%centp
         !    im = maxloc(abs(dd))
         !    cr%gradient(:,im(1)) = (qvf%v - qvr%v) / dd(im(1))
      case(2)
         if (nmeq > 0) qvr % v(imb:ime) = two * bc % vars(:ndim) - qvl % v(imb:ime)
         if (s_ivis == 2) call no_wall_function(cf, bc % vars(ndim + 2))
         tw = bc % vars(ndim + 1)
         !    if(tw < qvl%v(ien)) then
         !    qvr%v(ien) = max(two * tw - qvl%v(ien),0.6_rfp * qvl%v(ien))
         !    else
         !    qvr%v(ien) = min(two * tw - qvl%v(ien),1.5_rfp * qvl%v(ien))
         !    end if
         qvr % v(ien) = tw * tw / qvl % v(ien)
         if (s_iaux == 1) then ! swirl
            qvr % v(iaub) = two * bc % vars(ndim + 4) - qvl % v(iaub)
            !   dp/dr=vq^2/r
            qvr % v(ico) = qvl % v(ico) + rhof(qvl) * (cr % centp(i) - cl % centp(i)) * &
               bc % vars(ndim + 4)**2 /cf % centp(i)
         end if
      case(20)
         if (nmeq > 0) qvr % v(imb:ime) = -qvl % v(imb:ime)
         if (s_ivis == 2) call no_wall_function(cf, bc % vars(ndim + 2))
         tw = bc % vars(bc % n + ien)
         qvr % v(ien) = tw * tw / qvl % v(ien)
         if (s_iaux == 1) then ! swirl
            qvr % v(iaub) = two * bc % vars(ndim + 4) - qvl % v(iaub)
            !   dp/dr=vq^2/r
            qvr % v(ico) = qvl % v(ico) + rhof(qvl) * (cr % centp(i) - cl % centp(i)) * &
               bc % vars(ndim + 4)**2 /cf % centp(i)
         end if
         bc % n = bc % n + neq
         !
      case(3)
         !
         qvr % v(imb:ime) = two * bc % vars(:ndim) - qvl % v(imb:ime)
         if (s_ivis == 2) call turb_wall_function(cf)
         qvr % v(ien) = qvl % v(ien) + sum(cr % centp - cl % centp * cf % vecn) * bc % vars(ndim + 3) / lamda(qvl, cl)
      case(4)
         !
         qvr % v(imb:ime) = two * bc % vars(:ndim) - qvl % v(imb:ime)
         if (s_ivis == 2) call turb_wall_function(cf)
         qvr % v(ien) = two * bc % vars(ndim + 1) - qvl % v(ien)
         ! tw = bc % vars(ndim + 1)
         ! qvr % v(ien) = tw * tw / qvl % v(ien)
      case(5) ! rotation
         qvr % v(imb:ime) = bc % vars(3) * (cf % centp(1) - bc % vars(1))
         qvr % v(imb:ime) = two * qvr % v(imb:ime) - qvl % v(imb:ime)
      case(6) ! rotation
         qvr % v(imb) = -bc % vars(3) * (cf % centp(2) - bc % vars(2))
         qvr % v(imb + 1) = bc % vars(3) * (cf % centp(1) - bc % vars(1))
         qvr % v(imb:ime) = two * qvr % v(imb:ime) - qvl % v(imb:ime)
         qvr % v(ien) = two * bc % vars(4) - qvl % v(ien)
      case(7) ! moving body with rotation
         qvr % v(imb:ime) = two * mesh_vel(cf % centp) - qvl % v(imb:ime)
         ! Pressure BC in all rotating walls
         ! dp/dn = RhoOmega^2(yn_y+zn_z)
         dd = cr % centp - cl % centp
         dn = sqrt(sum(dd * 2))
         dd = dd * cf % centp !/ dn
         dd(s_iaxis) = zero
         qvr % v(ico) = qvl % v(ico) + rhof((qvl + qvr) * half) * s_omega**2 * sum(dd)
         !
         if (s_ivis == 2) call no_wall_function(cf, bc % vars(ndim + 2))
      end select
      !
      if (g_igrav > 0) then
         qvr % v(ico) = qvl % v(ico) + rhof(qvl) * sum((cr % centp(gravity_direction) - cl % centp(gravity_direction)) * &
            cf % vecn) * sum(g_gravity * cf % vecn)
      end if
      !
      ! EB boundary condition
      it = bc % itype / 1000
      if (neqm > 0) call mhd_boundary(cf, bc, it)
   end subroutine wall_boundary

   subroutine radiation_boundary(cf, bc)
      implicit none
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(bc_type), pointer :: bc
      type(vector), pointer :: qvr, qvl
      real(rfp) :: dn, tw, alpha
      integer :: it, i
      cr => cf % right_cell
      cl => cf % left_cell
      qvl => cl % qv
      qvr => cr % qv
      if (nrad == 0) return
      dn = sum(cf % vecn * (cr % centp - cl % centp))
      alpha = two * two / (three * absorption_coef(qvl) * dn)
      !    tw =  bc%vars(3)  !half * ( qvr%v(ien) + qvl%v(ien))
      tw = half * (qvr % v(ien) + qvl % v(ien))
      qvr % v(irad) = (8._rfp * stephan_boltzmann * tw**4 - &
         (one - alpha) * qvl % v(irad))/(one + alpha)
      select case(bc % igrp) ! boundary type
      case(inlet_type, outlet_type, farfield_type)
         qvr % v(irad) = qvl % v(irad)
         !    qvr%v(irad) = (8._rfp * stephan_boltzmann * tw**4 - &
         !                    (one - alpha) * qvl%v(irad))/(one + alpha)     
      case(wall_type)
         it = mod(bc % itype, 1000)
         select case(it)
         case(0, 7)
            qvr % v(irad) = qvl % v(irad)
            !      case(2)
            !    qvr%v(irad) = 4._rfp * stephan_boltzmann * tw**4 
         case default
            !    qvr%v(irad) = 8._rfp * stephan_boltzmann * tw**4 - qvl%v(irad)
            !    alpha =  two * two / (three * 100._rfp * dn)
            !    qvr%v(irad) = (8._rfp * stephan_boltzmann * tw**4 - &
            !                    (one - alpha) * qvl%v(irad))/(one + alpha)     

         end select
      end select
   end subroutine radiation_boundary

   function viscous_correction(cf, nodes) result(dpdx)
      type(node), pointer :: nodes(:)
      type(face), pointer :: cf, sf
      type(cell), pointer :: cc
      type(vector) :: qv, qvf, gqv(ndim)
      real(rfp) :: zmu, zk, vecn(ndim), vjacob(ndim, 2), dpdx, coef, zd(nspe), tauf, tauc
      integer :: i, k1, k2
      integer, pointer :: f2n(:)
      !
      dpdx = zero
      return
      cc => cf % left_cell
      zmu = mumix(cc % qv)
      zk = lamda(cc % qv, cc)
      qv % v = zero
      f2n => cf % f2n
      !   call face_gradient(cf%right_cell,cf%left_cell,nodes(f2n),gqv,cf%vecn)
      tauf = sum(qv % v(imb:ime) * cf % vecn)

      tauc = sum(qv % v(imb:ime) * cf % vecn)

      dpdx = two * (tauf - tauc)
   end function viscous_correction

   subroutine geom_boundary(cf, bc, cells)
      implicit none
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl, cells(:)
      type(vector), pointer :: qvr, qvl

      select case(mod(bc % itype, 1000))
      case(0)
         call symmetry(cf, bc)
      case(1)
         call axisymmetry(cf)
      case(2)
         call periodic(cf, cells)
      case(3)
         call centerline(cf)
      end select
   end subroutine geom_boundary

   subroutine symmetry(cf, bc)
      implicit none
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      real(rfp) :: un

      cl => cf % left_cell
      cr => cf % right_cell
      qvl => cl % qv
      qvr => cr % qv
      qvr = qvl
      un = sum(cf % vecn * qvl % v(imb:ime))
      qvr % v(imb:ime) = qvl % v(imb:ime) - two * un * cf % vecn(:)
      !    qvr%v(imb:ime) = qvl%v(imb:ime) -  un * cf%vecn(:)
      !
      ! EB boundary condition
      if (neqm > 0) call mhd_boundary(cf, bc, bc % itype / 1000)
   end subroutine symmetry

   subroutine centerline(cf)
      implicit none
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      real(rfp) :: un
      cl => cf % left_cell
      cr => cf % right_cell
      qvl => cl % qv
      qvr => cr % qv
      qvr = qvl
      un = sum(cf % vecn * qvl % v(imb:ime))
      qvr % v(imb:ime) = qvl % v(imb:ime) - two * un * cf % vecn(:)
      if (nrad > 0) qvr % v(irad) = qvl % v(irad)
      if (s_iaux > 0) qvr % v(iaub) = -qvl % v(iaub)
   end subroutine centerline

   subroutine axisymmetry(cf)
      implicit none
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      real(rfp) :: rr, rl, ur, ut, cos_sin(2)

      cl => cf % left_cell
      cr => cf % right_cell
      qvl => cl % qv
      qvr => cr % qv
      qvr = qvl
      !    cr%gradient = cl%gradient
      !
      rl = one / sqrt(cl % centp(2)**2 + cl % centp(3)**2)
      cos_sin = cl % centp(2:3) * rl
      ur = cos_sin(1) * qvl % v(imb + 1) + cos_sin(2) * qvl % v(ime)
      ut = -cos_sin(2) * qvl % v(imb + 1) + cos_sin(1) * qvl % v(ime)
      !
      rr = one / sqrt(cr % centp(2)**2 + cr % centp(3)**2)
      cos_sin = cr % centp(2:3) * rr
      qvr % v(imb + 1) = cos_sin(1) * ur - cos_sin(2) * ut
      qvr % v(ime) = cos_sin(2) * ur + cos_sin(1) * ut
      !
   end subroutine axisymmetry

   subroutine periodic(cf, cells)
      implicit none
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl, cells(:)
      return ! do nothing
   end subroutine periodic

   subroutine periodic_correction(cf, cells)
      implicit none
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl, cells(:)
      type(vector), pointer :: qvr
      real(rfp) :: angle, rr, u(2), p(2), cos_sin(2)
      integer :: i, j

      cr => cf % right_cell
      qvr => cr % qv
      if (s_iperiod > 3) return ! no need to modify the data
      !
      !    cr%gradient = cl%gradient
      !
      !    rl = one / sqrt(sum((cl%centp(:2)- bc(cf%itype)%vars(:2))**2))
      !    cos_sin = (cl%centp(1:2) - bc(cf%itype)%vars(:2)) * rl
      !    ur =  cos_sin(1)*qvl%v(imb) + cos_sin(2)*qvl%v(imb+1)
      !    ut = -cos_sin(2)*qvl%v(imb) + cos_sin(1)*qvl%v(imb+1)
      !
      j = 0
      do i = 1, ndim
         if (i == s_iperiod) cycle
         j = j + 1
         u(j) = qvr % v(imb + i - 1)
         p(j) = cr % centp(i) - s_pcenter(i)
      end do
      rr = one / sqrt(sum(p**2))
      cos_sin = p * rr
      ![ux]=[cos   -sin][ur]
      ![uy]=[sin    cos][ut]
      !
      !    angle = bc(cf%itype)%vars(4)
      !    cos_sin(1) = cos(angle)
      !    cos_sin(2) = sin(angle)
      p(1) = cos_sin(1) * u(1) - cos_sin(2) * u(2)
      p(2) = cos_sin(2) * u(1) + cos_sin(1) * u(2)
      j = 0
      do i = 1, ndim
         if (i == s_iperiod) cycle
         j = j + 1
         qvr % v(imb + i - 1) = p(j)
      end do
   end subroutine periodic_correction

   subroutine turb_wall_function(cf)
      type(face), pointer :: cf
      type(cell), pointer :: cl, cr
      type(vector) :: qv
      real(rfp), pointer :: vecn(:)
      real(rfp) :: vel(ndim), velw(ndim), yp, up, kp, ypn, upn, rho, zmu, dy, ut
      real(rfp) :: cons, karman, kinv, beta2, E, omega
      beta2 = sqrt(t_betas0) !**0.25_rfp
      karman = 0.4_rfp
      kinv = one / karman
      cons = 5.1_rfp
      E = 9.0_rfp
      vecn => cf % vecn
      cl => cf % left_cell
      qv = cl % qv
      !
      !   if(s_ialg /= 0 ) &
      !    qv = qv + qv_from_highorder(cl,cf)
      rho = rhof(qv)
      zmu = mumix(qv)
      kp = qv % v(ike)
      omega = qv % v(iom)
      vel = qv % v(imb:ime)
      !
      up = sqrt(dot_product(vel, vel) - dot_product(vel, vecn)**2)
      dy = dot_product(vecn, cf % centp - cl % centp)
      yp = rho * sqrt(beta2) * dy * sqrt(kp) / zmu
      !
      !
      ut = cf % right_cell % rp ! rp store ut
      if (ut <= zero) ut = sqrt(zmu * up / (rho * dy))
      call Reichardt_wall_law(up, upn, ut, dy, yp, rho, zmu, karman, kinv)

      !  if( yp > 3.0_rfp) then   !11.3
      !   ut = sqrt(sqrt(beta2) * karman * sqrt(kp) * up / log(yp*e))
      !   upn = 5.0 * ut
      velw = vel - dot_product(vel, vecn) * vecn ! decide tangential velocity
      velw = velw * upn / (sqrt(sum(velw * velw)) + mytiny) ! cal velocity at boundary
      qv % v(imb:ime) = two * velw - vel ! slip wall with velw
      !   qv%v(imb:ime) = velw   ! slip wall with velw
      !   cl%qv%v(ike) = ut * ut * beta2
      !   cl%qv%v(iom) = two * ut * kinv / (beta2 * dy)
      !  else
      !   ut = sqrt(zmu * up / (rho * dy))
      qv % v(ike) = zero
      !   qv%v(iom) = two * ut * kinv / (beta2 * dy)
      qv % v(iom) = 400.0_rfp * zmu / (rho * dy * dy) !- omega
      !
      cf % right_cell % rp = ut
      cf % right_cell % qv = qv
   end subroutine turb_wall_function

   subroutine no_wall_function(cf, ks)
      type(face), pointer :: cf
      type(vector), pointer :: qv
      real(rfp) :: vel(ndim), yp, up, rho, zmu, dy, ut, k, omega, ks, utaw, ksplus, sr
      qv => cf % left_cell % qv
      rho = rhof(qv)
      zmu = mumix(qv)
      vel = qv % v(imb:ime)
      k = qv % v(ike)
      omega = qv % v(iom)
      ut = sum(vel * vel) - sum(vel * cf % vecn)**2
      dy = dot_product(cf % vecn, cf % centp - cf % left_cell % centp)
      if (ks < mytiny .or. ut < mytiny) then
         ksplus = zero
      else
         utaw = zmu * sqrt(ut) / dy
         utaw = sqrt(utaw / rho)
         ks = max(dy, ks) ! correction for grid
         ksplus = ks * utaw * rho / zmu
      end if
      !  utaw = zmu * sqrt(ut) / dy
      !  utaw = sqrt(utaw / rho )
      !  qv => cf%left_cell%qv
      !  qv%v(ike) = utaw**2 / 0.3_rfp
      !  qv%v(iom) = utw / ( 0.3_rfp * 0.41_rfp * dy) 
      qv => cf % right_cell % qv
      qv % v(ike) = zero
      if (ksplus > mytiny) then
         if (ksplus < 25.0_rfp) then
            sr = (50.0_rfp / ksplus)**2
         else
            sr = 100.0 / ksplus
         end if
         !    qv%v(iom) = 2500.0_rfp * zmu / (rho * ks * ks)  !- omega
         qv % v(iom) = utaw**2 * rho * sr / zmu
      else
         qv % v(iom) = 400.0_rfp * zmu / (rho * dy * dy) !- omega
      end if
      !  try dk^2/dy^2=0
      !  qv%v(ike) = cf%left_cell%qv%v(ike) +  &
      !    sum(cf%left_cell%gradient(ike,:) * (cf%right_cell%centp - cf%left_cell%centp))
      !    if(s_nstep < 1000) then
      !       qv%v(iom) = min(10._rfp * omega, qv%v(iom))
      !       qv%v(ike) = k * 0.2_rfp
      !    end if
      !    qv%v(iom) = min(1.e6_rfp, qv%v(iom))
   end subroutine no_wall_function

   subroutine Reichardt_wall_law(u, un, ut, y, yp, rho, zmu, karman, kinv)
      real(rfp), intent(in) :: u, y, rho, zmu, karman, kinv
      real(rfp), intent(inout) :: ut, un, yp
      real(rfp) :: dypdut, yp11, eyp11, eyp33, f0, df, up
      integer :: i

      dypdut = rho * y / zmu
      do i = 1, 4
         yp = dypdut * ut !rho * ut * y / zmu
         yp11 = yp / 11.0_rfp
         eyp11 = exp(-yp11)
         eyp33 = exp(-0.33_rfp * yp)
         up = kinv * log(one + karman * yp)+ &
            7.8_rfp * (one - eyp11 - yp11 * eyp33)
         f0 = ut * up - u
         df = up + ut * dypdut * (one/(one + karman * yp) + &
            (7.8_rfp/11._rfp)*(eyp11 - (one - 0.33_rfp * yp) * eyp33))
         ut = ut - f0 / (df + mytiny)
      end do
      !
      yp = dypdut * ut
      yp11 = yp / 11.0_rfp
      eyp11 = exp(-yp11)
      eyp33 = exp(-0.33_rfp * yp)
      df = ut * ut * rho / zmu * (one/(one + karman * yp)+ &
         (7.8_rfp/11._rfp)*(eyp11 - (one - 0.33_rfp * yp) * eyp33))
      un = max(u - y * df, zero)
   end subroutine Reichardt_wall_law

   subroutine mhd_boundary(cf, bc, it)
      implicit none
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      real(rfp) :: r, r1, r2, ci, p(ndim), btheta, current, current0
      real(rfp) :: bj0, by0, bj1, by1, omega
      real(rfp) :: nv(3), mv(3), lv(3)
      integer :: it
      cl => cf % left_cell
      cr => cf % right_cell
      qvl => cl % qv
      qvr => cr % qv
      qvr % v(ibb:iee) = zero ! initial
      ci = g_mu0 * bc % vars(1) * half / pi
      p = cf % centp
      r = sqrt(sum(p(:2) * p(:2)))
      if (ndim == 2.and.s_iaxis > 0) r = cf % centp(s_iaxis)
      select case(it)
      case(0)
         r1 = bc % vars(2); r2 = bc % vars(3)
         current = bc % vars(1)
         btheta = bthetaf(r1, r2, current, r)
         qvr % v(ibb) = -p(2) * btheta / r
         qvr % v(ibb + 1) = p(1) * btheta / r
         qvr % v(iee) = e_resistivity(qvl, cl) * current / pi / (r2 * r2 - r1 * r1)
      case(1)
         btheta = ci / r
         qvr % v(ibb) = -p(2) * btheta / r
         qvr % v(ibb + 1) = p(1) * btheta / r
      case(2) ! free boundary condition
         qvr % v(ieb:iee) = qvl % v(ieb:iee)
      case(20) ! free boundary with zero e-field
         qvr % v(ibb:ibe) = qvl % v(ibb:ibe)
         qvr % v(ieb:iee) = zero
      case(3)
         r1 = bc % vars(2); r2 = bc % vars(3)
         current = current_vs_time(bc % vars(1))
         current0 = current_vs_time(bc % vars(ifb))
         if (r < r1) then
            ci = current0
         else if (r < r2) then
            ci = current * (r * r - r1 * r1)/(r2 * r2 - r1 * r1) + current0
         else
            ci = current + current0
         end if
         qvr % v(ibe) = g_mu0 * ci / (two * pi * (r + mytiny))
         qvr % v(ieb:iee) = qvl % v(ieb:iee)
         !     qvr%v(ieb:iee) = zero
         !     qvr%v(ieb) = bc%vars(1)/(pi *(r2*r2-r1*r1)) * e_resistivity(cl%qv,cl)
      case(4) ! r = 0  centrline
         qvr % v(ibb:ibe) = zero
         qvr % v(ieb:iee) = qvl % v(ieb:iee)
         !     qvr%v(ieb) = zero
         qvr % v(ieb + 1:iee) = zero
      case(5)
         qvr % v(ieb:iee) = bc % vars(ieb:iee)
         qvr % v(ibb:ibe) = qvl % v(ibb:ibe)
      case(6) ! dielectric boundaru condition
         qvr % v(ieb) = qvl % v(ieb)
         r1 = sqrt(sum((cl % centp - cr % centp)**2)) * half
         r2 = cf % centp(2) / cr % centp(2)
         qvr % v(ibe) = qvl % v(ibe) + g_mu0 / e_resistivity(cl % qv, cl) * qvr % v(ieb) * r1
         qvr % v(ibe) = r2 * qvr % v(ibe)
      case(7)
         r1 = bc % vars(2); r2 = bc % vars(3)
         r = cf % centp(1)
         current = bc % vars(1)
         qvr % v(ibe) = zero
         if (r >= r1.and.r <= r2) qvr % v(ibe) = g_mu0 * current
         qvr % v(ieb:iee) = qvl % v(ieb:iee)
      case(8)
         r1 = bc % vars(2); r2 = bc % vars(3)
         ci = bc % vars(1)
         qvr % v(ibe) = g_mu0 * ci / (two * pi * (r + mytiny))
         qvr % v(ieb:iee) = qvl % v(ieb:iee)
      case(9)
         qvr % v = zero
         qvr % v(ibe) = bc % vars(1)
         qvr % v(ieb:iee) = qvl % v(ieb:iee)
      case(11)
         r1 = bc % vars(ibb + 1); r2 = bc % vars(ibe)
         current = bc % vars(ibb) / (pi * (r2 **2 - r1 **2))
         if (r < r1 .or. r > r2) current = zero
         qvr % v(ibb:ibe) = qvl % v(ibb:ibe)
         qvr % v(ieb) = bc % vars(ibb) !current * e_resistivity(cl%qv,cl)
      case(12)
         r1 = bc % vars(ndim + 6); r2 = bc % vars(ndim + 7)
         current = bc % vars(ndim + 5)
         !     current = current_vs_time(bc%vars(ndim+5))
         qvr % v(ibe) = bthetaf(r1, r2, current, r)
         qvr % v(ieb:iee) = qvl % v(ieb:iee)
      case(13) ! 2d linear distribution
         r1 = bc % vars(2)
         r2 = bc % vars(3)
         current = bc % vars(1)
         current0 = bc % vars(ifb)
         r1 = (p(1) - r1) / (r2 - r1)
         r1 = min(max(zero, r1), one)
         current = r1 * current + current0
         qvr % v(ibe) = current * g_mu0 / (two * pi * r)
         qvr % v(ieb:iee) = qvl % v(ieb:iee)
      case(23) ! 3d linear distribution
         r = sqrt(sum(p(2:3)**2))
         r1 = bc % vars(2)
         r2 = bc % vars(3)
         current = bc % vars(1)
         current0 = bc % vars(ifb)
         r1 = (p(1) - r1) / (r2 - r1)
         r1 = min(max(zero, r1), one)
         current = r1 * current + current0
         omega = current * g_mu0 / (two * pi * r)
         qvr % v(ibb + 1) = -omega * p(3) / r
         qvr % v(ibe) = omega * p(2) / r
         qvr % v(ieb:iee) = qvl % v(ieb:iee)
      case(14)
         r1 = bc % vars(ndim + 6)
         r2 = bc % vars(ndim + 7)
         current = bc % vars(ndim + 5)
         if (p(1) < r1.or.p(1) > r2) current = zero
         r1 = (p(1) - r1) / (r2 - r1)
         r1 = min(max(zero, r1), one)
         qvr % v(ibe) = r1 * bc % vars(ndim + 5) * g_mu0
         qvr % v(ieb:iee) = qvl % v(ieb:iee)
         !     qvr%v(ieb+1) = -current * e_resistivity(cl%qv,cl)
      case(15) ! Mike's suggested condition
         !
         qvr % v(ibb:iee) = qvl % v(ibb:iee)
         qvr % v(ieb + 1) = zero
         r1 = sum(qvl % v(imb:ime) * cf % vecn)
         r2 = sum((cr % centp - cl % centp) * cf % vecn)
         qvr % v(ibe) = qvl % v(ibe) * exp(g_mu0 * r1 * r2/e_resistivity(qvl, cl))
      case(16) ! non-reflection bc
         nv = zero
         nv = avto3(cf % vecn)
         !  (E-cBxn)xn
         mv = qvl % v(ieb:iee) - s_cspeed * acrossb(qvl % v(ibb:ibe), nv)
         mv = acrossb(mv, nv) * half
         qvr % v(ibb:ibe) = sum(qvl % v(ibb:ibe) * nv) * nv + mv / s_cspeed
         qvr % v(ieb:iee) = qvl % v(ieb:iee) - acrossb(qvl % v(ieb:iee), nv) + mv
         qvr % v(ibe) = zero
         qvr % v(ieb + 1:iee) = zero
      case(17)
         r1 = bc % vars(2); r2 = bc % vars(3)
         omega = two * pi * 1.e4_rfp
         ci = omega * s_elapsed_time
         current = bc % vars(1) * pi * r2**2 * sin(ci)
         btheta = bthetaf(r1, r2, current, r)
         qvr % v(ibb) = -p(2) * btheta / r
         qvr % v(ibb + 1) = p(1) * btheta / r
         current = bc % vars(1) * pi * r2**2 * cos(ci)
         current = current * two * pi * 1.e4_rfp
         btheta = bthetaf(r1, r2, current, r)
         qvr % v(iee) = qvl % v(iee) + half * btheta * sqrt(sum((cr % centp - cl % centp)**2))
         qvr % v(ifb) = qvl % v(ifb)
         ! try analytical solution
         !     current = s_cspeed**2 * g_mu0 * bc%vars(1) / omega
         !     ci = omega * r2 / s_cspeed
         !     CALL JY01A(ci,BJ0,r1,BJ1,r1,BY0,r1,BY1,r1)
         !     current = current  * bj1 / (by1 * bj0 - by0 * bj1)
         !     ci = omega * r / s_cspeed
         !     CALL JY01A(ci,BJ0,r1,BJ1,r1,BY0,r1,BY1,r1)
         !     r1 = atan(by0/bj0)
         !     qvr%v(iee) = current * ( bj0 * sin(omega * s_elapsed_time+r1) - &
         !                 by0 * cos(omega * s_elapsed_time+r1))
      case(21) ! puff case  db/dn=0;E_T=0
         qvr % v(ibb:ibe) = qvl % v(ibb:ibe)
         nv = avto3(cf % vecn)
         qvl % v(ieb:iee) = sum(qvl % v(ieb:iee) * nv) * nv
      case(22)
         qvr % v(ieb:iee) = qvl % v(ieb:iee) ! dE/dn=0
         current = 4.e6_rfp * sin(1.16e7_rfp * s_elapsed_time)
         ci = g_mu0 * current * half / pi
         qvr % v(ibb:ibb + 1) = zero
         qvr % v(ibe) = -ci / r
      end select
      qvr % v(ifb) = zero
      qvr % v(ife) = zero
      if (it /= 16) qvr % v(ibb:iee) = qvr % v(ibb:iee) * two - qvl % v(ibb:iee)
      cr % jv = cl % jv
   end subroutine mhd_boundary

   function current_vs_time(c)result(c1)
      real(rfp) :: c, c1
      c1 = c
      return
      c1 = min(one, s_elapsed_time / 1.e-6_rfp) * c
   end function current_vs_time

   subroutine eb_userdefine(cf, bc, it)
      implicit none
      type(bc_type), pointer :: bc
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      real(rfp) :: r, r1, r2, ci, p(ndim), btheta, current, j0
      integer :: it
      cl => cf % left_cell
      cr => cf % right_cell
      qvl => cl % qv
      qvr => cr % qv
      p = cf % centp
      j0 = 1.e5_rfp
      r = sqrt(sum(p(:2) * p(:2)))
      if (ndim == 2) r = cf % centp(s_iaxis)
      if (r < 0.2_rfp) then
         current = j0 * r**2/0.2**2
         !    qvr%v(ieb) = j0 / (0.04_rfp * pi) * e_resistivity(cl%qv,cl)
      else
         current = j0
      end if
      btheta = g_mu0 * current * half / pi / r
      qvr % v(ibb) = zero
      qvr % v(ibb + 1) = zero
      qvr % v(ibe) = btheta
      qvr % v(ieb:iee) = qvl % v(ieb:iee)
   end subroutine eb_userdefine

   function ojacob(cf)result(a) !vecn,qv,itype)result(a)
      type(face), pointer :: cf
      type(matrix) :: a
      integer :: it, i
      it = cf % itype
      a % e = zero
      if (it >= partition_face_no) return
      select case(bc(it) % igrp) ! boundary type
      case(inlet_type)
         call inlet_jacob(cf, a)
      case(outlet_type)
         call outlet_jacob(cf, a)
      case(farfield_type)
         call farfield_jacob(cf, a)
      case(wall_type)
         call wall_jacob(cf, a)
      case(geom_topology)
         call geom_jacob(cf, a)
      case(MHD_type)
         call mhd_jacob(bc(it) % itype, cf, a)
      case default
         print *, ' Donot know the boundary type', bc(it) % itype, it, bc(it) % igrp
      end select
      !
      ! Radiation
      !   if(nrad > 0) call radiation_jacob(cf,a)

      ! EB boundary condition
      it = bc(it) % itype / 1000
      if (neqm > 0.and.neqf > 0) call mhd_jacob(it, cf, a)
   end function ojacob

   subroutine radiation_jacob(cf, a)
      implicit none
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      real(rfp) :: dn, tw, alpha
      type(matrix) :: a
      integer :: it, i
      cr => cf % right_cell
      cl => cf % left_cell
      qvl => cl % qv
      qvr => cr % qv
      if (nrad == 0) return
      dn = sum(cf % vecn * (cr % centp - cl % centp))
      alpha = two * two / (three * absorption_coef(qvl) * dn)
      tw = half * (qvr % v(ien) + qvl % v(ien))
      a % e(irad, irad) = -(one - alpha)/(one + alpha)
      !    a%e(irad,ien)  = 16._rfp * stephan_boltzmann * tw**3
      select case(bc(cf % itype) % igrp) ! boundary type
      case(inlet_type, outlet_type, farfield_type)
         a % e(irad, irad) = one
      case(wall_type)
         it = mod(bc(cf % itype) % itype, 1000)
         select case(it)
         case(0, 7)
            a % e(irad, irad) = one
         case default
            !    qvr%v(irad) = 8._rfp * stephan_boltzmann * tw**4 - qvl%v(irad)
            alpha = two * two / (three * 100._rfp * dn)
            a % e(irad, irad) = -(one - alpha)/(one + alpha)
            !    a%e(irad,ien)  = 16._rfp * stephan_boltzmann * tw**3

         end select
      end select
   end subroutine radiation_jacob

   subroutine mhd_jacob(it, cf, a)
      implicit none
      type(matrix) :: a
      type(face), pointer :: cf
      type(cell), pointer :: cr, cl
      type(vector), pointer :: qvr, qvl
      integer :: i, it
      if (neqf == 0) a % e = zero
      call put_diag(a % e, -one, ibb, iee)
      select case(it)
      case(0)
         call put_diag(a % e, one, ibb, iee)
         a % e(ibb, ibb) = -one
         a % e(ibb + 1, ibb + 1) = -one
         a % e(iee, iee) = -one
      case(1)
         call put_diag(a % e, one, ibb, iee)
         a % e(ibb, ibb) = -one
         a % e(ibb + 1, ibb + 1) = -one
      case(2) ! free boundary condition
         call put_diag(a % e, one, ibb, iee)
      case(20) ! free boundary condition
         call put_diag(a % e, one, ibb, ibe)
         call put_diag(a % e, zero, ieb, iee)
      case(3)
         call put_diag(a % e, one, ieb, iee)
         !     a%e(ibb,ibe) = zero
      case(4) ! r = 0  centrline
         a % e(ieb, ieb) = one
      case(5)
         call put_diag(a % e, one, ibb, ibe)
      case(6) ! dielectric boundaru condition
         !     a%e(ibe,ibe) = zero
      case(7)
         call put_diag(a % e, one, ieb, iee)
      case(8)
         call put_diag(a % e, one, ibb, iee)
         a % e(ibe, ibe) = -one
      case(9, 12, 13, 23)
         call put_diag(a % e, one, ieb, iee)
      case(11)
         call put_diag(a % e, one, ibb, ibe)
      case(16)
         call put_diag(a % e, one, ibb, iee)
         a % e(ibe, ibe) = zero
         a % e(ieb + 1:iee, ieb + 1:iee) = zero
      case(17)
         call put_diag(a % e, zero, ibb, iee)
         call put_diag(a % e, one, iee, iee)
         call put_diag(a % e, one, ifb, ifb)
      case(21)
         call put_diag(a % e, one, ibb, ibe)
      case(22)
         call put_diag(a % e, zero, ibb, ibe)
      end select
   end subroutine mhd_jacob

   subroutine inlet_jacob(cf, a)
      implicit none
      type(face), pointer :: cf
      type(matrix) :: a
      type(vector), pointer :: qv, qvr
      real(rfp) :: vel(ndim), vd(ndim), u, rho, vecn(ndim), v(3)
      real(rfp) :: vr(ndim), vr2, sos
      integer :: i, it
      qv => cf % left_cell % qv
      qvr => cf % right_cell % qv
      vr = qvr % v(imb:ime)
      vr2 = sum(vr**2)
      sos = sound_speed(qv)

      a % e = zero
      it = mod(bc(cf % itype) % itype, 1000)
      select case(it)
      case(0, 1, 8) ! Specify H^0 and S
         rho = rhof(qv)
         a % e(ico, imb:ime) = -rho * qv % v(imb:ime)
         a % e(ien, imb:ime) = -(one - rho * hpf(qv)) / htf(qv) * qv % v(imb:ime)
         !    
         vel = qv % v(imb:ime)
         vel = vel / (sqrt(sum(vel * vel)) + mytiny)
         vd = bc(cf % itype) % vars(neq4:neq3d)
         if (it == 1) vd = -cf % vecn
         do i = imb, ime
            a % e(i, imb:ime) = vel * vd(i - imb + 1)
         end do
      case(2, 20) ! specify velocity T
         if(vr2 <= sos) then
            a % e(ico, ico) = one
         end if
      case(5) ! specify velocity T
         if(vr2 <= sos) then
            a % e(ico, ico) = one
         end if
      case(3, 30) ! mass flow & T
         vecn = -cf % vecn
         rho = rhof(cf % right_cell % qv)
         u = -rhopf(qv) * sum(qv % v(imb:ime) * vecn) / rho
         vel = u * vecn
         a % e(ico, ico) = one
         a % e(imb:ime, ico) = vel
      case(6)
         vel = qv % v(imb:ime)
         rho = rhof(qv)
         a % e(ico, imb:ime) = -rho * vel
         a % e(ien, imb:ime) = -(one - rho * hpf(qv)) / htf(qv) * vel
         return
         if (bc(cf % itype) % n == size(bc(cf % itype) % vars)) then
            bc(cf % itype) % n = 0
         else
            bc(cf % itype) % n = bc(cf % itype) % n + neq
         end if
         i = bc(cf % itype) % n
         vd(:ndim) = bc(cf % itype) % vars(i + imb:i + ime)
         do i = imb, ime
            a % e(i, imb:ime) = vel * vd(i - imb + 1)
         end do
      end select
   end subroutine inlet_jacob

   subroutine outlet_jacob(cf, a)
      implicit none
      type(face), pointer :: cf
      type(matrix) :: a
      type(vector), pointer :: qv
      real(rfp) :: v2, rho0, c0, p0, alpha, beta, Tl, pl, gm
      integer :: i, it, nindex

      qv => cf % left_cell % qv
      a % e = zero
      do i = 1, neqf
         a % e(i, i) = one
      end do
      it = mod(bc(cf % itype) % itype, 1000)
      if (it == 1) return ! full developed flow 
      v2 = sum(qv % v(imb:ime)**2)
      if (v2 > sound_speed(qv)) then
         !    cf%right_cell%qv%v(ico) = qv%v(ico)
         return ! supersonic
      end if

      select case(it)
      case default
         a % e(ico, ico) = zero
      end select
   end subroutine outlet_jacob

   subroutine farfield_jacob(cf, a)
      implicit none
      type(face), pointer :: cf
      type(cell), pointer :: cl
      type(matrix) :: a
      !
      real(rfp) :: un
      cl => cf % left_cell
      un = sum(cf % vecn * cl % qv % v(imb:ime))
      if (un < zero) then ! inlet
         call inlet_jacob(cf, a)
      else
         call outlet_jacob(cf, a)
      end if
   end subroutine farfield_jacob

   subroutine wall_jacob(cf, a)
      implicit none
      type(face), pointer :: cf
      type(matrix) :: a
      integer :: i, it
      real(rfp) :: tw, t
      it = mod(bc(cf % itype) % itype, 1000)
      call put_diag(a % e, one, neqf, neqf)
      !
      select case(it)
      case(0)
         call symmetry_jacob(cf, a)
      case(10)
         call symmetry_jacob(cf, a)
         a % e(ien, ien) = zero
      case(1, 3, 5, 7)
         if (nmeq > 0) then
            a % e(ico, ico) = one
            a % e(ien, ien) = one
            do i = imb, ime
               a % e(i, i) = -one
            end do
         end if
         a % e(ien, ien) = one
         if (s_iaux == 1) a % e(iaub, iaub) = -one
         if (s_ivis == 2) call put_diag(a % e, zero, ike, iom)
      case(2, 4, 6, 20)
         if (nmeq > 0) then
            a % e(ico, ico) = one
            do i = imb, ime
               a % e(i, i) = -one
            end do
         end if
         tw = bc(cf % itype) % vars(ndim + 1)
         t = cf % left_cell % qv % v(ien)
         a % e(ien, ien) = -(tw/t)**2
         if (s_iaux == 1) a % e(iaub, iaub) = -one
         if (s_ivis == 2) call put_diag(a % e, zero, ike, iom)
      case(8)
         call periodic_jacob(cf, a)
      end select
   end subroutine wall_jacob

   subroutine geom_jacob(cf, a)
      implicit none
      type(face), pointer :: cf
      type(matrix) :: a
      integer :: it
      it = mod(bc(cf % itype) % itype, 1000)
      select case(it)
      case(0)
         call symmetry_jacob(cf, a)
      case(1)
         call axisymmetry_jacob(cf, a)
      case(2)
         call periodic_jacob(cf, a)
      end select
   end subroutine geom_jacob

   subroutine symmetry_jacob(cf, a)
      implicit none
      type(face), pointer :: cf
      type(matrix) :: a
      integer :: i
      !   a%e = zero
      !   return
      do i = 1, neqf
         a % e(i, i) = one
      end do
      do i = imb, ime
         a % e(i, imb:ime) = -two * cf % vecn * cf % vecn(i - imb + 1)
         a % e(i, i) = one + a % e(i, i)
      end do
   end subroutine symmetry_jacob

   subroutine axisymmetry_jacob(cf, a)
      implicit none
      type(face), pointer :: cf
      type(cell), pointer :: cl, cr
      type(matrix) :: a
      real(rfp) :: r, r1(2), r2(2), h1, h2
      integer :: i
      a % e = zero
      return
      do i = 1, neqf
         a % e(i, i) = one
      end do
      cl => cf % left_cell
      cr => cf % right_cell
      !
      r = one / sqrt(cl % centp(2)**2 + cl % centp(3)**2)
      r1 = cl % centp(2:3) * r
      r2 = cr % centp(2:3) * r
      h1 = r1(1) * r2(1) + r1(2) * r2(2)
      h2 = r1(2) * r2(1) - r1(1) * r2(2)
      a % e(imb + 1, imb + 1) = h1
      a % e(imb + 1, ime) = h2
      a % e(ime, imb + 1) = -h2
      a % e(ime, ime) = h1
   end subroutine axisymmetry_jacob

   subroutine periodic_jacob(cf, a)
      implicit none
      type(face), pointer :: cf
      type(matrix) :: a
      integer :: i
      a % e = zero
   end subroutine periodic_jacob
   !
   !  implicit boundary condition
   !*****************************************************************
   subroutine implicit_boundary_condition(faces)
      implicit none
      type(face), pointer :: faces(:), cf
      type(cell), pointer :: cr, cl
      type(matrix) :: omega
      integer :: i
      do i = 1, size(faces)
         cf => faces(i)
         if (cf % itype == 0) exit
         if (cf % itype >= partition_face_no) cycle
         cl => cf % left_cell
         cr => cf % right_cell
         omega = ojacob(cf)
         cl % dm = cl % dm + cf % ajr * omega !%vecn,cr%qv,cf%itype)
         cf % ajl = omega
      end do
   end subroutine implicit_boundary_condition

   subroutine update_boundary(faces, interf)
      type(face), pointer :: faces(:), cf
      type(vector), pointer :: qvr, dqv
      type(itf) :: interf
      real(rfp) :: dh, ds, rho, t, hp, ht, dp, dt
      integer :: i

      do i = interf % nitf + 1, size(faces)
         cf => faces(i)
         if (cf % itype == 0) exit
         qvr => cf % right_cell % qv
         dqv => cf % left_cell % dqv
         qvr = qvr + cf % ajl * dqv
      end do
   end subroutine update_boundary

   subroutine update_geom_pinterface(faces, nodes, cells, interf)
      implicit none
      type(itf) :: interf
      type(face), pointer :: faces(:), cf
      type(node), pointer :: nodes(:)
      type(cell), pointer :: cells(:), cc, cg
      real(rfp), allocatable :: rb(:), sb(:)
      real(rfp) :: dd, dm
      integer :: i, ist, n, ib, j, k, ip, nd, ierr, st(mpi_status_size), nr, nr1, ir, it
      logical :: ismatch

      if (s_iperiod <= 0) return ! no need
      nr = 1 + ndim
      nr1 = nr - 1
      !
      call mpi_comm_rank(mpi_comm_world, ip, ierr)
      !
      ib = 0; ir = 0
      do i = 1, interf % nc
         nd = nr * interf % ncp(i)
         allocate(sb(nd))
         nd = 1
         do j = 1, interf % ncp(i)
            ib = ib + 1
            cc => interf % scell(ib) % to_cell
            ismatch = .false.
            ! find periodic face
            do k = 1, size(cc % sface)
               cf => cc % sface(k) % to_face
               it = cf % itype
               if (it == 0 .or. it >= partition_face_no) cycle

               if(associated(interf % sface)) then
                  if(associated(cf, interf % sface(ib) % to_face)) then
                     ismatch = .true.; exit
                  end if
               else
                  if(bc(it) % igrp == geom_topology .and. &
                     mod(bc(it) % itype, 1000) == 2) then
                     ismatch = .true.; exit
                  end if
               end if
            end do

            if(.not.ismatch) then
               print *, 'cannot find matching face in update_geom_pinterface'; stop
            end if

            sb(nd) = cc % vol
            call xy2rq(cf % centp, cc % centp - cf % centp, sb(nd + 1:nd + nr1))
            !     sb(nd+1:nd+nr1) = cc%centp - cf%centp   ! relative distance to face
            nd = nd + nr
         end do
         !
         nd = nr * interf % nnp(i)
         allocate(rb(nd))
         call mpi_sendrecv(sb, size(sb) * rfp, MPI_BYTE, interf % cp(i), ip, &
            rb, size(rb) * rfp, MPI_BYTE, interf % np(i), interf % np(i), mpi_comm_world, st, ierr)
         !
         nd = 1
         do j = 1, interf % nnp(i)
            ir = ir + 1
            cc => interf % pcell(ir)
            !
            cc % vol = rb(nd)
            cc % centp = rb(nd + 1:nd + nr1)
            ! need to be added face coordinate
            nd = nd + nr
         end do
         deallocate(rb, sb)
      end do
      ! corrected
      !
      do i = 1, size(faces)
         if (faces(i) % itype == 0) exit
         it = faces(i) % itype
         if (it == 0 .or. it >= partition_face_no) cycle
         if (bc(it) % igrp == geom_topology.and.mod(bc(it) % itype, 1000) == 2) then
            cc => faces(i) % right_cell
            call rq2xy(faces(i) % centp, cc % centp, cc % centp)
            cc % centp = cc % centp + faces(i) % centp
         end if
      end do
      !
   end subroutine update_geom_pinterface

   subroutine xy2rq(v1, v2, v3)
      real(rfp) :: v1(ndim), v2(ndim), v3(ndim)
      integer :: k1, k2
      real(rfp) :: p(2), u(2), rl, cos_sin(2)
      if (s_iperiod > 3) then
         v3 = v2
         return
      end if
      k1 = 0
      do k2 = 1, ndim
         if (k2 == s_iperiod) cycle
         k1 = k1 + 1
         p(k1) = v1(k2)
         u(k1) = v2(k2)
      end do
      rl = one / sqrt(sum(p**2))
      cos_sin = p * rl
      p(1) = cos_sin(1) * u(1) + cos_sin(2) * u(2)
      p(2) = -cos_sin(2) * u(1) + cos_sin(1) * u(2)
      k1 = 0
      do k2 = 1, ndim
         if (k2 == s_iperiod) then
            v3(k2) = v2(k2)
         else
            k1 = k1 + 1
            v3(k2) = p(k1)
         end if
      end do
   end subroutine xy2rq
   !
   subroutine rq2xy(v1, v2, v3)
      real(rfp) :: v1(ndim), v2(ndim), v3(ndim)
      integer :: k1, k2
      real(rfp) :: p(2), u(2), rl, cos_sin(2)
      if (s_iperiod > 3) then
         v3 = v2
         return
      end if
      k1 = 0
      do k2 = 1, ndim
         if (k2 == s_iperiod) cycle
         k1 = k1 + 1
         p(k1) = v1(k2)
         u(k1) = v2(k2)
      end do
      rl = one / sqrt(sum(p**2))
      cos_sin = p * rl
      p(1) = cos_sin(1) * u(1) - cos_sin(2) * u(2)
      p(2) = cos_sin(2) * u(1) + cos_sin(1) * u(2)
      k1 = 0
      do k2 = 1, ndim
         if (k2 == s_iperiod) then
            v3(k2) = v2(k2)
         else
            k1 = k1 + 1
            v3(k2) = p(k1)
         end if
      end do
   end subroutine rq2xy

   subroutine update_pinterface(nodes, cells, faces, interf)
      implicit none
      type(itf) :: interf
      type(face), pointer :: faces(:), cf
      type(node), pointer :: nodes(:)
      type(cell), pointer :: cells(:), cc
      real(rfp), allocatable :: sb(:), rb(:)
      integer :: i, ist, n, ib, j, ip, nd, ierr, st(mpi_status_size), neq1, ir, k1, k2, it
      real(rfp) :: rl, cos_sin(2), u(2), p(2)
      type(vector) :: qv
      ! 
      if (s_iperiod <= 0) return
      neq1 = neq - 1
      !
      call mpi_comm_rank(mpi_comm_world, ip, ierr)
      !
      ib = 0; ir = 0
      do i = 1, interf % nc
         nd = neq * interf % ncp(i)
         allocate(sb(nd))
         nd = 1
         do j = 1, interf % ncp(i)
            ib = ib + 1
            cc => interf % scell(ib) % to_cell
            qv = cc % qv
            !
            if (s_iperiod <= 3) then ! for cylindar wedge case
               k1 = 0
               do k2 = 1, ndim
                  if (k2 == s_iperiod) cycle
                  k1 = k1 + 1
                  p(k1) = cc % centp(k2) - s_pcenter(k2)
                  u(k1) = qv % v(imb + k2 - 1)
               end do
               rl = one / sqrt(sum(p**2))
               cos_sin = p * rl
               p(1) = cos_sin(1) * u(1) + cos_sin(2) * u(2)
               p(2) = -cos_sin(2) * u(1) + cos_sin(1) * u(2)
               k1 = 0
               do k2 = 1, ndim
                  if (k2 == s_iperiod) cycle
                  k1 = k1 + 1
                  qv % v(imb + k2 - 1) = p(k1)
               end do
               !
            end if
            sb(nd:nd + neq1) = qv % v
            nd = nd + neq
         end do
         !
         nd = neq * interf % nnp(i)
         allocate(rb(nd))
         if (ip == interf % np(i)) then
            rb = sb
         else
            call mpi_sendrecv(sb, size(sb) * rfp, MPI_BYTE, interf % cp(i), ip, &
               rb, size(rb) * rfp, MPI_BYTE, interf % np(i), interf % np(i), mpi_comm_world, st, ierr)
         end if
         !     
         nd = 1
         do j = 1, interf % nnp(i)
            ir = ir + 1
            cc => interf % pcell(ir)
            cc % qv % v = rb(nd:nd + neq1)
            nd = nd + neq
         end do
         deallocate(rb, sb)

      end do
      ! correction
      do i = 1, size(faces)
         cf => faces(i)
         it = cf % itype
         if (it == 0) exit ! no boundary faces, exit
         if (it >= partition_face_no) cycle
         if (bc(it) % igrp /= geom_topology) cycle
         if (mod(bc(it) % itype, 1000) /= 2) cycle
         call periodic_correction(cf, cells)
      end do
      !  
   end subroutine update_pinterface

   subroutine update_interface(faces, nodes, cells, interf)
      implicit none
      type(itf) :: interf
      type(face), pointer :: faces(:)
      type(node), pointer :: nodes(:)
      type(cell), pointer :: cells(:), cc
      real(rfp), allocatable :: sb(:), rb(:)
      integer :: i, ist, n, ib, j, ip, nd, ierr, st(mpi_status_size), neq1, ir, meq, op
      !
      call mpi_comm_rank(mpi_comm_world, ip, ierr)
      !
      ib = 0; ir = 0
      do i = 1, interf % nc
         op = interf % np(i)
         meq = pz(op) % neq
         !     
         nd = meq * interf % ncp(i)
         neq1 = meq - 1
         allocate(sb(nd))
         nd = 1
         do j = 1, interf % ncp(i)
            ib = ib + 1

            !pack information from scell into nd
            cc => interf % scell(ib) % to_cell
            sb(nd:nd + neq1) = cc % qv % v(pz(op) % loc)
            nd = nd + meq
         end do
         !
         nd = meq * interf % nnp(i)
         allocate(rb(nd))
         call mpi_sendrecv(sb, size(sb) * rfp, MPI_BYTE, interf % cp(i), ip, &
            rb, size(rb) * rfp, MPI_BYTE, interf % np(i), interf % np(i), mpi_comm_world, st, ierr)
         nd = 1
         do j = 1, interf % nnp(i)
            ir = ir + 1

            !receive information into rb and transfer it to pcell
            cc => interf % pcell(ir)
            cc % qv % v(pz(op) % loc) = rb(nd:nd + neq1)
            nd = nd + meq
         end do
         deallocate(rb, sb)
      end do
      !
      !
      !         ib = 0; ir = 0
      !         do i = 1, interf % nc
      !             op = interf % np(i)
      !             meq = 2
      !             if (neqm > 0) then
      !                 if (pz(op) % imeth == 3.or.pz(op) % imeth == 4) meq = meq + 3
      !             end if
      !             nd = meq * interf % ncp(i)
      !             !      neq1 = meq - 1
      !             allocate(sb(nd))
      !             nd = 1
      !             do j = 1, interf % ncp(i)
      !                 ib = ib + 1
      !                 cc => interf % scell(ib) % to_cell
      !                 sb(nd) = lamda(cc % qv, cc)
      !                 nd = nd + 1
      !                 sb(nd) = e_resistivity(cc % qv, cc)
      !                 nd = nd + 1
      !                 if (neqm > 0) then
      !                     if (pz(op) % imeth == 3.or.pz(op) % imeth == 4) then
      !                         sb(nd:nd + 2) = cc % jv
      !                         nd = nd + 3
      !                     end if
      !                 end if
      !             end do
      !             !
      !             nd = meq * interf % nnp(i)
      !             allocate(rb(nd))
      !             call mpi_sendrecv(sb, size(sb) * rfp, MPI_BYTE, interf % cp(i), ip, &
      !             rb, size(rb) * rfp, MPI_BYTE, interf % np(i), interf % np(i), mpi_comm_world, st, ierr)
      !             nd = 1
      !             do j = 1, interf % nnp(i)
      !                 ir = ir + 1
      !                 cc => interf % pcell(ir)
      !                 cc % rp = rb(nd)
      !                 nd = nd + 1
      !                 cc % er = rb(nd)
      !                 nd = nd + 1
      !                 if (neqm > 0) then
      !                     if (pz(op) % imeth == 3.or.pz(op) % imeth == 4) then
      !                         cc % jv = rb(nd:nd + 2)
      !                         nd = nd + 3
      !                     end if
      !                 end if
      !             end do
      !             deallocate(rb, sb)
      !         end do
      !
      call interface_bc(faces, interf)
      !  
   end subroutine update_interface

   subroutine interface_bc(faces, interf)
      type(itf) :: interf
      type(face), pointer :: faces(:), cf
      type(cell), pointer :: cc, cg
      type(vector), pointer :: qvr, qvl
      type(vector) :: qvf
      real(rfp) :: u, zk1, zk2, t1, t2, t, r1, r2, dd(ndim), vn(3), br(3), bl(3)
      integer :: i, j, im(1), ip, ierr
      call mpi_comm_rank(mpi_comm_world, ip, ierr)
      do i = 1, interf % nitf
         cf => faces(i)
         qvl => cf % left_cell % qv
         qvr => cf % right_cell % qv
         select case(cf % itype - partition_face_no)
         case(0)
            cycle
         case(1) ! wall + heat transfer
            !
            zk1 = lamda(qvl, cf % left_cell)
            zk2 = cf % right_cell % rp
            r1 = abs(dot_product(cf % vecn, cf % left_cell % centp - cf % centp))
            r2 = abs(dot_product(cf % vecn, cf % right_cell % centp - cf % centp))
            t1 = qvl % v(ien)
            t2 = qvr % v(ien)
            zk1 = zk1 / r1
            zk2 = zk2 / r2
            t = (zk1 * t1 + zk2 * t2)/ (zk1 + zk2)
            !       qvr%v(ien) = t  +  (r2 / r1) * (t - t1)
            !       if(qvr%v(ien) < zero ) qvr%v(ien) = t * t / t1
            !        qvr%v(ien) = t * (t / t1)  !**(r2/r1)
            !
            qvf = qvl
            qvf % v(ien) = t
            !
            if (nmeq > 0) then
               if (s_ivis == 0) then
                  u = sum(cf % vecn * qvl % v(imb:ime))
                  qvf % v(imb:ime) = qvl % v(imb:ime) - u * cf % vecn(:)
               else
                  qvf % v(imb:ime) = zero
               end if
               !         r1 = r2 / r1
               !         qvr%v(imb:ime) =  (one + r1) * qvf%v(imb:ime) - r1 * qvl%v(imb:ime)
               qvr % v(ico) = qvl % v(ico)
               qvr % v(imb:ime) = two * qvf % v(imb:ime) - qvl % v(imb:ime)
               if (neqf > ien) qvr % v(ien + 1:neqf) = qvl % v(ien + 1:neqf)
               if (s_ivis == 2) call no_wall_function(cf, zero)
            end if
         case(2)
            if (s_ivis == 2) then ! k-omega inlet condition
               u = max(sqrt(sum(qvr % v(imb:ime)**2)), 1.e-10_rfp)
               qvr % v(iom) = t_omegainf * u / b_lenref
               qvr % v(ike) = t_keinf * mumix(qvr) * qvr % v(iom) / rhof(qvr)
            end if
         case(3) ! maxwell-plasma
            if (neqf > 0) then
               qvf % v(:neqf) = qvl % v(:neqf)
               if (nmeq > 0) then
                  if (s_ivis == 0) then
                     u = sum(cf % vecn * qvl % v(imb:ime))
                     qvf % v(imb:ime) = qvl % v(imb:ime) - u * cf % vecn(:)
                  else
                     qvf % v(imb:ime) = zero
                  end if
                  qvr % v(:neqf) = two * qvf % v(:neqf) - qvl % v(:neqf)
                  ! for isothermal wall BC
                  !
                  qvr % v(ien) = 1000._rfp **2 / qvl % v(ien)
                  if (s_ivis == 2) call no_wall_function(cf, zero)
               end if
            end if
            !  treat for interface between two media
            !
            !      vn = avto3(cf%vecn)
            ! E field  ||Et||=0
            !   zk1 = e_resistivity(qvl,cf%left_cell)
            !   zk2 = e_resistivity(qvr,cf%right_cell)
            !   u = zk1 / zk2
            !!      t = sum((qvl%v(ieb:iee) - qvr%v(ieb:iee)) * vn)
            !      t = sum(qvr%v(ieb:iee) * vn) * ( u - one)
            !      qvr%v(ieb:iee) = qvr%v(ieb:iee) + t * vn
            ! B field  ||Bn||=0
            !!     t = sum((qvl%v(ibb:ibe) - qvr%v(ibb:ibe)) * vn)
            !!     qvr%v(ibb:ibe) = qvl%v(ibb:ibe) - t * vn
            !     br = qvr%v(ibb:ibe) - sum(qvr%v(ibb:ibe) * vn) * vn
            !     bl = qvl%v(ibb:ibe) - sum(qvl%v(ibb:ibe) * vn) * vn
            !     br = bl + two * one / ( one + u) * ( br - bl)
            !!     br = (br + u * bl ) / ( one + u) 
            !     qvr%v(ibb:ibe) = br + sum(qvr%v(ibb:ibe)*vn) * vn
            !     cc => cf%right_cell
            !     cc%gradient(ibb:ife,:) = zero
         case(4) ! fluid - plasma
            !       qvr%v(ibb:ibe) = qvl%v(ibb:ibe)
            qvr % v(ibb:ibe) = zero
            qvr % v(ieb:iee) = qvl % v(ieb:iee)
            !       qvr%v(ieb)   = zero
            qvr % v(ifb:ife) = zero
         end select
         ! 
         !    if(neqf > 0) then
         !    cc => cf%right_cell
         !    dd = cf%centp-cc%centp
         !    im = maxloc(abs(dd))
         !    cc%gradient(:neqf,im(1)) = (qvf%v(:neqf) - qvr%v(:neqf)) / dd(im(1))
         !    if(s_ivis == 2) cc%gradient(ike:iom,:) = zero
         !    end if
         !
      end do
   end subroutine interface_bc

   subroutine update_interface_dqv(nodes, cells, interf)
      implicit none
      type(itf) :: interf
      type(node), pointer :: nodes(:)
      type(cell), pointer :: cells(:), cc
      real(rfp), allocatable :: sb(:), rb(:)
      integer :: i, ist, n, ib, j, ip, nd, ierr, st(mpi_status_size), neq1, ir, meq, op
      !
      call mpi_comm_rank(mpi_comm_world, ip, ierr)
      !
      !
      ib = 0; ir = 0
      do i = 1, interf % nc
         op = interf % np(i)
         meq = pz(op) % neq
         nd = meq * interf % ncp(i)
         neq1 = meq - 1
         allocate(sb(nd))
         nd = 1
         do j = 1, interf % ncp(i)
            ib = ib + 1
            cc => interf % scell(ib) % to_cell
            sb(nd:nd + neq1) = cc % dqv % v(pz(op) % loc)
            nd = nd + meq
         end do
         !
         nd = meq * interf % nnp(i)
         allocate(rb(nd))
         call mpi_sendrecv(sb, size(sb) * rfp, MPI_BYTE, interf % cp(i), ip, &
            rb, size(rb) * rfp, MPI_BYTE, interf % np(i), interf % np(i), mpi_comm_world, st, ierr)
         nd = 1
         do j = 1, interf % nnp(i)
            ir = ir + 1
            cc => interf % pcell(ir)
            cc % dqv % v(pz(op) % loc) = rb(nd:nd + neq1)
            nd = nd + meq
         end do
         deallocate(rb, sb)
      end do
   end subroutine update_interface_dqv

   subroutine update_geom_interface(faces, nodes, cells, interf, bnodes)
      implicit none
      type(itf) :: interf
      type(face), pointer :: faces(:), cf
      type(node), pointer :: nodes(:)
      type(bnode), pointer :: bnodes(:)
      type(cell), pointer :: cells(:), cc, cg
      real(rfp), allocatable :: rb(:), sb(:)
      real(rfp) :: dd, dm
      integer :: i, ist, n, ib, j, ip, nd, ierr, st(mpi_status_size), nr, nr1, ir
      nr = 1 + ndim
      nr1 = nr - 1
      !
      call mpi_comm_rank(mpi_comm_world, ip, ierr)
      !
      ib = 0; ir = 0
      do i = 1, interf % nc
         nd = nr * interf % ncp(i)
         allocate(sb(nd))
         nd = 1
         do j = 1, interf % ncp(i)
            ib = ib + 1
            cc => interf % scell(ib) % to_cell
            sb(nd) = cc % vol
            sb(nd + 1:nd + nr1) = cc % centp
            nd = nd + nr
         end do
         !
         nd = nr * interf % nnp(i)
         allocate(rb(nd))
         call mpi_sendrecv(sb, size(sb) * rfp, MPI_BYTE, interf % cp(i), ip, &
            rb, size(rb) * rfp, MPI_BYTE, interf % np(i), interf % np(i), mpi_comm_world, st, ierr)
         !
         nd = 1
         do j = 1, interf % nnp(i)
            ir = ir + 1
            cc => interf % pcell(ir)
            !
            cc % vol = rb(nd)
            cc % centp = rb(nd + 1:nd + nr1)
            cc % dqv % v = zero
            cc % qv % v = zero
            cc % itype = partition_face_no ! default
            if (s_idt > 0) allocate(cc % qvn(s_idt))
            !
            cc % srf = pz(interf % np(i)) % imeth ! from zone
            !
            nd = nd + nr
         end do
         deallocate(rb, sb)
      end do
      !
      !  correction for corner between different zones
      !
      !   do i = 1, size(interf%pcell)
      !    cc => interf%pcell(i)
      !    if(cc%itype /= 0) cycle  ! need to be worried
      !     dm = 1.e30_rfp
      !    do j = 1, interf%nitf
      !     cf => faces(j)
      !     dd = sum((cf%centp-cc%centp)**2)
      !     if(dd < dm) then
      !      cc%itype = -j
      !      dm = dd
      !     end if
      !    end do
      !   end do 
      !
      select case(pz(ip) % imeth)
      case(0) ! n-s
         do i = 1, interf % nitf
            cf => faces(i)
            ir = cf % right_cell % srf
            if (ir == 2) cf % itype = cf % itype + 1 ! heat
            if (ir == 4) cf % itype = cf % itype + 4 ! plasma
         end do
      case(1) ! k-omega
         do i = 1, interf % nitf
            cf => faces(i)
            ir = cf % right_cell % srf
            if (ir == 2) cf % itype = cf % itype + 1 ! heat
            if (ir == 0) cf % itype = cf % itype + 2
            if (ir == 4) cf % itype = cf % itype + 4 ! plasma
         end do
      case(2) ! heat
         do i = 1, interf % nitf
            cf => faces(i)
            ir = cf % right_cell % srf
            cf % itype = cf % itype + 1 ! heat
         end do
      case(3) ! maxwell
         do i = 1, interf % nitf
            cf => faces(i)
            ir = cf % right_cell % srf
            cf % itype = cf % itype + 3 ! maxwell-plasma
            if (ir /= 3.and.ir /= 4) print *, 'not apply this type interfacei', 3, ir
         end do
      case(4) ! maxwell
         do i = 1, interf % nitf
            cf => faces(i)
            ir = cf % right_cell % srf
            if (ir == 3) cf % itype = cf % itype + 3 ! maxwell-plasma
            if (ir == 0) cf % itype = cf % itype + 4 ! plasma-fluid
            if (ir == 1) cf % itype = cf % itype + 4 ! plasma-fluid
         end do
      end select
      !  collect the nodes on the interface
      !
      call collect_bnode(faces, cells, nodes, bnodes, interf)
      ! distinguish the corner cells
      do i = 1, size(interf % pcell)
         cc => interf % pcell(i)
         cc % srf = -one
      end do
      do i = 1, interf % nitf
         cf => faces(i)
         cf % right_cell % srf = one
      end do
   end subroutine update_geom_interface

   subroutine collect_bnode(faces, cells, nodes, bnodes, interf)
      implicit none
      type(itf) :: interf
      type(face), pointer :: faces(:), cf
      type(node), pointer :: nodes(:)
      type(bnode), pointer :: bnodes(:)
      type(cell), pointer :: cells(:), cc
      type(neighbour_cell) :: scell(100)
      integer :: nd(size(nodes))
      integer, pointer :: n2b, c2n(:)
      integer :: i, j, n
      !
      nd = 0
      do i = 1, interf % nitf
         cf => faces(i)
         select case(cf % itype - partition_face_no)
         case(1, 3) ! heat interface
            nd(cf % f2n) = 1 ! marked as interface nodes
         end select
      end do
      n = sum(nd)
      allocate(bnodes(n))
      j = 0
      do i = 1, size(nd)
         if (nd(i) > 0) then
            j = j + 1
            nd(i) = j
            bnodes(j) % tn => nodes(i)
            allocate(bnodes(j) % scell(100)) ! allocate large enough
         end if
      end do
      !
      do i = 1, interf % nitf
         cf => faces(i)
         select case(cf % itype - partition_face_no)
         case(1, 3) ! heat interface
            bnodes(nd(cf % f2n)) % itype = cf % itype - partition_face_no
         end select
      end do
      !
      bnodes % n = 0
      do i = 1, size(cells)
         cc => cells(i)
         c2n => cc % c2n
         do j = 1, size(c2n)
            n = nd(c2n(j))
            if (n > 0) then
               n2b => bnodes(n) % n
               n2b = n2b + 1
               bnodes(n) % scell(n2b) % to_cell => cells(i)
            end if
         end do
      end do

      do i = 1, size(interf % pcell)
         cc => interf % pcell(i)
         c2n => cc % c2n
         do j = 1, size(c2n)
            n = nd(c2n(j))
            if (n > 0) then
               n2b => bnodes(n) % n
               n2b = n2b + 1
               bnodes(n) % scell(n2b) % to_cell => interf % pcell(i)
            end if
         end do
      end do

      do i = interf % nitf + 1, size(faces)
         if (faces(i) % itype == 0) exit
         cf => faces(i)
         cc => cf % right_cell
         c2n => cf % f2n
         do j = 1, size(c2n)
            n = nd(c2n(j))
            if (n > 0) then
               n2b => bnodes(n) % n
               n2b = n2b + 1
               bnodes(n) % scell(n2b) % to_cell => cf % right_cell
            end if
         end do
      end do
      !
      do i = 1, size(bnodes)
         if (bnodes(i) % n > 100) print *, 'You need to increase number of scell in collect_node subroutine!!'
         scell = bnodes(i) % scell
         deallocate(bnodes(i) % scell)
         allocate(bnodes(i) % scell(bnodes(i) % n))
         bnodes(i) % scell = scell(:bnodes(i) % n)
      end do
   end subroutine collect_bnode

   subroutine update_gradient(nodes, cells, interf)
      implicit none
      type(itf) :: interf
      type(node), pointer :: nodes(:)
      type(cell), pointer :: cells(:), cc
      real(rfp), allocatable :: rb(:), sb(:)
      integer :: i, ist, n, ib, j, ip, nd, ierr, st(mpi_status_size), nblock, nrec, ir
      integer :: meq, op
      !  if(s_ivis < 2 .and. s_ialg == 0) return
      nrec = s_nrec
      !
      call mpi_comm_rank(mpi_comm_world, ip, ierr)
      !
      ib = 0; ir = 0
      do i = 1, interf % nc
         op = interf % np(i)
         meq = pz(op) % neq
         nblock = meq * nrec
         nd = nblock * interf % ncp(i)
         allocate(sb(nd))
         nd = 1
         do j = 1, interf % ncp(i)
            ib = ib + 1
            cc => interf % scell(ib) % to_cell
            sb(nd:nd + nblock - 1) = &
               reshape(cc % gradient(pz(op) % loc,:), shape = (/nblock/))
            nd = nd + nblock
         end do
         !
         nd = nblock * interf % nnp(i)
         allocate(rb(nd))

         call mpi_sendrecv(sb, size(sb) * rfp, MPI_BYTE, interf % cp(i), ip, &
            rb, size(rb) * rfp, MPI_BYTE, interf % np(i), interf % np(i), mpi_comm_world, st, ierr)
         nd = 1
         do j = 1, interf % nnp(i)
            ir = ir + 1
            cc => interf % pcell(ir)
            cc % gradient(pz(op) % loc,:) = reshape(rb(nd:nd + nblock - 1), (/meq, nrec/))
            nd = nd + nblock
         end do
         deallocate(rb, sb)
      end do
   end subroutine update_gradient

   subroutine setup_periodic_face(faces, cells, interf)
      implicit none
      type(itf) :: interf
      type(face), pointer :: faces(:), cf
      type(cell), pointer :: cells(:), cc, rc
      integer :: ip, ierr, ib, ir
      integer :: st(mpi_status_size)
      integer :: i, j, nd, nr, k, it, nscell
      real(rfp), allocatable :: rb(:), sb(:)
      real(rfp) :: vecn(ndim), cross(ndim)
      logical :: ismatch

      if(s_iperiod /= 4) return
      nscell = size(interf % scell)
      allocate(interf % sface(nscell))

      nr = ndim

      call mpi_comm_rank(mpi_comm_world, ip, ierr)

      ib = 0; ir = 0
      do i = 1, interf % nn
         nd = nr * interf % nnp(i)
         allocate(sb(nd))
         nd = 1
         do j = 1, interf % nnp(i)
            ib = ib + 1
            cc => interf % pcell(ib)

            ismatch = .false.
            ! find periodic face
            do k = 1, size(faces)
               cf => faces(k)
               it = cf % itype
               if(it == 0) exit
               if(it >= partition_face_no) cycle
               if(bc(it) % igrp == geom_topology .and. &
                  mod(bc(it) % itype, 1000) == 2) then
                  rc => cf % right_cell
                  if(associated(rc, cc)) then
                     ismatch = .true.
                     exit
                  end if
               end if
            end do
            if(.not.ismatch) then
               print *, 'cannot find matching face for pcell'; stop
            end if

            sb(nd:nd+nr-1) = cf % vecn(:)
            nd = nd + nr
         end do
         !
         nd = nr * interf % ncp(i)
         allocate(rb(nd))
         call mpi_sendrecv(sb, size(sb) * rfp, MPI_BYTE, interf % np(i), ip, &
            rb, size(rb) * rfp, MPI_BYTE, interf % cp(i), interf % cp(i), mpi_comm_world, st, ierr)
         !
         nd = 1
         do j = 1, interf % ncp(i)
            ir = ir + 1
            cc => interf % scell(ir) % to_cell
            !
            vecn(:) = rb(nd:nd+nr-1)
            nd = nd + nr

            ! find matching face
            ismatch = .false.
            do k = 1, size(cc % sface)
               cf => cc % sface(k) % to_face
               it = cf % itype
               if(it == 0 .or. it >= partition_face_no) cycle
               if(bc(it) % igrp == geom_topology .and. &
                  mod(bc(it) % itype, 1000) == 2) then
                  if(check_parallel(cf % vecn, vecn)) then
                     ismatch = .true.; exit
                  end if
               end if
            end do

            if(ismatch) then
               interf % sface(ir) % to_face => cf
            else
               print *, 'cannot find matching face for scell'; stop
            end if
         end do
         deallocate(rb, sb)
      end do
   end subroutine setup_periodic_face

end module gems_bound
