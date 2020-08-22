module gems_fv
   use gems_library
   use gems_state
   use gems_jacob
   use gems_disc
   use gems_precon
   use gems_data
   use gems_turb
   implicit none
contains
   !
   !********************************************************
   subroutine residual(cells, faces, nodes)
      !********************************************************
      type(node), intent(in) :: nodes(:)
      type(cell), pointer :: cells(:), dbcell
      type(face), pointer :: faces(:), current_face, mon_face
      ! local variables
      !
      type(matrix), pointer :: dml, dmr
      ! variables at face
      type(vector) :: faceflux
      real(kind = rfp) :: vecn(ndim)
      type(matrix) :: aav
      !  
      ! variables at left and right cell  
      type(cell), pointer :: cl, cr
      type(vector) :: qvl, qvr, dqvlr
      type(vector) :: evl, evr
      type(matrix) :: ajl, ajr
      !  type(matrix),pointer::ajl,ajr
      real(kind = rfp) :: tl, tr, rrr, dp
      integer, pointer :: c2n(:), f2n(:)
      !
      ! variables for computation
      integer :: i, it, n

      !only heat conduction
      !still need to modify for level-set
      if (nmeq == 0.and.neqm == 0) then
         call vis_spectrum(faces)
         return
      end if

      do i = 1, size(faces)

         current_face => faces(i)
         it = current_face % itype
         cl => current_face % left_cell
         cr => current_face % right_cell
         qvl = cl % qv
         qvr = cr % qv
         vecn = current_face % vecn

         !
         if (abs(s_ialg) >= 1) then
            !     if(it > 0 .and. it < partition_face_no) then
            !     qvr = (qvr + qvl) * half
            !     else	
            qvr = qvr + qv_from_highorder(cr, current_face, nodes)
            qvl = qvl + qv_from_highorder(cl, current_face, nodes)
            !     if(qvl%v(ien) <= zero) qvl = cl%qv
            !     if(qvr%v(ien) <= zero) qvr = cr%qv
         end if

         ajl % e = zero; ajr % e = zero;
         !  force wall boundary conditions
         !    if(it /= 0 .and. it < partition_face_no .and. bc(it)%igrp == wall_type) then
         !      qvl = (qvl + qvr) * half
         !      faceflux = fluxv(vecn,qvl,cl)
         !      ajl = ajacobian(vecn * half,qvl,cl)
         !      ajr = ajacobian(vecn * half,qvr,cr)
         !    else
         !     aav = a_average_jacobian(current_face,cl,cr,qvl,qvr,aav1)
         select case(s_ischeme)
         case(1)
            call hlle_flux(current_face, cl, cr, qvl, qvr, faceflux, aav, tl, tr, ajl, ajr)
         case(2)
            call hllc_flux(current_face, cl, cr, qvl, qvr, faceflux, aav, tl, tr, ajl, ajr)
         case(4)
            call rusanov_flux(current_face, cl, cr, qvl, qvr, faceflux, aav, tl, tr, ajl, ajr)
         case(6)
            call roe_flux(current_face, cl, cr, qvl, qvr, faceflux, aav, tl, tr, ajl, ajr)
         case default
            call facefluxscheme(current_face, cl, cr, qvl, qvr, faceflux, aav, tl, tr, ajl, ajr)

            if (s_ischeme /= 5) then
               !calculate the convective jacobian
               ajl = ajacobian(vecn * tl, qvl, cl, current_face % centp)
               ajr = ajacobian(vecn * tr, qvr, cr, current_face % centp)
               ! add the Roe's average matrix
               ajl = ajl + aav
               ajr = ajr - aav
            end if
         end select
         ! positive normal as pointer to right cell       
         !  LHS
         !
         dml => cl % dm
         dmr => cr % dm

         !assemble res, dm, ajl, and ajr
         cl % res = cl % res - faceflux ! -Res
         if (it == 0) then
            cr % res = cr % res + faceflux !-Res
         end if
         dml = dml + ajl
         ajl % e = -ajl % e
         dmr = dmr - ajr
         current_face % ajl = current_face % ajl + ajl
         current_face % ajr = current_face % ajr + ajr
      end do

      nullify(cl, cr, dmr, dml, current_face)
   end subroutine residual
   !
   subroutine vis_spectrum(faces)
      type(face) :: faces(:)
      type(cell), pointer :: cl, cr
      integer :: i
      real(rfp) :: a, cmax, zkl, zkr
      do i = 1, size(faces)
         cl => faces(i) % left_cell
         cr => faces(i) % right_cell
         zkl = lamda(cl % qv, cl)
         zkr = lamda(cr % qv, cr)
         a = faces(i) % area
         cmax = max(zkl, zkr) / abs(sum(faces(i) % vecn * (cr % centp - cl % centp)))
         cmax = cmax / rhof(cl % qv) / htf(cl % qv) ! thermal diffusivity
         if (s_cfl > zero) then
            cl % srf = cl % srf + cmax * a
            cr % srf = cr % srf + cmax * a
         else
            cl % srf = max(cl % srf, cmax)
            cr % srf = max(cr % srf, cmax)
         end if
         !
         !   cl%srf = cl%srf + cmax * a
         !   cr%srf = cr%srf + cmax * a
      end do
   end subroutine vis_spectrum

   subroutine roe_flux(cf, cl, cr, qvl, qvr, flux, aa, tl, tr, ajl, ajr)
      implicit none
      type(cell), pointer :: cl, cr
      type(face) :: cf
      !
      !  local variables
      real(kind = rfp) :: rhol, rhor
      real(kind = rfp) :: vecn(ndim), tl, tr, a, sl, sr, sm
      type(vector) :: qvl, qvr, qva, evl, evr, qcl, qcr, flux, dsmdq
      real(rfp) :: sqrhol, sqrhor, c1, c2, h0l, h0r, cmax, lamda_l, lamda_r, delta
      real(rfp) :: fac, fac1, h0av, rhoav, rhoppav, alpha
      real(rfp) :: velav(ndim), unav, dl, unl, unr, pl, pr, rhos, ps
      type(matrix) :: aa, zm, zminv, gj, ajl, ajr
      integer :: ic, i
      !
      vecn = cf % vecn
      a = cf % area
      evl = fluxv(vecn, qvl, cl, cf % centp)
      evr = fluxv(vecn, qvr, cr, cf % centp)
      qcl = conservative_qv(qvl)
      qcr = conservative_qv(qvr)
      !   faceflux = evl * tl + evr * tr  - aav * ( qvr - qvl)
      !
      if (neqf > 0) then
         ! special treat for inlet boundary condition
         !
         !  if(bc(cf%itype)%igrp == inlet_type) then
         !    tl = zero; tr = a; aa%e = zero;
         !    flux = evr
         !    return
         !  end if

         rhol = rhof(qvl)
         rhor = rhof(qvr)
         sqrhol = sqrt(rhol)
         sqrhor = sqrt(rhor)
         rhoav = sqrhol * sqrhor
         fac = sqrhol / (sqrhol + sqrhor)
         fac1 = one - fac
         !   qva%v(ico) = fac * qvl%v(ico) + fac1 * qvr%v(ico)
         !   fac = (rhol / rhoav) * fac
         !   fac1 = (rhor / rhoav) * fac1
         qva % v(:neqf) = qvl % v(:neqf) * fac + qvr % v(:neqf) * fac1
         !   qva%v(imb:ime) = (qvl%v(imb:ime) * rhol * fac + &
         !                     qvr%v(imb:ime) * rhor * fac1 ) / rhoav
         if(neeq /= 0) then
            if (qva % v(ien) <= zero) then
               print *, qvl, qvr
               print *, cf % itype, cf % area
               print *, rhol, rhor, fac
               do i = 1, size(cl % sface)
                  print *, cl % sface(i) % to_face % itype
               end do
            end if
         end if
         !
         h0l = h0f(qvl)
         h0r = h0f(qvr)
         !h0av = (fac * h0l /rhol + fac1 * h0r / rhor) * rhoav
         h0av = h0l * fac + h0r * fac1          
         call qvfromh0(rhoav, h0av, qva, qvl, qvr)

         !
         if (s_ischeme == -1) qva = (qvl + qvr) * (half)
         unav = dot_product(vecn, qva % v(imb:ime) - mesh_vel(cf % centp))
         !
         call precondition_face(s_ipre, qva, cf, cl, cr, rhoppav)
         call eigenvalues(s_ipre, c1, c2, vecn, qva, rhoppav, cf % centp)
         cmax = max(abs(c1), abs(c2))
         gj = gamma(s_ipre, qva, rhoppav)
      end if

      !
      !  if(neqm > 0.and.neqf == 0) cmax = s_cspeed
      !
      if (s_cfl > zero) then
         cl % srf = cl % srf + cmax * a
         cr % srf = cr % srf + cmax * a
      else
         cl % srf = max(cl % srf, cmax)
         cr % srf = max(cr % srf, cmax)
      end if
      !
      aa % e = zero
      if (neqf > 0) then
         call eigenvector(c1, c2, zm, zminv, cf % vecn, qva, rhoppav, cf % centp)
         !     zminv%e(:ndim,:neqf) = zminv%e(:ndim,:neqf) * abs(unav)
         !     zminv%e(  ime,:neqf) = zminv%e(  ime,:neqf) * abs(c1)
         !     zminv%e(  ien,:neqf) = zminv%e(  ien,:neqf) * abs(c2)

         !Wenda: replace the above 3 lines with the following 3 lines
         zminv % e(:ndim,:neqf) = zminv % e(:ndim,:neqf) * abs(unav)
         zminv % e(imb:ime,:neqf) = zminv % e(imb:ime,:neqf) * abs(c1)
         zminv % e(ien,:neqf) = zminv % e(ien,:neqf) * abs(c2)



         if (neqf > ien) zminv % e(ien + 1:neqf,:neqf) = zminv % e(ien + 1:neqf,:neqf) * abs(unav)
         !     if(s_ivis == 2) then ! special treat for k-omega
         !       zminv%e(ike:iom ,:neqf) = zminv%e(ike:iom,:neqf) * sqrt(sum(qva%v(imb:ime)**2))
         !       zminv%e(ike:iom ,:neqf) = zminv%e(ike:iom,:neqf) * cmax
         !       if(neqf > iom) zminv%e(iom+1:neqf ,:neqf) = zminv%e(iom+1:neqf,:neqf) * abs(unav)
         !     else if(neqf > ien) then
         !       zminv%e(ien+1:neqf ,:neqf) = zminv%e(ien+1:neqf,:neqf) * abs(unav)
         !     end if
         zm % e(:neqf,:neqf) = matmul(zm % e(:neqf,:neqf), zminv % e(:neqf,:neqf))
         aa % e(:neqf,:neqf) = matmul(gj % e(:neqf,:neqf), zm % e(:neqf,:neqf))
      end if
      !
      if (neqm > 0) then ! Maxwell
         call aaws(tl, sl, cl)
         call aaws(tr, sr, cr)
         !     tl = half * ( tl + tr)  !min(tl,sl)
         !     tr = half * ( sl + sr)  !min(tl,sl)
         tl = tl * tr / (tl + tr) !min(tl,sl)
         tr = sl * sr / (sl + sr) !min(tl,sl)
         sl = tl * tr
         if (s_cfl > zero) then
            cl % srfm = cl % srfm + sl * a
            cr % srfm = cr % srfm + sl * a
         else
            cl % srfm = max(cl % srfm, sl)
            cr % srfm = max(cr % srfm, sl)
         end if
         call meigenmatrix(tl, tr, cf % vecn, aa % e)
         !     qvr%v(ieb:iee) = qvr%v(ieb:iee) / e_resistivity(qvr,cr)
         !     qvl%v(ieb:iee) = qvl%v(ieb:iee) / e_resistivity(qvl,cl)
      end if
      aa = aa * (a * half)
      !     sl = abs(sum((cl%centp - cf%centp)*cf%vecn))
      !     sr = abs(sum((cr%centp - cf%centp)*cf%vecn))
      !     tl = a * sr /(sl+sr); tr = a * sl / ( sl + sr)
      tl = a * half; tr = a * half
      !    if(s_irealgas == 10) aa%e = zero
      flux = evl * tl + evr * tr - aa * (qvr - qvl)
      ajl = ajacobian(vecn * tl, qvl, cl, cf % centp)
      ajr = ajacobian(vecn * tr, qvr, cr, cf % centp)
      ajl = ajl + aa
      ajr = ajr - aa
      !
      return
      if (neqf > ien) then
         unl = sum(vecn * (qvl % v(imb:ime) - mesh_vel(cf % centp)))
         unr = sum(vecn * (qvr % v(imb:ime) - mesh_vel(cf % centp)))
         if (unav > 0) then ! pure upwind
            qva = rhoxf(qvl)
            do i = ien + 1, neqf
               flux % v(i) = rhol * unl * qvl % v(i) * a
               ajl % e(i,:neqf) = qva % v(:neqf) * unl * qvl % v(i) * a
               ajl % e(i, imb:ime) = ajl % e(i, imb:ime) + rhol * vecn * qvl % v(i) * a
               ajl % e(i, i) = ajl % e(i, i) + rhol * unl * a
               ajr % e(i,:neqf) = zero
            end do
         else
            qva = rhoxf(qvr)
            do i = ien + 1, neqf
               flux % v(i) = rhor * unr * qvr % v(i) * a
               ajr % e(i,:neqf) = qva % v(:neqf) * unr * qvr % v(i) * a
               ajr % e(i, imb:ime) = ajr % e(i, imb:ime) + rhor * vecn * qvr % v(i) * a
               ajr % e(i, i) = ajr % e(i, i) + rhor * unr * a
               ajl % e(i,:neqf) = zero
            end do
         end if
      end if
      return
      if (neqf > ien) then
         unl = sum(vecn * (qvl % v(imb:ime) - mesh_vel(cf % centp)))
         unr = sum(vecn * (qvr % v(imb:ime) - mesh_vel(cf % centp)))
         qcl = rhoxf(qvl) * ((unl - unav) / rhol * half * fac * a)
         qcl % v(imb:ime) = qcl % v(imb:ime) + vecn * fac * a
         qcr = rhoxf(qvr) * ((unr - unav) / rhor * half * fac1 * a)
         qcr % v(imb:ime) = qcr % v(imb:ime) + vecn * fac1 * a
         do i = ien + 1, neqf
            ajl % e(i,:neqf) = ajl % e(i,:neqf) - half * sum(gj % e(i,:neqf) * (qvr % v(:neqf) - qvl % v(:neqf))) * &
               abs(qcl % v(:neqf))
            ajr % e(i,:neqf) = ajr % e(i,:neqf) - half * sum(gj % e(i,:neqf) * (qvr % v(:neqf) - qvl % v(:neqf))) * &
               abs(qcr % v(:neqf))
         end do
      end if

   end subroutine roe_flux

   subroutine Rusanov_flux(cf, cl, cr, qvl, qvr, flux, aa, tl, tr, ajl, ajr)
      implicit none
      type(cell), pointer :: cl, cr
      type(face) :: cf
      !
      !  local variables
      real(kind = rfp) :: rhol, rhor
      real(kind = rfp) :: vecn(ndim), tl, tr, a, sl, sr, sm
      type(vector) :: qvl, qvr, qva, evl, evr, qcl, qcr, flux
      real(rfp) :: sqrhol, sqrhor, c1, c2, h0l, h0r, cmax, lamda_l, lamda_r, delta
      real(rfp) :: fac, fac1, h0av, rhoav, rhoppav, alpha
      real(rfp) :: velav(ndim), unav, dl, unl, unr, pl, pr, rhos, ps
      type(matrix) :: aa, zm, zminv, gj, ajl, ajr
      integer :: ic, i
      !
      vecn = cf % vecn
      a = cf % area
      evl = fluxv(vecn, qvl, cl, cf % centp)
      evr = fluxv(vecn, qvr, cr, cf % centp)
      qcl = conservative_qv(qvl)
      qcr = conservative_qv(qvr)
      !   faceflux = evl * tl + evr * tr  - aav * ( qvr - qvl)
      !
      if (neqf > 0) then
         ! special treat for inlet boundary condition
         !
         !  if(bc(cf%itype)%igrp == inlet_type) then
         !    tl = zero; tr = a; aa%e = zero;
         !    flux = evr
         !    return
         !  end if

         rhol = rhof(qvl)
         rhor = rhof(qvr)
         sqrhol = sqrt(rhol)
         sqrhor = sqrt(rhor)
         rhoav = sqrhol * sqrhor
         fac = sqrhol / (sqrhol + sqrhor)
         fac1 = one - fac
         qva % v(:neqf) = qvl % v(:neqf) * fac + qvr % v(:neqf) * fac1
         h0l = h0f(qvl)
         h0r = h0f(qvr)
         h0av = (fac * h0l /rhol + fac1 * h0r / rhor) * rhoav
         call qvfromh0(rhoav, h0av, qva, qvl, qvr)
         unav = dot_product(vecn, qva % v(imb:ime) - mesh_vel(cf % centp))
         call precondition_face(s_ipre, qva, cf, cl, cr, rhoppav)
         call eigenvalues(s_ipre, c1, c2, vecn, qva, rhoppav, cf % centp)
         cmax = max(abs(c1), abs(c2))
         gj = gamma(s_ipre, qva, rhoppav)
      end if
      !
      if (s_cfl > zero) then
         cl % srf = cl % srf + cmax * a
         cr % srf = cr % srf + cmax * a
      else
         cl % srf = max(cl % srf, cmax)
         cr % srf = max(cr % srf, cmax)
      end if
      ajl = ajacobian(vecn, qvl, cl, cf % centp)
      ajr = ajacobian(vecn, qvr, cr, cf % centp)
      if (neqf > 0) then
         tl = half * a; tr = half * a
         cmax = cmax * a * half
         aa % e(:ien,:neqf) = gj % e(:ien,:neqf) * cmax
         flux % v(:ien) = evl % v(:ien) * tl + evr % v(:ien) * tr - &
            matmul(aa % e(:ien,:neqf), qvr % v(:neqf) - qvl % v(:neqf))
         ajl % e(:ien,:neqf) = ajl % e(:ien,:neqf) * tl + aa % e(:ien,:neqf)
         ajr % e(:ien,:neqf) = ajr % e(:ien,:neqf) * tr - aa % e(:ien,:neqf)
         if (neqf > ien) then
            cmax = abs(unav * a) * half
            do i = ien + 1, neqf
               aa % e(i,:neqf) = gj % e(i,:neqf) * cmax
               flux % v(i) = evl % v(i) * tl + evr % v(i) * tr - &
                  sum(aa % e(i,:neqf)*(qvr % v(:neqf) - qvl % v(:neqf)))
               ajl % e(i,:neqf) = ajl % e(i,:neqf) * tl + aa % e(i,:neqf)
               ajr % e(i,:neqf) = ajr % e(i,:neqf) * tr - aa % e(i,:neqf)
            end do
         end if
      end if

      if (neqm > 0) then ! Maxwell
         call aaws(tl, sl, cl)
         call aaws(tr, sr, cr)
         !     tl = half * ( tl + tr)  !min(tl,sl)
         !     tr = half * ( sl + sr)  !min(tl,sl)
         tl = tl * tr / (tl + tr) !min(tl,sl)
         tr = sl * sr / (sl + sr) !min(tl,sl)
         sl = tl * tr
         if (s_cfl > zero) then
            cl % srfm = cl % srfm + sl * a
            cr % srfm = cr % srfm + sl * a
         else
            cl % srfm = max(cl % srfm, sl)
            cr % srfm = max(cr % srfm, sl)
         end if
         !     cl%Srfm = max(cl%srfm , sl)
         !     cr%srfm = max(cr%srfm , sl)
         call meigenmatrix(tl, tr, cf % vecn * (a * half), aa % e)
         tl = half * a; tr = half * a
         flux % v(ibb:ife) = evl % v(ibb:ife) * tl + evr % v(ibb:ife) * tr - &
            matmul(aa % e(ibb:ife, ibb:ife), qvr % v(ibb:ife) - qvl % v(ibb:ife))
         ajl % e(ibb:ife, ibb:ife) = ajl % e(ibb:ife, ibb:ife) * tl + aa % e(ibb:ife, ibb:ife)
         ajr % e(ibb:ife, ibb:ife) = ajr % e(ibb:ife, ibb:ife) * tr - aa % e(ibb:ife, ibb:ife)
      end if
      !
   end subroutine Rusanov_flux

   subroutine hlle_flux(cf, cl, cr, qvl, qvr, flux, aa, tl, tr, ajl, ajr)
      implicit none
      type(cell), pointer :: cl, cr
      type(face) :: cf
      !
      !  local variables
      real(kind = rfp) :: rhol, rhor
      real(kind = rfp) :: vecn(ndim), tl, tr, a, sl, sr, sm
      type(vector) :: qvl, qvr, qva, evl, evr, qcl, qcr, flux
      real(rfp) :: sqrhol, sqrhor, c1, c2, h0l, h0r, cmax, lamda_l, lamda_r, delta
      real(rfp) :: fac, fac1, h0av, rhoav, rhoppav, alpha
      real(rfp) :: velav(ndim), unav, dl, unl, unr, pl, pr, rhos, ps
      type(matrix) :: aa, zm, zminv, gj, ajl, ajr
      integer :: ic, i
      !
      vecn = cf % vecn
      a = cf % area
      evl = fluxv(vecn, qvl, cl, cf % centp)
      evr = fluxv(vecn, qvr, cr, cf % centp)
      qcl = conservative_qv(qvl)
      qcr = conservative_qv(qvr)
      !   faceflux = evl * tl + evr * tr  - aav * ( qvr - qvl)
      !
      if (neqf > 0) then
         ! special treat for inlet boundary condition
         !
         !  if(bc(cf%itype)%igrp == inlet_type) then
         !    tl = zero; tr = a; aa%e = zero;
         !    flux = evr
         !    return
         !  end if

         rhol = rhof(qvl)
         rhor = rhof(qvr)
         sqrhol = sqrt(rhol)
         sqrhor = sqrt(rhor)
         rhoav = sqrhol * sqrhor
         fac = sqrhol / (sqrhol + sqrhor)
         fac1 = one - fac
         qva % v(:neqf) = qvl % v(:neqf) * fac + qvr % v(:neqf) * fac1
         !   h0l    = h0f(qvl)
         !   h0r    = h0f(qvr)
         !   h0av   = (fac * h0l /rhol + fac1 * h0r / rhor ) * rhoav
         !   call qvfromh0(rhoav,h0av,qva,qvl,qvr)
         unav = dot_product(vecn, qva % v(imb:ime) - mesh_vel(cf % centp))
         call precondition_face(s_ipre, qva, cf, cl, cr, rhoppav)
         call eigenvalues(s_ipre, c1, c2, vecn, qva, rhoppav, cf % centp)
         cmax = max(abs(c1), abs(c2))
         gj = gamma(s_ipre, qva, rhoppav)
      end if
      !
      if (s_cfl > zero) then
         cl % srf = cl % srf + cmax * a
         cr % srf = cr % srf + cmax * a
      else
         cl % srf = max(cl % srf, cmax)
         cr % srf = max(cr % srf, cmax)
      end if
      ajl = ajacobian(vecn, qvl, cl, cf % centp)
      ajr = ajacobian(vecn, qvr, cr, cf % centp)
      aa % e = zero
      if (neqf > 0) then
         !      unl = sum(vecn * (qvl%v(imb:ime) - mesh_vel(cf%centp)))
         !      unr = sum(vecn * (qvr%v(imb:ime) - mesh_vel(cf%centp)))
         call eigenvalues(s_ipre, tl, sl, vecn, qvl, rhoppav, cf % centp)
         call eigenvalues(s_ipre, sr, tr, vecn, qvr, rhoppav, cf % centp)
         sl = min(sl, c2)
         sr = max(sr, c1)
         tr = (min(sr, zero) - min(zero, sl)) / (sr - sl) * a
         tl = (a - tr)
         cmax = (sr * abs(sl) - sl * abs(sr)) * half / (sr - sl) * a
         aa % e(:neqf,:neqf) = gj % e(:neqf,:neqf) * cmax
         flux % v(:neqf) = evl % v(:neqf) * tl + evr % v(:neqf) * tr - &
            matmul(aa % e(:neqf,:neqf), qvr % v(:neqf) - qvl % v(:neqf))
         ajl % e(:neqf,:neqf) = ajl % e(:neqf,:neqf) * tl + aa % e(:neqf,:neqf)
         ajr % e(:neqf,:neqf) = ajr % e(:neqf,:neqf) * tr - aa % e(:neqf,:neqf)

         !     aa%e(:ien,:neqf) = gj%e(:ien,:neqf) * cmax
         !     flux%v(:ien) = evl%v(:ien) * tl + evr%v(:ien) * tr  - &
         !           matmul(aa%e(:ien,:neqf), qvr%v(:neqf) - qvl%v(:neqf)) 
         !     ajl%e(:ien,:neqf) = ajl%e(:ien,:neqf) * tl + aa%e(:ien,:neqf)
         !     ajr%e(:ien,:neqf) = ajr%e(:ien,:neqf) * tr - aa%e(:ien,:neqf)
         !
         !
         !    if(neqf > ien) then 
         !      unl = sum(vecn * (qvl%v(imb:ime) - mesh_vel(cf%centp)))
         !      unr = sum(vecn * (qvr%v(imb:ime) - mesh_vel(cf%centp)))
         !      if(unav > 0) then   ! pure upwind
         !      qva = rhoxf(qvl)
         !      do i = ien + 1, neqf
         !        flux%v(i) = rhol * unl * qvl%v(i) * a
         !       ajl%e(i,:neqf) = qva%v(:neqf) * unl * qvl%v(i) * a
         !       ajl%e(i,imb:ime) = ajl%e(i,imb:ime) + rhol * vecn * qvl%v(i) * a
         !       ajl%e(i,i) = ajl%e(i,i) + rhol * unl * a
         !       ajr%e(i,:neqf) = zero
         !      end do
         !      else
         !      qva = rhoxf(qvr)
         !      do i = ien + 1, neqf
         !        flux%v(i) = rhor * unr * qvr%v(i) * a
         !       ajr%e(i,:neqf) = qva%v(:neqf) * unr * qvr%v(i) * a
         !       ajr%e(i,imb:ime) = ajr%e(i,imb:ime) + rhor * vecn * qvr%v(i) * a
         !       ajr%e(i,i) = ajr%e(i,i) + rhor * unr * a
         !       ajl%e(i,:neqf) = zero
         !      end do
         !      end if
         !     end if
         !---

         !     if(neqf > ien) then
         !      tl = half * a; tr = half * a; cmax = abs(unav * a * half)
         !      do i = ien+1,neqf
         !       aa%e(i,:neqf) = gj%e(i,:neqf) * cmax
         !       flux%v(i) = evl%v(i) * tl + evr%v(i) * tr  - &
         !           sum(aa%e(i,:neqf)*(qvr%v(:neqf) - qvl%v(:neqf)))
         !       ajl%e(i,:neqf) = ajl%e(i,:neqf) * tl + aa%e(i,:neqf)
         !       ajr%e(i,:neqf) = ajr%e(i,:neqf) * tr - aa%e(i,:neqf)
         !      end do
         !     end if
      end if
      !
      if (neqm > 0) then ! Maxwell
         call aaws(tl, sl, cl)
         call aaws(tr, sr, cr)
         !     tl = half * ( tl + tr)  !min(tl,sl)
         !     tr = half * ( sl + sr)  !min(tl,sl)
         tl = tl * tr / (tl + tr) !min(tl,sl)
         tr = sl * sr / (sl + sr) !min(tl,sl)
         sl = tl * tr
         if (s_cfl > zero) then
            cl % srfm = cl % srfm + sl * a
            cr % srfm = cr % srfm + sl * a
         else
            cl % srfm = max(cl % srfm, sl)
            cr % srfm = max(cr % srfm, sl)
         end if
         call meigenmatrix(tl, tr, cf % vecn, aa % e)
         aa % e(ibb:ife, ibb:ife) = aa % e(ibb:ife, ibb:ife) * (a * half)
         tl = a * half; tr = a * half
         !     flux = evl * tl + evr * tr  - aa * ( qvr - qvl)
         flux % v(ibb:ife) = evl % v(ibb:ife) * tl + evr % v(ibb:ife) * tr - &
            matmul(aa % e(ibb:ife, ibb:ife), qvr % v(ibb:ife) - qvl % v(ibb:ife))
         ajl % e(ibb:ife, ibb:ife) = ajl % e(ibb:ife, ibb:ife) * tl + aa % e(ibb:ife, ibb:ife)
         ajr % e(ibb:ife, ibb:ife) = ajr % e(ibb:ife, ibb:ife) * tr - aa % e(ibb:ife, ibb:ife)
      end if
      !     ajl = ajl + aa
      !     ajr = ajr - aa
      !
   end subroutine hlle_flux

   subroutine facefluxscheme(cf, cl, cr, qvl, qvr, flux, aa, tl, tr, ajl, ajr)
      implicit none
      type(cell), pointer :: cl, cr
      type(face) :: cf
      !
      !  local variables
      real(kind = rfp) :: rhol, rhor
      real(kind = rfp) :: vecn(ndim), tl, tr, a, sl, sr, sm
      type(vector) :: qvl, qvr, qva, evl, evr, qcl, qcr, flux
      real(rfp) :: sqrhol, sqrhor, c1, c2, h0l, h0r, cmax, lamda_l, lamda_r, delta
      real(rfp) :: fac, fac1, h0av, rhoav, rhoppav, alpha
      real(rfp) :: velav(ndim), unav, dl, unl, unr, pl, pr, rhos, ps
      type(matrix) :: aa, zm, zminv, gj, ajl, ajr

      integer :: ic, i
      type(matrix) :: aa2, zm2, zminv2
      real(rfp) :: alpha_blend

      alpha_blend = 0.5
      !
      vecn = cf % vecn
      a = cf % area

      !convective flux and conservative variables
      evl = fluxv(vecn, qvl, cl, cf % centp)
      evr = fluxv(vecn, qvr, cr, cf % centp)
      qcl = conservative_qv(qvl)
      qcr = conservative_qv(qvr)

      !   faceflux = evl * tl + evr * tr  - aav * ( qvr - qvl)
      !
      if (neqf > 0) then
         ! special treat for inlet boundary condition
         !
         !  if(bc(cf%itype)%igrp == inlet_type) then
         !    tl = zero; tr = a; aa%e = zero;
         !    flux = evr
         !    return
         !  end if

         ! calculate Roe's average
         rhol = rhof(qvl)
         rhor = rhof(qvr)
         h0l = h0f(qvl)
         h0r = h0f(qvr)

         sqrhol = sqrt(rhol)
         sqrhor = sqrt(rhor)
         rhoav = sqrhol * sqrhor
         fac = sqrhol / (sqrhol + sqrhor)
         fac1 = one - fac
         !   qva%v(ico) = fac * qvl%v(ico) + fac1 * qvr%v(ico)
         !   fac = (rhol / rhoav) * fac
         !   fac1 = (rhor / rhoav) * fac1
         qva % v(:neqf) = qvl % v(:neqf) * fac + qvr % v(:neqf) * fac1
         !   qva%v(imb:ime) = (qvl%v(imb:ime) * rhol * fac + &
         !                     qvr%v(imb:ime) * rhor * fac1 ) / rhoav
         if(neeq /= 0) then
            if (qva % v(ien) <= zero) then
               print *, qvl, qvr
               print *, cf % itype, cf % area
               print *, rhol, rhor, fac
               do i = 1, size(cl % sface)
                  print *, cl % sface(i) % to_face % itype
               end do
            end if
         end if
         !
         ! correct primary variables of Roe's average from total enthalpy 
         h0av = h0l * fac + h0r * (1-fac)
         !h0av = h0l*(1-fac) + h0r*fac

         call qvfromh0(rhoav, h0av, qva, qvl, qvr)
         !
         if (s_ischeme == -1) qva = (qvl + qvr) * (half)
         unav = dot_product(vecn, qva % v(imb:ime) - mesh_vel(cf % centp))

         ! calculate the preconditioning factor rhopp according to Roe's average
         call precondition_face(s_ipre, qva, cf, cl, cr, rhoppav)

         ! calculate the eigenvalues of the preconditioned system
         ! calculate the preconditioned primary jacobian
         call eigenvalues(s_ipre, c1, c2, vecn, qva, rhoppav, cf % centp)
         gj = gamma(s_ipre, qva, rhoppav)

         cmax = max(abs(c1), abs(c2))
      end if

      !
      !  if(neqm > 0.and.neqf == 0) cmax = s_cspeed
      !
      ! spectrum radius
      if (s_cfl > zero) then
         cl % srf = cl % srf + cmax * a
         cr % srf = cr % srf + cmax * a
      else
         cl % srf = max(cl % srf, cmax)
         cr % srf = max(cr % srf, cmax)
      end if
      !
      select case(s_ischeme)
      case(0, -1) ! Roe's
         ! entropy fix
         ! c1=c+, c2=c-
         !
         !if(unav > zero) then   ! flow from left to right
         !   if(abs(c2) < 0.125_rfp) then ! check entropy violation
         !      call eigenvalues(0,delta,lamda_l,vecn,qvl,rhoppav,cf%centp)
         !      call eigenvalues(0,delta,lamda_r,vecn,qvr,rhoppav,cf%centp)
         !      delta = max(mytiny, lamda_r - lamda_l)
         !      if(abs(c2) < delta) then
         !         c2 = (c2*c2 / delta + delta) * half
         !      end if
         !   end if
         !else
         !   if(abs(c1) < 0.125_rfp) then ! check entropy violation
         !      call eigenvalues(0,delta,lamda_l,vecn,qvl,rhoppav,cf%centp)
         !      call eigenvalues(0,delta,lamda_r,vecn,qvr,rhoppav,cf%centp)
         !      delta = max(mytiny, lamda_l - lamda_r)
         !      if(abs(c1) < delta) then
         !         c1 = (c1*c1 / delta + delta) * half
         !      end if
         !   end if
         !end if    
         !!
         aa % e = zero
         if (neqf > 0) then
            ! calculate the eigen vectors of the preconditioned system
            call eigenvector(c1, c2, zm, zminv, cf % vecn, qva, rhoppav, cf % centp)

            zminv % e(:ndim,:neqf) = zminv % e(:ndim,:neqf) * abs(unav)
            zminv % e(ime,:neqf) = zminv % e(ime,:neqf) * abs(c1)
            zminv % e(ien,:neqf) = zminv % e(ien,:neqf) * abs(c2)
            if (neqf > ien) zminv % e(ien + 1:neqf,:neqf) = zminv % e(ien + 1:neqf,:neqf) * abs(unav)
            !     if(s_ivis == 2) then ! special treat for k-omega
            !       zminv%e(ike:iom ,:neqf) = zminv%e(ike:iom,:neqf) * sqrt(sum(qva%v(imb:ime)**2))
            !       zminv%e(ike:iom ,:neqf) = zminv%e(ike:iom,:neqf) * cmax
            !       if(neqf > iom) zminv%e(iom+1:neqf ,:neqf) = zminv%e(iom+1:neqf,:neqf) * abs(unav)
            !     else if(neqf > ien) then
            !       zminv%e(ien+1:neqf ,:neqf) = zminv%e(ien+1:neqf,:neqf) * abs(unav)
            !     end if
            zm % e(:neqf,:neqf) = matmul(zm % e(:neqf,:neqf), zminv % e(:neqf,:neqf))

            ! aa is Roe's average matrix
            aa % e(:neqf,:neqf) = matmul(gj % e(:neqf,:neqf), zm % e(:neqf,:neqf))
         end if
         !
         if (neqm > 0) then ! Maxwell
            call aaws(tl, sl, cl)
            call aaws(tr, sr, cr)
            !     tl = half * ( tl + tr)  !min(tl,sl)
            !     tr = half * ( sl + sr)  !min(tl,sl)
            tl = tl * tr / (tl + tr) !min(tl,sl)
            tr = sl * sr / (sl + sr) !min(tl,sl)
            sl = tl * tr
            if (s_cfl > zero) then
               cl % srfm = cl % srfm + sl * a
               cr % srfm = cr % srfm + sl * a
            else
               cl % srfm = max(cl % srfm, sl)
               cr % srfm = max(cr % srfm, sl)
            end if
            call meigenmatrix(tl, tr, cf % vecn, aa % e)
            !     qvr%v(ieb:iee) = qvr%v(ieb:iee) / e_resistivity(qvr,cr)
            !     qvl%v(ieb:iee) = qvl%v(ieb:iee) / e_resistivity(qvl,cl)
         end if
         ! multiply the area of face
         aa = aa * (a * half)
         !     sl = abs(sum((cl%centp - cf%centp)*cf%vecn))
         !     sr = abs(sum((cr%centp - cf%centp)*cf%vecn))
         !     tl = a * sr /(sl+sr); tr = a * sl / ( sl + sr)
         tl = a * half; tr = a * half
         !    if(s_irealgas == 10) aa%e = zero
         ! flux is the convective flux
         flux = evl * tl + evr * tr - aa * (qvr - qvl)
         !flux = evl * tl + evr * tr
      case(1) ! HLLE's
         aa % e = zero

         if (neqf > 0) then
            call eigenvalues(s_ipre, tl, sl, vecn, qvl, rhoppav, cf % centp)
            call eigenvalues(s_ipre, sr, tr, vecn, qvr, rhoppav, cf % centp)
            sl = min(sl, c2)
            sr = max(sr, c1)
            tr = (min(sr, zero) - min(zero, sl)) / (sr - sl) * a
            tl = (a - tr)
            cmax = (sr * abs(sl) - sl * abs(sr)) * half / (sr - sl) * a
            aa % e(:neqf,:neqf) = gj % e(:neqf,:neqf) * cmax
            if (neqf > ien) aa % e(ien + 1:neqf,:neqf) = gj % e(ien + 1:neqf,:neqf) * abs(unav)
         end if
         if (neqm > 0) then ! Maxwell
            call aaws(tl, sl, cl)
            call aaws(tr, sr, cr)
            !     tl = half * ( tl + tr)  !min(tl,sl)
            !     tr = half * ( sl + sr)  !min(tl,sl)
            tl = tl * tr / (tl + tr) !min(tl,sl)
            tr = sl * sr / (sl + sr) !min(tl,sl)
            sl = tl * tr
            if (s_cfl > zero) then
               cl % srfm = cl % srfm + sl * a
               cr % srfm = cr % srfm + sl * a
            else
               cl % srfm = max(cl % srfm, sl)
               cr % srfm = max(cr % srfm, sl)
            end if
            call meigenmatrix(tl, tr, cf % vecn * (a * half), aa % e)
            !
            !     tl = half * (e_resistivity(qvl,cl)+e_resistivity(qvr,cr))
            !     aa%e(ife,ife) = aa%e(ife,ife) * tl
         end if
         !
         flux = evl * tl + evr * tr - aa * (qvr - qvl)
      case(2) ! HLLC's
         call eigenvalues(s_ipre, tl, sl, vecn, qvl, rhoppav, cf % centp)
         call eigenvalues(s_ipre, sr, tr, vecn, qvr, rhoppav, cf % centp)
         sl = min(sl, c2)
         sr = max(sr, c1)
         if (sl >= zero) then
            tr = zero; tl = a; aa % e = zero; flux = evl * a
            !      ajl = ajacobian(vecn*a,qvl);ajr%e = zero
            return
         end if
         if (sr <= zero) then
            tr = a; tl = zero; aa % e = zero; flux = evr * a
            !      ajr = ajacobian(vecn*a,qvr);ajl%e = zero
            return
         end if

         unl = sum(vecn * qvl % v(imb:ime))
         unr = sum(vecn * qvr % v(imb:ime))
         pl = pf(qvl)
         pr = pf(qvr)
         sm = rhor * unr * (sr - unr) - rhol * unl * (sl - unl) + pl - pr
         sm = sm / (rhor * (sr - unr) - rhol * (sl - unl))
         if (sm > zero) then
            !  F*=(Sm(S_L*Q_L-F_L) +S_L*Q_p*)/(S_L-Sm)
            !  Q_p = (0,pn,p*sm,0)^T
            !
            ps = pl + rhol * (unl - sl) * (unl - sm)
            qva % v = zero
            qva % v(imb:ime) = vecn !ps * vecn
            qva % v(ien) = sm !ps * sm
            lamda_l = sm * a / (sl - sm)
            lamda_r = sl * a / (sl - sm)
            qva = qva * lamda_r
            !      flux = ((qcl * sl - evl)* sm + qva * (sl*ps)) * ( a/ (sl-sm))
            flux = (qcl * sl - evl) * lamda_l + qva * ps
            !      ajl = (gamma(0,qvl,tl) * sl - ajacobian(vecn,qvl))* lamda_l
            !      ajl%e(:,ico) = ajl%e(:,ico) + qva%v;ajr%e = zero
         else
            ps = pr + rhor * (unr - sr) * (unr - sm)
            qva % v = zero
            qva % v(imb:ime) = vecn !ps * vecn
            qva % v(ien) = sm !ps * sm
            lamda_l = sm * a / (sr - sm)
            lamda_r = sr * a / (sr - sm)
            qva = qva * lamda_r
            !      flux = ((qcr * sr - evr)* sm + qva * (sr*ps)) * ( a/ (sr-sm))
            flux = (qcr * sr - evr) * lamda_l + qva * ps
            !      ajr = (gamma(0,qvr,tl) * sr - ajacobian(vecn,qvr))* lamda_l
            !      ajr%e(:,ico) = ajr%e(:,ico) + qva%v;ajl%e = zero
         end if
         tr = (min(sr, zero) - min(zero, sl)) / (sr - sl) * a
         tl = (a - tr)
         cmax = (sr * abs(sl) - sl * abs(sr)) * half / (sr - sl) * a
         aa = gj * cmax
      case(3) ! Linde's HLLC
         call eigenvalues(s_ipre, tl, sl, vecn, qvl, rhoppav, cf % centp)
         call eigenvalues(s_ipre, sr, tr, vecn, qvr, rhoppav, cf % centp)

         sl = min(sl, c2)
         sr = max(sr, c1)
         !
         if (sl > zero) then
            tl = a; tr = zero; aa % e = zero
            flux = evl * a
            return
         end if

         if (sr < zero) then
            tl = zero; tr = a; aa % e = zero
            flux = evr * a
            return
         end if
         tl = sr / (sr - sl) * a
         tr = -sl / (sr - sl) * a
         unl = sum(vecn * qvl % v(imb:ime))
         unr = sum(vecn * qvr % v(imb:ime))
         pl = pf(qvl)
         pr = pf(qvr)
         sm = rhor * unr * (sr - unr) - rhol * unl * (sl - unl) + pl - pr
         sm = sm / (rhor * (sr - unr) - rhol * (sl - unl))
         qva = evr - evl - (qcr - qcl) * sm
         alpha = sum(abs(qva % v)) / sum(abs(qcr % v - qcl % v)) / (c1 - c2)
         alpha = max(zero, one - alpha)
         if (sm > zero) then
            cmax = sl * ((one - alpha) * sr + alpha * sm) / (sr - sl) * a
         else
            cmax = sr * ((one - alpha) * sl + alpha * sm) / (sr - sl) * a
         end if
         aa = gj * (-cmax)
         flux = evl * tl + evr * tr - aa * (qvr - qvl) ! * cmax
      case(4) ! Rusanov's (TVD Lax_Friedrich)
         aa % e = zero
         if (neqf > 0) then
            aa % e(:neqf,:neqf) = gj % e(:neqf,:neqf) * cmax
         end if
         !
         if (neqm > 0) then ! Maxwell
            call aaws(tl, sl, cl)
            call aaws(tr, sr, cr)
            !     tl = half * ( tl + tr)  !min(tl,sl)
            !     tr = half * ( sl + sr)  !min(tl,sl)
            tl = tl * tr / (tl + tr) !min(tl,sl)
            tr = sl * sr / (sl + sr) !min(tl,sl)
            sl = tl * tr
            if (s_cfl > zero) then
               cl % srfm = cl % srfm + sl * a
               cr % srfm = cr % srfm + sl * a
            else
               cl % srfm = max(cl % srfm, sl)
               cr % srfm = max(cr % srfm, sl)
            end if
            call meigenmatrix(tl, tr, cf % vecn, aa % e)
         end if
         aa = aa * (a * half)
         tl = a * half; tr = a * half
         !    if(s_irealgas == 10) aa%e = zero
         flux = evl * tl + evr * tr - aa * (qvr - qvl)
      case(5) ! AUSM+
         !call eigenvalues(s_ipre,tl,sl,vecn,qvl,rhoppav,cf%centp)
         !call eigenvalues(s_ipre,sr,tr,vecn,qvr,rhoppav,cf%centp)
         !sl = min(sl,c2)
         !sr = max(sr,c1)
         !tr = (min(sr,zero) - min(zero,sl)) / (sr - sl) * a
         !tl = (a - tr )
         !cmax = (sr * abs(sl) - sl * abs(sr) ) * half / (sr - sl) * a

         sm = half * (c1 - c2) !sqrt(sound_speed(qva))   ! average sound speed
         unl = sum(vecn * qvl % v(imb:ime))
         unr = sum(vecn * qvr % v(imb:ime))
         sl = unl / sm ! M_L  left mach number
         sr = unr / sm ! M_R  right mach number
         c1 = m4fun(sl, tl, 1); c2 = m4fun(sr, tr, -1) ! M+ M-
         qcl % v(ien) = qcl % v(ien) + pf(qvl) ! rho_L*Psi_L
         qcr % v(ien) = qcr % v(ien) + pf(qvr) ! rho_R*Psi_R
         cmax = (c1 + c2) * a
         if (cmax > zero) then
            flux = qcl * (sm * cmax)
            ajl = gamma(0, qvl, rhoppav)
            ajl % e(ien, ico) = ajl % e(ien, ico) + one
            ajl = ajl * (sm * cmax)
            do i = imb, ime
               ajl % e(:, i) = ajl % e(:, i) + qcl % v * (a * tl * vecn(i - imb + 1))
               ajr % e(:, i) = ajr % e(:, i) + qcl % v * (a * tr * vecn(i - imb + 1))
            end do
         else
            flux = qcr * (sm * cmax)
            ajr = gamma(0, qvr, rhoppav)
            ajr % e(ien, ico) = ajr % e(ien, ico) + one
            ajr = ajr * (sm * cmax)
            do i = imb, ime
               ajl % e(:, i) = ajl % e(:, i) + qcr % v * (a * tl * vecn(i - imb + 1))
               ajr % e(:, i) = ajr % e(:, i) + qcr % v * (a * tr * vecn(i - imb + 1))
            end do
         end if
         lamda_l = p5fun(sl, tl, 1); lamda_r = p5fun(sr, tr, -1) !p+ p-
         flux % v(imb:ime) = flux % v(imb:ime) + &
            a * (lamda_l * qvl % v(ico) + lamda_r * qvr % v(ico)) * vecn
         ajl % e(imb:ime, ico) = ajl % e(imb:ime, ico) + lamda_l * a * vecn
         ajr % e(imb:ime, ico) = ajr % e(imb:ime, ico) + lamda_r * a * vecn
         do i = imb, ime
            ajl % e(imb:ime, i) = ajl % e(imb:ime, i) + vecn * qvl % v(ico) * (a * tl * vecn(i - imb + 1)) / sm
            ajr % e(imb:ime, i) = ajr % e(imb:ime, i) + vecn * qvr % v(ico) * (a * tr * vecn(i - imb + 1)) / sm
         end do

      case(11) ! Blended Scalar Scheme to fix excessive artifial dissipation
         aa % e = zero; aa2 % e = zero
         if (neqf > 0) then
            ! calculate the eigen vectors of the preconditioned system
            call eigenvector(c1, c2, zm, zminv, cf % vecn, qva, rhoppav, cf % centp)
            zminv2 % e = zminv % e; zm2 % e= zm % e 

            ! conventional artificial dissipation
            zminv % e(:ndim,:neqf) = zminv % e(:ndim,:neqf) * abs(unav)
            zminv % e(ime,:neqf) = zminv % e(ime,:neqf) * abs(c1)
            zminv % e(ien,:neqf) = zminv % e(ien,:neqf) * abs(c2)
            if (neqf > ien) zminv % e(ien + 1:neqf,:neqf) = zminv % e(ien + 1:neqf,:neqf) * abs(unav)
            zm % e(:neqf,:neqf) = matmul(zm % e(:neqf,:neqf), zminv % e(:neqf,:neqf))
            aa % e(:neqf,:neqf) = matmul(gj % e(:neqf,:neqf), zm % e(:neqf,:neqf))

            ! scalar modification of artificial dissipation
            call put_diag(aa2 % e, max(abs(c1),abs(c2)), ico, ico)
            call put_diag(aa2 % e, abs(unav), imb, neqf)
            aa2 % e(:neqf, :neqf) = matmul(gj % e(:neqf, :neqf), aa2 % e(:neqf, :neqf))
         end if
         !
         ! multiply the area of face
         aa = aa * (a * half)
         aa2 = aa2 * (a * half * alpha_blend)

         ! modify the convective flux
         tl = a * half; tr = a * half

         flux = evl * tl + evr * tr - aa2 * (qvr - qvl)
      end select

   end subroutine facefluxscheme

   subroutine hllc_flux(cf, cl, cr, qvl, qvr, flux, aa, tl, tr, ajl, ajr)
      implicit none
      type(cell), pointer :: cl, cr
      type(face) :: cf
      !
      !  local variables
      real(kind = rfp) :: rhol, rhor
      real(kind = rfp) :: vecn(ndim), tl, tr, a, sl, sr, sm
      type(vector) :: qvl, qvr, qva, evl, evr, qcl, qcr, qsl, qsr, flux, dsmdq, dwdq
      real(rfp) :: sqrhol, sqrhor, c1, c2, h0l, h0r, cmax, lamda_l, lamda_r, delta
      real(rfp) :: fac, fac1, h0av, rhoav, rhoppav, alpha
      real(rfp) :: velav(ndim), unav, dl, unl, unr, pl, pr, rhos, ps, ug
      type(matrix) :: aa, zm, zminv, gj, ajl, ajr
      integer :: ic, i
      !
      vecn = cf % vecn
      a = cf % area
      evl = fluxv(vecn, qvl, cl, cf % centp)
      evr = fluxv(vecn, qvr, cr, cf % centp)
      qcl = conservative_qv(qvl)
      qcr = conservative_qv(qvr)
      if (neqf > 0) then
         rhol = rhof(qvl)
         rhor = rhof(qvr)
         sqrhol = sqrt(rhol)
         sqrhor = sqrt(rhor)
         rhoav = sqrhol * sqrhor
         fac = sqrhol / (sqrhol + sqrhor)
         fac1 = one - fac
         qva % v(:neqf) = qvl % v(:neqf) * fac + qvr % v(:neqf) * fac1
         !   h0l    = h0f(qvl)
         !   h0r    = h0f(qvr)
         !   h0av   = (fac * h0l /rhol + fac1 * h0r / rhor ) * rhoav
         !   call qvfromh0(rhoav,h0av,qva,qvl,qvr)
         ug = dot_product(vecn, mesh_vel(cf % centp))
         unav = dot_product(vecn, qva % v(imb:ime)) - ug
         call precondition_face(s_ipre, qva, cf, cl, cr, rhoppav)
         call eigenvalues(s_ipre, c1, c2, vecn, qva, rhoppav, cf % centp)
         cmax = max(abs(c1), abs(c2))
         gj = gamma(s_ipre, qva, rhoppav)
         if (s_cfl > zero) then
            cl % srf = cl % srf + cmax * a
            cr % srf = cr % srf + cmax * a
         else
            cl % srf = max(cl % srf, cmax)
            cr % srf = max(cr % srf, cmax)
         end if
         !
         call eigenvalues(s_ipre, tl, sl, vecn, qvl, rhoppav, cf % centp)
         call eigenvalues(s_ipre, sr, tr, vecn, qvr, rhoppav, cf % centp)
         sl = min(sl, c2)
         sr = max(sr, c1)

         if (sl >= zero) then
            flux = evl * a
            ajl = ajacobian(vecn * a, qvl, cl, cl % centp); ajr % e = zero
         else if (sr <= zero) then
            flux = evr * a
            ajr = ajacobian(vecn * a, qvr, cr, cr % centp); ajl % e = zero
         else

            unl = sum(vecn * (qvl % v(imb:ime) - mesh_vel(cf % centp)))
            unr = sum(vecn * (qvr % v(imb:ime) - mesh_vel(cf % centp)))
            pl = pf(qvl) !qvl%v(ico)
            pr = pf(qvr) !qvr%v(ico)
            sm = rhor * unr * (sr - unr) - rhol * unl * (sl - unl) + pl - pr
            sm = sm / (rhor * (sr - unr) - rhol * (sl - unl))
            if (sm > zero) then
               !  F*=(Sm(S_L*Q_L-F_L) +S_L*Q_p*)/(S_L-Sm)
               !  Q_p = (0,pn,p*sm,0)^T
               !
               ps = pl + rhol * (unl - sl) * (unl - sm)
               qva % v = zero
               qva % v(imb:ime) = vecn !ps * vecn
               qva % v(ien) = sm !ps * sm
               qsl = qcl * ((sl - unl)/(sl - sm)) ! + qva * (( ps - pl) / (sl - sm))
               qsl % v(imb:ime) = qsl % v(imb:ime) + vecn * ((ps - pl) / (sl - sm))
               qsl % v(ien) = qsl % v(ien) + ((sm * ps - unl * pl) / (sl - sm))
               flux = qsl * sm + qva * ps
               flux % v(ien) = flux % v(ien) + ug * ps

               ! Left
               rhoav = rhor * (sr - unr) - rhol * (sl - unl)
               dsmdq = rhoxf(qvl) * ((sl - unl)*(sm - unl))
               dsmdq % v(imb:ime) = dsmdq % v(imb:ime) + vecn * (rhol * (two * unl - sm - sl))
               dsmdq % v(ico) = dsmdq % v(ico) + one
               dsmdq = dsmdq * (one / rhoav)
               !
               ajl = vec2mat(dsmdq, qsl)
               !  dP*/dQp
               dwdq = qva * (rhor * (sr - unr))
               dwdq % v(ien) = dwdq % v(ien) + ps
               ajl = ajl + vec2mat(dsmdq, dwdq)
               rhoav = sm / (sl - sm)
               ajl = ajl * (one + rhoav) + gamma(0, qvl, zero) * ((sl - unl) * rhoav)
               dwdq % v = zero; dwdq % v(imb:ime) = vecn ! dU/dq
               ajl = ajl - vec2mat(dwdq, qcl * rhoav)
               !
               ajl % e(imb:ime, ico) = ajl % e(imb:ime, ico) - vecn * rhoav
               ajl % e(ien, ico) = ajl % e(ien, ico) - unl * rhoav
               ajl % e(ien, imb:ime) = ajl % e(ien, imb:ime) - vecn * (rhoav * pl)
               ajl % e(ien,:) = ajl % e(ien,:) + ug * rhor * (sr - unr) * dsmdq % v
               ! Right
               rhoav = rhor * (sr - unr) - rhol * (sl - unl)
               dsmdq = rhoxf(qvr) * ((sr - unr)*(sm - unr))
               dsmdq % v(imb:ime) = dsmdq % v(imb:ime) + vecn * (rhor * (two * unr - sm - sr))
               dsmdq % v(ico) = dsmdq % v(ico) + one
               dsmdq = dsmdq * (-one / rhoav)
               ajr = vec2mat(dsmdq, qsl)
               !  dP*/dQp
               dwdq = qva * (rhol * (sl - unl))
               dwdq % v(ien) = dwdq % v(ien) + ps
               ajr = ajr + vec2mat(dsmdq, dwdq)
               rhoav = sl / (sl - sm)
               ajr = ajr * rhoav
               ajr % e(ien,:) = ajr % e(ien,:) + ug * rhol * (sl - unl) * dsmdq % v
            else
               ps = pr + rhor * (unr - sr) * (unr - sm)
               qva % v = zero
               qva % v(imb:ime) = vecn !ps * vecn
               qva % v(ien) = sm !ps * sm
               qsr = qcr * ((sr - unr)/(sr - sm)) ! + qva * (( ps - pr) / (sr - sm))
               qsr % v(imb:ime) = qsr % v(imb:ime) + vecn * ((ps - pr) / (sr - sm))
               qsr % v(ien) = qsr % v(ien) + ((sm * ps - unr * pr) / (sr - sm))
               flux = qsr * sm + qva * ps
               flux % v(ien) = flux % v(ien) + ug * ps
               ! Left
               rhoav = rhor * (sr - unr) - rhol * (sl - unl)
               dsmdq = rhoxf(qvl) * ((sl - unl)*(sm - unl))
               dsmdq % v(imb:ime) = dsmdq % v(imb:ime) + vecn * (rhol * (two * unl - sm - sl))
               dsmdq % v(ico) = dsmdq % v(ico) + one
               dsmdq = dsmdq * (one / rhoav)
               !
               ajl = vec2mat(dsmdq, qsr)
               !  dP*/dQp
               dwdq = qva * (rhor * (sr - unr))
               dwdq % v(ien) = dwdq % v(ien) + ps
               ajl = ajl + vec2mat(dsmdq, dwdq)
               rhoav = sr / (sr - sm)
               ajl = ajl * rhoav
               ajl % e(ien,:) = ajl % e(ien,:) + ug * rhor * (sr - unr) * dsmdq % v
               !
               ! Right
               rhoav = rhor * (sr - unr) - rhol * (sl - unl)
               dsmdq = rhoxf(qvr) * ((sr - unr)*(sm - unr))
               dsmdq % v(imb:ime) = dsmdq % v(imb:ime) + vecn * (rhor * (two * unr - sm - sr))
               dsmdq % v(ico) = dsmdq % v(ico) + one
               dsmdq = dsmdq * (-one / rhoav)
               ajr = vec2mat(dsmdq, qsr)
               !  dP*/dQp
               dwdq = qva * (rhol * (sl - unl))
               dwdq % v(ien) = dwdq % v(ien) + ps
               ajr = ajr + vec2mat(dsmdq, dwdq)
               rhoav = sm / (sr - sm)
               ajr = ajr * (one + rhoav) + gamma(0, qvr, zero) * ((sr - unr) * rhoav)
               dwdq % v = zero; dwdq % v(imb:ime) = vecn ! dU/dq
               ajr = ajr - vec2mat(dwdq, qcr * rhoav)
               !
               ajr % e(imb:ime, ico) = ajr % e(imb:ime, ico) - vecn * rhoav
               ajr % e(ien, ico) = ajr % e(ien, ico) - unr * rhoav
               ajr % e(ien, imb:ime) = ajr % e(ien, imb:ime) - vecn * (rhoav * pr)
               ajr % e(ien,:) = ajr % e(ien,:) + ug * rhol * (sl - unl) * dsmdq % v
               !
            end if
            flux % v(:neqf) = flux % v(:neqf) * a
            ajl = ajl * a
            ajr = ajr * a
         end if
      end if
      !
      if (neqm > 0) then ! Maxwell
         call aaws(tl, sl, cl)
         call aaws(tr, sr, cr)
         !     tl = half * ( tl + tr)  !min(tl,sl)
         !     tr = half * ( sl + sr)  !min(tl,sl)
         tl = tl * tr / (tl + tr) !min(tl,sl)
         tr = sl * sr / (sl + sr) !min(tl,sl)
         sl = tl * tr
         if (s_cfl > zero) then
            cl % srfm = cl % srfm + sl * a
            cr % srfm = cr % srfm + sl * a
         else
            cl % srfm = max(cl % srfm, sl)
            cr % srfm = max(cr % srfm, sl)
         end if
         call meigenmatrix(tl, tr, cf % vecn, aa % e)
         aa % e(ibb:ife, ibb:ife) = aa % e(ibb:ife, ibb:ife) * (a * half)
         tl = a * half; tr = a * half
         !     flux = evl * tl + evr * tr  - aa * ( qvr - qvl)
         flux % v(ibb:ife) = evl % v(ibb:ife) * tl + evr % v(ibb:ife) * tr - &
            matmul(aa % e(ibb:ife, ibb:ife), qvr % v(ibb:ife) - qvl % v(ibb:ife))
         ajl % e(ibb:ife, ibb:ife) = ajl % e(ibb:ife, ibb:ife) * tl + aa % e(ibb:ife, ibb:ife)
         ajr % e(ibb:ife, ibb:ife) = ajr % e(ibb:ife, ibb:ife) * tr - aa % e(ibb:ife, ibb:ife)
      end if
      !
   end subroutine hllc_flux


   function m4fun(m, dm, plus)result(mf)
      real(rfp) :: m, mf, dm
      integer, intent(in) :: plus
      if (plus > 0) then ! plus part
         if (abs(m) >= one) then
            mf = max(m, zero)
            dm = half * (one + sign(1._rfp, m))
         else
            mf = 0.25_rfp * (m + one)**2 + 0.125_rfp * (m * m - one)**2
            dm = half * ((m + one) + m * (m * m - one))
         end if
      else
         if (abs(m) >= one) then
            mf = min(m, zero)
            dm = half * (sign(1._rfp, m) - one)
         else
            mf = -0.25_rfp * (m - one)**2 - 0.125_rfp * (m * m - one)**2
            dm = -half * ((m - one) + m * (m * m - one))
         end if
      end if
   end function m4fun

   function p5fun(m, dm, plus)result(mf)
      real(rfp) :: m, mf, dm
      integer, intent(in) :: plus
      if (plus > 0) then ! plus part
         if (abs(m) >= one) then
            mf = max(m, zero) / (m + mytiny)
            dm = zero
         else
            mf = 0.25 * (m + one)**2 * (two - m) + 0.1875_rfp * m * (m * m - one)**2
            dm = half * (m + one) * (two - m) - 0.25 * (m + one)**2 - 0.1875_rfp * (m * m - one)**2 &
               -0.75_rfp * m * m * (m * m - one)
         end if
      else
         if (abs(m) >= one) then
            mf = min(m, zero) / (m - mytiny)
            dm = zero
         else
            mf = 0.25 * (m - one)**2 * (two + m) - 0.1875_rfp * m * (m * m - one)**2
            dm = half * (m - one) * (two + m) + 0.25 * (m - one)**2 - 0.1875_rfp * (m * m - one)**2 &
               -0.75_rfp * m * m * (m * m - one)
         end if
      end if
   end function p5fun

   function m2fun(m, plus)result(mf)
      real(rfp) :: m, mf
      integer, intent(in) :: plus
      if (plus > 0) then ! plus part
         mf = 0.25 * (m + one) **2
      else
         mf = -0.25 * (m - one)**2
      end if
   end function m2fun

   function a_average_jacobian(cf, cl, cr, qvl, qvr, aa1)result(aa)
      implicit none
      type(cell), pointer :: cl, cr
      type(face) :: cf
      !
      !  local variables
      real(kind = rfp) :: rhol, rhor
      real(kind = rfp) :: vecn(ndim)
      type(vector) :: qvl, qvr, qva
      real(rfp) :: sqrhol, sqrhor, c1, c2, h0l, h0r, cmax, lamda_l, lamda_r, delta
      real(rfp) :: fac, fac1, h0av, rhoav, rhoppav
      real(rfp) :: velav(ndim), unav, dl
      type(matrix) :: aa, zm, zminv, gj, aa1
      integer :: ic
      if (neqf == 0) return
      !
      !   qvl = cl%qv
      !   qvr = cr%qv
      vecn = cf % vecn
      !   
      rhol = rhof(qvl)
      rhor = rhof(qvr)
      sqrhol = sqrt(rhol)
      sqrhor = sqrt(rhor)
      rhoav = sqrhol * sqrhor
      fac = sqrhol / (sqrhol + sqrhor)
      fac1 = one - fac
      qva = qvl * fac + qvr * fac1
      qva % v(ico) = qvl % v(ico) * fac1 + qvr % v(ico) * fac

      !   h0l    = h0f(qvl)
      !   h0r    = h0f(qvr)
      !   h0av   = fac * h0l + fac1 * h0r  
      !   call qvfromh0(rhoav,h0av,qva)
      !
      !
      !  qva = (qvl + qvr) * half

      unav = dot_product(vecn, qva % v(imb:ime) - mesh_vel(cf % centp))
      !
      call precondition_face(s_ipre, qva, cf, cl, cr, rhoppav)
      call eigenvalues(s_ipre, c1, c2, vecn, qva, rhoppav, cf % centp)
      ! entropy fix
      ! c1=c+, c2=c-
      !
      if (unav > zero) then ! flow from left to right
         if (abs(c2) < 1.e-3) then ! check entropy violation
            call eigenvalues(0, c1, lamda_l, vecn, qvl, rhoppav, cf % centp)
            call eigenvalues(0, c1, lamda_r, vecn, qvr, rhoppav, cf % centp)
            delta = max(mytiny, lamda_r - lamda_l)
            if (abs(c2) < delta) then
               c2 = (c2 * c2 / delta + delta) * half
            end if
         end if
      else
         if (abs(c1) < 1.e-3) then ! check entropy violation
            call eigenvalues(0, c2, lamda_l, vecn, qvl, rhoppav, cf % centp)
            call eigenvalues(0, c2, lamda_r, vecn, qvr, rhoppav, cf % centp)
            delta = max(mytiny, lamda_l - lamda_r)
            if (abs(c1) < delta) then
               c1 = (c1 * c1 / delta + delta) * half
            end if
         end if
      end if
      !!
      cmax = max(abs(c1), abs(c2))
      !
      !   if(s_ivis > 0) then
      !   fac = s_cfl/s_vnn / dot_product(vecn,cl%centp-cr%centp)
      !   cmax = max(cmax,mu_suth(qva)/rhof(qva) * abs(fac))
      !    cmax = cmax + abs(g_zmu) / (b_rhoinf &
      !           * abs(dot_product(vecn,cl%centp-cr%centp)))
      !   end if
      !   cmax = dot_product(qva%v(imb:ime),vecn)
      !   dl = dot_product(cr%centp-cl%centp,vecn)
      !   cmax = abs(cmax / dl)


      ! try to use inverse average to increase convergence
      !   cmax = one /cmax
      !
      if (s_cfl > zero) then
         !   cmax = cmax * cf%area 
         cl % srf = cl % srf + cmax * cf % area
         cr % srf = cr % srf + cmax * cf % area
      else
         cl % srf = max(cl % srf, cmax)
         cr % srf = max(cr % srf, cmax)
         !   cmax = cmax * cf%area 
      end if

      gj = gamma(s_ipre, qva, rhoppav)
      if (s_irec == 0) then ! large artificial dissipation
         aa = gj * cmax
      else
         call eigenvector(c1, c2, zm, zminv, cf % vecn, qva, rhoppav, cf % centp)
         zminv % e(:ndim,:) = zminv % e(:ndim,:) * abs(unav)
         zminv % e(ime,:) = zminv % e(ime,:) * abs(c1)
         zminv % e(ien,:) = zminv % e(ien,:) * abs(c2)
         if (neq > ien) zminv % e(ien + 1:,:) = zminv % e(ien + 1:,:) * abs(unav)
         !   print *,'aj=',ajacobian(vecn,qva)
         !   aa = zm * zminv   !  zm = Lamda * M
         aa = gj * zm * zminv !  zm = Lamda * M
         !   if(s_alpha < one) then
         !   aa = aa * s_alpha + gj * (cmax * (one -s_alpha))
         !   end if
         !   aa1 =  gj * cmax 
      end if
      aa1 = aa
      !   aa1%e(:,ien) = zero
      !   print *,'av=',aa1
   end function a_average_jacobian

   !********************************************************
   subroutine lhs(cells, faces, nodes)
      !********************************************************
      implicit none
      !
      type(node) :: nodes(:)
      type(cell), pointer :: cells(:)
      type(face) :: faces(:)
      type(face), pointer :: cf
      ! local variables
      type(cell), pointer :: current_cell, cn
      type(vector), pointer :: qv, res
      type(vector) :: qvc
      real(rfp) :: dtau, dtaum, dtinv, coef, rhopp, dl, gear(10)
      type(matrix) :: gm
      type(matrix), pointer :: dm
      integer :: i, j, ndt1
      real(rfp) :: theat = 0.5_rfp, bspeed, cspeed
      !
      ndt1 = s_idt + 1
      !
      do i = 1, size(cells)
         ! share left cell variables
         current_cell => cells(i)

         qv => current_cell % qv
         !    call precondition_cell(s_ipre,current_cell,rhopp)
         if (neqm > 0) call aaws(bspeed, cspeed, current_cell)

         gm = gamma(s_ipre, qv, cells(i) % rp, current_cell, bspeed, cspeed)
         !
         if (s_idt == 0.and.s_dt_inv < zero) then
            dtau = current_cell % vol * abs(s_dt_inv)
            dtau = max(real(s_damping, rfp)/real(s_nstep, rfp), one) * dtau
            dtaum = dtau
         else if (s_cfl > zero) then
            dtau = current_cell % srf / abs(s_cfl)
            !       dtaum   = current_cell%srfm / abs(s_cfl)
            dtaum = current_cell % srfm / abs(s_cflmm)
         else
            dl = length_for_cfl(current_cell, nodes)
            dtau = current_cell % vol * current_cell % srf / abs(s_cfl) / dl
            !       dtaum = current_cell%vol * current_cell%srfm / abs(s_cfl) / dl 
            dtaum = current_cell % vol * current_cell % srfm / abs(s_cflmm) / dl
         end if
         !       dl = length_for_cfl(current_cell,nodes)
         !       dtaum = current_cell%vol * current_cell%srfm / abs(s_cfl) / dl 

         dm => current_cell % dm
         if (s_imeth == 0) then ! explicit
            dm = gm * dtau !* s_dt_inv
            call matrixinv(dm % e)
            current_cell % dqv = dm * current_cell % res
            cycle
         end if
         !    do j = 1, size(current_cell%scell)
         !     cn => current_cell%scell(j)%to_cell
         !     dm = dm + cn%dm
         !    end do
         !     dm = dm *(one / real(size(current_cell%scell),rfp))

         !     dm  = dm + gm * dtau

         ! add preconditioned primary jacobian times 1/dtau
         if (neqf > 0) dm % e(:neqf,:) = dm % e(:neqf,:) + gm % e(:neqf,:) * dtau
         if (neqm > 0) then
            dm % e(ibb:ife,:) = dm % e(ibb:ife,:) + gm % e(ibb:ife,:) * dtaum
         end if
         !
         !    if(s_idt /= 0) then  ! unsteady case
         !     dtinv = s_dt_inv * current_cell%vol
         !     gm   = gamma(0,qv,rhopp)   ! without precondition
         !
         !     res => current_cell%res
         !
         !     if(s_idt == 1) then ! first order
         !      res = res - (gm * dtinv - dm * theat) * (qv - current_cell%qvn)
         !     else if(s_idt >= 2) then
         !!      res = res - gm *(qv * (dtinv * three_second) - &
         !!                       current_cell%qvn * (dtinv * two) + &
         !!                       current_cell%qvm * (dtinv * half) )
         !      res = res - (conservative_qv(qv) * (dtinv * three_second) - &
         !                       conservative_qv(current_cell%qvn) * (dtinv * two) + &
         !                       conservative_qv(current_cell%qvm) * (dtinv * half) )
         !     end if
         !!
         !     if (s_idt == 1) then
         !      coef = dtinv  ! first order time marching
         !      dm = dm * theat + gm * coef
         !     else if (s_idt == 2) then
         !      coef = dtinv * three_second ! first order time marching
         !      dm = dm + gm * coef
         !     end if
         !!     
         !    end if

         ! unsteady
         if (s_idt > 0) then 
            gear(:ndt1) = current_cell % vol * gear_c

            ! without precondition
            gm = gamma(0, qv, rhopp, cells(i), one, s_cspeed) 

            res => current_cell % res
            !conservative variables
            qvc = conservative_qv(qv)

            res = res - qvc * gear(1)
            do j = 1, s_idt
               qvc = conservative_qv(current_cell % qvn(gear_in(j)))
               res = res - qvc * gear(j + 1)
            end do
            dm = dm + gm * gear(1)
            !  
         end if
         !    call matrixinv(dm%e)
         !  initial DeltQ
         !    current_cell%dqv = dm * current_cell%res
      end do

      !    if(s_idt == 1) then
      !     do i = 1, size(faces)
      !      faces(i)%ajl = faces(i)%ajl * theat
      !      faces(i)%ajr = faces(i)%ajr * theat
      !     end do
      !    end if

   end subroutine lhs

   function length_for_cfl(cc, nodes) result(dl)
      type(node) :: nodes(:)
      type(cell), pointer :: cc
      real(rfp) :: dl, p(size(cc % c2n)), dd(ndim), vel(ndim)
      integer :: i
      !
      !  dl = cc%vol
      !  if(ndim==2.and.s_iaxis > 0) dl = dl / cc%centp(s_iaxis)
      !  dl = cc%vol**(one/real(ndim,rfp))
      !  return

      ! if(nmeq > 0) then
      ! vel = abs(cc%qv%v(imb:ime))
      ! else
      vel = one
      ! end if
      do i = 1, ndim
         p = nodes(cc % c2n) % xyz(i)
         dd(i) = maxval(p) - minval(p)
      end do
      dl = sqrt(sum(vel * vel))
      if (dl < mytiny) then
         dl = dd(1)
      else
         dl = dot_product(dd, vel) / dl !sqrt(sum(vel*vel))
      end if

   end function length_for_cfl
   !
   !********************************************************
   subroutine viscous(faces, nodes, cells)
      !********************************************************
      implicit none
      !
      type(face), pointer :: faces(:), cf, mon_face
      type(cell), pointer :: cells(:), lcell, rcell
      type(node), pointer :: nodes(:)
      ! local variables
      type(vector) :: qv, qvf, gqv(ndim), dqv !(ndim) !,qvl,qvr
      type(matrix) :: rR, rL
      real(rfp) :: vecn(ndim), coef, t, diffc(neq)
      real(rfp) :: vjacob(ndim, 2)
      real(rfp) :: dl, dr, vel(ndim), un, p, rho, volume_inv, epsi, alpha, ls_face_sharp, ls_face_smear
      integer, pointer :: f2n(:)
      integer :: i, j, it, k, k1, k2
      !
      !  initial all cell's variables to zero
      epsi = 1.e-20_rfp
      alpha = two
      !  call set_zero(cells,faces)
      if (neqf == 0) return ! Maxwell only
      if (s_ivis == 0 .and. s_ialg == 0) return
      if (s_ivis == 0 .and. s_irec /= 2) return ! Euler eqn

      !
      do i = 1, size(faces)
         !
         cf => faces(i)
         vecn = cf % vecn * cf % area
         it = cf % itype
         lcell => faces(i) % left_cell
         rcell => faces(i) % right_cell

         dl = one !/ sqrt(sum((lcell%centp - cf%centp)**2))
         if (cf % itype == 0) dr = one !/ sqrt(sum((rcell%centp - cf%centp)**2))
         !
         ! try to use Mitchell's stencil to calculate viscous flux
         f2n => cf % f2n
         !
         call face_gradient(cf, rcell, lcell, nodes(f2n), gqv, vjacob)
         !
         ! save gqv for reconstructure
         !
         if (s_irec == 2) then
            k = ndim
            dqv % v = zero
            do k1 = 1, ndim
               dqv % v = dqv % v + abs(gqv(k1) % v) ! * gqv(k1)%v
            end do
            dqv % v = one / (dqv % v + epsi)
            !
            lcell % dqv = lcell % dqv + dqv * dl
            if (cf % itype == 0) rcell % dqv = rcell % dqv + dqv * dr

            do k1 = 1, ndim
               lcell % gradient(:, k1) = lcell % gradient(:, k1) + gqv(k1) % v * dqv % v * dl
               if (cf % itype == 0) rcell % gradient(:, k1) = rcell % gradient(:, k1) + gqv(k1) % v * dqv % v * dr
               if (abs(s_ialg) <= 1) cycle
               do k2 = k1, ndim
                  k = k + 1
                  lcell % gradient(:, k) = lcell % gradient(:, k) + gqv(k1) % v * vecn(k2) !d^2q/dk1dk2
                  if (cf % itype == 0) &
                     rcell % gradient(:, k) = rcell % gradient(:, k) - gqv(k1) % v * vecn(k2) !d^2q/dk1dk2
               end do
            end do
         end if

         if (s_ivis == 0.and.nmeq > 0) cycle
         coef = one / real(size(f2n), rfp)
         !calculate face average value
         qvf = sum_qv_node(nodes(cf % f2n)) * coef

         !calculate diffusion related coefficient
         diffc = diff_coef(qvf, cf)

         !calculate viscous flux as qv
         qv = vis_t(qvf, diffc, gqv, vecn)

         !calculate viscous jacobian as rL and rR
         rL = vis_tc(qvf, diffc, vjacob(:, 1), vecn)
         rR = vis_tc(qvf, diffc, vjacob(:, 2), vecn)

         if (nmeq > 0) then
            rR % e(ien, imb:ime) = rR % e(ien, imb:ime) + qv % v(imb:ime)
            rL % e(ien, imb:ime) = rL % e(ien, imb:ime) + qv % v(imb:ime)
         end if
         !
         if (s_ivis == 2) then ! k-omega
            rho = rhof(qvf)
            !    dqv%v(imb:ime) = qv%v(imb:ime) / diffc(imb) *  rho / qvf%v(iom)
            !    rR%e(imb:ime,ike) = rR%e(imb:ime,ike) + dqv%v(imb:ime)
            !    rL%e(imb:ime,ike) = rL%e(imb:ime,ike) + dqv%v(imb:ime)
            !!
            !    dqv%v(imb:ime) =  - qv%v(imb:ime) / diffc(imb) * mu_turb(qvf) / qvf%v(iom)
            !    rR%e(imb:ime,iom) = rR%e(imb:ime,iom) + dqv%v(imb:ime)
            !    rL%e(imb:ime,iom) = rL%e(imb:ime,iom) + dqv%v(imb:ime)
            !!
            dqv % v(ike:iom) = qv % v(ike:iom) / diffc(ike:iom) * rho / qvf % v(iom) * t_psi
            rR % e(ike:iom, ike) = rR % e(ike:iom, ike) + dqv % v(ike:iom)
            rL % e(ike:iom, ike) = rL % e(ike:iom, ike) + dqv % v(ike:iom)
            !!
            !    dqv%v(ike:iom) =  - qv%v(ike:iom) / diffc(ike:iom) * &
            !                        mu_turb(qvf) / qvf%v(iom) * t_psis
            !    rR%e(ike:iom,iom) = rR%e(ike:iom,iom) + dqv%v(ike:iom)
            !    rL%e(ike:iom,iom) = rL%e(ike:iom,iom) + dqv%v(ike:iom)
         end if

         !assemble res, dm, ajl, and ajr	      
         lcell % res = lcell % res + qv
         rcell % res = rcell % res - qv
         !
         !
         lcell % dm = lcell % dm - rL
         cf % ajr = cf % ajr - rR
         !
         rcell % dm = rcell % dm + rR
         cf % ajl = cf % ajl + rL

      end do ! end of face loop

      if (s_irec == 2) then
         do i = 1, size(cells)
            do k = 1, neq
               cells(i) % gradient(k,:) = cells(i) % gradient(k,:) / cells(i) % dqv % v(k)
            end do
         end do
      end if

      nullify(cf, lcell, rcell, f2n)

   end subroutine viscous
   !
   subroutine set_zero(cells, faces)
      type(cell), pointer :: cells(:), cc
      type(face), pointer :: faces(:)
      integer :: i
      do i = 1, size(cells)
         cells(i) % srf = zero
         cells(i) % srfm = zero
         cells(i) % dm % e = zero
         cells(i) % res % v = zero
         cells(i) % gradient = zero
         cells(i) % dqv % v = zero
      end do
      do i = 1, size(faces)
         faces(i) % ajl % e = zero
         faces(i) % ajr % e = zero
         if (faces(i) % itype == 0) cycle
         cc => faces(i) % right_cell
         cc % srf = zero
         cc % srfm = zero
         cc % dm % e = zero
         cc % res % v = zero
         cc % gradient = zero
         cc % dqv % v = zero
      end do
   end subroutine set_zero

   function vis_t(qvf, diffc, gqv, vecn) result(qv)
      type(vector) :: qv, qvf, gqv(:)
      real(rfp), intent(in) :: diffc(neq), vecn(:)
      real(rfp) :: vel(ndim), hi(nspe), dyi(nspe), zd(nspe), sigma(ndim, ndim), sq
      integer :: k
      qv % v = zero
      vel = zero
      if (nmeq > 0) then
         sq = zero
         do k = 1, ndim
            sq = sq + gqv(k) % v(imb + k - 1)
         end do
         vel = qvf % v(imb:ime)
         do k = 1, ndim
            sigma(k,:) = diffc(imb) * (gqv(k) % v(imb:ime) + gqv(:) % v(imb + k - 1))
            sigma(k, k) = sigma(k, k) + diffc(ime) * sq
         end do

         vel(:) = matmul(sigma, vel)
         qv % v(imb:ime) = matmul(vecn, sigma)
      end if
      if (neeq > 0) then
         vel = vel + diffc(ien) * gqv(:) % v(ien)
         qv % v(ien) = dot_product(vecn, vel(:ndim))
      end if

      if (s_ivis == 2) then
         qv % v(ike) = diffc(ike) * sum(vecn * gqv % v(ike))
         qv % v(iom) = diffc(iom) * sum(vecn * gqv % v(iom))
      end if

      if (nspe > 1) then ! mass diffusion
         do k = isb, ise
            dyi(k - isb + 1) = sum(vecn * gqv % v(k))
            qv % v(k) = diffc(k) * sum(vecn * gqv % v(k))
         end do
         dyi(nspe) = -sum(dyi(:nspm1))
         hi = hif(qvf)
         zd = (/diffc(isb:ise), diffc(ico)/)
         qv % v(ien) = qv % v(ien) + sum(hi * zd * dyi)
      end if

      if (naux > 0) then
         qv % v(iaub) = diffc(iaub) * (sum(vecn * gqv % v(iaub)))
         qv % v(ien) = qv % v(ien) + qvf % v(iaub) * qv % v(iaub)
      end if
      if (nrad > 0) then
         qv % v(irad) = diffc(irad) * (sum(vecn * gqv % v(irad)))
         qv % v(ien) = qv % v(ien) - qv % v(irad) ! radiation flux added to energy eq.
         !   qv%v(ien) = qv%v(ien) + qv%v(irad)  ! radiation flux added to energy eq.
      end if
   end function vis_t

   function vis_tc(qv, diffc, v, vecn) result(r)
      type(vector) :: qv
      type(matrix) :: r
      real(rfp), intent(in) :: diffc(neq), v(:), vecn(:)
      real(rfp) :: sq, zmu2o3, vel(ndim), zdhn, hi(nspe), zd(nspe)
      integer :: l, k
      sq = sum(vecn * v(:))
      r % e = zero
      !
      if (nmeq > 0) then
         vel = qv % v(imb:ime)
         do k = imb, ime
            l = k - imb + 1
            r % e(k, imb:ime) = diffc(imb) * vecn * v(l) + diffc(ime) * vecn(l) * v
            r % e(k, k) = r % e(k, k) + diffc(imb) * sq
         end do
         r % e(ien, imb:ime) = diffc(imb) * (vecn * (sum(v * vel)) + vel * sq) + &
            diffc(ime) * sum(vecn * vel) * v
      end if
      if (neeq > 0) then
         r % e(ien, ien) = diffc(ien) * sq
      end if

      if (s_ivis == 2) then
         r % e(ike, ike) = diffc(ike) * sq
         r % e(iom, iom) = diffc(iom) * sq
      end if

      if (nspe > 1) then ! mass diffusion
         do k = isb, ise
            r % e(k, k) = diffc(k) * sq
         end do
         hi = hif(qv)
         zd = (/diffc(isb:ise), diffc(ico)/)
         zdhn = zd(nspe) * hi(nspe)
         r % e(ien, isb:ise) = (zd(:nspm1) * hi(:nspm1) - zdhn) * sq
      end if

      if (naux > 0) r % e(iaub, iaub) = diffc(iaub) * sq
      if (nrad > 0) then
         r % e(irad, irad) = diffc(irad) * sq
         r % e(ien, irad) = -r % e(irad, irad)
         !    r%e(ien,irad) = r%e(irad,irad)
      end if

   end function vis_tc

   function diff_coef(qv, cf)result(c)
      type(face) :: cf
      type(vector) :: qv
      real(rfp) :: zmu, zk, zmuk, zmuw, zmut, zkt, zdt
      real(rfp) :: zmu2o3, zd(nspe)
      real(rfp) :: c(neq)
      c = zero

      ! viscosity
      zmu = mumix(qv)

      !thermal conductivity
      zk = lamda_at_face(cf, qv)
      zmut = zero; zkt = zero; zdt = zero

      !mass diffusivity
      if(nspe > 1) then
         zd = mass_diffusivity(qv)
      end if

      ! turbulent model
      if (s_ivis >= 2) then
         if (s_ivis == 2) then !K-Omega
            zmut = mu_turb(qv)
            zmuk = zmu + zmut * t_psi
            zmuw = zmu + zmut * t_psis
         else
            if (cf % itype == 0) then
               zmut = (cf % left_cell % zmut + cf % right_cell % zmut) * half
            else
               zmut = cf % left_cell % zmut
            end if
         end if
         if (nspe > 1) then
            zdt = zmut / g_scht
            !     zd = zd * ( one  + zmut / zmu)    !zdt
            zd = zd + zdt
         end if
         zmu = zmu + zmut
         !    zkt = zk * zmut / zmu * g_prt     !* g_pr / g_prt
         zkt = zmut * htf(qv) / g_prt
         zk = zk + zkt
         c(ike) = zmuk
         c(iom) = zmuw
      end if

      if (nmeq > 0) then
         c(imb) = zmu
         c(ime) = -two_third * zmu
      end if
      if (neeq > 0) c(ien) = zk
      if (nspe > 1) then
         c(isb:ise) = zd(1:nspm1)
         c(ico) = zd(nspe)
      end if
      if (naux > 0) c(iaub) = zmu
      !  if(nrad> 0) c(irad) = -one_third / absorption_coef(cf%left_cell,qv)
      if (nrad > 0) c(irad) = -one_third / absorption_coef(qv)
      ! flux-limited radiation diffusion
      !
      !  if(nrad> 0) then
      !   zk = sqrt(sum(cf%left_cell%gradient(irad,:)**2)) / abs(qv%v(irad))
      !   c(irad) = -one_third / (absorption_coef(qv)+zk)
      !  end if

   end function diff_coef
   !
   function lamda_at_face(cf, qvf)result(zk)
      type(vector) :: qvf
      type(face) :: cf
      type(cell), pointer :: lc, rc
      real(rfp) :: zk1, zk2, zk, r1, r2, lsl, lsr
      lc => cf % left_cell
      rc => cf % right_cell

      !   if(cf%itype == 0.or.cf%itype < partition_face_no) then   ! inside faces
      zk1 = lamda(lc % qv, lc)
      zk2 = lamda(rc % qv, rc)

      !   zk = lamda(qvf,lc)
      !   else
      !    zk1 = lamda(lc%qv,lc)
      !    zk2 = rc%rp
      !   end if
      r1 = abs(dot_product(cf % vecn, lc % centp - cf % centp))
      r2 = abs(dot_product(cf % vecn, rc % centp - cf % centp))
      zk = (r1 + r2) * zk1 * zk2 / (zk1 * r2 + zk2 * r1)
   end function lamda_at_face
end module gems_fv
!
