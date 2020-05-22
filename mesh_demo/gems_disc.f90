module gems_disc
   use gems_data
   use gems_library
   use gems_state
   implicit none

contains
   !
   subroutine face_reconstruction(faces, nodes)
      type(face), pointer :: faces(:)
      type(cell), pointer :: rc, lc
      type(face), pointer :: cf
      type(node), pointer :: nodes(:)
      real(rfp) :: a(ndim, ndim)
      integer :: iloc(1), i, k
      integer, pointer :: f2n(:)
      do i = 1, size(faces)
         cf => faces(i)
         lc => cf % left_cell
         rc => cf % right_cell
         f2n => cf % f2n
         a(1,:) = rc % centp - lc % centp
         if(size(f2n)/=1) then
            select case(size(f2n))
            case(2)
               a(2,:) = nodes(f2n(2)) % xyz - nodes(f2n(1)) % xyz
            case(3)
               a(2,:) = nodes(f2n(2)) % xyz - nodes(f2n(1)) % xyz
               a(3,:) = nodes(f2n(3)) % xyz - nodes(f2n(2)) % xyz
            case default
               a(2,:) = nodes(f2n(3)) % xyz - nodes(f2n(1)) % xyz
               a(3,:) = nodes(f2n(4)) % xyz - nodes(f2n(2)) % xyz
            end select
         end if
         !
         call direct_inv(a)
         allocate(cf % avec(ndim, ndim))
         cf % avec = a
      end do

   end subroutine face_reconstruction

   subroutine face_gradient(cf, rcell, lcell, nodes, gqv, vjacob)
      type(face), pointer :: cf
      type(cell) :: rcell, lcell
      type(node) :: nodes(:)
      real(rfp), optional :: vjacob(ndim, 2)
      type(vector) :: gqv(ndim)
      !  if(s_ifmt == 1) then
      call face_gradient_dli2(cf, rcell, lcell, nodes, gqv, vjacob)
   end subroutine face_gradient

   subroutine face_gradient_dli2(cf, rcell, lcell, nodes, gqv, vjacob)
      type(face), pointer :: cf
      type(cell) :: rcell, lcell
      type(node) :: nodes(:)
      real(rfp), optional :: vjacob(ndim, 2)
      type(vector) :: gqv(ndim), dqv(ndim)
      integer :: k, n, i
      type(vector) :: qvl, qvr, qv_node(size(nodes))
      real(rfp) :: qv_db(ndim)

      qvr = rcell%qv
      qvl = lcell%qv

      do i=1,size(nodes)
         qv_node(i) = nodes(i) % qv
      end do

      dqv(1) = qvr - qvl
      select case(size(nodes))
      case(2)
         dqv(2) = qv_node(2) - qv_node(1)
      case(3)
         dqv(2) =  qv_node(2) - qv_node(1)
         dqv(3) =  qv_node(3) - qv_node(2)
      case default
         dqv(2) =  qv_node(3) - qv_node(1)
         dqv(3) =  qv_node(4) - qv_node(2)
      end select
      !
      do k = 1, neq
         qv_db = dqv(:)%v(k)
         qv_db = matmul(cf%avec,dqv(:)%v(k))
         gqv(:)%v(k) = qv_db
         !gqv(:) % v(k) = matmul(cf % avec, dqv(:) % v(k)) !cramer(a,dqv(:)%v(k))
      end do

      if (present(vjacob)) then
         vjacob(:, 1) = -cf % avec(:, 1) !left cell
         vjacob(:, 2) = cf % avec(:, 1) !right cell
      end if
   end subroutine face_gradient_dli2

   subroutine face_gradient_dli1(cf, rcell, lcell, nodes, gqv, vjacob)
      type(face), pointer :: cf
      type(cell) :: rcell, lcell
      type(node) :: nodes(:)
      real(rfp), optional :: vjacob(ndim, 2)
      type(vector) :: gqv(:)
      integer :: i, j

      do i = 1, ndim
         gqv(i) % v = cf % avec(i, 1) * lcell % qv % v + cf % avec(i, 2) * rcell % qv % v
      end do
      do j = 1, size(nodes)
         do i = 1, ndim
            gqv(i) % v = gqv(i) % v + cf % avec(i, j + 2) * nodes(j) % qv % v
         end do
      end do

      if (present(vjacob)) vjacob = cf % avec(:, 1:2)
   end subroutine face_gradient_dli1

   subroutine cell_gradient(cells, faces, nodes)
      type(cell), pointer :: cells(:), cl, cr
      type(face), pointer :: faces(:), cf
      type(node) :: nodes(:)
      type(vector) :: qv
      real(rfp) :: rl, rr, vecn(ndim)
      integer :: i, j
      !
      do i = 1, size(cells)
         cells(i) % gradient = zero
      end do
      !
      do i = 1, size(faces)
         cf => faces(i)
         cl => cf % left_cell
         cr => cf % right_cell

         vecn = cf % vecn * cf % area
         !   rl = one / sqrt(sum((cf%centp - cl%centp)**2))
         !   rr = one / sqrt(sum((cf%centp - cr%centp)**2))
         !   rl = rl / (rl + rr)
         !   rr = one - rl
         !   qv = cl%qv * rl + cr%qv * rr
         qv = sum_qv_node(nodes(cf % f2n)) *(one / real(size(cf % f2n), rfp))

         do j = 1, ndim
            cl % gradient(:, j) = cl % gradient(:, j) + qv % v * vecn(j)
            if(cf % itype == 0) &
               cr % gradient(:, j) = cr % gradient(:, j) - qv % v * vecn(j)
         end do
      end do
      !
      do i = 1, size(cells)
         cells(i) % gradient = cells(i) % gradient / cells(i) % vol
         if (s_iaxis > 0 .and. ndim == 2) &
            cells(i) % gradient(:, s_iaxis) = cells(i) % gradient(:, s_iaxis) - &
            cells(i) % qv % v / cells(i) % centp(s_iaxis)
      end do
   end subroutine cell_gradient

   subroutine face_gradient_orth(cf, cr, cl, nodes, gqv, vjacob)
      type(face), pointer :: cf
      type(cell) :: cl, cr
      type(node) :: nodes(:)
      real(rfp), optional :: vjacob(ndim, 2)
      type(vector) :: gqv(:)
      real(rfp) :: fac, fac1, vecn(ndim)
      integer :: i, j

      if (cf % itype == 0 .or.cf % itype >= partition_face_no) then
         fac = cl % vol / (cl % vol + cr % vol)
         fac1 = one - fac
         do i = 1, ndim
            gqv(i) % v = cl % gradient(:, i) * fac + cr % gradient(:, i) * fac1
         end do
      else
         !   fac = one;fac1=zero
         do i = 1, ndim
            gqv(i) % v = cl % gradient(:, i)
         end do
      end if
      fac = one / sqrt(sum((cr % centp - cl % centp)**2))
      vecn = (cr % centp - cl % centp) * fac
      !
      do i = 1, ndim
         gqv(i) = (cr % qv - cl % qv) * (fac * vecn(i)) + gqv(i) * (cf % vecn(i) - vecn(i))
      end do

      if (present(vjacob)) then
         vecn = vecn * fac
         vjacob(:, 1) = -vecn ! left cell
         vjacob(:, 2) = vecn ! right cell
         !      vjacob(:,1) = cl%vj * fac   ! left cell
         !      vjacob(:,2) = cr%vj * fac1   ! right cell
      end if
   end subroutine face_gradient_orth

   !
   !********************************************************
   ! Weighted Averaging Procedure 
   ! to determined the solution at each node
   ! qn = sum(Wci.qci)/sum(Wci)
   ! Wci = 1-Cx(Xci-Xn)-Cy(Yci-Yn)-Cz(Zci-Zn)
   ! [C] = - [A]^-1.[R]
   ! [R] =   Rci-Rn            R = X,Y,Z
   ! [A] =   (Rci-Rn).(Rcj-Rn) R = X,Y,Z,i = 1..3,j = 1..3

   subroutine weighted_average(cells, nodes, faces, interf)
      !********************************************************
      implicit none
      !
      type(face), pointer :: faces(:)
      type(cell), pointer :: cells(:), cc
      type(node), pointer :: nodes(:)
      type(itf) :: interf
      !  
      ! local variables
      real(kind = rfp) :: a(ndim, ndim, size(nodes)), dr(ndim, size(nodes))
      real(kind = rfp) :: dn(size(nodes))
      real(kind = rfp) :: rn(ndim), rc(ndim), dd(ndim), ai(ndim, ndim), di, weight
      integer, pointer :: c2n(:)
      integer :: i, j, n, k1, k2
      !  
      ! search max distance surrounding the nodes
      dn = zero
      do i = 1, size(cells)
         c2n => cells(i) % c2n
         rc = cells(i) % centp
         !
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dn(n) = max(dn(n), sum((rc - rn)*(rc - rn)))
         end do ! end loop nodes of cell
      end do ! end loop cells
      !
      ! for interface cells
      !
      do i = 1, size(interf % pcell)
         c2n => interf % pcell(i) % c2n
         rc = interf % pcell(i) % centp
         !
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dn(n) = max(dn(n), sum((rc - rn)*(rc - rn)))
         end do ! end loop nodes of cell
      end do ! end loop cells

      ! for boundary nodes
      do i = interf % nitf + 1, size(faces)
         if (faces(i) % itype == 0) exit
         c2n => faces(i) % f2n
         cc => faces(i) % right_cell
         rc = cc % centp
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dn(n) = max(dn(n), sum((rc - rn)*(rc - rn)))
            !
         end do ! end loop nodes of cell
         !
      end do ! end loop faces

      dn = sqrt(dn)
      a = zero
      dr = zero
      !
      do i = 1, size(cells)
         c2n => cells(i) % c2n
         rc = cells(i) % centp
         !
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dd = (rc - rn) / dn(n)
            dr(:, n) = dr(:, n) + dd
            do k1 = 1, ndim
               do k2 = 1, ndim
                  a(k1, k2, n) = a(k1, k2, n) + dd(k1) * dd(k2)
               end do
            end do
         end do ! end loop nodes of cell
      end do ! end loop cells
      !
      ! for interface cells
      !
      do i = 1, size(interf % pcell)
         c2n => interf % pcell(i) % c2n
         rc = interf % pcell(i) % centp
         !
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dd = (rc - rn) / dn(n)
            dr(:, n) = dr(:, n) + dd
            do k1 = 1, ndim
               do k2 = 1, ndim
                  a(k1, k2, n) = a(k1, k2, n) + dd(k1) * dd(k2)
               end do
            end do
         end do ! end loop nodes of cell
      end do ! end loop cells

      ! for boundary nodes
      do i = interf % nitf + 1, size(faces)
         if (faces(i) % itype == 0) exit
         c2n => faces(i) % f2n
         cc => faces(i) % right_cell
         rc = cc % centp
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dd = (rc - rn)/dn(n)
            dr(:, n) = dr(:, n) + dd
            do k1 = 1, ndim
               do k2 = 1, ndim
                  a(k1, k2, n) = a(k1, k2, n) + dd(k1) * dd(k2)
               end do
            end do
         end do ! end loop nodes of cell
         !
      end do ! end loop faces

      !   
      do i = 1, size(nodes)
         !
         ai = a(:,:, i)
         call direct_inv(ai)
         dr(:, i) = matmul(ai, dr(:, i)) ! A^-1.R
      end do
      !
      ! use a(1,1,:) to store sum of weights
      a(1, 1,:) = zero
      do i = 1, size(cells)
         c2n => cells(i) % c2n
         rc = cells(i) % centp
         !
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dd = (rc - rn)/dn(n)
            weight = one - dot_product(dr(:, n), dd)
            if (weight < mytiny) then
               weight = mytiny
            else if (weight > two) then
               weight = two
            end if
            cells(i) % weight(j) = weight
            a(1, 1, n) = a(1, 1, n) + weight
         end do
      end do

      do i = 1, size(interf % pcell)
         cc => interf % pcell(i)
         c2n => cc % c2n
         rc = cc % centp
         !
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dd = (rc - rn)/dn(n)
            weight = one - dot_product(dr(:, n), dd)
            if (weight < mytiny) then
               weight = mytiny
            else if (weight > two) then
               weight = two
            end if
            cc % weight(j) = weight
            a(1, 1, n) = a(1, 1, n) + weight
         end do
      end do

      do i = interf % nitf + 1, size(faces)
         if (faces(i) % itype == 0) exit
         c2n => faces(i) % f2n
         cc => faces(i) % right_cell
         rc = cc % centp
         !
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dd = (rc - rn)/dn(n)
            weight = one - dot_product(dr(:, n), dd)
            if (weight < mytiny) then
               weight = mytiny
            else if (weight > two) then
               weight = two
            end if
            cc % weight(j) = weight
            a(1, 1, n) = a(1, 1, n) + weight
         end do
      end do
      !
      ! normalize weights
      !
      do i = 1, size(cells)
         c2n => cells(i) % c2n
         do j = 1, size(c2n)
            n = c2n(j)
            cells(i) % weight(j) = cells(i) % weight(j) / a(1, 1, n)
         end do
      end do
      !
      do i = 1, size(interf % pcell)
         cc => interf % pcell(i)
         c2n => cc % c2n
         do j = 1, size(c2n)
            n = c2n(j)
            cc % weight(j) = cc % weight(j) / a(1, 1, n)
         end do
      end do

      ! boundary nodes
      ! 
      do i = interf % nitf + 1, size(faces)
         if (faces(i) % itype == 0) exit
         c2n => faces(i) % f2n
         cc => faces(i) % right_cell
         do j = 1, size(c2n)
            n = c2n(j)
            cc % weight(j) = cc % weight(j) / a(1, 1, n)
         end do
      end do

   end subroutine weighted_average

   subroutine inverse_average(cells, nodes, faces, interf)
      !********************************************************
      implicit none
      !
      type(cell), pointer :: cells(:), cc
      type(node), pointer :: nodes(:)
      type(face), pointer :: faces(:)
      type(itf) :: interf
      ! local variables
      real(kind = rfp) :: dr(size(nodes))
      real(kind = rfp) :: rn(ndim), rc(ndim), dd
      integer, pointer :: c2n(:)
      integer :: i, j, n
      !  
      dr = zero
      !
      do i = 1, size(cells)
         c2n => cells(i) % c2n
         rc = cells(i) % centp
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dd = one / sqrt(dot_product(rc - rn, rc - rn))
            cells(i) % weight(j) = dd
            dr(n) = dr(n) + dd
         end do ! end loop nodes of cell

      end do ! end loop cells
      !
      do i = 1, size(interf % pcell)
         cc => interf % pcell(i)
         c2n => cc % c2n
         rc = cc % centp
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dd = one / sqrt(dot_product(rc - rn, rc - rn))
            cc % weight(j) = dd
            dr(n) = dr(n) + dd
         end do ! end loop nodes of cell

      end do ! end loop cells
      !
      ! boundary nodes

      do i = interf % nitf + 1, size(faces)
         if (faces(i) % itype == 0) exit
         c2n => faces(i) % f2n
         rc = faces(i) % right_cell % centp
         do j = 1, size(c2n)
            n = c2n(j)
            rn = nodes(n) % xyz
            dd = one / sqrt(dot_product(rc - rn, rc - rn))
            faces(i) % right_cell % weight(j) = dd ! save weight in ghost cell
            dr(n) = dr(n) + dd
         end do ! end loop nodes of face
      end do ! end loop faces
      !   
      ! normalize weights
      !
      do i = 1, size(cells)
         c2n => cells(i) % c2n
         do j = 1, size(c2n)
            n = c2n(j)
            cells(i) % weight(j) = cells(i) % weight(j) / dr(n)
         end do
      end do
      !
      do i = 1, size(interf % pcell)
         cc => interf % pcell(i)
         c2n => cc % c2n
         do j = 1, size(c2n)
            n = c2n(j)
            cc % weight(j) = cc % weight(j) / dr(n)
         end do
      end do
      !   
      ! normalize weights of face
      !
      do i = interf % nitf + 1, size(faces)
         if (faces(i) % itype == 0) exit
         c2n => faces(i) % f2n
         do j = 1, size(c2n)
            n = c2n(j)
            faces(i) % right_cell % weight(j) = faces(i) % right_cell % weight(j) / dr(n)
         end do
      end do

   end subroutine inverse_average
   !

   !
   !********************************************************
   ! Linear least square reconstruction
   !
   subroutine reconstruction(cells, faces)
      !********************************************************
      implicit none
      !
      type(cell), pointer :: cells(:)
      type(face), intent(in) :: faces(:)
      ! local variables
      type(cell), pointer :: current_cell, ghost_cell
      type(face), pointer :: current_face
      ! variables at face
      real(kind = rfp) :: distance(ndim), current_center(ndim), wi(20), srate
      integer :: i, j, k1, k2, l, nrec, ierr
      !  
      ! Use SVD algorithm to solve least square equation
      ! A(MxN)x= b   A = U.w.v^T  A^-1 = v.1/w.U^T  
      !
      if (s_irec == 2) return
      !
      nrec = size(cells(1) % gradient, dim = 2)
      do i = 1, size(cells)
         current_cell => cells(i)
         !  initial gradient to zero
         !   
         current_cell % gradient = zero

         current_center = current_cell % centp
         ierr = 0
         !   
         10 do j = 1, size(current_cell % scell)
            !
            ! search for neighbour cell
            !
            ghost_cell => current_cell % scell(j) % to_cell
            distance = ghost_cell % centp - current_center
            !   if(neqf == 0.and.neqm > 0) then
            !    srate = e_resistivity(current_cell%qv,current_cell) - &
            !            e_resistivity(ghost_cell%qv,ghost_cell)
            !!    srate  = min(srate,one/srate)
            !!    wi(j) = srate
            !    wi(j) = one
            !    if(abs(srate) > 1.e-10_rfp) wi(j) = 1.e-30_rfp
            !   else
            wi(j) = one
            !   end if

            !    if(ierr == 0.and.ndim == 2) then
            !      wi(j) = one / sqrt(sum(distance**2))
            !    else
            !      wi(j) = one
            !    end if
            !     wi(j) = wi(j) * srate
            current_cell % svd(j,:ndim) = distance * wi(j)
            cycle
            if (nrec <= ndim) cycle
            l = ndim
            do k1 = 1, ndim
               do k2 = k1, ndim
                  l = l + 1
                  current_cell % svd(j, l) = distance(k1) * distance(k2)
               end do
            end do
         end do ! end loop neighbour cells
         ! use SVD to inverse SVD matrix
         !
         call svdcmp_inv(current_cell % svd, ierr)
         if (ierr == 1) then
            goto 10 ! redo the matrix
         else if (ierr > 1) then
            print *, ' Sorry, I can not handle your grid in reconstruction'
            stop
         end if
         !attantion:
         !svd store A(not A^T) when interpolate, use
         ! q = matmul(b,A) (not matmul(A,b))
         !
         !  modified to weighted least square
         !
         do j = 1, size(current_cell % scell)
            current_cell % svd(j,:) = current_cell % svd(j,:) * wi(j)
         end do
         !
      end do ! end loop cells

   end subroutine reconstruction

   !********************************************************
   ! Calculate qv at nodes
   !
   subroutine cal_node_qv(cells, nodes, faces, interf, bnodes)
      !********************************************************
      implicit none
      !
      type(cell) :: cells(:)
      type(cell), pointer :: cc
      type(node) :: nodes(:)
      type(bnode), pointer, optional :: bnodes(:)
      type(face) :: faces(:)
      type(vector) :: qv, gqv
      type(itf) :: interf
      ! local variables
      integer, pointer :: c2n(:)
      integer :: i, k, j, n
      real(rfp) :: t, r, w, wq, cls

      do k = 1, neq
         nodes % qv % v(k) = zero
      end do

      ! contribution from  cells
      do i = 1, size(cells)

         c2n => cells(i) % c2n
         qv = cells(i) % qv

         do j = 1, size(c2n)
            n = c2n(j)        
            nodes(n) % qv = nodes(n) % qv + qv * cells(i) % weight(j)
         end do  
      end do

      ! contribution from pcell
      do i = 1, size(interf % pcell)
         cc => interf % pcell(i)
         c2n => cc % c2n
         qv = cc % qv
         do j = 1, size(c2n)
            n = c2n(j)
            nodes(n) % qv = nodes(n) % qv + qv * cc % weight(j)
         end do
      end do

      ! contribution from rc
      do i = interf % nitf + 1, size(faces)
         if (faces(i) % itype == 0) exit
         cc => faces(i) % right_cell
         qv = cc % qv
         c2n => faces(i) % f2n
         do j = 1, size(c2n)
            n = c2n(j)
            nodes(n) % qv = nodes(n) % qv + qv * cc % weight(j)
         end do

         ! correct for no-slip wall
         !
         !    if(bc(faces(i)%itype)%igrp == wall_type) then
         !    if(bc(faces(i)%itype)%itype == 5) cycle  ! rotation
         !    if(bc(faces(i)%itype)%itype > 0 ) then  ! no-slip wall
         !     c2n => faces(i)%f2n
         !     do j = 1, size(c2n)
         !      n = c2n(j)
         !      nodes(n)%qv%v(imb:ime) = bc(faces(i)%itype)%vars(:ndim)
         !     end do
         !    end if
         !   end if
         !  
      end do

      !   Special treatment for the nodes at the interface
      !
      if(present(bnodes)) then
         do i = 1, size(bnodes)
            T = zero; w = zero; wq = zero
            qv % v = zero
            do j = 1, bnodes(i) % n
               cc => bnodes(i) % scell(j) % to_cell
               r = one / sqrt(sum((cc % centp - bnodes(i) % tn % xyz)**2))
               !     if(nmeq > 0.and.cc%itype /= partition_face_no) then
               if (neqf > 0) then
                  if (cc % itype /= partition_face_no.or.cc % srf > zero) then
                     qv % v(:neqf) = qv % v(:neqf) + cc % qv % v(:neqf) * r
                     wq = wq + r
                  end if
               end if
               if (bnodes(i) % itype == 1) then
                  r = r * lamda(cc % qv, cc)
                  T = T + r * cc % qv % v(ien)
                  w = w + r
               end if
            end do
            if (bnodes(i) % itype == 1) T = T * (one / w)
            if (neqf > 0) then
               qv % v(:neqf) = qv % v(:neqf) * (one / wq)
               if (nmeq > 0.and.s_ivis > 0) qv % v(imb:ime) = zero
            end if
            if (bnodes(i) % itype == 1) qv % v(ien) = T
            bnodes(i) % tn % qv % v(:neqf) = qv % v(:neqf)
         end do
      end if
   end subroutine cal_node_qv

   !
   !********************************************************
   ! Calculate gradient matrix from 
   ! linear least square reconstruction
   !
   subroutine grad_recon(cells, faces, nodes)
      !********************************************************
      implicit none
      !
      type(node), intent(in) :: nodes(:)
      type(cell), pointer :: cells(:)
      type(face), pointer :: faces(:), fs
      ! local variables
      type(cell), pointer :: current_cell, ghost_cell
      type(face), pointer :: current_face, ghost_face
      type(vector) :: qv, qvn, qvmax, qvmin
      real(rfp) :: cc(ndim), cf(ndim), phi(neq), phid(neq), epsi
      !
      ! variables for computation
      integer :: i, j, k, k1, k2

      !  
      if (s_irec == 2) return
      if (s_irec == 3) then ! donot use it
         !
         do i = 1, size(cells)
            current_cell => cells(i)
            current_cell % gradient(:,:ndim) = zero
         end do

         do i = 1, size(faces)
            fs => faces(i)
            cc = fs % vecn * fs % area
            qv % v = zero
            do k = 1, neq
               qv % v(k) = sum(nodes(fs % f2n) % qv % v(k))
            end do
            qv % v = qv % v / real(size(fs % f2n), rfp)
            do k = 1, ndim
               current_cell % gradient(:, k) = current_cell % gradient(:, k) + &
                  qv % v * cc(k)
            end do
         end do
         do i = 1, size(cells)
            current_cell => cells(i)
            epsi = one / current_cell % vol
            current_cell % gradient(:,:) = current_cell % gradient(:,:) * epsi
         end do
      else if (s_irec == 2) then
         !
         do i = 1, 0 !size(cells)
            current_cell => cells(i)
            current_cell % gradient(:,:ndim) = &
               current_cell % gradient(:,:ndim) / real(size(current_cell % sface), rfp)
            if (abs(s_ialg) <= 1) cycle
            j = ndim
            do k1 = 1, ndim
               do k2 = k1, ndim
                  j = j + 1
                  if (k1 == k2) then
                     current_cell % gradient(:, j) = half * current_cell % gradient(:, j) &
                        / current_cell % vol
                  else
                     current_cell % gradient(:, j) = current_cell % gradient(:, j) &
                        / current_cell % vol
                  end if
               end do
            end do
         end do ! end loop cells
      else
         do i = 1, size(cells)
            current_cell => cells(i)
            qv = current_cell % qv
            current_cell % gradient(:,:) = zero


            do j = 1, size(current_cell % scell)
               !
               ! search for neighbour cell
               !
               ghost_cell => current_cell % scell(j) % to_cell
               qvn = ghost_cell % qv

               do k = 1, s_nrec
                  current_cell % gradient(:, k) = current_cell % gradient(:, k) + &
                     (qvn % v - qv % v) * current_cell % svd(j, k)
               end do

            end do ! end loop neighbor
            !
         end do ! end loop cells

      end if

      !  return
      if (s_ialg >= 0.or.s_irec /= 1) return
      !
      !   dli's limiter
      !   call try_new_limiter(cells,faces,nodes)
      !   return
      !
      ! Barth Limiter
      !   
      call barth_limiter(cells, faces, nodes)

   end subroutine grad_recon

   subroutine dli_limiter(cc, keq, qmax, qmin)
      type(cell), pointer :: cc, gc
      real(rfp) :: gp(ndim), y0, p0(ndim), qmax, qmin
      real(rfp) :: a(ndim, ndim), f(ndim), sigma, d, f1, f2, lamda, dmin
      real(rfp), pointer :: pn(:,:), yn(:), en(:), q(:)
      integer :: keq, i, j, n, k
      !
      y0 = cc % qv % v(keq)
      p0 = cc % centp
      gp = cc % gradient(keq,:)
      !
      n = size(cc % scell)
      allocate(pn(ndim, n), yn(n), en(n), q(n))
      !
      lamda = 1.0e5_rfp
      sigma = zero
      do i = 1, n
         gc => cc % scell(i) % to_cell
         pn(:, i) = gc % centp - p0
         yn(i) = gc % qv % v(keq) - y0
         q(i) = sum(gp * pn(:, i))
         en(i) = q(i) - yn(i)
         sigma = sigma + en(i) * en(i)
      end do
      sigma = sqrt(sigma/real(n, rfp))
      en = en / sigma
      pn = pn / sigma
      !
      f = zero; a = zero
      do j = 1, ndim
         do i = 1, n
            !    d = one / (one + en(i)*en(i))
            !    f(j) = f(j) + two * en(i) * pn(j,i) * d
            !    a(j,:) = a(j,:) + two * d * pn(j,i) * pn(:,i) *(one - two * (en(i) * d)**2)

            f(j) = f(j) + en(i) * pn(j, i)
            a(j,:) = a(j,:) + pn(j, i) * pn(:, i)
            !
            if (q(i) > qmax) then
               f(j) = f(j) + lamda * (q(i) - qmax)/ sigma * pn(j, i)
               a(j,:) = a(j,:) + lamda * pn(j, i) * pn(:, i)
            else if (q(i) < qmin) then
               f(j) = f(j) + lamda * (q(i) - qmin)/ sigma * pn(j, i)
               a(j,:) = a(j,:) + lamda * pn(j, i) * pn(:, i)
            end if
         end do
      end do
      p0 = maxval(abs(a), dim = 1)
      if (any(p0 == zero)) then
         print *, 'singular matrix', a
      else
         do i = 1, ndim
            a(:, i) = a(:, i) / p0(i)
         end do
         call direct_inv(a)
         gp = gp - matmul(a, f) / p0
         cc % gradient(keq,:) = gp
      end if
      !
      deallocate(pn, yn, en, q)

   end subroutine dli_limiter

   !
   function psif(dp, dm)
      real(rfp) :: y, y1, psif, dp, dm
      !  r = 10.0 * (vol)**(one/real(ndim,rfp))
      psif = dp/dm
      return
      y = dp / dm
      psif = (y * y + two * y) / (y * y + y + two) ! (y^2+2y)/(y^2+y+2)
   end function psif

   function qv_from_highorder(cc, cf, nodes)result(qv)
      type(node), intent(in) :: nodes(:)
      type(cell), pointer :: cc
      type(face), intent(in) :: cf
      type(face), pointer :: fs
      type(vector) :: qv, qvc !,qvn,qvmax,qvmin
      real(rfp) :: distance(ndim), pn, pmax, pmin
      real(rfp) :: dc(size(cc % gradient, dim = 2))
      real(rfp) :: c1, c2, epsi
      integer :: l, k1, k2, nrec, j
      nrec = size(cc % gradient, dim = 2)
      if (s_ialg == 0) then
         qv % v = zero
         return
      end if
      !  call optimum_coef(ndim,size(cc%c2n),size(cf%f2n),c1,c2) 
      !  qvc = sum_qv_node(nodes(cc%c2n))
      !  qvf = sum_qv_node(nodes(cf%f2n))
      !  qv = qvf * c1 - qvc * c2
      !  return

      distance = cf % centp - cc % centp
      if (nrec <= ndim) then
         qv % v = matmul(cc % gradient, distance)
      else
         dc(:ndim) = distance
         l = ndim
         do k1 = 1, ndim
            do k2 = k1, ndim
               l = l + 1
               dc(l) = distance(k1) * distance(k2)
            end do
         end do
         qv % v = matmul(cc % gradient, dc)
      end if
      if (s_ivis == 2) qv % v(ike:iom) = zero
      return
      !try some limiter
      qv = qv + cc % qv
      do k1 = 1, neq
         qv % v(k1) = min(qv % v(k1), &
            max(cf % left_cell % qv % v(k1), cf % right_cell % qv % v(k1)))
         qv % v(k1) = max(qv % v(k1), &
            min(cf % left_cell % qv % v(k1), cf % right_cell % qv % v(k1)))
      end do
      qv = qv - cc % qv
      !  
      ! for k-w equation, always using first order
      ! if(s_ivis == 2 .or. s_ialg < 0) then
      if (s_ivis == 2) then
         if (s_ivis == 2) qv % v(ike:iom) = zero
         !try some limiter for pressure
         !   if(abs(qv%v(ico)) >      &
         !     abs(cf%left_cell%qv%v(ico) - cf%right_cell%qv%v(ico))) &
         !   qv%v(ico) = zero
         ! search for neighbour cell and search max and min values
         !
         qvc = cc % qv
         pmax = qvc % v(ico)
         pmin = qvc % v(ico)
         do j = 1, size(cc % sface)
            fs => cc % sface(j) % to_face
            !
            if (associated(cc, fs % left_cell)) then
               pn = fs % right_cell % qv % v(ico)
            else
               pn = fs % left_cell % qv % v(ico)
            end if
            if (pmax < pn) pmax = pn
            if (pmin > pn) pmin = pn
            !
         end do


         !     do l = 1, neq
         !     if(l == ike .or. l == iom) cycle
         !     epsi = 0.001_rfp * qvc%v(l)
         !     if(qv%v(l) > epsi) then
         !      qv%v(l) = min(qvmax%v(l) - qvc%v(l),qv%v(l))
         !     else if(qv%v(l) < -epsi) then
         !      qv%v(l) = max(qvmin%v(l) - qvc%v(l),qv%v(l))
         !     end if
         !     end do
         !  switch off turbulent parameters
         pn = qvc % v(ico) + qv % v(ico)
         !    if(any(qvn%v > qvmax%v) .or. any(qvn%v < qvmin%v)) then
         if (pn > pmax .or. pn < pmin) then
            qv % v(ico) = zero
         end if
         !
      end if
   end function qv_from_highorder

   subroutine face_qv(qvl, qvr, cl, cr, cf)
      type(cell), pointer :: cl, cr
      type(face), intent(in) :: cf
      type(vector) :: qvl, qvr
      real(rfp) :: rl(ndim), rr(ndim)
      real(rfp) :: gn(neq), gnl(neq), gnr(neq), alpha(neq), dl, dr

      if (cf % itype /= 0) return !cr%gradient = zero
      rl = cf % centp + (cl % centp - cf % centp) * cf % vecn
      rr = cf % centp + (cr % centp - cf % centp) * cf % vecn
      !
      qvl % v = cl % qv % v + matmul(cl % gradient, rl - cl % centp)
      qvr % v = cr % qv % v + matmul(cr % gradient, rr - cr % centp)
      ! qvl%v = cl%qv%v + matmul(half*(cl%gradient+cr%gradient),cf%centp - cl%centp)
      ! qvr%v = cr%qv%v + matmul(half*(cl%gradient+cr%gradient),cf%centp - cr%centp)
      !return


      gnl = matmul(cl % gradient, cf % vecn)
      gnr = matmul(cr % gradient, cf % vecn)
      !
      alpha = one - abs(gnr - gnl)/(abs(gnl) + abs(gnr) + 1.e-30_rfp)
      !
      gn = half * (gnl + gnr) * alpha

      dl = sqrt(sum((rl - cf % centp)**2))
      dr = sqrt(sum((rr - cf % centp)**2))

      qvl % v = qvl % v + gn * dl
      qvr % v = qvr % v - gn * dr
      !

   end subroutine face_qv

   subroutine optimum_coef(nd, nv, nf, cf, cv)
      integer, intent(in) :: nd, nv, nf
      real(rfp), intent(out) :: cv, cf
      ! cv = 1/((nd+1)(nv-nf))
      ! cf = cv*nv/nf
      !
      if (nd == 2) then
         select case(nv)
         case(3)
            cv = one_third
            cf = half
         case(4)
            cv = 0.1666666666666666667_rfp
            cf = one_third
         case default
            cv = one / (real(nd + 1, rfp) * real(nv - nf, rfp))
            cf = cv * real(nv, rfp) / real(nf, rfp)
         end select
      else if (nd == 3) then
         ! Tets
         if (nv == 4 .and. nf == 3) then
            cv = 0.25_rfp
            cf = one_third
            ! prism
         else if (nv == 6 .and. nf == 3) then
            cv = 0.083333333333333_rfp
            cf = 0.166666666666667_rfp
         else if (nv == 6 .and. nf == 4) then
            cv = 0.125_rfp
            cf = 0.1875_rfp
         else if (nv == 8 .and. nf == 4) then
            cv = 0.0625_rfp
            cf = 0.125_rfp
         else
            cv = one / (real(nd + 1, rfp) * real(nv - nf, rfp))
            cf = cv * real(nv, rfp) / real(nf, rfp)
         end if
      end if

   end subroutine optimum_coef

   function sum_qv_node(nodes)result(qv)
      type(node), intent(in) :: nodes(:)
      type(vector) :: qv
      integer :: i
      qv % v = zero
      do i = 1, size(nodes)
         qv = qv + nodes(i) % qv
      end do
   end function sum_qv_node

   subroutine barth_limiter(cells, faces, nodes)
      !********************************************************
      implicit none
      !
      type(node), intent(in) :: nodes(:)
      type(cell), pointer :: cells(:)
      type(face), pointer :: faces(:), fs, ghost_face
      ! local variables
      type(cell), pointer :: current_cell, ghost_cell
      type(face), pointer :: current_face
      type(vector) :: qv, qvn, qvmax, qvmin
      real(rfp) :: cc(ndim), cf(ndim), phi(neq), phid(neq), epsi
      !
      ! variables for computation
      integer :: i, j, k, n


      ! search for max and min value surounding a cell
      !

      10 epsi = 1.e-30_rfp
      do i = 1, size(cells)
         !
         current_cell => cells(i)

         cc = current_cell % centp
         qv = current_cell % qv
         qvmax = qv
         qvmin = qv

         ! search for neighbour cell and search max and min values
         !
         do j = 1, size(current_cell % scell)

            ghost_cell => current_cell % scell(j) % to_cell
            qvn = ghost_cell % qv

            where(qvmax % v < qvn % v) 
               qvmax % v = qvn % v
            endwhere
            where(qvmin % v > qvn % v)
               qvmin % v = qvn % v
            endwhere
         end do
         !
         phi = one
         do j = 1, size(current_cell % sface)
            fs => current_cell % sface(j) % to_face !%centp
            qvn = qv_from_highorder(current_cell, fs, nodes)

            do k = 1, neq
               if (abs(qvn % v(k)) <= epsi) then
                  phid(k) = one
               else if (qvn % v(k) > epsi) then
                  phid(k) = psif((qvmax % v(k) - qv % v(k)), qvn % v(k))
               else
                  phid(k) = psif((qvmin % v(k) - qv % v(k)), qvn % v(k))
               end if
               phi(k) = min(phi(k), phid(k))
            end do
            !
         end do ! loop faces

         do k = 1, neq
            current_cell % gradient(k,:) = current_cell % gradient(k,:) * phi(k)
         end do

      end do ! end loop cells

   end subroutine barth_limiter


   subroutine try_new_limiter(cells, faces, nodes)
      !********************************************************
      implicit none
      !
      type(node), intent(in) :: nodes(:)
      type(cell), pointer :: cells(:)
      type(face), pointer :: faces(:), fs
      ! local variables
      type(cell), pointer :: current_cell, ghost_cell
      type(face), pointer :: current_face
      type(vector) :: qv, qvn, qvmax, qvmin
      real(rfp) :: cc(ndim), cf(ndim), phi(neq), phid(neq), epsi
      !
      ! variables for computation
      integer :: i, j, k

      do i = 1, size(cells)
         !
         current_cell => cells(i)
         cc = current_cell % centp
         qv = current_cell % qv
         qvmax = qv
         qvmin = qv

         do j = 1, size(current_cell % scell)
            ghost_cell => current_cell % scell(j) % to_cell
            qvn = ghost_cell % qv
            where(qvmax % v < qvn % v)
               qvmax % v = qvn % v
            end where
            where(qvmin % v > qvn % v)
               qvmin % v = qvn % v
            end where
         end do
         qvmax = qvmax - qv
         qvmin = qvmin - qv
         !
         phi = zero
         do j = 1, size(current_cell % scell)
            ghost_cell => current_cell % scell(j) % to_cell
            qvn % v = matmul(current_cell % gradient, ghost_cell % centp - cc)
            where(qvn % v > qvmax % v)
               phi = two
            end where
            where(qvn % v < qvmin % v) 
               phi = two
            end where
         end do
         !         if(any(phi > one)) k1=k1+1
         !         if(any(phi > one)) current_cell%gradient = zero
         !	 return

         do k = 1, neq
            if (phi(k) > one) then ! need limiter
               call dli_limiter(current_cell, k, qvmax % v(k), qvmin % v(k))
            end if
         end do

      end do ! end loop cells

   end subroutine try_new_limiter
end module gems_disc

