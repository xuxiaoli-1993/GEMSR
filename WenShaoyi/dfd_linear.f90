MODULE dfd_linear
   use dfd_library
   use dfd_data
   use dfd_bound
contains
   !
   subroutine linear_solve(cells, faces, nodes, interf, raw)
      !********************************************************
      implicit none
      !
      type(cell), pointer :: cells(:)
      type(node), pointer :: nodes(:)
      type(face), pointer :: faces(:)
      type(itf) :: interf
      type(raw_chain) :: raw(ndim)
      ! local variables
      type(cell), pointer :: current_cell, ghost_cell
      type(face), pointer :: current_face
      ! variables at face
      type(matrix), pointer :: dm
      type(vector), pointer :: qv, dqv, res
      type(vector) :: qvs, deltv, sumdq
      integer :: i, j, k, istatus
      !  
      !*** loop all colors to calculate LHS
      ! do icolor = 1, size(colors)
      !  current_color = colors(icolor)
      !
      !  call normalize_system

      select case(s_imeth)
      case(0) ! explicit
         return
      case(1) ! Block (point) Gauss seidel
         call point_gauss_seidel
      case(2) ! Block DLU  
         call block_DLU
      case(3) ! Block GMRES with (point) Gauss seidel preconditioning
         call gmres_dilu
      case(4) ! CG with DLU preconditioning
         call bicgstab_dilu
      case(5)
         call dl_lgs ! Ding's LGS
      case(6)
         call gmres_lgs
      case(7)
         !    call smooth_residual
         call dl_lgs ! Ding's LGS
      end select

      nullify(current_cell, ghost_cell, current_face, res, dqv)

   contains

      subroutine normalize_system
         integer :: i
         type(cell), pointer :: cc
         type(face), pointer :: cf
         real(rfp) :: volinv
         do i = 1, size(cells)
            cc => cells(i)
            volinv = one / cc % vol
            cc % dm = cc % dm * volinv
            cc % res = cc % res * volinv
         end do
         do i = 1, size(faces)
            cf => faces(i)
            cc => cf % left_cell
            volinv = one / cc % vol
            cf % ajl = cf % ajl * volinv
            !   if(cf%itype /= 0 ) cycle
            cc => cf % right_cell
            volinv = one / cc % vol
            cf % ajr = cf % ajr * volinv
         end do
      end subroutine normalize_system

      subroutine smooth_residual
         integer :: i, j, k, nc
         type(cell), pointer :: cc, cn
         type(face), pointer :: cf
         do i = 1, size(cells)
            cc => cells(i)
            nc = size(cc % sface)
            k = 0
            do j = 1, nc
               cf => cc % sface(j) % to_face
               if (cf % itype /= 0) cycle
               k = k + 1
               cn => cf % left_cell
               if (associated(cc, cn)) cn => cf % right_cell
               cc % res = cc % res + cn % res
            end do
            cc % res = cc % res * (one / real(k, rfp))
         end do
      end subroutine smooth_residual

      subroutine dl_lgs
         implicit none
         !
         type(chain), pointer :: bchain
         type(matrix) :: dm
         type(matrix), pointer :: cm(:)
         type(vector), pointer :: res(:)
         real(rfp) :: ar
         integer :: l, nfb, n
         !
         ! initial dqv
         !
         !   do i = 1, size(cells)
         !!    dm = cells(i)%dm
         !!    cells(i)%dqv = cells(i)%res
         !!    call normalize(dm%e,cells(i)%dqv%v)
         !    cells(i)%dqv%v = zero
         !   end do

         !   do i = 1, interf%nitf
         !    ghost_cell => faces(i)%right_cell
         !    ghost_cell%dqv%v = zero
         !   end do
         !
         do l = 1, s_isub

            do k = 1, ndim
               nfb = size(raw(k) % bchain)
               cells % srf = zero
               do i = 1, nfb
                  n = raw(k) % nchain(i)
                  bchain => raw(k) % bchain(i)
                  if (n <= 1) then
                     call pgsforlgs(bchain % cc)
                     cycle
                  end if
                  allocate(cm(n), res(n))
                  call forward_sweep(cm, res, bchain, n)
                  deallocate(cm, res)
               end do
               !
               !    call pgsforlgs(1,1,k) 
               !    
               !    if(s_isub == 1) cycle
               !    call update_interface_dqv(nodes,cells,interf)
               !    cycle   ! skip reverse loop
               !    
               cells % srf = zero
               do i = nfb, 1, -1
                  n = raw(k) % nchain(i)
                  bchain => raw(k) % bchain(i)
                  if (n <= 1) then
                     call pgsforlgs(bchain % cc)
                     cycle
                  end if
                  allocate(cm(n), res(n))
                  bchain => raw(k) % bchain(i)
                  call forward_sweep(cm, res, bchain, n)
                  deallocate(cm, res)
               end do
               !
               !     call update_interface_dqv(nodes,cells,interf)
               !
               !   call pgsforlgs(1,1,k) 

            end do ! end dimension
            !
         end do ! end subiterations

         return

         istatus = 0
         do i = 1, size(cells)
            ! check the high aspact ratio
            cells(i) % srf = one
            current_face => current_cell % sface(1) % to_face
            !   if(ndim == 2) then
            ar = current_face % area**ndim / cells(i) % vol
            !   else
            !    ar = rt(current_face%area) / (cells(i)%vol**one_third)
            !   end if
            if (ar > 500.0_rfp) then
               call matrixinv(cells(i) % dm % e)
               !     istatus = istatus + 1
               cells(i) % srf = 0
            end if
         end do
         !   print *,'there are',istatus, 'cells AR > 100'

         !   call pgsforlgs(s_isub,0,1)

      end subroutine dl_lgs


      subroutine pgsforlgs(cc)
         implicit none
         type(cell), pointer :: cc
         integer :: i, j, k, m
         type(matrix) :: dm
         current_cell => cc
         deltv % v = zero
         do j = 1, size(current_cell % sface) ! loop faces surrounding the cell
            current_face => current_cell % sface(j) % to_face
            if (current_face % itype /= 0) cycle
            ghost_cell => current_face % left_cell
            if (associated(ghost_cell, current_cell)) then
               ghost_cell => current_face % right_cell
               dqv => ghost_cell % dqv
               deltv = deltv + current_face % ajr * dqv
            else
               dqv => ghost_cell % dqv
               deltv = deltv + current_face % ajl * dqv
            end if
         end do ! end loop faces
         res => current_cell % res
         dqv => current_cell % dqv
         dm = current_cell % dm
         dqv = res - deltv
         call normalize(dm % e, dqv % v)
         !
      end subroutine pgsforlgs

      subroutine forward_sweep(chain_cm, chain_res, bchain, nc)
         implicit none
         !
         type(chain), pointer :: bchain, cchain
         type(cell), pointer :: cc, cn
         type(face), pointer :: cfl, cfr, cff
         type(neighbour), pointer :: sf(:)
         ! variables at face
         type(matrix) :: chain_cm(:), dm, al, ar
         type(vector) :: chain_res(:), res
         integer :: nc, i, il, ir, ichain
         !
         cchain => bchain
         do ichain = 1, nc
            cc => cchain % cc
            !    cc%srf = one
            cfl => cchain % cf
            cfr => cchain % next_chain % cf
            !  
            cn => cfl % left_cell
            if (associated(cc, cn)) then
               al = cfl % ajr
            else
               al = cfl % ajl
            end if
            cn => cfr % left_cell
            if (associated(cc, cn)) then
               cn => cfr % right_cell
               ar = cfr % ajr
            else
               ar = cfr % ajl
            end if
            res = cc % res
            sf => cc % sface
            do i = 1, size(sf)
               cff => sf(i) % to_face
               if (associated(cff, cfl)) cycle
               if (associated(cff, cfr)) cycle
               if (cff % itype /= 0.and.cff % itype >= partition_face_no) cycle !skip boundary face
               !    if(cff%itype /= 0 ) cycle  !skip boundary face
               cn => cff % left_cell
               if (associated(cc, cn)) then
                  cn => cff % right_cell
                  res = res - cff % ajr * cn % dqv
               else
                  res = res - cff % ajl * cn % dqv
               end if
            end do

            if (ichain == 1) then ! the begin cell for sweep
               dm = cc % dm
               call matrixinv(dm % e)
               chain_cm(ichain) = dm * ar
               chain_res(ichain) = dm * res
            else if (ichain == nc) then
               dm = cc % dm - al * chain_cm(nc - 1)
               res = res - al * chain_res(nc - 1)
               call matrixinv(dm % e)
               chain_res(nc) = dm * res
            else
               dm = cc % dm - al * chain_cm(ichain - 1)
               res = res - al * chain_res(ichain - 1)
               call matrixinv(dm % e)
               chain_cm(ichain) = dm * ar
               chain_res(ichain) = dm * res
            end if
            cchain => cchain % next_chain
         end do
         cchain => cchain % up_chain
         cc => cchain % cc
         cc % dqv = chain_res(nc)
         do ichain = nc - 1, 1, -1
            cchain => cchain % up_chain
            cc => cchain % cc
            chain_res(ichain) = chain_res(ichain) - &
               chain_cm(ichain) * chain_res(ichain + 1)
            cc % dqv = chain_res(ichain)
         end do
      end subroutine forward_sweep

      subroutine point_gauss_seidel
         ! initial dqv 
         do i = 1, size(cells)
            current_cell => cells(i)
            dm => current_cell % dm
            call matrixinv(dm % e)
            current_cell % dqv = dm * current_cell % res
            !    current_cell%dqv%v = zero


         end do
         !


         do k = 1, s_isub
            do i = 1, size(cells)
               current_cell => cells(i)
               deltv % v = zero
               do j = 1, size(current_cell % sface) ! loop faces surrounding the cell
                  current_face => current_cell % sface(j) % to_face
                  if (current_face % itype /= 0) cycle
                  ghost_cell => current_face % left_cell
                  if (associated(ghost_cell, current_cell)) then
                     ghost_cell => current_face % right_cell
                     dqv => ghost_cell % dqv
                     deltv = deltv + current_face % ajr * dqv
                     !print *,j,'right'
                     !print *,'dqv',dqv
                     !print *,'dv',deltv
                  else
                     dqv => ghost_cell % dqv
                     deltv = deltv + current_face % ajl * dqv
                     !print *,j,'left'
                     !print *,'dqv',dqv
                     !print *,'dv',deltv
                  end if


               end do ! end loop faces

               res => current_cell % res
               dqv => current_cell % dqv
               dqv = current_cell % dm * (res - deltv)

               if (isnan(dqv % v(1))) then
                  print *, 'id', i
                  print *, dqv % v(:)
                  stop


               end if
            end do ! end loop cells
            !       call update_interface_dqv(nodes,cells,interf)
         end do ! end iterations

      end subroutine point_gauss_seidel
      !
      subroutine block_dlu
         !
         !  call my_pivot
         !
         do i = 1, size(cells)
            call matrixinv(cells(i) % dm % e)
         end do
         !  
         call my_solve_lu(cells % res)

      end subroutine block_DLU

      !  subroutine block_dilup
      !
      !  integer::ip = 3

      !  call my_dilup(ip)
      !
      !  call my_solve_lu(cells%res)

      !  end subroutine block_DLU

      subroutine gmres_lgs
         real(rfp) :: norm, toler
         integer :: k, j

         ! normalize linear eqs
         !
         !  do k = 1, neq
         !  norm = maxval(abs(cells%res%v(k)))
         !  if(norm /= zero) then
         !   norm = one / norm
         !   cells%res%v(k) = cells%res%v(k) * norm
         !   do j = 1, neq
         !    cells%dm%e(k,j) = cells%dm%e(k,j) * norm
         !    faces%ajl%e(k,j) = faces%ajl%e(k,j) * norm
         !    faces%ajr%e(k,j) = faces%ajr%e(k,j) * norm   
         !   end do
         !  end if
         !  end do  

         toler = sqrt(my_dot_product(cells % res, cells % res))
         toler = toler * 1.0e - 5_rfp
         !
         call my_solve_lgs(cells % res)
         !
         call my_gmres_lgs(toler, cells % dm)
         !
      end subroutine gmres_lgs

      subroutine my_solve_lgs(res)
         implicit none
         type(vector), intent(in) :: res(:)
         type(chain), pointer :: bchain
         type(matrix) :: dm
         type(matrix), pointer :: c_cm(:)
         type(vector), pointer :: c_res(:)
         real(rfp) :: ar
         integer :: l, nfb, n
         ! initial dqv
         cells % res = res
         do i = 1, size(cells)
            dm = cells(i) % dm
            cells(i) % dqv = cells(i) % res
            call normalize(dm % e, cells(i) % dqv % v)
         end do
         !
         do k = 1, ndim
            nfb = size(raw(k) % bchain)
            do i = 1, nfb
               n = raw(k) % nchain(i)
               if (n <= 1) cycle
               allocate(c_cm(n), c_res(n))
               bchain => raw(k) % bchain(i)
               call forward_sweep(c_cm, c_res, bchain, n)
               deallocate(c_cm, c_res)
            end do
            !
            !    call pgsforlgs(1,1,1) 
            !
         end do
         !    
      end subroutine my_solve_lgs

      subroutine my_gmres_lgs(toler, pdm)
         type(matrix), intent(in) :: pdm(:)
         integer :: ic, min_iter, max_iter, i, i1, j, k1, k, ci, cii, nc, n, m
         real(rfp) :: toler, norm, t, gam
         real(rfp), allocatable :: h(:,:), c(:), s(:), rs(:)
         type(vector), pointer :: v(:,:)
         type(vector) :: dqv0(size(cells))
         type(cell), pointer :: cc
         !
         m = s_isub
         min_iter = 1
         max_iter = m
         nc = size(cells)
         allocate(v(nc, m + 1))
         allocate(h(m + 1, m), c(m), s(m), rs(m + 1))

         do ci = 1, nc
            cc => cells(ci)
            v(ci, 1) = cells(ci) % res - my_matmul(cc, pdm(ci)) ! b-Ax^0
            dqv0(ci) = cells(ci) % dqv ! save dqv^0 in res 
         end do
         ic = 0
         do
            norm = my_norm2(v(:, 1))

            !     ** initialize 1-st term  of rhs of hessenberg system..            
            rs(1) = norm
            i = 0
            do
               !        print *, ic, norm
               if (i == m .or. norm <= toler .or. IC >= max_iter) EXIT
               ic = ic + 1; i = i + 1; i1 = i + 1
               call my_solve_lgs(v(:, i))
               do ci = 1, nc
                  cc => cells(ci)
                  v(ci, i1) = my_matmul(cc, pdm(ci))
               end do
               !-----------------------------------------                              
               !     modified gram - schmidt...                                        
               !-----------------------------------------                              
               do j = 1, i
                  t = my_dot_product(v(:, j), v(:, i1))
                  h(j, i) = t
                  v(:, i1) = v(:, i1) - v(:, j) * t
               end do
               t = sqrt(my_dot_product(v(:, i1), v(:, i1)))
               h(i1, i) = t
               if (t /= zero) v(:, i1) = v(:, i1) * (one / t)
               !                                                                       
               !     done with modified gram schimd and arnoldi step..                 
               !     now  update factorization of h                                   
               !                                                                       
               !--------perfrom previous transformations  on i-th column of h          
               do k = 2, i
                  k1 = k - 1
                  t = h(k1, i)
                  h(k1, i) = c(k1) * t + s(k1) * h(k, i)
                  h(k, i) = -s(k1) * t + c(k1) * h(k, i)
               end do
               gam = sqrt(h(i, i)**2 + h(i1, i)**2)
               !  
               !     if gamma is zero then any small value will do...                  
               !     will affect only residual estimate                                
               !                                                                       
               IF (gam .eq. zero) gam = EPSILON(gam)
               !                                                                        
               !     get  next plane rotation 
               c(i) = h(i, i)/gam
               s(i) = h(i1, i)/gam
               rs(i1) = -s(i) * rs(i)
               rs(i) = c(i) * rs(i)
               !    
               !     determine residual norm and test for convergence- 
               ! 
               h(i, i) = c(i) * h(i, i) + s(i) * h(i1, i)
               norm = ABS(rs(i1))
            end do
            !                                                                       
            !     now compute solution. first solve upper triangular system.        
            !                                                                       
            DO k = i, 1, -1
               rs(k) = (rs(k) - dot_product(h(k, k + 1:i), rs(k + 1:i))) / h(k, k)
            END DO
            !
            ! form linear combination of v(*,i)'s to get solution and call preconditioner
            !
            cells % dqv = v(:, 1) * rs(1)
            do k = 2, i
               cells % dqv = cells % dqv + v(:, k) * rs(k)
            end do

            call my_solve_lgs(cells % dqv)
            dqv0 = dqv0 + cells % dqv
            cells % dqv = dqv0 ! save dqv in res
            !                                                                       
            !     restart outer loop  when necessary                                
            !                                                                       
            IF (norm <= toler .or. IC >= max_iter) EXIT

            !                                                                       
            !     else compute residual vector and continue..                       
            !                                                                       
            DO j = i + 1, 2, -1
               rs(j - 1) = -s(j - 1) * rs(j)
            END DO

            rs(2:i + 1) = c(1:i) * rs(2:i + 1)
            v(:, 1) = v(:, 1) * rs(1)
            do k = 2, i + 1
               v(:, 1) = v(:, 1) + v(:, k) * rs(k)
            end do

         end do
         !      sord%subit = ic
         deallocate(v, h, c, s, rs)
      end subroutine my_gmres_lgs


      subroutine gmres_dilu
         type(matrix) :: pdm(size(cells))
         real(rfp) :: norm, toler
         integer :: k, j

         ! normalize linear eqs
         !
         do k = 1, neq
            norm = maxval(abs(cells % res % v(k)))
            if (norm > 1.0e - 10_rfp) then
               norm = one / norm
               cells % res % v(k) = cells % res % v(k) * norm
               do j = 1, neq
                  cells % dm % e(k, j) = cells % dm % e(k, j) * norm
                  faces % ajl % e(k, j) = faces % ajl % e(k, j) * norm
                  faces % ajr % e(k, j) = faces % ajr % e(k, j) * norm
               end do
            end if
         end do

         toler = sqrt(my_dot_product(cells % res, cells % res))
         !   toler = sqrt(my_dot_product(cells%dqv,cells%dqv))
         toler = toler * 1.0e - 2_rfp
         pdm = cells % dm
         !
         call my_pivot
         !
         call my_solve_lu(cells % res)
         !
         call my_gmres(toler, pdm)
         !
      end subroutine gmres_dilu

      subroutine my_gmres(toler, pdm)
         type(matrix), intent(in) :: pdm(:)
         integer :: ic, min_iter, max_iter, i, i1, j, k1, k, ci, cii, nc, n, m
         real(rfp) :: toler, norm, t, gam
         real(rfp), allocatable :: h(:,:), c(:), s(:), rs(:)
         type(vector), pointer :: v(:,:)
         type(cell), pointer :: cc
         !
         m = s_isub
         min_iter = 1
         max_iter = m
         nc = size(cells)
         allocate(v(nc, m + 1))
         allocate(h(m + 1, m), c(m), s(m), rs(m + 1))

         do ci = 1, nc
            cc => cells(ci)
            v(ci, 1) = cells(ci) % res - my_matmul(cc, pdm(ci)) ! b-Ax^0
            cells(ci) % res = cells(ci) % dqv ! save dqv^0 in res 
         end do
         ic = 0
         do
            norm = my_norm2(v(:, 1))

            !     ** initialize 1-st term  of rhs of hessenberg system..            
            rs(1) = norm
            i = 0
            do
               !        print *, ic, norm
               if (i == m .or. norm <= toler .or. IC >= max_iter) EXIT
               ic = ic + 1; i = i + 1; i1 = i + 1
               call my_solve_lu(v(:, i))
               do ci = 1, nc
                  cc => cells(ci)
                  v(ci, i1) = my_matmul(cc, pdm(ci))
               end do
               !-----------------------------------------                              
               !     modified gram - schmidt...                                        
               !-----------------------------------------                              
               do j = 1, i
                  t = my_dot_product(v(:, j), v(:, i1))
                  h(j, i) = t
                  v(:, i1) = v(:, i1) - v(:, j) * t
               end do
               t = sqrt(my_dot_product(v(:, i1), v(:, i1)))
               h(i1, i) = t
               if (t /= zero) v(:, i1) = v(:, i1) * (one / t)
               !                                                                       
               !     done with modified gram schimd and arnoldi step..                 
               !     now  update factorization of h                                   
               !                                                                       
               !--------perfrom previous transformations  on i-th column of h          
               do k = 2, i
                  k1 = k - 1
                  t = h(k1, i)
                  h(k1, i) = c(k1) * t + s(k1) * h(k, i)
                  h(k, i) = -s(k1) * t + c(k1) * h(k, i)
               end do
               gam = sqrt(h(i, i)**2 + h(i1, i)**2)
               !  
               !     if gamma is zero then any small value will do...                  
               !     will affect only residual estimate                                
               !                                                                       
               IF (gam .eq. zero) gam = EPSILON(gam)
               !                                                                        
               !     get  next plane rotation 
               c(i) = h(i, i)/gam
               s(i) = h(i1, i)/gam
               rs(i1) = -s(i) * rs(i)
               rs(i) = c(i) * rs(i)
               !    
               !     determine residual norm and test for convergence- 
               ! 
               h(i, i) = c(i) * h(i, i) + s(i) * h(i1, i)
               norm = ABS(rs(i1))
            end do
            !                                                                       
            !     now compute solution. first solve upper triangular system.        
            !                                                                       
            DO k = i, 1, -1
               rs(k) = (rs(k) - dot_product(h(k, k + 1:i), rs(k + 1:i))) / h(k, k)
            END DO
            !
            ! form linear combination of v(*,i)'s to get solution and call preconditioner
            !
            cells % dqv = v(:, 1) * rs(1)
            do k = 2, i
               cells % dqv = cells % dqv + v(:, k) * rs(k)
            end do

            call my_solve_lu(cells % dqv)
            cells % res = cells % res + cells % dqv
            cells % dqv = cells % res ! save dqv in res
            !                                                                       
            !     restart outer loop  when necessary                                
            !                                                                       
            IF (norm <= toler .or. IC >= max_iter) EXIT

            !                                                                       
            !     else compute residual vector and continue..                       
            !                                                                       
            DO j = i + 1, 2, -1
               rs(j - 1) = -s(j - 1) * rs(j)
            END DO

            rs(2:i + 1) = c(1:i) * rs(2:i + 1)
            v(:, 1) = v(:, 1) * rs(1)
            do k = 2, i + 1
               v(:, 1) = v(:, 1) + v(:, k) * rs(k)
            end do

         end do
         !      sord%subit = ic
         deallocate(v, h, c, s, rs)
      end subroutine my_gmres

      function my_norm2(qv) result(tnorm)
         type(vector), intent(inout) :: qv(:)
         real(rfp) :: tnorm, norm_inv
         integer :: k
         tnorm = zero
         do k = 1, neq
            tnorm = tnorm + dot_product(qv % v(k), qv % v(k))
         end do
         tnorm = sqrt(tnorm)
         if (tnorm <= zero) return
         norm_inv = one / tnorm
         qv(:) = qv(:) * norm_inv
      end function my_norm2

      function my_dot_product(qv1, qv2) result(dot)
         type(vector), intent(in) :: qv1(:), qv2(:)
         real(rfp) :: dot
         integer :: k
         dot = zero
         do k = 1, neq
            dot = dot + dot_product(qv1 % v(k), qv2 % v(k))
         end do
      end function my_dot_product

      subroutine bicgstab_dilu
         type(matrix) :: pdm(size(cells))
         real(rfp) :: norm, toler
         integer :: k, j

         ! normalize linear eqs
         !
         do k = 1, neq
            norm = maxval(abs(cells % res % v(k)))
            if (norm /= zero) then
               norm = one / norm
               cells % res % v(k) = cells % res % v(k) * norm
               do j = 1, neq
                  cells % dm % e(k, j) = cells % dm % e(k, j) * norm
                  faces % ajl % e(k, j) = faces % ajl % e(k, j) * norm
                  faces % ajr % e(k, j) = faces % ajr % e(k, j) * norm
               end do
            end if
         end do

         toler = sqrt(my_dot_product(cells % res, cells % res))
         toler = toler * 0.01_rfp
         !   toler = 0.001_rfp

         pdm = cells % dm
         !
         call block_dlu ! initial dqv using D-ILU

         call my_bicgstab(toler, pdm)

      end subroutine bicgstab_dilu

      subroutine my_bicgstab(toler, pdm)
         type(matrix), intent(in) :: pdm(:)
         integer :: ic, min_iter, max_iter, i, i1, j, k1, k, ci, nc
         real(rfp) :: toler, norm, gam
         REAL(rfp) :: alpha, beta, rho, omega, temp, c, norm_s, norm_t
         type(vector), dimension(size(pdm)) :: p, r, r_bar, s, t, v, y, z
         type(vector) :: qv
         type(cell), pointer :: cc
         !
         min_iter = 1
         max_iter = s_isub
         nc = size(cells)
         ic = 0
         do ci = 1, nc
            cc => cells(ci)
            r(ci) = cells(ci) % res - my_matmul(cc, pdm(ci))
         end do

         cells % res = cells % dqv ! save dqv^0 in res 

         ic = 0; r_bar = r; rho = one; alpha = one; omega = one;
         do k = 1, neq
            p % v(k) = zero
            v % v(k) = zero
         end do

         DO
            norm = sqrt(my_dot_product(r, r))
            !       Print *, ic, norm
            IF (norm <= toler .or. ic >= max_iter) EXIT
            ic = ic + 1
            temp = my_dot_product(r_bar, r)
            beta = alpha * (temp/(rho + SIGN(EPSILON(1.0_rfp), rho)))/omega
            rho = temp
            p = r + (p - v * omega) * beta
            call my_solve_lu(p)
            y = cells % dqv
            do ci = 1, nc
               cc => cells(ci)
               v(ci) = my_matmul(cc, pdm(ci))
            end do

            alpha = rho / (temp + SIGN(EPSILON(1.0_rfp), temp))
            s = r - v * alpha
            call my_solve_lu(s)
            z = cells % dqv
            do ci = 1, nc
               cc => cells(ci)
               t(ci) = my_matmul(cc, pdm(ci))
            end do

            !  Stabilizing BiCGSTAB by alternative evaluation of omega

            norm_s = sqrt(my_dot_product(s, s))
            norm_t = sqrt(my_dot_product(t, t))
            c = my_dot_product(s, t)/(norm_s * norm_t)
            omega = sign(max(abs(c), 0.7_rfp), c) * norm_s/norm_t
            cells % res = cells % res + y * alpha + z * omega
            r = s - t * omega
         END DO

         cells % dqv = cells % res
      end subroutine my_bicgstab

      function my_matmul(cc, dm)result(xv)
         type(cell), pointer :: cc
         type(matrix), intent(in) :: dm
         type(vector) :: xv
         integer :: j
         xv = dm * cc % dqv
         do j = 1, size(cc % sface) ! loop faces surrounding the cell
            current_face => cc % sface(j) % to_face
            if (current_face % itype /= 0) cycle
            ghost_cell => current_face % left_cell
            if (associated(ghost_cell, cc)) then
               ghost_cell => current_face % right_cell
               dqv => ghost_cell % dqv
               xv = xv + current_face % ajr * dqv
            else
               dqv => ghost_cell % dqv
               xv = xv + current_face % ajl * dqv
            end if
         end do ! end loop faces

      end function my_matmul

      subroutine my_pivot
         !
         do i = 1, size(cells)
            current_cell => cells(i)
            dm => current_cell % dm
            call matrixinv(dm % e)
            cycle
            do j = 1, size(current_cell % sface) ! loop faces surrounding the cell
               current_face => current_cell % sface(j) % to_face
               if (current_face % itype /= 0) cycle
               ghost_cell => current_face % right_cell
               if (associated(ghost_cell, current_cell)) cycle
               ghost_cell % dm = ghost_cell % dm - &
                  current_face % ajr * dm * current_face % ajl
            end do ! end loop faces
         end do
      end subroutine my_pivot


      subroutine my_solve_lu(rhs)
         type(vector) :: rhs(:)
         integer :: i, j
         !
         ! solve (L+D)z = y
         do i = 1, size(cells)
            current_cell => cells(i)
            deltv % v = zero
            do j = 1, size(current_cell % sface) ! loop faces surrounding the cell
               current_face => current_cell % sface(j) % to_face
               if (current_face % itype /= 0) cycle
               ghost_cell => current_face % left_cell
               if (associated(ghost_cell, current_cell)) cycle
               dqv => ghost_cell % dqv
               deltv = deltv + current_face % ajl * dqv
            end do ! end loop faces
            dqv => current_cell % dqv
            dqv = current_cell % dm * (rhs(i) - deltv)
         end do ! end loop cells
         ! solve (I+D^-1U)x = z
         do i = size(cells) - 1, 1, -1
            current_cell => cells(i)
            deltv % v = zero
            do j = 1, size(current_cell % sface) ! loop faces surrounding the cell
               current_face => current_cell % sface(j) % to_face
               if (current_face % itype /= 0) cycle
               ghost_cell => current_face % right_cell
               if (associated(ghost_cell, current_cell)) cycle
               dqv => ghost_cell % dqv
               deltv = deltv + current_face % ajr * dqv
            end do ! end loop faces
            dqv => current_cell % dqv
            dqv = dqv - current_cell % dm * deltv
         end do ! end loop cells

      end subroutine my_solve_lu

   end subroutine linear_solve

end module dfd_linear
