!***********************************************************************
!   This program is for
!   General Flow Calculation Code in hybrid unstructured grid
!   using procondition's method designed in Fortran90
!
!               by Dr. Ding Li
!
!
!This code is designed to solve general equations:
!
! Gamma.DQ/Dt + D.dot.E = D.dot.Ev + S
!where Q=[p,u,v,w,t,k,....]^T  non-conservative variables
!      E is flux vector   Ev is viscous vector  S is source vector 
!
! If your problems can be described in this form governing equations,
! you can use this code to solve your problems
!
!
!1.  March, 1999      Version 0.1
!  only for Eluer Equation with 0-order upwind finite volume schemes
!2.  Sept. 1999       Version 0.2
!  Add high-order upwind finite volume schemes and viscous terms for 
!  Navier-Stokes equations
!3.  Oct. 1999        Version 0.5
!  Directly read structured grid in plot3d format and changed viscous term's treatment
!  Added DILU, GMRES and BICGSTAB linear solve
!4. Oct., 1999   
!  Developed a modified Line Gauss Seidel method for hybrid grid
!  Use new method for reconstructure replace least square, improve stability
!5. Nov., 1999
!  Add k-omega turbulent model
!
!6. July 2002
!  Rewrite the boundary conditions and add B-L turbulence model
!7. Aug. 2002
! Version 1.0 released
!8. Oct. 2005
!  Add Conjugated Heat Transfer
!  Add Radiation model
!  Add DES model
!  Add Maxwell Equations
!  ADD Multi-physics zone mehtods
!  Add job schedule control
!
!------------------makefile--------------------
!
! Makefile_s or Makefile_mpi
!
!----- input file
!
! gems.inp
!
!-----Output file
! gems.res.dat  residual history records
! gems.plt.x datafile for tecplot post-processor
!
!------------------------------------------------
!
!report bugs to dli@purdue.edu or ding_li@yahoo.com
!
!
!**************************************************************

PROGRAM dfd_main

    USE dfd_data
    USE dfd_disc
    USE dfd_linear
    USE dfd_fv
    USE dfd_input
    USE dfd_output
    USE dfd_bound
    USE dfd_react
    USE dfd_source
    USE mpi

    implicit none
    !
    type(face), pointer, save :: faces(:)
    type(cell), pointer, save :: cells(:)
    type(node), pointer, save :: nodes(:)
    type(bnode), pointer, save :: bnodes(:)
    type(itf) :: interf, pinterf
    type(raw_chain) :: raw(ndim)
    !
    real(rfp) :: cpu_time(6), tcpu(6), tb !---   record escaped time in seconds
    logical :: iexist
    integer :: nadv, k, istop, ierr, np, id, ischedule !---   Iterative number
    integer :: lsonly

    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, id, ierr)
    call mpi_comm_size(mpi_comm_world, np, ierr)
    if (id == 0) then
        iexist = .false.
        call unlink('stopfile')
    end if
    !    
    !***  read the input data and allocate neccesary variables
    !
    call read_data(nodes, cells, faces, interf, pinterf, raw, bnodes)
    !
    !*** set initial parameters and initial conditions
    !
    !calculate the pseudo inverse of the least-square-fit problem, save it in cells%svd
    call reconstruction(cells, faces)
    !calculate cells%weight for nodal value calculation
    call inverse_average(cells, nodes, faces, interf) ! pass for multi-domain parallel process

    call cal_node_qv(cells, nodes, faces, interf, bnodes)
    
    !level-set preparation
    if(usels>0) then
      call prepare_levelset(cells,interf)
      !debug
      !call levelset_check(cells,faces,interf)
      !debug
    end if

    nadv = 0
    call wenda_output(cells,faces,np,id,nadv) 
    !call output_result(cells, nodes, faces, interf, id, nadv)
    !call output_to_plot(cells, nodes, faces, interf, id)

    lsonly = 0
    if (laser1 % laser_on < 0.5) then
        lsonly = 1
    end if

    print *, 'lsonly', lsonly
    call mpi_barrier(mpi_comm_world, ierr)
    ! initial cpu_time clock
    tb = mpi_wtime()

    open(79, file='./output/mass.dat')
    do ischedule = 1, schedule % n
        do nadv = s_nstart, s_nend
            !
!test level-set
	  call levelset_advection(cells,interf,nadv)
	  if(mod(nadv,10)==0) call levelset_reinitialization(cells,interf,nadv)
	  if(mod(nadv,50)==0) call wenda_output(cells,faces,np,id,nadv)
	  if(nadv==300) call flip_velocity(cells)
	  if(mod(nadv,10)==0) then
	    write(79,*) nadv,first_order_volume(cells)
	  end if
	  call update_levelset(cells)

            if (id == 0 .and. lsonly == 0) print *, 'ct and time', nadv, s_elapsed_time

            if (lsonly == 0) call decide_cfl(nadv)

            if (lsonly == 0) call set_zero(cells, faces)

            if (lsonly == 0) call grad_recon(cells, faces, nodes)
            if (lsonly == 0) call update_gradient(nodes, cells, interf)

            if (lsonly == 0) call precondition_cell(s_ipre, cells, nodes, nadv)
            if (lsonly == 0) call residual(cells, faces, nodes)

            if (lsonly == 0) call cell_gradient(cells, faces, nodes)
            if (lsonly == 0) call update_gradient(nodes, cells, interf)

            !   call turb_model_1998(cells)
            if (lsonly == 0) call turb_model_2006(cells)
            if (lsonly == 0) call viscous(faces, nodes, cells)
            if (lsonly == 0) call source(cells, nodes)
            if (lsonly == 0) call chemical_reaction(cells)
            if (lsonly == 0) call implicit_boundary_condition(faces)
            if (lsonly == 0) call lhs(cells, faces, nodes)

            if (lsonly == 0) call linear_solve(cells, faces, nodes, interf, raw)
            !
            !*** update variables

            if (lsonly == 0) call update_cells(cells, faces, nadv)
            !
            if (s_iperiod > 0) call update_pinterface(nodes, cells, faces, pinterf)
            !
            if (lsonly == 0)call boundary_condition(faces, nodes, cells, interf)

            if (lsonly == 0)call update_interface(faces, nodes, cells, interf)
            !
            if (s_ivis == 3) call turb_baldwin_lomax(cells, faces, interf)
            !
            if (lsonly == 0)call cal_node_qv(cells, nodes, faces, interf, bnodes)

            !if (mod(nadv, s_imonit * s_nt) == 0) call wenda_output(cells, faces, np, id, nadv, s_nplot)
	    !call output_result(cells,nodes,faces,interf,id,nadv)
            !   if(s_idt /= 0 ) call output_result(cells,nodes,faces,interf,id,nadv)
            !
            !**** Check convergency

            if (lsonly == 0) call check_conv(cells, faces, nadv)

            if (mod(nadv, 100) == 0) call save_binary(nadv, cells, faces, nodes)
            if (lsonly == 0) call update_unsteady_data(cells, faces, interf, nadv)
	    call mpi_barrier(mpi_comm_world, ierr)
        end do
        !
        s_nend = min(s_nend, nadv)
        call set_schedule(ischedule)
    end do
    !*** output elapsed total time in seconds
    !   if(id == 0) print *,' Total wall time =',tcpu(1),(tcpu(k)-tcpu(k-1),k=2,6),' Sec.' 
    tb = mpi_wtime() - tb
    if (id == 0) print *, ' Total CPU time =', tb, 'Sec'
    call report_cpu_usage(id, tb)
    !
    call boundary_condition(faces, nodes, cells, interf)

    call update_interface(faces, nodes, cells, interf)

    call cal_node_qv(cells, nodes, faces, interf, bnodes)

    call save_binary(nadv - 1, cells, faces, nodes)
    call deformation_error(cells)
    close(79)

    !call output_to_plot(cells, nodes, faces, interf, id)

    !call user_define_output(cells, nodes, faces, interf, id)

    call mpi_finalize(ierr)

contains
    !
    subroutine update_cells(cells, faces, nit)
        integer, intent(in) :: nit
        type(face), pointer :: faces(:), sf
        type(cell), pointer :: cells(:), cc
        type(vector) :: qv
        integer :: i, j
        real(rfp) :: c1, cinv, alpha, sumy
        real(rfp) :: vel(ndim), dv(ndim), kmin, ke, vecn(ndim), mpbase
        mpbase = -g_pbase + 1.e-8_rfp
        !
        cinv = one / 0.2_rfp
        do i = 1, size(cells)

            !cells(i)%dqv%v(ils) = zero !levelset is updated in levelset_calculation (in source.f90)

            cells(i) % qv = cells(i) % qv + cells(i) % dqv

            if (nspe > 1) then
                cells(i) % qv % v(isb:ise) = min(max(zero, cells(i) % qv % v(isb:ise)), one)
                sumy = sum(cells(i) % qv % v(isb:ise))
                if (sumy > one) cells(i) % qv % v(isb:ise) = cells(i) % qv % v(isb:ise) * (one / sumy)
            end if
            ! pressure correction
            if (nmeq > 0) then
                call correct_value(cells, i, ico, mpbase)
            end if

            ! temperature correction
            if (neqf > 0) then
                call correct_value(cells, i, ien, zero)
            end if
            !
            if (s_ivis == 2) then
                call correct_value(cells, i, ike, zero)
                call correct_value(cells, i, iom, zero)
            end if
        end do

        !   if(s_irec == 0)  then   ! smooth results
        !    do i = 1, size(cells)
        !      qv%v = zero
        !     do j = 1, size(cells(i)%scell)
        !       cc => cells(i)%scell(j)%to_cell 
        !       qv = qv + cc%qv - cells(i)%qv
        !     end do
        !       qv%v = qv%v / real(size(cells(i)%scell),rfp)
        !       cells(i)%qv = cells(i)%qv + qv
        !     end do
        !   end if
        !
        if (neqm > 0) cells % qv % v(ife) = zero
        if (neqm > 0) cells % dqv % v(ife) = zero
        !
    end subroutine update_cells
    !
    subroutine update_unsteady_data(cells, faces, intf, nit)
        type(face), pointer :: faces(:)
        type(cell), pointer :: cells(:), cc
        type(itf) :: intf
        integer :: i, j, nit

        if (s_idt == 0 .or. mod(nit, s_nt) /= 0) return
        s_elapsed_time = s_elapsed_time + one / s_dt_inv
        if (npos > 0) then
            if (nit/s_nt == 1) then
                open(8, file = 'gems.points.' // i2s(id))
                write(8, *) 'variables=t'
                do i = 1, neq * npos
                    write(8, '(a)') "v" // i2s(i)
                end do
            else
                open(8, file = 'gems.points.' // i2s(id), position = 'append')
            end if
            write(8, *) s_elapsed_time, (cells(s_ipos(j)) % qv % v, j = 1, npos)
            close(8)
        end if
        gear_in = cshift(gear_in, -1)
        do i = 1, size(cells)
            cells(i) % qvn(gear_in(1)) = cells(i) % qv
        end do
        !
        do i = intf % nitf + 1, size(faces)
            if (faces(i) % itype /= 0) exit
            cc => faces(i) % right_cell
            cc % qvn(gear_in(1)) = cc % qv
        end do

        do i = 1, size(intf % pcell)
            cc => intf % pcell(i)
            cc % qvn(gear_in(1)) = cc % qv
        end do
    end subroutine update_unsteady_data


    subroutine correct_value(cells, i, k, c)
        integer :: k, j, i
        type(face), pointer :: sf
        type(cell), pointer :: cells(:), cc
        real(rfp) :: c, kmin, kmax, ke
        if (cells(i) % qv % v(k) < c) then
            kmin = myhuge
            do j = 1, size(cells(i) % sface)
                sf => cells(i) % sface(j) % to_face
                ke = sf % left_cell % qv % v(k)
                if (ke > c) kmin = min(kmin, ke)
                ke = sf % right_cell % qv % v(k)
                if (ke > c) kmin = min(kmin, ke)
            end do
            cells(i) % qv % v(k) = max(kmin, c)
        end if
    end subroutine correct_value

    subroutine decide_cfl(step1)
        integer :: step1, step
        real(rfp) :: alpha

        step = step1
        !  if(s_idt > 0 ) step = mod(step1,s_nt) + 1
        !   s_cfl = s_cflmin  + (abs(s_cflmax) - s_cflmin)* (atan(8.0 *(step/s_damping - 0.8_rfp)) + 1.4) * 0.4
        alpha = min(max(1.e-10_rfp, step/s_damping), one)
        s_cfl = s_cflmin + (abs(s_cflmax) - s_cflmin) * alpha
        !   s_cfl = max(s_cflmin,s_cfl)
        !   s_cfl = min(abs(s_cflmax),s_cfl)
        if (s_cflm < mytiny) then
            s_cflmm = s_cfl
        else
            s_cflmm = alpha * s_cflm
        end if
        s_cfl = sign(s_cfl, s_cflmax)
        s_omega = min(one, real(step, rfp) / s_damping_o) * s_omega0
        !   if(step > two * s_damping) s_ischeme = 0
        s_nstep = step
    end subroutine decide_cfl
    !
    !*************************************************************
    subroutine check_conv(cells, faces, nadv)
        !*************************************************************
        implicit none
        type(cell), pointer :: cells(:)
        type(face), pointer :: faces(:)
        type(cell), pointer :: celli
        !  real(rfp)::errk(neq),sumq(neq)
        real(rfp) :: rate, l2q
        real(rfp), allocatable :: sumqt(:), errt(:), errk(:), sumq(:)
        integer :: nadv, ifile, i, neqt, k, id, neqx
        integer, allocatable :: displs(:), rcount(:)
        !*** Compute nomiralized residuals
        call mpi_comm_rank(mpi_comm_world, id, ierr)
        neqx = maxval(pz % neq)
        call mpi_comm_size(mpi_comm_world, k, ierr)
        !    call mpi_bcast(neqx,1,mpi_integer,0,mpi_comm_world,ierr)
        neqt = neqx * k
        allocate(sumqt(neqt), errt(neqt), displs(k), rcount(k))
        allocate(sumq(neqx), errk(neqx))
        sumq = zero !one    !mytiny 
        errk = zero
        do i = 1, neq
            sumq(i) = sumq(i) + snrm2(cells % qv % v(i))
            errk(i) = snrm2(cells % dqv % v(i))
            if (sumq(i) < 1.e-20_rfp) sumq(i) = one
        end do
        if (nmeq > 0) sumq(imb:ime) = sqrt(sum(sumq(imb:ime)**2))
        if (neqm > 0) sumq(ibb:ibe) = sqrt(sum(sumq(ibb:ibe)**2))
        if (neqm > 0) sumq(ieb:iee) = sqrt(sum(sumq(ieb:iee)**2))
        call mpi_gather(sumq, rfp * neqx, mpi_byte, sumqt, neqx * rfp, &
        mpi_byte, 0, mpi_comm_world, ierr)
        call mpi_gather(errk, rfp * neqx, mpi_byte, errt, neqx * rfp, &
        mpi_byte, 0, mpi_comm_world, ierr)
        s_istop = 0 ! not use err criteria condition
        if (id == 0) then
            !    k = pz(0)%neq + 1
            k = neqx + 1
            do i = 1, size(pz) - 1
                sumq(pz(i) % loc) = sumq(pz(i) % loc) + sumqt(k:k + pz(i) % neq - 1)
                errk(pz(i) % loc) = errk(pz(i) % loc) + errt(k:k + pz(i) % neq - 1)
                !    k = k + pz(i)%neq
                k = k + neqx
            end do
            !
            sumq = errk / sumq + 1.e-20_rfp
            rate = log10(sum(sumq(:neqx))/real(neqx, rfp))
            !    s_alpha = min(max(-rate / 5.0_rfp,zero),one)
            !    s_alpha = min(one,real(nadv,rfp)/1000.0_rfp)
            !
            if (rate > s_errm.or.nadv < 10) s_istop = 0
            if (s_iscr == 1) write (*, 501) nadv, s_cfl, rate
            if (s_iscr == 2) write (*, 501) nadv, s_cfl, rate, errk
            if (s_iscr == 3) write (*, 501) nadv, s_cfl, rate, log10(sumq)
            if (s_ifile == 1) write (7, 501) nadv, rate, log10(sumq)
            if (s_ifile == 2) write (7, 501) nadv, rate, errk
            if (s_ifile == 3) then
                do i = 1, neq
                    errk(i) = snrm2(cells % res % v(i)/cells % vol)
                end do
                errk = errk / real(size(cells), rfp)
                write (7, 501) nadv, rate, log10(errk)
            end if
            if (mod(nadv, 100) == 0) call flush(7)
            deallocate(sumqt, errt, displs, rcount)
        end if ! end if of id
        deallocate(sumq, errk)
        !   if(s_ialg == 0) then
        !   if(rate < -5 .or. nadv > 400) s_ialg = -1
        !   end if
        return
        !
        501 format (1x, i5, 1x, 20e15.6)
    end subroutine check_conv

    subroutine report_cpu_usage(id, cpu)
        integer, intent(in) :: id
        real(rfp) :: cpu, it

        !    open(10,file='~/gems/survey/gems_cputime.dat',position='append')
        !    if(id == 0) then
        !    write(10,*)' ------ survey of CPU usage ------'
        !    write(10,*)' total cell number=',size(cells)
        !    write(10,*)' total node number=',size(nodes)
        !    write(10,*)' total face number=',size(faces)
        !
        !    write(10,*)' ndim=',ndim
        !    write(10,*)' neq =',neq
        !    write(10,*)' meth =',s_imeth
        !    write(10,*)' ialg =',s_ialg
        !    write(10,*)' isub =',s_isub
        !    write(10,*)' irea =', s_irealgas
        !    write(10,*)' isource =',s_isource
        !    write(10,*)' total iteration =',s_nend - s_nstart + 1
        !    end if
        !!  
        !    write(10,*)' *** Processor #',id,'*** time in second'
        !    write(10,*)' Total_CPU_times =',cpu
        it = s_nend - s_nstart + 1
        !    write(10,*)' CPUtime per iteration=',cpu/real(s_nend-s_nstart+1,rfp)
        !    write(10,*)' CPUtime per iteration per cell=',cpu/real(s_nend-s_nstart+1,rfp) / &
        !                 real(size(cells),rfp)

        !   write(10,*)size(cells),neq,s_imeth,s_ialg,it,cpu,cpu/it,cpu/it/real(size(cells),rfp)
        !    close(10)
    end subroutine report_cpu_usage
    !
    !*************************************************************
    subroutine release_var(cells)
        !*************************************************************
        implicit none
        integer :: istatus, i, j
        type(cell), pointer :: cells(:)
        !
        do i = 1, size(cells)
            do j = 1, size(cells(i) % sface)
                deallocate(cells(i) % sface(j) % to_face)
            end do
            if (s_idt > 0) deallocate(cells(i) % qvn)
            deallocate(cells(i) % c2n, cells(i) % sface, cells(i) % svd)
        end do
        deallocate(cells)
    end subroutine release_var

end program dfd_main
