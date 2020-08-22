module gems_input
   use gems_data
   use gems_geom
   use gems_bound
   use gems_turb
   use gems_react
   USE RPTBDB_setup, ONLY: ReadDatabase
   use mpi
   implicit none
   character(len = 80) :: gridfile, bfile, gridfile0
   character(len = 80) :: fnames(0:500)
contains

   !
   !********************************************************
   !*** read input data and allocate neccesary variables
   !********************************************************
   subroutine read_data(nodes, cells, faces, interf, pinterf, raw, bnodes)
      type(cell), pointer :: cells(:)
      type(face), pointer :: faces(:)
      type(node), pointer :: nodes(:)
      type(bnode), pointer :: bnodes(:)
      type(itf) :: interf, pinterf
      type(raw_chain) :: raw(ndim)
      logical :: yes
      character(len = 40) :: fname
      integer :: id, ierr, np, id1, id2, itype

      call mpi_comm_rank(mpi_comm_world, id, ierr)
      call mpi_comm_size(mpi_comm_world, np, ierr)
      if (np == 1) then ! single processor
         !   open (unit=11, file='gems.inp',status = 'old')
         fnames(0) = 'gems.inp'
      else
         inquire(file = 'gems.par', exist = yes)
         if (.not.yes) then
            print *, 'Sorry, could not find the gems.par file for multi-zone problems'
            call mpi_abort(mpi_comm_world, ierr, ierr)
            stop
         end if
         open (unit = 12, file = 'gems.par', status = 'old')
         read(12, *)
         do
            read(12, *) id1, id2, itype, fname
            fnames(id1:id2) = fname
            if (id2 == np - 1) exit
         end do
         close(12)
      end if
      inquire(file = fnames(id), exist = yes)
      if (.not.yes) then
         print *, 'Sorry, could not find the input file', fnames(id)
         print *, 'This file is need for partition', id
         call mpi_abort(mpi_comm_world, ierr, ierr)
         stop
      end if

      !read from 11 which is gems.inp
      !write whatever is read from 11 to 32 which is DL_gems.inp
      open (unit = 11, file = fnames(id), status = 'old')
      !open (unit = 32, file = 'DL_' // fnames(id))

      !set up physics, i.e. number of equations for each PE
      !store in pz which is of the type phy
      call setup_problem

      !read the namelist geom from gems.inp
      call grid_information(nodes, cells, faces)

      !read the namelist solution from gems.inp
      call solution_information
      call read_schedule

      !read the namelist unsteady from gems.inp
      call unsteady_parameter

      !set up fluid and solid properties
      call fluid_property

      if(s_irealgas > 110) call solid_property
      !read the namelist turb from gems.inp
      call turbulence_model_paramter
      !
      call read_grid(gridfile, cells, faces, nodes, interf, pinterf)
      print *, 'number of cells=', size(cells), id
      print *, 'number of faces=', size(faces), id
      print *, 'number of nodes=', size(nodes), id

      ! calculate geometry information
      call geometry(nodes, cells, faces) 
      if (s_imeth == 5) call setup_chain(cells, faces, raw)
      if (s_imeth == 7) call setup_chain(cells, faces, raw)

      !read the namelist on boundary condition from gems.inp
      call read_boundary(cells, faces, nodes, interf)

      !creat the right cells from boundary faces
      call set_boundary(faces, nodes, cells, interf)

      ! set debug monitor
      call set_monitor(cells, id)

      call initial(cells, faces, nodes, bnodes,interf)
      call output_option
      call read_external_bfield(cells)
      if (s_ibfield == 2) call output_cell_tecplot_v('bfiled', cells, nodes, faces, id, e_b0)
      ! just for MHD channel
      !
      if (m_input == 1) then
         !  if(neqm > 0.and.vc(cells(1)%itype)%vars(1) < zero) then
         open(10, file = 'sigma.dbs.' // i2s(id))
         read(10, *) cells % er
         do ierr = 1, size(faces)
            if (faces(ierr) % itype == 0) exit
            read(10, *) faces(ierr) % right_cell % er
         end do
         close(10)
      end if
      !for periodic boundary condition
      call setup_periodic_face(faces, cells, pinterf)
      call update_geom_pinterface(faces, nodes, cells, pinterf)
      call update_pinterface(nodes, cells, faces, pinterf)

      !enforce boundary condition
      if (s_nstart <= 1) call boundary_condition(faces, nodes, cells, interf)

      call update_geom_interface(faces, nodes, cells, interf, bnodes)
      call update_interface(faces, nodes, cells, interf)
      !plot wall
      call woutput_option(faces)

      !read debug info
      call read_debug()
      if (id == 0) print *, ' reading input data has been finished'
      if (s_irec /= 2) & ! for least-square reconstruction
         call search_neighbor(nodes, cells, faces)
      if (s_ivis == 3) then
         call init_baldwin_lomax(cells, faces, interf)
      end if

      ! add axisymmetric treatment
      !
      if (ndim == 2 .and. s_iaxis > 0) then
         faces % area = faces % area * abs(faces % centp(s_iaxis))
         cells % vol = cells % vol * abs(cells % centp(s_iaxis))
      end if
      !
      ! 
      !  if(s_ifmt /= 1) call face_geom(nodes,cells,faces)  

      !calcualte face%avec which stores the 3 normal vectors of a face 
      call face_reconstruction(faces, nodes) 

      close(11)
      close(32)

   end subroutine read_data
   !
   subroutine setup_problem
      integer :: id, ierr, np, imeth, j, i
      integer, pointer :: neq_np(:), imeth_np(:), ien_np(:)
      !  0:N+S
      !  1:N+S+Turbulent Equations
      !  2:Heat Transfer
      !  3:MHD
      !  4:Plasma flow
      !
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      call mpi_comm_size(mpi_comm_world, np, ierr)
      np = np - 1
      allocate(neq_np(0:np), imeth_np(0:np), ien_np(0:np))
      allocate(pz(0:np))

      !see module gems_constant for constants
      print *, '-------', id, neq

      imeth = 0 ! N+S Eq.
      if (nteq > 0) imeth = 1 ! N+S+Turbulent Eq.
      if (neq == 1) imeth = 2 ! Heat Transfer Eq.
      if (neqf == 0.and.neqm > 0) imeth = 3 ! Pure Maxwell
      if (neqf > 0 .and.neqm > 0) imeth = 4 ! Plasma flow
      call mpi_allgather(imeth, 1, MPI_integer, imeth_np, 1, MPI_integer, mpi_comm_world, ierr)
      call mpi_allgather(neq, 1, MPI_integer, neq_np, 1, MPI_integer, mpi_comm_world, ierr)
      call mpi_allgather(ien, 1, MPI_integer, ien_np, 1, MPI_integer, mpi_comm_world, ierr)
      if (id == 0) print *, 'number of equations at partitions =', neq_np
      if (id == 0) print *, 'Equation type at partitions =', imeth_np
      pz % imeth = imeth_np
      select case(imeth)
      case(0) ! Inviscid or Laminar flow
         do j = 0, np
            select case(imeth_np(j))
            case(0, 1, 4)
               pz(j) % neq = neq
               allocate(pz(j) % loc(neq))
               pz(j) % loc = (/(i, i = 1, neq)/)
            case(2)
               pz(j) % neq = 1
               allocate(pz(j) % loc(1))
               pz(j) % loc = ien
            case(3)
               pz(j) % neq = 0 ! nothing need to exchange
            case default
               print *, 'not apply', imeth, imeth_np(j)
            end select
         end do
      case(1) ! Turburlent flow
         do j = 0, np
            select case(imeth_np(j))
            case(0)
               pz(j) % neq = neq - nteq
               allocate(pz(j) % loc(pz(j) % neq))
               pz(j) % loc(:ien) = (/(i, i = 1, ien)/)
               if (nspe > 1) pz(j) % loc(ien + 1:) = (/(i, i = isb, ise)/)
            case(1, 4)
               pz(j) % neq = neq
               allocate(pz(j) % loc(neq))
               pz(j) % loc = (/(i, i = 1, neq)/)
            case(2)
               pz(j) % neq = 1
               allocate(pz(j) % loc(1))
               pz(j) % loc = ien
            case(3)
               pz(j) % neq = 0 ! nothing need to exchange
            case default
               print *, 'not apply', imeth, imeth_np(j)
            end select
         end do
      case(2) !   Heat transfer
         do j = 0, np
            select case(imeth_np(j))
            case(0, 1, 2)
               pz(j) % neq = 1
               allocate(pz(j) % loc(1))
               pz(j) % loc = ien
            case default
               print *, 'not apply', imeth, imeth_np(j)
            end select
         end do
      case(3) ! Pure Maxwell  neqf ==0, neqm > 0
         do j = 0, np
            select case(imeth_np(j))
            case(0, 1)
               pz(j) % neq = 0 ! nothing need to exchange
            case(3)
               pz(j) % neq = neq
               allocate(pz(j) % loc(neq))
               pz(j) % loc(:) = (/(i, i = 1, neq)/)
            case(4)
               pz(j) % neq = neqm
               allocate(pz(j) % loc(neqm))
               pz(j) % loc(:) = (/(i, i = ibb, ife)/)
            case default
               print *, 'not apply', imeth, imeth_np(j)
            end select
         end do
      case(4) ! plasma flow neqf >0, neqm > 0
         do j = 0, np
            select case(imeth_np(j))
            case(0, 1)
               pz(j) % neq = neqf
               allocate(pz(j) % loc(neqf))
               pz(j) % loc(:) = (/(i, i = 1, neqf)/)
            case(4)
               pz(j) % neq = neq
               allocate(pz(j) % loc(neq))
               pz(j) % loc(:) = (/(i, i = 1, neq)/)
            case(3)
               pz(j) % neq = neqm
               allocate(pz(j) % loc(neqm))
               pz(j) % loc(:) = (/(i, i = ibb, ife)/)
            case default
               print *, 'not apply', imeth, imeth_np(j)
            end select
         end do
      end select
      deallocate(neq_np, imeth_np, ien_np)
      ! if(id == 0) print *,(pz(i)%neq,i=0,np)
      ! if(id == 0) print *,(pz(i)%loc,i=0,np)
   end subroutine setup_problem

   subroutine read_external_bfield(cells)
      type(cell), pointer :: cells(:)
      integer :: i, j, k, l, ni, nj, nk
      real(rfp) :: p(3), fp, ft, fc
      real(rfp), allocatable :: x(:), y(:), z(:), f(:,:,:)
      if (s_ibfield /= 2) return
      allocate(e_b0(size(cells)))
      open(11, file = 'external_bfield.dat')
      read(11, *) ni, nj, nk
      allocate(x(ni), y(nj), z(nk), f(ni, nj, nk))
      read(11, *) x, y, z
      read(11, *) f
      close(11)
      x = x * b_lenref
      y = y * b_lenref
      z = z * b_lenref
      do l = 1, size(cells)
         ! bi-linear interpolation
         p = avto3(cells(l) % centp)
         i = location(p(1), x(1), x(2) - x(1), ni, fp)
         j = location(p(2), y(1), y(2) - y(1), nj, ft)
         k = location(p(3), z(1), z(2) - z(1), nk, fc)
         e_b0(l) = &
            f(i, j, k) * (one - fp) * (one - ft) * (one - fc) + &
            f(i + 1, j, k) * fp * (one - ft) * (one - fc) + &
            f(i, j + 1, k) * (one - fp) * ft * (one - fc) + &
            f(i, j, k + 1) * (one - fp) * (one - ft) * fc + &
            f(i + 1, j + 1, k) * fp * (ft) * (one - fc) + &
            f(i, j + 1, k + 1) * (one - fp) * (ft) * (fc) + &
            f(i + 1, j, k + 1) * (fp) * (one - ft) * (fc) + &
            f(i + 1, j + 1, k + 1) * (fp) * (ft) * (fc)
      end do
      deallocate(x, y, z, f)
   end subroutine read_external_bfield

   subroutine interpolate_from_coarsegrid(cells, faces, nodes, bnodes)
      type(cell), pointer :: cells(:), cells_c(:), cc
      type(face), pointer :: faces(:), faces_c(:)
      type(node), pointer :: nodes(:), nodes_c(:)
      type(bnode), pointer :: bnodes(:)
      type(itf) :: interf, pinterf
      logical :: inside
      integer :: np1, np2, n_o, n_d, ivar_o(2 * neq), ivar_d(neq)
      character(len = 80) :: gridfile, binfile_location
      integer :: i, j, id, ierr
      namelist /interpolation/gridfile, binfile_location, np1, np2, n_o, n_d, ivar_o, ivar_d
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      open(12, file = fnames(id))
      n_o = neq; n_d = neq;
      np1 = 0; np2 = 0;
      ivar_o(:neq) = (/(i, i = 1, neq)/)
      ivar_d(:neq) = (/(i, i = 1, neq)/)
      binfile_location = ''
      read(12, interpolation)
      !xxxxxxxxxxxxxx
      !write(32, interpolation)


      close(12)
      !  j = 1;last=.false.
      !  do while(.not.last)
      do j = np1, np2
         call read_grid_interp(gridfile, cells_c, faces_c, nodes_c, interf, pinterf, j + 1)
         !  call read_grid(gridfile,cells_c,faces_c,nodes_c,interf,pinterf)
         call geometry(nodes_c, cells_c, faces_c) ! cal geometry information
         call set_boundary(faces_c, nodes_c, cells_c, interf)
         call weighted_average(cells_c, nodes_c, faces_c, interf)
         !  call read_binary(cells_c,faces_c,j-1,neq)
         call read_binary_interp(binfile_location, cells_c, faces_c, j, n_o, ivar_o(:n_d), ivar_d(:n_d))
         call cal_node_qv(cells_c, nodes_c, faces_c, interf, bnodes)
         !  do i = 1,size(cells_c)
         !  cells_c(i)%itype = i
         !  end do
         do i = 1, size(cells)
            call search_location_at_cell(cells(i) % centp, cells_c, faces_c, nodes_c, cc, inside)
            if (inside)call interpolateqv(cells(i) % qv, cc, nodes_c, cells(i) % centp)
         end do
         call release(cells_c, nodes_c, faces_c)
         !  j = j + 1
      end do

   end subroutine interpolate_from_coarsegrid
   !
   subroutine interpolateqv(qv, cc, nodes, p)
      type(cell), pointer :: cc
      type(node), pointer :: nodes(:)
      type(vector) :: qv
      real(rfp) :: p(ndim), r, rt
      integer :: i, n
      n = size(cc % c2n)
      r = sqrt(sum((p - cc % centp)**2))
      if (r < mytiny) then
         qv = cc % qv
         return
      end if
      r = one / r
      qv = cc % qv * r
      rt = r
      do i = 1, n
         r = sqrt(sum((p - nodes(cc % c2n(i)) % xyz)**2))
         if (r < mytiny) then
            qv = nodes(cc % c2n(i)) % qv
            return
         end if
         r = one / r
         qv = qv + nodes(cc % c2n(i)) % qv * r
         rt = rt + r
      end do
      qv = qv * (one/ rt)
   end subroutine interpolateqv

   subroutine release(cells, nodes, faces)
      type(cell), pointer :: cells(:)
      type(face), pointer :: faces(:)
      type(node), pointer :: nodes(:)
      integer :: i, j
      deallocate(nodes)
      do i = 1, size(cells)
         do j = 1, size(cells(i) % sface)
            nullify(cells(i) % sface(j) % to_face)
         end do
         do j = 1, size(cells(i) % scell)
            nullify(cells(i) % scell(j) % to_cell)
         end do
         if (size(cells(i) % sface) > 0) deallocate(cells(i) % sface)
         if (size(cells(i) % scell) > 0) deallocate(cells(i) % scell)
         if (size(cells(i) % svd) > 0) deallocate(cells(i) % svd)
         if (size(cells(i) % weight) > 0) deallocate(cells(i) % weight)
         if (size(cells(i) % gradient) > 0) deallocate(cells(i) % gradient)
         if (size(cells(i) % c2n) > 0) deallocate(cells(i) % c2n)
         if (size(cells(i) % qvn) > 0) deallocate(cells(i) % qvn)
         if (size(cells(i) % jv) > 0) deallocate(cells(i) % jv)
      end do
      deallocate(cells)
      do i = 1, size(faces)
         nullify(faces(i) % left_cell)
         nullify(faces(i) % right_cell)
         if (size(faces(i) % f2n) > 0) deallocate(faces(i) % f2n)
         if (size(faces(i) % avec) > 0) deallocate(faces(i) % avec)
      end do
      deallocate(faces)

   end subroutine release


   subroutine search_location_at_cell(p, cells, faces, nodes, cc, inside)
      implicit none
      type(cell), pointer :: cells(:), cc
      type(face), pointer :: faces(:), cf
      type(node), pointer :: nodes(:)
      real(rfp) :: p(2), dn, dmin, p1(2), p2(2), x
      logical :: inside
      !
      integer :: i, imin
      ! only avariable for 2D now!!!!!!!!!
      dmin = myhuge
      inside = .false.
      do i = 1, size(faces)
         if (faces(i) % itype == 0) exit
         p1 = nodes(faces(i) % f2n(1)) % xyz(:2)
         p2 = nodes(faces(i) % f2n(2)) % xyz(:2)
         if (p(2) <= min(p1(2), p2(2))) cycle
         if (p(2) > max(p1(2), p2(2))) cycle
         if (abs(p2(2) - p1(2)) < mytiny) then
            if (p(1) <= min(p1(1), p2(1))) cycle
            if (p(1) > max(p1(1), p2(1))) cycle
            inside = .true.
            imin = i
            exit
         else
            x = p1(1) + (p2(1) - p1(1)) * (p(2) - p1(2)) / (p2(2) - p1(2))
            if (x < p(1)) cycle
            if (abs(x - p(1)) < mytiny) then
               inside = .true.
               imin = i
               exit
            end if
            inside = .not.inside
            dn = sum((p - faces(i) % centp(:2))**2)
            if (dmin > dn) then
               imin = i
               dmin = dn
            end if
         end if
      end do
      if (inside) then
         cf => faces(imin)
         call check_at_cell(cf, faces(imin) % left_cell, nodes, p, cc) !cc,cells(i)%centp)
      end if
      !
   end subroutine search_location_at_cell
   !
   function insection(p, cf, nodes)result(in)
      implicit none
      type(face), pointer :: cf
      type(node), pointer :: nodes(:)
      real(rfp) :: p(2), dn, dmin, p1(2), p2(2), x
      logical :: in
      !
      in = .false.
      p1 = nodes(cf % f2n(1)) % xyz(:2)
      p2 = nodes(cf % f2n(2)) % xyz(:2)
      if (p(2) <= min(p1(2), p2(2))) return
      if (p(2) > max(p1(2), p2(2))) return
      if (abs(p2(2) - p1(2)) < mytiny) then
         if (p(1) <= min(p1(1), p2(1))) return
         if (p(1) > max(p1(1), p2(1))) return
         in = .true.
         return
      else
         x = p1(1) + (p2(1) - p1(1)) * (p(2) - p1(2)) / (p2(2) - p1(2))
         if (x < p(1)) return
         in = .true.
      end if
   end function insection

   recursive subroutine check_at_cell(cfr, cc, nodes, p, ch)
      type(cell), pointer :: cc, cg, ch
      type(face), pointer :: cf, cfr
      type(node), pointer :: nodes(:)
      real(rfp) :: vecn(ndim), dn, p(2)
      integer :: i, istop
      do i = 1, size(cc % sface)
         cf => cc % sface(i) % to_face
         if (associated(cfr, cf)) cycle
         cg => cf % right_cell
         if (associated(cc, cg)) cg => cf % left_cell
         if (insection(p, cf, nodes)) then
            call check_at_cell(cf, cg, nodes, p, ch)
            return
         end if
      end do
      ch => cc
   end subroutine check_at_cell
   !
   subroutine grid_information(nodes, cells, faces)
      implicit none
      type(cell), pointer :: cells(:)
      type(face), pointer :: faces(:)
      type(node), pointer :: nodes(:)
      !
      integer :: ifmt, iaxis
      real(rfp) :: lenref, nblade

      real(rfp) :: rbuff(2)
      integer :: ibuff(2)
      integer :: i, j, nd
      namelist /geom/gridfile, bfile, ifmt, iaxis, lenref, nd, gridfile0, nblade

      equivalence(ibuff(1), ifmt), (ibuff(2), iaxis)
      equivalence(rbuff(1), lenref), (rbuff(2), nblade)
      integer :: id, ierr
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      ibuff = 0; rbuff = zero
      !
      !  if(id == 0) then
      iaxis = 0; ifmt = 2
      nd = 2
      nblade = one

      !read from gems.inp
      !geom is a namelist
      read (11, geom)
      !xxxxxxxxxxxxxx
      !write(32, geom)
      if (nd /= ndim) then
         print *, 'please change the ndim in als_constant.f90 to ', nd
         print *, 'Re-compile before you run the case'
         call mpi_abort(mpi_comm_world, ierr, ierr)
      end if
      !  end if
      !
      !  call mpi_bcast(rbuff,size(rbuff)*rfp,mpi_byte,0,mpi_comm_world,ierr)
      !  call mpi_bcast(ibuff,size(ibuff),mpi_integer,0,mpi_comm_world,ierr)
      !  call mpi_bcast(gridfile0,80,mpi_character,0,mpi_comm_world,ierr)

      s_ifmt = ifmt
      s_iaxis = iaxis
      b_lenref = lenref
      b_nblade = nblade
      !
      !
   end subroutine grid_information

   subroutine solution_information
      implicit none
      !
      ! Solution control
      !
      integer :: irestart, nsteps, init, ialg, iscr, imeth, irec, iplot, nplot, nsave
      integer :: isub, ivis, irea, isource, iperiodic, iaux, ibfield
      real :: errm, sor, dt, cflmin, cflmax, cflm, vnn, pbase, cspeed, damping
      real :: bspeed, espeed, fspeed, b0(3), rpm, damping_o, bodyvel(ndim)
      integer :: istatus, icolor, ifile, ipre, ischeme, ischedule, minput
      integer :: i, j
      !
      real(rfp) :: rbuff(20), vpre_uns
      integer :: ibuff(30)
      !**** namelis
      !
      namelist /solution/irestart, nsteps, init, &
         ialg, irec, ischeme, ischedule, &
         imeth, isub, ibfield, &
         iscr, ifile, bodyvel, &
         cflmin, cflmax, vnn, damping, rpm, damping_o, cflm, &
         errm, sor, cspeed, bspeed, espeed, fspeed, &
         ipre, pbase, b0, nplot, &
         ivis, irea, iplot, isource, iperiodic, iaux, minput, vpre_uns, nsave
      equivalence(ibuff(1), irestart), (ibuff(2), nsteps), (ibuff(3), init)
      equivalence(ibuff(4), ialg), (ibuff(5), irec), (ibuff(6), imeth)
      equivalence(ibuff(7), iscr), (ibuff(8), ifile), (ibuff(9), isub)
      equivalence(ibuff(10), ipre), (ibuff(11), ivis), (ibuff(12), irea)
      equivalence(ibuff(13), iplot), (ibuff(14), isource), (ibuff(15), iperiodic)
      equivalence(ibuff(16), iaux), (ibuff(17), ischeme), (ibuff(18), ischedule)
      equivalence(ibuff(19), nplot), (ibuff(20), minput), (ibuff(21), ibfield)
      equivalence(ibuff(22), nsave)
      !
      ! for real
      equivalence(rbuff(1), cflmin), (rbuff(2), cflmax), (rbuff(3), vnn)
      equivalence(rbuff(4), errm), (rbuff(5), sor), (rbuff(6), dt)
      equivalence(rbuff(7), pbase), (rbuff(8), cspeed), (rbuff(9), damping)
      equivalence(rbuff(10), bspeed), (rbuff(11), espeed), (rbuff(12), fspeed)
      equivalence(rbuff(13), rpm), (rbuff(14), damping_o)
      equivalence(rbuff(15), cflm), (rbuff(16), bodyvel(1)), (rbuff(17), vpre_uns)
      integer :: id, ierr
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      !  defaults
      !  if(id == 0) then
      ibuff = 0
      rbuff = zero
      pbase = g_patm !1atm = 101.325 kPa
      iscr = 1; ialg = 0; ischeme = 0; damping = 100.; damping_o = 1.
      ischedule = 0; ibfield = 0
      iperiodic = 0; b0 = 0.0
      nplot = 100
      cspeed = s_cspeed
      rpm = 0.0_rfp; bodyvel = 0.0_rfp
      nsave = 4000000
      fspeed = zero
      read(11, solution)
      !xxxxxxxxxxxxxx
      !write(32, solution)
      !  end if
      !
      !  call mpi_bcast(rbuff,size(rbuff)*rfp,mpi_byte,0,mpi_comm_world,ierr)
      !  call mpi_bcast(ibuff,size(ibuff),mpi_integer,0,mpi_comm_world,ierr)
      ! use first input file to control iterations
      call mpi_bcast(ibuff, 2, mpi_integer, 0, mpi_comm_world, ierr)
      !
      s_nstart = irestart
      s_nend = nsteps
      s_init = init
      !
      s_ialg = ialg
      s_irec = irec
      s_ischeme = ischeme
      s_ischedule = ischedule
      !
      s_imeth = imeth
      s_isub = isub
      !
      s_ipre = ipre; s_vpre_uns = vpre_uns
      s_ivis = ivis
      s_irealgas = irea
      s_isource = isource
      s_iperiod = iperiodic ! <3 for cylindar =4 for planar
      s_iaux = iaux
      !
      s_iscr = iscr
      s_ifile = ifile
      s_iplot = iplot
      s_nplot = nplot; s_nsave = nsave
      !
      s_vnn = vnn
      s_cfl = cflmin
      s_cflmin = cflmin
      s_cflmax = cflmax
      s_cflm = cflm
      s_cspeed = cspeed
      s_damping = damping
      s_damping_o = damping_o
      s_errm = errm
      s_bspeed = bspeed
      s_espeed = espeed
      s_fspeed = fspeed
      s_omega0 = rpm * two * pi / 60.0_rfp ! convert to 1/s
      s_body_vel = bodyvel

      g_pbase = pbase
      s_ibfield = ibfield
      g_b0 = b0
      !
      m_input = minput ! 1=read conductivity from fluid solution

      if (ialg == 0) ialg = 1
      s_nrec = ((abs(ialg) - 1) * ndim * ndim + (abs(ialg) + 1) * ndim)/2
      ! 
      call mpi_bcast(s_ischedule, 1, mpi_integer, 0, mpi_comm_world, ierr)

      if (s_ivis == 2.and.id == 0) then
         print *, 'You choose the K-Omega turbulence model'
         if (nteq /= 2) then
            print *, 'please set nteq = 2 in constant.f90'
            print *, 'Re-compile before you run the case'
            call mpi_abort(mpi_comm_world, ierr, ierr)
         end if
      end if


      if (s_ivis == 3.and.id == 0) then
         print *, 'You choose the Baldwin-Lomax turbulence model'
         if (nteq > 0) then
            print *, 'You do not need to set turbulence equation'
            print *, 'Please set nteq = 0'
            call mpi_abort(mpi_comm_world, ierr, ierr)
            stop
         end if
      end if
      !
   end subroutine solution_information
   !
   subroutine unsteady_parameter
      implicit none
      !
      integer :: idt, nt, imonit, ipos(20), i
      real(rfp) :: rbuff(2)
      integer :: ibuff(4)
      real(rfp) :: dt, elapsedtime
      equivalence(ibuff(1), idt), (ibuff(2), nt)
      equivalence(ibuff(3), imonit) !,(ibuff(4),ipos)
      equivalence(rbuff(1), dt), (rbuff(2), elapsedtime)

      namelist /unsteady/idt, nt, imonit, dt, ipos, elapsedtime, npos
      integer :: id, ierr
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      !
      ! if(id == 0) then
      ibuff = 0; rbuff = zero;
      ipos = 0; npos = 0
      elapsedtime = 0
      read(11, unsteady)
      !xxxxxxxxxxxxxx
      !write(32, unsteady)
      ! end if
      ! Bcast dt
      call mpi_bcast(dt, 1 * rfp, mpi_byte, 0, mpi_comm_world, ierr)
      ! call mpi_bcast(rbuff,size(rbuff)*rfp,mpi_byte,0,mpi_comm_world,ierr)
      ! call mpi_bcast(ibuff,size(ibuff),mpi_integer,0,mpi_comm_world,ierr)

      if (idt > 4) idt = 4 ! max is fourth-order
      s_idt = idt
      s_nt = nt
      s_imonit = imonit
      s_dt_inv = one/dt
      s_dt = dt
      if (npos /= 0) then
         allocate(s_ipos(npos))
         s_ipos = ipos(:npos)
      end if
      s_elapsed_time = elapsedtime

      if (idt > 0) then
         !   s_damping = s_nt
         allocate(gear_c(idt + 1), gear_in(idt))
         select case(idt)
         case(1)
            gear_c(1) = one !1
            gear_c(2) = -one !1
         case(2)
            gear_c(1) = three_second !3/2
            gear_c(2) = -two !4/2
            gear_c(3) = half !1/2
         case(3)
            gear_c(1) = 11._rfp/ 6._rfp !11/6
            gear_c(2) = -18._rfp/ 6._rfp !18/6
            gear_c(3) = 9._rfp/ 6._rfp !19/6
            gear_c(4) = -2._rfp/ 6._rfp ! 2/6
         case(4)
            gear_c(1) = 25._rfp/ 12._rfp !25/12
            gear_c(2) = -48._rfp/ 12._rfp !48/12
            gear_c(3) = 36._rfp/ 12._rfp !36/12
            gear_c(4) = -16._rfp/ 12._rfp !16/12
            gear_c(5) = 3._rfp/ 12._rfp ! 3/12   
         end select
         gear_c = gear_c * s_dt_inv
         do i = 1, idt
            gear_in(i) = i
         end do
      end if
      !
   end subroutine unsteady_parameter

   subroutine fluid_property
      implicit none
      ! fluid property
      real :: mwi(nspe), cpi(nspe), zi(nspe), rho0(nspe), hrefi(nspe), htref(nspe)
      real :: zmu(nspe), tref(nspe), sref(nspe), pr(nspe), gravity(ndim), absorptivity
      real :: sigma(nspe), ne(nspe), cdest, cprod
      character*80 :: pfile(nspe), species_name(nspe)
      !
      integer :: i, j, n, ifmt, itype(nspe)
      real(rfp) :: rbuff(12 * nspe + ndim + 1)
      namelist /gas_prop/pfile, mwi, species_name, sigma, ne, ifmt, itype, cdest, cprod, &
         cpi, zi, rho0, hrefi, htref, zmu, tref, sref, pr, gravity, absorptivity
      ! for real
      equivalence(rbuff(1), zmu)
      equivalence(rbuff(nspe + 1), htref)
      equivalence(rbuff(2 * nspe + 1), tref)
      equivalence(rbuff(3 * nspe + 1), sref)
      equivalence(rbuff(4 * nspe + 1), pr)
      equivalence(rbuff(5 * nspe + 1), mwi), (rbuff(6 * nspe + 1), cpi)
      equivalence(rbuff(7 * nspe + 1), hrefi), (rbuff(8 * nspe + 1), zi)
      equivalence(rbuff(9 * nspe + 1), rho0)
      equivalence(rbuff(10 * nspe + 1), sigma)
      equivalence(rbuff(11 * nspe + 1), ne)
      equivalence(rbuff(12 * nspe + 1), gravity)
      equivalence(rbuff(12 * nspe + ndim + 1), absorptivity)
      !
      integer :: id, ierr
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      !
      !  if(id == 0) then
      rbuff = zero;
      pfile = 'gems.property'
      ifmt = 0 ! format of the file
      !
      gravity = zero
      mwi = 28.9_rfp
      itype = 0 ! perfect gas
      zi = one
      rho0 = zero
      zmu = 1.e-5_rfp
      pr = 0.72_rfp
      tref = 273._rfp
      sref = 110.5_rfp
      sigma = one; ne = 1.e30_rfp;
      cdest = one; cprod = one;
      do i = 1, nspe
         species_name(i) = ''
         write(species_name(i), '("Gas",a)') trim(i2s(i))
      end do
      read (11, gas_prop)
      !xxxxxxxxxxxxxx
      !write(32, gas_prop)
      !  end if

      !   call mpi_bcast(rbuff,size(rbuff)*rfp,mpi_byte,0,mpi_comm_world,ierr)
      !   call mpi_bcast(species_name,size(species_name)*80,mpi_character,0,mpi_comm_world,ierr)
      !
      ! assign fluid property(perfect gas)
      !
      g_mwi = mwi
      g_cpi = cpi !g_r * gamma / (gamma - one)
      g_zi = zi
      g_rho0 = rho0
      g_hrefi = hrefi
      g_htref = htref
      g_itype = itype
      !
      g_sigma = sigma
      g_ne = ne
      g_zmu = zmu
      g_tref = tref
      if (s_isource /= 5) then
         g_sref = sref / tref
      else
         g_sref = sref
      end if
      g_pr = pr
      if (all(abs(gravity) < mytiny)) then
         g_igrav = 0
      else
         g_igrav = 1
         g_gravity = gravity
      end if
      g_absorptivity = absorptivity
      species % name = species_name
      g_cdest = cdest; g_cprod = cprod;

      g_gm = g_cpi(1) / (g_cpi(1)-g_runi/g_mwi(1))
      g_r = g_runi / g_mwi(1)
      !
      !  g_prt     = 0.9_rfp
      !
      select case(s_irealgas)
      case(1)
         call read_fluid_table(pfile, ifmt, id) ! real gas
      case(2, 3)
         call init_reaction(pfile(1)) ! NASA gas table
         if (any(g_mwi < mytiny)) g_mwi = species % mw
      case(7)
         call read_RPTBDB(pfile)
      end select
      !
   end subroutine fluid_property
   !
   subroutine solid_property
      implicit none
      ! fluid property
      real :: rhol, rhos, cpl, cps, lheat_melting, lheat_evaporation, kl, ks, mul, beta, sigma, sigmat, k0, dl, alpha, emi
      real :: tl, ts, gv(ndim), n, k, Mw, pb, Tb
      integer :: Darcy_on, Thermal_Buoyancy_on, Gravity_on, Enthalpy_Correction_on 
      namelist /solid_prop/rhol, rhos, cpl, cps, lheat_melting, lheat_evaporation, kl, ks, mul, beta, sigma, sigmat, k0, dl, alpha, emi, tl, ts, gv, n, k, Mw, pb, Tb, Darcy_on, Thermal_Buoyancy_on, Gravity_on, Enthalpy_Correction_on
      integer :: id, ierr

      call mpi_comm_rank(mpi_comm_world, id, ierr)
      rhol = 6518._rfp; rhos = 7870._rfp
      cpl = 800._rfp; cps = 640._rfp
      lheat_melting = 2.7196e5_rfp; lheat_evaporation = 6e6_rfp
      kl = 43.99_rfp; ks = 40.96_rfp
      mul = 0.0050_rfp
      beta = 1.45e-4_rfp
      sigma = 1.6_rfp
      sigmat = -0.00049_rfp
      k0 = 1.e-10_rfp
      dl = 2.e-18_rfp
      alpha = 0.1_rfp
      emi = 0.3_rfp
      tl = 1795._rfp; ts = 1750._rfp
      gv = g_gravity; n=3_rfp; k=4_rfp; Mw=56_rfp 
      pb=1.01e5; Tb=3135
      Darcy_on = 0
      Thermal_Buoyancy_on = 0
      Gravity_on = 0
      Enthalpy_Correction_on = 1
      read (11, solid_prop)
      !xxxxxxxxxxxxxx
      !write(32, solid_prop)
      metal % rhol = rhol
      metal % rhos = rhos
      metal % cpl = cpl
      metal % cps = cps
      metal % lheat_melting = lheat_melting
      metal % lheat_evaporation = lheat_evaporation
      metal % kl = kl
      metal % ks = ks
      metal % mul = mul
      metal % beta = beta
      metal % sigma = sigma
      metal % sigmat = sigmat
      metal % k0 = k0
      metal % dl = dl
      metal % alpha = alpha
      metal % emi = emi
      metal % tl = tl
      metal % ts = ts
      metal % gravity = gv
      metal % n = n
      metal % k = k
      metal % Mw = Mw
      metal % pb = pb; metal % Tb = Tb
      metal % Darcy_on = Darcy_on
      metal % Gravity_on = Gravity_on
      metal % Thermal_Buoyancy_on = Thermal_Buoyancy_on
      metal % Enthalpy_Correction_on = Enthalpy_Correction_on
      !
   end subroutine solid_property

   subroutine read_fluid_table(pfile, ifmt, id)
      character(*) :: pfile(:)
      integer, intent(in) :: id
      integer :: np, nt, nc, ierr, n, ifmt
      integer :: ibuff(3), i
      logical :: yes
      equivalence(ibuff(1), np)
      equivalence(ibuff(2), nt)
      equivalence(ibuff(3), nc)

      do i = 1, nspe
         !  if(id == 0) then
         ibuff = 0
         inquire(file = pfile(i), exist = yes)
         select case(ifmt)
         case(0) ! separate files
            open (unit = 12, file = pfile(i), status = 'old') !form='binary',
         case(1) ! mhd property 
            open (unit = 12, file = pfile(i), form = 'unformatted', status = 'old')
         case(2) ! single file for multi-species
            if (i > 1) then
               yes = .true.
            else
               open (unit = 12, file = pfile(i), status = 'old') !form='binary',
            end if
         case(3) ! single file +absorption no entropy
            open (unit = 12, file = pfile(i), form = 'unformatted', status = 'old')
         end select
         if (.not.yes) then
            print *, 'Sorry, could not find the property file --->>>', pfile(i)
            call mpi_abort(mpi_comm_world, ierr, ierr)
            stop
         end if
         !  end if

         !   if(id == 0) then
         read(12) np, fluid(i) % p0, fluid(i) % dp
         read(12) nt, fluid(i) % t0, fluid(i) % dt
         !   read(12)fluid(i)%nc,fluid(i)%c0,fluid(i)%dc
         nc = 1 ! just mixture
         if (id == 0) then
            print *, ' Using real gas and look up from property table'
            print *, ' The filename is  ', trim(pfile(i))
            print *, 'Pmin=', exp(fluid(i) % p0), 'Pmax=', exp(fluid(i) % p0 + np * fluid(i) % dp)
            print *, 'Tmin=', fluid(i) % t0, 'Tmax=', fluid(i) % t0 + nt * fluid(i) % dt
         end if
         !  call mpi_bcast(ibuff,size(ibuff),mpi_integer,0,mpi_comm_world,ierr)   
         fluid(i) % np = np
         fluid(i) % nt = nt
         fluid(i) % nc = nc
         !  species(i)%name='fluid'//trim(i2s(i))
         !
         allocate(fluid(i) % rho(np, nt, nc))
         allocate(fluid(i) % rhop(np, nt, nc))
         allocate(fluid(i) % rhot(np, nt, nc))
         !   allocate(fluid(i)%rhoc(np,nt,nc))
         allocate(fluid(i) % h(np, nt, nc))
         allocate(fluid(i) % hp(np, nt, nc))
         allocate(fluid(i) % ht(np, nt, nc))
         !   allocate(fluid(i)%hc(np,nt,nc))
         allocate(fluid(i) % s(np, nt, nc))
         allocate(fluid(i) % mu(np, nt, nc))
         allocate(fluid(i) % lamd(np, nt, nc))
         allocate(fluid(i) % sigma(np, nt, nc))
         allocate(fluid(i) % ne(np, nt, nc))
         allocate(fluid(i) % alpha(np, nt, nc))

         !  if(id == 0) then
         read(12) fluid(i) % rho
         read(12) fluid(i) % rhop
         read(12) fluid(i) % rhot
         !   read(12)fluid%rhoc
         read(12) fluid(i) % h
         read(12) fluid(i) % hp
         read(12) fluid(i) % ht
         !   read(12)fluid(i)%hc  
         read(12) fluid(i) % s
         read(12) fluid(i) % mu
         read(12) fluid(i) % lamd
         if (ifmt == 1) then
            read(12) fluid(i) % sigma
            read(12) fluid(i) % ne
         else if (ifmt == 3) then
            read(12) fluid(i) % alpha
            fluid(i) % sigma = fluid(i) % lamd
            fluid(i) % lamd = fluid(i) % mu
            fluid(i) % mu = fluid(i) % s
            fluid(i) % ne = 1.e22_rfp
         end if

         !   fluid%p = fluid%p*1.e6_rfp
         !   fluid%mu = fluid%mu * 1.e-7_rfp
         !   fluid%lamd = fluid%lamd * 1.e-3_rfp
         !  end if

         n = ibuff(1) * ibuff(2) * ibuff(3)
         !  call mpi_bcast(fluid(i)%p0,rfp,mpi_byte,0,mpi_comm_world,ierr)
         !  call mpi_bcast(fluid(i)%dp,rfp,mpi_byte,0,mpi_comm_world,ierr)
         !  call mpi_bcast(fluid(i)%t0,rfp,mpi_byte,0,mpi_comm_world,ierr)
         !  call mpi_bcast(fluid(i)%dt,rfp,mpi_byte,0,mpi_comm_world,ierr)
         !!

         !  call mpi_bcast(fluid(i)%rho,n*rfp,mpi_byte,0,mpi_comm_world,ierr)
         !  call mpi_bcast(fluid(i)%rhop,n*rfp,mpi_byte,0,mpi_comm_world,ierr)
         !  call mpi_bcast(fluid(i)%rhot,n*rfp,mpi_byte,0,mpi_comm_world,ierr)  
         !!
         !  call mpi_bcast(fluid(i)%h,n*rfp,mpi_byte,0,mpi_comm_world,ierr)
         !  call mpi_bcast(fluid(i)%hp,n*rfp,mpi_byte,0,mpi_comm_world,ierr)
         !  call mpi_bcast(fluid(i)%ht,n*rfp,mpi_byte,0,mpi_comm_world,ierr)
         !!
         !  call mpi_bcast(fluid(i)%s,n*rfp,mpi_byte,0,mpi_comm_world,ierr)
         !  call mpi_bcast(fluid(i)%mu,n*rfp,mpi_byte,0,mpi_comm_world,ierr)
         !  call mpi_bcast(fluid(i)%lamd,n*rfp,mpi_byte,0,mpi_comm_world,ierr)
         if (ifmt /= 2) close(12)

      end do ! end loop species
      if (ifmt == 2) close(12)

      !  if(s_isource == 5) then
      !    allocate(fluid(1)%satp(nt))
      !!  if(id == 0)  then
      !     read(12)fluid(1)%satp
      !     close(12)
      !!  end if
      !!   call mpi_bcast(fluid(1)%satp,nt*rfp,mpi_byte,0,mpi_comm_world,ierr)
      !  end if
   end subroutine read_fluid_table

   subroutine read_RPTBDB(pfile)
      character(*) :: pfile(:)
      integer :: ibuff(3), i
      character*80 :: pfil
      logical :: yes

      do i = 1, nspe
         pfil = pfile(i)
         pfil = pfil(:len_trim(pfil))
         inquire(file = pfil, exist = yes)
         if (.not.yes) then
            print *, 'Sorry, could not find the property file --->>>', pfil
            stop
         end if
         RB(i) % p => ReadDatabase(pfil)
      enddo

   end subroutine read_RPTBDB

   subroutine turbulence_model_paramter
      implicit none
      ! Turbulence model
      real :: lamda, cmu, ksplus, cdes
      real(rfp) :: rbuff(4)
      equivalence(rbuff(1), lamda)
      equivalence(rbuff(2), cmu), (rbuff(3), ksplus)
      equivalence(rbuff(4), cdes)
      namelist /turb/lamda, cmu, ksplus, cdes
      integer :: id, ierr
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      !
      !   if(id == 0) then
      rbuff = zero
      cdes = -one
      read (11, turb)
      !xxxxxxxxxxxxxx
      !write(32, turb)
      !   end if
      !
      ! call mpi_bcast(rbuff,size(rbuff)*rfp,mpi_byte,0,mpi_comm_world,ierr)
      !
      !  
      t_omegainf = lamda ! temperally store this constant
      t_keinf = cmu
      t_cdes = cdes ! ->0.86 to use DES model
      !
      if (ksplus < 25.0) then
         t_sr = (50.0_rfp / ksplus)**2
      else
         t_sr = 100.0_rfp / ksplus
      end if
      !
   end subroutine turbulence_model_paramter

   !
   !************************************************************
   subroutine initial(cells, faces, nodes, bnodes, interf)
      !************************************************************
      implicit none
      type(vector) :: qv
      type(cell), pointer :: cells(:), cc
      type(face), pointer :: faces(:), cf
      type(node), pointer :: nodes(:)
      type(bnode), pointer :: bnodes(:)
      type(itf) :: interf
      integer :: i, k, islast, j, ne(10), ii
      integer :: ivar_o(2 * neq), ivar_d(neq), n_o, n_d, label
      real(rfp) :: p, t, u, kin, omega, yi(nspe), zmu, mach, alpha, beta, v(ndim)
      real(rfp) :: xl, yl, zl, xh, yh, zh, lh(6), b(3), e(3), fi(2), l
      character(len = 80) :: location
      !
      namelist /initial_condition/p, t, u, yi, kin, omega, mach, alpha, beta, &
         xl, yl, zl, xh, yh, zh, islast, b, e, fi, &
         ivar_o, ivar_d, n_o, n_d, location, label
      integer :: id, ierr, n, np
      logical :: yes
      real(rfp)::ymax,ls,p_gravity,r,xc,yc,x,y,ra,vortex_strength,vtheta
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      call mpi_comm_size(mpi_comm_world, np, ierr)

      islast = 0
      k = 0
      do while (islast == 0)
         k = k + 1
         !
         ! if(id == 0) then

         ! Defaults
         p = g_patm !1atm = 101.325 kPa
         t = g_troom ! room temperature
         u = 0.0
         mach = 0.0
         yi = 0.0; yi(1) = 1.0;
         alpha = 0.0; beta = 0.0;
         omega = 0.0;
         b = zero; e = zero; fi = zero
         xl = -myhuge; yl = -myhuge; zl = -myhuge;
         xh = myhuge; yh = myhuge; zh = myhuge; ! set domain to (-infinite, infinite)
         label = 0
         islast = 1
         location = ''
         n_o = neq; n_d = neq;
         ivar_o(:neq) = (/(k, k = 1, neq)/)
         ivar_d(:neq) = (/(k, k = 1, neq)/)
         !
         ! Read in initial condition 
         read(11, initial_condition)
         !xxxxxxxxxxxxxx
         !write(32, initial_condition)
         !
         lh(:) = (/xl, yl, zl, xh, yh, zh/)
         if (neqf > 0) then
            if (nmeq == 0) then
               qv % v(ien) = t
            else
               qv % v(ico) = p - g_pbase
               qv % v(ien) = t
               !
               u = max(mach * sqrt(sound_speed(qv)), u)
               alpha = alpha * pi / 180.0_rfp
               beta = beta * pi / 180.0_rfp
               v(1) = cos(alpha) * cos(beta)
               if(ndim/=1) v(2) = sin(alpha) * cos(beta)
               if (ndim == 3) v(ndim) = sin(beta)
               !
               qv % v(imb:ime) = v * u
               if (nspe > 1) qv % v(isb:ise) = yi(:nspm1)
               if (nrad > 0) qv % v(irad) = 4. * stephan_boltzmann * t**4
               !
               !   zmu = mumix(qv)
            end if

            if (s_ivis == 2) then
               if (omega < mytiny) then
                  zmu = mumix(qv)
                  omega = t_omegainf * max(u, one) / b_lenref
                  kin = t_keinf * zmu * omega / rhof(qv)
               end if
               qv % v(iom) = omega
               qv % v(ike) = kin
            end if
         end if
         !
         !  Maxwell's Equ
         if (neqm > 0) then
            qv % v(ibb:ibe) = b
            qv % v(ieb:iee) = e
            qv % v(ifb:ife) = fi
         end if
         !
         ! call mpi_bcast(islast,1,mpi_integer,0,mpi_comm_world,ierr) 
         ! call mpi_bcast(qv%v,neq*rfp,mpi_byte,0,mpi_comm_world,ierr) 
         ! call mpi_bcast(lh,size(lh)*rfp,mpi_byte,0,mpi_comm_world,ierr) 
         if (k == 1) then
            if (s_init == 3.and.s_init == 6) cycle
            do i = 1, size(cells)
               cells(i) % qv = qv ! first
               !   cells(i)%qv%v(ibe) = bthetaf(0._rfp,0.1_rfp,100._rfp,cells(i)%centp(2))
            end do
         else
            do i = 1, size(cells)
               if (label /= 0) then
                  if (vc(cells(i) % itype) % label == label) then
                     cells(i) % qv = qv
                  end if
               else
                  if (any(cells(i) % centp/b_lenref < lh(:ndim))) cycle
                  if (any(cells(i) % centp/b_lenref > lh(4:3 + ndim))) cycle
                  cells(i) % qv = qv
               end if
            end do
         end if

      end do
      ! ghost cell
      do i = 1, size(faces)
         cf => faces(i)
         if (cf % itype == 0) exit
         cf % right_cell % qv = cf % left_cell % qv
      end do

      do j = 1, size(bc)
         if (bc(j) % igrp /= inlet_type) cycle
         if (bc(j) % itype == 6) bc(j) % n = 0
         if (bc(j) % itype == 40) bc(j) % n = 0 ! inlet profile
         do i = 1, size(faces)
            cf => faces(i)
            if (cf % itype == 0) exit
            if (cf % itype /= j) cycle
            if (bc(j) % itype == 6) then
               n = bc(cf % itype) % n
               cf % right_cell % qv % v(:ien) = bc(cf % itype) % vars(n + 1:n + ien)
               b = zero; b(s_iaxis) = s_omega
               b = acrossb(b, avto3(cf % centp))
               cf % right_cell % qv % v(imb:ime) = -b(:ndim)
               bc(j) % n = bc(j) % n + neq
            else if (bc(j) % itype == 40) then
               n = bc(cf % itype) % n
               cf % right_cell % qv % v(:neqf) = bc(cf % itype) % vars(n + ico:n + neqf)
               bc(j) % n = bc(j) % n + neq
            else
               cf % right_cell % qv % v(:neqf) = bc(cf % itype) % vars(:neqf)
            end if
         end do
      end do
      !
      if (s_nstart > 0) then
         call read_binary(cells, faces, id, neq)
      else

         select case(s_init)
         case(1)
            call read_binary(cells, faces, id, neq)
            !    do i = 1, size(cells)
            !     if(cells(i)%centp(1) < 0.008) cycle
            !     if(cells(i)%centp(2) < 0.005) cycle
            !     if(cells(i)%centp(2) < 0.002) cycle
            !     if(cells(i)%centp(2) > 0.005) cycle
            !     cells(i)%qv%v(ien) = 400.
            !     cells(i)%qv%v(isb:ise) = (/0.,0.,0.,0.0,0.,1.,0./)
            !    end do
         case(20) ! user choice which variables need to be used from old solution
            call read_binary_interp(location, cells, faces, id, n_o, ivar_o(:n_d), ivar_d(:n_d))
         case(2) ! from laminar -> turbulence
            call read_binary(cells, faces, id, neq - 2)
            if (neq - 2 > ien) then
               do i = 1, size(cells)
                  cells(i) % qv % v(isb:ise) = cells(i) % qv % v(isb - 2:ise - 2)
               end do
               do i = 1, size(faces)
                  if (faces(i) % itype == 0) exit
                  faces(i) % right_cell % qv % v(isb:ise) = faces(i) % right_cell % qv % v(isb - 2:ise - 2)
               end do
            end if
            do i = 1, size(cells)
               p = sqrt(sum(cells(i) % qv % v(imb:ime)**2)) !/ b_vinf
               cells(i) % qv % v(ike) = 10.0_rfp !p * t_keinf
               cells(i) % qv % v(iom) = 1.e6 !sqrt(p) * t_omegainf / b_lenref
            end do
         case(3) ! initial from tecplot interpolate
            inquire(file = 'gems.init.' // i2s(id), exist = yes)
            if (.not.yes) then
               print *, 'Sorry, could not find the gems.init file'
               call mpi_abort(mpi_comm_world, ierr, ierr)
               stop
            end if
            open(5, file = 'gems.init.' // i2s(id))
            do i = 1, neq + 4
               read(5, *)
            end do
            read(5, *) (nodes(i) % qv % v(:), i = 1, size(nodes))
            do i = 1, size(cells)
               do k = 1, neq
                  cells(i) % qv % v(k) = sum(nodes(cells(i) % c2n(:)) % qv % v(k))
               end do
               cells(i) % qv % v = cells(i) % qv % v / real(size(cells(i) % c2n), rfp)
            end do
            !   close(5)
         case(30) ! initial from tecplot interpolate
            inquire(file = 'gems.init.' // i2s(id), exist = yes)
            if (.not.yes) then
               print *, 'Sorry, could not find the gems.init file'
               call mpi_abort(mpi_comm_world, ierr, ierr)
               stop
            end if
            open(5, file = 'gems.init.' // i2s(id))
            read(5, *) n, ne(:n)
            read(5, *) (nodes(i) % qv % v(ne(:n)), i = 1, size(nodes))
            do i = 1, size(cells)
               do k = 1, n
                  cells(i) % qv % v(ne(k)) = sum(nodes(cells(i) % c2n(:)) % qv % v(ne(k)))
               end do
               cells(i) % qv % v(ne(:)) = cells(i) % qv % v(ne(:)) / real(size(cells(i) % c2n), rfp)
            end do
            !   close(5)
            !
         case(6) ! correct initial condition
            open(5, file = 'gems.bin.' // i2s(id), form = 'unformatted')
            read(5, end = 100) k !s_nstart !,s_cfl
            read(5, end = 100) (cells(i) % qv, i = 1, size(cells))
            100 continue
         case(10)
            call read_binary(cells, faces, id, neq)
            call user_defined_init(cells, faces)
         case(11)
            call interpolate_from_coarsegrid(cells, faces, nodes, bnodes)

         case(9)
            !1D Tube Test
            do i = 1, size(cells)
               if(cells(i)%centp(1)<-one_third .or. cells(i)%centp(1)>one_third) then
                  !cells(i)%qv%v(ico)=1000*9.8*(2.-cells(i)%centp(1))
                  !cells(i)%qv%v(imb)=0.
                  !cells(i)%qv%v(ien)=300.
                  cells(i)%qv%v(ico)=0.
                  cells(i)%qv%v(imb)=0.
                  cells(i)%qv%v(ien)=300.
               else
                  !cells(i)%qv%v(ico)=19600-1.1613689*9.8*cells(i)%centp(1)
                  !cells(i)%qv%v(ico)=19600
                  !cells(i)%qv%v(imb)=0.
                  !cells(i)%qv%v(ien)=300.
                  cells(i)%qv%v(ico)=0.
                  cells(i)%qv%v(imb)=0.
                  cells(i)%qv%v(ien)=300.
               end if
            end do
         case(91)
            ! lamb vortex
            xc = 0.2_rfp
            yc = 0.1_rfp
            ra = 0.03_rfp
            vortex_strength = 1e-4
            do i=1,size(cells)
               x=cells(i)%centp(1); y=cells(i)%centp(2)
               r=sqrt((x-xc)**2+(y-yc)**2)
               if(r<=1e-10) r=1e-10
               vtheta=vortex_strength*(1-exp(-r**2/ra**2))/r
               cells(i)%qv%v(imb)=-vtheta*(y-yc)/r; cells(i)%qv%v(imb+1)=vtheta*(x-xc)/r
            end do
         case(92)
            ! inviscid duct
            l = 1e10
            do i = 1, size(cells)
               if(norm2(cells(i)%centp) < l) then
                  l = norm2(cells(i)%centp)
                  ii = i
               end if
            end do

            cells(ii) % qv % v(ico) = cells(ii) % qv % v(ico) + &
               0.1 * norm2(cells(ii) % qv % v(imb:ime)) ** 2
         case(93)
            ! viscous duct
            do i = 1, size(cells)
               yc = cells(i) % centp(2)
               cells(i) % qv % v(imb) = 1.5 * 3.14 * 4 * (yc / 0.5 - (yc / 0.5)**2)
               cells(i) % qv % v(ime) = zero
            end do
         case(94)
            ! 2D Riemann Problem
            do i = 1, size(cells)
               cc => cells(i)
               xc = cc % centp(1)
               yc = cc % centp(2)
               if(xc <= 0.5 .and. yc <= 0.5) then
                  cc % qv % v = (/0.029d0, 1.206d0, 1.206d0, 7.33004919238d-4/)
               else if(xc > 0.5 .and. yc <= 0.5) then
                  cc % qv % v = (/0.3d0, 0.0d0, 1.206d0, 0.00196574259d0/)
               else if(xc > 0.5 .and. yc > 0.5) then
                  cc % qv % v = (/1.5d0, 0.0d0, 0.0d0, 0.00348788261d0/)
               else
                  cc % qv % v = (/0.3d0, 1.206d0, 0.0d0, 0.00196574259d0/)
               end if
            end do

         end select
      end if
      !
      ! if(s_init == 1) then
      !  lh(:3) = (/-82.,5.69,-100000.0/)
      !  lh(4:6) = (/-81.9,5.75,100000.0/)
      !  do i = 1, size(cells)
      !   if(any(cells(i)%centp/b_lenref < lh(:ndim))) cycle
      !   if(any(cells(i)%centp/b_lenref > lh(4:3+ndim))) cycle
      !   cells(i)%qv%v(ien) = qv%v(ien)
      !  end do
      ! end if

      ! if(s_idt /= 0) then
      !   s_elapsed_time = zero
      !  do i = 1, size(cells)
      !   cells(i)%qvn = cells(i)%qv
      !   if(s_idt > 1) cells(i)%qvm = cells(i)%qv
      !  end do
      !  if(s_nstart > 0 ) then
      !   read(5,end=10)(cells(i)%qvn%v,i=1,size(cells))
      !  if(s_idt > 1) read(5,end=10)(cells(i)%qvm%v,i=1,size(cells))
      !  end if
      ! end if  
      !10 close(5)

      if (s_idt > 0) then
         do j = 1, s_idt
            do i = 1, size(cells)
               cells(i) % qvn(j) = cells(i) % qv
            end do
         end do

         !read qvn from binary files
         if (s_nstart > 0) then
            do j=1, s_idt
               do i=1, size(cells)
                  read(5, end = 10) cells(i) % qvn(j)
               end do
            end do
         end if
      end if

      10 close(5)

      if (id == 0) then
         if (s_nstart == 0) call unlink('gems.res.dat')
         open(7, file = 'gems.res.dat')

         if (s_nstart > 0) then
            do i = 1, s_nstart
               read(7, *, end = 22)
            end do
            22 backspace(7) !ifile)
            s_nstart = s_nstart + 1
            s_nend = s_nstart + s_nend - 1
         else
            s_nstart = 1
         end if
      end if
      call mpi_bcast(s_nstart, 1, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(s_nend, 1, mpi_integer, 0, mpi_comm_world, ierr)
      ! if(s_nstart > 1000) s_alpha = one
      s_omega = min(one, real(s_nstart, rfp)/s_damping_o) * s_omega0

   end subroutine initial
   !
   subroutine read_binary_interp(location, cells, faces, id, n_o, ivar_o, ivar_d)
      type(cell), pointer :: cells(:)
      type(face), pointer :: faces(:)
      integer :: i, id, n_o, k, n, ierr
      real(rfp) :: vars(n_o), skip
      integer :: ivar_o(:), ivar_d(:)
      character(len = *) :: location
      logical :: yes
      n = min(n_o, neq)
      if (size(ivar_o) /= size(ivar_d)) print *, ' interpolate size not match', &
         size(ivar_o), size(ivar_d)
      if (size(ivar_o) > n) print *, 'interpolate variables large than data file'
      inquire(file = trim(location) // 'gems.bin.' // i2s(id), exist = yes)
      if (.not.yes) then
         print *, 'Sorry, could not find-->', trim(location) // 'gems.bin.' // i2s(id)
         call mpi_abort(mpi_comm_world, ierr, ierr)
      end if
      open(5, file = trim(location) // 'gems.bin.' // i2s(id), form = 'unformatted')
      read(5) k
      if (s_nstart >= 1) s_nstart = k
      do i = 1, size(cells)
         read(5) vars
         cells(i) % qv % v(ivar_d) = vars(ivar_o)
         !   (cells(i)%qv%v(:n),skip,skip,i=1,size(cells))
      end do
      do i = 1, size(faces)
         if (faces(i) % itype == 0) exit
         !    read(5,end=10)faces(i)%right_cell%qv%v(:n),skip,skip
         read(5) vars
         faces(i) % right_cell % qv % v(ivar_d) = vars(ivar_o)
      end do
   end subroutine read_binary_interp
   !
   subroutine read_binary(cells, faces, id, n)
      type(cell), pointer :: cells(:)
      type(face), pointer :: faces(:)
      integer :: i, id, n, k

      open(5, file = 'gems.bin.' // i2s(id), form = 'unformatted')

      read(5, end = 10) k ! k=nadv-1 from last time 
      if (s_nstart >= 1) s_nstart = k
      read(5, end = 10) s_elapsed_time

      if(id==0) then
         print*,'reading from binary'
         print*,'s_nstart=',s_nstart
         print*, 's_elapsed_time=',s_elapsed_time
      end if

      do i=1, size(cells)
         read(5, end = 10) cells(i) % qv % v(:n)
      end do

      do i=1, size(faces)
         if (faces(i) % itype == 0) exit
         read(5, end = 10) faces(i) % right_cell % qv % v(:n)
      end do

      return
      10 print *, ' Nothing read in from gems.bin'
   end subroutine read_binary

   !************************************************************
   subroutine user_defined_init(cells, faces)
      type(cell), pointer :: cells(:)
      type(face), pointer :: faces(:)
      type(vector), pointer :: qv
      real(rfp) :: r, r1, r2, ci, p(ndim), btheta, x(12), cu(11), x0, current
      integer :: i, k
      x = (/231.783, 247.658, 250.916, 266.791, 270.05, 285.925, &
         480.522, 496.41, 499.535, 515.517, 518.882, 534.503/)
      cu = (/17.32, 0.0, 8.678, 0.0, 2.255, 0.0, -7.729, 0.0, -9.494, 0.0, -10.93/)
      x = x * b_lenref
      !
      !    ci = g_mu0*bc(1)%vars(1) * half / pi
      do i = 1, size(cells)
         qv => cells(i) % qv
         !     qv%v = zero
         p = cells(i) % centp
         r = sqrt(sum(p(2:3)**2))
         x0 = p(1)
         if (x0 <= x(1).or.x0 >= x(12)) cycle
         current = zero
         do k = 2, 12
            if (x0 < x(k)) exit
            current = current + cu(k - 1)
         end do
         r1 = (x0 - x(k - 1)) / (x(k) - x(k - 1))
         r1 = min(max(zero, r1), one)
         current = r1 * cu(k - 1) + current
         ci = current * g_mu0 / (two * pi * r)
         qv % v(ibb + 1) = -ci * p(3) / r
         qv % v(ibe) = ci * p(2) / r
         qv % v(ieb) = zero
         !      ci = cu(k-1) / (( x(k) - x(k-1)) * two * pi * r) * e_resistivity(qv,cells(i))
         !      qv%v(ieb+1) =   ci * p(2) / r
         !      qv%v(iee  ) =   ci * p(3) / r
      end do
   end subroutine user_defined_init

   !************************************************************
   subroutine read_boundary(cells, faces, nodes, interf)
      !************************************************************
      implicit none
      type(cell), pointer :: cells(:)
      type(face), pointer :: faces(:)
      type(node), pointer :: nodes(:)
      type(itf) :: interf
      integer :: n_total, n_in, n_out, n_far, n_wall, n_geom, n, n_vol, n_mhd
      integer :: i, j, k
      namelist /number_of_boundary/n_in, n_out, n_far, n_wall, n_geom, n_vol, n_mhd
      integer :: ibuff(7)
      equivalence(ibuff(1), n_in), (ibuff(2), n_out), (ibuff(3), n_far)
      equivalence(ibuff(4), n_wall), (ibuff(5), n_geom), (ibuff(6), n_vol)
      equivalence(ibuff(7), n_mhd)
      real(rfp) :: total_area
      integer :: id, ierr
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      ! 
      ! if(id == 0) then
      ibuff = 0
      !
      read(11, number_of_boundary)
      !xxxxxxxxxxxxxx
      !write(32, number_of_boundary)
      ! end if

      ! call mpi_bcast(ibuff,size(ibuff),mpi_integer,0,mpi_comm_world,ierr)
      !
      n_total = n_in + n_out + n_far + n_wall + n_geom + n_mhd

      allocate(bc(n_total))
      ! if(id == 0) then
      k = 0
      call read_inlet(n_in, k)
      !
      !print *,'inlet',id
      call read_outlet(n_out, k)
      !print *,'outlet',id
      !
      call read_farfield(n_far, k)
      !print *,'farfield',id
      !
      call read_wall(n_wall, k)
      ! print *,'wall',id
      !
      call read_geom_boundary(n_geom, k)
      !print *,'geom',id
      !
      call read_mhd_boundary(n_mhd, k)
      !print *,'MHD',id

      if (k /= n_total) print *, 'boundary number is not matched', k, n_total
      if (id == 0) print *, 'reading boundary conditions are finished'
      !end if
      ! 
      ! do i = 1, n_total
      !  if(id == 0) n = size(bc(i)%vars)
      !  call mpi_bcast(bc(i)%label,1,mpi_integer,0,mpi_comm_world,ierr)
      !  call mpi_bcast(bc(i)%itype,1,mpi_integer,0,mpi_comm_world,ierr)
      !  call mpi_bcast(bc(i)%igrp,1,mpi_integer,0,mpi_comm_world,ierr)
      !  call mpi_bcast(n,1,mpi_integer,0,mpi_comm_world,ierr)
      !  if(id /= 0) allocate(bc(i)%vars(n))
      !  call mpi_bcast(bc(i)%vars,n*rfp,mpi_byte,0,mpi_comm_world,ierr) 
      ! end do
      !! some reference values
      ! call mpi_bcast(b_p0,rfp,mpi_byte,0,mpi_comm_world,ierr)
      ! call mpi_bcast(b_pinf,rfp,mpi_byte,0,mpi_comm_world,ierr)
      ! call mpi_bcast(b_rhoinf,rfp,mpi_byte,0,mpi_comm_world,ierr)
      ! call mpi_bcast(b_vinf,rfp,mpi_byte,0,mpi_comm_world,ierr)
      ! call mpi_bcast(b_re,rfp,mpi_byte,0,mpi_comm_world,ierr)

      ! make face%itype equal to bc%label
      ! count no. of boundary faces for each boundary patch
      ! cal. boundary area for each boundary patch
      bc % area = zero
      bc % n = 0
      do i = interf % nitf + 1, size(faces)
         if (faces(i) % itype == 0) exit
         do j = 1, n_total
            if (faces(i) % itype == bc(j) % label) then
               faces(i) % itype = j
               bc(j) % n = bc(j) % n + 1
               if (ndim == 2.and.s_iaxis > 0) then
                  bc(j) % area = bc(j) % area + faces(i) % area * &
                     pi * abs(nodes(faces(i) % f2n(1)) % xyz(s_iaxis) + nodes(faces(i) % f2n(2)) % xyz(s_iaxis))
               else
                  bc(j) % area = bc(j) % area + faces(i) % area
               end if

               exit
            end if
         end do
         if (j > n_total) then
            print *, 'Can not match the boundary type', faces(i) % itype
            print *, 'Please check the boundary label in your grid and gems.inp'
            call mpi_abort(mpi_comm_world, ierr, ierr)
            stop
         end if
      end do

      call read_boundary_data_file(faces, bc)

      ! do i = 1, n_total
      !  total_area = zero
      !  call mpi_allreduce(bc(i)%area,total_area,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
      !  if(total_area < mytiny) total_area = one
      !  bc(i)%area = total_area
      !  if(bc(i)%igrp == inlet_type) then
      !    if(bc(i)%itype == 3 .or. bc(i)%itype == 30) &
      !      bc(i)%vars(neq+3) = bc(i)%vars(neq+3) / total_area
      !  end if
      ! end do


      ! s_isource = n_vol
      ! if(n_vol == 0) then
      !  return
      ! else
      !  if(s_isource ==0) s_isource = 1
      ! end if
      if (n_vol == 0) return
      allocate(vc(n_vol))
      ! if( id == 0) then
      k = 0
      call read_volume_conditions(n_vol, k)
      ! end if

      ! do i = 1, n_vol
      !  if(id == 0) n = size(vc(i)%vars)
      !  call mpi_bcast(vc(i)%label,1,mpi_integer,0,mpi_comm_world,ierr)
      !  call mpi_bcast(vc(i)%itype,1,mpi_integer,0,mpi_comm_world,ierr)
      !  call mpi_bcast(vc(i)%igrp,1,mpi_integer,0,mpi_comm_world,ierr)
      !  call mpi_bcast(n,1,mpi_integer,0,mpi_comm_world,ierr)
      !  if(id /= 0) allocate(vc(i)%vars(n))
      !  call mpi_bcast(vc(i)%vars,n*rfp,mpi_byte,0,mpi_comm_world,ierr) 
      ! end do
      !
      do i = 1, size(cells)
         do j = 1, n_vol
            if (cells(i) % itype == vc(j) % label) then
               cells(i) % itype = j
               exit
            end if
         end do
         if (j > n_vol) then
            print *, 'Can not match the boundary type', cells(i) % itype
            print *, 'Please check the cell label in your grid and dfd.input'
            call mpi_abort(mpi_comm_world, ierr, ierr)
            stop
         end if
      end do

   end subroutine read_boundary

   subroutine read_inlet(n, k)
      integer, intent(in) :: n
      integer :: i, k, label, itype
      real :: p, t, u, mach, v(ndim), yi(nspe), zmu, alpha, beta, vd(ndim), omega, kin, rhou, mass
      real :: on_time, off_time, vswirl, current, r1, r2, bthickness, ycenter
      type(vector) :: qv
      namelist /inlet/label, itype, p, t, mach, u, yi, alpha, beta, omega, kin, rhou, mass, &
         on_time, off_time, vswirl, current, r1, r2, bthickness, ycenter
      if (n <= 0) return
      do i = 1, n
         ! Defaults
         itype = 0
         p = g_patm !1atm = 101.325 kPa
         t = g_troom ! room temperature
         mach = .1 ! m represent mach number
         u = -1.e30_rfp ! magnitude of velocity
         yi = 0.0; yi(1) = 1.0
         alpha = 0.0; beta = 0.0
         omega = zero
         on_time = one; off_time = one
         vswirl = zero
         bthickness = one
         ycenter = zero
         mass = zero
         !
         k = k + 1
         read(11, inlet)
         !xxxxxxxxxxxxxx
         !write(32, inlet)
         !
         allocate(bc(k) % vars(neq + 12 + ndim))
         !
         qv % v(ico) = p - g_pbase
         qv % v(ien) = t
         if (nspe > 1) then
            qv % v(isb:ise) = yi(:nspm1)
         end if
         !
         alpha = alpha * pi / 180.0_rfp
         beta = beta * pi / 180.0_rfp

         vd(1) = cos(alpha) * cos(beta)
         if(ndim/=1) vd(2) = sin(alpha) * cos(beta)
         if (ndim == 3) vd(ndim) = sin(beta)

         if (u < -1.e20) u = mach * sqrt(sound_speed(qv))
         qv % v(imb:ime) = vd * u
         !
         if (itype == 3.or.itype == 30) then
            if (mass <= zero) print *, 'please specify mass flow (Kg/m^2/s)'
            !     u = rhou / rhof(qv)   ! use mass flow
            !     qv%v(imb:ime) = vd * rhou 
         end if
         !
         zmu = mumix(qv)
         if (s_ivis == 2) then
            s_visb = 0
            if (omega < mytiny) then
               omega = t_omegainf * max(u, 1.e-2_rfp) / b_lenref
               kin = t_keinf * zmu * omega / rhof(qv)
               s_visb = 1
            end if
            qv % v(iom) = omega
            qv % v(ike) = kin
            !    print *,'k and omega in inlet =',kin,omega
         end if

         if (s_iaux == 1) then
            qv % v(iaub) = vswirl
         end if
         ! Maxwell
         if (neqm > 0) then
            qv % v(ibb) = current
            qv % v(ibb + 1) = r1 * b_lenref
            qv % v(ibe) = r2 * b_lenref
         end if

         ! some reference parameters
         !
         b_p0 = entropyf(qv)
         b_pinf = p
         b_rhoinf = rhof(qv)
         b_vinf = u
         b_re = b_rhoinf * u / zmu

         !   store to bc database
         !
         bc(k) % label = label
         bc(k) % itype = itype
         bc(k) % igrp = inlet_type

         !    if(itype == 6) then    ! inlet fan
         !     call read_fan_inlet(bc(k),qv)
         !     return
         !    end if

         bc(k) % vars(:neq) = qv % v
         !
         bc(k) % vars(neq + 1) = h0f(qv) !h+1/2v^2=CpT+1/2v^2
         bc(k) % vars(neq + 2) = entropyf(qv) ! entropy
         bc(k) % vars(neq + 3) = u
         bc(k) % vars(neq + 4) = mass
         select case(mod(itype, 1000))
         case(3) ! mass flow
            bc(k) % vars(neq + 3) = mass
            !    case(4)   ! swirl flow
            !     bc(k)%vars(neq + 1) = vswirl
         case(21) ! pipe flow velocity profile
            bc(k) % vars(neq + 1) = bthickness * b_lenref
            bc(k) % vars(neq + 2) = ycenter * b_lenref
            bc(k) % vars(neq + 4 + ndim) = r1 * b_lenref
         case(30) ! unsteady flow with on/off switch
            bc(k) % vars(neq + 3) = mass
            bc(k) % vars(neq + 1) = on_time
            off_time = on_time + off_time
            bc(k) % vars(neq + 2) = off_time
         end select
         bc(k) % vars(neq4:neq3d) = vd ! direction    
         !
      end do
      ! 
   end subroutine read_inlet

   subroutine read_boundary_data_file(faces, bc)
      type(face), pointer :: faces(:)
      type(bc_type) :: bc(:)
      integer :: i, it
      do i = 1, size(bc)
         it = mod(bc(i) % itype, 1000)
         select case(bc(i) % igrp)
         case(inlet_type)
            select case(it)
            case(6)
               call read_inlet_fan(faces, i, bc(i))
            case(40)
               call read_inlet_profile_2dto3d(faces, i, bc(i))
            end select
         case(wall_type)
            if (it == 20) &
               call read_wall_temperature(faces, i, bc(i))
         case(outlet_type)
            if (it == 5) &
               call read_back_pressure(faces, i, bc(i))
         end select
      end do
   end subroutine read_boundary_data_file

   subroutine read_back_pressure(faces, it, bc)
      type(face), pointer :: faces(:), cf
      type(vector) :: qv
      type(bc_type) :: bc
      integer :: it, i, j, n, nf, iaxis
      real(rfp), pointer :: r(:), p(:), b(:), c(:), d(:)
      real(rfp) :: pref, p0, r0
      if (bc % n == 0) return ! no such boundary
      !
      open(8, file = 'back_pressure.dat')
      read(8, *) n, pref
      allocate(r(n), p(n), b(n), c(n), d(n))
      read(8, *)
      do i = 1, n
         ! t
         read(8, *) r(i), p(i)
         r(i) = b_lenref * r(i)
         p(i) = pref * p(i)
      end do
      call spline(r, p, b, c, d)
      close(8)
      nf = bc % n
      deallocate(bc % vars)
      allocate(bc % vars(nf)) ! store variables on bface
      j = 0
      do i = 1, size(faces)
         cf => faces(i)
         if (cf % itype == 0) exit
         if (cf % itype /= it) cycle
         r0 = sqrt(sum(cf % centp**2) - cf % centp(s_iaxis)**2)
         p0 = eval_spline(r0, r, p, b, c, d)
         j = j + 1
         bc % vars(j) = p0 - g_pbase
      end do ! end loop of faces
      deallocate(r, p, b, c, d)
   end subroutine read_back_pressure

   subroutine read_wall_temperature(faces, it, bc)
      type(face), pointer :: faces(:), cf
      type(vector) :: qv
      type(bc_type) :: bc
      integer :: it, i, j, n, nf, iaxis
      real(rfp), pointer :: r(:), t(:), b(:), c(:), d(:)
      real(rfp) :: tref, t0, r0
      if (bc % n == 0) return ! no such boundary
      !
      open(8, file = 'wall_temp.dat')
      read(8, *) n, iaxis, tref
      allocate(r(n), t(n), b(n), c(n), d(n))
      read(8, *)
      do i = 1, n
         ! t
         read(8, *) r(i), t(i)
         r(i) = b_lenref * r(i)
         t(i) = tref * t(i)
      end do
      call spline(r, t, b, c, d)
      close(8)
      nf = bc % n
      deallocate(bc % vars)
      allocate(bc % vars(neq * nf)) ! store variables on bface
      j = 0
      do i = 1, size(faces)
         cf => faces(i)
         if (cf % itype == 0) exit
         if (cf % itype /= it) cycle
         r0 = cf % centp(iaxis) ! just for 2D now
         !      sqrt(sum(cf%centp**2) - cf%centp(s_iaxis)**2)
         t0 = eval_spline(r0, r, t, b, c, d)
         qv % v = zero
         !     qv%v(ico) = p - g_pbase
         qv % v(ien) = t0
         !     qv%v(imb:ime) = vd(:ndim)
         bc % vars(j + ico:j + ien) = qv % v(:ien)
         j = j + neq
      end do ! end loop of faces
      deallocate(r, t, b, c, d)
   end subroutine read_wall_temperature

   subroutine read_inlet_profile_2dto3d(faces, it, bc)
      type(face), pointer :: faces(:), cf
      type(vector) :: qv
      type(bc_type) :: bc
      integer :: it, i, n, j, nf, m
      real(rfp), pointer :: r(:), pt(:,:), b(:,:), c(:,:), d(:,:)
      real(rfp) :: xcos, xsin, vr(3), vd(3), p, t, ux, uy, uz, ke, om, ur, r0
      nf = bc % n
      if (nf == 0) return
      open(8, file = 'inlet.dat')
      read(8, *) n, m ! n : no. of points, m: no. of variables
      allocate(r(n), pt(n, m), b(n, m), c(n, m), d(n, m))
      read(8, *)
      do i = 1, n
         ! p0,t0,ar,at
         read(8, *) r(i), (pt(i, j), j = 1, m)
         r(i) = b_lenref * r(i)
      end do
      !
      do j = 1, m
         call spline(r, pt(:, j), b(:, j), c(:, j), d(:, j))
      end do
      close(8)
      deallocate(bc % vars)
      allocate(bc % vars(neq * nf)) ! store variables on bface
      j = 0
      do i = 1, size(faces)
         cf => faces(i)
         if (cf % itype == 0) exit
         if (cf % itype /= it) cycle
         r0 = sqrt(sum(cf % centp(2:3)**2))
         p = eval_spline(r0, r, pt(:, 1), b(:, 1), c(:, 1), d(:, 1))
         ux = eval_spline(r0, r, pt(:, 2), b(:, 2), c(:, 2), d(:, 2))
         ur = eval_spline(r0, r, pt(:, 3), b(:, 3), c(:, 3), d(:, 3))
         t = eval_spline(r0, r, pt(:, 4), b(:, 4), c(:, 4), d(:, 4))
         ke = eval_spline(r0, r, pt(:, 5), b(:, 5), c(:, 5), d(:, 5))
         om = eval_spline(r0, r, pt(:, 6), b(:, 6), c(:, 6), d(:, 6))
         !
         xcos = cf % centp(2) / r0
         xsin = cf % centp(3) / r0
         ! convert to Cartesian Coord.
         vd(1) = ux
         vd(2) = ur * xcos ! vr*cos + vq*sin  clockwise
         vd(3) = ur * xsin ! vr*sin - vq*cos  clockwise
         !
         qv % v(ico) = p - g_pbase
         qv % v(imb:ime) = vd(:ndim)
         qv % v(ien) = t
         qv % v(ike) = ke
         qv % v(iom) = om
         !
         bc % vars(j + ico:j + neqf) = qv % v(:neqf)
         j = j + neq
      end do ! end loop of faces
      deallocate(r, pt, b, c, d)

   end subroutine read_inlet_profile_2dto3d

   subroutine read_inlet_fan(faces, it, bc)
      type(face), pointer :: faces(:), cf
      type(vector) :: qv
      type(bc_type) :: bc
      integer :: it, i, n, j, nf
      real(rfp), pointer :: r(:), pt(:,:), b(:,:), c(:,:), d(:,:)
      real(rfp) :: xcos, xsin, vr(3), vd(3), p, t, beta, alpha, r0
      open(8, file = 'inlet.dat')
      read(8, *) n
      allocate(r(n), pt(n, 4), b(n, 4), c(n, 4), d(n, 4))
      read(8, *)
      do i = 1, n
         ! p0,t0,ar,at
         read(8, *) r(i), (pt(i, j), j = 1, 4)
         r(i) = b_lenref * r(i)
      end do
      ! dimensional values
      pt(:, 1) = pt(:, 1) * (bc % vars(ico) + g_pbase)
      pt(:, 2) = pt(:, 2) * bc % vars(ien)
      !
      do j = 1, 4
         call spline(r, pt(:, j), b(:, j), c(:, j), d(:, j))
      end do
      close(8)
      nf = bc % n
      deallocate(bc % vars)
      allocate(bc % vars(neq * nf)) ! store variables on bface
      j = 0
      do i = 1, size(faces)
         cf => faces(i)
         if (cf % itype == 0) exit
         if (cf % itype /= it) cycle
         r0 = sqrt(sum(cf % centp**2) - cf % centp(s_iaxis)**2)
         p = eval_spline(r0, r, pt(:, 1), b(:, 1), c(:, 1), d(:, 1)) !* g_patm
         t = eval_spline(r0, r, pt(:, 2), b(:, 2), c(:, 2), d(:, 2)) !* g_troom
         beta = eval_spline(r0, r, pt(:, 3), b(:, 3), c(:, 3), d(:, 3)) * pi / 180 ! Meridional
         alpha = eval_spline(r0, r, pt(:, 4), b(:, 4), c(:, 4), d(:, 4)) * pi / 180 ! Tangential
         xcos = cf % centp(2) / r0
         xsin = cf % centp(3) / r0
         !   Cylindrical coordinates
         vr(1) = cos(alpha) * cos(beta) ! x-
         vr(2) = sin(alpha) * cos(beta) ! theta
         vr(3) = sin(beta) ! r
         ! convert to Cartesian Coord.
         vd(1) = vr(1)
         vd(2) = vr(3) * xcos + vr(2) * xsin ! vr*cos + vq*sin  clockwise
         vd(3) = vr(3) * xsin - vr(2) * xcos ! vr*sin - vq*cos  clockwise
         !
         qv % v(ico) = p - g_pbase
         qv % v(ien) = t
         qv % v(imb:ime) = vd(:ndim)
         bc % vars(j + ico:j + ien) = qv % v(:ien)
         j = j + neq
      end do ! end loop of faces

      deallocate(r, pt, b, c, d)

   end subroutine read_inlet_fan

   subroutine read_farfield(n, k)
      integer, intent(in) :: n
      integer :: i, k, label, itype
      real :: p, t, u, mach, v(ndim), yi(nspe), zmu, alpha, beta, vd(ndim), omega, kin
      type(vector) :: qv
      namelist /farfield/label, itype, p, t, mach, u, yi, alpha, beta, omega, kin
      if (n <= 0) return
      do i = 1, n
         ! Defaults
         itype = 0
         p = g_patm !1atm = 101.325 kPa
         t = g_troom ! room temperature
         mach = 0. ! m represent mach number
         u = -1.e30_rfp ! magnitude of velocity
         yi = 0.0; yi(1) = 1.0
         alpha = 0.0; beta = 0.0
         omega = zero
         !
         k = k + 1
         read(11, farfield)
         !xxxxxxxxxxxxxx
         !write(32, farfield)
         allocate(bc(k) % vars(neq + 5 + ndim))
         !
         qv % v(ico) = p - g_pbase
         qv % v(ien) = t
         if (nspe > 1) then
            qv % v(isb:ise) = yi(:nspm1)
         end if
         !    if( itype == 0) then  ! use p,t,m
         if (u < -1.e20) u = mach * sqrt(sound_speed(qv))
         !    else if(itype == 1) then ! use p,t,u
         !    else if(itype == 2) then ! use mass flow
         !    end if
         !
         alpha = alpha * pi / 180.0_rfp
         beta = beta * pi / 180.0_rfp
         vd(1) = cos(alpha) * cos(beta)
         vd(2) = sin(alpha) * cos(beta)
         if (ndim == 3) vd(ndim) = sin(beta)
         qv % v(imb:ime) = vd * u
         !
         if (s_ivis == 2) then
            if (omega < mytiny) then
               zmu = mumix(qv)
               omega = t_omegainf * u / b_lenref
               kin = t_keinf * zmu * omega / rhof(qv)
            end if
            qv % v(ike) = max(1.e-3_rfp, kin)
            qv % v(iom) = max(1.e3_rfp, omega)
         end if

         !   store to bc database
         !
         bc(k) % label = label
         bc(k) % itype = itype
         bc(k) % igrp = farfield_type
         !

         bc(k) % vars(:neq) = qv % v
         !
         bc(k) % vars(neq + 1) = h0f(qv) !h+1/2v^2=CpT+1/2v^2
         bc(k) % vars(neq + 2) = entropyf(qv) ! entropy
         bc(k) % vars(neq + 3) = u
         bc(k) % vars(neq4:neq3d) = vd ! direction    
         !
      end do
      ! 
   end subroutine read_farfield


   subroutine read_outlet(n, k)
      integer, intent(in) :: n
      integer :: i, label, itype, k, ni
      real :: pback, mass, r1, r2
      namelist /outlet/label, itype, pback, ni, r1, r2
      if (n <= 0) return
      itype = 0
      pback = g_patm !1atm = 101.325 kPa
      ni = 1
      label = 1
      mass = zero
      r1 = zero
      r2 = zero

      do i = 1, n
         k = k + 1
         read(11, outlet)
         !xxxxxxxxxxxxxx
         !write(32, outlet)
         allocate(bc(k) % vars(10))
         bc(k) % vars(4:7) = -19621011._rfp ! initialize
         !
         !   store to bc database
         !
         bc(k) % label = label
         bc(k) % itype = itype
         bc(k) % igrp = outlet_type
         !    
         bc(k) % vars(1) = pback - g_pbase
         bc(k) % vars(2) = ni + 0.1
         bc(k) % vars(3) = mass
         bc(k) % vars(4) = r1 * b_lenref
         bc(k) % vars(5) = r2 * b_lenref
      end do
      ! 
   end subroutine read_outlet
   !   
   subroutine read_wall(n, k)
      integer, intent(in) :: n
      integer :: i, k, label, itype
      real :: wvel(ndim), wtemp, ks, qw, rcenter(2), omega, vswirl, current, r1, r2
      namelist /wall/label, itype, wvel, wtemp, ks, qw, rcenter, omega, vswirl, current, r1, r2
      if (n <= 0) return
      do i = 1, n
         !  Default values
         itype = 0
         wvel = 0.0; wtemp = 297.0_rfp; ks = 0.0
         vswirl = 0.0
         current = zero
         qw = zero ! qw = 0 means adiabatic wall
         r1 = zero
         r2 = zero
         !
         k = k + 1
         read(11, wall)
         !xxxxxxxxxxxxxx
         !write(32, wall)
         allocate(bc(k) % vars(ndim + 10))
         !
         !   store to bc database
         !
         bc(k) % label = label
         bc(k) % itype = itype
         bc(k) % igrp = wall_type

         bc(k) % vars(:ndim) = wvel
         bc(k) % vars(ndim + 1) = wtemp
         bc(k) % vars(ndim + 2) = ks
         bc(k) % vars(ndim + 3) = qw
         bc(k) % vars(ndim + 4) = vswirl
         bc(k) % vars(ndim + 5) = current
         bc(k) % vars(ndim + 6) = r1 * b_lenref
         bc(k) % vars(ndim + 7) = r2 * b_lenref
         if (itype == 5.or.itype == 6) then !rotation
            bc(k) % vars(:2) = rcenter
            bc(k) % vars(3) = omega * two * pi
            bc(k) % vars(4) = wtemp
         end if
      end do
      ! 
   end subroutine read_wall

   subroutine read_geom_boundary(n, k)
      integer, intent(in) :: n
      integer :: i, k, label, itype, ni, id
      real(rfp) :: rcenter(ndim), angle
      namelist /geom_type/label, itype, rcenter, angle, ni, id
      if (n <= 0) return
      do i = 1, n
         !  Default values
         itype = 0; ni = 0; id = 0
         rcenter = zero
         angle = zero
         !
         k = k + 1
         read(11, geom_type)
         !xxxxxxxxxxxxxx
         !write(32, geom_type)
         allocate(bc(k) % vars(4))
         !
         !   store to bc database
         !
         bc(k) % label = label
         bc(k) % itype = itype ! = 0 symmetry 1 Axisymmetry 2 periodic
         bc(k) % igrp = geom_topology
         select case(itype)
         case(2)
            s_pcenter = rcenter
            bc(k) % vars(1:ndim) = rcenter
            bc(k) % vars(4) = angle * pi / 180.0_rfp
         case(4) ! 2d to 3d interface
            bc(k) % vars(1) = ni
            bc(k) % vars(2) = id
         end select
      end do
      ! 
   end subroutine read_geom_boundary

   subroutine read_mhd_boundary(n, k)
      integer, intent(in) :: n
      integer :: i, k, label, itype
      real(rfp) :: rin, rout, current, current0, ev(3)
      namelist /mhd_btype/label, itype, rin, rout, current, ev, current0
      if (n <= 0) return
      do i = 1, n
         !  Default values
         itype = 0
         current = zero; current0 = zero
         rin = zero; rout = zero;
         !
         k = k + 1
         read(11, mhd_btype)
         !xxxxxxxxxxxxxx
         !write(32, mhd_btype)
         allocate(bc(k) % vars(10))
         !
         !   store to bc database
         !
         bc(k) % label = label
         bc(k) % itype = itype ! = 0 the end of dommain 1 outside of domain
         bc(k) % igrp = mhd_type
         bc(k) % vars(1) = current
         bc(k) % vars(ifb) = current0
         bc(k) % vars(2) = rin * b_lenref
         bc(k) % vars(3) = rout * b_lenref
         bc(k) % vars(ieb:iee) = ev
      end do
      ! 
   end subroutine read_mhd_boundary

   subroutine read_volume_conditions(n, k)
      integer, intent(in) :: n
      integer :: i, k, label, itype
      real(rfp) :: f(ndim), hx, sigma, ne, speed, x0, dx
      namelist /volume_con/label, itype, f, hx, sigma, ne, speed, x0, dx
      if (n <= 0) return
      do i = 1, n
         !  Default values
         itype = 0
         f = zero
         hx = zero
         ne = 1.e12_rfp; sigma = one
         !
         k = k + 1
         read(11, volume_con)
         !xxxxxxxxxxxxxx
         !write(32, volume_con)
         allocate(vc(k) % vars(10)) !ndim+1))
         !
         !   store to bc database
         !
         vc(k) % label = label
         vc(k) % itype = itype !0 non 1 only f 2 only hx 3 both
         vc(k) % igrp = 1
         vc(k) % vars(:ndim) = f
         vc(k) % vars(ndim + 1) = hx
         !    vc(k)%vars(1) = sigma
         !    vc(k)%vars(2) = ne
         !    vc(k)%vars(3) = speed
         !    vc(k)%vars(4) = x0
         !    vc(k)%vars(5) = dx
      end do
      ! 
   end subroutine read_volume_conditions
   !
   subroutine output_option
      implicit none
      !
      integer :: Perturbation_Pressure, Static_Pressure, Total_Pressure
      integer :: pitot_pressure, Velx, Vely, Velz
      integer :: Temperature, Total_Temperature
      integer :: Species_mass_fraction, Species_volume_fraction, Species_Concentration
      integer :: Turbulent_kinetic, Turbulent_dissipation, Turbulent_Viscosity
      integer :: Mach_Number, Sound_speed
      INteger :: Density, Enthalpy, Entropy
      integer :: Total_Vorticity, Vorticity
      integer :: Viscosity, Thermo_Conductivity
      integer :: intensity, shadowgraph
      integer :: bx, by, bz, ex, ey, ez, fib, fie, jx, jy, jz, electrical_conductivity
      integer :: lorentz_x, lorentz_y, lorentz_z, je, absorptivity
      integer :: RHeatFlux_x, RHeatFlux_y, RHeatFlux_z, DivE, dBxdT, dBydt, dBzdt, dExdT, dEydt, dEzdt
      integer :: reaction_rate
      integer :: ibuff(100)
      equivalence(ibuff(1), Perturbation_Pressure), &
         (ibuff(2), Static_Pressure), &
         (ibuff(3), Total_Pressure), &
         (ibuff(4), pitot_pressure), &
         (ibuff(5), Velx), &
         (ibuff(6), Vely), &
         (ibuff(7), Velz), & !7
         (ibuff(8), Temperature), &
         (ibuff(9), Total_Temperature), & !9
         (ibuff(10), Turbulent_kinetic), &
         (ibuff(11), Turbulent_dissipation), &
         (ibuff(12), Turbulent_Viscosity), & !12
         (ibuff(13), Mach_Number), &
         (ibuff(14), Sound_speed), & ! 14
         (ibuff(15), Density), &
         (ibuff(16), Enthalpy), &
         (ibuff(17), Entropy), & ! 17
         (ibuff(18), Total_Vorticity), &
         (ibuff(19), Vorticity), & ! 19
         (ibuff(20), Species_mass_fraction), &
         (ibuff(21), Species_volume_fraction), &
         (ibuff(22), Species_Concentration), &
         (ibuff(23), Viscosity), &
         (ibuff(24), Thermo_Conductivity), &
         (ibuff(25), Intensity), (ibuff(26), Shadowgraph), &
         (ibuff(27), bx), (ibuff(28), by), (ibuff(29), bz), &
         (ibuff(30), ex), (ibuff(31), ey), (ibuff(32), ez), &
         (ibuff(33), fib), (ibuff(34), fie), &
         (ibuff(35), jx), (ibuff(36), jy), (ibuff(37), jz), &
         (ibuff(38), Electrical_Conductivity), &
         (ibuff(39), lorentz_x), (ibuff(40), lorentz_y), &
         (ibuff(41), lorentz_z), (ibuff(42), je), &
         (ibuff(43), absorptivity), &
         (ibuff(44), RHeatFlux_x), &
         (ibuff(45), RHeatFlux_y), &
         (ibuff(46), RHeatFlux_z), &
         (ibuff(47), divE), &
         (ibuff(48), dBxdT), &
         (ibuff(49), dBydT), &
         (ibuff(50), dBzdT), &
         (ibuff(51), dExdT), &
         (ibuff(52), dEydT), &
         (ibuff(53), dEzdT), &
         (ibuff(54), reaction_rate)
      namelist /output_vars/&
         Perturbation_Pressure, Static_Pressure, Total_Pressure, &
         pitot_pressure, Velx, Vely, Velz, &
         Temperature, Total_Temperature, &
         Species_mass_fraction, Species_volume_fraction, Species_Concentration, &
         Turbulent_kinetic, Turbulent_dissipation, Turbulent_Viscosity, &
         Mach_Number, Sound_speed, &
         Density, Enthalpy, Entropy, &
         Total_Vorticity, Vorticity, &
         Viscosity, Thermo_Conductivity, &
         Intensity, shadowgraph, bx, by, bz, ex, ey, ez, fib, fie, jx, jy, jz, &
         electrical_conductivity, lorentz_x, lorentz_y, lorentz_z, je, nwplot, ivplot, &
         absorptivity, RHeatFlux_x, RHeatFLux_y, RHeatFlux_z, DivE, &
         dBxdt, dBydt, dBzdt, dExdt, dEydt, dEzdt, reaction_rate
      !
      integer :: id, ierr, i
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      vname(:54) = (/&
         'Perturbation_Pressure  ', &
         'Static_Pressure        ', &
         'Total_Pressure         ', &
         'Pitot_pressure         ', &
         'Velx                   ', &
         'Vely                   ', &
         'Velz                   ', &
         'Temperature            ', &
         'Total_Temperature      ', &
         'Turbulent_kinetic      ', &
         'Turbulent_dissipation  ', &
         'Turbulent_Viscosity    ', &
         'Mach_Number            ', &
         'Sound_speed            ', &
         'Density                ', &
         'Enthalpy               ', &
         'Entropy                ', &
         'Total_Vorticity        ', &
         'Vorticity              ', &
         'Species_mass_fraction  ', &
         'Species_volume_fraction', &
         'Species_Concentration  ', &
         'Viscosity              ', &
         'Thermo_Conductivity    ', &
         'Intensity              ', &
         'Shadowgraph            ', &
         'bx                     ', &
         'by                     ', &
         'bz                     ', &
         'ex                     ', &
         'ey                     ', &
         'ez                     ', &
         'fib                    ', &
         'fie                    ', &
         'jx                     ', &
         'jy                     ', &
         'jz                     ', &
         'Electrical_Conductivity', &
         'Lorentz_x              ', &
         'Lorentz_y              ', &
         'Lorentz_z              ', &
         'je                     ', &
         'Absorptivity           ', &
         'RadiationHeatFlux_X    ', &
         'RadiationHeatFlux_Y    ', &
         'RadiationHeatFlux_Z    ', &
         'DivE                   ', &
         'dBxdt                  ', &
         'dBydt                  ', &
         'dBzdt                  ', &
         'dExdt                  ', &
         'dEydt                  ', &
         'dEzdt                  ', &
         'Reaction_Rate          '/)
      !  if(id == 0) then
      !  defaults
      ibuff = 0
      nwplot = 0
      ivplot = 0
      if (neqf > 0) then
         !  Static_pressure = 1
         !Velx = 1; Vely = 1; if (ndim == 3) Velz = 1;
         if (nmeq == 0) ibuff = 0
         !Temperature = 1
      end if
      if (neqm > 0) then
         bz = 1; ex = 1; ey = 1;
      end if
      read(11, output_vars)
      !xxxxxxxxxxxxxx
      !write(32, output_vars)
      !
      !  end if
      !
      !  call mpi_bcast(ibuff,size(ibuff),mpi_integer,0,mpi_comm_world,ierr)
      output_selection = .false.
      if (id == 0) print *, ' You choose following variables as output'
      if (s_ivis /= 2) ibuff(12) = 0
      if (s_idt == 0) ibuff(48:53) = 0

      do i = 1, size(ibuff)
         if (ibuff(i) == 1) then
            output_selection(i) = .true.
            if (id == 0) print *, trim(vname(i))
         end if
      end do
      !
   end subroutine output_option

   subroutine woutput_option(faces)
      implicit none
      type(face), pointer :: faces(:)
      !
      integer :: pressure, temperature, density, heat_flux, taw, utaw, yplus, label
      integer :: vx, vy, vz, rheatflux, ex, ey, ez, enthalpy, stagnation_pressure
      integer :: stagnation_temperature, levelset
      integer :: ibuff(18)
      equivalence(ibuff(1), Pressure), &
         (ibuff(2), Temperature), &
         (ibuff(3), Density), &
         (ibuff(4), taw), &
         (ibuff(5), utaw), & ! 17
         (ibuff(6), yplus), &
         (ibuff(7), heat_flux), &
         (ibuff(8), Vx), &
         (ibuff(9), Vy), &
         (ibuff(10), Vz), &
         (ibuff(11), rheatflux), &
         (ibuff(12), ex), &
         (ibuff(13), ey), &
         (ibuff(14), ez), &
         (ibuff(15), enthalpy), &
         (ibuff(16), stagnation_pressure), &
         (ibuff(17), stagnation_temperature), &
         (ibuff(18), levelset)

      namelist /woutput_vars/&
         label, pressure, temperature, density, heat_flux, taw, utaw, yplus, &
         vx, vy, vz, rheatflux, ex, ey, ez, enthalpy, &
         stagnation_pressure, stagnation_temperature, levelset
      !
      integer :: id, ierr, i, n, iw
      if (nwplot == 0) return
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      wname(:18) = (/&
         'Static_Pressure        ', &
         'Temperature            ', &
         'Density                ', &
         'Taw                    ', &
         'Utaw                   ', &
         'Yplus                  ', &
         'Heat_Flux              ', &
         'Velx                   ', &
         'Vely                   ', &
         'Velz                   ', &
         'RadiationHeatFlux      ', &
         'Ex                     ', &
         'Ey                     ', &
         'Ez                     ', &
         'Enthalpy               ', &
         'Stagnation_pressure    ', &
         'Stagnation_temperature ', &
         'Levelset '&
         /)

      allocate(wlabel(nwplot), woutput_selection(100, nwplot))
      do iw = 1, nwplot
         !  defaults
         label = 0
         ibuff = 0
         !if (nmeq > 0) ibuff(1:2) = 1
         read(11, woutput_vars)
         !xxxxxxxxxxxxxx
         !write(32, woutput_vars)
         !  call mpi_bcast(ibuff,size(ibuff),mpi_integer,0,mpi_comm_world,ierr)
         woutput_selection(:, iw) = .false.
         if (id == 0) print *, ' You cheese following variables as wall plot for label', label
         do i = 1, size(ibuff)
            if (ibuff(i) == 1) then
               woutput_selection(i, iw) = .true.
               if (id == 0) print *, trim(wname(i))
            end if
         end do
         if (label /= 0) then
            n = 0
            do i = 1, size(faces)
               if (faces(i) % itype == 0) exit
               if (label == -1) then
                  if (faces(i) % itype == partition_face_no + 1) n = n + 1
               else
                  if (faces(i) % itype >= partition_face_no) cycle
                  if (bc(faces(i) % itype) % label == label) then
                     n = n + 1
                  end if
               end if
            end do
            if (n == 0) label = 0
         end if
         wlabel(iw) = label
      end do ! end loop of woutput
      !
   end subroutine woutput_option

   subroutine read_schedule
      integer :: i, n, ischeme(5), ialg(5), steps(5), ipre(5)
      real(rfp) :: cfl(5), rbuff(5)
      integer :: id, ierr, ibuff(4 * 5)
      equivalence(ibuff(1), steps(1)), &
         (ibuff(6), ialg(1)), &
         (ibuff(11), ischeme(1)), (ibuff(16), ipre(1))
      equivalence(rbuff(1), cfl(1))
      namelist /job_schedule/ischeme, ialg, steps, cfl, ipre
      if (s_ischedule == 0) then
         schedule % n = 1
         return
      end if
      call mpi_comm_rank(mpi_comm_world, id, ierr)
      !  if(id == 0) then
      steps = 100
      cfl = s_cflmax
      ischeme = s_ischeme
      ialg = s_ialg
      ipre = s_ipre
      read(11, job_schedule)
      !xxxxxxxxxxxxxx
      !write(32, job_schedule)
      !  end if
      !  open(3,file='gems.job')
      !  read(3,*)n
      call mpi_bcast(ibuff, size(ibuff), mpi_integer, 0, mpi_comm_world, ierr)
      !  call mpi_bcast(rbuff,size(rbuff)*rfp,mpi_byte,0,mpi_comm_world,ierr)
      n = min(5, s_ischedule)
      schedule % n = n + 1
      allocate(schedule % steps(n), schedule % ischeme(n), schedule % ialg(n), schedule % cfl(n))
      allocate(schedule % iprecon(n))
      !  do i = 1, n
      !  read(3,*)schedule%steps(i),schedule%cfl(i),schedule%ischeme(i),schedule%ialg(i)
      schedule % steps = steps(:n)
      schedule % ischeme = ischeme(:n)
      schedule % ialg = ialg(:n)
      schedule % iprecon = ipre(:n)
      schedule % cfl = cfl(:n)
      !  end do
   end subroutine read_schedule

   subroutine set_schedule(i)
      integer :: i
      if (i == schedule % n) return
      s_nstart = s_nend + 1
      s_nend = s_nend + schedule % steps(i)
      s_ischeme = schedule % ischeme(i)
      s_ialg = schedule % ialg(i)
      s_ipre = schedule % iprecon(i)
      s_cflmax = schedule % cfl(i)
   end subroutine set_schedule

   subroutine output_cell_tecplot_v(fn, cells, nodes, faces, it, var)
      !********************************************************
      implicit none
      !
      type(cell), pointer :: cells(:), cc
      type(node), pointer :: nodes(:)
      type(face), pointer :: faces(:)

      ! local variables
      type(vector) :: qv
      type(vector), allocatable :: qvn(:)
      real(rfp) :: rc(ndim), rn(ndim), mach, qinf, var(:)
      integer, pointer :: c2n(:), nps(:)
      integer :: i, j, k, n, it
      character(len = *) :: fn
      !  
      open(2, file = trim(fn) // '.' // trim(i2s(it)) // '.dat')
      write(2, *) 'variables ='
      if (ndim == 1) write(2, *) 'x'
      if (ndim == 2) write(2, *) 'x y'
      if (ndim == 3) write(2, *) 'x y z'
      n = ndim + 1
      write(2, *) 'var'
      ! head   
      if (ndim == 2) then
         write(2, *) 'zone n=', size(nodes), ',e=', size(cells), &
            !            ',varlocation=([3-'//trim(i2s(n))//']=cellcentered)',      &
         ',varlocation=([' // trim(i2s(n)) // ']=cellcentered)', &
            ', f=feblock,et=quadrilateral'
      else
         write(2, *) 'zone n=', size(nodes), ',e=', size(cells), ',f=feblock', &
            !            ',varlocation=([4-'//trim(i2s(n))//']=cellcentered),et=brick'
         ',varlocation=([' // trim(i2s(n)) // ']=cellcentered),et=brick'
      end if
      ! coordinates
      do i = 1, ndim
         write(2, 20) nodes(:) % xyz(i) / b_lenref
      end do
      !
      write(2, 20) var
      !  
      20 format(5e20.8)
      do i = 1, size(cells)
         c2n => cells(i) % c2n
         n = size(c2n)
         select case(ndim)
         case(2)
            if (n == 3) then
               write(2, 10) c2n, c2n(3)
            else
               write(2, 10) c2n
            end if
         case(3)
            if (n == 4) then
               write(2, 10) c2n(:3), c2n(3), c2n(4), c2n(4), c2n(4), c2n(4)
            else if (n == 5) then
               write(2, 10) c2n(:4), c2n(5), c2n(5), c2n(5), c2n(5)
            else if (n == 6) then
               write(2, 10) c2n(:3), c2n(3), c2n(4:6), c2n(6)
            else
               write(2, 10) c2n
            end if
         end select
      end do

      close(2)
      10 format(10I8)
   end subroutine output_cell_tecplot_v

   subroutine read_debug()
      implicit none
      integer :: pretime, nadv
      real(rfp) :: xc, yc, zc

      pretime = 0
      xc = zero; yc = zero; zc = zero
      nadv = 1000000

      namelist /debug/pretime, xc, yc, zc, nadv
      read(11, debug)

      db_pretime = pretime
      db_xc = xc
      db_yc = yc
      db_zc = zc
      db_nadv = nadv
   end subroutine read_debug

   subroutine set_monitor(cells, id)
      implicit none
      type(cell), pointer :: cells(:), cc
      type(face), pointer :: cf
      integer :: i, id, ierr, cellid(1), j
      real(rfp) :: dist(1), monp(ndim), d_in(2), d_out(2), area, dx
      real(rfp), allocatable :: dist_array(:), dx_array(:)

      ! find the cell that is closest to the monitor point
      allocate(dist_array(size(cells)))
      allocate(dx_array(size(cells)))

      do i = 1, size(cells)
         if(ndim == 1) then
            dist_array(i) = abs(cells(i)%centp(1) - db_xc)
         else if(ndim == 2) then
            dist_array(i) = sqrt((cells(i)%centp(1)-db_xc)**2 + &
               (cells(i)%centp(2)-db_yc)**2)
         else
            dist_array(i) = sqrt((cells(i)%centp(1)-db_xc)**2 + &
               (cells(i)%centp(2)-db_yc)**2 + (cells(i)%centp(3)-db_zc)**2)
         end if
      end do

      dist = minval(dist_array)

      d_in = (/dist(1), real(id, 8)/)
      call mpi_allreduce(d_in, d_out, 2, MPI_2DOUBLE_PRECISION, MPI_MINLOC, &
         mpi_comm_world, ierr)

      db_pid = IDNINT(d_out(2))
      if(id == db_pid) then
         cellid = minloc(dist_array)
         db_cellid = cellid(1)
         db_moncell => cells(db_cellid)
         print *, 'monitor cell is in partition ', id
         print *, 'monitor cell id = ', db_cellid
         print *, 'cell center = ', db_moncell % centp
      else
         db_cellid = -1
      end if

      ! find the charecteristic mesh cell, save to s_dx
      do i = 1, size(cells)
         cc => cells(i)
         area = -one
         do j = 1, size(cc % sface)
            cf => cc % sface(j) % to_face
            if(cf % area > area) area = cf % area
         end do
         dx_array(i) = cc % vol / area
      end do

      dx = minval(dx_array)
      call mpi_allreduce(dx, s_dx, 1, MPI_DOUBLE, MPI_MIN, mpi_comm_world, ierr)

      if(id == 0) print *, 'characteristic mesh size = ', s_dx
      deallocate(dist_array, dx_array)
   end subroutine set_monitor
end module gems_input

