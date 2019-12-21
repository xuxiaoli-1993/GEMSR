MODULE dfd_levelset
  USE dfd_data
  USE dfd_state
  USE mpi
implicit none
 contains
subroutine read_levelset_parameter
        implicit none
        ! level-set parameters
        real :: thickness, cfl_adv, cfl_RI
        integer :: material, intf_treatment, init, adv_freq, RI_freq, RI_quit, adv_nRK, RI_nRK
        namelist /levelset_parameter/thickness, material, intf_treatment, cfl_adv, cfl_RI, init, adv_freq, RI_freq, RI_quit, &
        adv_nRK, RI_nRK
        
        thickness = 1.5_rfp
	material = 110	
        intf_treatment = 2
	cfl_adv = 0.5_rfp
	cfl_RI = 0.5_rfp
	init = 1
	adv_freq = 10
	RI_freq = 10
	RI_quit = 1
	adv_nRK = 2
	RI_nRK = 2
	
        read(11, levelset_parameter)
        write(32, levelset_parameter)
        
        ls_thickness = thickness
	ls_material = material	
        ls_intf_treatment = intf_treatment
	ls_cfl_adv = cfl_adv
	ls_cfl_RI = cfl_RI
	ls_init = init
	ls_adv_freq = adv_freq
	ls_RI_freq = RI_freq
	ls_RI_quit = RI_quit
	ls_adv_nRK = adv_nRK
	ls_RI_nRK = RI_nRK

end subroutine read_levelset_parameter

subroutine levelset_allocation(cells,faces,interf)
  implicit none
  integer :: ierr, i
  type(cell), pointer :: cells(:), rc
  type(face), pointer :: faces(:)
  type(itf) :: interf
  
  do i = 1,size(cells)
    allocate(cells(i)%ls, cells(i)%ls0, cells(i)%ls1, cells(i)%ls2, cells(i)%gradls(ndim), stat=ierr)
    allocate(cells(i)%lsn, stat=ierr)
    allocate(cells(i)%itypels, stat=ierr)
    allocate(cells(i)%inb(2*ndim), stat=ierr)
    cells(i)%ls=0
    cells(i)%ls0=0
    cells(i)%ls1=0
    cells(i)%ls2=0
    cells(i)%lsn=0
    cells(i)%gradls=0
    cells(i)%itypels=0
    cells(i)%inb=0
  end do
  
  do i=1,size(faces)
    if(faces(i)%itype>0 .and. faces(i)%itype .le. size(BC)) then
      rc => faces(i)%right_cell
      allocate(rc%ls, rc%ls0, rc%ls1, rc%ls2, rc%lsn, stat=ierr)
      allocate(rc%itypels, stat=ierr)
      rc%ls=0
      rc%ls0=0
      rc%ls1=0
      rc%ls2=0
      rc%lsn=0
      rc%itypels=0
    end if
  end do
  
  do i=1,size(interf%scell)
    allocate(interf%scell(i)%inbmore, stat=ierr)
    interf%scell(i)%inbmore=0
  end do
  
  do i=1,size(interf%pcell)
    allocate(interf%pcell(i)%ls, interf%pcell(i)%ls0, interf%pcell(i)%ls1, interf%pcell(i)%ls2, stat=ierr)
    allocate(interf%pcell(i)%morex, interf%pcell(i)%morels, interf%pcell(i)%morels0, stat=ierr)
    allocate(interf%pcell(i)%lsn, stat=ierr)
    allocate(interf%pcell(i)%itypels, stat=ierr)
    interf%pcell(i)%ls=0
    interf%pcell(i)%ls0=0
    interf%pcell(i)%ls1=0
    interf%pcell(i)%ls2=0
    interf%pcell(i)%lsn=0
    interf%pcell(i)%morex=zero
    interf%pcell(i)%morels=zero
    interf%pcell(i)%morels0=zero
    interf%pcell(i)%itypels=0
  end do
  
end subroutine levelset_allocation

subroutine levelset_initialization(cells,faces)
  implicit none
  integer :: id, ierr, i
  type(cell), pointer :: cells(:), rc, lc
  type(face), pointer :: faces(:)
  real(rfp) :: x, y, z

  call mpi_comm_rank(mpi_comm_world, id, ierr)
  do i = 1,size(cells)
    select case(ls_init)
    case default
      cells(i)%ls = zero
    case(1)
      cells(i)%ls = cells(i)%centp(2)
    case(2)
      cells(i)%ls = real(id,rfp)
    case(3)
    !sphere deformation test
    !a sphere at (0.35,0.35,0.35) with a radius of 0.15
    !a deformation velocity field with a period of 3 second
    !simulate 1/5 of the period to test the deformation
      x=cells(i)%centp(1)
      y=cells(i)%centp(2)
      z=cells(i)%centp(3)
      cells(i)%ls = sqrt((x-0.35_rfp)**2+(y-0.35_rfp)**2+(z-0.35_rfp)**2)-0.15_rfp
      cells(i)%qv%v(imb) = 2*dsin(pi*x)*dsin(pi*x)*dsin(2*pi*y)*dsin(2*pi*z)
      cells(i)%qv%v(imb+1) = -dsin(2*pi*x)*dsin(pi*y)*dsin(pi*y)*dsin(2*pi*z)
      cells(i)%qv%v(imb+2) = -dsin(2*pi*x)*dsin(2*pi*y)*dsin(pi*z)*dsin(pi*z)
    end select
    
    cells(i)%ls0=cells(i)%ls
    cells(i)%ls1=cells(i)%ls
    cells(i)%ls2=cells(i)%ls
    cells(i)%lsn=cells(i)%ls
  end do
    
  do i=1,size(faces)
    if(faces(i)%itype>0 .and. faces(i)%itype .le. size(BC)) then
      rc => faces(i)%right_cell
      lc => faces(i)%left_cell
      rc%ls=lc%ls
      rc%ls0=lc%ls0
      rc%ls1=lc%ls1
      rc%ls2=lc%ls2
      rc%lsn=lc%lsn
    end if
  end do 

end subroutine levelset_initialization

subroutine prepare_levelset(cells,interf)
  implicit none
  type(cell), pointer :: cells(:)
  type(itf) :: interf
  integer :: i, ierr
  
  !give value to inb of each cell
  call find_neighbor_direction(cells)
  call tag_inner_cell(cells)
  call gather_more_neighbors(cells,interf)
  call mpi_barrier(mpi_comm_world, ierr)
  call levelset_transfer(cells,interf,1)
  
  ! for diffuse interface only
  if(ls_intf_treatment==1) then
    call levelset_thickness(cells)
  end if
  
end subroutine prepare_levelset

subroutine find_neighbor_direction(cells)
  implicit none
  type(cell), pointer :: cells(:),current_cell,neighbor_cell
  integer :: i,j,k,m,tag
  
  do i=1,size(cells)
    current_cell => cells(i)
    do j=1,size(current_cell%scell)
      neighbor_cell => current_cell%scell(j)%to_cell
      
      if(current_cell%centp(1)-neighbor_cell%centp(1)>1e-12) current_cell%inb(1)=j
      if(current_cell%centp(1)-neighbor_cell%centp(1)<-1e-12) current_cell%inb(2)=j
      
      if(ndim>1) then
	if(current_cell%centp(2)-neighbor_cell%centp(2)>1e-12) current_cell%inb(3)=j
	if(current_cell%centp(2)-neighbor_cell%centp(2)<-1e-12) current_cell%inb(4)=j
      end if
      
      if(ndim>2) then
	if(current_cell%centp(3)-neighbor_cell%centp(3)>1e-12) current_cell%inb(5)=j
	if(current_cell%centp(3)-neighbor_cell%centp(3)<-1e-12) current_cell%inb(6)=j
      end if
    end do
    
    !debug
    do k=1,2*ndim
      tag=0
      do j=1,size(current_cell%inb)
	if(current_cell%inb(j)==k) tag=1
      end do
      if(tag==0) then
	print*,'can not find neighbor in direction: ',k
	print*, current_cell%inb
      end if
    end do
    !debug
    
  end do
end subroutine find_neighbor_direction

subroutine tag_inner_cell(cells)
  implicit none
  type(cell), pointer :: cells(:), neighbor_cell
  type(face), pointer :: neighbor_face
  integer :: i,j,no_BC,id,ierr,ct
  
  !borrow rp
  cells%rp = 0. 
  call mpi_comm_rank(mpi_comm_world, id, ierr)
  no_BC = size(bc)
 
  !tag the outer layer
  do i=1,size(cells)
    do j=1,size(cells(i)%sface)
      neighbor_face => cells(i)%sface(j)%to_face
      !boundary face
      if(neighbor_face%itype > 0 .and. neighbor_face%itype .le. no_BC) then
	cells(i)%itypels = 1
	exit
      end if
    end do
  end do
  
  !tag the next layer using rp
  do i=1,size(cells)
    if(cells(i)%itypels == 1) cycle 
    cells(i)%itypels = 2
  end do

  !debug
  ct=0
  do i=1,size(cells)
    if(cells(i)%itypels == 1) then
      ct = ct+1
    end if
  end do
  
  print*, 'id=',id,' ,levelset outer cell=',ct
  !debug

end subroutine tag_inner_cell

subroutine gather_more_neighbors(cells,interf)

  implicit none
  type(cell), pointer :: cells(:),current_cell,neighbor_cell
  type(itf) :: interf
  integer :: i,j,neighbor_index,this_index,this_coordinate

  do i=1,size(interf%scell)
    current_cell => interf%scell(i)%to_cell
    if(current_cell%itypels == 0) cycle
    
    this_index=0
    
    !find which pcell(i) is associated with interf%scell(i)%to_cell
    !save the direction in this_index
    do j=1,size(current_cell%inb)
      neighbor_index = current_cell%inb(j)
      neighbor_cell => current_cell%scell(neighbor_index)%to_cell
      if(associated(neighbor_cell,interf%pcell(i))) then
	this_index = j
	exit
      end if
    end do

    !find the opposite cell
    if(this_index /= 0) then
      if(mod(this_index,2)==0) then
	interf%scell(i)%inbmore = this_index-1
      else
	interf%scell(i)%inbmore = this_index+1
      end if
    end if

  end do
end subroutine gather_more_neighbors

subroutine levelset_transfer(cells,interf,tag)
  !tag = 1: preparation transer, transfer ls and morex
  !tag = 2: transfer ls
  !tag = 3: transfer ls0
  implicit none
  type(itf) :: interf
  type(cell), pointer :: cells(:), cc, nc
  type(neighbour_cell), pointer :: sc
  real(rfp), allocatable :: sb(:), rb(:)
  integer :: i, n, ib, j, ip, nd, ierr, st(mpi_status_size), ir, meq, tag
  integer :: inbmore1, ndim1, more_index
  
  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_comm_rank(mpi_comm_world, ip, ierr)

  ib = 0; ir = 0
  do i = 1, interf % nc
    
    select case(tag)
    case(1)
      meq = 3
    case(2,3)
      meq = 2
    end select
    
    nd = meq * interf % ncp(i)
    allocate(sb(nd))
    sb = 0
    nd = 1
    do j = 1, interf % ncp(i)
      ib = ib + 1
      !pack information from scell into nd
      cc => interf % scell(ib) % to_cell
      sc => interf % scell(ib)
      inbmore1 = sc%inbmore
      ndim1=0
      if(inbmore1 /= 0) then
	if(inbmore1 < 2.5) then
	  ndim1=1
	else if(inbmore1 < 4.5) then
	  ndim1=2
	else
	  ndim1=3
	end if
      end if

      !intieror cells
      if(cc%itypels /= 0) then
	select case(tag)
	case(1)
	  sb(nd) = cc % ls
	  if(inbmore1 /= 0) then
	    more_index = cc%inb(inbmore1)
	    nc => cc%scell(more_index)%to_cell
	    sb(nd+1)=nc%ls
	    sb(nd+2)=nc%centp(ndim1)
	  end if
	  
	case(2)
	  sb(nd) = cc % ls
	  if(inbmore1 /= 0) then
	    more_index = cc%inb(inbmore1)
	    nc => cc%scell(more_index)%to_cell
	    sb(nd+1)=nc%ls
	  end if
	  
	case(3)
	  sb(nd) = cc % ls0
	  if(inbmore1 /= 0) then
	    more_index = cc%inb(inbmore1)
	    nc => cc%scell(more_index)%to_cell
	    sb(nd+1)=nc%ls0
	  end if
	  
	end select
      end if
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
      select case(tag)
      case(1)
	cc % ls = rb(nd)
	cc % morels = rb(nd+1)
	cc % morex = rb(nd+2)
      case(2)
	cc % ls = rb(nd)
	cc % morels = rb(nd+1)
      case(3)
	cc % ls0 = rb(nd)
	cc % morels0 = rb(nd+1)  
      end select
      nd = nd + meq
    end do
    deallocate(rb, sb)
  end do

  call mpi_barrier(mpi_comm_world, ierr)
end subroutine levelset_transfer

subroutine levelset_advection(cells, interf, nadv)
  implicit none
  type(cell), pointer :: cells(:), cc, nc
  type(cell), pointer :: xlc, xrc, ylc, yrc,zlc, zrc, xllc, xrrc, yllc, yrrc, zllc, zrrc
  type(itf) :: interf
  real(rfp) :: dt_adv, dx, maxu, maxv, maxw, maxu_all, maxv_all, maxw_all, templs
  real(rfp) :: xsten(3,5), fsten(3,5), dd1(3,4), dd2(3,3), dphi(3,2), HJ
  real(rfp) :: u, v, w
  integer :: i, j, ii, jj, n, ind, nadv, id, np, ierr
  integer :: xll, xl, xrr, xr, yll, yl, yrr, yr, zll, zl, zrr, zr
  
  call mpi_comm_rank(mpi_comm_world, id, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)
  
  cc => cells(1)
  ind = cc%inb(1)
  nc => cc%scell(ind)%to_cell
  dx = cc%centp(1) - nc%centp(1)
  
  maxu=maxval(cells%qv%v(imb))
  maxv=maxval(cells%qv%v(imb+1))
  maxw=maxval(cells%qv%v(imb+2))
  
  call mpi_allreduce(maxu, maxu_all, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_world, ierr)
  call mpi_allreduce(maxv, maxv_all, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_world, ierr)
  call mpi_allreduce(maxw, maxw_all, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_world, ierr)
  
  dt_adv=ls_cfl_adv*dx/sqrt((dabs(maxu_all))**2+(dabs(maxv_all))**2+(dabs(maxw_all))**2)
  !grab level-set from the previous real time step
  do i=1,size(cells)
    cells(i)%ls=cells(i)%lsn
  end do
  call levelset_transfer(cells,interf,2)
  
  do n=1, ls_adv_nRK
    do i=1,size(cells)
      if(cells(i)%itypels /= 2) cycle
      
      !find stencils
      cc => cells(i)
      xl = cc%inb(1)
      xr = cc%inb(2)
      yl = cc%inb(3)
      yr = cc%inb(4)
      zl = cc%inb(5)
      zr = cc%inb(6)
      u = cc%qv%v(imb)
      v = cc%qv%v(imb+1)
      w = cc%qv%v(imb+2)
      
      xlc => cc%scell(xl)%to_cell
      xrc => cc%scell(xr)%to_cell
      ylc => cc%scell(yl)%to_cell
      yrc => cc%scell(yr)%to_cell
      zlc => cc%scell(zl)%to_cell
      zrc => cc%scell(zr)%to_cell
      
      xsten(1,2)=xlc%centp(1)
      xsten(1,3)=cc%centp(1)
      xsten(1,4)=xrc%centp(1)
      xsten(2,2)=ylc%centp(2)
      xsten(2,3)=cc%centp(2)
      xsten(2,4)=yrc%centp(2)
      xsten(3,2)=zlc%centp(3)
      xsten(3,3)=cc%centp(3)
      xsten(3,4)=zrc%centp(3)
      
      fsten(1,2)=xlc%ls
      fsten(1,3)=cc%ls
      fsten(1,4)=xrc%ls
      fsten(2,2)=ylc%ls
      fsten(2,3)=cc%ls
      fsten(2,4)=yrc%ls
      fsten(3,2)=zlc%ls
      fsten(3,3)=cc%ls
      fsten(3,4)=zrc%ls
      
      if(xlc%itypels /= 0) then
	xll = xlc%inb(1)
	xllc => xlc%scell(xll)%to_cell
	xsten(1,1)=xllc%centp(1)
	fsten(1,1)=xllc%ls
      else if(associated(xlc%morex)) then
	xsten(1,1) = xlc%morex
	fsten(1,1) = xlc%morels
      else
	print*, 'can not find x left left'
      end if
      
      if(xrc%itypels /= 0) then
	xrr = xrc%inb(2)
	xrrc => xrc%scell(xrr)%to_cell
	xsten(1,5)=xrrc%centp(1)
	fsten(1,5)=xrrc%ls
      else if(associated(xrc%morex)) then
	xsten(1,5) = xrc%morex
	fsten(1,5) = xrc%morels
      else
	print*, 'can not find x right right'
      end if
      
      if(ylc%itypels /= 0) then
	yll = ylc%inb(3)
	yllc => ylc%scell(yll)%to_cell
	xsten(2,1)=yllc%centp(2)
	fsten(2,1)=yllc%ls
      else if(associated(ylc%morex)) then
	xsten(2,1) = ylc%morex
	fsten(2,1) = ylc%morels
      else
	print*, 'can not find y left left'
      end if
	     
      if(yrc%itypels /= 0) then
	yrr = yrc%inb(4)
	yrrc => yrc%scell(yrr)%to_cell
	xsten(2,5)=yrrc%centp(2)
	fsten(2,5)=yrrc%ls
      else if(associated(yrc%morex)) then
	xsten(2,5) = yrc%morex
	fsten(2,5) = yrc%morels
      else
	print*, 'can not find y right right'
      end if
	     
      if(zlc%itypels /= 0) then
	zll = zlc%inb(5)
	zllc => zlc%scell(zll)%to_cell
	xsten(3,1)=zllc%centp(3)
	fsten(3,1)=zllc%ls
      else if(associated(zlc%morex)) then
	xsten(3,1) = zlc%morex
	fsten(3,1) = zlc%morels
      else
	print*, 'can not find z left left'
      end if
	     
      if(zrc%itypels /= 0) then
	zrr = zrc%inb(6)
	zrrc => zrc%scell(zrr)%to_cell
	xsten(3,5)=zrrc%centp(3)
	fsten(3,5)=zrrc%ls
      else if(associated(zrc%morex)) then
	xsten(3,5) = zrc%morex
	fsten(3,5) = zrc%morels
      else
	print*, 'can not find z right right'
      end if
	     
      !calculate divided difference
      do ii=1,3
	do jj=1,4
	  dd1(ii,jj)=(fsten(ii,jj+1)-fsten(ii,jj))/(xsten(ii,jj+1)-xsten(ii,jj))
	end do
      end do
	  
      do ii=1,3
	do jj=1,3
	  dd2(ii,jj)=(dd1(ii,jj+1)-dd1(ii,jj))/(xsten(ii,jj+2)-xsten(ii,jj))
	end do
      end do
  
      do ii=1,3
	dphi(ii,1)=dd1(ii,2)+minmod(dd2(ii,1),dd2(ii,2))*(xsten(ii,3)-xsten(ii,2))
	dphi(ii,2)=dd1(ii,3)+minmod(dd2(ii,2),dd2(ii,3))*(xsten(ii,3)-xsten(ii,4))
      end do
      
      !calcualte the Hamilton-Jacobian by upwind scheme
      HJ=0	    
      if(u>0) then
	HJ=HJ+u*dphi(1,1)
      else
	HJ=HJ+u*dphi(1,2)
      end if
	      
      if(v>0) then
	HJ=HJ+v*dphi(2,1)
      else
	HJ=HJ+v*dphi(2,2)
      end if
	       
      if(w>0) then
	HJ=HJ+w*dphi(3,1)
      else
	HJ=HJ+w*dphi(3,2)
      end if
	     
      !time evolution
      if(n==1) then
	cc%ls1=cc%ls-HJ*dt_adv
      else
	cc%ls2=cc%ls-HJ*dt_adv
      end if
  
    end do
    
    !swap ls and ls1
    if(n==1) then
      do i=1,size(cells)
	templs=cells(i)%ls
	cells(i)%ls=cells(i)%ls1
	cells(i)%ls1=templs
      end do
	
      call levelset_transfer(cells,interf,2)
    end if
	   
    !add time n+2 and n
    if(n==2) then
      do i=1,size(cells)
	cells(i)%ls = 0.5*(cells(i)%ls2 + cells(i)%ls1)
      end do
	
      call levelset_transfer(cells,interf,2)
    end if
  
  end do
  if(id==0) print*, 'ADV', nadv
  
end subroutine levelset_advection

subroutine levelset_reinitialization(cells, interf, nadv)
  implicit none
  type(cell), pointer :: cells(:), cc, nc
  type(cell), pointer :: xlc, xrc, ylc, yrc,zlc, zrc, xllc, xrrc, yllc, yrrc, zllc, zrrc
  type(itf) :: interf
  real(rfp) :: dt_RI, dx, templs,  it_err, it_err0, it_err_all, it_err0_all, bandsize
  real(rfp) :: xsten(3,5), fsten(3,5), fsten0(3,5), dd1(3,4), dd2(3,3), dphi(3,2), HJ
  real(rfp) :: intpx, a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2
  integer :: i, j, ii, jj, n, ind, nadv, id, np, ierr, qt_flag
  integer :: xll, xl, xrr, xr, yll, yl, yrr, yr, zll, zl, zrr, zr
  integer :: nit, ct_bd, ct_bd_all

  call mpi_comm_rank(mpi_comm_world, id, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)
  
  cc => cells(1)
  ind = cc%inb(1)
  nc => cc%scell(ind)%to_cell
  dx = cc%centp(1) - nc%centp(1)
  
  dt_RI=ls_cfl_RI*dx/1.
  qt_flag=0
  nit=0
  bandsize=5*dx
  
  do i=1,size(cells)
    cells(i)%ls0=cells(i)%ls
  end do
    
  call levelset_transfer(cells,interf,3)
  call levelset_transfer(cells,interf,2)
  
  do while(qt_flag .ne. 1)
    do n=1,ls_RI_nRK
      
      do i=1,size(cells)
	if(cells(i)%itypels /= 2) cycle
	
	!find stencils
	cc => cells(i)
	xl = cc%inb(1)
	xr = cc%inb(2)
	yl = cc%inb(3)
	yr = cc%inb(4)
	zl = cc%inb(5)
	zr = cc%inb(6)
	
	xlc => cc%scell(xl)%to_cell
	xrc => cc%scell(xr)%to_cell
	ylc => cc%scell(yl)%to_cell
	yrc => cc%scell(yr)%to_cell
	zlc => cc%scell(zl)%to_cell
	zrc => cc%scell(zr)%to_cell
	
	xsten(1,2)=xlc%centp(1)
	xsten(1,3)=cc%centp(1)
	xsten(1,4)=xrc%centp(1)
	xsten(2,2)=ylc%centp(2)
	xsten(2,3)=cc%centp(2)
	xsten(2,4)=yrc%centp(2)
	xsten(3,2)=zlc%centp(3)
	xsten(3,3)=cc%centp(3)
	xsten(3,4)=zrc%centp(3)
	
	fsten(1,2)=xlc%ls
	fsten(1,3)=cc%ls
	fsten(1,4)=xrc%ls
	fsten(2,2)=ylc%ls
	fsten(2,3)=cc%ls
	fsten(2,4)=yrc%ls
	fsten(3,2)=zlc%ls
	fsten(3,3)=cc%ls
	fsten(3,4)=zrc%ls
	
	fsten0(1,2)=xlc%ls0
	fsten0(1,3)=cc%ls0
	fsten0(1,4)=xrc%ls0
	fsten0(2,2)=ylc%ls0
	fsten0(2,3)=cc%ls0
	fsten0(2,4)=yrc%ls0
	fsten0(3,2)=zlc%ls0
	fsten0(3,3)=cc%ls0
	fsten0(3,4)=zrc%ls0
	
	if(xlc%itypels /= 0) then
	  xll = xlc%inb(1)
	  xllc => xlc%scell(xll)%to_cell
	  xsten(1,1)=xllc%centp(1)
	  fsten(1,1)=xllc%ls
	  fsten0(1,1)=xllc%ls0
	else if(associated(xlc%morex)) then
	  xsten(1,1) = xlc%morex
	  fsten(1,1) = xlc%morels
	  fsten0(1,1) = xlc%morels0
	else
	  print*, 'can not find x left left'
	end if
	
	if(xrc%itypels /= 0) then
	  xrr = xrc%inb(2)
	  xrrc => xrc%scell(xrr)%to_cell
	  xsten(1,5)=xrrc%centp(1)
	  fsten(1,5)=xrrc%ls
	  fsten0(1,5)=xrrc%ls0
	else if(associated(xrc%morex)) then
	  xsten(1,5) = xrc%morex
	  fsten(1,5) = xrc%morels
	  fsten0(1,5) = xrc%morels0
	else
	  print*, 'can not find x right right'
	end if
	
	if(ylc%itypels /= 0) then
	  yll = ylc%inb(3)
	  yllc => ylc%scell(yll)%to_cell
	  xsten(2,1)=yllc%centp(2)
	  fsten(2,1)=yllc%ls
	  fsten0(2,1)=yllc%ls0
	else if(associated(ylc%morex)) then
	  xsten(2,1) = ylc%morex
	  fsten(2,1) = ylc%morels
	  fsten0(2,1) = ylc%morels0
	else
	  print*, 'can not find y left left'
	end if
	      
	if(yrc%itypels /= 0) then
	  yrr = yrc%inb(4)
	  yrrc => yrc%scell(yrr)%to_cell
	  xsten(2,5)=yrrc%centp(2)
	  fsten(2,5)=yrrc%ls
	  fsten0(2,5)=yrrc%ls0
	else if(associated(yrc%morex)) then
	  xsten(2,5) = yrc%morex
	  fsten(2,5) = yrc%morels
	  fsten0(2,5) = yrc%morels0
	else
	  print*, 'can not find y right right'
	end if
	      
	if(zlc%itypels /= 0) then
	  zll = zlc%inb(5)
	  zllc => zlc%scell(zll)%to_cell
	  xsten(3,1)=zllc%centp(3)
	  fsten(3,1)=zllc%ls
	  fsten0(3,1)=zllc%ls0
	else if(associated(zlc%morex)) then
	  xsten(3,1) = zlc%morex
	  fsten(3,1) = zlc%morels
	  fsten0(3,1) = zlc%morels0
	else
	  print*, 'can not find z left left'
	end if
	      
	if(zrc%itypels /= 0) then
	  zrr = zrc%inb(6)
	  zrrc => zrc%scell(zrr)%to_cell
	  xsten(3,5)=zrrc%centp(3)
	  fsten(3,5)=zrrc%ls
	  fsten0(3,5)=zrrc%ls0
	else if(associated(zrc%morex)) then
	  xsten(3,5) = zrc%morex
	  fsten(3,5) = zrc%morels
	  fsten0(3,5) = zrc%morels0
	else
	  print*, 'can not find z right right'
	end if
  
	!sub-cell fix
	!re-define the stencils
	if(fsten0(1,2)*fsten0(1,3)<0) then
	  intpx=cubic_intp(xsten(1,1:4),fsten0(1,1:4))
	  xsten(1,1)=xsten(1,2)
	  xsten(1,2)=intpx
	  fsten(1,1)=fsten(1,2)
	  fsten(1,2)=0.
	end if
	    
	if(fsten0(1,3)*fsten0(1,4)<0) then
	  intpx=cubic_intp(xsten(1,2:5),fsten0(1,2:5))
	  xsten(1,5)=xsten(1,4)
	  xsten(1,4)=intpx
	  fsten(1,5)=fsten(1,4)
	  fsten(1,4)=0.
	end if
	    
	if(fsten0(2,2)*fsten0(2,3)<0) then
	  intpx=cubic_intp(xsten(2,1:4),fsten0(2,1:4))
	  xsten(2,1)=xsten(2,2)
	  xsten(2,2)=intpx
	  fsten(2,1)=fsten(2,2)
	  fsten(2,2)=0.
	end if
	      
	if(fsten0(2,3)*fsten0(2,4)<0) then
	  intpx=cubic_intp(xsten(2,2:5),fsten0(2,2:5))
	  xsten(2,5)=xsten(2,4)
	  xsten(2,4)=intpx
	  fsten(2,5)=fsten(2,4)
	  fsten(2,4)=0.	  
	end if
	      
	if(fsten0(3,2)*fsten0(3,3)<0) then
	  intpx=cubic_intp(xsten(3,1:4),fsten0(3,1:4))
	  xsten(3,1)=xsten(3,2)
	  xsten(3,2)=intpx
	  fsten(3,1)=fsten(3,2)
	  fsten(3,2)=0.
	end if
	      
	if(fsten0(3,3)*fsten0(3,4)<0) then
	  intpx=cubic_intp(xsten(3,2:5),fsten0(3,2:5))
	  xsten(3,5)=xsten(3,4)
	  xsten(3,4)=intpx
	  fsten(3,5)=fsten(3,4)
	  fsten(3,4)=0.	  
	end if
  
	!calculate divided difference
	do ii=1,3
	  do jj=1,4
	    dd1(ii,jj)=(fsten(ii,jj+1)-fsten(ii,jj))/(xsten(ii,jj+1)-xsten(ii,jj))
	  end do
	end do
	    
	do ii=1,3
	  do jj=1,3
	    dd2(ii,jj)=(dd1(ii,jj+1)-dd1(ii,jj))/(xsten(ii,jj+2)-xsten(ii,jj))
	  end do
	end do
  
	do ii=1,3
	  dphi(ii,1)=dd1(ii,2)+minmod(dd2(ii,1),dd2(ii,2))*(xsten(ii,3)-xsten(ii,2))
	  dphi(ii,2)=dd1(ii,3)+minmod(dd2(ii,2),dd2(ii,3))*(xsten(ii,3)-xsten(ii,4))
	end do
  
	!calcualte the Hamilton-Jacobian by Godunov's scheme
	if(cc%ls0>0) then
	  a1=max(dphi(1,1),0.)
	  b2=min(dphi(1,2),0.)
	  c1=max(dphi(2,1),0.)
	  d2=min(dphi(2,2),0.)
	  e1=max(dphi(3,1),0.)
	  f2=min(dphi(3,2),0.)
	  HJ=sqrt(max(a1**2,b2**2)+max(c1**2,d2**2)+max(e1**2,f2**2))-1
	else
	  a2=min(dphi(1,1),0.)
	  b1=max(dphi(1,2),0.)
	  c2=min(dphi(2,1),0.)
	  d1=max(dphi(2,2),0.)
	  e2=min(dphi(3,1),0.)
	  f1=max(dphi(3,2),0.)
	  HJ=sqrt(max(a2**2,b1**2)+max(c2**2,d1**2)+max(e2**2,f1**2))-1
	end if
  
	!use a smoothed sign distance function	  	  
	if(fsten0(1,2)*fsten0(1,3)<0 .or. fsten0(1,3)*fsten0(1,4)<0 .or. &
	  fsten0(2,2)*fsten0(2,3)<0 .or. fsten0(2,3)*fsten0(2,4)<0 .or. &
	  fsten0(3,2)*fsten0(3,3)<0 .or. fsten0(3,3)*fsten0(3,4)<0) then
	      
	  HJ=HJ*smoothed_sign(fsten0,dx)
	else
	  HJ=HJ*Ssign(cc%ls0)
	end if
	
	!time evolution
	if(n==1) then
	  cc%ls1=cc%ls-HJ*dt_RI
	else
	  cc%ls2=cc%ls-HJ*dt_RI
	end if
	
      end do
   
      !swap ls and ls1
      if(n==1) then
	do i=1,size(cells)
	  templs=cells(i)%ls
	  cells(i)%ls=cells(i)%ls1
	  cells(i)%ls1=templs
	end do
	
	call levelset_transfer(cells,interf,2)
      end if
	   
      !add time n+2 and n
      if(n==2) then
	do i=1,size(cells)
	  cells(i)%ls = 0.5*(cells(i)%ls2 + cells(i)%ls1)
	end do
	
	call levelset_transfer(cells,interf,2)
      end if
	     
    end do
      
    !calculate quit flag of the while loop
    if(nit==1) then 
      it_err0_all=it_err_all
    end if
      
    it_err=0
    ct_bd=0
    do i=1,size(cells)
      if(dabs(cells(i)%ls1) .le. bandsize/2) then
	it_err=it_err+dabs(cells(i)%ls-cells(i)%ls1)
	ct_bd=ct_bd+1
      end if
    end do
    
    call mpi_allreduce(it_err, it_err_all, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_world, ierr)
    call mpi_allreduce(ct_bd, ct_bd_all, 1, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
    
    if(ct_bd_all /= 0) then
      it_err_all=it_err_all/ct_bd_all
    else
      it_err_all=0
      print*, 'interface is not in calculation domain'
    end if
      
    if(it_err_all/dx/dt_RI < 0.1) then
      qt_flag=1
    end if
     
    if(id==0) print*, 'RI', nit, ' error=', it_err_all
    nit=nit+1
    
    !if can not convergence
    if(mod(nit,10)==0) then
      if(it_err_all>it_err0_all .or. isnan(it_err_all)) then
	!not advisable to do re-initialization anymore
	if(id==0) print*,'re-initialization failed'
	do i=1,size(cells)
	  cells(i)%ls=cells(i)%ls0
	end do
	call levelset_transfer(cells,interf,2)
	goto 200
      end if	
    end if

    !force to quit after 100 iteration
    if(nit==100) then
      do i=1,size(cells)
	cells(i)%ls=cells(i)%ls0
      end do
      call levelset_transfer(cells,interf,2)
      if(id==0) print*,'re-initialization failed'
      goto 200
    end if

  end do

  200 continue
  if(id==0) print*, 'reinitialization finished'
  
end subroutine levelset_reinitialization

subroutine deformation_error(cells)
  implicit none
  type(cell), pointer :: cells(:), cc, nc
  integer :: i, ct, ct_all, ind, ierr, id, np
  real(rfp) :: ave_err, max_err, err, dx, true_ls, x, y, z
  real(rfp) :: ave_err_all, max_err_all
  
  call mpi_comm_rank(mpi_comm_world, id, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)
  
  ave_err=0.
  max_err=0.
  ct=0
  
  cc => cells(1)
  ind = cc%inb(1)
  nc => cc%scell(ind)%to_cell
  dx = cc%centp(1) - nc%centp(1)
  
  do i=1,size(cells)
    if(dabs(cells(i)%ls) .le. 5*dx) then
      x=cells(i)%centp(1)
      y=cells(i)%centp(2)
      z=cells(i)%centp(3)
      true_ls = sqrt((x-0.35_rfp)**2+(y-0.35_rfp)**2+(z-0.35_rfp)**2)-0.15_rfp
      err = dabs(true_ls-cells(i)%ls)
      if(err>max_err) max_err=err
      ave_err=ave_err+err
      ct=ct+1
    end if
  end do
  
  call mpi_allreduce(ave_err, ave_err_all, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(max_err, max_err_all, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_world, ierr)
  call mpi_allreduce(ct, ct_all, 1, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
  
  ave_err_all=ave_err_all/ct_all
  
  if(id == 0) then
    open(78, file='./output/deformation_error.dat')
    write(78,*) 'average error=',ave_err_all/dx
    write(78,*) 'maximum error=',max_err_all/dx     
    close(78)
  end if
  call mpi_barrier(mpi_comm_world, ierr)
  
end subroutine deformation_error

function first_order_volume(cells) result(V_all)
  implicit none
  type(cell), pointer :: cells(:), cc, nc
  integer :: i, ind, ierr, id, np
  real(rfp) :: dx, eps, dV, V, V_all
  
  call mpi_comm_rank(mpi_comm_world, id, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)
  
  cc => cells(1)
  ind = cc%inb(1)
  nc => cc%scell(ind)%to_cell
  dx = cc%centp(1) - nc%centp(1)
  eps = 1.5*dx
  dV = dx**3
  V=0.
  
  do i=1,size(cells)
    if(dabs(cells(i)%ls) .le. 6*eps) then
      V=V+(1-heaviside(cells(i)%ls,eps))*dV
    end if
  end do
    
  call mpi_allreduce(V, V_all, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_world, ierr)
  call mpi_barrier(mpi_comm_world, ierr)
end function first_order_volume

!debug
subroutine levelset_check(cells,faces,interf)
  implicit none
  type(itf) :: interf
  type(cell), pointer :: cells(:), current_cell, more_cell, left_cell, right_cell
  type(face), pointer :: faces(:), current_face
  integer :: i,j,more_index,ierr,id
  call mpi_comm_rank(mpi_comm_world, id, ierr)

!   if(id == 4) then
!     open(77, file='./output/interface.dat')
!     write(77,*) interf%nitf
!     write(77,*) interf%nn, interf%nns
!     write(77,*) interf%nc, interf%ncs
!     do i=1,size(interf%np)
!       write(77,*) interf%np(i), interf%nnp(i)
!     end do
!       
!     do i=1,size(interf%cp)
!       write(77,*) interf%cp(i), interf%ncp(i)
!     end do
!     
!     do i=1,size(interf%pcell)
!       write(77,*) i, interf%pcell(i)%ls, interf%pcell(i)%morex, interf%pcell(i)%morels
!     end do
!       
!     close(77)
!   end if
   
  call mpi_barrier(mpi_comm_world, ierr)
end subroutine levelset_check  
!debug

subroutine levelset_thickness(cells)
  implicit none
  type(cell), pointer :: cells(:), cc
  real(rfp) :: vol, dl, fi0, fi1, temp
  integer, pointer :: c2n(:)
  integer :: i, j
  dl = 1.e-1
  do i = 1, size(cells)
    cc => cells(i)
    c2n => cc % c2n
    temp = sqrt(cc % centp(2)**2) !to find the minimum grid size
    if (dl > 0 .and. dl > temp) then
      dl = temp
    endif
  end do
  ls_xi = dl * two * ls_thickness ! two time grid size
  print *, 'level set thickness =', ls_xi
end subroutine levelset_thickness

function minmod(a,b) result(r)
    implicit none  
    real(rfp)::a,b,r
    
    if(dabs(a) .le. dabs(b) .and. a*b>0) then
      r=a
    else if(dabs(a) > dabs(b) .and. a*b>0) then
      r=b
    else
      r=0
    end if
end function minmod

function cubic_intp(xs,ps) result(xA)
    implicit none
    real(rfp) ::xs(4),ps(4),xA,x1,x2,y1,y2,k1,k2,t,a,b
    
    x1=ps(2)
    x2=ps(3)
    !x1*x2<0
    
    y1=xs(2)
    y2=xs(3)
    k1=(xs(3)-xs(1))/(ps(3)-ps(1))
    k2=(xs(4)-xs(2))/(ps(4)-ps(2))
    
    !x=0
    t=-x1/(x2-x1)
    a=k1*(x2-x1)-(y2-y1)
    b=-k2*(x2-x1)+(y2-y1)
    xA=(1-t)*y1+t*y2+t*(1-t)*(a*(1-t)+b*t)
    
end function cubic_intp

function smoothed_sign(ps0,dx) result(s)
    implicit none
    real(rfp) :: ps0(3,5),dx,s
    real(rfp) :: dpx,dpy,dpz,eps

    eps=0.1*dx
    
    dpx=(ps0(1,4)-ps0(1,2))/2.
    dpy=(ps0(2,4)-ps0(2,2))/2.
    dpz=(ps0(3,4)-ps0(3,2))/2.
    s=ps0(3,3)/sqrt(dpx**2+dpy**2+dpz**2)
    
end function smoothed_sign

!sign function
function Ssign(x) result(r)
    implicit none
    real(rfp) :: x,r
    if(x>0) then
      r=1.
    else if(x<0) then
      r=-1.
    else
      r=0.
    end if
end function Ssign

subroutine flip_velocity(cells)
  implicit none
  type(cell), pointer :: cells(:)
  integer :: i
  
  do i=1,size(cells)
    cells(i)%qv%v(imb:ime)=-cells(i)%qv%v(imb:ime)
  end do
end subroutine flip_velocity

subroutine update_levelset(cells)
  implicit none
  type(cell), pointer :: cells(:)
  integer :: i
  
  do i=1,size(cells)
    cells(i)%lsn=cells(i)%ls
  end do
end subroutine update_levelset

end MODULE dfd_levelset
