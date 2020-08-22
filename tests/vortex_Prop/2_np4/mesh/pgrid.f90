program pgrid
   implicit none
   integer, parameter :: rfp  = selected_real_kind(8)  ! resolution of floating point
   type node
      real(kind=rfp),pointer::xyz(:)   ! nodal coordinates
      !   integer,pointer::n2n(:)
      integer::nw  ! temp var
   end type node

   type cell
      integer,pointer::c2n(:) ! point to nodes
      !   integer,pointer::c2c(:) ! point to nodes
      integer::c2c(6) ! point to cells
      integer::nw  ! temp var
      integer::ip
      integer::itype  ! mutli-domain
   end type cell

   type face
      integer::il,ir  !  indices of left and right cell
      integer,pointer::f2n(:) ! pointer to nodes
      integer::nw  ! temp var
      integer::key(3)  ! idenfity unique face
   end type face

   type itf
      integer::n,m
      integer,pointer::cl(:)
   end type

   type(cell),pointer :: cells(:)
   type(face),pointer :: faces(:),new_faces(:),cf
   type(face)::fswap
   type(node),pointer :: nodes(:)
   type(itf), pointer :: mypinterface(:,:)
   type(itf), pointer :: myinterface(:,:)
   !  
   character(len=80)::gridfile,newfile,line
   integer::istatus,nc,f2n(50),c2n(50),nnodes,ncells,nfaces,itype
   integer :: i,ii,j,jj,k,l,fb,im,edgecut,k1,k2,k3,nbsets,nperiod,lab_period,nperiod_p
   integer ::nparts,nct,ntf,nn,nf,ipl,ipr,ndim,nsend,nrecv,intf
   integer ::nvc,vc(50),vnparts(0:50),vncells
   integer,allocatable::epart(:),npart(:),adj(:),n2n(:),esize(:),ind(:)

   ! Xuxiao
   ! for division in all X, Y, and Z direction 
   integer ::npartsX,npartsY,npartsZ,ngeom,write_format
   integer,allocatable::epartX(:),epartY(:),epartZ(:),npartX(:),npartY(:),npartZ(:)

   integer,allocatable::in(:,:,:),ic(:,:,:),ibc(:,:,:),vepart(:)
   integer::iftm,ijk(3),nijk(3),ni,nj,nk,f234(3),fn234(3),max_c2c,nedge,idum
   real(rfp)::tdum,x,y,z
   real(rfp),allocatable::xl(:),fac(:)
   character(len=80)::arg1, arg2, arg3
   ! for multiple periodic boundary condition
   integer :: npbc, ipbc
   integer, allocatable :: npbcs(:), lpbcs(:)

   write_format = 0   ! 0 = binary; 1 = unformatted
   npbc = 2
   ipbc = 0
   if(npbc /= 0) allocate(npbcs(npbc), lpbcs(npbc))

   print *,'input grid file'
   print *,'input file format'
   print *,' 0 = Uniform'
   print *,' 1 = Structure'
   print *,' 2 = Unstructure hybrid (GEMS-DFD)'
   print *,' 21 = Unstructure hybrid (GEMS-DFD) with volume index'
   print *,' 3 = VGRID'
   print *,' 4 = CFDRC Mixed format'
   print *,' 5 = GAMBIT neutral file format'
   print *,' 6 = FIELDVIEW unstructured format(3D only***)'
   print *,' 7 = su2 format, created from GMSH'
   print *,'***  if you choose 2D for FV format, I will try to give you z= 0 surface grid'
   ! gridfile = "hybridMesh"
   ! gridfile = "periodicMesh"
   ! gridfile = "simpleMesh"
   gridfile = "vortexPropMesh"
   ! gridfile = "2D_mesh.dat"
   !read(*,'(a)')gridfile

   print *, 'input file name is ',gridfile

   iftm = 7
   print *, 'input file format is ',iftm
   !read(*,*)iftm

   print *,'input the dimension of problem (1=1D,2=2D,3=3D):'
   ndim = 2
   !read(*,*)ndim
   print *, 'input dimension is ',ndim
   !
   max_c2c = 6    !2 * ndim**3
   !
   ! Read grid data file--------------------------------

   open (2, file = gridfile)
   rewind(2)
   !
   select case(iftm)  
   case(0,1)   ! structured grid(1) or uniform grid(0)
      !  open (2, file = gridfile)
      read(2,*)ni,nj,nk
      if(ndim == 2) nk = 1
      if(ndim == 1) then
         nj = 1;nk = 1
      end if
      nijk(:) = (/ni,nj,nk/)
      ! boundary information
      allocate(ibc(ni,nj,nk))
      call boundary_info(ni,nj,nk,ibc)
      !
      nnodes = product(nijk(:ndim))
      ncells = product(nijk(:ndim) - 1)

      if(ndim == 1)nfaces = ndim * nnodes
      if(ndim == 2)nfaces = ndim * nnodes - sum(nijk(:ndim))
      if(ndim == 3)nfaces = ndim * nnodes + sum(nijk(:ndim)) - &
         2 *(nijk(1)*nijk(2) +nijk(1)*nijk(3) +nijk(3)*nijk(2))
      allocate(nodes(nnodes),cells(ncells),faces(nfaces),stat = istatus)
      if(istatus /= 0 ) print *, 'error in allocate cells and faces'
      ! node's information
      if(iftm == 1) then
         do i = 1, nnodes
            allocate(nodes(i)%xyz(ndim))
         end do
         do i = 1, ndim    
            read(2,*)(nodes(j)%xyz(i),j=1,nnodes)
         end do

      else
         allocate(xl(ndim),fac(ndim))
         read(2,*)xl
         read(2,*)fac
         do i = 1, ndim
            tdum = 0.0_rfp
            do j = 1, nijk(i) - 1
               tdum =  tdum +  fac(i)**j    !1.0_rfp + tdum * fac(i)
            end do
            xl(i) = xl(i) / tdum
         end do
         l = 0
         z = 0.0_rfp
         do k = 1, nk
            y = 0.0_rfp
            do j = 1, nj
               x = 0.0_rfp
               do i = 1, ni
                  l = l + 1
                  allocate(nodes(l)%xyz(ndim))
                  nodes(l)%xyz(:2) = (/x,y/)
                  if(ndim == 3) nodes(l)%xyz(3) = z
                  x = x + fac(1)**i * xl(1)
               end do
               y = y + fac(2)**j * xl(2)
            end do
            z = z + fac(3)**k * xl(3)
         end do
         deallocate(xl,fac)
      end if

      select case(ndim)
      case(1)
         j = 1
         k = 1
         allocate(in(nijk(1),1,1),ic(nijk(1)-1,1,1))
         in = reshape((/(i,i=1,nnodes)/),shape=(/nijk(1),1,1/))
         ic = reshape((/(i,i=1,ncells)/),shape=(/nijk(1)-1,1,1/))
         do i = 1, ni - 1
            l = ic(i,j,k)
            nc = 2
            allocate(cells(l)%c2n(nc))
            c2n(1) = in(i,  j,  k)
            c2n(2) = in(i+1,j,  k)
            cells(l)%c2n = c2n(:nc)
         end do

      case(2)
         k = 1
         allocate(in(nijk(1),nijk(2),1),ic(nijk(1)-1,nijk(2)-1,1))
         in = reshape((/(i,i=1,nnodes)/),shape=(/nijk(1),nijk(2),1/))
         ic = reshape((/(i,i=1,ncells)/),shape=(/nijk(1)-1,nijk(2)-1,1/))
         !  cell's information
         do j = 1, nj -1
            do i = 1, ni - 1
               l = ic(i,j,k)
               nc = 4
               allocate(cells(l)%c2n(nc))
               c2n(1) = in(i,  j,  k)
               c2n(2) = in(i+1,j,  k)
               c2n(3) = in(i+1,j+1,k)
               c2n(4) = in(i,  j+1,k)
               cells(l)%c2n = c2n(:nc)
            end do
         end do
      case(3)
         allocate(in(nijk(1),nijk(2),nijk(3)),ic(nijk(1)-1,nijk(2)-1,nijk(3)-1))
         in = reshape((/(i,i=1,nnodes)/),shape=(/nijk(1),nijk(2),nijk(3)/))
         ic = reshape((/(i,i=1,ncells)/),shape=(/nijk(1)-1,nijk(2)-1,nijk(3)-1/))
         !  cell's information
         do k = 1, nk - 1 
            do j = 1, nj -1
               do i = 1, ni - 1
                  l = ic(i,j,k)
                  nc = 8
                  allocate(cells(l)%c2n(nc))
                  c2n(1) = in(i,  j,  k)
                  c2n(2) = in(i,  j+1,k)
                  c2n(3) = in(i+1,j+1,k)
                  c2n(4) = in(i+1,j,  k)
                  c2n(5) = in(i,  j,  k+1)
                  c2n(6) = in(i,  j+1,k+1)
                  c2n(7) = in(i+1,j+1,k+1)
                  c2n(8) = in(i+1,j,  k+1)
                  cells(l)%c2n = c2n(:nc)
                  !
               end do
            end do
         end do
      end select 
      !
   case(2,21)
      read(2,*)nnodes,ncells,nfaces
      !
      allocate(nodes(nnodes),cells(ncells),faces(nfaces),stat = istatus)
      if(istatus /= 0 ) print *, 'error in allocate cells and faces'
      !   
      do i = 1, nnodes
         allocate(nodes(i)%xyz(ndim))
         read(2,*)nodes(i)%xyz 
      end do
      do i = 1, ncells
         if(iftm == 2) then
            read(2,*)nc,c2n(1:nc) 
            itype = 0
         else
            read(2,*)nc,c2n(1:nc),itype
         end if
         !     if(nc /= nct) then 
         !     print *, 'so far, can not handle hybrid grid',nc,nct
         !     stop
         !     end if
         allocate(cells(i)%c2n(nc))
         cells(i)%c2n = c2n(:nc)
         cells(i)%itype = itype
      end do
      !
   case(3)  ! unstructured grid output by vgrid and neal's preprocessor
      l = len_trim(gridfile)
      newfile = gridfile(:l)//'.cogsg'   
      open (2, file = newfile,form='unformatted')
      newfile = gridfile(:l)//'.iface'   
      open (3, file = newfile,form='unformatted')
      read(2)idum,ncells,nnodes,idum,idum,idum,tdum
      read(3)nfaces,fb
      rewind(2)
      rewind(3)
      print *,nnodes,ncells,nfaces
      allocate(nodes(nnodes),cells(ncells),faces(nfaces),stat = istatus)
      if(istatus /= 0 ) print *, 'error in allocate cells and faces'

      do i = 1, nnodes
         allocate(nodes(i)%xyz(ndim),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate gradient'
      end do

      nc = 4
      do i = 1, ncells
         allocate(cells(i)%c2n(nc),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate gradient'
      end do
      !
      read(2)idum,ncells,nnodes,idum,idum,idum,tdum,&
         ((cells(i)%c2n(j),i=1,ncells),j=1,nc)
      read(2)((nodes(i)%xyz(j),i=1,nnodes),j=1,ndim)
      close(2)
      !
   case(4)  ! CFDRC-GEOM's mixed unstructured grid
      do i = 1, 5
         read(2,*)    ! skip first 5 headlines
      end do
      read(2,*)nnodes
      read(2,*)f234   ! boundary faces two, thress or four nodes's face
      fb = sum(f234)
      nfaces = fb
      fn234 = 0
      read(2,*)ijk   ! interface faces
      nfaces = nfaces + sum(ijk) 
      fn234 = fn234 + ijk
      read(2,*)ijk   ! interior faces
      fn234 = fn234 + ijk
      nfaces = nfaces + sum(ijk) 
      read(2,*)ijk,nijk  !quad, tri, hex,tet,prism,and pyramids
      if(ndim == 2) then
         ncells = sum(ijk(:2))
      else
         ncells = ijk(3) + sum(nijk)
      end if
      allocate(nodes(nnodes),cells(ncells),faces(nfaces),stat = istatus)
      if(istatus /= 0 ) print *, 'error in allocate cells and faces'
      read(2,*)ni,nj,nk   ! read number of boundary, interface and cell patches
      read(2,*)  ! skip node comment's line
      do i = 1, nnodes
         allocate(nodes(i)%xyz(ndim))
         read(2,*)nodes(i)%xyz
      end do
      read(2,*)  ! skip cell comment's line
      nc = 4
      i = 0
      do l = 1, ijk(1) ! quad's cells
         i = i + 1
         allocate(cells(i)%c2n(nc))
         read(2,*)cells(i)%c2n
      end do
      !
      nc = 3          ! tri's cells
      do l = 1, ijk(2)
         i = i + 1
         allocate(cells(i)%c2n(nc))
         read(2,*)cells(i)%c2n
      end do
      !
      nc = 8          ! Hex's cells
      do l = 1, ijk(3)
         i = i + 1
         allocate(cells(i)%c2n(nc))
         read(2,*)cells(i)%c2n
      end do
      !
      nc = 4          ! tet's cells
      do l = 1, nijk(1)
         i = i + 1
         allocate(cells(i)%c2n(nc))
         read(2,*)cells(i)%c2n
      end do

      nc = 6          ! Prism's cells
      do l = 1, nijk(2)
         i = i + 1
         allocate(cells(i)%c2n(nc))
         read(2,*)cells(i)%c2n
      end do

      if(i /= ncells) print *,'error in grid files',i,ncells
      !
      ! GAMBIT Neutral File Format ----------------
   case(5)
      call gambit2dfd
   case(6)
      call fieldview2dfd
   case(7)
      call su2dfd
   end select
   ! 
   print *,'nnodes =',nnodes
   print *,'ncells =',ncells
   print *,'nfaces =',nfaces

   print *,'How many partitions? (if input <=1, just for convert to GEMS-DFD format)'
   nparts = 4
   ! call getarg(1, arg1)
   ! read(arg1, *) nparts
   ! read(*,*)nparts
   print *, 'partition number is ',nparts

   if (nparts < 0) then  !for plot only
      open(12,file='geom.plt')
      write(12,'("variables = ")')
      do j = 1,ndim
         write(12,'(1x,1hx,i1)')j
      end do

      ! need to be modified for your grid
      !
      ! write(*,*)'which type mesh 0=quadrilateral 1 = Hex?'
      ! read(*,*)k
      k = 2

      if(k == 0) then
         write(12,*)'zone t="',i,'",n=',size(nodes),',e=',size(cells),  &
            ',f=feblock,et=quadrilateral'
      else if(k == 1) then
         write(12,*)'zone t="',i,'",n=',size(nodes),',e=',size(cells),  &
            ',f=feblock,et=brick'
      else if(k == 2) then
         write(12,*)'zone t="',i,'",n=',size(nodes),',e=',size(cells),  &
            ',f=feblock,et=tetrahedron'
      end if
      !
      do j = 1, ndim
         write(12,20)(nodes(i)%xyz(j),i=1,size(nodes))
      end do

      do j = 1, ncells
         l = size(cells(j)%c2n)
         c2n(:l) = cells(j)%c2n
         ! this is to replenish the artifial nodes for tecplot
         select case(k)
         case(0)
            if(l == 3)  c2n(4) = c2n(3)
            l = 4
         case(1)
            if(l == 4) then
               c2n(4) = c2n(3)
               c2n(5:8) = c2n(4)
            else if(l == 5) then
               c2n(5:8) = c2n(5)
            else if(l == 6) then
               c2n(8) = c2n(6)
               c2n(5:7) = c2n(4:6)
               c2n(4) = c2n(3)
            end if
            l = 8
         end select
         write(12,10)c2n(:l)
      end do
      close(12)
      stop
   end if
   !
   ! End of Ploting


   !Read face information------------------------------------------------
   !
   select case(iftm)
   case(0,1)
      select case(ndim)
      case(1)
         l = 0; fb = 0; j=1; k=1
         do i=1,ni
            l=l+1
            if(i==1) then
               faces(l)%il = ic(i,j,k)
               faces(l)%ir = - max(ibc(i+1,j,k),ibc(i,j,k))
               fb=fb+1
            else if(i==ni) then
               faces(l)%il = ic(i-1,j,k)
               faces(l)%ir = - max(ibc(i-1,j,k),ibc(i,j,k))
               fb = fb + 1
            else
               faces(l)%il  = ic(i-1,j,k)
               faces(l)%ir  = ic(i,j,k)
            end if

            nc=1
            allocate(faces(l)%f2n(nc))
            faces(l)%f2n(1) = in(i,j,k)
         end do

      case(2)
         !
         l = 0; fb = 0
         do j = 1, nj
            do i = 1, ni - 1
               l = l + 1
               if(j == 1) then
                  faces(l)%il = ic(i,j,k)
                  faces(l)%ir = - max(ibc(i+1,j,k),ibc(i,j,k))
                  fb = fb + 1
               else if(j == nj) then
                  faces(l)%il = ic(i,j-1,k)
                  faces(l)%ir = - max(ibc(i+1,j,k),ibc(i,j,k))
                  fb = fb + 1
               else
                  faces(l)%il  = ic(i,j-1,k)
                  faces(l)%ir  = ic(i,j,  k)
               end if
               nc = 2
               allocate(faces(l)%f2n(nc))
               faces(l)%f2n(1) = in(i,  j,k)
               faces(l)%f2n(2) = in(i+1,j,k)
            end do
         end do

         do i = 1, ni
            do j = 1, nj - 1
               l = l +1
               if(i == 1) then
                  fb = fb + 1
                  faces(l)%il = ic(i,j,k)
                  faces(l)%ir = -max(ibc(i,j,k),ibc(i,j+1,k))
               else if(i == ni) then
                  fb = fb + 1
                  faces(l)%il = ic(i-1,j,k)
                  faces(l)%ir = -max(ibc(i,j,k),ibc(i,j+1,k))
               else
                  faces(l)%il  = ic(i-1,j,k)
                  faces(l)%ir  = ic(i,  j,k)
               end if
               nc = 2
               allocate(faces(l)%f2n(nc))
               faces(l)%f2n(1) = in(i,j,  k)
               faces(l)%f2n(2) = in(i,j+1,k)
            end do
         end do
      case(3)
         !
         l = 0;fb=0

         ! i face
         do i = 1, ni
            do k = 1, nk-1
               do j = 1, nj - 1
                  l = l +1
                  if(i == 1) then
                     fb = fb + 1
                     faces(l)%il = ic(i,j,k)
                     faces(l)%ir = -maxval(ibc(i:i,j:j+1,k:k+1))
                  else if(i == ni) then
                     faces(l)%il = ic(i-1,j,k)
                     faces(l)%ir = -maxval(ibc(i:i,j:j+1,k:k+1))
                  else
                     faces(l)%il  = ic(i-1,j,k)
                     faces(l)%ir = ic(i,  j,k)
                  end if
                  nc = 4
                  allocate(faces(l)%f2n(nc))
                  faces(l)%f2n(1) = in(i,j,  k)
                  faces(l)%f2n(2) = in(i,j+1,k)
                  faces(l)%f2n(3) = in(i,j+1,k+1)
                  faces(l)%f2n(4) = in(i,j,  k+1)
                  if(faces(l)%ir < 0) fb = fb + 1
               end do
            end do
         end do
         ! j face
         do j = 1, nj
            do k = 1, nk - 1
               do i = 1, ni - 1
                  l = l + 1
                  if(j == 1) then
                     faces(l)%il = ic(i,j,k)
                     faces(l)%ir = -maxval(ibc(i:i+1,j:j,k:k+1))
                  else if(j == nj) then
                     faces(l)%il = ic(i,j-1,k)
                     faces(l)%ir = -maxval(ibc(i:i+1,j:j,k:k+1))
                  else
                     faces(l)%il  = ic(i,j-1,k)
                     faces(l)%ir = ic(i,j,k  )
                  end if
                  nc = 4
                  allocate(faces(l)%f2n(nc))
                  faces(l)%f2n(1) = in(i,  j,  k)
                  faces(l)%f2n(2) = in(i+1,j,  k)
                  faces(l)%f2n(3) = in(i+1,j,  k+1)
                  faces(l)%f2n(4) = in(i,  j,  k+1)
                  if(faces(l)%ir< 0) fb = fb + 1
               end do
            end do
         end do

         ! k face
         do k = 1,nk
            do j = 1, nj - 1
               do i = 1, ni - 1
                  l = l + 1
                  if(k == 1) then
                     faces(l)%il = ic(i,j,k)
                     faces(l)%ir = -maxval(ibc(i:i+1,j:j+1,k:k))
                  else if(k == nk) then
                     faces(l)%il = ic(i,j,k-1)
                     faces(l)%ir = -maxval(ibc(i:i+1,j:j+1,k:k))
                  else
                     faces(l)%il  = ic(i,j,k-1)
                     faces(l)%ir  = ic(i,j,k  )
                  end if
                  nc = 4
                  allocate(faces(l)%f2n(nc))
                  faces(l)%f2n(1) = in(i,  j,  k)
                  faces(l)%f2n(2) = in(i+1,j,  k)
                  faces(l)%f2n(3) = in(i+1,j+1,k)
                  faces(l)%f2n(4) = in(i,  j+1,k)
                  if(faces(l)%ir < 0) fb = fb + 1
               end do
            end do
         end do
      end select
      deallocate(in,ic,ibc)

   case(2,21)
      fb = 0
      do i = 1, nfaces
         cf => faces(i)
         read(2,*)cf%il,cf%ir,nc,f2n(1:nc) 
         allocate(cf%f2n(nc))
         cf%f2n = f2n(1:nc)
         if(cf%ir < 0) fb = fb + 1
      end do

   case(3)
      nc = 3
      do i = 1, nfaces
         allocate(faces(i)%f2n(nc),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate face array'
      end do
      read(3)nfaces,fb,(faces(i)%il,i=1,nfaces),&
         (faces(i)%ir,i=1,nfaces),&
         ((faces(i)%f2n(j),i=1,nfaces),j=1,nc)
      close(3)
      k = maxval(faces(:fb)%ir)
      allocate(ind(k),stat=istatus)  !k is max patch number
      print *,'max patch number is',k

      newfile = gridfile(:l)//'.mapbc'
      open(2,file=newfile)
      do i = 1, 4
         read(2,'(a)')   ! skip first 4 lines
      end do
      !
      !   do while(.true.)
      do i = 1, k
         read(2,*)j,l
         ind(j) = neal2dli(l)
      end do
      !
      do i = 1, fb
         faces(i)%ir = -ind(faces(i)%ir)
      end do
      deallocate(ind)
      close(2)

   case(4)
      read(2,*)  ! skip face comment's line

      i = 0
      k = 0
      nc = 2
      do l = 1, f234(1) !fb              !boundary faces
         i = i + 1
         allocate(faces(i)%f2n(nc),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate face array'
         read(2,*)faces(i)%f2n,faces(i)%il,faces(i)%ir,itype
         faces(i)%ir = itype
         k = max(k,itype)
      end do
      !
      do l = 1, fn234(1)  !fb+1,nfaces  !line face 
         i = i + 1
         allocate(faces(i)%f2n(nc),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate face array'
         read(2,*)faces(i)%f2n,faces(i)%il,faces(i)%ir
      end do

      nc = 3
      do l = 1, f234(2) !fb              !boundary faces
         i = i + 1
         allocate(faces(i)%f2n(nc),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate face array'
         read(2,*)faces(i)%f2n,faces(i)%il,faces(i)%ir,itype
         faces(i)%ir = itype
         k = max(k,itype)
      end do

      do l = 1, fn234(2)  !fb+1,nfaces  !tri's faces
         i = i + 1
         allocate(faces(i)%f2n(nc),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate face array'
         read(2,*)faces(i)%f2n,faces(i)%il,faces(i)%ir
      end do

      nc = 4
      do l = 1, f234(3)  !fb             !boundary faces
         i = i + 1
         allocate(faces(i)%f2n(nc),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate face array'
         read(2,*)faces(i)%f2n,faces(i)%il,faces(i)%ir,itype
         faces(i)%ir = itype
         k = max(k,itype)
      end do
      !
      do l = 1, fn234(3)  !fb+1,nfaces  !quad faces
         i = i + 1
         allocate(faces(i)%f2n(nc),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate face array'
         read(2,*)faces(i)%f2n,faces(i)%il,faces(i)%ir
      end do
      !
      read(2,*)   ! skip boundary condition comment's line
      allocate(ind(k),stat=istatus)  !k is max key number
      print *,'max key number is',k
      if(istatus /= 0 ) print *, 'error in allocate in array'
      !      ni = ni + nj
      do i = 1, ni
         read(2,'(a)')line
         if(nj > k) cycle
         call lookup_key(line,ind)  !
      end do
      !      
      do i = 1, fb
         faces(i)%ir = -ind(faces(i)%ir)   ! convert bounary condition
      end do
      deallocate(ind)
   case(5, 7)   ! FE mesh: gambit and su2
      call fe2dfd_face
   end select
   !
   ! do swap: move the boundary faces to the front of faces array
   fb = 0
   do i = 1, nfaces
      if(faces(i)%ir <= 0) then
         fb = fb + 1
         if (i == fb) cycle
         fswap = faces(fb)
         faces(fb) = faces(i)
         faces(i) = fswap
      end if
   end do
   if(iftm == 5) call read_gambit_boundary
   if(iftm == 7) call read_su_boundary
   print *,'number of boundary faces=',fb
   !
   nperiod = -1
   close(2)

   ! do i = 1, nfaces
   !    print *, i, faces(i) % il, faces(i) % ir
   ! end do
   do ii = 1, npbc
      ipbc = ipbc + 1
      call treatment_periodic_bc
   end do
   ! do i = 1, nfaces
   !    print *, i, faces(i) % il, faces(i) % ir
   ! end do

   if(nparts <= 1) then
      if(write_format == 0) then
         open(2,file=gridfile(:len_trim(gridfile))//'.1',form='binary')
      else
         open(2,file=gridfile(:len_trim(gridfile))//'.1.annotate')
      end if
      print *,gridfile(:len_trim(gridfile))//'.1 is generated'
      if(write_format == 0) then
         write(2)nnodes,ncells,nfaces,nparts
      else
         write(2,*)nnodes,ncells,nfaces,nparts
      end if

      do k = 1, nnodes
         if(write_format == 0) then
            write(2)nodes(k)%xyz
         else
            write(2,*)nodes(k)%xyz
         end if
      end do
      !
      do j = 1, ncells
         l = size(cells(j)%c2n)
         if(write_format == 0) then
            write(2)l,cells(j)%c2n,cells(j)%itype
         else
            write(2,*)l,cells(j)%c2n,cells(j)%itype
         end if
      end do
      !
      do j = 1, nfaces
         if(write_format == 0) then
            write(2)faces(j)%il,faces(j)%ir,size(faces(j)%f2n),faces(j)%f2n
         else
            write(2,*)faces(j)%il,faces(j)%ir,size(faces(j)%f2n),faces(j)%f2n
         end if
      end do
      !
      if(nperiod >= 0) then
         if(write_format == 0) then
            write(2)nperiod,lab_period
            write(2)1,fb - nperiod
            write(2)0,fb-nperiod,-faces(nperiod+1:fb)%ir
            write(2)1,fb - nperiod
            write(2)0,fb-nperiod,(i,i=nperiod+1,fb)
         else
            write(2,*)nperiod,lab_period
            write(2,*)1,fb - nperiod
            write(2,*)0,fb-nperiod,-faces(nperiod+1:fb)%ir
            write(2,*)1,fb - nperiod
            write(2,*)0,fb-nperiod,(i,i=nperiod+1,fb)
         end if
      end if
      !
      !
      stop
   end if
   !
   if(nparts > 1) then

      print *,'Which partitioning algorithm(0=user defined,1=metis edgecut,2=metis Comvolume,3=x,4=y,5=z,6=multi-domain,7=xyz)'
      im = 7
      print *, 'partitioning algorithm is ',im
      !read(*,*)im

      !
      ! Do partitioning ---------------------------------------
      !
      print *,'Begin partition............'   
      open(12,file='gems.par')
      write(12,*)'id1    id2    itype    file_name'
      select case(im)
      case(0)
         allocate(epart(ncells),stat=istatus)
         epart = cells%itype
         call check_parts
         write(*,*)'you have',nvc,' partitions with itype=',vc(:nvc)
         do i = 1, nvc
            write(12,*)i-1,i-1,vc(i),'  "gems.'//trim(i2s(i))//'.inp"'
            do j = 1, size(cells)
               if(cells(j)%itype /= vc(i)) cycle
               epart(j) = i
            end do
         end do
         nparts = nvc
      case(1,2)
      ! Missing METIS, do not use
         write(12,*)0,nparts-1,1,'  "gems.inp"'
         do i = 1, ncells
            !     allocate(cells(i)%c2c(max_c2c))
            cells(i)%nw = 0
         end do
         do i = 1, nfaces
            call face_2_c2c(faces(i)%il,faces(i)%ir,cells)
         end do
         nedge = sum(cells%nw)
         print *,nedge,nfaces

         allocate(adj(ncells+1),n2n(nedge),epart(ncells),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate face array'
         adj(1)=1
         do i = 2, ncells+1
            adj(i) = adj(i-1) + cells(i-1)%nw
            n2n(adj(i-1):adj(i)-1) = cells(i-1)%c2c(:cells(i-1)%nw)
         end do
         if(im==1) then
            ! call METIS_PartGraphKway(ncells,adj,n2n,0,0,0,1,nparts,0,edgecut,epart)
            print *,'total edges cut:',edgecut
         else
            ! call METIS_PartGraphvKway(ncells,adj,n2n,0,0,0,1,nparts,0,edgecut,epart)
            print *,'total communication volume:',edgecut
         end if
         deallocate(adj,n2n)
      case(3,4,5)
         write(12,*)0,nparts-1,1,'  "gems.inp"'
         allocate(npart(size(nodes)))
         allocate(epart(ncells),stat=istatus)

         call Geom_Partxyz(npart, nparts, im-2)
         ! call Geom_Partijk(npart,nparts,im-2)

         do i = 1, ncells
            epart(i) = minval(npart(cells(i)%c2n))
         end do
         deallocate(npart)
      case(6)
         allocate(epart(ncells),stat=istatus)
         call check_parts
         print *, ' You have',nvc, 'domains with itype=',vc(:nvc)
         print *, ' please input number of partitions for each domain'
         read(*,*)vnparts(1:nvc)
         vnparts(0) = 0
         do k = 1, nvc
            itype = vc(k)
            vncells = 0
            do i = 1, ncells
               cells(i)%nw = 0
               if(cells(i)%itype /= itype) cycle
               vncells = vncells + 1
               cells(i)%ip = vncells
            end do

            do i = 1, nfaces
               call vface_2_c2c(faces(i)%il,faces(i)%ir,cells,itype)
            end do

            nedge = sum(cells%nw)
            print *,k,nedge,nfaces,vncells
            allocate(adj(vncells+1),n2n(nedge),vepart(vncells),stat=istatus)
            if(istatus /= 0 ) print *, 'error in allocate face array'
            adj(1)=1
            j = 1
            do i = 1, ncells
               if(cells(i)%itype /= itype) cycle
               j = j + 1
               adj(j) = adj(j-1) + cells(i)%nw
               n2n(adj(j-1):adj(j)-1) = cells(i)%c2c(:cells(i)%nw)
            end do
            !     if(im==1) then
            !     call METIS_PartGraphKway(ncells,adj,n2n,0,0,0,1,&
            !                            nparts,0,edgecut,epart)
            !     print *,'total edges cut:',edgecut
            !     else
            if(vnparts(k) == 1) then
               vepart = 1
            else
               !      call METIS_PartGraphvKway(vncells,adj,n2n,0,0,0,1,&
               !			    vnparts(k),0,edgecut,vepart)
            end if
            !      print *,'total communication volume:',edgecut
            !    end if
            j = 0
            do i = 1, ncells
               if(cells(i)%itype /= itype) cycle
               j = j + 1
               epart(i) = vepart(j) + vnparts(k-1)
            end do
            vnparts(k) = vnparts(k) + vnparts(k-1)

            deallocate(adj,n2n,vepart)
         end do
         nparts = vnparts(nvc)
         do i = 1, nvc
            write(12,*)vnparts(i-1),vnparts(i)-1,vc(i),'  "gems.'//trim(i2s(i))//'.inp"'
         end do

      case(7)
         if(nparts == 1) then
            npartsX = 1; npartsY = 1; npartsZ = 1
         else if(nparts == 2) then
            npartsX = 2; npartsY = 1; npartsZ = 1
         else if(nparts == 4) then
            npartsX = 2; npartsY = 2; npartsZ = 1
         else if(nparts == 8) then
            npartsX = 2; npartsY = 2; npartsZ = 2
         else if(nparts == 16) then
            npartsX = 2; npartsY = 2; npartsZ = 4
         else if(nparts == 20) then
            npartsX = 2; npartsY = 5; npartsZ = 2
         else if(nparts == 32) then
            npartsX = 4; npartsY = 4; npartsZ = 2
         else if(nparts == 40) then
            npartsX = 4; npartsY = 5; npartsZ = 2
         else if(nparts == 80) then
            npartsX = 4; npartsY = 5; npartsZ = 4
         else
            print *, 'partioning the domain in all X, Y, and Z directions'
            print *, 'How many partitions along X direction'
            read(*, *) npartsX
            print *, 'How many partitions along Y direction'
            read(*, *) npartsY
            print *, 'How many partitions along Z direction'
            read(*, *) npartsZ
         end if

         nparts = npartsX * npartsY * npartsZ

         write(12, *)0, nparts - 1, 1, '  "gems.inp"'
         allocate(npartX(size(nodes)), npartY(size(nodes)), npartZ(size(nodes)))
         allocate(epartX(ncells), stat = istatus)
         allocate(epartY(ncells), stat = istatus)
         allocate(epartZ(ncells), stat = istatus)
         allocate(epart(ncells), stat = istatus)

         ! call Geom_Partijk(npartX, npartsX, 1)
         ! call Geom_Partijk(npartY, npartsY, 2)
         ! call Geom_Partijk(npartZ, npartsZ, 3)
         call Geom_Partxyz(npartX, npartsX, 1)
         call Geom_Partxyz(npartY, npartsY, 2)
         call Geom_Partxyz(npartz, npartsZ, 3)

         do i = 1, ncells
            epartX(i) = minval(npartX(cells(i) % c2n))
            epartY(i) = minval(npartY(cells(i) % c2n))
            epartZ(i) = minval(npartZ(cells(i) % c2n))
         end do

         !*********** check how many cells are in each partition ***********
         print *, 'X direction'
         !call countcells(npartX, npartsX)
         call countcells(epartX, npartsX)
         print *, 'Y direction'
         !call countcells(npartY, npartsY)
         call countcells(epartY, npartsY)
         print *, 'Z direction'
         !call countcells(npartZ, npartsZ)
         call countcells(epartZ, npartsZ)
         !print *, 'press enter to continue'
         !read(*, *)
         !*********** end of counting cells *****************

         do i = 1, ncells
            epart(i) = (epartX(i)-1)*npartsY*npartsZ + (epartY(i)-1)*npartsZ + epartZ(i)
         end do

         print *, 'number of cells in each partition'
         call countcells(epart, nparts)
         ! print *, 'press enter to continue'
         ! read(*, *)
         deallocate(npartX, npartY, npartZ, epartX, epartY, epartZ)

      end select
      close(12)

   end if
   !

   !Remeshing ----------------------------------------------------------
   !
   print *,'Begin remesh............'   
   allocate(esize(nparts))
   !
   call setup_interface
   !
   if(nperiod >= 0) call setup_periodic_interface

   !
   esize = 0
   do i = 1,ncells
      esize(epart(i)) = esize(epart(i)) + 1
      cells(i)%ip = epart(i)
      !   cells(i)%nw = esize(epart(i))
   end do

   do i = 1, ncells
      print *, i, cells(i) % ip
   end do
   !   
   print *,esize
   open(12,file='dfd.plt')
   !
   do i = 1,nparts
      ! find cells, nodes, faces for this partition
      ! must re-index for the cells, nodes and faces for this partition
      if(write_format == 0) then
         open(2,file=gridfile(:len_trim(gridfile))//'.'//i2s(i),form='binary')
      else
         open(2,file=gridfile(:len_trim(gridfile))//'.'//trim(i2s(i))//'.annotate')
      end if
      
      print *,'Partition:',i
      write(12,'("variables = ")')
      do j = 1,ndim
         write(12,'(1x,1hx,i1)')j
      end do

      nodes%nw = 0
      cells%nw = 0
      nn = 0
      nc = 0
      do j = 1, ncells
         if(epart(j) /= i) cycle
         nc = nc + 1
         cells(j)%nw = nc
         !    if( cells(j)%nw /= nc) print *,'something wrong!!!',cells(j)%nw,nc
         nct = size(cells(j)%c2n)
         do k = 1, nct
            l = cells(j)%c2n(k)
            if(nodes(l)%nw == 0) then
               nn = nn + 1
               nodes(l)%nw = nn
            end if
            !     cells(j)%c2n(k) = nodes(l)%nw
         end do
      end do
      allocate(npart(nn))
      !
      !  Interface
      nsend = sum(myinterface(i,:)%n)
      nrecv = 0
      do j = 1, nparts   
         nf = myinterface(j,i)%n
         !    if(nf <= 0) cycle
         do l = 1, nf
            !     nrecv = nrecv + 1
            k = myinterface(j,i)%cl(l)
            if(cells(k)%nw /= 0) print *,'error in interface'
            cells(k)%nw = nrecv + l
         end do
         nrecv = nrecv + nf + myinterface(j,i)%m
      end do
      !
      do j = 1, nnodes
         if(nodes(j)%nw == 0) cycle
         npart(nodes(j)%nw) = j
      end do
      !
      !  treat face's information
      nf = 0  ! no. of totall faces
      intf = 0  ! no. of partitioning faces
      fb = 0  ! no. of boundary faces
      nperiod_p = 0   ! no. of periodic boundary faces
      do j = 1, nfaces
         cf => faces(j)
         if(cf%ir < 0 ) then
            ipl = cells(cf%il)%ip
            if( ipl /= i) cycle
            fb = fb + 1
            nf = nf + 1
            if(nperiod >= 0.and.j > nperiod) then
               !      cf%ir = -cells(abs(cf%ir))%nw !+ lab_period
               nperiod_p = nperiod_p + 1
               cf%ir = -nperiod_p
               !      cf%ir = - lab_period
               !      cf%ir = -cf%ir
            end if

         else
            !
            ipl = cells(cf%il)%ip
            ipr = cells(cf%ir)%ip
            if(ipl /= i .and. ipr /= i) cycle
            nf = nf + 1
            if(ipl == i .and. ipr == i) cycle
            intf = intf + 1
         end if
      end do   
      print *,'ok',nf,intf
      allocate(new_faces(nf))
      !
      k1 = 0
      k2 = 0
      k3 = 0
      do j = 1, nfaces
         cf => faces(j)
         if(cf%ir < 0 ) then
            ipl = cells(cf%il)%ip
            if( ipl /= i) cycle
            k1 = k1 + 1
            l = intf + k1  ! count from intf for boundary faces
            new_faces(l) = faces(j)
            faces(j)%nw = l
         else
            !
            ipl = cells(cf%il)%ip
            ipr = cells(cf%ir)%ip
            if(ipl /= i .and. ipr /= i) cycle
            if(ipl == i .and. ipr == i) then
               k2 = k2 + 1
               l = intf + fb + k2  ! count from fb + intf for interior faces
               new_faces(l) = faces(j)
            else if(ipl == i) then  ! partitioning faces
               k3 = k3 + 1
               l = k3
               new_faces(l) = faces(j)
            else  ! partitioning faces
               k3 = k3 + 1
               l = k3
               new_faces(l) = faces(j)
               ! make sure the left cell stays inside the current partition
               new_faces(l)%il = faces(j)%ir  
               new_faces(l)%ir = faces(j)%il
            end if
         end if
         !      new_faces(l)%f2n => faces(j)%f2n
      end do   

      !   print *,sum(intf)+k+nf,size(new_faces)
      !   if(i < 10) then
      !   write(newfile,'(a,1h.,i1)')gridfile(:len_trim(gridfile)),i   
      !   else if(i < 100) then
      !   write(newfile,'(a,1h.,i2)')gridfile(:len_trim(gridfile)),i   
      !   else
      !   write(newfile,'(a,1h.,i3)')gridfile(:len_trim(gridfile)),i   
      !   end if

      ! do ii = 1, size(cells)
      !    print *, ii, cells(ii) % nw
      ! end do
      ! print *, '-------------'

      ! do ii = 1, size(nodes)
      !    print *, ii, nodes(ii) % nw
      ! end do
      ! print *, '-------------'

      ! do ii = 1, size(new_faces)
      !    print *, ii, new_faces(ii) % il, new_faces(ii) % ir
      ! end do
      ! print *, '-------------'

      ! do ii = 1, nfaces
      !    print *, ii, faces(ii) % il, faces(ii) % ir
      ! end do
      ! print *, '-------------'

      ! do ii = 1, nparts
      !    do jj = 1, nparts
      !       print *, ii, jj, myinterface(ii,jj) % n, myinterface(ii,jj) % m
      !       print *, myinterface(ii,jj) % cl
      !    end do
      ! end do
      ! print *, '-------------'

      ! do ii = 1, nparts
      !    do jj = 1, nparts
      !       print *, ii, jj, mypinterface(ii,jj) % n, mypinterface(ii,jj) % m
      !       print *, mypinterface(ii,jj) % cl
      !    end do
      ! end do
      ! stop

      ! Just for ploting
      !
      !print *,'-------'
      !   print *,cells(new_faces(:ind(2))%il)%nw 
      !   cells(new_faces(:ind(2))%ir)%nw
      !print *,'-------'

      if(write_format == 0) then
         write(2)nn,esize(i),size(new_faces),nparts
      else
         write(2,*)nn,esize(i),size(new_faces),nparts
      end if
      ngeom = 1
      !
      !
      select case(ngeom)
      case(1)
         write(12,*)'zone t="',i,'",n=',nn,',e=',esize(i), ',f=feblock,et=triangle'
      case(2)
         write(12,*)'zone t="',i,'",n=',nn,',e=',esize(i), ',f=feblock,et=quadrilateral'
      case(3)
         write(12,*)'zone t="',i,'",n=',nn,',e=',esize(i), ',f=feblock,et=tetrahedron'
      case(4)
         write(12,*)'zone t="',i,'",n=',nn,',e=',esize(i), ',f=feblock,et=brick'
      case default
         if(ndim == 2) then
            write(12,*)'zone t="',i,'",n=',nn,',e=',esize(i), ',f=feblock,et=quadrilateral'
         else if(ndim == 3) then
            write(12,*)'zone t="',i,'",n=',nn,',e=',esize(i), ',f=feblock,et=brick'
         else
            print *, 'cannot plot 1D mesh'; stop
         end if
      end select

      do k = 1, nn
         if(write_format == 0) then
            write(2)nodes(npart(k))%xyz
         else
            write(2,*)nodes(npart(k))%xyz
         end if
      end do

      do j = 1, ncells
         if(epart(j) == i) then
            l = size(cells(j)%c2n)
            if(write_format == 0) then
               write(2)l,nodes(cells(j)%c2n)%nw,cells(j)%itype !cells(j)%c2n
            else
               write(2,*)l,nodes(cells(j)%c2n)%nw,cells(j)%itype !cells(j)%c2n
            end if
         end if
      end do

      do j = 1, ndim
         write(12,20)(nodes(npart(k))%xyz(j),k=1,nn)
      end do

      do j = 1, ncells
         if(epart(j) == i) then
            l = size(cells(j)%c2n)
            c2n(:l) = cells(j)%c2n
            if(ngeom > 4) then  ! mixed mesh
               select case(ndim)
               case(2)
                  if(l == 3)  c2n(4) = c2n(3)
                  l = 4
               case(3)
                  if(l == 4) then
                     c2n(4) = c2n(3)
                     c2n(5:8) = c2n(4)
                  else if(l == 5) then
                     c2n(5:8) = c2n(5)
                  else if(l == 6) then
                     c2n(8) = c2n(6)
                     c2n(5:7) = c2n(4:6)
                     c2n(4) = c2n(3)
                  end if
                  l = 8
               end select
            end if
            write(12,10)nodes(c2n(:l))%nw
         end if
      end do
      !
      ! End of Ploting
      do j = 1, size(new_faces)
         ipl = new_faces(j)%il
         ipr = new_faces(j)%ir
         if(ipr < 0) then
            if(write_format == 0) then
               write(2)cells(ipl)%nw,ipr,&
                  size(new_faces(j)%f2n),nodes(new_faces(j)%f2n)%nw
            else
               write(2,*)cells(ipl)%nw,ipr,&
                  size(new_faces(j)%f2n),nodes(new_faces(j)%f2n)%nw
            end if
         else
            if(write_format == 0) then
               write(2)cells(ipl)%nw,cells(ipr)%nw,&
                  size(new_faces(j)%f2n),nodes(new_faces(j)%f2n)%nw
            else
               write(2,*)cells(ipl)%nw,cells(ipr)%nw,&
                  size(new_faces(j)%f2n),nodes(new_faces(j)%f2n)%nw
            end if
         end if
      end do
      !
      deallocate(npart)
      !
      call collection_period
      !
      if(write_format == 0) then
         write(2) intf   ! number of partitioning faces
      else
         write(2,*) intf
      end if
      nn = 0
      nc = 0
      do j = 1, nparts
         if(myinterface(i,j)%n + myinterface(i,j)%m /= 0) nn = nn + 1
         if(myinterface(j,i)%n + myinterface(j,i)%m /= 0) nc = nc + 1
      end do
      if(write_format == 0) then
         write(2)nn,sum(myinterface(i,:)%n) + sum(myinterface(i,:)%m)    ! send
      else
         write(2,*)nn,sum(myinterface(i,:)%n) + sum(myinterface(i,:)%m)    ! send
      end if
      do j = 1, nparts
         k1 = myinterface(i,j)%n
         k2 = myinterface(i,j)%m
         if(k1+k2 /= 0) then
            if(write_format == 0) then
               write(2)j-1,k1,k2,cells(myinterface(i,j)%cl(:k1))%nw,&
                  faces(myinterface(i,j)%cl(k1+1:))%nw
            else
               write(2,*)j-1,k1,k2,cells(myinterface(i,j)%cl(:k1))%nw,&
                  faces(myinterface(i,j)%cl(k1+1:))%nw
            end if
         end if
      end do

      if(write_format == 0) then
         write(2)nc,sum(myinterface(:,i)%n) + sum(myinterface(:,i)%m)    ! receive      
      else
         write(2,*)nc,sum(myinterface(:,i)%n) + sum(myinterface(:,i)%m)    ! receive      
      end if

      do j = 1, nparts
         k1 = myinterface(j,i)%n
         k2 = myinterface(j,i)%m
         if(k1 + k2 /= 0) then
            if(write_format == 0) then
               write(2)j-1,k1+k2
            else
               write(2,*)j-1,k1+k2
            end if
         end if
      end do
      !
      do j = 1, nparts
         do l = 1,myinterface(j,i)%n
            k1 = myinterface(j,i)%cl(l)
            k3 = 0
            do k2 = 1, size(cells(k1)%c2n)
               if(nodes(cells(k1)%c2n(k2))%nw == 0) cycle
               k3 = k3 + 1
               c2n(k3) = nodes(cells(k1)%c2n(k2))%nw
            end do
            if(k3 == 0) print *,'something wrong ',nodes(cells(k1)%c2n)%nw
            if(write_format == 0) then
               write(2)k3,c2n(:k3)
            else
               write(2,*)k3,c2n(:k3)
            end if
            !
         end do

         do k = 1,myinterface(j,i)%m
            l = k + myinterface(j,i)%n
            k1 = myinterface(j,i)%cl(l)
            k3 = 0
            do k2 = 1, size(faces(k1)%f2n)
               if(nodes(faces(k1)%f2n(k2))%nw == 0) cycle
               k3 = k3 + 1
               c2n(k3) = nodes(faces(k1)%f2n(k2))%nw
            end do
            if(k3 == 0) print *,'something wrong ',nodes(faces(k1)%f2n)%nw
            if(write_format == 0) then
               write(2)k3,c2n(:k3)
            else
               write(2,*)k3,c2n(:k3)
            end if
            !
         end do
      end do
      close(2)
      !
   end do   ! partitioning
   10 format(10I8)
   20  format(5e20.8)

contains

   subroutine check_parts
      nvc = 1
      vc(1) = cells(1)%itype
      do i = 2, ncells
         if(any((vc(:nvc) - cells(i)%itype) == 0)) cycle
         nvc = nvc + 1
         vc(nvc) = cells(i)%itype
      end do
      call sort(vc(:nvc))
   end subroutine check_parts 

   subroutine collection_period
      if(nperiod < 0) return
      if(write_format == 0) then
         write(2)fb + intf - nperiod_p,lab_period
      else
         write(2,*)fb + intf - nperiod_p,lab_period
      end if

      nc = 0
      nn = 0
      do j = 1, nparts
         if(mypinterface(i,j)%m /= 0) nn = nn + 1
         if(mypinterface(j,i)%m /= 0) nc = nc + 1
      end do

      if(write_format == 0) then
         write(2)nn,sum(mypinterface(i,:)%m)    ! send
      else
         write(2,*)nn,sum(mypinterface(i,:)%m)    ! send
      end if
      do j = 1, nparts
         k2 = mypinterface(i,j)%m
         if(write_format == 0) then
            if(k2 /= 0) write(2)j-1,k2,cells(mypinterface(i,j)%cl(:k2))%nw
         else
            if(k2 /= 0) write(2,*)j-1,k2,cells(mypinterface(i,j)%cl(:k2))%nw
         end if
      end do
      !
      if(write_format == 0) then
         write(2)nc,sum(mypinterface(:,i)%m)    ! receive      
      else
         write(2,*)nc,sum(mypinterface(:,i)%m)    ! receive      
      end if
      do j = 1, nparts
         k2 = mypinterface(j,i)%m
         if(write_format == 0) then
            if(k2 /= 0) write(2)j-1,k2,faces(mypinterface(j,i)%cl(k2+1:))%nw
         else
            if(k2 /= 0) write(2,*)j-1,k2,faces(mypinterface(j,i)%cl(k2+1:))%nw
         end if
      end do
      !
   end subroutine collection_period

   subroutine setup_periodic_interface
      integer::i,j,k
      integer,pointer::n,m
      !
      allocate(mypinterface(nparts,nparts))
      !
      !
      mypinterface%m = 0
      do k = nperiod+1, fb    !size(faces)
         i = epart(faces(k)%il)
         j = epart(abs(faces(k)%ir))
         m => mypinterface(j,i)%m
         m = m + 1
      end do
      !
      do i = 1, nparts
         do j = 1, nparts
            m => mypinterface(i,j)%m
            if(m == 0) cycle
            allocate(mypinterface(i,j)%cl(m*2))
         end do
      end do
      mypinterface%n = mypinterface%m
      mypinterface%m = 0
      do k = nperiod+1, fb    !size(faces)
         i = epart(faces(k)%il)
         j = epart(abs(faces(k)%ir))
         m => mypinterface(j,i)%m
         n => mypinterface(j,i)%n
         m = m + 1
         mypinterface(j,i)%cl(m) = abs(faces(k)%ir)   ! cell
         mypinterface(j,i)%cl(n+m) = k    !face	
      end do

   end subroutine setup_periodic_interface

   ! gather interface info to myinterface
   ! this includes interface cells as interface boundary faces
   subroutine setup_interface
      integer::i,j,k,ii
      integer,pointer::n,m
      !
      allocate(myinterface(nparts,nparts))
      !
      myinterface%n = 0   ! cell
      myinterface%m = 0   ! face
      do i = 1, nparts
         nodes%nw = 0
         do j = 1, ncells
            if(epart(j) /= i) cycle
            nodes(cells(j)%c2n)%nw = 1
         end do

         ! do ii = 1, size(nodes)
         !    print *, ii, nodes(ii) % nw
         ! end do
         !
         cells%nw = 0   
         faces%nw = 0
         do j = 1, ncells
            if(epart(j) == i) cycle 
            !     if(cells(j)%nw /= 0) cycle
            if(any(nodes(cells(j)%c2n)%nw /= 0)) then
               !      cells(j)%nw = 1
               myinterface(epart(j),i)%n = myinterface(epart(j),i)%n + 1
            end if
         end do
         do k = 1, fb   !size(faces)  !fb
            !     if(faces(k)%ir > 0) cycle
            !     if(faces(k)%nw /= 0) cycle
            j = faces(k)%il
            if(epart(j) == i) cycle 
            if(any(nodes(faces(k)%f2n)%nw /= 0)) then
               !      faces(j)%nw = 1
               myinterface(epart(j),i)%m = myinterface(epart(j),i)%m + 1
            end if
         end do

      end do
      !
      ! 
      do i = 1, nparts
         do j = 1, nparts
            n => myinterface(i,j)%n
            m => myinterface(i,j)%m
            if(n+m == 0) cycle
            allocate(myinterface(i,j)%cl(n+m))
         end do
      end do
      myinterface%n = 0
      myinterface%m = 0
      do i = 1, nparts
         nodes%nw = 0
         do j = 1, ncells
            if(epart(j) /= i) cycle
            nodes(cells(j)%c2n)%nw = 1
         end do
         !
         cells%nw = 0   
         do j = 1, ncells
            if(epart(j) == i) cycle 
            !     if(cells(j)%nw /= 0) cycle
            if(any(nodes(cells(j)%c2n)%nw /= 0)) then
               !      cells(j)%nw = 1
               n => myinterface(epart(j),i)%n
               n = n + 1
               myinterface(epart(j),i)%cl(n) = j
            end if
         end do
         !
         faces%nw = 0
         do k = 1, fb    !size(faces)
            !      if(faces(k)%ir > 0) cycle
            !      if(faces(k)%nw /= 0) cycle
            j = faces(k)%il
            if(epart(j) == i) cycle 
            if(any(nodes(faces(k)%f2n)%nw /= 0)) then
               !      faces(j)%nw = 1
               m => myinterface(epart(j),i)%m
               n => myinterface(epart(j),i)%n
               m = m + 1
               myinterface(epart(j),i)%cl(m+n) = k
            end if
         end do

      end do

   end subroutine setup_interface


   subroutine Geom_Partxyz(p,nps,ix)
      integer,intent(out)::p(:)
      integer::nps,ix,i
      real(rfp)::xmax,xmin,dx
      xmin = nodes(1)%xyz(ix)
      xmax = xmin
      do i = 1, size(nodes)
         xmax = max(xmax,nodes(i)%xyz(ix))
         xmin = min(xmin,nodes(i)%xyz(ix))
      end do
      dx = (xmax - xmin) / real(nps,rfp)
      do i = 1, size(nodes)
         p(i) = (nodes(i)%xyz(ix) - xmin) / dx + 1
         p(i) = min(nps,p(i))
         p(i) = max(1,p(i))
      end do
   end subroutine geom_partxyz

   !Xuxiao
   !partition along xyz not by distance but by number
   subroutine Geom_Partijk(p,nps,ix)
      integer,intent(out)::p(:)
      integer::nps,ix,i,ctxx
      real(rfp),allocatable::xcoord(:),xcd(:)

      allocate(xcoord(10000))
      call init_array(xcoord)

      call append_array(xcoord,nodes(1)%xyz(ix))
      ctxx = 1
      do i = 2, size(nodes)
         if(ishere(xcoord,ctxx,nodes(i)%xyz(ix)) == 0) then
            call append_array(xcoord,nodes(i)%xyz(ix))
            ctxx = ctxx + 1
         end if     
      end do

      allocate(xcd(ctxx))
      do i = 1, ctxx
         xcd(i) = xcoord(i)
      end do

      call sort_array(xcd)
      do i = 1, size(nodes)
         p(i) = findindex(xcd,nodes(i)%xyz(ix),nps)
         if(p(i)>nps .or. p(i)<1) then
            print *, 'wrong epart(i)', p(i), i
         end if
         p(i) = min(nps,p(i))
         p(i) = max(1,p(i))
      end do 
      deallocate(xcd, xcoord)
   end subroutine Geom_Partijk

   !some subroutines to operate an array written by Xuxiao
   subroutine init_array(a)
      real(rfp),allocatable::a(:)
      integer i
      do i = 1, size(a)
         a(i) = 1e20
      end do
   end subroutine init_array

   integer function issame(a,b)
      real(rfp)::a,b,e
      e = 1e-12
      if(abs(a-b) < e) then
         issame = 1
      else
         issame = 0
      end if
      return  
   end function issame

   subroutine append_array(a,v)
      real(rfp),allocatable::a(:)
      real(rfp) v
      integer i
      do i = 1, size(a)
         if(a(i) == 1e20) then
            a(i) = v
            exit
         end if
      end do  
   end subroutine append_array

   integer function ishere(a,n,v)
      real(rfp),allocatable::a(:)
      real(rfp)::v
      integer::n,i
      ishere = 0
      do i = 1, n
         if(issame(a(i),v) == 1) then
            ishere = 1
            exit
         end if
      end do
   end function ishere

   subroutine sort_array(ARR)
      real(rfp),allocatable::ARR(:)
      real(rfp) a
      do j=2, size(ARR)
         a=ARR(j)
         do i=j-1,1,-1
            if (ARR(i)<=a) goto 10
            ARR(i+1)=ARR(i)
         end do
         i=0
         10    ARR(i+1)=a
      end do
   end subroutine sort_array

   integer function findindex(a,v,nps)
      real(rfp),allocatable::a(:)
      integer,allocatable::iinp(:),sinp(:)
      real(rfp)::v
      integer i,inum,it,nps,j

      allocate(iinp(nps))
      allocate(sinp(nps+1))
      do i = 1, nps
         iinp(i) = size(a)/nps
         if(i == nps) then
            iinp(i) = size(a)/nps + mod(size(a),nps)
         end if
      end do

      if(mod(size(a),nps) > 1) then
         do i = 1, mod(size(a),nps)-1
            iinp(i) = iinp(i)+1
         end do
      end if
      iinp(nps) = size(a)/nps + 1

      sinp(1) = 0
      do i = 1, nps
         sinp(i+1) = iinp(1)
         if(i>1) then
            do j = 2, i
               sinp(i+1) = sinp(i+1) + iinp(j)
            end do
         end if
      end do

      do i = 1, size(a)
         if(issame(a(i),v) == 1) then
            inum = i
            exit
         end if
      end do

      do i = 1, nps
         if(inum > sinp(i) .and. inum .le. sinp(i+1)) then
            findindex = i
         end if
      end do

   end function findindex

   subroutine countcells(ep, np)
      integer,allocatable::ep(:),ctcells(:)
      integer:: np, i, whichp

      allocate(ctcells(np))
      ctcells = 0
      do i = 1, size(ep)
         whichp = ep(i)
         if(whichp > np) then
            print *, 'in epart there is a partition index larger than nparts'
         else
            ctcells(whichp) = ctcells(whichp) + 1
         end if
      end do

      do i = 1, np
         print *, 'partition ', i, ctcells(i)
      end do

   end subroutine countcells
   !************** end of xuxiao's subroutines ***********

   subroutine lookup_key(line,ikey)
      character(*)::line
      character(len=40)::terms(4)
      integer::ikey(:),ic,i,nterms

      line = adjustl(line)
      ic = 0
      nterms = 0
      do i = 1, len_trim(line)
         if(line(i:i) /= ' ' .and. line(i:i) /= '    ') then
            if(ic == 0) then   ! new terms
               nterms = nterms + 1
               terms(nterms) = ''
               ic = 1
               if(nterms > 4) print *,'error in MFG data file'
            else
               ic = ic + 1
            end if
            terms(nterms)(ic:ic) = line(i:i)
         else
            if(ic /= 0) ic = 0
         end if
      end do
      read(terms(1),*)ic
      read(terms(4),*)ikey(ic)
   end subroutine lookup_key

   subroutine vface_2_c2c(il,ir,cells,itype)
      integer,intent(in)::il,ir,itype
      type(cell),pointer::cells(:)
      integer,pointer::n
      if( ir < 0) return
      if(cells(il)%itype /= itype) return  ! interface between two domains
      if(cells(ir)%itype /= itype) return
      n => cells(il)%nw
      n = n + 1
      cells(il)%c2c(n) = cells(ir)%ip
      !
      n => cells(ir)%nw
      n = n + 1
      cells(ir)%c2c(n) = cells(il)%ip

   end subroutine vface_2_c2c

   subroutine face_2_c2c(il,ir,cells)
      integer,intent(in)::il,ir
      type(cell),pointer::cells(:)
      integer,pointer::n
      if( ir < 0) return
      n => cells(il)%nw
      n = n + 1
      cells(il)%c2c(n) = ir
      !
      n => cells(ir)%nw
      n = n + 1
      cells(ir)%c2c(n) = il

   end subroutine face_2_c2c

   subroutine c2n_2_n2n(c2n,nodes)
      integer,intent(in)::c2n(:)
      type(node),pointer::nodes(:)
      integer::i,n,ip,it
      n = size(c2n)
      if(ndim == 2) then
         select case(n)
         case(:2)
            print *,'wrong in c2n',n
         case(3:)  !Poly
            do i = 1, n
               ip = i + 1
               if(ip > n) ip = 1
               !    call add_edge(c2n(i),c2n(ip),nodes)
            end do
         end select
      else
         select case(n)
         case(:3)
            print *,'wrong in c2n',n
         case(4)  !Tet
            it = n
            do i = 1, 3
               ip = i + 1
               if(ip > 3) ip = 1
               !    call add_edge(c2n(i),c2n(ip),nodes)
               !    call add_edge(c2n(i),c2n(it),nodes)
            end do
         case(5)  ! Pysarim
            it = n
            do i = 1, 4
               ip = i + 1
               if(ip > 4) ip = 1
               !    call add_edge(c2n(i),c2n(ip),nodes)
               !    call add_edge(c2n(i),c2n(it),nodes)
            end do
         case(6,8)  !Prisam,Hex
            it = n/2
            do i = 1, it
               ip = i + 1
               if(ip > it) ip = 1
               !    call add_edge(c2n(i),c2n(ip),nodes)
               !    call add_edge(c2n(i+it),c2n(ip+it),nodes)
               !    call add_edge(c2n(i),c2n(i+it),nodes)
            end do
         case default
            print *,'not support yet',n
         end select
      end if
   end subroutine c2n_2_n2n

   !  subroutine add_edge(i1,i2,nodes)
   !  integer,intent(in)::i1,i2
   !  type(node),pointer::nodes(:)
   !  integer,pointer::n2n(:),n
   !! 
   !  n2n => nodes(i1)%n2n
   !  n => nodes(i1)%nw
   !  if(any(n2n(:n) == i2)) goto 10
   !    n = n + 1
   !    if(n > max_n2n) print *,'please increase max_n2n',n
   !    n2n(n) = i2
   !!
   !10   n2n => nodes(i2)%n2n
   !  n => nodes(i2)%nw
   !  if(any(n2n(:n) == i1)) return
   !    n = n + 1
   !   if(n > max_n2n) print *,'please increase max_n2n',n
   !    n2n(n) = i1
   !   end subroutine add_edge

   function neal2dli(ib)result(id)
      integer::ib,id
      select case(ib)
      case(1)
         id = 2
      case(2)
         id = 2
      case(3)
         id = 1
      case(4)
         id = 3
      end select
   end function neal2dli

   subroutine boundary_info(ni,nj,nk,ibc)
      character(len=40)::boundfile
      integer,intent(in)::ni,nj,nk
      integer,intent(inout)::ibc(:,:,:)
      integer::i1,i2,j1,j2,k1,k2,n,i,it,l1,l2,l3
      ibc = 0 
      !print *,'input boundary data file'
      !read(*,'(a)')boundfile
      !boundfile='3D_mesh_boundary.dat'
      boundfile='2D_mesh_boundary.dat'
      print*,'input boundary data file is ',boundfile
      open(3,file=boundfile)
      read(3,*)n
      do i = 1, n
         read(3,*)i1,i2,j1,j2,k1,k2,it
         !
         if(i1 == 0) i1 = ni
         if(j1 == 0) j1 = nj
         if(k1 == 0) k1 = nk
         if(i2 == 0) i2 = ni
         if(j2 == 0) j2 = nj
         if(k2 == 0) k2 = nk
         !
         ! special treating for joint points
         !   ibc(i1:i2,j1:j2,k1:k2) = it - ibc(i1:i2,j1:j2,k1:k2)
         do l1 = i1,i2
            do l2 = j1,j2
               do l3 = k1,k2
                  if(ibc(l1,l2,l3) /= 0) then
                     ibc(l1,l2,l3) = -1
                  else
                     ibc(l1,l2,l3) = it
                  end if
               end do
            end do
         end do
         !
      end do
      close(3)
   end subroutine boundary_info

   function i2s(i)
      character(len=4)::i2s
      integer,intent(in)::i
      select case(i)
      case(0:9)
         write(i2s,'(i1)')i
      case(10:99)
         write(i2s,'(i2)')i
      case(100:999)
         write(i2s,'(i3)')i
      case(1000:9999)
         write(i2s,'(i4)')i
      case default
         write(i2s,'(4h****)')
      end select
   end function i2s

   subroutine su2dfd
      integer::nd
      integer,pointer::cl_grp(:)
      character(len=50)::str

      ! dimension
      read(2,*) str, nd
      if(nd /= ndim) then
         print *, 'SU mesh dimension is inconsistent with specified dimension'
         stop
      end if

      ! connectivity
      read(2,*) str, ncells
      allocate(cells(ncells))
      do i = 1, ncells
         read(2,*) l, c2n(:get_suNodeNumber(l)), k
         nc = get_suNodeNumber(l)  ! number of nodes per cell
         c2n(:nc) = c2n(:nc) + 1  ! su2 mesh number nodes from 0
         allocate(cells(i) % c2n(nc))
         cells(i) % c2n = c2n(:nc)
         nfaces = nfaces + Nface_of_cell(nc, ndim)
      end do

      ! nodal coordinates
      read(2,*) str, nnodes
      allocate(nodes(nnodes))
      do i = 1, nnodes
         allocate(nodes(i) % xyz(ndim))
         read(2,*) nodes(i) % xyz, l
      end do
   end subroutine su2dfd

   function get_suNodeNumber(id) result(n)
      implicit none
      integer :: id, n

      select case(id)
      case(3)
         n = 2  ! Line
      case(5)
         n = 3  ! Triangle
      case(9)
         n = 4  ! Quadrilateral
      case(10)
         n = 4  ! Tetrahedral
      case(12)
         n = 8  ! Hexahedral
      case(13)
         n = 6  ! Prism
      case(14)
         n = 5  ! Pyramid
      end select
   end function get_suNodeNumber

   subroutine gambit2dfd
      integer::ngrps,ndfcd,ndfvl,nelgp,iskip
      integer,pointer::cl_grp(:)
      do i = 1, 4
         read(2,*)
      end do
      !Problem size section
      read(2,*)
      read(2,*)
      read(2,*)nnodes,ncells,ngrps,nbsets,iskip,ndfvl
      read(2,*)
      !Nodal Coordinates section
      read(2,*)
      allocate(nodes(nnodes),cells(ncells),stat = istatus)
      if(istatus /= 0 ) print *, 'error in allocate cells and faces'
      do i = 1, nnodes
         allocate(nodes(i)%xyz(ndim))
         read(2,*)l,nodes(i)%xyz
      end do
      read(2,*)   ! end section
      !Element/Cell Connectivity  
      read(2,*)
      nfaces = 0
      do i = 1, ncells  
         read(2,*)l,k,nc,c2n(:nc)
         c2n(:nc) = c2n(:nc)   
         allocate(cells(i)%c2n(nc))
         cells(i)%c2n = c2n(:nc)
         ! rearrange the local connectivity to fit DFD convention
         if(nc == 8) then  ! Brick
            cells(i)%c2n(3) = c2n(4)
            cells(i)%c2n(4) = c2n(3)
            cells(i)%c2n(7) = c2n(8)
            cells(i)%c2n(8) = c2n(7)
         else if(nc == 5) then  ! Pyramid
            cells(i)%c2n(3) = c2n(4)
            cells(i)%c2n(4) = c2n(3)
         end if
         nfaces = nfaces + Nface_of_cell(nc,ndim)
      end do
      read(2,*)   ! skip head lines
      ! Groups
      do i = 1, ngrps
         read(2,*)
         read(2,'(28x,I10)')nelgp
         allocate(cl_grp(nelgp))
         read(2,*,err=10)itype
         10  read(2,*)
         read(2,*)(cl_grp(k),k = 1, nelgp) ! skip
         print *,'volume type=',itype
         cells(cl_grp)%itype = itype
         read(2,*)
         deallocate(cl_grp)
      end do

   end subroutine gambit2dfd

   subroutine fe2dfd_face

      allocate(faces(nfaces))

      faces%il = 0
      faces%ir = 0
      k = 0
      ! double-count all the faces based on cells
      do i = 1, ncells  
         call cell2face(cells(i)%c2n,k,i)
      end do
      if(k /= nfaces) print *,'error in faces list',k,nfaces

      ! do i = 1, nfaces
      !    print *, faces(i) % il
      ! end do

      call sort_face

   end subroutine fe2dfd_face

   ! extract c2n in cells to f2n in faces
   ! c2n must be aligned such that f2n makes sense
   ! if c2n is not aligned as the convention in this subroutine
   ! must rearrange c2n when reading connectivity data
   subroutine cell2face(c2n,k,ic)
      integer,intent(in)::c2n(:),ic
      integer::k,nc
      nc = size(c2n)
      select case(ndim)
      case(2)  ! 2-dimension
         if(nc < 3 ) then
            print *, ' no cell could be formed',nc
            stop
         end if
         do j = 1, nc - 1
            k = k + 1
            allocate(faces(k)%f2n(2))
            faces(k)%f2n = c2n(j:j+1)
            faces(k)%il = ic
         end do
         k = k + 1
         allocate(faces(k)%f2n(2))
         faces(k)%f2n(1) = c2n(nc)
         faces(k)%f2n(2) = c2n(1)     
         faces(k)%il = ic
      case(3)
         if(nc == 4) then  ! Tet
            allocate(faces(k+1)%f2n(3))
            faces(k+1)%f2n = c2n((/1,2,3/))
            allocate(faces(k+2)%f2n(3))
            faces(k+2)%f2n = c2n((/2,3,4/))
            allocate(faces(k+3)%f2n(3))
            faces(k+3)%f2n = c2n((/3,4,1/))
            allocate(faces(k+4)%f2n(3))
            faces(k+4)%f2n = c2n((/4,1,2/))
            faces(k+1:k+4)%il = ic
            k = k + 4
         else if(nc == 5) then  ! Pyramid
            allocate(faces(k+1)%f2n(4))
            faces(k+1)%f2n = c2n((/1,2,3,4/))
            allocate(faces(k+2)%f2n(3))
            faces(k+2)%f2n = c2n((/1,2,5/))
            allocate(faces(k+3)%f2n(3))
            faces(k+3)%f2n = c2n((/2,3,5/))
            allocate(faces(k+4)%f2n(3))
            faces(k+4)%f2n = c2n((/3,4,5/))
            allocate(faces(k+5)%f2n(3))
            faces(k+5)%f2n = c2n((/4,1,5/))
            faces(k+1:k+5)%il = ic
            k = k + 5
         else if(nc == 6) then  ! Prism
            allocate(faces(k+1)%f2n(3))
            faces(k+1)%f2n = c2n((/1,2,3/))
            allocate(faces(k+2)%f2n(3))
            faces(k+2)%f2n = c2n((/4,5,6/))
            allocate(faces(k+3)%f2n(4))
            faces(k+3)%f2n = c2n((/1,2,5,4/))
            allocate(faces(k+4)%f2n(4))
            faces(k+4)%f2n = c2n((/2,3,6,5/))
            allocate(faces(k+5)%f2n(4))
            faces(k+5)%f2n = c2n((/3,1,4,6/))
            faces(k+1:k+5)%il = ic
            k = k + 5
         else if(nc == 8) then  ! Brick
            allocate(faces(k+1)%f2n(4))
            faces(k+1)%f2n = c2n((/1,2,3,4/))
            allocate(faces(k+2)%f2n(4))
            faces(k+2)%f2n = c2n((/5,6,7,8/))
            allocate(faces(k+3)%f2n(4))
            faces(k+3)%f2n = c2n((/1,2,6,5/))
            allocate(faces(k+4)%f2n(4))
            faces(k+4)%f2n = c2n((/2,3,7,6/))
            allocate(faces(k+5)%f2n(4))
            faces(k+5)%f2n = c2n((/3,4,8,7/))
            allocate(faces(k+6)%f2n(4))
            faces(k+6)%f2n = c2n((/4,1,5,8/))
            faces(k+1:k+6)%il = ic
            k = k + 6
         else
            print *, ' no cell could be formed',nc
            stop
         end if
      end select
   end subroutine cell2face

   ! move the non-repeating to the front of faces array
   ! re-define nfaces to be the number of non-repeating faces
   subroutine sort_face
      type(face)::fswep
      integer::incr,n,i,j,k
      logical::check

      ! assign a key to each face
      ! the key is the sorted array of the nodal indices
      do i = 1, nfaces
         !    allocate(faces(i)%key(ndim))
         faces(i)%key = n2key(faces(i)%f2n)   
         ! print *, i, faces(i) % key, faces(i) % il, faces(i) % ir
      end do
      ! print *, '----------------------------'
      !
      n = nfaces
      incr = 10000
      !
      !    Loop : calculate the increment
      !
      do while(incr <= n)
         incr = 3 * incr + 1
      end do
      !
      !   Loop : Shell-Metzner sort
      !   Sort by face keys, the sorted array starts from the smallest key to largest
      !   When face key is multi-dimensional, start sorting from the first dim., then second, third
      !
      20 INCR = INCR / 3
      I = INCR + 1
      30 IF (I .GT. N) GOTO 60
      fswep = faces(I)
      J = I
      40 do k = 1, ndim
         if(faces(j-incr)%key(k) < fswep%key(k)) goto 50     
         if(faces(j-incr)%key(k) > fswep%key(k)) exit
      end do     
      !      
      !      IF (faces(J - INCR)%nw < fswep%nw) GOTO 50
      !
      faces(J) = faces(J - INCR)
      J = J - INCR
      IF (J .GT. INCR) GOTO 40
      50 faces(J) = fswep
      I = I + 1
      GOTO 30
      60 IF (INCR .GT. 1) GOTO 20

      ! do i = 1, nfaces
      !    print *, i, faces(i) % key, faces(i) % il, faces(i) % ir
      ! end do
      ! print *, '----------------------------'

      ! move non-repeating faces to front
      k = 1
      fswep = faces(1)
      do i = 2, nfaces
         check = .true.
         do j = 1, ndim
            if(faces(i)%key(j) /= fswep%key(j)) then
               check = .false.
               exit
            end if
         end do
         !   if(faces(i)%nw == fswep%nw) then
         if(check) then
            faces(k)%ir = faces(i)%il
         else
            k = k + 1
            fswep = faces(i)
            faces(k) = fswep
         end if
      end do 

      ! do i = 1, nfaces
      !    print *, i, faces(i) % key, faces(i) % il, faces(i) % ir
      ! end do
      ! print *, '----------------------------'

      nfaces = k
      !
   end subroutine sort_face
   !
   subroutine read_gambit_boundary
      integer::itype,ic,nelem,nv,it,iface,iskip,nt
      nt = 0
      do i = 1, nbsets
         read(2,*)
         read(2,*)itype,ic,nelem,iskip,iskip
         if(ic /= 1) then
            print *,'only accept cell centered boundary condition',ic
            stop
         end if
         nt = nt + nelem
         !
         print *,fb,nt
         do j = 1, nelem
            read(2,*)ic,it,iface
            call map2face_gambit(ic,it,iface,itype)
         end do
         read(2,*)
      end do
      if(fb /= nt) print *,'error, your boundary conditions are not set correctly'
      ! close(2)
   end subroutine read_gambit_boundary

   subroutine read_su_boundary
      character(len=50)::str
      integer::itype,nelem,it,iface,nt,nbounds
      integer::ff2n(10)

      read(2,*) str, nbounds
      nt = 0
      do i = 1, nbounds
         read(2,*)  ! label of the boundary
         read(2,*) str, nelem  ! number of cells at boundary
         nt = nt + nelem
         !
         print *,fb,nt
         itype = i  ! itype = the ith boundary
         do j = 1, nelem
            read(2, *) l, ff2n(:get_suNodeNumber(l))
            it = get_suNodeNumber(l) ! number of nodes at face
            ff2n(:it) = ff2n(:it) + 1  ! su2 mesh index from 0
            call map2face_su(ff2n, it, itype)
         end do
      end do
      if(fb /= nt) print *,'error, your boundary conditions are not set correctly'
      ! close(2)
   end subroutine read_su_boundary

   subroutine map2face_su(f2n,it,itype)
      integer::itype,ic,it,iface,n,nc,j
      integer::f2n(10),key(3),jb,je
      logical::check

      n = it  ! number of nodes per element at boundary
      key(:ndim) = n2key(f2n(:n))

      jb = 0
      je = fb + 1
      out_loop: do  
         j = (jb + je)  / 2
         do k = 1, ndim
            if(key(k) < faces(j)%key(k))  then
               if(je - jb == 1) then
                  print *,'could not find the boundary face',ic,key,itype
                  exit out_loop
               end if
               je = j
               cycle out_loop
            else if(key(k) > faces(j)%key(k) )  then
               if(je - jb == 1) then
                  print *,'could not find the boundary face',ic,key,itype
                  exit out_loop
               end if
               jb = j
               cycle out_loop
            end if
         end do
         exit out_loop
      end do out_loop

      faces(j)%ir = -itype
   end subroutine map2face_su

   subroutine map2face_gambit(ic,it,iface,itype)
      integer::itype,ic,it,iface,n,nc,j
      integer::f2n(10),key(3),jb,je
      logical::check
      select case(it)
      case(2:3)  ! triangle and Quad
         n = 2
         nc = size(cells(ic)%c2n)
         if(iface == nc) then
            f2n(1) = cells(ic)%c2n(nc)
            f2n(2) = cells(ic)%c2n(1)   
         else
            f2n(1:2) = cells(ic)%c2n(iface:iface+1)
         end if
      case(4)  ! Brick
         n = 4
         if(iface == 1) then
            f2n(1:4) = cells(ic)%c2n((/1,2,5,6/))
         else if(iface == 2) then
            f2n(1:4) = cells(ic)%c2n((/2,3,7,6/))
         else if(iface == 3) then
            f2n(1:4) = cells(ic)%c2n((/3,4,8,7/))
         else if(iface == 4) then
            f2n(1:4) = cells(ic)%c2n((/4,1,5,8/))
         else if(iface == 5) then
            f2n(1:4) = cells(ic)%c2n((/2,1,4,3/))
         else if(iface == 6) then
            f2n(1:4) = cells(ic)%c2n((/5,6,7,8/))
         end if
      case(5)  ! Prisam
         n = 4
         if(iface == 1) then
            f2n(1:4) = cells(ic)%c2n((/1,2,5,4/))
         else if(iface == 2) then
            f2n(1:4) = cells(ic)%c2n((/2,3,6,5/))
         else if(iface == 3) then
            f2n(1:4) = cells(ic)%c2n((/3,1,4,6/))
         else if(iface == 4) then
            n = 3
            f2n(1:3) = cells(ic)%c2n((/1,2,3/))
         else if(iface == 5) then
            n = 3
            f2n(1:3) = cells(ic)%c2n((/4,5,6/))
         end if
      case(6)  ! Tet
         n = 3
         if(iface == 1) then
            f2n(1:n) = cells(ic)%c2n((/1,2,3/))
         else if(iface == 2) then
            f2n(1:n) = cells(ic)%c2n((/1,2,4/))
         else if(iface == 3) then
            f2n(1:n) = cells(ic)%c2n((/2,3,4/))
         else if(iface == 4) then
            f2n(1:n) = cells(ic)%c2n((/4,3,1/))
         end if
      case(7)  ! Pyramid
         n = 3
         if(iface == 1) then
            n = 4
            f2n(1:n) = cells(ic)%c2n((/1,2,3,4/))
         else if(iface == 2) then
            f2n(1:n) = cells(ic)%c2n((/1,2,5/))
         else if(iface == 3) then
            f2n(1:n) = cells(ic)%c2n((/2,3,5/))
         else if(iface == 4) then
            f2n(1:n) = cells(ic)%c2n((/3,4,5/))
         else if(iface == 5) then
            f2n(1:n) = cells(ic)%c2n((/4,1,5/))
         end if     
      end select
      !
      key(:ndim) = n2key(f2n(:n))


      jb = 0
      je = fb + 1
      out_loop: do  
         j = (jb + je)  / 2
         do k = 1, ndim
            if(key(k) < faces(j)%key(k))  then
               if(je - jb == 1) then
                  print *,'could not find the boundary face',ic,key,itype
                  exit out_loop
               end if
               je = j
               cycle out_loop
            else if(key(k) > faces(j)%key(k) )  then
               if(je - jb == 1) then
                  print *,'could not find the boundary face',ic,key,itype
                  exit out_loop
               end if
               jb = j
               cycle out_loop
            end if
         end do
         exit out_loop
      end do out_loop

      faces(j)%ir = -itype

      !   do i = 1, nfaces
      !   if(faces(i)%ir > 0) exit
      !   check = .true.
      !   do j = 1, ndim
      !   if(faces(i)%key(j) /= key(j))  then
      !    check = .false.
      !    exit
      !   end if
      !   end do
      !    if(check) then
      !      faces(i)%ir = -itype
      !      exit
      !    end if
      !   end do
      !   if(i > fb) then
      !    print *,'could not find the boundary face'
      !    print *,f2n(:n),key
      !   end if
   end subroutine map2face_gambit

   function nface_of_cell(nc,ndim)result(nf)
      integer,intent(in)::ndim,nc
      integer::nf
      select case(ndim)
      case(2)  ! 2-dimension
         if(nc < 3 ) then
            print *, ' no cell could be formed',nc
            stop
         end if
         nf = nc
      case(3)
         if(nc == 4) then  ! Tet
            nf = 4
         else if(nc == 5) then  ! Pyramid
            nf = 5
         else if(nc == 6) then  ! Prism
            nf = 5
         else if(nc == 8) then  ! Brick
            nf = 6
         else
            print *, ' no cell could be formed',nc
            stop
         end if
      end select
   end function nface_of_cell  

   function n2key(f2n)result(key)
      implicit none
      integer,intent(in)::f2n(:)
      integer::ip(size(f2n)),n
      integer::key(3),shift,k,k0,i,j
      ip = f2n
      n = size(ip)
      call sort(ip)
      key = 0
      key = ip(:ndim)
      return    
      !     key = ip(1)
      !     do i = 2,n
      !      key = ieor(key,ip(i))
      !       k = ip(i) !- ip(i-1)
      !       k0 = k
      !      j = 0
      !     do while(k > 0)
      !      k = ishft(k,-1)
      !      j = j + 1
      !     end do
      !     print *,'j=',j
      !      key = ishft(key,j+1) + k0
      !      key = ishft(key,16) + k0
      !      key = key*10 + k0
      !     end do
      !
   end function n2key

   subroutine sort(x)
      integer:: x(:)
      integer:: temp
      !
      INTEGER I, J, INCR,n
      !
      n = size(x)
      if(n == 2) then
         if(x(1) > x(2) ) x(1:2) = x(2:1:-1)
         return
      end if
      INCR = 1
      !
      !        Loop : calculate the increment
      !
      do while(incr <= n)
         incr = 3 * incr + 1
      end do
      !  10 INCR = 3 * INCR + 1
      !     IF (INCR .LE. N) GOTO 10

      !
      !   Loop : Shell-Metzner sort
      !
      20 INCR = INCR / 3
      I = INCR + 1
      30 IF (I .GT. N) GOTO 60
      TEMP = X(I)
      J = I
      40 IF (X(J - INCR) .LT. TEMP) GOTO 50
      X(J) = X(J - INCR)
      J = J - INCR
      IF (J .GT. INCR) GOTO 40
      50 X(J) = TEMP
      I = I + 1
      GOTO 30
      60 IF (INCR .GT. 1) GOTO 20
      !
   end subroutine sort

   subroutine fieldview2dfd
      integer,pointer::itype(:),nf(:),f2n(:,:),bc(:)
      integer::key(3),nv,nb,nnc(4)=(/4,8,6,5/),mytype,myf2n(4)
      integer::jb,je,iv
      real(rfp)::xyz(3)
      logical::check
      read(2,'(a)')line
      iv = 0
      if(index(line,'2 4') >0) iv = 2
      if(index(line,'3 0') >0) iv = 3
      !Boundary Table
      if(iv == 0) then
         print *, 'Not recognized fieldview data format, please contact dli@purdue.edu'
         stop
      end if
      do
         read(2,'(a)')line
         if(line(:5) == 'Bound') exit
      end do
      !
      if(iv == 3) then
         read(line(15:),*)nb
      else
         read(2,*)nb
      end if
      allocate(bc(nb))
      print *,'Pgrid assigns labels for boundary patchs!!!!!!!!!!!'
      Print *,' Make sure the labels in gems.inp is matched to these labels'
      do i = 1, nb
         read(2,'(a)')line
         bc(i) = i
         print *,'Label =',bc(i),'--------->',line(:20)
         !    read(2,*)bc(i)
      end do
      !Variable Name
      !    read(2,*)
      !    read(2,*)nv
      !    do i = 1, nv
      !     read(2,*)
      !    end do
      !Nodal Coordinates section
      do 
         read(2,'(a)')line
         if(line(:5) == 'Nodes') exit
      end do
      if(iv == 3) then
         read(line(6:),*)nnodes
      else
         read(2,*)nnodes
      end if
      allocate(nodes(nnodes),stat = istatus)
      !  allocate(nodes(nnodes),cells(ncells),stat = istatus)
      if(istatus /= 0 ) print *, 'error in allocate cells and faces'
      k = 0
      do i = 1, nnodes
         allocate(nodes(i)%xyz(ndim))
         read(2,*)xyz
         nodes(i)%nw = 0
         if(ndim == 2.and. abs(xyz(3))>1.e-10_rfp) cycle
         k = k + 1
         nodes(k)%xyz = xyz(:ndim)
         nodes(i)%nw = k
      end do
      nnodes = k
      !
      ! Boundary Faces
      do i = 1, i+1
         read(2,'(a)')line
         if(line(:5) == 'Bound') exit
      end do
      if(iv == 3) then
         read(line(15:),*)nb
      else
         read(2,*)nb
      end if
      allocate(itype(nb),nf(nb),f2n(4,nb),stat = istatus)
      k = 0
      do i = 1, nb
         read(2,*)mytype,l,myf2n(:l)
         if(ndim == 2) then
            if(all(nodes(myf2n(:l))%nw ==0))cycle 
            if(all(nodes(myf2n(:l))%nw > 0))cycle
            k = k + 1
            k1 = 0
            do j = 1, l
               if(nodes(myf2n(j))%nw == 0) cycle
               k1 = k1 + 1
               myf2n(k1) = nodes(myf2n(j))%nw
            end do
            itype(k) = mytype
            nf(k) = k1
            f2n(:k1,k) = myf2n(:k1)
         else
            k = k + 1
            nf(i) = l
            itype(i) = mytype
            f2n(:l,i) = myf2n(:l)  
         end if
      end do
      nb = k
      !
      !Element/Cell Connectivity  
      read(2,*)
      !
      ! Decide the cells number
      ncells = 0
      do 
         read(2,'(a)',end=10)line
         if(line(1:3) =='Var') exit
         ncells = ncells + 1
      end do
      10 continue
      do i = 1, ncells + 1
         backspace(2)
      end do
      print *,'total number of cells =',ncells
      print *,'please waiting for face connectivity searching'
      !   
      allocate(cells(ncells),stat = istatus)
      if(istatus /= 0 ) print *, 'error in allocate cells and faces'
      nfaces = 0
      do i = 1, ncells  
         read(2,*)l,k,c2n(:nnc(l))
         nc = nnc(l)
         !
         if(ndim == 2) then
            k = 0
            do j = 1, nc
               if(nodes(c2n(j))%nw == 0) cycle
               k = k + 1
               c2n(k) = nodes(c2n(j))%nw
            end do
            nc = k
         end if

         allocate(cells(i)%c2n(nc))
         cells(i)%c2n = c2n(:nc)
         if(nc == 4) then
            cells(i)%c2n((/1,2,3,4/)) = c2n((/1,2,4,3/))
         else if(nc == 8) then  ! Brick --- stupid connectivity
            cells(i)%c2n((/1,2,3,4/)) = c2n((/1,2,4,3/))
            cells(i)%c2n((/5,6,7,8/)) = c2n((/5,6,8,7/))
         else if(nc == 6) then  ! Prism
            cells(i)%c2n((/1,2,3/)) = c2n((/1,4,6/))
            cells(i)%c2n((/4,5,6/)) = c2n((/2,3,5/))
         end if
         nfaces = nfaces + Nface_of_cell(nc,ndim)
      end do
      !  read(2,*)   ! skip head lines
      close(2)
      !
      allocate(faces(nfaces))
      faces%il = 0
      faces%ir = 0
      k = 0
      do i = 1, ncells  
         call cell2face(cells(i)%c2n,k,i)
      end do
      if(k /= nfaces) print *,'error in faces list',k,nfaces

      call sort_face

      ! do swap
      fb = 0
      do i = 1, nfaces
         if(faces(i)%ir <= 0) then
            fb = fb + 1
            if (i == fb) cycle
            fswap = faces(fb)
            faces(fb) = faces(i)
            faces(i) = fswap
         end if
      end do
      !
      do i = 1, nb
         !
         key = n2key(f2n(:nf(i),i))
         jb = 0
         je = fb + 1
         out_loop: do  
            j = (jb + je)  / 2
            do k = 1, ndim
               if(key(k) < faces(j)%key(k))  then
                  if(je - jb == 1) then
                     print *,'could not find the boundary face',i,key,itype(i)
                     exit out_loop
                  end if
                  je = j
                  cycle out_loop
               else if(key(k) > faces(j)%key(k) )  then
                  if(je - jb == 1) then
                     print *,'could not find the boundary face',i,key,itype(i)
                     exit out_loop
                  end if
                  jb = j
                  cycle out_loop
               end if
            end do
            exit out_loop
         end do out_loop

         !   if(check) then
         faces(j)%ir = -bc(itype(i))
         !     exit
         !    end if
         !    end do
         !   if(j > fb) print *,'could not find the boundary face',i,key,itype(i)
         ! if(mod(i,10) == 0) print *,i
      end do
      !
      ! 
      deallocate(itype,nf,f2n)

   end subroutine fieldview2dfd

   subroutine treatment_periodic_bc
      integer::iy,nb,iax,id,npb
      real(rfp)::x0(2),rci,rcj,zci,zcj,eps = 1.e-5_rfp
      write(*,*)'Do you have periodic boundary condition 0=no,1=cylindrical(3D),2=planar(3D),3=2D case'
      !read(*,*)id
      id = 3
      print *, 'periodic BC = ', id

      if(id ==0) return
      write(*,*)'Please input the labels for Periodic BC'
      !read(*,*)lab_period
      if(ipbc == 1) then
         lab_period = 1
      else if(ipbc == 2) then
         lab_period = 2
      end if
      print *, 'periodic label = ', lab_period

      if(id == 1) then
         write(*,*)'please tell me which axis is cylindrical axis,1=x,2=y,3=z?'
      else if(id == 2) then
         write(*,*)'which coordinates are used to match,1=y,z 2=x,z 3=x,y?'
      else
         write(*,*)'which coordinate is used to match,1=x,2=y?'
      end if
      !read(*,*)iax
      if(ipbc == 1) then
         iax = 2
      else if(ipbc == 2) then
         iax = 1
      end if
      print *, 'axis ', iax, 'is chosen'

      if(id == 1) then
         write(*,*)'please input the center of periodic surface?(0.0,0.0)'
         read(*,*)x0
      end if

      nb = 0
      do i = 1,nfaces 
         if(faces(i)%ir > 0) exit
         nb = nb + 1
      end do

      do i = 1, nb
         if(faces(i)%ir > 0 .or. abs(faces(i)%ir) /= lab_period) cycle
         iy = 0
         ! find matching coordinates save to rci and zci
         call rc_zc(i,x0,rci,zci,iax,id)
         do j = i+1,nb
            if(faces(j)%ir > 0 .or. abs(faces(j)%ir) /= lab_period) cycle
            call rc_zc(j,x0,rcj,zcj,iax,id)
            if(ndim == 3.and.abs(zci-zcj) > eps)cycle
            if(abs(rci - rcj) > eps)cycle
            ! they are matched
            faces(i)%ir = faces(j)%il
            faces(j)%ir = faces(i)%il
            iy = 1
            exit
         end do
         if(iy == 0) then
            print *,'Couldnot find matched face',i,faces(i)%f2n
         end if
      end do

      ! do swap
      iy = 0
      if(ipbc == 1) then  ! first periodic BC
         npb = nb
      else  ! build on the first periodic BC
         npb = nperiod
      end if

      do i = 1, npb
         if(faces(i)%ir <= 0) then
            iy = iy + 1
            if (i == iy) cycle
            fswap = faces(iy)
            faces(iy) = faces(i)
            faces(i) = fswap
         end if
      end do

      ! So periodic boundary faces located between iy+1:nb
      if(ipbc == 1) then
         faces(iy+1:nb)%ir = - faces(iy+1:nb)%ir
         npbcs(ipbc) = nb - iy
      else
         faces(iy+1:nperiod)%ir = - faces(iy+1:nperiod)%ir
         npbcs(ipbc) = nperiod - iy
      end if
      nperiod = iy
      lpbcs(ipbc) = lab_period
   end subroutine treatment_periodic_bc

   subroutine correction_periodic_bc
      integer::il,ir
      real(rfp)::x0(2),rci,rcj,zci,zcj,eps = 1.e-10_rfp

      if(nperiod < 0) return

      do i = nperiod + 1, fb
         il = faces(i)%il
         ir = faces(i)%ir

         faces(i)%ir = - ir  !lab_period
         if(epart(il) == epart(ir)) cycle
         epart(ir) = epart(il)
      end do
      !    
   end subroutine correction_periodic_bc


   subroutine rc_zc(i,x0,rc,zc,iax,id)
      integer,intent(in)::i,iax
      integer::n,k,k1,id
      real(rfp),intent(in)::x0(2)
      real(rfp),intent(out)::rc,zc
      real(rfp)::pc(ndim),pc1(3)

      pc = 0.0_rfp
      n = size(faces(i)%f2n)
      !
      do k = 1, n
         pc = pc + nodes(faces(i)%f2n(k))%xyz
      end do
      pc = pc / real(n,rfp)
      pc1 = 0.0
      pc1(:ndim) = pc
      !
      select case(id)
      case(1)
         zc = pc(iax)
         k1 = 0
         rc = 0.0_rfp
         do k = 1, ndim
            if(k == iax) cycle
            k1 = k1 + 1
            rc = rc + (pc(k) - x0(k1))**2
         end do
         rc = sqrt(rc)
      case(2)
         select case(iax)
         case(1)
            rc = pc1(2);zc=pc1(3)
         case(2)
            rc = pc1(1);zc=pc1(3)
         case(3)
            rc = pc1(1);zc=pc1(2)
         case default
            print *,'You input wrong coordinate',iax
         end select
      case(3)
         rc = pc(iax)
      case default
         print *,'Sorry, I can not handle this option:',id
      end select
      !
   end subroutine rc_zc

end program pgrid
