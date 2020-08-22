module gems_geom
   use gems_data
   use gems_library
   use mpi
   implicit none
contains
   !
   subroutine read_grid(gridfile,cells,faces,nodes,interf,pinterf)
      character(len=*)::gridfile
      type(cell),pointer::cells(:),cc
      type(face),pointer::faces(:)
      type(node),pointer::nodes(:)
      type(itf)::interf,pinterf
      logical::yes
      !
      integer::mb(16)
      integer::i,j,n,nn,nc,nf,nparts,nt,ntf,nitf,itype,lf,rf,id,ierr,st(mpi_status_size)
      integer,pointer::c2n(:,:),f2n(:,:),nlist(:),clist(:),pnlist(:),pclist(:)
      integer::tnn,tnns,tnc,tncs,np1,np,ip,m,nperiod,lab_p
      integer::pnn,pnns,pnc,pncs
      real(rfp),pointer::xyz(:,:)
      equivalence (mb(1),nn),(mb(2),nc),(mb(3),nf),(mb(4),nt),(mb(5),ntf),(mb(6),nitf)
      equivalence (mb(7),tnn),(mb(8),tnns),(mb(9),tnc),(mb(10),tncs)
      equivalence (mb(11),nperiod),(mb(12),lab_p)
      equivalence (mb(13),pnn),(mb(14),pnns),(mb(15),pnc),(mb(16),pncs)
      !
      call mpi_comm_size(mpi_comm_world,n,ierr)
      call mpi_comm_rank(mpi_comm_world,id,ierr)
      nperiod = -1; tnns = 0; tncs = 0;

      !only use one PE to read all grid data
      if(id == 0) then
         nt = 2**ndim + 2
         ntf = 2**ndim + 1
         if(ndim==1) ntf=ntf+1
         do i = n - 1, 0, -1
            inquire(file=gridfile(:len_trim(gridfile))//'.'//i2s(i+1),exist=yes)
            if(.not.yes) then
               print *,'Sorry, could not find the grid file ---->>>',trim(gridfile)//'.'//i2s(i+1)
               call mpi_abort(mpi_comm_world,ierr,ierr)
               stop
            end if

            open(2,file=gridfile(:len_trim(gridfile))//'.'//i2s(i+1),form='binary',status='old')

            !# of nodes, cells, faces and parts
            read(2)nn,nc,nf,nparts
            if( nparts /= n) then
               print *,'grid partitions:',nparts,' /= number of process ',n
            end if
            !
            allocate(xyz(nn,ndim),c2n(nc,nt),f2n(nf,ntf))
            !read coordinates of nodes 
            do j = 1, nn
               read(2)xyz(j,:)
            end do
            !       xyz(:,1) = xyz(:,1)*two  
            xyz = xyz * b_lenref

            !read: # of nodes per cell
            !read: indices of neighbor nodes per cell
            !read: cell itype
            do j = 1, nc
               read(2)c2n(j,1),c2n(j,2:c2n(j,1)+2)
            end do

            !read face left&right cell
            !read # of nodes per face
            !read indices of neighbor nodes per face
            do j = 1, nf
               read(2)f2n(j,:2),f2n(j,3),f2n(j,4:3+f2n(j,3))
            end do
            !
            if(s_iperiod > 0) then
               read(2)nperiod,lab_p
               read(2)pnc,pncs      ! cell information
               allocate(pclist(2*pnc + pncs))
               np1=2*pnc + 1
               do j = 1, pnc
                  read(2)ip,m,pclist(np1:np1+m-1)
                  pclist(j) = ip
                  pclist(j+pnc) = m
                  np1 = np1+m
               end do
               !
               read(2)pnn,pnns      ! face information
               allocate(pnlist(2 * pnn + pnns))
               np1 = 2 * pnn + 1 
               do j = 1, pnn
                  read(2)ip,m,pnlist(np1:np1+m-1)
                  pnlist(j) = ip
                  pnlist(j+pnn) = m
                  np1 = np1 + m
               end do

            end if

            !
            !      read(2)ni
            !      allocate(exf(2*ni))
            !      read(2)exf
            !
            if(n == 1) then
               nitf = 0
               close(2)
               cycle
            end if 

            !read # of cells per interface
            read(2)nitf

            !send
            !read # of interface related to this part
            !read # of cell+face related to this part's interface
            read(2)tnc,tncs
            allocate(clist(3 * tnc + tncs))
            np1=3*tnc + 1
            do j = 1, tnc

               !read index of the part
               !read # of cells
               !read # of faces
               !read indices of cells and faces 
               !save all of this information to clist for sending (that is clist)
               !the first few indices of clist are destination id, # of cells to send for each destination and # of face to send for each destination
               read(2)ip,np,m,clist(np1:np1+np+m-1)
               clist(j) = ip
               clist(j+tnc) = np
               clist(j+2*tnc) = m
               np1 = np1+np+m
            end do

            !receives
            !save receiving information to nlist
            !the first few indices of nlist are origin id, # of cells and faces to receive for each origin
            read(2)tnn,tnns
            allocate(nlist(2 * tnn + tnns * ntf))
            nlist = 0
            !
            do j = 1, tnn
               read(2)nlist(j),nlist(tnn + j)
            end do
            !
            np1 = 2 * tnn + 1
            !
            do j =1,tnns
               read(2)np,nlist(np1+1:np1+np)
               nlist(np1) = np
               np1 = np1 + ntf
            end do
            !
            close(2)

            !send the information (mb,xyz,c2n,f2n,nlist,clist) from id=0 to all PE's
            if(i == 0) cycle   ! do not need to send

            call mpi_send(mb,size(mb),MPI_INTEGER,i,0,mpi_comm_world,ierr)
            call mpi_send(xyz,nn*ndim,MPI_DOUBLE_PRECISION,i,1,mpi_comm_world,ierr)
            call mpi_send(c2n,nc*nt,MPI_INTEGER,i,2,mpi_comm_world,ierr)
            call mpi_send(f2n,nf*ntf,MPI_INTEGER,i,3,mpi_comm_world,ierr)
            !      call mpi_send(exf,2*ni,MPI_INTEGER,i,4,mpi_comm_world,ierr)
            call mpi_send(nlist,size(nlist),MPI_INTEGER,i,5,mpi_comm_world,ierr)
            call mpi_send(clist,size(clist),MPI_INTEGER,i,6,mpi_comm_world,ierr)
            !
            if(s_iperiod > 0) then
               call mpi_send(pnlist,size(pnlist),MPI_INTEGER,i,7,mpi_comm_world,ierr)
               call mpi_send(pclist,size(pclist),MPI_INTEGER,i,8,mpi_comm_world,ierr)
               deallocate(pnlist,pclist)
            end  if
            !
            deallocate(xyz,c2n,f2n,nlist,clist)

         end do

      else
         !receive the information sent from id=0  
         call mpi_recv(mb,size(mb),MPI_INTEGER,0,0,mpi_comm_world,st,ierr)
         !
         allocate(xyz(nn,ndim),c2n(nc,nt),f2n(nf,ntf),nlist(2*tnn+tnns*ntf),clist(3*tnc+tncs))

         !
         call mpi_recv(xyz,nn*ndim,MPI_DOUBLE_PRECISION,0,1,mpi_comm_world,st,ierr)
         call mpi_recv(c2n,nc*nt,MPI_INTEGER,0,2,mpi_comm_world,st,ierr)
         call mpi_recv(f2n,nf*ntf,MPI_INTEGER,0,3,mpi_comm_world,st,ierr)
         !     call mpi_recv(exf,2*ni,MPI_INTEGER,0,4,mpi_comm_world,st,ierr)
         call mpi_recv(nlist,size(nlist),MPI_INTEGER,0,5,mpi_comm_world,st,ierr)
         call mpi_recv(clist,size(clist),MPI_INTEGER,0,6,mpi_comm_world,st,ierr)
         !
         if(s_iperiod > 0) then
            allocate(pclist(2 * pnc + pncs),pnlist(2 * pnn + pnns))
            call mpi_recv(pnlist,size(pnlist),MPI_INTEGER,0,7,mpi_comm_world,st,ierr)
            call mpi_recv(pclist,size(pclist),MPI_INTEGER,0,8,mpi_comm_world,st,ierr)
         end if
         !
      end if
      !  
      allocate(nodes(nn),cells(nc),faces(nf),interf%pcell(tnns))

      do i = 1, nn
         nodes(i)%xyz = xyz(i,:)
      end do
      deallocate(xyz)
      !
      do i = 1, nc
         nt = c2n(i,1)
         allocate(cells(i)%c2n(nt),cells(i)%weight(nt),cells(i)%gradient(neq,s_nrec),stat=ierr)
         cells(i)%gradient = zero   ! initial
         cells(i)%weight   = zero
         if(neqm > 0) allocate(cells(i)%jv(3))
         if(ierr /= 0 ) print *, 'error in allocate gradient'
         if(s_idt >= 1) allocate(cells(i)%qvn(s_idt))
         cells(i)%c2n = c2n(i,2:nt+1)
         cells(i)%itype = c2n(i,nt+2)
         if(cells(i)%itype /= 0) print *, cells(i)%itype
      end do

      !
      do i = 1, nf
         lf = f2n(i,1)
         rf = f2n(i,2)
         ntf = f2n(i,3)
         allocate(faces(i)%f2n(ntf))
         if(i <= nitf) then   ! interface face
            itype = partition_face_no
            faces(i)%left_cell => cells(lf)
            faces(i)%right_cell => interf%pcell(rf)
         else
            if(rf < 0) then ! boundary face
               faces(i)%left_cell => cells(lf)
               if(nperiod >= 0.and.i > nperiod) then
                  itype = lab_p
                  !            faces(i)%right_cell => pinterf%pcell(abs(rf))
               else
                  allocate(faces(i)%right_cell,stat=ierr)
                  if(ierr /= 0 ) print *, 'error in allocate boundary faces'
                  itype = abs(rf)
                  !           if(s_ivis == 0.and.itype == 3) itype = 2
               end if
            else ! interior face
               itype = 0
               if(s_ifmt == 1) then
                  faces(i)%left_cell => cells(lf)
                  faces(i)%right_cell => cells(rf)
               else
                  ! unstructured grid, make sure lf < rf
                  faces(i)%left_cell => cells(min(lf,rf))
                  faces(i)%right_cell => cells(max(lf,rf))
               end if
            end if
         end if
         !	
         faces(i)%itype = itype
         faces(i)%f2n = f2n(i,4:ntf+3)
      end do
      deallocate(c2n,f2n)
      !
      !  Periodic Boundary Condition
      !
      if(s_iperiod > 0) then
         ! Receiving.............
         allocate(pinterf%pcell(pnns),pinterf%np(pnn),pinterf%nnp(pnn))

         do i = 1, pnns
            nt = pnlist(2*pnn+i)
            cc => pinterf%pcell(i)
            faces(nt)%right_cell => cc  ! no need to allocate face right cell
         end do

         pinterf%nitf = nperiod
         pinterf%nn = pnn
         pinterf%np = pnlist(:pnn)
         pinterf%nnp = pnlist(pnn+1:2*pnn)
         if(sum(pinterf%nnp) /= pnns) print *,'error.............'
         !
         ! Sending...........
         !
         allocate(pinterf%cp(pnc),pinterf%ncp(pnc))
         allocate(pinterf%scell(pncs))
         !
         ip = 0
         do i = 1, pnc
            m = pclist(pnc + i)
            do j = 1, m
               ip = ip + 1
               pinterf%scell(ip)%to_cell => cells(pclist(2*pnc+ip))
            end do
         end do
         !
         if(ip /= pncs) print *,'err......',np1,tncs
         pinterf%nc = pnc
         pinterf%cp = pclist(:pnc)
         pinterf%ncp = pclist(pnc+1:2*pnc)
         if(sum(pinterf%ncp) /= pncs) print *,'error.............'
         deallocate(pnlist,pclist)
         !
      end if

      if(n == 1) then
         interf%nitf = nitf  ! if single core, nitf = 0
         return
      end if
      !
      allocate(interf%np(tnn),interf%nnp(tnn))
      np1 = 2 * tnn + 1
      ntf = 2**ndim + 1
      do i = 1, tnns
         nt = nlist(np1)
         cc => interf%pcell(i)
         allocate(cc%c2n(nt),cc%weight(nt),cc%gradient(neq,s_nrec),stat=ierr)
         if(neqm > 0) allocate(cc%jv(3))
         if(ierr /= 0 ) print *, 'error in allocate gradient'
         cc%c2n = nlist(np1+1:np1+nt)
         np1 = np1 + ntf
      end do
      !save in interf%nn: number of batches to receive      
      !save in interf%np: origin id for each batch to receive
      !save in interf%nnp: how many cells&faces in each batch
      interf%nitf = nitf
      interf%nn = tnn
      interf%np = nlist(:tnn)
      interf%nnp = nlist(tnn+1:2*tnn)
      if(sum(interf%nnp) /= tnns) print *,'error.............'
      !
      allocate(interf%cp(tnc),interf%ncp(tnc))
      allocate(interf%scell(tncs))
      np1 = 0
      ip = 0
      do i = 1, tnc
         np = clist(tnc + i)
         m  = clist(2*tnc + i)
         do j = 1, np
            np1 = np1 + 1
            ip  = ip + 1
            interf%scell(np1)%to_cell => cells(clist(3*tnc+ip))
         end do

         do j = 1, m
            np1 = np1 + 1
            ip  = ip + 1
            interf%scell(np1)%to_cell => faces(clist(3*tnc+ip))%right_cell
         end do
      end do
      if(np1 /= tncs) print *,'err......',np1,tncs

      !save in interf%nc: number of batches to send
      !save in interf%cp: destination id for each batch to send
      !save in interf%ncp: how many cells&faces in each batch
      interf%nc = tnc
      interf%cp = clist(:tnc)
      interf%ncp = clist(tnc+1:2*tnc) + clist(2*tnc+1:3*tnc)
      if(sum(interf%ncp) /= tncs) print *,'error.............'
      deallocate(nlist,clist)
      !
   end subroutine read_grid

   subroutine read_grid_interp(gridfile,cells,faces,nodes,interf,pinterf,ipart)  
      character(len=*)::gridfile
      type(cell),pointer::cells(:),cc
      type(face),pointer::faces(:)
      type(node),pointer::nodes(:)
      type(itf)::interf,pinterf
      logical::yes
      !
      integer::mb(16),ipart
      integer::i,j,n,nn,nc,nf,nparts,nt,ntf,nitf,itype,lf,rf,id,ierr,st(mpi_status_size)
      integer,pointer::c2n(:,:),f2n(:,:),nlist(:),clist(:),pnlist(:),pclist(:)
      integer::tnn,tnns,tnc,tncs,np1,np,ip,m,nperiod,lab_p
      integer::pnn,pnns,pnc,pncs
      real(rfp),pointer::xyz(:,:)
      equivalence (mb(1),nn),(mb(2),nc),(mb(3),nf),(mb(4),nt),(mb(5),ntf),(mb(6),nitf)
      equivalence (mb(7),tnn),(mb(8),tnns),(mb(9),tnc),(mb(10),tncs)
      equivalence (mb(11),nperiod),(mb(12),lab_p)
      equivalence (mb(13),pnn),(mb(14),pnns),(mb(15),pnc),(mb(16),pncs)
      !
      nperiod = 0; tnns = 0; tncs = 0;
      !
      nt = 2**ndim + 2
      ntf = 2**ndim + 1
      inquire(file=trim(gridfile)//'.'//i2s(ipart),exist=yes)
      if(.not.yes) then
         print *,'Sorry, could not find the grid file ---->>>',gridfile//'.'//i2s(ipart)
         !     call mpi_abort(mpi_comm_world,ierr,ierr)
         return
      end if
      !
      open(2,file=trim(gridfile)//'.'//i2s(ipart),form='binary',status='old')

      read(2)nn,nc,nf,nparts
      !
      ! see if reach to last partition?
      !
      !    if(ipart == nparts) last=.true.
      !
      allocate(xyz(nn,ndim),c2n(nc,nt),f2n(nf,ntf))
      do j = 1, nn
         read(2)xyz(j,:)
      end do
      xyz = xyz * b_lenref
      do j = 1, nc
         read(2)c2n(j,1),c2n(j,2:c2n(j,1)+2)
      end do
      !
      do j = 1, nf
         read(2)f2n(j,:2),f2n(j,3),f2n(j,4:3+f2n(j,3))
      end do
      !
      !
      if(nparts == 1) then
         nitf = 0
      else
         read(2)nitf
      end if
      close(2)
      !
      allocate(nodes(nn),cells(nc),faces(nf),interf%pcell(tnns))

      do i = 1, nn
         nodes(i)%xyz = xyz(i,:)
      end do
      deallocate(xyz)
      !
      do i = 1, nc
         nt = c2n(i,1)
         allocate(cells(i)%c2n(nt),cells(i)%weight(nt),cells(i)%gradient(neq,s_nrec),stat=ierr)
         if(neqm > 0) allocate(cells(i)%jv(3))
         if(ierr /= 0 ) print *, 'error in allocate gradient'
         if(s_idt >= 1) allocate(cells(i)%qvn(s_idt))
         cells(i)%c2n = c2n(i,2:nt+1)
         cells(i)%itype = c2n(i,nt+2)
      end do
      !
      do i = 1, nf
         lf = f2n(i,1)
         rf = f2n(i,2)
         ntf = f2n(i,3)
         allocate(faces(i)%f2n(ntf))
         if(i <= nitf) then   ! interface
            itype = partition_face_no
            faces(i)%left_cell => cells(lf)
            allocate(faces(i)%right_cell,stat=ierr)
         else
            if(rf < 0) then
               faces(i)%left_cell => cells(lf)
               allocate(faces(i)%right_cell,stat=ierr)
               if(ierr /= 0 ) print *, 'error in allocate boundary faces'
               itype = abs(rf)
            else
               itype = 0
               faces(i)%left_cell => cells(min(lf,rf))
               faces(i)%right_cell => cells(max(lf,rf))
            end if
         end if
         !	
         faces(i)%itype = itype
         faces(i)%f2n = f2n(i,4:ntf+3)
      end do
      deallocate(c2n,f2n)
      !
      interf%nitf = 0
      !
   end subroutine read_grid_interp

   !
   subroutine geometry(nodes,cells,faces)
      implicit none
      type(node),pointer::nodes(:)
      type(cell),pointer::cells(:),cc
      type(face),pointer::faces(:),cf,cfmax
      integer,pointer::f2n(:)
      real(rfp)::vol,vecn(ndim),tvol
      integer::i,j,istatus,n,id
      !
      call mpi_comm_rank(mpi_comm_world,id,n)

      !
      !  if(s_nstart > 0) then
      ! read the grid information
      !   open(2,file='dfd.grd.'//i2s(id),form='unformatted')
      !   read(2)cells%vol,(cells%centp(i),i=1,ndim)
      !   read(2)faces%area,(faces%vecn(i),faces%centp(i),i=1,ndim)
      !   close(2)
      !   goto 10
      !  end if

      ! center of cells
      !
      tvol = zero
      do i = 1, size(cells)
         call element_center(nodes(cells(i)%c2n),cells(i)%vol,cells(i)%centp)
         tvol = tvol + cells(i)%vol
      end do
      !  tvol = zero
      !  do i = 1, size(cells)
      !   n = size(cells(i)%c2n)
      !   call element_center(nodes(cells(i)%c2n),cells(i)%vol,cells(i)%centp)
      !   tvol = tvol + cells(i)%vol
      !   do j = 1, ndim
      !   cells(i)%centp(j) = sum(nodes(cells(i)%c2n)%xyz(j)) / real(n,rfp)
      !   end do
      !  end do
      ! center of faces
      !
      do i = 1, size(faces)
         call face_center(nodes(faces(i)%f2n),faces(i)%area,faces(i)%centp)
      end do

      !  do i = 1, size(faces)
      !     n = size(faces(i)%f2n)
      !    call face_center(nodes(faces(i)%f2n),faces(i)%area,faces(i)%centp)
      !   do j = 1, ndim
      !    print *,faces(i)%centp(j),sum(nodes(faces(i)%f2n)%xyz(j)) / real(n,rfp)
      !   end do
      !  end do
      !
      ! calculate volume and normal vector
      !
      !  cells%vol = zero
      !
      do i = 1, size(faces)
         cf => faces(i)
         f2n => cf%f2n
         if(ndim == 2) then
          !  vecn(:ndim) = area2d(nodes(f2n(1))%xyz,nodes(f2n(2))%xyz)
          !  cc => cf%left_cell
          !  vol = volume2d(cc%centp,nodes(f2n(1))%xyz,nodes(f2n(2))%xyz)
          !  if( vol < zero ) vecn = -vecn
          !  cf%area = sqrt(sum(vecn*vecn))
          !  vecn = vecn / cf%area
          !  !
          !  cf%vecn = vecn
          !  !    cc%vol = cc%vol + abs(vol)
          !  if(cf%itype /= 0) cycle
          !  cc => cf%right_cell
          !  vol = volume2d(cc%centp,nodes(f2n(1))%xyz,nodes(f2n(2))%xyz)
          !  !    cc%vol = cc%vol + abs(vol)
         else if(ndim==3) then
          !  !
          !  cc => cf%left_cell
          !  vecn = zero
          !  vol = zero
          !  n = size(f2n)
          !  !   do j = 1, n - 2
          !  !    vecn=vecn + area3d(nodes(f2n(1))%xyz,  &
          !  !         nodes(f2n(j+1))%xyz,nodes(f2n(j+2))%xyz)
          !  !    vol= vol+volume3d(cc%centp,nodes(f2n(1))%xyz, &
          !  !                nodes(f2n(j+1))%xyz,nodes(f2n(j+2))%xyz)
          !  !   end do
          !  do j = 1, n - 1
          !     vecn=vecn + area3d(cf%centp,  &
          !        nodes(f2n(j))%xyz,nodes(f2n(j+1))%xyz)
          !     vol= vol+volume3d(cc%centp,cf%centp, &
          !        nodes(f2n(j))%xyz,nodes(f2n(j+1))%xyz)
          !  end do
          !  vecn=vecn + area3d(cf%centp,  &
          !     nodes(f2n(n))%xyz,nodes(f2n(1))%xyz)
          !  vol= vol+volume3d(cc%centp,cf%centp, &
          !     nodes(f2n(n))%xyz,nodes(f2n(1))%xyz)
          !  !
          !  ! decide the normal vector's direction
          !  !
          !  !   cc%vol = cc%vol + abs(vol)
          !  if(cf%itype ==0) then
          !     vol = dot_product(vecn,cf%right_cell%centp - cc%centp)
          !  else
          !     vol = dot_product(vecn,cf%centp - cc%centp)
          !  end if
          !  if( vol < zero) then
          !     if(s_ifmt /= 1) then
          !        vecn = -vecn
          !     end if
          !  end if
          !  !
          !  cf%area = sqrt(sum(vecn*vecn))
          !  if(cf%area < 1.e-15_rfp) then
          !     print *,'too small face',cf%f2n
          !     print *,nodes(cf%f2n(1))%xyz
          !     print *,nodes(cf%f2n(2))%xyz
          !  end if
          !  vecn = vecn / cf%area
          !  !
          !  cf%vecn = vecn
          !  !
          !  if(cf%itype /= 0) cycle   ! boundary cells
          !  cc => cf%right_cell
          !  vol=zero
          !  do j = 1, n - 1
          !     vol= vol+volume3d(cc%centp,cf%centp, &
          !        nodes(f2n(j))%xyz,nodes(f2n(j+1))%xyz)
          !  end do
          !  vol= vol+volume3d(cc%centp,cf%centp, &
          !     nodes(f2n(n))%xyz,nodes(f2n(1))%xyz)
          !  !    cc%vol = cc%vol + abs(vol)
         else
            cc => cf%left_cell
            vecn=-(cc%centp(1)-cf%centp(1))/abs(cc%centp(1)-cf%centp(1))
            cf%vecn=vecn
         end if
      end do
      !
      ! save the grid information
      !  open(2,file='dfd.grd.'//i2s(id),form='unformatted')
      !  write(2)cells%vol,(cells%centp(i),i=1,ndim)
      !  write(2)faces%area,(faces%vecn(i),faces%centp(i),i=1,ndim)
      !  close(2)

      10  cells%srf = zero
      ! borrow srf to store the number of faces surounding the cell  
      do i = 1, size(faces)
         cf => faces(i)
         cc => cf%left_cell
         cc%srf = cc%srf + one
         if(cf%itype /= 0) cycle
         cc => cf%right_cell
         cc%srf = cc%srf + one
      end do
      !
      vol=zero
      do i = 1, size(cells)
         cc => cells(i)
         vol=vol+cc%vol
         if(cc%vol < 1.0e-20_rfp) print *,'No.',i,cc%vol,'is too small'
         n = int(cc%srf + 0.1)   ! worry about the error( + 0.5)
         allocate(cells(i)%sface(n),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate sface and svd'
      end do
      !
      if(s_nstart == 0.and. id == 0) then
         print *, 'checking the volume'
         print *, 'total volume calculated by two methods =',tvol,vol,tvol-vol
      end if
      !
      cells%srf = zero
      !
      vecn = zero
      !find sface of cells
      do i = 1, size(faces)
         cf => faces(i)

         !calculate the boundary area
         if(cf%itype /= 0) &
            vecn = vecn + cf%vecn * cf%area
         !
         cc => cf%left_cell
         cc%srf = cc%srf + one
         n = int(cc%srf + 0.1)
         !    cc%sface(n)%to_face => cf
         if(n <= 1) then
            cc%sface(n)%to_face => cf
         else
            cfmax => cc%sface(1)%to_face
            if(cf%area > cfmax%area) then
               cc%sface(1)%to_face => cf
               cc%sface(n)%to_face => cfmax
            else
               cc%sface(n)%to_face => cf
            end if
         end if
         ! right cell
         if(cf%itype /= 0) cycle
         cc => cf%right_cell
         cc%srf = cc%srf + one
         n = int(cc%srf + 0.1)
         !   cc%sface(n)%to_face => cf
         if(n <= 1) then
            cc%sface(n)%to_face => cf
         else
            cfmax => cc%sface(1)%to_face
            if(cf%area > cfmax%area) then
               cc%sface(1)%to_face => cf
               cc%sface(n)%to_face => cfmax
            else
               cc%sface(n)%to_face => cf
            end if
         end if
      end do  

      if(sum(vecn*vecn) > 1.0e-16_rfp) then
         print *,'computational domain is not closed'
         print *, 'total surface area='
         write(*,'(3e16.6)')vecn
      end if

      !
   end subroutine geometry

   subroutine setup_chain(cells,faces,raw)
      implicit none
      !
      type(raw_chain)::raw(ndim)
      type(cell),pointer::cells(:),cc
      type(face),pointer::faces(:)
      type(chain),pointer::cchain,next_chain
      integer::k,i,j,nf,index(1),is,im,ih,ic
      integer::nc(size(cells)),ics(size(cells)),iend(size(cells))
      !
      do k = 1, ndim
         !
         im = ndim - k + 1
         if(im == k) im = 3
         !
         do i=1,size(cells)
            cells(i)%srf = real(i,rfp)   ! assigen srf = index of cells
         end do
         nf = 0
         is = 1
         do while(any(cells%srf > zero))
            nf = nf + 1
            ! select min coordination 
            index = minloc(cells%centp(im),cells%srf > zero)
            !
            call search_chain(index(1),cells,k,nc(nf),ics(is:),iend(nf))
            is = is + nc(nf)
         end do

         !
         if(nf > 0) allocate(raw(k)%bchain(nf),raw(k)%nchain(nf))
         is = 0
         do i = 1, nf
            nullify(raw(k)%bchain(i)%up_chain)
            cchain => raw(k)%bchain(i)
            raw(k)%nchain(i) = nc(i)
            do j = 1, nc(i)
               is = is + 1
               ic = ics(is)
               cchain%cc => cells(ic)
               ih = abs(cells(ic)%srf) + 0.1_rfp
               cchain%cf => cells(ic)%sface(ih)%to_face
               allocate(next_chain)
               cchain%next_chain => next_chain
               next_chain%up_chain => cchain
               cchain => next_chain
            end do
            cchain%cf => cells(ic)%sface(iend(i))%to_face
         end do
         !
         do i = 1, size(faces)
            if(faces(i)%itype < 0) faces(i)%itype = 0
         end do
         !
      end do  ! end dimension
      !
   end subroutine setup_chain

   subroutine search_chain(istart,cells,idir,nc,ics,iend)
      implicit none
      type(cell),pointer::cells(:),cc,cn
      type(face),pointer::cf,cfr
      type(neighbour),pointer::sf(:)
      integer,intent(in)::idir,istart
      integer::nf,i,j,im(1),in(1),nc,k1,k2,iend
      integer::ics(:),icf(size(ics)),icb(size(ics))
      real(kind=rfp)::vn(20)
      cc => cells(istart)
      sf => cc%sface
      nf = size(sf)
      do i = 1, nf
         cfr => sf(i)%to_face
         vn(i) = cfr%vecn(idir)
         if(associated(cfr%right_cell,cc)) vn(i) = - vn(i)
      end do
      in = maxloc(vn(:nf))
      im = minloc(vn(:nf))
      !  forward search
      cf => sf(in(1))%to_face
      k1 = 1
      icf(k1) = cc%srf + 0.1_rfp
      cc%srf = -im(1)    ! record left face no.
      do while(cf%itype == 0)
         cf%itype = -1
         k1 = k1 + 1
         cn => cf%left_cell
         if(associated(cn,cc)) cn => cf%right_cell
         icf(k1) = cn%srf + 0.1_rfp
         sf => cn%sface
         nf = size(sf)
         do i = 1, nf
            cfr => sf(i)%to_face
            if(associated(cf,cfr)) cn%srf = -i
            vn(i) = cfr%vecn(idir)
            if(associated(cfr%right_cell,cn)) vn(i) = - vn(i)
         end do
         in = maxloc(vn(:nf))
         cf => sf(in(1))%to_face
         cc => cn
      end do
      iend = in(1)
      !
      !  backward search   
      !
      cc => cells(istart)
      sf => cc%sface
      cfr => sf(im(1))%to_face
      k2 = 0
      do while(cfr%itype == 0)
         cfr%itype = -1
         k2 = k2 + 1
         cn => cfr%left_cell
         if(associated(cn,cc)) cn => cfr%right_cell
         icb(k2) = cn%srf + 0.1_rfp
         sf => cn%sface
         nf = size(sf)
         do i = 1, nf
            cfr => sf(i)%to_face
            vn(i) = cfr%vecn(idir)
            if(associated(cfr%right_cell,cn)) vn(i) = - vn(i)
         end do
         im = minloc(vn(:nf))
         cn%srf = - im(1)
         cfr => sf(im(1))%to_face
         cc => cn
      end do
      !
      ! join forward and backward chains
      !
      nc = 0   
      do i =k2, 1, -1
         nc = nc + 1
         ics(nc) = icb(i)
         sf => cells(icb(i))%sface
         do j = 1, size(sf)
            cfr => sf(j)%to_face
            if(cfr%itype == 0) cfr%itype = -1
         end do
      end do
      do i =1, k1
         nc = nc + 1
         ics(nc) = icf(i)
         sf => cells(icf(i))%sface
         do j = 1, size(sf)
            cfr => sf(j)%to_face
            if(cfr%itype == 0) cfr%itype = -1
         end do
      end do
   end subroutine search_chain

   function average_center(mynodes)result(c)
      real(rfp)::c(ndim)
      type(node)::mynodes(:)
      integer::i
      real(rfp)::inv_n
      inv_n = one / real(size(mynodes),rfp)
      do i = 1, ndim
         c(i) = sum(mynodes%xyz(i)) * inv_n
      end do
   end function average_center

   recursive subroutine face_center(mynodes,s,c)
      real(rfp)::c(ndim),c1(ndim),c2(ndim),s,s1,s2,sv(ndim)
      type(node)::mynodes(:)
      integer::n
      n = size(mynodes)
      if(ndim == 2) then
         c = average_center(mynodes)
      else
         if(n <= 3) then
            if(n > 1) then
               sv= area3d(mynodes(1)%xyz,    &   
                  mynodes(2)%xyz,    &
                  mynodes(3)%xyz)
               s = sqrt(sum(sv*sv))
               c = average_center(mynodes)
            else
               s = 1.
               c = mynodes(1)%xyz
            end if
         else
            call face_center(mynodes(:n-1),s1,c1)
            call face_center(mynodes((/1,n-1,n/)),s2,c2)
            s = s1 + s2
            if(s > mytiny) then
               c = (c1 * s1 + c2 * s2)/s
            else
               c = c1
            end if
         end if
      end if
   end subroutine face_center

   recursive subroutine element_center(mynodes,vol,c)
      real(rfp)::c(ndim),c1(ndim),c2(ndim),vol,vol1,vol2
      type(node)::mynodes(:)
      integer::n
      n = size(mynodes)
      if(ndim == 2) then
       !  select case(n)
       !  case(3)  ! triangle
       !     c = average_center(mynodes)
       !     vol = volume2d(mynodes(1)%xyz,  &
       !        mynodes(2)%xyz,  &
       !        mynodes(3)%xyz)
       !     vol = abs(vol)
       !  case (4:)
       !     call element_center(mynodes(:n-1),vol1,c1)
       !     call element_center(mynodes((/1,n-1,n/)),vol2,c2)
       !     vol = vol1 + vol2
       !     c = (vol1 * c1 + vol2 * c2 ) / vol
       !  end select
      else if(ndim == 3) then
       !  select case(n)
       !  case(4)   ! tet
       !     c = average_center(mynodes)
       !     vol = volume3d(mynodes(1)%xyz,  &
       !        mynodes(2)%xyz,  &
       !        mynodes(3)%xyz,  &
       !        mynodes(4)%xyz)
       !     vol = abs(vol)
       !  case(5)   ! pysarim
       !     call element_center(mynodes((/1,2,3,5/)),vol1,c1)
       !     call element_center(mynodes((/1,3,4,5/)),vol2,c2)
       !     vol = vol1 + vol2
       !     c = (vol1 * c1 + vol2 * c2 ) / vol
       !  case(6)   ! prism
       !     call element_center(mynodes((/1,3,6,4,5/)),vol1,c1)
       !     call element_center(mynodes((/1,2,3,5/)),vol2,c2)
       !     vol = vol1 + vol2
       !     c = (vol1 * c1 + vol2 * c2 ) / vol
       !  case(8)   ! hex
       !     call element_center(mynodes((/1,2,3,5,6,7/)),vol1,c1)
       !     call element_center(mynodes((/1,3,4,5,7,8/)),vol2,c2)
       !     vol = vol1 + vol2
       !     c = (vol1 * c1 + vol2 * c2 ) / vol
       !  case default
       !     print *,'something wrong'
       !  end select
       !  !   vol = vol1 + vol2
       !  !  if(vol > mytiny) then
       !  !   c = (vol1 * c1 + vol2 * c2 ) / vol
       !  !  else
       !  !   c = c1
       !  !  end if

       !  !ndim=1, 1D case
      else
         vol=abs(mynodes(2)%xyz(1)-mynodes(1)%xyz(1))
         c=(mynodes(2)%xyz(1)+mynodes(1)%xyz(1))/2

      end if   
   end subroutine element_center

   subroutine search_neighbor(nodes,cells,faces)
      type(node),pointer::nodes(:)
      type(cell),pointer::cells(:),cc,nc
      type(face),pointer::faces(:),cf
      !  type(neighbour_cell)::neighbor(size(nodes),50),mnc(50)
      integer,pointer::c2n(:)
      integer::nn(size(nodes))
      integer::i,istatus,nl,nrec,j,k,l,id,ierr
      real(rfp)::vecn(ndim)
      ! initial   
      ! borrow srf to store the number of faces surounding the cell  
      !
      call mpi_comm_rank(mpi_comm_world,id,ierr)
      nrec = size(cells(1)%gradient,dim=2)
      do i = 1, size(cells)
         cc => cells(i)
         c2n => cc%c2n
         vecn = zero
         nl = size(cc%sface)
         allocate(cc%scell(nl),cc%svd(nl,nrec),stat=istatus)
         if(istatus /= 0 ) print *, 'error in allocate sface and svd'
         do j = 1,nl
            cf => cc%sface(j)%to_face
            if(associated(cc,cf%left_cell)) then
               cc%scell(j)%to_cell => cf%right_cell
               vecn = vecn + cf%vecn * cf%area
            else
               cc%scell(j)%to_cell => cf%left_cell
               vecn = vecn - cf%vecn * cf%area
            end if
         end do
         if(sum(vecn*vecn) > 1.0e-16)then
            print *,'this cell is not closed', i,sum(vecn*vecn),c2n,' at',id
         end if
      end do
      !
      return

      !  nn = 0
      !  do i = 1, size(cells)
      !   cc => cells(i)
      !   c2n => cc%c2n
      !   do j = 1, size(c2n)
      !     nn(c2n(j)) = nn(c2n(j)) + 1
      !     neighbor(c2n(j),nn(c2n(j)))%to_cell => cc
      !   end do
      !  end do
      !  
      !  do i = 1, size(faces)       ! for ghost cell
      !   cf => faces(i)
      !   if(cf%itype == 0) exit
      !    c2n => cf%f2n
      !    do j = 1, size(c2n)
      !     nn(c2n(j)) = nn(c2n(j)) + 1
      !     neighbor(c2n(j),nn(c2n(j)))%to_cell => cf%right_cell   
      !    end do
      !  end do
      !   
      !  do i = 1, size(cells)
      !   cc => cells(i)   
      !   c2n => cc%c2n
      !   
      !   nl = 0
      !   do j = 1, size(c2n)
      !kloop:  do k = 1, nn(c2n(j))
      !     nc => neighbor(c2n(j),k)%to_cell
      !     if(associated(cc,nc)) cycle kloop
      !     do l = 1, nl
      !      if(associated(nc,mnc(l)%to_cell)) cycle kloop
      !     end do
      !     nl = nl + 1
      !     mnc(nl)%to_cell => nc
      !    end do kloop
      !   end do
      !   if(nl > 50) print *,'please increase array in search sub'
      !   allocate(cc%scell(nl),cc%svd(nl,nrec),stat=istatus)
      !   if(istatus /= 0 ) print *, 'error in allocate sface and svd'
      !    cc%scell(:) = mnc(:nl)
      ! end do
      !

      !   nn = 0
      !   nl = 0
      !   cells%srf = -1.0_rfp
      !   cc%srf = 1.0_rfp
      !   call cell_neighbor(cc,nn,neighbor,1)
      !   allocate(cc%scell(nn),cc%svd(nn,nrec),stat=istatus)
      !   if(istatus /= 0 ) print *, 'error in allocate sface and svd'
      !   cc%scell(:) = neighbor(1:nn)
      !  end do

   end subroutine search_neighbor

   !  recursive subroutine cell_neighbor(cc,nn,neighbor,nl)
   !   type(cell),pointer::cc,nc
   !   type(face),pointer::cf
   !   type(neighbour_cell)::neighbor(:)
   !   integer::nn,nl
   !   integer::j,is
   !   is = abs(s_ialg)
   !   if(is == 0) is = 1
   !  if(nl > is) return
   !   do j = 1, size(cc%sface)
   !    cf => cc%sface(j)%to_face
   !	if(cf%itype /= 0) then
   !	 nn = nn + 1
   !	 neighbor(nn)%to_cell => cf%right_cell
   !	else
   !     nc => cf%left_cell
   !	 if(nc%srf < zero) then
   !	  nn = nn + 1
   !	  neighbor(nn)%to_cell => nc
   !	  nc%srf = 1.0_rfp
   !      call cell_neighbor(nc,nn,neighbor,nl+1)
   !	 end if
   !     nc => cf%right_cell
   !	 if(nc%srf < zero) then
   !	  nn = nn + 1
   !	  neighbor(nn)%to_cell => nc
   !	  nc%srf = 1.0_rfp
   ! 	  call cell_neighbor(nc,nn,neighbor,nl+1)
   ! 	 end if
   !	end if
   !  end do  
   ! end subroutine cell_neighbor

   !subroutine face_geom(nodes,cells,faces)
   !   implicit none
   !   type(node),pointer::nodes(:)
   !   type(cell),pointer::cells(:),cl,cr
   !   type(face),pointer::faces(:),cf
   !   integer,pointer::f2n(:)
   !   real(rfp)::vol,vecn(ndim)
   !   integer::i,j,j1,n,ia
   !   ia = s_iaxis
   !   !
   !   do i = 1, size(faces)
   !      cf => faces(i)
   !      f2n => cf%f2n
   !      n = size(f2n)
   !      cl => cf%left_cell; cr=> cf%right_cell
   !      allocate(cf%avec(ndim,n+2))
   !      cf%avec = zero
   !      !
   !      if(ndim == 2) then
   !         vol = volume2d(cl%centp,nodes(f2n(1))%xyz,nodes(f2n(2))%xyz)
   !         cf%vl = abs(vol)
   !         !   if(s_iaxis > 0) cf%vl = cf%vl * one_third *abs(sum(nodes(f2n(:))%xyz(ia)) + cl%centp(ia))
   !         vecn(:ndim) = area2d(nodes(f2n(2))%xyz,cl%centp)
   !         vecn = sign(one,vol) * vecn
   !         !   if(s_iaxis > 0) vecn = vecn * half *abs(nodes(f2n(2))%xyz(ia) + cl%centp(ia))
   !         cf%avec(:,1) = cf%avec(:,1) + vecn(:ndim)
   !         cf%avec(:,4) = cf%avec(:,4) + vecn(:ndim)
   !         !
   !         vecn(:ndim) = area2d(cl%centp,nodes(f2n(1))%xyz)
   !         vecn = sign(one,vol) * vecn
   !         !   if(s_iaxis > 0) vecn = vecn * half *abs(nodes(f2n(1))%xyz(ia) + cl%centp(ia))
   !         cf%avec(:,1) = cf%avec(:,1) + vecn(:ndim)
   !         cf%avec(:,3) = cf%avec(:,3) + vecn(:ndim)
   !         vol = volume2d(cr%centp,nodes(f2n(1))%xyz,nodes(f2n(2))%xyz)
   !         cf%vr = abs(vol)
   !         !   if(s_iaxis > 0) cf%vr = cf%vr * one_third *abs(sum(nodes(f2n(:))%xyz(ia)) + cr%centp(ia))
   !         vecn(:ndim) = area2d(nodes(f2n(2))%xyz,cr%centp)
   !         vecn = sign(one,vol) * vecn
   !         !   if(s_iaxis > 0) vecn = vecn * half *abs(nodes(f2n(2))%xyz(ia) + cr%centp(ia))
   !         cf%avec(:,2) = cf%avec(:,2) + vecn(:ndim)
   !         cf%avec(:,4) = cf%avec(:,4) + vecn(:ndim)
   !         !
   !         vecn(:ndim) = area2d(cr%centp,nodes(f2n(1))%xyz)
   !         vecn = sign(one,vol) * vecn
   !         !   if(s_iaxis > 0) vecn = vecn * half *abs(nodes(f2n(1))%xyz(ia) + cr%centp(ia))
   !         cf%avec(:,2) = cf%avec(:,2) + vecn(:ndim)
   !         cf%avec(:,3) = cf%avec(:,3) + vecn(:ndim)
   !      else
   !         vecn = zero
   !         cf%vl = zero;  cf%vr = zero;
   !         do j = 1, n
   !            j1 = j + 1
   !            if(j1 > n) j1 = 1
   !            vecn = area3d(cl%centp,  &
   !               nodes(f2n(j1))%xyz,nodes(f2n(j))%xyz)
   !            vol= volume3d(cl%centp, &
   !               nodes(f2n(j))%xyz,nodes(f2n(j1))%xyz,cf%centp)
   !            if(sum(vecn * ( cl%centp - cf%centp)) < zero) vecn = -vecn
   !            cf%vl = cf%vl + abs(vol)
   !            cf%avec(:,1) = cf%avec(:,1) +  vecn
   !            cf%avec(:,j+2) = cf%avec(:,j+2) +  vecn
   !            cf%avec(:,j1+2) = cf%avec(:,j1+2) +  vecn

   !            vecn = area3d(cr%centp,  &
   !               nodes(f2n(j))%xyz,nodes(f2n(j1))%xyz)
   !            vol= volume3d(cr%centp, &
   !               nodes(f2n(j1))%xyz,nodes(f2n(j))%xyz,cf%centp)
   !            if(sum(vecn * ( cr%centp - cf%centp)) < zero) vecn = -vecn
   !            cf%vr = cf%vr + abs(vol)
   !            cf%avec(:,2) = cf%avec(:,2) +  vecn
   !            cf%avec(:,j+2) = cf%avec(:,j+2) +  vecn
   !            cf%avec(:,j1+2) = cf%avec(:,j1+2) +  vecn
   !         end do
   !      end if
   !      cf%avec = cf%avec / (real(ndim,rfp) * (cf%vl + cf%vr))
   !   end do
   !end subroutine face_geom


end module gems_geom
