MODULE dfd_turb
  USE dfd_library
  USE dfd_state
  USE dfd_data
  USE mpi
implicit none
 contains

!********************************************************
 subroutine turb_model_2006(cells)
!********************************************************
implicit none
!
  type(cell),pointer::cells(:),cc
  type(face),pointer::cf
! local variables
  type(vector)::qv,gqv(ndim),rhox
  real(rfp)::gvel(ndim,ndim),gk(ndim),go(ndim)
  real(rfp)::st,so,zmut,vk,vo,ta,xwi,xw,xk,w3inv,vecn(ndim),o1,o2
  real(rfp)::beta,betas,su,fb,fbs,resk,reso,rho,vol,rt,alpha
  real(rfp)::ksi=1.5_rfp,m_t0=0.0625_rfp,m_t,a   ! compressible correction
  real(rfp)::lt     ! length scale for DES
  real(rfp)::sigmad,dkdw

  integer::i,j,it,k,k1,k2
!
  
  if(s_ivis /= 2.or.neqf == 0) return  ! no k-omega turbulent model
!
  do i = 1, size(cells)
!
  cc => cells(i)
  vol = cc%vol
!  do k = 1, ndim
!   gqv(k)%v = zero
!  end do
!  do k = 1, size(cc%sface)
!   cf => cc%sface(k)%to_face
!   qv = sum_qv_node(nodes(cf%f2n)) / real(size(cf%f2n),rfp)
!   vecn = cf%vecn * cf%area
!   do k2 = 1,neq
!   gqv%v(k2) = gqv%v(k2) + qv%v(k2) * vecn
!   end do 
!  end do
!  do k = 1, ndim
!   gqv(k)%v = gqv(k) / cc%vol
!  end do
!
  gvel = cc%gradient(imb:ime,1:ndim)  !du/dx
  gk   = cc%gradient(ike,1:ndim)  !dk/dxi
  go   = cc%gradient(iom,1:ndim)  !do/dxi
  qv   = cc%qv
  vk = qv%v(ike)
  vo = qv%v(iom)
!
  dkdw = sum(gk * go)   ! Wilcox(2006)
  ta = zero
  do k1 = 1, ndim
   do k2 = 1, ndim
    ta = ta + (gvel(k1,k2) + gvel(k2,k1))**2 * half 
   end do
  end do
  ta = sqrt(ta / t_betas0) * 7._rfp/ 8._rfp
  qv%v(iom) = max(vo,ta)
  cc%qv%v(iom) = qv%v(iom)
!
  rho = rhof(qv)
  zmut = mu_turb(qv)
!
  su = zero
  do k = 1, ndim
  su = su + gvel(k,k)
  end do
  st = zero
  do k1 = 1, ndim
   do k2 = 1, ndim
    ta = zmut * (gvel(k1,k2) + gvel(k2,k1))
    if(k1 == k2) ta = ta - two_third * (zmut * su + rho * vk)
    st = st + ta * (gvel(k1,k2) + gvel(k2,k1)) * half
   end do
   end do
  xw = zero  
  xk = zero
  do k = 1, ndim
   do k1 = 1, ndim
    if(k == k1) cycle
    o1 = (gvel(k,k1) - gvel(k1,k))
! rotating frame
!    if(abs(s_omega) > mytiny) o1 = o1 + permutation(two,k,s_iaxis,k1) * s_omega
    do k2 = 1, ndim
    if(k1 == k2) cycle
    o2 = (gvel(k1,k2) - gvel(k2,k1))
! rotating frame
!    if(abs(s_omega) > mytiny) o2 = o2 + permutation(two,k1,s_iaxis,k2) * s_omega
    ta = o1 * o2 * (gvel(k2,k ) + gvel(k ,k2))
    xw = xw + ta
   end do
  end do
    xk = xk + gk(k) * go(k)
 end do
  w3inv = one / vo**3
  xw = 0.125_rfp * abs(xw * w3inv) / t_betas3
!  xk = xk * w3inv
! Wilcox(2006)
!  fb = (one + 70.0_rfp * xw) / (one + 80.0_rfp * xw)
  fb = (one + 85.0_rfp * xw) / (one + 100.0_rfp * xw)
!  xk = max(zero,xk)
!  if(xk > 1.e10_rfp) then
!   fbs = 1.7_rfp
!  else
!  xk = xk * xk
!  fbs = (one + 680.0_rfp * xk) / (one + 400.0_rfp * xk)
!  end if

! try Low-Reynolds-number version
!
!  rt = zmut / mumix(qv)
  
!  alpha = (0.025_rfp + rt / 6.0_rfp)/(one + rt / 6.0_rfp)
!  alpha = (0.1_rfp + rt/2.7_rfp)/(one + rt/2.7_rfp)/alpha
!  rt = (rt/8.0_rfp)**4
!  fbs = fbs *  (5.0_rfp/18.0_rfp + rt)/(one + rt)
!
!
!
  beta  = t_beta0  * fb
!  betas = t_betas0 * fbs
  betas = t_betas0    ! Wilcox(2006)
!
! Compressible correction
!
  m_t = two * vk / sound_speed(qv)
  if(m_t > m_t0) then
   m_t = ksi * (m_t - m_t0)   ! ksi*(M_t^2-M_t0^2)
   beta = beta - betas * m_t
   betas = betas *( one + m_t)
  end if
!
  so = t_alpha * vo / vk * st
!  so = t_alpha * vo / vk * st * alpha  ! low-RE
! DES correction
 if(t_cdes <= zero) then   ! no change  
  resk = st - betas * rho * vk * vo
 else
  if(ndim == 2.and.s_iaxis > 0) then
  lt = (vol/cc%centp(s_iaxis))**(one/real(ndim,rfp)) * t_cdes
  else
  lt = (vol)**(one/real(ndim,rfp)) * t_cdes
  end if
  lt = min(lt, sqrt(vk) / (betas * vo))
  resk = st - rho * vk**1.5_rfp / lt
 end if

!  reso = so - beta * rho * vo * vo
! Wilcox 2006
  sigmad = max(zero,dkdw * t_psido / vo )
  reso = so - beta * rho * vo * vo + rho * sigmad

  cc%res%v(ike) = cc%res%v(ike) + resk * vol
  cc%res%v(iom) = cc%res%v(iom) + reso * vol
!  LHS
!
  rhox = rhoxf(qv)
!
!  st = st * vol
!  so = so * vol
!  rho = rho * vol
!
!  su =  two_third * rho * su
!  cc%dm%e(ike,ike) = cc%dm%e(ike,ike) + betas * rho * vo   + st / vk
!!       
!  cc%dm%e(ike,iom) = cc%dm%e(ike,iom) + betas * rho * vk   - st / vo
!!
!  cc%dm%e(iom,iom) = cc%dm%e(iom,iom) + two * beta * vo * rho + su * t_alpha 

!cycle
  cc%dm%e(iom,:neqf) = cc%dm%e(iom,:neqf) + ((two_third * t_alpha * ( su + vo) * su + &
                        beta * vo * vo ) * vol) * rhox%v(:neqf)
  cc%dm%e(iom,iom) = cc%dm%e(iom,iom) + rho * vol * ( two_third * t_alpha * su + two * beta * vo)
!                                         
!   dmut/dQp
!
!
  rhox%v(:neqf) = (zmut / rho) * rhox%v(:neqf)
  rhox%v(ike)   = rhox%v(ike) + rho / vo
  rhox%v(iom)   = rhox%v(iom) - zmut / vo 
!
  cc%dm%e(ike,:neqf) = cc%dm%e(ike,:neqf) + (vol * (two_third * ( su + vo) * su + &
                                            betas * vo * vo ) ) * rhox%v(:neqf)
  cc%dm%e(ike,iom) = cc%dm%e(ike,iom) + two_third * zmut * ( su + three * betas * vo) * vol
!  cc%dm%e(ike,ike) = cc%dm%e(ike,ike) +  rho * betas * vo * vol
!  cc%dm%e(ike,iom) = cc%dm%e(ike,iom) +  rho * betas * vk * vol
!

 end do
!
 end subroutine turb_model_2006
!
!********************************************************
 subroutine turb_model_1998(cells)
!********************************************************
implicit none
!
  type(cell),pointer::cells(:),cc
  type(face),pointer::cf
! local variables
  type(vector)::qv,gqv(ndim)
  real(rfp)::gvel(ndim,ndim),gk(ndim),go(ndim)
  real(rfp)::st,so,zmut,vk,vo,ta,xwi,xw,xk,w3inv,vecn(ndim)
  real(rfp)::beta,betas,su,fb,fbs,resk,reso,rho,vol,rt,alpha
  real(rfp)::ksi=1.5_rfp,m_t0=0.0625_rfp,m_t,a   ! compressible correction

  integer::i,j,it,k,k1,k2
!
  
  if(s_ivis /= 2) return  ! no k-omega turbulent model
!
  do i = 1, size(cells)
!
  cc => cells(i)
  vol = cc%vol
!  do k = 1, ndim
!   gqv(k)%v = zero
!  end do
!  do k = 1, size(cc%sface)
!   cf => cc%sface(k)%to_face
!   qv = sum_qv_node(nodes(cf%f2n)) / real(size(cf%f2n),rfp)
!   vecn = cf%vecn * cf%area
!   do k2 = 1,neq
!   gqv%v(k2) = gqv%v(k2) + qv%v(k2) * vecn
!   end do 
!  end do
!  do k = 1, ndim
!   gqv(k)%v = gqv(k) / cc%vol
!  end do
!
  gvel = cc%gradient(imb:ime,1:ndim)  !du/dx
  gk   = cc%gradient(ike,1:ndim)  !dk/dxi
  go   = cc%gradient(iom,1:ndim)  !do/dxi
  qv   = cc%qv
  vk = abs(qv%v(ike))
  vo = abs(qv%v(iom))
!
  rho = rhof(qv)
  zmut = mu_turb(qv)
!
  su = zero
  do k = 1, ndim
  su = su + gvel(k,k)
  end do
  st = zero
  do k1 = 1, ndim
   do k2 = 1, ndim
    ta = zmut * (gvel(k1,k2) + gvel(k2,k1))
    if(k1 == k2) ta = ta - two_third * (zmut * su + rho * vk)
    st = st + ta * (gvel(k1,k2) + gvel(k2,k1)) * half
   end do
   end do
  xw = zero  
  xk = zero
  do k = 1, ndim
   do k1 = 1, ndim
    do k2 = 1, ndim
    ta = (gvel(k,k1) - gvel(k1,k))*(gvel(k1,k2) - gvel(k2,k1))
    ta = ta * (gvel(k2,k ) + gvel(k ,k2))
    xw = xw + ta
   end do
  end do
    xk = xk + gk(k) * go(k)
 end do
  w3inv = one / vo**3
  xw = 0.125_rfp * abs(xw * w3inv) / t_betas3
  xk = xk * w3inv
    
  fb = (one + 70.0_rfp * xw) / (one + 80.0_rfp * xw)

  xk = max(zero,xk)
  if(xk > 1.e10_rfp) then
   fbs = 1.7_rfp
  else
  xk = xk * xk
  fbs = (one + 680.0_rfp * xk) / (one + 400.0_rfp * xk)
  end if

! try Low-Reynolds-number version
!
!  rt = zmut / mu_suth(qv)
  
!  alpha = (0.025_rfp + rt / 6.0_rfp)/(one + rt / 6.0_rfp)
!  alpha = (0.1_rfp + rt/2.7_rfp)/(one + rt/2.7_rfp)/alpha
!  rt = (rt/8.0_rfp)**4
!  fbs = fbs *  (5.0_rfp/18.0_rfp + rt)/(one + rt)
!
!
!
  beta  = t_beta0  * fb
  betas = t_betas0 * fbs
!
! Compressible correction
!
  m_t = two * vk / sound_speed(qv)
  if(m_t > m_t0) then
   m_t = ksi * (m_t - m_t0)   ! ksi*(M_t^2-M_t0^2)
   beta = beta - betas * m_t
   betas = betas *( one + m_t)
  end if
!
  so = t_alpha * vo / vk * st
!  so = t_alpha * vo / vk * st * alpha  ! low-RE
  
  resk = st - betas * rho * vk * vo
  reso = so - beta * rho * vo * vo

  cc%res%v(ike) = cc%res%v(ike) + resk * vol
  cc%res%v(iom) = cc%res%v(iom) + reso * vol
!  LHS
!
  st = st * vol
!  so = so * vol
  rho = rho * vol

  su =  two_third * rho * su

  cc%dm%e(ike,ike) = cc%dm%e(ike,ike) + betas * rho * vo + st / vk
       
  cc%dm%e(ike,iom) = cc%dm%e(ike,iom) + betas * rho * vk - st / vo
!
  cc%dm%e(iom,iom) = cc%dm%e(iom,iom) + two * beta * vo * rho + su * t_alpha 

 end do
!
 end subroutine turb_model_1998


 function permutation(c,i,j,k)result(c1)
 integer::in(3,3,3),i,j,k
 real(rfp)::c,c1
 in = 0
 in(1,2,3) = 1;in(2,3,1) = 1; in(3,1,2) = 1;
 in(1,3,2) = -1;in(2,1,3) = -1; in(3,2,1) = -1;
 c1 = zero
 if(in(i,j,k) == 1) c1 = c
 if(in(i,j,k) == -1) c1 = -c
 end function permutation

!********************************************************
 subroutine turb_baldwin_lomax(cells,faces,interf)
!********************************************************
implicit none
!
  type(cell),pointer::cells(:),cc
  type(face),pointer::faces(:),cf
  type(wallface),pointer::wf
  type(itf)::interf
! local variables
  type(vector)::qv
  real(rfp)::gvel(ndim,ndim)
  real(rfp)::omega,y,yplus,dn,nu
  real(rfp)::rho,vel(ndim),u,lmix,f
  real(rfp)::fwake,delta,fkleb
  integer::i,j,k,k1,k2,i2bf,nf,n
  integer::id,ierr,nproc,myrank,loc
!  Closure Coefficients
  real(rfp)::kc = 0.40,alpha=0.0168,aplus0=26.0
  real(rfp)::ccp=1.6,ckleb=0.3,cwk=0.25   !1.0
!  work array
!
  real(rfp)::inbuff(2),outbuff(2)
  real(rfp),pointer::rbuff(:)
!
  if(s_ivis /= 3) return  ! no B-L turbulent model
!
  call mpi_comm_size(mpi_comm_world,nproc,ierr)
  call mpi_comm_rank(mpi_comm_world,id,ierr)

  do k = 1, nproc
   if(wallbc(k)%n == 0) cycle
   wallbc(k)%faces%ymax = mytiny
   wallbc(k)%faces%umax = mytiny
   wallbc(k)%faces%fmax = mytiny
  end do
!
  id = id + 1
    k = 0
  do i = 1, size(faces)
   cf => faces(i)
   if(cf%itype == 0) exit
   if(cf%itype >= partition_face_no) cycle
   if(bc(cf%itype)%igrp /= wall_type) cycle   ! no-wall
   if(bc(cf%itype)%itype == 0) cycle   ! no-wall
!
   cc => cf%left_cell
   vel = cc%qv%v(imb:ime)
   dn = sum(cf%vecn*(cf%centp - cc%centp))
   nu = mumix(cc%qv) / rhof(cc%qv)
   k = k + 1
   wallbc(id)%faces(k)%utawonu = sqrt(sqrt(sum(vel*vel)) / dn / nu)
  end do
! send to all processors
!
  do i = 1, nproc
    n = wallbc(i)%n
    if(n == 0) cycle
    allocate(rbuff(n))
    if(i == id ) rbuff = wallbc(i)%faces%utawonu
    call mpi_bcast(rbuff,n*rfp,mpi_byte,i - 1,mpi_comm_world,ierr)
    if(i /= id) wallbc(i)%faces%utawonu = rbuff
    deallocate(rbuff)
   end do
!
  do i = 1, size(cells)
!
   cc => cells(i)
   call find_rank_and_location(p2f_bl(i),nproc,myrank,loc)
   wf => wallbc(myrank)%faces(loc)
   qv = cc%qv
   vel = qv%v(imb:ime)
   u = sum(vel*vel)
   gvel = cc%gradient(imb:ime,1:ndim)  !du/dx
   rho = rhof(qv)
!
!Magnitude of the vorticity
   omega = zero
  do k1 = 1, ndim - 1
   do k2 = k1+1, ndim
    omega = omega + (gvel(k2,k1) - gvel(k1,k2))**2
   end do
  end do
   omega = sqrt(omega)
   y =  sum(wf%vecn * (wf%centp - cc%centp))
   if(y <= zero) y = 1.e10_rfp
   yplus = y * wf%utawonu
   lmix = kc * y * ( one - exp(-yplus / aplus0))
   zmut_bl(i) = rho * lmix * lmix * omega  
!
   F = lmix * omega / kc
! Detect the ymax
!
   if( f > wf%fmax) then
    wf%ymax = y
    wf%fmax = f
   end if
  if(u > wf%umax) wf%umax = u
!
 end do
!
!seek global maximum
!
 if(nproc < 0) then   ! skip this part, too cost for time
! if(nproc > 1) then
  do k = 1, nproc
   n = wallbc(k)%n
   if(n == 0) cycle
   do i = 1, n
   wf =>wallbc(k)%faces(i)
   call mpi_reduce(wf%umax,outbuff,1,mpi_double_precision,mpi_max,0,mpi_comm_world,ierr)
   call mpi_bcast(outbuff,rfp,mpi_byte,0,mpi_comm_world,ierr)
   wf%umax = outbuff(1)
!
   inbuff(1) = wf%fmax
   inbuff(2) = id
   call mpi_reduce(inbuff,outbuff,1,mpi_2double_precision,mpi_maxloc,0,mpi_comm_world,ierr)
   if(id == 1) j = int(outbuff(2))   ! location of max fmax
   call mpi_bcast(j,1,mpi_integer,0,mpi_comm_world,ierr)
! the max location
!
   if( id == j) then
    outbuff(1) = wf%ymax
    outbuff(2) = wf%fmax
   end if
   call mpi_bcast(outbuff,2*rfp,mpi_byte,j-1,mpi_comm_world,ierr)
   wf%ymax = outbuff(1)
   wf%fmax = outbuff(2)   
   end do
  end do
 end if
!
  do i = 1, size(cells)
!
  cc => cells(i)
   call find_rank_and_location(p2f_bl(i),nproc,myrank,loc)
   wf => wallbc(myrank)%faces(loc)
!
  y =  sum(wf%vecn * (wf%centp - cc%centp))
!
!Outer Layer
!
  fwake = wf%ymax * min(wf%fmax, cwk * wf%umax / wf%fmax)
  delta = y / wf%ymax
  if( delta > 1.e5_rfp) then
    fkleb = zero
  else
    fkleb = one / (one + 5.5_rfp *(ckleb * delta)**6)
  end if
  nu = rhof(cc%qv) * alpha * ccp * fwake * fkleb
  if(nu < zmut_bl(i) ) zmut_bl(i) = nu
 end do

 end subroutine turb_baldwin_lomax
 
  subroutine init_baldwin_lomax(cells,faces,interf)
  type(cell),pointer::cells(:),cc
  type(face),pointer::faces(:),cf
  type(wallface),pointer::wf
  type(itf)::interf
  integer::i,j,k,jmin,kmin,id,ierr,nb,nproc,n
  real(rfp)::dmin,dl
  real(rfp),pointer::rbuff(:)
  if(s_ivis /= 3) return
  ! search for no-slip wall boundary faces
  !
   call mpi_comm_size(mpi_comm_world,nproc,ierr)
   call mpi_comm_rank(mpi_comm_world,id,ierr)
   allocate(wallbc(nproc))
    nb = 0
  do i = interf%nitf + 1, size(faces)
    cf => faces(i)
    if(cf%itype == 0) exit
    if(cf%itype >= partition_face_no) cycle
    if(bc(cf%itype)%igrp /= wall_type) cycle
    if(bc(cf%itype)%itype == 0) cycle
    nb = nb + 1
  end do
   id = id + 1
    wallbc(id)%n = nb
  if(nb >  0) then
    allocate(wallbc(id)%faces(nb))
    j = 0
   do i = interf%nitf + 1, size(faces)
    cf => faces(i)
    if(cf%itype == 0) exit
    if(cf%itype >= partition_face_no) cycle
    if(bc(cf%itype)%igrp /= wall_type) cycle
    if(bc(cf%itype)%itype == 0) cycle
    j = j + 1
    wallbc(id)%faces(j)%centp = faces(i)%centp
    wallbc(id)%faces(j)%vecn  = faces(i)%vecn
   end do
  end if
  
  do i = 1, nproc
   if(i == id) n = wallbc(i)%n
   call mpi_bcast(n,1,mpi_integer,i - 1,mpi_comm_world,ierr)
    if(i /= id) then
     wallbc(i)%n = n
     if(n /= 0) allocate(wallbc(i)%faces(n))
    end if
    if(n == 0) cycle
    nb = ndim * n * 2
    allocate(rbuff(nb))
    if(i == id) then
     do j = 1, n
      k = (j -1)*ndim*2    
      rbuff(k+1:k+ndim) = wallbc(i)%faces(j)%centp
      rbuff(k+1+ndim:k+2*ndim) = wallbc(i)%faces(j)%vecn
     end do
    end if
    call mpi_bcast(rbuff,size(rbuff)*rfp,mpi_byte,i - 1,mpi_comm_world,ierr)

    if(i /= id) then
     do j = 1, n
      k = (j -1)*ndim*2    
      wallbc(i)%faces(j)%centp = rbuff(k+1:k+ndim) 
      wallbc(i)%faces(j)%vecn  = rbuff(k+1+ndim:k+2*ndim)
     end do
    end if
    deallocate(rbuff)
  end do   

  allocate(zmut_bl(size(cells)),p2f_bl(size(cells)))
  zmut_bl = zero
  p2f_bl  = 0
  
  do i = 1, size(cells)
    cc => cells(i)
    cc%zmut => zmut_bl(i)
    dmin = myhuge
    do j = 1, nproc
    do k = 1, wallbc(j)%n
     wf => wallbc(j)%faces(k)
     dl = sum((cc%centp - wf%centp)**2)
     if(dl < dmin) then
      jmin = j
      kmin = k
      dmin = dl
     end if
    end do
    end do
    p2f_bl(i) = nproc * kmin + jmin - 1
  end do
 end subroutine init_baldwin_lomax

 subroutine find_rank_and_location(is,np,myrank,loc)
  integer,intent(in)::is,np
  integer,intent(out)::myrank,loc
   myrank = mod(is,np) + 1
   loc    = is / np 
 end subroutine find_rank_and_location

end MODULE dfd_turb
