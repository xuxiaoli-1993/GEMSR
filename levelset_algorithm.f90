!**********************************************
!This program is for the stage 1 test of the development of LSGEMS
!Test the level-set advection equation and the re-initialization equation
!level-set advection (AD) equation: dphi/dt + u*dphi/dx + v*dphi/dy + w*dphi/dz = 0
!re-initialization (RI) equation: dphi/dt + sgn(phi0)*(abs(grad(phi))-1) = 0
!2nd-order ENO + 1st~2nd order TVD-RK
!major reference: "A Remark on Computing Distance Functions", "A review of level-set methods and some recent applications"
!01/08/2018
!**********************************************
Program Levelset
    implicit none
    real*8::PI=4*atan (1.0_8)
    integer::i,j,k,ii,jj,kk,n,ijksize(3),nRK_RI,nRK_AD,nit,Nt
    real*8::dx(3),dt_RI,cfl_RI,dt_AD,cfl_AD,xyzsize(6)
    real*8,dimension(:,:,:),allocatable::u,v,w,phi,phi0,phi_n,phi_n1,phi_n2,gradphi,ex_phi1,ex_phi2
    real*8,dimension(:),allocatable::xc,yc,zc
    real*8::bandsize_RI(2),it_err,it_err0,init_trial1_para(7),init_trial2_para(2),init_trial3_para(2,4),final_para(2,4)
    real*8::xsten(3,5),fsten(3,5),dd1(3,4),dd2(3,3),dphi(3,2),HJ
    character::plot_ads*100, charplot*5
    
    !*********************************************************
    !parameters
    ijksize(1)=100
    ijksize(2)=100
    ijksize(3)=100
    
    xyzsize(1)=0.
    xyzsize(2)=1.0
    xyzsize(3)=0.
    xyzsize(4)=1.0
    xyzsize(5)=0.
    xyzsize(6)=1.0
    
    cfl_RI=0.5
    cfl_AD=0.5
    plot_ads='./init_ls.dat'
    nRK_RI=2
    nRK_AD=2

    !initialization trial 1: zero level-set as sphere or ellipsoid
    !from parameter 1 - 7:
    !A,B,C,x0,y0,z0,epsilon
    init_trial1_para(1)=1.
    init_trial1_para(2)=1.
    init_trial1_para(3)=1.
    init_trial1_para(4)=1.
    init_trial1_para(5)=1.
    init_trial1_para(6)=1.
    init_trial1_para(7)=0.1
    
    !from parameter 1-2:
    !edge length of the cube, initial level-set value
    init_trial2_para(1)=0.75
    init_trial2_para(2)=5.

    do i=1,3
      dx(i)=(xyzsize(2*i)-xyzsize(2*i-1))/ijksize(i)
    end do
    
    bandsize_RI(1)=1000*maxval(dx)
    bandsize_RI(2)=5*maxval(dx)
    !*********************************************************

    allocate(u(ijksize(1),ijksize(2),ijksize(3)))
    allocate(v(ijksize(1),ijksize(2),ijksize(3)))
    allocate(w(ijksize(1),ijksize(2),ijksize(3)))
    allocate(phi0(ijksize(1),ijksize(2),ijksize(3)))
    allocate(phi(ijksize(1),ijksize(2),ijksize(3)))
    allocate(phi_n(ijksize(1),ijksize(2),ijksize(3)))
    allocate(phi_n1(ijksize(1),ijksize(2),ijksize(3)))
    allocate(phi_n2(ijksize(1),ijksize(2),ijksize(3)))
    allocate(gradphi(ijksize(1),ijksize(2),ijksize(3)))
    allocate(ex_phi1(ijksize(1),ijksize(2),ijksize(3)))
    allocate(ex_phi2(ijksize(1),ijksize(2),ijksize(3)))
    allocate(xc(ijksize(1)))
    allocate(yc(ijksize(2)))
    allocate(zc(ijksize(3)))
    
    print*,'xc:'
    do i=1,ijksize(1)
      xc(i)=xyzsize(1)+dx(1)/2+(i-1)*dx(1)
      print*,xc(i)
    end do
    
    print*,'yc:'
    do j=1,ijksize(2)
      yc(j)=xyzsize(3)+dx(2)/2+(j-1)*dx(2)
      print*,yc(j)
    end do
    
    print*,'zc:'
    do k=1,ijksize(3)
      zc(k)=xyzsize(5)+dx(3)/2+(k-1)*dx(3)
      print*,zc(k)
    end do
    
    !para(1,1:4): position and radius of the first sphere
    !para(2,1:4): position and radius of the second sphere
    !para(1,1:4)=para(2,1:4) means single sphere
    init_trial3_para(1,1)=0.35
    init_trial3_para(1,2)=0.35
    init_trial3_para(1,3)=0.35
    init_trial3_para(1,4)=0.15
    !init_trial3_para(1,4)=7.5*maxval(dx)
    init_trial3_para(2,1:4)=init_trial3_para(1,1:4)

    !initialization
    !call init_trial1(init_trial1_para,phi0,xc,yc,zc,ijksize)
    !call init_trial2(init_trial2_para,phi0,xc,yc,zc,ijksize,dx)
    call init_trial3(init_trial3_para,phi0,xc,yc,zc,ijksize)
    call plot_ls(phi0,xc,yc,zc,ijksize,plot_ads,0)
    stop


    call deformation_field(u,v,w,xc,yc,zc,ijksize)
    !call rotational_field(u,v,w,xc,yc,zc,ijksize)
    dt_AD=cfl_AD*minval(dx)/sqrt(maxval(dabs(u))**2+maxval(dabs(v))**2+maxval(dabs(w))**2)
    !0.6 is 1/5 of a period of the deformation field
    Nt=0.6/dt_AD
    print*,'dt_AD=',dt_AD
    print*,'Nt=',Nt
    read(*,*)
    
    final_para(1,1)=xc(25)
    final_para(1,2)=yc(13)-0.5*maxval(dx)
    final_para(1,3)=zc(25)
    final_para(1,4)=init_trial3_para(1,4)
    final_para(2,1:4)=init_trial3_para(1,1:4)
    
    call calculate_exact_ls1(final_para(1,1:4),ex_phi1,xc,yc,zc,ijksize) 
    call calculate_exact_ls1(final_para(2,1:4),ex_phi2,xc,yc,zc,ijksize) 

    phi=phi0
    phi_n1=phi
    open(2,file='./mass.dat')
    !call reinitialization(phi,xc,yc,zc,dx,ijksize,cfl_RI,nRK_RI,bandsize)

    !time evolution main loop
    do nit=1,2*Nt
	!reverse motion backwards 
      if(mod(nit,Nt)==1 .and. nit .ne. 1) then
	u=-u
	v=-v
	w=-w
      end if

      !Runge-Kutta loop
      do n=1,nRK_AD
	
	if(n==1) then
	  phi_n=phi
	end if
	
	!loop to update phi
	do k=3,ijksize(3)-2; do j=3,ijksize(2)-2; do i=3,ijksize(1)-2 
	
	    !find stencil
	    xsten(1,1)=xc(i-2)
	    xsten(1,2)=xc(i-1)
	    xsten(1,3)=xc(i)
	    xsten(1,4)=xc(i+1)
	    xsten(1,5)=xc(i+2)
	    xsten(2,1)=yc(j-2)
	    xsten(2,2)=yc(j-1)
	    xsten(2,3)=yc(j)
	    xsten(2,4)=yc(j+1)
	    xsten(2,5)=yc(j+2)
	    xsten(3,1)=zc(k-2)
	    xsten(3,2)=zc(k-1)
	    xsten(3,3)=zc(k)
	    xsten(3,4)=zc(k+1)
	    xsten(3,5)=zc(k+2)
	
	    fsten(1,1)=phi(i-2,j,k)
	    fsten(1,2)=phi(i-1,j,k)
	    fsten(1,3)=phi(i,j,k)
	    fsten(1,4)=phi(i+1,j,k)
	    fsten(1,5)=phi(i+2,j,k)
	    fsten(2,1)=phi(i,j-2,k)
	    fsten(2,2)=phi(i,j-1,k)
	    fsten(2,3)=phi(i,j,k)
	    fsten(2,4)=phi(i,j+1,k)
	    fsten(2,5)=phi(i,j+2,k)
	    fsten(3,1)=phi(i,j,k-2)
	    fsten(3,2)=phi(i,j,k-1)
	    fsten(3,3)=phi(i,j,k)
	    fsten(3,4)=phi(i,j,k+1)
	    fsten(3,5)=phi(i,j,k+2)

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
	    
	    !calcualte the upwind and downwind difference
	    !notations:
	    !dphi(1,1)=dxminus_phi
	    !dphi(1,2)=dxplus_phi
	    !dphi(2,1)=dyminus_phi
	    !dphi(2,2)=dyplus_phi
	    !dphi(3,1)=dzminus_phi
	    !dphi(3,2)=dzplus_phi
	    do ii=1,3
	      dphi(ii,1)=dd1(ii,2)+minmod(dd2(ii,1),dd2(ii,2))*(xsten(ii,3)-xsten(ii,2))
	      dphi(ii,2)=dd1(ii,3)+minmod(dd2(ii,2),dd2(ii,3))*(xsten(ii,3)-xsten(ii,4))
	    end do

	    !calcualte the Hamilton-Jacobian by upwind scheme
	    HJ=0	    
	    if(u(i,j,k)>0) then
	      HJ=HJ+u(i,j,k)*dphi(1,1)
	    else
	      HJ=HJ+u(i,j,k)*dphi(1,2)
	    end if
	    
	    if(v(i,j,k)>0) then
	      HJ=HJ+v(i,j,k)*dphi(2,1)
	    else
	      HJ=HJ+v(i,j,k)*dphi(2,2)
	    end if
	    
	    if(w(i,j,k)>0) then
	      HJ=HJ+w(i,j,k)*dphi(3,1)
	    else
	      HJ=HJ+w(i,j,k)*dphi(3,2)
	    end if
	      
	    !time evolution
	    if(n==1) then
	      phi_n1(i,j,k)=phi(i,j,k)-HJ*dt_AD
	    else if(n==2) then
	      phi_n2(i,j,k)=phi(i,j,k)-HJ*dt_AD
	      phi_n1(i,j,k)=0.5*(phi_n(i,j,k)+phi_n2(i,j,k))
	    end if
	
	!end of loop to update phi
	enddo; enddo; enddo
	
	phi=phi_n1
      !end of Runge-Kutta loop
      end do
      
      print*,'ADV ',nit
      if(mod(nit,10)==0) call reinitialization(phi,xc,yc,zc,dx,ijksize,cfl_RI,nRK_RI,bandsize_RI)
      if(mod(nit,10)==0) then
	write(2,*)nit,first_order_volume(phi,ijksize,dx)
      end if
      
      !monitor level-set 
      if(mod(nit,25)==0) then
	write(charplot,'(I5)') nit
	plot_ads='./plot2/ls_'//trim(adjustl(charplot))//'.dat'
	call plot_ls(phi,xc,yc,zc,ijksize,plot_ads,nit)
	!print*,'exact volume=',4*PI*init_trial3_para(1,4)**3/3.
	!print*,'actual volume=',first_order_volume(phi,ijksize,dx)
! 	if(mod(nit,2*Nt)==0) then
! 	  call calcualte_error(phi,ex_phi2,ijksize,5*maxval(dx),dx)
! 	else
! 	  call calcualte_error(phi,ex_phi1,ijksize,5*maxval(dx),dx)
! 	end if
      end if
    !end of time evolution main loop 
    enddo

!     write(charplot,'(I5)') nit
!     plot_ads='./plot/final_ls_adv_'//trim(adjustl(charplot))//'.dat'
!     call plot_ls(phi,xc,yc,zc,ijksize,plot_ads,nit)
! 
!     call calculate_exact_ls1(final_para,ex_phi,xc,yc,zc,ijksize) 
!     write(charplot,'(I5)') nit+1
!     plot_ads='./plot/final_ex_ls_'//trim(adjustl(charplot))//'.dat'
!     call plot_ls(ex_phi,xc,yc,zc,ijksize,plot_ads,nit+1)

    !estimate the error
    close(2)
    call calcualte_error(phi,ex_phi2,ijksize,10*maxval(dx),dx)
    

contains

!generate a rotating velocity field
subroutine rotational_field(u,v,w,x,y,z,sz)
    implicit none
    real*8,dimension(:,:,:),allocatable::u,v,w
    real*8,dimension(:),allocatable::x,y,z
    integer::sz(3),i,j,k

    do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
      u(i,j,k)=PI/314.*(yc(sz(2)/2)-yc(j))
      v(i,j,k)=PI/314.*(xc(i)-xc(sz(1)/2))
      w(i,j,k)=0.
    enddo; enddo; enddo
    
end subroutine rotational_field

!generate a deformation velocity field
subroutine deformation_field(u,v,w,x,y,z,sz)
    implicit none
    real*8,dimension(:,:,:),allocatable::u,v,w
    real*8,dimension(:),allocatable::x,y,z
    integer::sz(3),i,j,k

    do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
      u(i,j,k)=2*dsin(PI*x(i))*dsin(PI*x(i))*dsin(2*PI*y(j))*dsin(2*PI*z(k))
      v(i,j,k)=-dsin(2*PI*x(i))*dsin(PI*y(j))*dsin(PI*y(j))*dsin(2*PI*z(k))
      w(i,j,k)=-dsin(2*PI*x(i))*dsin(2*PI*y(j))*dsin(PI*z(k))*dsin(PI*z(k))
    enddo; enddo; enddo
    
end subroutine deformation_field

!make re-initialization as a subroutine
subroutine reinitialization(phi,xc,yc,zc,dx,ijksize,cfl_RI,nRK,bandsize)
    implicit none
    
    real*8,dimension(:,:,:),allocatable::phi,phi0,phi_n,phi_n1,phi_n2,gradphi,ex_phi
    real*8,dimension(:),allocatable::xc,yc,zc
    integer::i,j,k,ii,jj,kk,n,ijksize(3),nRK,nit
    real*8::xsten(3,5),fsten(3,5),intpx,dd1(3,4),dd2(3,3),dphi(3,2),fsten0(3,5),HJ
    real*8::dx(3),dt_RI,cfl_RI,bandsize(2)
    integer::qt_flag,ct_bd
    character::plot_ads*100, charplot*5
    real*8::a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2
    
    allocate(phi0(ijksize(1),ijksize(2),ijksize(3)))
    allocate(phi_n(ijksize(1),ijksize(2),ijksize(3)))
    allocate(phi_n1(ijksize(1),ijksize(2),ijksize(3)))
    allocate(phi_n2(ijksize(1),ijksize(2),ijksize(3)))
    allocate(gradphi(ijksize(1),ijksize(2),ijksize(3)))
    allocate(ex_phi(ijksize(1),ijksize(2),ijksize(3)))
    
    phi0=phi
    dt_RI=cfl_RI*minval(dx)/1.
    qt_flag=0
    nit=0
    phi_n1=phi
    !main while loop
    !do while(qt_flag .ne. 1)
    do while(qt_flag .ne. 1)
    !do while(qt_flag .ne. 1 .or. nit .le. 10)
      !Runge-Kutta loop
      do n=1,nRK
	
	if(n==1) then
	  phi_n=phi
	end if
	
	!loop to update phi
	!include all of the interior domain since the initial level-set is pretty distorted
	!only need to calcualte within the band if the level-set from the previous advection time step is already re-initialized
	do k=3,ijksize(3)-2; do j=3,ijksize(2)-2; do i=3,ijksize(1)-2 
	  
	  if(phi0(i,j,k) .le. 0.5*bandsize(1)) then
	    !find stencil
	    xsten(1,1)=xc(i-2)
	    xsten(1,2)=xc(i-1)
	    xsten(1,3)=xc(i)
	    xsten(1,4)=xc(i+1)
	    xsten(1,5)=xc(i+2)
	    xsten(2,1)=yc(j-2)
	    xsten(2,2)=yc(j-1)
	    xsten(2,3)=yc(j)
	    xsten(2,4)=yc(j+1)
	    xsten(2,5)=yc(j+2)
	    xsten(3,1)=zc(k-2)
	    xsten(3,2)=zc(k-1)
	    xsten(3,3)=zc(k)
	    xsten(3,4)=zc(k+1)
	    xsten(3,5)=zc(k+2)
	
	    fsten(1,1)=phi(i-2,j,k)
	    fsten(1,2)=phi(i-1,j,k)
	    fsten(1,3)=phi(i,j,k)
	    fsten(1,4)=phi(i+1,j,k)
	    fsten(1,5)=phi(i+2,j,k)
	    fsten(2,1)=phi(i,j-2,k)
	    fsten(2,2)=phi(i,j-1,k)
	    fsten(2,3)=phi(i,j,k)
	    fsten(2,4)=phi(i,j+1,k)
	    fsten(2,5)=phi(i,j+2,k)
	    fsten(3,1)=phi(i,j,k-2)
	    fsten(3,2)=phi(i,j,k-1)
	    fsten(3,3)=phi(i,j,k)
	    fsten(3,4)=phi(i,j,k+1)
	    fsten(3,5)=phi(i,j,k+2)

	    fsten0(1,1)=phi0(i-2,j,k)
	    fsten0(1,2)=phi0(i-1,j,k)
	    fsten0(1,3)=phi0(i,j,k)
	    fsten0(1,4)=phi0(i+1,j,k)
	    fsten0(1,5)=phi0(i+2,j,k)
	    fsten0(2,1)=phi0(i,j-2,k)
	    fsten0(2,2)=phi0(i,j-1,k)
	    fsten0(2,3)=phi0(i,j,k)
	    fsten0(2,4)=phi0(i,j+1,k)
	    fsten0(2,5)=phi0(i,j+2,k)
	    fsten0(3,1)=phi0(i,j,k-2)
	    fsten0(3,2)=phi0(i,j,k-1)
	    fsten0(3,3)=phi0(i,j,k)
	    fsten0(3,4)=phi0(i,j,k+1)
	    fsten0(3,5)=phi0(i,j,k+2)
	    
	    !if not using sub-cell fix
!   	   if(fsten0(1,2)*fsten0(1,3)<0 .or. fsten0(1,3)*fsten0(1,4)<0 .or. &
!   	     fsten0(2,2)*fsten0(2,3)<0 .or. fsten0(2,3)*fsten0(2,4)<0 .or. &
!   	     fsten0(3,2)*fsten0(3,3)<0 .or. fsten0(3,3)*fsten0(3,4)<0) then
!   
!   	    HJ=first_order_upwind_HJ(phi(i,j,k),fsten0,dx)
!   	    goto 100
!   	  end if
!   	  goto 50
	    !end of not using sub-cell
	    
	    !sub-cell fix
	    if(fsten0(1,2)*fsten0(1,3)<0) then
	      intpx=cubic_intp(xsten(1,1:4),fsten0(1,1:4))
	      xsten(1,1)=xc(i-1)
	      xsten(1,2)=intpx
	      fsten(1,1)=phi(i-1,j,k)
	      fsten(1,2)=0.
	    end if
	    
	    if(fsten0(1,3)*fsten0(1,4)<0) then
	      intpx=cubic_intp(xsten(1,2:5),fsten0(1,2:5))
	      xsten(1,4)=intpx
	      xsten(1,5)=xc(i+1)
	      fsten(1,4)=0.
	      fsten(1,5)=phi(i+1,j,k)
	    end if
	    
	    if(fsten0(2,2)*fsten0(2,3)<0) then
	      intpx=cubic_intp(xsten(2,1:4),fsten0(2,1:4))
	      xsten(2,1)=yc(j-1)
	      xsten(2,2)=intpx
	      fsten(2,1)=phi(i,j-1,k)
	      fsten(2,2)=0.
	    end if
	      
	    if(fsten0(2,3)*fsten0(2,4)<0) then
	      intpx=cubic_intp(xsten(2,2:5),fsten0(2,2:5))
	      xsten(2,4)=intpx
	      xsten(2,5)=yc(j+1)
	      fsten(2,4)=0.
	      fsten(2,5)=phi(i,j+1,k)
	    end if
	      
	    if(fsten0(3,2)*fsten0(3,3)<0) then
	      intpx=cubic_intp(xsten(3,1:4),fsten0(3,1:4))
	      xsten(3,1)=zc(k-1)
	      xsten(3,2)=intpx
	      fsten(3,1)=phi(i,j,k-1)
	      fsten(3,2)=0.
	    end if
	      
	    if(fsten0(3,3)*fsten0(3,4)<0) then
	      intpx=cubic_intp(xsten(3,2:5),fsten0(3,2:5))
	      xsten(3,4)=intpx
	      xsten(3,5)=zc(k+1)
	      fsten(3,4)=0.
	      fsten(3,5)=phi(i,j,k+1)
	    end if
	    50 continue
	    !end of sub-cell fix

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
	    
	    !calcualte the upwind and downwind difference
	    !notations:
	    !dphi(1,1)=dxminus_phi=a
	    !dphi(1,2)=dxplus_phi=b
	    !dphi(2,1)=dyminus_phi=c
	    !dphi(2,2)=dyplus_phi=d
	    !dphi(3,1)=dzminus_phi=e
	    !dphi(3,2)=dzplus_phi=f
	    do ii=1,3
	      dphi(ii,1)=dd1(ii,2)+minmod(dd2(ii,1),dd2(ii,2))*(xsten(ii,3)-xsten(ii,2))
	      dphi(ii,2)=dd1(ii,3)+minmod(dd2(ii,2),dd2(ii,3))*(xsten(ii,3)-xsten(ii,4))
	    end do

	    !calcualte the Hamilton-Jacobian by Godunov's scheme
	    !notations:
	    !a1=aplus=max(a,0)
	    !a2=aminus=min(a,0)
	    !b1=bplus=max(b,0)
	    !b2=bminus=min(b,0)
	    !c1=cplus=max(c,0)
	    !c2=cminus=min(c,0)
	    if(phi0(i,j,k)>0) then
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

	    !ill-advised sign distance function
	    !HJ=HJ*(phi0(i,j,k)/sqrt(phi0(i,j,k)**2+maxval(dx)**2))
	      
	    !use a smoothed sign distance function	  	  
	    if(fsten0(1,2)*fsten0(1,3)<0 .or. fsten0(1,3)*fsten0(1,4)<0 .or. &
	      fsten0(2,2)*fsten0(2,3)<0 .or. fsten0(2,3)*fsten0(2,4)<0 .or. &
	      fsten0(3,2)*fsten0(3,3)<0 .or. fsten0(3,3)*fsten0(3,4)<0) then
	      
	      HJ=HJ*smoothed_sign(fsten0,dx)
	    else
	      HJ=HJ*Ssign(phi0(i,j,k))
	    end if
	    
	    !time evolution
	    100 continue
	    if(n==1) then
	      phi_n1(i,j,k)=phi(i,j,k)-HJ*dt_RI
	    else if(n==2) then
	      phi_n2(i,j,k)=phi(i,j,k)-HJ*dt_RI
	      phi_n1(i,j,k)=0.5*(phi_n(i,j,k)+phi_n2(i,j,k))
	    end if
	  end if
	
	!end of loop to update phi
	enddo; enddo; enddo
	
	phi=phi_n1
      !end of Runge-Kutta loop
      end do
      
      !calculate quit flag of the while loop
      if(nit==1) it_err0=it_err

      it_err=0
      ct_bd=0
      do k=1,ijksize(3); do j=1,ijksize(2); do i=1,ijksize(1)
	if(dabs(phi_n(i,j,k)) .le. bandsize(2)/2) then
	  it_err=it_err+dabs(phi_n1(i,j,k)-phi_n(i,j,k))
	  ct_bd=ct_bd+1
	end if
      enddo; enddo; enddo
      it_err=it_err/ct_bd
	  
      !if(it_err/minval(dx) < 0.02) then
      !if(it_err < dt_RI*minval(dx)**2) then
      if(it_err/minval(dx)/dt_RI < 0.1) then
	qt_flag=1
      end if
      print*,'RI',nit,' error=',it_err
      nit=nit+1
      
      !monitor level-set 
!       if(mod(nit,1)==0) then
! 	write(charplot,'(I5)') nit
!	plot_ads='./plot/ls_'//trim(adjustl(charplot))//'.dat'
!	call plot_ls(phi,xc,yc,zc,ijksize,plot_ads,nit)
! 	if(nit==10) stop
!       end if
    !end of main while loop  

      !if can not convergence
      if(mod(nit,10)==0) then
	if(it_err>it_err0 .or. isnan(it_err) .or. dabs(maxval(phi))>10*(maxval(xc)-minval(xc))) then
	  !not advisable to do re-initialization anymore
	  print*,'re-initialization failed'
	  phi=phi0
	  goto 200
	end if	
      end if

      !force to quit
      if(nit==100) then
	phi=phi0
	print*,'re-initialization failed'
	goto 200
      end if
    end do
    
    200 continue
!     call fillin_levelset(phi,ijksize,2)
!     call calculate_gradient(gradphi,phi,xc,yc,zc,ijksize)
!     call calculate_exact_ls3(init_trial3_para,ex_phi,xc,yc,zc,ijksize)
!     write(charplot,'(I5)') nit
!     plot_ads='./plot/afterRI_ls_'//trim(adjustl(charplot))//'.dat'
!     call plot_ls_more1(phi,gradphi,ex_phi,xc,yc,zc,ijksize,plot_ads,nit)

    print*,'reinitialization finished'
!     write(charplot,'(I5)') nit
!     plot_ads='./plot/afterRI_ls_'//trim(adjustl(charplot))//'.dat'
!     call plot_ls(phi,xc,yc,zc,ijksize,plot_ads,nit)
    
    deallocate(phi0)
    deallocate(phi_n)
    deallocate(phi_n1)
    deallocate(phi_n2)
    deallocate(gradphi)
    deallocate(ex_phi)
    
end subroutine reinitialization

!calculate the maximum and average error betweent level-set and the exact level-set
!the error is in terms of dx
subroutine calcualte_error(p,ep,sz,bandsize,dx)
    implicit none
    integer::i,j,k,sz(3),ct
    real*8,dimension(:,:,:),allocatable::p,ep
    real*8::bandsize,dx(3),err,ave_err,max_err

    ave_err=0.
    max_err=0.
    ct=0

    do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
      if(dabs(ep(i,j,k)) .le. 0.5*bandsize) then
	err=dabs(p(i,j,k)-ep(i,j,k))
	if(err>max_err) max_err=err
	ave_err=ave_err+err
	ct=ct+1
      end if
    enddo; enddo; enddo
    ave_err=ave_err/ct

    print*,'average error=',ave_err/minval(dx)
    print*,'maximum error=',max_err/minval(dx)

end subroutine calcualte_error

!paranoid programing
subroutine fillin_levelset(p,sz,l)
    implicit none
    integer::i,j,k,ii,jj,kk,l,sz(3)
    real*8,dimension(:,:,:),allocatable::p
    
    do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
      ii=i
      jj=j
      kk=k
      if(i<1+l) ii=1+l
      if(i>sz(1)-l) ii=sz(1)-l
      if(j<1+l) jj=1+l
      if(j>sz(2)-l) jj=sz(2)-l
      if(k<1+l) kk=1+l
      if(k>sz(3)-l) kk=sz(3)-l
      
      p(i,j,k)=p(ii,jj,kk)
    enddo; enddo; enddo
end subroutine fillin_levelset

!*****************************************************************************************************************
!description of subroutine init_trial1:
!try to initialize phi as phi(x,y,z,0)=(epsilon + (x-x0)^2 + (y-y0)^2 + (z-z0)^2) * ((sqrt(x^2/A^2 + y^2/B^2 + z^2/C^2)) -1) 
!reference: "A Remark on Computing Distance Functions"
!input: para(7)=A,B,C,x0,y0,z0,epsilon, p=phi
!output: p
!turn out to be a good way to test re-initialization algorithms
!*****************************************************************************************************************
subroutine init_trial1(para,p,x,y,z,sz)
    implicit none    
    real*8,dimension(:,:,:),allocatable::p
    real*8,dimension(:),allocatable::x,y,z
    real*8::para(7)
    integer::i,j,k,sz(3)
    
    do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
      p(i,j,k)=(para(7)+(x(i)-para(4))**2+(y(j)-para(5))**2+(z(k)-para(6))**2)*(sqrt((x(i)/para(1))**2+(y(j)/para(2))**2+&
      (z(k)/para(3))**2)-1)      
    enddo; enddo; enddo
    
end subroutine init_trial1

!*****************************************************************************************************************
!description of subroutine calculate_exact_ls1:
!calculate the exact level set of a sphere
!input: para(4)=x0,y0,z0,r,ep=exact level-set,x=xc,y=yc,z=zc,sz=ijksize
!*****************************************************************************************************************
subroutine calculate_exact_ls1(para,ep,x,y,z,sz) 
    implicit none
    real*8,dimension(:,:,:),allocatable::ep
    real*8,dimension(:),allocatable::x,y,z
    real*8::para(4)
    integer::i,j,k,sz(3)
    
    do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
      ep(i,j,k)=sqrt((x(i)-para(1))**2+(y(j)-para(2))**2+(z(k)-para(3))**2)-para(4)
    enddo; enddo; enddo
end subroutine calculate_exact_ls1

!*****************************************************************************************************************
!description of subroutine init_trial2:
!try to initialize phi as a cube 
!reference: "Second-Order Accurate Computation of Curvatures in a Level Set Framework Using Novel High-Order Reinitialization Schemes"
!input: para(2)=(edge length,jump), p=phi
!in order to test sharp changes around zero level-set
!*****************************************************************************************************************
subroutine init_trial2(para,p,x,y,z,sz,dx)
    implicit none    
    real*8,dimension(:,:,:),allocatable::p
    real*8,dimension(:),allocatable::x,y,z
    real*8::para(2),dx(3),x0,y0,z0,eps
    integer::i,j,k,sz(3)
    
    do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
      if(dabs(x(i))>para(1)/2+dx(1)/3 .or. dabs(y(j))>para(1)/2+dx(2)/3 .or. dabs(z(k))>para(1)/2+dx(3)/3) then
	p(i,j,k)=para(2)
      else 
	p(i,j,k)=-para(2)
      end if
    enddo; enddo; enddo
    
    !try a distortion
!     x0=0.
!     y0=0.
!     z0=0.
!     eps=0.0
!     do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
!       p(i,j,k)=p(i,j,k)*((x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2+eps)
!     enddo; enddo; enddo
    
end subroutine init_trial2

!two spheres
subroutine init_trial3(para,p,x,y,z,sz)
    implicit none
    real*8,dimension(:,:,:),allocatable::p
    real*8,dimension(:),allocatable::x,y,z
    real*8::para(2,4),d(2),x0,y0,z0,eps
    integer::i,ii,j,k,sz(3)

!     do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
!       do ii=1,2
! 	d(ii)=sqrt((x(i)-para(ii,1))**2+(y(j)-para(ii,2))**2+(z(k)-para(ii,3))**2)-para(ii,4)
!       enddo
!       p(i,j,k)=minval(d)
!     enddo; enddo; enddo

    do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
      p(i,j,k)=sqrt((x(i)-0.35)**2+(y(j)-0.35)**2+(z(k)-0.35)**2)-0.15
    enddo; enddo; enddo

    !try a distortion
!      x0=1.
!      y0=1.
!      z0=1.
!      eps=0.1
!      do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
!        p(i,j,k)=p(i,j,k)*((x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2+eps)
!      enddo; enddo; enddo

end subroutine init_trial3

subroutine calculate_exact_ls3(para,ep,x,y,z,sz) 
    implicit none
    real*8,dimension(:,:,:),allocatable::ep
    real*8,dimension(:),allocatable::x,y,z
    real*8::para(2,4),d(2)
    integer::i,ii,j,k,sz(3)
    
    do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
      do ii=1,2
	d(ii)=sqrt((x(i)-para(ii,1))**2+(y(j)-para(ii,2))**2+(z(k)-para(ii,3))**2)-para(ii,4)
      enddo
      ep(i,j,k)=minval(d)
    enddo; enddo; enddo

end subroutine calculate_exact_ls3

!*****************************************************************************************************************
!description of subroutine plot_ls:
!plot the level-set using tecplot
!input: p=level set, sz=ijksize, ads=plot address, nt=time frame
!output: write the file to ads
!*****************************************************************************************************************
subroutine plot_ls(p,x,y,z,sz,ads,nt)
    implicit none
    real*8,dimension(:,:,:),allocatable::p
    real*8,dimension(:),allocatable::x,y,z
    integer::i,j,k,sz(3),nt
    character(100)::ads
    
    open(1,file=trim(ads))
    write(1,*) 'TITLE="vector"'
    write(1,*) 'VARIABLES="x" "y" "z" "levelset"'
    write(1,*) 'ZONE T="',nt,'" I=',sz(1),' J=',sz(2),' K=',sz(3),' F=POINT'
    print*, 'plotting levelset time = ',nt
    
    do k=1,sz(3)
        do j=1,sz(2)
            do i=1,sz(1)
                write(1,*)xc(i),yc(j),zc(k),p(i,j,k)
		!write(1,*)i,j,k,p(i,j,k)
            end do
        end do
   end do
   close(1)
   
end subroutine plot_ls

!plot also the gradient and relative error of the levels-set
subroutine plot_ls_more1(p,gp,ep,x,y,z,sz,ads,nt)
    implicit none
    real*8,dimension(:,:,:),allocatable::p,gp,ep
    real*8,dimension(:),allocatable::x,y,z
    integer::i,j,k,sz(3),nt
    character(100)::ads
    
    open(1,file=trim(ads))
    write(1,*) 'TITLE="vector"'
    write(1,*) 'VARIABLES="x" "y" "z" "levelset" "gradient" "error"'
    write(1,*) 'ZONE T="',nt,'" I=',sz(1),' J=',sz(2),' K=',sz(3),' F=POINT'
    print*, 'plotting levelset time = ',nt
    
    do k=1,sz(3)
        do j=1,sz(2)
            do i=1,sz(1)
                write(1,*)i,j,k,p(i,j,k),gp(i,j,k),dabs((p(i,j,k)-ep(i,j,k))/ep(i,j,k))
            end do
        end do
   end do
   close(1)
end subroutine plot_ls_more1

!*****************************************************************************************************************
!description of function ct_in_band:
!count the number of cells that have a level-set less than the band size
!to calculate the quit criterion
!input: p=level set, sz=ijksize, bdsz=band size
!output: M
!NOT IN USE
!*****************************************************************************************************************
function ct_in_band(p,sz,bdsz) result(M)
    implicit none
    real*8,dimension(:,:,:),allocatable::p
    integer::i,j,k,sz(3),M
    real*8::bdsz
    
    M=0
    do k=1,sz(3); do j=1,sz(2); do i=1,sz(1)
      if(dabs(p(i,j,k)) .le. bdsz/2) then
	M=M+1
      end if
    enddo; enddo; enddo
    
end function ct_in_band

!smoothed sign function calculated as level-set divided by the distance to the interfacec
!reference: "A Remark on Computing Distance Functions"
function smoothed_sign(ps0,dx) result(s)
    implicit none
    real*8::ps0(3,5),dx(3),s
    real*8::dpx,dpy,dpz,eps

    eps=0.1*minval(dx)
    
    !dpx=max(dabs(ps0(1,4)-ps0(1,2))/2.,dabs(ps0(1,3)-ps0(1,2)),dabs(ps0(1,4)-ps0(1,3)),eps)
    !dpy=max(dabs(ps0(2,4)-ps0(2,2))/2.,dabs(ps0(2,3)-ps0(2,2)),dabs(ps0(2,4)-ps0(2,3)),eps)
    !dpz=max(dabs(ps0(3,4)-ps0(3,2))/2.,dabs(ps0(3,3)-ps0(1,2)),dabs(ps0(3,4)-ps0(3,3)),eps)
    
    dpx=(ps0(1,4)-ps0(1,2))/2.
    dpy=(ps0(2,4)-ps0(2,2))/2.
    dpz=(ps0(3,4)-ps0(3,2))/2.
    s=ps0(3,3)/sqrt(dpx**2+dpy**2+dpz**2)
    
end function smoothed_sign

!minmod function
function minmod(a,b) result(r)
    implicit none  
    real*8::a,b,r
    
    if(dabs(a) .le. dabs(b) .and. a*b>0) then
      r=a
    else if(dabs(a) > dabs(b) .and. a*b>0) then
      r=b
    else
      r=0
    end if
end function minmod

!sign function
function Ssign(x) result(r)
    implicit none
    real*8::x,r
    if(x>0) then
      r=1.
    else if(x<0) then
      r=-1.
    else
      r=0.
    end if
end function Ssign

!*****************************************************************************************************************
!description of function cubic_intp:
!given p1,p2,p3 and p4, x(p1)=x1,x(p2)=x2,x(p3)=x3, and x(p4)=x4, and p2*p3<0
!find xA=x(0)
!reference: https://en.wikipedia.org/wiki/Spline_interpolation
!input: xs=(x1,x2,x3,x4),ps=(p1,p2,p3,p4)
!output: xA
!*****************************************************************************************************************
function cubic_intp(xs,ps) result(xA)
    implicit none
    real*8::xs(4),ps(4),xA,x1,x2,y1,y2,k1,k2,t,a,b
    
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

!linear interpolation
!not in use
function linear_intp(xs,ps) result(xA)
    implicit none
    real*8::xs(2),ps(2),xA
    
    xA=xs(1)+(xs(2)-xs(1))/(ps(2)-ps(1))*(-ps(1))
    
end function linear_intp

!*****************************************************************************************************************
!description of function first_order_upwind_HJ:
!use this if not using 2nd order sub-cell fix near the interface
!calculate the Hamilton-Jacobian S*(grad(phi)-1)
!reference: "A Remark on Computing Distance Functions"
!input: p=phi(i,j,k),ps0=stencil of phi0, dx=mesh size
!output: HJ
!*****************************************************************************************************************  
function first_order_upwind_HJ(p,ps0,dx) result(H)
    implicit none
    real*8::ps0(3,5),dx(3),p,H
    real*8::dphi(3),eps,D
    integer::i
    
    eps=0.05*minval(dx)
    do i=1,3      
      dphi(i)=max(dabs(ps0(i,4)-ps0(i,2))/2.,dabs(ps0(i,3)-ps0(i,2)),dabs(ps0(i,4)-ps0(i,3)),eps)/dx(i)
    end do
    
    D=ps0(3,3)/sqrt(dphi(1)**2+dphi(2)**2+dphi(3)**2)
    H=smoothed_sign(ps0,dx)*(dabs(p)/dabs(D)-1.)
end function first_order_upwind_HJ

subroutine calculate_gradient(gp,p,x,y,z,sz)
    implicit none
    real*8,dimension(:,:,:),allocatable::gp,p
    real*8,dimension(:),allocatable::x,y,z
    real*8::gx,gy,gz
    integer::i,j,k,sz(3)
    
    gp=1.
    !second-order central difference
    do k=2,sz(3)-1; do j=2,sz(2)-1; do i=2,sz(1)-1
      gx=(p(i+1,j,k)-p(i-1,j,k))/(x(i+1)-x(i-1))
      gy=(p(i,j+1,k)-p(i,j-1,k))/(y(j+1)-y(j-1))
      gz=(p(i,j,k+1)-p(i,j,k-1))/(z(k+1)-z(k-1))
      gp(i,j,k)=sqrt(gx**2+gy**2+gz**2)
    enddo; enddo; enddo

end subroutine calculate_gradient

!smeared heaviside function
function heaviside(p,eps) result(H)
    implicit none
    real*8::p,H,eps

    if(p<-eps) then
      H=0.
    else if(dabs(p) .le. eps) then
      H=0.5+0.5*p/eps+0.5*dsin(PI*p/eps)/PI
    else
      H=1.
    end if
end function heaviside      

End Program
