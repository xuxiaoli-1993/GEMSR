Program meshproducer
   implicit none
   integer::i,j,k,ijksize(1),totalijksize(1)
   integer::morei1,morei2,morej1,morej2
   real*8::dx(1),xyzsize(2),ratiox1,ratiox2,ratioy1,ratioy2
   real*8::maxdxi1,maxdxi2,maxdxj1,maxdxj2,dd
   real*8,dimension(:),allocatable::xc,yc,zc,xn,yn,zn
   character::ads*100

   !*********************************************************
   !parameters
   ijksize(1)=400
   morei1=0
   morei2=0
   totalijksize(1)=ijksize(1)+morei1+morei2

   xyzsize(1)=0.0
   xyzsize(2)=1.0

   ratiox1=1.15
   ratiox2=1.15

   allocate(xc(totalijksize(1)))
   allocate(xn(totalijksize(1)+1))

   do i=1,1
      dx(i)=(xyzsize(2*i)-xyzsize(2*i-1))/ijksize(i)
   end do

   maxdxi1=100*dx(1)
   maxdxi2=100*dx(1)

   xn(morei1+1)=xyzsize(1)
   do i=2,ijksize(1)+1
      xn(morei1+i)=xyzsize(1)+(i-1)*dx(1)
   end do

   dd=dx(1)
   do i=1,morei1
      dd=dd*ratiox1
      if(dd > maxdxi1) dd = maxdxi1
      xn(morei1+1-i)=xn(morei1+2-i)-dd
   end do

   dd=dx(1)
   do i=2,morei2+1
      dd=dd*ratiox2
      if(dd > maxdxi2) dd = maxdxi2
      xn(morei1+ijksize(1)+i)=xn(morei1+ijksize(1)+i-1)+dd
   end do

   ! change the boundary cells to normal size
   !xn(1) = xn(2)-dx(1)
   !xn(totalijksize(1)+1) = xn(totalijksize(1))+dx(1)

   print*,'xn:'
   do i=1,totalijksize(1)+1
      print*,xn(i)
   end do

   print*,'xc:'
   do i=1,totalijksize(1)
      xc(i)=(xn(i)+xn(i+1))/2
      print*,xc(i)
   end do

   ads='./1D_mesh.dat'
   open(1,file=trim(ads))
   write(1,*)totalijksize(1)+1,0,0
   do i=1,totalijksize(1)+1
      write(1,*) xn(i)
   end do

   close(1)

   !define boundary condition
   i=1
   j=1
   k=1
   ads='./1D_mesh_boundary.dat'
   open(1,file=trim(ads))
   write(1,*) size(xyzsize)
   write(1,*) i,i,0,0,0,0,1
   write(1,*) totalijksize(1)+1,totalijksize(1)+1,0,0,0,0,2
   close(1)

   ads='./tecplot_mesh.dat'
   open(1,file=trim(ads))
   write(1,*) 'TITLE="vector"'
   write(1,*) 'VARIABLES="x"'
   write(1,*) 'ZONE T="0" I=',totalijksize(1)+1,' F=POINT'
   do i=1,totalijksize(1)+1
      write(1,*) xn(i)
   end do
   close(1)

end Program meshproducer
