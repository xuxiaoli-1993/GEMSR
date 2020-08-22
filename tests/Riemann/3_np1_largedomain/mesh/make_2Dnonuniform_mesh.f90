Program meshproducer
   implicit none
   integer::i,j,k,ijksize(2),totalijksize(2)
   integer::morei1,morei2,morej1,morej2
   real*8::dx(2),xyzsize(4),ratiox1,ratiox2,ratioy1,ratioy2
   real*8::maxdxi1,maxdxi2,maxdxj1,maxdxj2,dd
   real*8,dimension(:),allocatable::xc,yc,zc,xn,yn,zn
   character::ads*100

   !*********************************************************
   !parameters
   ijksize(1)=400
   ijksize(2)=400
   morei1=20
   morei2=20
   morej1=20
   morej2=20
   totalijksize(1)=ijksize(1)+morei1+morei2
   totalijksize(2)=ijksize(2)+morej1+morej2

   xyzsize(1)=0.0
   xyzsize(2)=1.0
   xyzsize(3)=0.0
   xyzsize(4)=1.0

   ratiox1=1.2
   ratiox2=1.2
   ratioy1=1.2
   ratioy2=1.2

   allocate(xc(totalijksize(1)))
   allocate(yc(totalijksize(2)))
   allocate(xn(totalijksize(1)+1))
   allocate(yn(totalijksize(2)+1))

   do i=1,2
      dx(i)=(xyzsize(2*i)-xyzsize(2*i-1))/ijksize(i)
   end do

   maxdxi1=100*dx(1)
   maxdxi2=100*dx(1)
   maxdxj1=100*dx(2)
   maxdxj2=100*dx(2)

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

   yn(morej1+1)=xyzsize(3)
   do i=2,ijksize(2)+1
      yn(morej1+i)=xyzsize(3)+(i-1)*dx(2)
   end do

   dd=dx(2)
   do i=1,morej1
      dd=dd*ratioy1
      if(dd > maxdxj1) dd = maxdxj1
      yn(morej1+1-i)=yn(morej1+2-i)-dd
   end do

   dd=dx(2)
   do i=2,morej2+1
      dd=dd*ratioy2
      if(dd > maxdxj2) dd = maxdxj2
      yn(morej1+ijksize(2)+i)=yn(morej1+ijksize(2)+i-1)+dd
   end do

   ! change the boundary cells to normal size
   !yn(1) = yn(2)-dx(2)
   !yn(totalijksize(2)+1) = yn(totalijksize(2))+dx(2)

   print*,'yn:'
   do i=1,totalijksize(2)+1
      print*,yn(i)
   end do

   print*,'yc:'
   do i=1,totalijksize(2)
      yc(i)=(yn(i)+yn(i+1))/2
      print*,yc(i)
   end do

   ads='./2D_mesh.dat'
   open(1,file=trim(ads))
   write(1,*)totalijksize(1)+1,totalijksize(2)+1,0
   do j=1,totalijksize(2)+1; do i=1,totalijksize(1)+1
      write(1,*) xn(i)
   enddo; enddo

   do j=1,totalijksize(2)+1; do i=1,totalijksize(1)+1
      write(1,*) yn(j)
   enddo; enddo

   close(1)

   !define boundary condition
   i=1
   j=1
   k=1
   ads='./2D_mesh_boundary.dat'
   open(1,file=trim(ads))
   write(1,*) size(xyzsize)
   write(1,*) i,totalijksize(1)+1,j,j,0,0,1
   write(1,*) i,totalijksize(1)+1,totalijksize(2)+1,totalijksize(2)+1,0,0,2
   write(1,*) i,i,j,totalijksize(2)+1,0,0,3
   write(1,*) totalijksize(1)+1,totalijksize(1)+1,j,totalijksize(2)+1,0,0,4
   close(1)

   ads='./tecplot_mesh.dat'
   open(1,file=trim(ads))
   write(1,*) 'TITLE="vector"'
   write(1,*) 'VARIABLES="x(m)" "y(m)"'
   write(1,*) 'ZONE T="0" I=',totalijksize(1)+1,' J=',totalijksize(2)+1,' F=POINT'
   do j=1,totalijksize(2)+1
      do i=1,totalijksize(1)+1
         write(1,*) xn(i),yn(j)
      end do
   end do
   close(1)

end Program meshproducer
