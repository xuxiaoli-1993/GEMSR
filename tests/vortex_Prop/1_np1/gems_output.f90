module gems_output
   use gems_fv
   use gems_react
   use mpi
   implicit none
contains

   !*** write output in tecplot format
   !*************************************************************
   subroutine output_result(cells,nodes,faces,interf,id,it)
      !*************************************************************

      type(cell),pointer::cells(:)
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:)
      type(itf)::interf
      real(rfp)::t
      integer,pointer::c2n(:)
      integer::i,n,id,j,it,ierr
      character(len=40)::fn,fn2

      if(s_imonit == 0) return
      n = s_imonit*s_nt
      if(mod(it,n) /= 0) return

      j = it / n
      t = j * s_imonit / s_dt_inv

      fn = "./output/gems_mon_"//trim(i2s(j))//".plt"

      print *,id,'printing files'

      call output_cell_tecplot(fn,cells,nodes,faces,id)

      ! output the boundary data
      ! fn = "./output/gemsb."//trim(i2s(j))//".label"
      ! call wall_plot(fn,cells,nodes,faces,interf,id)
   end subroutine output_result

   subroutine output_cell_tecplot(fn,cells,nodes,faces,it)
      !********************************************************
      implicit none
      !
      type(cell),pointer::cells(:),cc,lcc,rcc,ucc,dcc
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:)

      ! local variables
      type(vector)::qv
      type(vector),allocatable::qvn(:)
      real(rfp)::rc(ndim),rn(ndim),mach,qinf,duc,dvc,uc,vc
      integer,pointer::c2n(:),nps(:)
      integer::i,j,k,n,it,nb
      character(len=*)::fn
      !  
      open(2,file=trim(fn)//'.'//trim(i2s(it))//'.dat')
      write(2,*)'variables ='
      if(ndim == 1) write(2,*)'x'
      if(ndim == 2) write(2,*)'x y'
      if(ndim == 3) write(2,*)'x y z'
      n = ndim
      do i = 1, size(output_selection)
         if(.not.output_selection(i)) cycle
         select case(i)
         case(20)
            do j = 1, nspe
               n = n + 1
               write(2,*)trim(species(j)%name)//'_mf'
            end do
         case(21)
            do j = 1, nspe
               n = n + 1
               write(2,*)trim(species(j)%name)//'_vf'
            end do
         case default
            n = n + 1
            write(2,*)trim(vname(i))
         end select
      end do

      ! head   
      if(ndim == 1) then
         write(2,*)'zone t="',s_elapsed_time,'",n=',size(nodes),',e=',size(cells),',ZONETYPE=FELINESEG',  &
            ',varlocation=([2-'//trim(i2s(n))//']=cellcentered)'
      else if(ndim == 2) then
         write(2,*)'zone t="',s_elapsed_time,'",n=',size(nodes),',e=',size(cells), &
            ',varlocation=([3-'//trim(i2s(n))//']=cellcentered)',      &
            ', f=feblock,et=quadrilateral'
      else if(ndim==3) then
         write(2,*)'zone t="',s_elapsed_time,'",n=',size(nodes),',e=',size(cells), &
            ',varlocation=([4-'//trim(i2s(n))//']=cellcentered)',      &
            ', f=feblock,et=quadrilateral'
      else
         print *, 'wrong dimension in output_cell_tecplot'; stop
      end if

      ! coordinates
      do i = 1, ndim
         write(2,20)nodes(:)%xyz(i)
      end do
      !
      n = size(cells)
      do j = 1, size(output_selection)
         if(.not.output_selection(j)) cycle
         select case(j)
         case(1)
            write(2,20)cells(:)%qv%v(ico)
         case(2)
            write(2,20)(pf(cells(i)%qv),i=1,n)
         case(3)
            write(2,20)(stagnation_pressure(cells(i)%qv,cells(i)%centp),i=1,n)
         case(5)
            write(2,20)cells(:)%qv%v(imb)
         case(6)
            write(2,20)cells(:)%qv%v(imb+1)
         case(7)
            if(naux > 0) then
               write(2,20)cells(:)%qv%v(iaub)
            else
               write(2,20)cells(:)%qv%v(imb+2)
            end if
         case(8)
            write(2,20)cells(:)%qv%v(ien)
         case(9)
            write(2,20)(stagnation_temperature(cells(i)%qv,cells(i)%centp),i=1,n)
         case(10)
            write(2,20)cells(:)%qv%v(ike)
         case(11)
            !    write(2,20)cells(:)%dqv%v(iom)
         case(12)
            write(2,20)(mu_turb(cells(i)%qv),i=1,n)
         case(13)
            write(2,20)(machf(cells(i)%qv),i=1,n)
         case(14)
            write(2,20)(sqrt(sound_speed(cells(i)%qv)),i=1,n)
         case(15)
            write(2,20)(rhof(cells(i)%qv),i=1,n)
         case(16)
            write(2,20)(hf(cells(i)%qv),i=1,n)
         case(17)
            write(2,20)(entropyf(cells(i)%qv),i=1,n)
         case(19) ! Vorticity
            do i=1,size(cells)
               cc => cells(i)
               write(2,20) fun_qv(cc,4)
            end do
         case(20)
            do i = isb,ise
               write(2,20)cells(:)%qv%v(i)
            end do
            write(2,20)((one - sum(cells(i)%qv%v(isb:ise))),i=1,n)
         case(21)
            do k = 1,nspe
               write(2,20)(volume_fraction(cells(i)%qv,k),i=1,n)
            end do
         case(23)
            write(2,20) (mumix(cells(i)%qv),i=1,n)
         case(24)
            write(2,20) (lamda(cells(i)%qv,cells(i)),i=1,n)
         case(25)
            write(2,20)cells(:)%qv%v(irad)
            !   case(26)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,1)  ! shadowgraph
         case(27)
            write(2,20)cells(:)%qv%v(ibb)
         case(28)
            write(2,20)cells(:)%qv%v(ibb+1)
         case(29)
            write(2,20)cells(:)%qv%v(ibe)
         case(30)
            write(2,20)cells(:)%qv%v(ieb)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,20+ieb) 
         case(31)
            write(2,20)cells(:)%qv%v(ieb+1)
         case(32)
            write(2,20)cells(:)%qv%v(iee)
         case(33)
            write(2,20)cells(:)%qv%v(ifb)
         case(34)
            write(2,20)cells(:)%qv%v(ife)
            !   case(35)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,5)  
            !   case(36)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,6) 
            !   case(37)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,7)
            !   case(38)  ! electrical conductivity
            !    write(2,20)nodevalue(acells,nodes,faces,interf,8)
            !   case(39)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,9)
            !   case(40)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,10)
            !   case(41)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,11)
            !   case(42)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,12)

         end select
      end do

      20  format(5e20.8)
      do i = 1, size(cells)
         c2n=>cells(i)%c2n
         n = size(c2n)
         select case(ndim)
         case(1)
            write(2,10)c2n
         case(2)
            if(n == 3) then
               write(2,10)c2n,c2n(3)
            else
               write(2,10)c2n
            end if
         case(3)
            if(n == 4) then
               write(2,10)c2n(:3),c2n(3),c2n(4),c2n(4),c2n(4),c2n(4)
            else if(n == 5) then
               write(2,10)c2n(:4),c2n(5),c2n(5),c2n(5),c2n(5)
            else if(n == 6) then
               write(2,10)c2n(:3),c2n(3),c2n(4:6),c2n(6)
            else
               write(2,10)c2n
            end if
         end select   
      end do

      close(2)
      10 format(10I8)
   end subroutine output_cell_tecplot

   subroutine output_to_plot(cells,nodes,faces,interf,it,iit)
      implicit none
      !
      type(cell),pointer::cells(:)
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:)
      type(itf)::interf
      real(rfp) :: t
      integer,intent(in)::it
      integer :: j,n,iit
      character(len=30)::fn

      n = s_imonit*s_nt
      j = iit / n
      t = j * s_imonit / s_dt_inv

      fn = './output/gems2plot_'//trim(i2s(j))

      select case(s_iplot)
      case(1) ! node based data
         call output_to_tecplot(fn,cells,nodes,faces,interf,it)
      case(2) 
         call output_to_fieldview(fn,cells,nodes,faces,it)
      case(3)
         call output_axisy_tecplot(fn,cells,nodes,faces,interf,it)
      case(4)  ! cell centered data
         call output_cell_tecplot(fn,cells,nodes,faces,it)
      case(5)  ! 2d to 3d only boundary surfaces output
         call output_axisy_surface(fn,cells,nodes,faces,interf,it)
      end select
      fn = "gemsb.label"
      call wall_plot(fn,cells,nodes,faces,interf,it)
      if(ivplot > 0) call integralvalue(faces,nodes,interf)

   end subroutine output_to_plot

   subroutine output_to_tecplot(fn,cells,nodes,faces,interf,it)
      !********************************************************
      implicit none
      !
      type(cell),pointer::cells(:),cc
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:)
      type(itf)::interf

      ! local variables
      type(vector)::qv
      type(vector),allocatable::qvn(:)
      real(rfp)::rc(ndim),rn(ndim),mach,qinf
      integer,pointer::c2n(:),nps(:)
      integer::i,j,k,n,it
      character(len=*)::fn
      !  
      open(2,file=trim(fn)//'.'//trim(i2s(it))//'.dat')
      write(2,*)'variables ='
      if(ndim == 1) write(2,*)'x'
      if(ndim == 2) write(2,*)'x y'
      if(ndim == 3) write(2,*)'x y z'
      do i = 1, size(output_selection)
         if(.not.output_selection(i)) cycle
         select case(i)
         case(20)
            do j = 1, nspe
               write(2,*)trim(species(j)%name)//'_mf'
            end do
         case(21)
            do j = 1, nspe
               write(2,*)trim(species(j)%name)//'_vf'
            end do
         case(22)
            do j = 1, nspe
               write(2,*)trim(species(j)%name)//'_c'
            end do
         case(54)
            do j = 1, nspe
               write(2,*)trim(species(j)%name)//'_rate'
            end do
         case default
            write(2,*)trim(vname(i))
         end select
      end do

      ! head   
      if(ndim == 2) then
         write(2,*)'zone t="',s_elapsed_time,'",n=',size(nodes),',e=',size(cells), &
            ', f=feblock,et=quadrilateral'
      else
         write(2,*) 'zone t="',s_elapsed_time,'",n=',size(nodes),',e=',size(cells),',f=feblock,et=brick'
         !write(2,*)'zone n=',size(nodes),',e=',size(cells),',f=feblock,et=brick'
      end if

      ! coordinates
      do i = 1, ndim
         write(2,20)nodes(:)%xyz(i) / b_lenref
      end do
      !
      n = size(nodes)
      do j = 1, size(output_selection)
         if(.not.output_selection(j)) cycle
         !write (2,*) 'printing',j

         select case(j)
         case(1)
            write(2,20)nodes(:)%qv%v(ico)
         case(2)
            write(2,20)(pf(nodes(i)%qv),i=1,n)
         case(3)
            write(2,20)(stagnation_pressure(nodes(i)%qv,nodes(i)%xyz),i=1,n)
         case(4)
            write(2,20)(pitot_pressure(nodes(i)%qv),i=1,n)
         case(5)
            write(2,20)nodes(:)%qv%v(imb)
         case(6)
            write(2,20)nodes(:)%qv%v(imb+1)
         case(7)
            if(naux > 0) then
               write(2,20)nodes(:)%qv%v(iaub)
            else
               write(2,20)nodes(:)%qv%v(imb+2)
            end if
         case(8)
            write(2,20)nodes(:)%qv%v(ien)
         case(9)
            write(2,20)(stagnation_temperature(nodes(i)%qv,nodes(i)%xyz),i=1,n)
         case(10)
            write(2,20)nodes(:)%qv%v(ike)
         case(11)
            !    write(2,20)nodes(:)%qv%v(iom)
         case(12)
            write(2,20)(mu_turb(nodes(i)%qv),i=1,n)
         case(13)
            write(2,20)(machf(nodes(i)%qv),i=1,n)
         case(14)
            write(2,20)(sqrt(sound_speed(nodes(i)%qv)),i=1,n)
         case(15)
            write(2,20)(rhof(nodes(i)%qv),i=1,n)
         case(16)
            write(2,20)(hf(nodes(i)%qv),i=1,n)
         case(17)
            write(2,20)(entropyf(nodes(i)%qv),i=1,n)
         case(18)
            write(2,20)nodevalue(cells,nodes,faces,interf,4)  ! vorticity
         case(19)
            write(2,20)nodevalue(cells,nodes,faces,interf,4)  ! vorticity
         case(20)
            do i = isb,ise
               write(2,20)nodes(:)%qv%v(i)
            end do
            write(2,20)((one - sum(nodes(i)%qv%v(isb:ise))),i=1,n)
         case(21)
            do k = 1,nspe
               write(2,20)(volume_fraction(nodes(i)%qv,k),i=1,n)
            end do
         case(22)
            do k = 1,nspe
               write(2,20)(concentration(nodes(i)%qv,k),i=1,n)
            end do
         case(23)
            write(2,20)(mumix(nodes(i)%qv),i=1,n)
         case(24)
            write(2,20)nodevalue(cells,nodes,faces,interf,13)  ! Thermal Conductivity
         case(25)
            write(2,20)nodes(:)%qv%v(irad)
         case(26)
            write(2,20)nodevalue(cells,nodes,faces,interf,1)  ! shadowgraph
         case(27)
            write(2,20)nodes(:)%qv%v(ibb)
         case(28)
            write(2,20)nodes(:)%qv%v(ibb+1)
         case(29)
            write(2,20)nodes(:)%qv%v(ibe)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,20+ibe) 
         case(30)
            write(2,20)nodes(:)%qv%v(ieb)
            !    write(2,20)nodevalue(cells,nodes,faces,interf,20+ieb) 
         case(31)
            write(2,20)nodes(:)%qv%v(ieb+1)
         case(32)
            write(2,20)nodes(:)%qv%v(iee)
         case(33)
            write(2,20)nodes(:)%qv%v(ifb)
         case(34)
            write(2,20)nodes(:)%qv%v(ife)
         case(35)
            write(2,20)nodevalue(cells,nodes,faces,interf,5)  
         case(36)
            write(2,20)nodevalue(cells,nodes,faces,interf,6) 
         case(37)
            write(2,20)nodevalue(cells,nodes,faces,interf,7)
         case(38)  ! electrical conductivity
            write(2,20)nodevalue(cells,nodes,faces,interf,8)
         case(39)
            write(2,20)nodevalue(cells,nodes,faces,interf,9)
         case(40)
            write(2,20)nodevalue(cells,nodes,faces,interf,10)
         case(41)
            write(2,20)nodevalue(cells,nodes,faces,interf,11)
         case(42)
            write(2,20)nodevalue(cells,nodes,faces,interf,12)
         case(43)
            write(2,20)(absorption_coef(nodes(i)%qv),i=1,n)
         case(44)  ! radiation heat flux 44-46
            write(2,20)nodevalue(cells,nodes,faces,interf,14)
         case(45)
            write(2,20)nodevalue(cells,nodes,faces,interf,15)
         case(46)
            write(2,20)nodevalue(cells,nodes,faces,interf,16)
         case(47)
            write(2,20)nodevalue(cells,nodes,faces,interf,17)
         case(48)  ! dBx/dt
            write(2,20)nodevalue(cells,nodes,faces,interf,18)
         case(49)  ! dBy/dt
            write(2,20)nodevalue(cells,nodes,faces,interf,19)
         case(50)  ! dBz/dt
            write(2,20)nodevalue(cells,nodes,faces,interf,20)
         case(51)  ! dEx/dt
            write(2,20)nodevalue(cells,nodes,faces,interf,21)
         case(52)  ! dEy/dt
            write(2,20)nodevalue(cells,nodes,faces,interf,22)
         case(53)  ! dEz/dt
            write(2,20)nodevalue(cells,nodes,faces,interf,23)
         case(54)
            do k = 1,nspe
               write(2,20)(write_reaction_rate(nodes(i)%qv,k),i=1,n)
            end do
         end select
      end do

      !write (2,*) 'second part'
      !  
      20  format(5e20.8)
      do i = 1, size(cells)
         c2n=>cells(i)%c2n
         n = size(c2n)
         select case(ndim)
         case(2)
            if(n == 3) then
               write(2,10)c2n,c2n(3)
            else
               write(2,10)c2n
            end if
         case(3)
            if(n == 4) then
               write(2,10)c2n(:3),c2n(3),c2n(4),c2n(4),c2n(4),c2n(4)
            else if(n == 5) then
               write(2,10)c2n(:4),c2n(5),c2n(5),c2n(5),c2n(5)
            else if(n == 6) then
               write(2,10)c2n(:3),c2n(3),c2n(4:6),c2n(6)
            else
               write(2,10)c2n
            end if
         end select   
      end do

      close(2)
      10 format(10I8)
   end subroutine output_to_tecplot

   !
   function nodevalue(cells,nodes,faces,interf,iv)result(p)
      type(cell),pointer::cells(:),cc
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:)
      type(itf)::interf
      real(rfp)::p(size(nodes)),f
      integer,pointer::c2n(:)
      integer::iv,n,i,j
      p = zero
      do i = 1, size(cells)
         c2n => cells(i)%c2n
         cc => cells(i)
         f = fun_qv(cc,iv)
         do j = 1, size(c2n)
            n = c2n(j)
            !   if(cc%weight(j) < 1.e-30_rfp) cycle
            p(n) = p(n) + f * cc%weight(j)
         end do
      end do
      !
      do i = 1, size(interf%pcell)
         cc => interf%pcell(i)
         c2n => cc%c2n
         f = fun_qv(cc,iv)
         do j = 1, size(c2n)
            n = c2n(j)
            !   if(cc%weight(j) < 1.e-30_rfp) cycle
            p(n) = p(n) + f * cc%weight(j)
         end do
      end do
      !
      do i = interf%nitf + 1, size(faces)
         if(faces(i)%itype == 0) exit
         cc => faces(i)%right_cell
         cc%gradient = faces(i)%left_cell%gradient
         f = fun_qv(cc,iv)
         c2n => faces(i)%f2n
         do j = 1, size(c2n)
            n = c2n(j)
            !   if(cc%weight(j) < 1.e-30_rfp) cycle
            p(n) = p(n) + f * cc%weight(j)
         end do
      end do

   end function nodevalue

   function fun_qv(cc,iv)result(f)
      type(cell),pointer::cc
      real(rfp)::dpdx(ndim),dtdx(ndim),drhodp,drhodt
      real(rfp)::f,dfdx(ndim),jb(3)
      integer::iv,i
      f = zero
      select case(iv) 
      case(1)    ! gradient of density
         drhodp = rhopf(cc%qv)
         drhodt = rhotf(cc%qv)
         dpdx = cc%gradient(ico,:)   
         dtdx = cc%gradient(ien,:)   
         dfdx = dpdx * drhodp + dtdx * drhodt
         f = sqrt(sum(dfdx**2))    ! Gradrho   Schlieren's image
      case(2)
         dpdx = cc%gradient(ico,:)   
         f = sqrt(sum(dpdx**2))    ! GradP   Schlieren's image
      case(3)
         dpdx = cc%gradient(ien,:)   
         f = sqrt(sum(dpdx**2))    ! GradT   Schlieren's image
      case(4)
         f = cc%gradient(ime,1) - cc%gradient(imb,2) ! vortex
      case(5)
         f = cc%jv(1)  !cc%qv%v(ieb) / e_resistivity(cc%qv,cc)
      case(6)
         f = cc%jv(2)  !cc%qv%v(ieb+1) / e_resistivity(cc%qv,cc)
      case(7)
         f = cc%jv(3)  !cc%qv%v(iee) / e_resistivity(cc%qv,cc)
      case(8)
         f = one / e_resistivity(cc%qv,cc)
      case(9)
         jb = acrossb(cc%jv,cc%qv%v(ibb:ibe))
         f = jb(1)
      case(10)
         jb = acrossb(cc%jv,cc%qv%v(ibb:ibe))
         f = jb(2)
      case(11)
         jb = acrossb(cc%jv,cc%qv%v(ibb:ibe))
         f = jb(3)
      case(12)
         f = sum(cc%jv * cc%qv%v(ieb:iee))
      case(13)
         f = lamda(cc%qv,cc)
      case(14)  ! Radiation Heat Flux:x
         f = - one_third / absorption_coef(cc%qv) * cc%gradient(irad,1)
      case(15)  ! Radiation Heat Flux:y
         f = - one_third / absorption_coef(cc%qv) * cc%gradient(irad,2)
      case(16)  ! Radiation Heat Flux:z
         f = - one_third / absorption_coef(cc%qv) * cc%gradient(irad,3)
      case(17)  ! Divergence E field
         f = diagsum(cc%gradient(ieb:iee,:))
      case(18)  ! dBx/dt
         do i = 2, s_idt + 1
            f = f + gear_c(i) * cc%qvn(gear_in(i-1))%v(ibb) 
         end do
         f = f + gear_c(1) * cc%qv%v(ibb)
      case(19)  ! dBy/dt
         do i = 2, s_idt + 1
            f = f + gear_c(i) * cc%qvn(gear_in(i-1))%v(ibb+1) 
         end do
         f = f + gear_c(1) * cc%qv%v(ibb+1)
      case(20)  ! dBz/dt
         do i = 2, s_idt + 1
            f = f + gear_c(i) * cc%qvn(gear_in(i-1))%v(ibe) 
         end do
         f = f + gear_c(1) * cc%qv%v(ibe)
      case(21)  ! dEx/dt
         do i = 2, s_idt + 1
            f = f + gear_c(i) * cc%qvn(gear_in(i-1))%v(ieb+1) 
         end do
         f = f + gear_c(1) * cc%qv%v(ieb)
      case(22)  ! dEy/dt
         do i = 2, s_idt + 1
            f = f + gear_c(i) * cc%qvn(gear_in(i-1))%v(ieb+1) 
         end do
         f = f + gear_c(1) * cc%qv%v(ieb+1)
      case(23)  ! dEz/dt
         do i = 2, s_idt + 1
            f = f + gear_c(i) * cc%qvn(gear_in(i-1))%v(iee) 
         end do
         f = f + gear_c(1) * cc%qv%v(iee)
      case(31:30+neq)
         f = cc%res%v(iv - 30)
      end select
   end function fun_qv


   subroutine output_to_fieldview(fn,cells,nodes,faces,it)
      !********************************************************
      implicit none
      !
      type(cell),pointer::cells(:),cc
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:)
      ! local variables
      type(vector)::qv
      type(vector),allocatable::qvn(:)
      real(rfp)::rc(ndim),rn(ndim),mach,qinf
      integer,pointer::c2n(:),nps(:)
      integer::i,j,n,it,nnodes
      character(len=*)::fn
      ! 
      !  treat no-slip wall
      !   n = 0
      !   do i = 1, size(faces)
      !    if(faces(i)%itype == 0) exit
      !    n = n + 1
      !    if(bc(faces(i)%itype)%igrp /= wall_type) cycle
      !    if(bc(faces(i)%itype)%itype /= 3) cycle
      !    do j = imb,ime
      !    nodes(faces(i)%f2n)%qv%v(j) = bc(faces(i)%itype)%vars(:ndim)
      !    end do
      !   end do

      !  open(2,file='dfd.fv')
      !  open(2,file='gems.fv.'//i2s(it))
      open(2,file=trim(fn)//'.'//trim(i2s(it))//'.fv')
      write(2,'(a)')'FIELDVIEW 2 5'
      write(2,'(a)')'Constants'
      write(2,*)0.0   !'Time'
      write(2,*)0.0  !FSMACH
      write(2,*)0.0   !ALPHA
      write(2,*)b_re  !Re

      write(2,'(a)')'GRIDS'
      write(2,*)1

      write(2,'(a)')'Boundary Table'
      write(2,*)7
      write(2,*) ' 0 0 1 Farfield'
      write(2,*) ' 0 0 1 Slip Wall'
      write(2,*) ' 0 0 1 Wall'
      write(2,*) ' 0 0 1 Moving Wall'
      write(2,*) ' 0 0 1 Inlet'
      write(2,*) ' 0 0 1 Outelet'
      write(2,*) ' 0 0 1 Heatwall'

      ! Houshang -begin

      write(2,'(a)')'Variable Names'
      write(2,*)neq + 3 - ndim
      !   write(2,*) 5            
      !   do i = 1,neq
      !   write(2,'(a1,i1)')'v',i
      !   end do
      write(2,'(a)')'pressure'
      write(2,'(a)')'u; velocity'
      write(2,'(a)')'v'
      if(ndim == 3) write(2,'(a)')'w'
      write(2,'(a)')'temperature'
      !   write(2,*)'sum_velocity'
      if(s_ivis == 2) then
         write(2,'(a)')'kin'
         write(2,'(a)')'omega'
      end if
      if(nspe > 1) then
         select case(s_irealgas)
         case(0)
            do i = 1, nspe
               write(2,'(1x,1hc,i1)')i
            end do
         case(2,3)
            do i = 1, nspe
               n = len_trim(species(i)%name)
               write(2,*)species(i)%name(:n)
            end do
         end select
      end if
      if(ndim == 2) write(2,'(a)')'dumy'
      !    write(2,*)'Mach Number'
      !    write(2,*)'Density'
      !    write(2,*)'Total pressure'

      write(2,'(a)')'Boundary variable Names'
      write(2,*)0


      write(2,'(a)')'Nodes'
      if(ndim == 2) then
         write(2,*)size(nodes) * 2
      else
         write(2,*)size(nodes)
      end if

      if(ndim ==2) then
         do i = 1, size(nodes)
            write(2,20)nodes(i)%xyz / b_lenref,zero
         end do
         do i = 1, size(nodes)
            write(2,20)nodes(i)%xyz / b_lenref,one
         end do
      else
         do i = 1, size(nodes)
            write(2,20)nodes(i)%xyz / b_lenref
         end do
      end if

      write(2,'(a)')'Boundary Faces'
      write(2,*)0
      !  write(2,*)n
      !  do i = 1, n
      !   j = faces(i)%itype
      !   if(j > 10) j = 5
      !   write(2,*)j,size(faces(i)%f2n),faces(i)%f2n
      !  end do

      write(2,'(a)')'Elements'
      !  write(2,*)size(cells)
      nnodes = size(nodes)
      do i = 1, size(cells)
         c2n=>cells(i)%c2n
         n = size(c2n)
         if(ndim == 2) then
            ! Houshang
            !     write(2,*)element_type(n,ndim),1,c2n,c2n + nnodes
            !    else
            !     write(2,*)element_type(n,ndim),1,c2n
            !SF -- begin
            if (size(c2n) == 4) then
               write(2,'(11i7)')element_type(n,ndim),1,c2n(1),c2n(2),c2n(4),c2n(3), &
                  c2n(1)+nnodes,c2n(2)+nnodes,c2n(4)+nnodes,c2n(3)+nnodes
            else if (size(c2n) == 3) then
               write(2,'(11i7)')element_type(n,ndim),1,c2n(1)+nnodes,c2n(1),c2n(3),&
                  c2n(3)+nnodes,c2n(2),c2n(2)+nnodes
            else
               write(2,'(11i7)')element_type(n,ndim),1,c2n,c2n + nnodes
            endif
         else
            if (size(c2n) == 8) then
               write(2,'(11i7)')element_type(n,ndim),1,c2n(1),c2n(5),c2n(2),c2n(6), &
                  c2n(4),c2n(8),c2n(3),c2n(7)
            else
               write(2,'(11i7)')element_type(n,ndim),1,c2n
            endif
            !SF - end
         end if
         !
      end do

      write(2,'(a)')'Variables'
      do i = 1, neq
         do j = 1, nnodes
            write(2,20)nodes(j)%qv%v(i)
         end do
         if(ndim == 2) then
            do j = 1, nnodes
               write(2,20)nodes(j)%qv%v(i)
            end do
         end if
         if(i == ime.and.ndim ==2) then
            do j = 1, 2*nnodes
               write(2,20) 0.0
            end do
         end if
      end do
      !   write(2,20)(machf(nodes(i)%qv),i=1,size(nodes))
      !   write(2,20)(rhof(nodes(i)%qv),i=1,size(nodes))

      !  if(it == 1) then
      !   write(2,20)(rhof(nodes(i)%qv),i=1,size(nodes))
      !   write(2,20)(p0f(nodes(i)%qv),i=1,size(nodes))
      !  end if
      !  qinf = half * b_rhoinf * b_velref**2
      !   write(2,20)((pf(nodes(i)%qv) - g_pbase) / qinf,i=1,size(nodes))


      write(2,'(a)')'Boundary Variables'

      20  format(5e20.8)
      close(2)
      10 format(10I8)
   end subroutine output_to_fieldview

   subroutine output_axisy_tecplot(fn,cells,nodes,faces,interf,it)
      !********************************************************
      implicit none
      !
      type(cell),pointer::cells(:),cc
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:)
      type(itf)::interf

      ! local variables
      type(vector)::qv
      type(vector),allocatable::qvn(:)
      real(rfp)::rc(ndim),rn(ndim),mach,qinf,theta,alpha
      integer,pointer::c2n(:),nps(:)
      integer::i,j,n,it,m,na,cn(4),k,numa,l
      character(len=*)::fn
      character(len=80)::str
      !  
      alpha = s_fspeed * pi / 180._rfp  ! borrow fspeed as rotating angls
      if(s_fspeed <= zero) alpha = two * pi
      open(2,file=trim(fn)//'.'//trim(i2s(it))//'.dat')
      write(2,*)'variables ='
      write(2,*)'x y z'
      numa = 3
      do i = 1, size(output_selection)
         if(.not.output_selection(i)) cycle
         numa = numa + 1
         select case(i)
         case(20)
            do j = 1, nspe
               write(2,*)trim(species(j)%name)//'_mf'
            end do
            numa = numa + nspm1
         case(21)
            do j = 1, nspe
               write(2,*)trim(species(j)%name)//'_vf'
            end do
            numa = numa + nspm1
         case(22)
            do j = 1, nspe
               write(2,*)trim(species(j)%name)//'_c'
            end do
            numa = numa + nspm1
         case default
            write(2,*)trim(vname(i))
         end select
      end do
      ! head   
      na = 30
      n = size(nodes)
      m = size(cells)
      write(2,*)'zone n=',n*na,',e=',m*(na-1),',f=feblock,et=brick'
      !
      ! coordinates

      do j = 1, na
         write(2,20)nodes(:)%xyz(1) / b_lenref
      end do
      do j = 1, na
         theta = real(j-1,rfp) * alpha / real(na-1,rfp)
         write(2,20)nodes(:)%xyz(2)  * cos(theta)/ b_lenref
      end do
      do j = 1, na
         theta = real(j-1,rfp) * alpha / real(na-1,rfp)
         write(2,20)nodes(:)%xyz(2)  * sin(theta)/ b_lenref
      end do
      !
      n = size(nodes)
      do j = 1, size(output_selection)
         if(.not.output_selection(j)) cycle
         do k = 1,na
            select case(j)
            case(1)
               write(2,20)nodes(:)%qv%v(ico)
            case(2)
               write(2,20)(pf(nodes(i)%qv),i=1,n)
            case(5)
               write(2,20)nodes(:)%qv%v(imb)
            case(6)  ! vy
               theta = real(k-1,rfp) * alpha / real(na-1,rfp)
               write(2,20)nodes(:)%qv%v(imb+1)*cos(theta) -nodes(:)%qv%v(iaub) * sin(theta)
            case(7)  ! vz
               theta = real(k-1,rfp) * alpha / real(na-1,rfp)
               write(2,20)nodes(:)%qv%v(imb+1)*sin(theta) +nodes(:)%qv%v(iaub) * cos(theta)
            case(8)
               write(2,20)nodes(:)%qv%v(ien)
            case(9)
               write(2,20)(stagnation_temperature(nodes(i)%qv,nodes(i)%xyz),i=1,n)
            case(10)
               write(2,20)nodes(:)%qv%v(ike)
            case(11)
               !    write(2,20)nodes(:)%qv%v(iom)
            case(12)
               write(2,20)(mu_turb(nodes(i)%qv),i=1,n)
            case(13)
               write(2,20)(machf(nodes(i)%qv),i=1,n)
            case(14)
               write(2,20)(sqrt(sound_speed(nodes(i)%qv)),i=1,n)
            case(15)
               write(2,20)(rhof(nodes(i)%qv),i=1,n)
            case(16)
               write(2,20)(hf(nodes(i)%qv),i=1,n)
            case(17)
               write(2,20)(entropyf(nodes(i)%qv),i=1,n)
            case(18)
               write(2,20)nodevalue(cells,nodes,faces,interf,4)  ! vorticity
            case(19)
               write(2,20)nodevalue(cells,nodes,faces,interf,4)  ! vorticity
               !   case(20)
               !    do i = isb,ise
               !    write(2,20)nodes(:)%qv%v(i)
               !    end do
               !    write(2,20)((one - sum(nodes(i)%qv%v(isb:ise))),i=1,n)
               !   case(21)
               !    do l = 1,nspe
               !    write(2,20)(volume_fraction(nodes(i)%qv,l),i=1,n)
               !    end do
               !   case(22)
               !    do l = 1,nspe
               !    write(2,20)(concentration(nodes(i)%qv,l),i=1,n)
               !    end do
            case(23)
               write(2,20)(mumix(nodes(i)%qv),i=1,n)
            case(24)
               write(2,20)nodevalue(cells,nodes,faces,interf,13)  ! Thermal Conductivity
            case(25)
               write(2,20)nodes(:)%qv%v(irad)
            case(26)
               write(2,20)nodevalue(cells,nodes,faces,interf,1)  ! shadowgraph
            case(27)
               write(2,20)nodes(:)%qv%v(ibb)
            case(28) ! By
               theta = real(k-1,rfp) * alpha / real(na-1,rfp)
               write(2,20)nodes(:)%qv%v(ibb+1)*cos(theta) -nodes(:)%qv%v(ibe) * sin(theta)
            case(29)  ! Bz
               theta = real(k-1,rfp) * alpha / real(na-1,rfp)
               write(2,20)nodes(:)%qv%v(ibb+1)*sin(theta) +nodes(:)%qv%v(ibe) * cos(theta)
            case(30)
               write(2,20)nodes(:)%qv%v(ieb)
            case(31)
               theta = real(k-1,rfp) * alpha / real(na-1,rfp)
               write(2,20)nodes(:)%qv%v(ieb+1)*cos(theta) -nodes(:)%qv%v(iee) * sin(theta)
            case(32)
               theta = real(k-1,rfp) * alpha / real(na-1,rfp)
               write(2,20)nodes(:)%qv%v(ieb+1)*sin(theta) +nodes(:)%qv%v(iee) * cos(theta)
            case(33)
               write(2,20)nodes(:)%qv%v(ifb)
            case(34)
               write(2,20)nodes(:)%qv%v(ife)
            case(35)
               write(2,20)nodevalue(cells,nodes,faces,interf,5)  
            case(36)
               write(2,20)nodevalue(cells,nodes,faces,interf,6) 
            case(37)
               write(2,20)nodevalue(cells,nodes,faces,interf,7)
            case(38)  ! electrical conductivity
               write(2,20)nodevalue(cells,nodes,faces,interf,8)
            case(39)
               write(2,20)nodevalue(cells,nodes,faces,interf,9)
            case(40)
               write(2,20)nodevalue(cells,nodes,faces,interf,10)
            case(41)
               write(2,20)nodevalue(cells,nodes,faces,interf,11)
            case(42)
               write(2,20)nodevalue(cells,nodes,faces,interf,12)
            case(43)
               write(2,20)(absorption_coef(nodes(i)%qv),i=1,n)
            case(44)  ! radiation heat flux 44-46
               write(2,20)nodevalue(cells,nodes,faces,interf,14)
            case(45)
               write(2,20)nodevalue(cells,nodes,faces,interf,15)
            case(46)
               write(2,20)nodevalue(cells,nodes,faces,interf,16)
            case(47)
               write(2,20)nodevalue(cells,nodes,faces,interf,17)
            end select
         end do
         !
         select case(j)
         case(20)
            do i = isb,ise
               do k = 1, na
                  write(2,20)nodes(:)%qv%v(i)
               end do
            end do
            do k = 1, na
               write(2,20)((one - sum(nodes(i)%qv%v(isb:ise))),i=1,n)
            end do
         case(21)
            do l = 1,nspe
               do k = 1,na
                  write(2,20)(volume_fraction(nodes(i)%qv,l),i=1,n)
               end do
            end do
         case(22)
            do l = 1,nspe
               do k = 1, na
                  write(2,20)(concentration(nodes(i)%qv,l),i=1,n)
               end do
            end do
         end select
      end do
      !  
      20  format(5e20.8)
      do i = 1, size(cells)
         c2n=>cells(i)%c2n
         m = size(c2n)
         cn(:m) = c2n
         if(m==3) cn(4) = c2n(3)
         do j = 1, na - 1
            write(2,10)cn,cn+n
            cn = cn + n
         end do
      end do
      str='d=('
      do i = 1, numa
         if(i > 1) str = trim(str)//','
         str = trim(str)//trim(i2s(i))
      end do
      str=trim(str)//')'

      do i = 1, size(bc)
         m = 0
         do j = 1, size(faces)
            if(faces(j)%itype == 0) exit
            if(faces(j)%itype == i) m = m + 1
         end do
         if(m == 0) cycle
         !  write(2,*)'zone n=',n*na,',e=',m*(na-1),',d=(1,2,3,4,5,6,7,8),f=feblock,et=quadrilateral'
         write(2,*)'zone t="','w'//trim(i2s(bc(i)%label))//'-'//trim(i2s(it)),'",'
         !  write(2,*)'zone n=',n*na,',e=',m*(na-1),trim(str)//',f=feblock,et=quadrilateral'
         write(2,*)'n=',n*na,',e=',m*(na-1),','
         write(2,*)trim(str)//',f=feblock,et=quadrilateral'
         do j = 1, size(faces)
            if(faces(j)%itype == 0) exit
            if(faces(j)%itype /= i) cycle
            c2n=>faces(j)%f2n
            cn(:2) = c2n
            do k = 1, na - 1
               write(2,10)cn(:2),cn(2:1:-1)+n
               cn(:2) = cn(:2) + n
            end do
         end do
      end do  ! end boundary
      close(2)
      10 format(10I8)
   end subroutine output_axisy_tecplot
   !  --------------------
   !
   subroutine output_axisy_surface(fn,cells,nodes,faces,interf,it)
      !********************************************************
      implicit none
      !
      type(cell),pointer::cells(:),cc
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:)
      type(itf)::interf

      ! local variables
      type(vector)::qv
      type(vector),allocatable::qvn(:)
      real(rfp)::rc(ndim),rn(ndim),mach,qinf,theta,alpha,alpha0
      integer,pointer::c2n(:),nps(:)
      integer::i,j,n,it,m,na,cn(4),k,numa,l
      character(len=*)::fn
      character(len=80)::str
      !  
      alpha0 = s_espeed * pi / 180._rfp
      alpha = (s_fspeed - s_espeed )* pi / 180._rfp  ! borrow fspeed as rotating angls
      if(s_fspeed <= zero) alpha = two * pi
      fn =trim(fn)//'.'//trim(i2s(it))//'.dat'
      open(2,file=trim(fn))
      write(2,*)'variables ='
      write(2,*)'x y z'
      numa = 3
      output_selection(17:) = .false.
      output_selection(14) = .false.
      do i = 1, size(output_selection)
         if(.not.output_selection(i)) cycle
         numa = numa + 1
         select case(i)
         case(20)
            do j = 1, nspe
               write(2,*)trim(species(j)%name)//'_mf'
            end do
            numa = numa + nspm1
         case(21)
            do j = 1, nspe
               write(2,*)trim(species(j)%name)//'_vf'
            end do
            numa = numa + nspm1
         case(22)
            do j = 1, nspe
               write(2,*)trim(species(j)%name)//'_c'
            end do
            numa = numa + nspm1
         case default
            write(2,*)trim(vname(i))
         end select
      end do
      ! head   
      na = 40
      n = size(nodes)
      m = size(cells)

      allocate(nps(n))
      nps = 0
      m = 0
      do i = 1, size(bc)
         do j = 1, size(faces)
            if(faces(j)%itype == 0) exit
            if(faces(j)%itype == i) then
               m = m + 1
               nps(faces(j)%f2n) = 1
            end if
         end do
      end do
      !
      n = 0
      do i = 1, size(nodes)
         if(nps(i) /= 0) then
            n = n + 1
            nps(n) = i
         end if
      end do
      if(n == 0) then
         close(2)
         call unlink(trim(fn))
         return
      end if
      write(2,*)'zone I=',n*na,',f=block'
      ! coordinates
      do j = 1, na
         write(2,20)nodes(nps(:n))%xyz(1) / b_lenref
      end do
      do j = 1, na
         theta = real(j-1,rfp) * alpha / real(na-1,rfp) + alpha0
         write(2,20)nodes(nps(:n))%xyz(2)  * cos(theta)/ b_lenref
      end do
      do j = 1, na
         theta = real(j-1,rfp) * alpha / real(na-1,rfp) + alpha0
         write(2,20)nodes(nps(:n))%xyz(2)  * sin(theta)/ b_lenref
      end do
      !
      !  n = size(nodes)
      do j = 1, size(output_selection)
         if(.not.output_selection(j)) cycle
         do k = 1,na
            select case(j)
            case(1)
               write(2,20)nodes(nps(:n))%qv%v(ico)
            case(2)
               write(2,20)(pf(nodes(nps(i))%qv),i=1,n)
            case(5)
               !    write(2,20)nodes(nps(:n))%qv%v(imb)
            case(6)  ! vy
               !    theta = real(k-1,rfp) * alpha / real(na-1,rfp) + alpha0
               !    write(2,20)nodes(nps(:n))%qv%v(imb+1)*cos(theta) -nodes(nps(:n))%qv%v(iaub) * sin(theta)
            case(7)  ! vz
               !    theta = real(k-1,rfp) * alpha / real(na-1,rfp) + alpha0
               !    write(2,20)nodes(nps(:n))%qv%v(imb+1)*sin(theta) +nodes(nps(:n))%qv%v(iaub) * cos(theta)
            case(8)
               write(2,20)nodes(nps(:n))%qv%v(ien)
            case(9)
               !    write(2,20)(stagnation_temperature(nodes(nps(i))%qv,nodes(i)%xyz),i=1,n)
            case(10)
               !    write(2,20)nodes(nps(:n))%qv%v(ike)
            case(11)
               !    write(2,20)nodes(nps(:n))%qv%v(iom)
            case(12)
               !    write(2,20)(mu_turb(nodes(nps(i))%qv),i=1,n)
            case(13)
               !    write(2,20)(machf(nodes(nps(i))%qv),i=1,n)
            case(15)
               write(2,20)(rhof(nodes(nps(i))%qv),i=1,n)
            case(16)
               write(2,20)(hf(nodes(nps(i))%qv),i=1,n)
            end select
         end do
      end do
      !
      20  format(5e20.8)
      !
      nps(n+1:) = 0
      do i=n,1,-1
         j = nps(i)
         if(j < i) print *,'something wrong'
         if(i/=j)then
            nps(j) = i
            nps(i) = 0
         end if
      end do
      !
      str='d=('
      do i = 1, numa
         if(i > 1) str = trim(str)//','
         str = trim(str)//trim(i2s(i))
      end do
      str=trim(str)//')'

      do i = 1, size(bc)
         m = 0
         do j = 1, size(faces)
            if(faces(j)%itype == 0) exit
            if(faces(j)%itype == i) m = m + 1
         end do
         if(m == 0) cycle
         write(2,*)'zone t="','w'//trim(i2s(bc(i)%label))//'-'//trim(i2s(it)),'",'
         write(2,*)'n=',n*na,',e=',m*(na-1),','
         write(2,*)trim(str)//',f=feblock,et=quadrilateral'
         do j = 1, size(faces)
            if(faces(j)%itype == 0) exit
            if(faces(j)%itype /= i) cycle
            c2n=>faces(j)%f2n
            cn(:2) = nps(c2n)
            do k = 1, na - 1
               write(2,10)cn(:2),cn(2:1:-1)+n
               cn(:2) = cn(:2) + n
            end do
         end do
      end do  ! end boundary
      close(2)
      deallocate(nps)
      10 format(10I8)
   end subroutine output_axisy_surface

   function element_type(n,nd)result(i)
      integer,intent(in)::n,nd
      integer::i
      select case(nd)
      case(2)  ! 2D
         if(n == 3) then
            i = 3
         else if(n == 4) then
            i = 2
         end if
      case(3)
         if(n == 4) then
            i = 1
         else if(n == 5) then
            i = 4
         else if(n == 6) then
            i = 3
         else if(n == 8) then
            i = 2
         end if
      end select
   end function element_type

   subroutine save_binary(n,cells,faces,nodes)
      type(cell),intent(in)::cells(:)
      type(face),intent(in)::faces(:)
      type(node),intent(in)::nodes(:)
      integer::n,k,id,j

      call mpi_comm_rank(mpi_comm_world,id,k)

      if(id==0) print*, 'save binary'

      open(5,file='gems.bin.'//i2s(id),form='unformatted')

      rewind(5)
      write(5) n ! n=nadv-1
      write(5) s_elapsed_time

      do k=1, size(cells)
         write(5) cells(k)%qv
      end do

      ! output values at boundary
      do k=1, size(faces)
         if(faces(k)%itype == 0) exit
         write(5) faces(k)%right_cell%qv
      end do
      !
      do j=1, s_idt
         do k=1, size(cells)
            write(5) cells(k)%qvn(gear_in(j))
         end do
      end do

      close(5)
   end subroutine save_binary

   subroutine integralvalue(faces,nodes,interf)
      type(face),pointer::faces(:),cf
      type(node),pointer::nodes(:)
      real(rfp)::vecn(ndim),un,rhou,area,c1,c2,c3,g,un1,g1
      real(rfp)::t0inlet,tinlet,pinlet,p0inlet
      real(rfp)::t0outlet,toutlet,poutlet,p0outlet
      real(rfp),pointer::vars(:,:),rbuff(:),tarea(:)
      type(vector)::qv,qvw
      type(itf)::interf
      integer::i,nb,it,n,j,ierr,id,nv,is,ndim1,md,inlet,outlet
      character(len=20)::vars_name(16),bname(6)
      vars_name(:16) = (/&
         'Mass (kg/s)          ',  &
         'Average P0 (N/m^2s)  ',  &
         'Average T0 (K)       ',  &
         'Momentum x(N/s)      ',  &
         'Momentum y(N/s)      ',  &
         'Momentum z(N/s)      ',  &
         'Force x (N)          ',  &
         'Force y (N)          ',  &
         'Force z (N)          ',  &
         'Power (W)            ',  &
         'SkinFriction x (N)   ',  &
         'SkinFriction y (N)   ',  &
         'SkinFriction z (N)   ',  &
         'Heat Flux (W)        ',  &
         'ElectricCurrent (A)  ',  &
         'RadiationFlux (W)    '   &
         /)
      bname   =(/'Inlet','Outlet','Farfield','Wall','GeometryWall','ElecMagBound'/)
      md = 3
      nv = 3*md + 6  ! total number of output
      if(nrad > 0) nv = nv + 1
      !  if(neqm > 0) nv = nv + 1
      ndim1 = ndim - 1
      call mpi_comm_rank(mpi_comm_world,id,ierr)
      nb = size(bc)
      allocate(rbuff(nb),vars(nv,nb),tarea(nb))
      vars = zero; tarea = zero   ! total area
      do i = interf%nitf + 1, size(faces)
         cf => faces(i)
         it = faces(i)%itype
         if(it == 0) exit
         if(it >= partition_face_no) cycle
         area = faces(i)%area
         if(s_iaxis > 0.and.ndim==2) area = two * pi * area
         tarea(it) = tarea(it) + area
         vecn = faces(i)%vecn * area
         ! face value
         call face_qv_tauw(cf,nodes,qv,qvw)
         g = htf(qv) / cvf(qv); g1 = ( g - one) / g
         un   = dot_product(vecn,qv%v(imb:ime))
         rhou = rhof(qv) * un
         !  Mass flux   1
         is = 1
         vars(is,it) = vars(is,it) + rhou
         !  Stagnation Pressure
         is = is + 1
         vars(is,it) = vars(is,it) + stagnation_pressure(qv,cf%centp)**g1 * rhou
         !    vars(is,it) = vars(is,it) + rhou * un / area + pf(qv) * area
         !  Stagnation Temperature
         is = is + 1
         vars(is,it) = vars(is,it) + stagnation_temperature(qv,cf%centp) * rhou
         !    vars(is,it) = vars(is,it) + rhou * h0f(qv) 
         !  Momentum flux !
         is = is + 1
         vars(is:is+ndim1,it) = vars(is:is+ndim1,it) + rhou * qv%v(imb:ime)
         !  Force on boundary surface
         is = is + md
         vars(is:is+ndim1,it) = vars(is:is+ndim1,it) + pf(qv) * vecn
         !  Energy
         is = is + md
         vars(is,it) = vars(is,it) +  rhou * h0f(qv)
         !  Skin friction
         is = is + 1
         vars(is:is+ndim1,it) = vars(is:is+ndim1,it) + qvw%v(imb:ime) * area
         !  Heat FLux
         is = is + md
         vars(is,it) = vars(is,it) +  face_value(cf,qv,qvw,7) * area
         ! Electric Current
         is = is + 1
         if(neqm > 0) then
            vars(is,it) = vars(is,it) +  sum(cf%left_cell%jv(:ndim) * vecn)
         end if
         ! Radiation Heat Flux
         if(nrad > 0) then
            is = is + 1
            vars(is,it) = vars(is,it) +  face_value(cf,qv,qvw,11) * area
         end if
      end do
      ! total area
      rbuff = zero 
      call mpi_reduce(tarea,rbuff,nb,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
      tarea = rbuff
      !
      do i = 1, nv
         rbuff = zero 
         call mpi_reduce(vars(i,:),rbuff,nb,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
         vars(i,:) = rbuff
      end do
      if(id == 0 ) then
         open(10,file='gems_integration.dat')
         write(10,*)'Integration Values on the Bounddary Surface'
         write(10,*)'Total has ',nb, ' surface patchs'
         do i = 1,nb
            if(tarea(i) <= mytiny) cycle
            write(10,*)'------ Label =',bc(i)%label, ' and Btype=',bname(bc(i)%igrp),'-------'
            write(10,*)'Surface Area (m^2)   ', tarea(i)
            vars(2,i) = vars(2,i) / vars(1,i)   ! average value
            vars(3,i) = vars(3,i) / vars(1,i)   ! average value
            !    vars(3,i) = vars(3,i) / tarea(i)   ! average value
            if(bc(i)%igrp == inlet_type) inlet = i
            if(bc(i)%igrp == outlet_type) outlet = i
            !
            do j = 1, nv
               if(j == 2.or.j==3) cycle
               write(10,'(a,e20.6)')vars_name(j)//'= ',vars(j,i)
            end do
         end do
         write(10,*) ' ================================================'
         !  inlet
         !   c1 = abs(vars(1,inlet)) / tarea(inlet)
         !   c2 = abs(vars(2,inlet))
         !   c3 = abs(vars(3,inlet))
         !   g = htf(qv) / cvf(qv)
         !   un = (c2 - sqrt(c2*c2 - 2._rfp * (g*g-one)*c3*c1/(g*g))) / &
         !         ( ( g + one) * c1 / g)
         !   t0inlet = c3 / c1 / htf(qv)
         !!   tinlet  = t0inlet - half * un * un / htf(qv)   wrong!!!
         !   pinlet  = c2 - c1 * un
         !   p0inlet = pinlet * ( t0inlet / tinlet )**(g/(g-one))
         p0inlet = abs(vars(2,inlet))** (one / g1)
         t0inlet = vars(3,inlet)
         !!
         !   tinlet = pinlet / (c1/un * g_runi / g_mwi(1))  ! t=p/Rrho
         !   un1 = un
         !! outlet
         !   c1 = abs(vars(1,outlet)) / tarea(outlet)
         !   c2 = abs(vars(2,outlet))
         !   c3 = abs(vars(3,outlet))
         !!   g = htf(qv) / cvf(qv)
         !   un = (c2 - sqrt(c2*c2 - 2._rfp * (g*g-one)*c3*c1/(g*g))) / &
         !         ( ( g + one) * c1 / g)
         !   t0outlet = c3 / c1 / htf(qv)
         !!   toutlet  = t0outlet - half * un * un / htf(qv)
         !   poutlet  = c2 - c1 * un
         !   p0outlet = poutlet * ( t0outlet / toutlet )**(g/(g-one))
         p0outlet = abs(vars(2,outlet))**(one/g1)
         t0outlet = vars(3,outlet)
         !   toutlet = poutlet / (c1/un * g_runi / g_mwi(1))  ! t=p/Rrho
         !!
         !!   un =  vars(2,outlet) / vars(2,inlet)
         write(10,*)'Ratio of Specific Heats',htf(qv) / cvf(qv)
         write(10,*)'Rotating Speed: rpm',s_omega / pi * 30._rfp
         write(10,*)'                   Inlet              Outlet'
         !   write(10,*)' Velocity      :', un1, un
         write(10,*)' Mass Flow Rate:', vars(1,inlet)*b_nblade,vars(1,outlet)*b_nblade
         write(10,*)' Total Pressure:', p0inlet,p0outlet
         write(10,*)' Total Temperature:', t0inlet,t0outlet
         !   write(10,*)' Static Pressure:',pinlet,poutlet
         !   write(10,*)' Static Temperature:',tinlet,toutlet
         un = p0outlet / p0inlet
         write(10,*) 'Pressure ratio =', un
         un = un**g1
         !   rhou = vars(3,outlet) / vars(3,inlet)
         rhou = t0outlet / t0inlet
         write(10,*) 'Efficiency =', (un - one ) / ( rhou - one)
         write(10,*)' Corrected mass flow rate:WRT/P(W: lb/s, T: K P: Psi)'
         write(10,*)'Inlet:',abs(vars(1,inlet)) * kg2lb * sqrt(t0inlet) / (p0inlet / psi2pa) * b_nblade
         write(10,*)'Outlet:',abs(vars(1,outlet)) * kg2lb * sqrt(t0inlet) / (p0inlet / psi2pa) * b_nblade
         write(10,*)' N/sqrt(T) =', s_omega * 30._rfp / pi / sqrt(t0inlet)
         close(10)
      end if
      !
      deallocate(vars,rbuff,tarea)
   end subroutine integralvalue

   subroutine massflow(faces,nodes,id,interf)
      type(face)::faces(:)
      type(node)::nodes(:)
      real(rfp)::massin,massout,vecn(ndim),p0,pn,fi,mtin,mtout
      type(vector)::qv
      type(itf)::interf
      integer::i,nin,nout,it,n,j,ierr,id
      massin = zero
      massout = zero
      p0 = zero
      pn = zero
      fi = zero
      nin = 0
      nout = 0
      !  do i = 1, size(faces)
      do i = interf%nitf + 1, size(faces)
         it = faces(i)%itype
         if(it == 0) exit
         if(it >= partition_face_no) cycle
         select case(bc(it)%igrp)
         case(inlet_type)
            vecn = faces(i)%vecn * faces(i)%area
            qv%v = zero
            n = size(faces(i)%f2n)
            do j = 1, n
               qv%v = qv%v + nodes(faces(i)%f2n(j))%qv%v
            end do
            qv%v = qv%v / real(n,rfp)
            nin = nin + 1
            p0 = p0 + pf(qv)
            fi = fi + sqrt(rhof(qv) * pf(qv)) 
            massin = massin - rhof(qv) * dot_product(vecn,qv%v(imb:ime))
            !
         case(outlet_type)
            vecn = faces(i)%vecn * faces(i)%area
            qv%v = zero
            n = size(faces(i)%f2n)
            do j = 1, n
               qv%v = qv%v + nodes(faces(i)%f2n(j))%qv%v
            end do
            qv%v = qv%v / real(n,rfp)
            nout = nout + 1
            pn = pn + pf(qv)
            fi = fi + sqrt(rhof(qv) * pf(qv)) 
            massout = massout + rhof(qv) * dot_product(vecn,qv%v(imb:ime))
         end select
      end do

      !   if(nin == 0.or.nout==0) return   ! nothing
      call mpi_reduce(massin,mtin,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
      call mpi_reduce(massout,mtout,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
      if(id == 0) then
         open(1,file='gems.mass')
         !   p0 = p0 / real(nin,rfp)
         !   pn = pn / real(nout,rfp)
         !   fi = fi / real(nin,rfp) * 0.001016**2 * 0.1 !0.36 * 0.001  !0.44 * 0.0254
         write(1,*)mtin,mtout
         close(1)
      end if
   end subroutine massflow

   subroutine force_on_wall(force,faces,interf)
      type(face)::faces(:)
      real(rfp)::force(ndim),vecn(ndim),coef
      type(cell),pointer::cr
      type(itf)::interf
      integer::i,it
      force = zero
      !  do i = 1, size(faces)
      do i = interf%nitf + 1, size(faces)
         it = faces(i)%itype
         if(bc(it)%igrp /= wall_type) cycle
         vecn = faces(i)%vecn * faces(i)%area

         cr => faces(i)%right_cell
         force = force + cr%qv%v(ico) * vecn
      end do

      !  coef = half * b_rhoinf * sum(b_velinf * b_velinf)
      force = force   !/ coef
   end subroutine force_on_wall


   subroutine thickness_bl(faces,nodes,ip,interf)
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:),ff
      real(rfp)::s,su,s1,s2,s3,g,dn
      real(rfp)::u1,u2,u0,rho,rho1,rho2
      real(rfp)::p(ndim),p1(ndim),p2(ndim),vecn(ndim)
      type(itf)::interf
      type(cell),pointer::cc
      type(vector)::qv
      integer::i,ip,j
      if(ndim /= 2) return  ! only for 2D
      if(ip == 0) then
         open(8,file='wthickness.dat')
         write(8,*)'variables ='
         write(8,*)'x y s0 su s1 s2 s3 x0 y0 x1 y1'
         write(8,*)'zone t="thickness",f="point"'
      else
         open(8,file='wthickness.'//i2s(ip))
      end if
      do i = interf%nitf + 1, size(faces)
         if(faces(i)%itype == 0) exit
         if(bc(faces(i)%itype)%igrp /= wall_type) cycle
         if(bc(faces(i)%itype)%itype == 0) cycle
         ff => faces(i)
         vecn = -ff%vecn
         cc => faces(i)%left_cell
         p1 = ff%centp
         qv = (nodes(ff%f2n(1))%qv + nodes(ff%f2n(2))%qv) * half
         rho1 = rhof(qv)
         u1 = zero
         p2 = p1 
         s=zero;s1=zero;s2=zero;s3=zero;su=zero
         do 
            call check_cross_face(p2,vecn,qv,ff,cc,nodes)
            dn = sum(vecn * (p2-p1))
            rho2 = rhof(qv)
            u2 = sqrt(sum(qv%v(imb:ime)**2))
            !
            u0 = half * ( u1 + u2)
            rho = half * (rho1 + rho2)
            s = s + dn
            g = u0 * dn
            su = su + g
            g = g * rho
            s1 = s1 + g
            g = g * u0
            s2 = s2 + g
            g = g * u0
            s3 = s3 + g
            if((u2-u1)/u2 < 0.005_rfp) exit
            if(ff%itype /= 0) exit   ! reach another side boundary
            rho1 = rho2; u1 = u2; p1 = p2 
         end do

         su = s - su / u2
         g = one / (rho2 * u2)
         s3 = g * (s1 - s3 / u2**2)
         s2 = g * (s1 - s2 / u2)
         s1 = s - s1 * g
         p = faces(i)%centp
         write(8,*)p,s,su,s1,s2,s3,p+vecn*s,p+vecn*s1 
      end do
      close(8)
   end subroutine thickness_bl
   !
   subroutine check_cross_face(p,vecn,qv,cfr,cc,nodes)
      type(cell),pointer::cc,cg,ch
      type(face),pointer::cf,cfr
      type(node),pointer::nodes(:)
      type(vector)::qv
      real(rfp)::vecn(ndim),p(ndim)
      integer::i,istop
      do i = 1, size(cc%sface)
         cf => cc%sface(i)%to_face
         if(associated(cfr,cf)) cycle
         cg => cf%right_cell
         if(associated(cc,cg)) cg => cf%left_cell
         if(insection(p,vecn,cf,nodes,qv)) then
            cc => cg
            cfr => cf
            return
         end if
      end do
   end subroutine check_cross_face
   !
   function insection(p,vecn,cf,nodes,qv)result(in)
      implicit none
      type(vector)::qv
      type(face),pointer::cf
      type(node),pointer::nodes(:)
      real(rfp)::p(ndim),p1(ndim),p2(ndim),f1,f2,vecn(:)
      logical::in
      !
      in =.false.
      p1 = nodes(cf%f2n(1))%xyz - p
      p2 = nodes(cf%f2n(2))%xyz - p
      f1 = p1(1) * vecn(2) - p1(2) * vecn(1)
      f2 = p2(1) * vecn(2) - p2(2) * vecn(1)
      if(f1 * f2 > zero) return
      in = .true.
      f1 = abs(f2) / (abs(f1) + abs(f2))
      f2 = one - f1
      qv = nodes(cf%f2n(1))%qv * f1 + nodes(cf%f2n(2))%qv * f2
      p = p1 * f1 + p2 * f2 + p
   end function insection

   subroutine wall_plot(fn,cells,nodes,faces,interf,it)
      !********************************************************
      implicit none
      !
      type(cell),pointer::cells(:)
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:),cf
      type(itf)::interf

      ! local variables
      type(vector)::qv,qvw
      real(rfp)::rc(ndim),rn(ndim),mach,qinf
      real(rfp),pointer::vars(:,:)
      integer::i,j,k,n,it,iw,nf,nn,f2n(4)
      integer,pointer::in(:),f2f(:)
      character(len=*)::fn
      ! 
      do iw = 1, nwplot 
         if(wlabel(iw) == 0) cycle    ! 
         select case(wlabel(iw))
         case(0)
            cycle  ! nothing for this boundary
         case(-1)      ! fluid/solid
            open(2,file=trim(fn)//'.sf.'//trim(i2s(it))//'.dat')
         case(-2)
            open(2,file=trim(fn)//'.tl.'//trim(i2s(it))//'.dat')
         case(-3)
            open(2,file=trim(fn)//'.mp.'//trim(i2s(it))//'.dat')
         case default 
            open(2,file=trim(fn)//'.'//trim(i2s(wlabel(iw)))//'.'//trim(i2s(it))//'.dat')
         end select
         !
         write(2,*)'variables ='
         if(ndim == 1) write(2,*)'x'
         if(ndim == 2) write(2,*)'x y'
         if(ndim == 3) write(2,*)'x y z'
         n = 0
         do i = 1, size(woutput_selection,1)
            if(.not.woutput_selection(i,iw)) cycle
            write(2,*)trim(wname(i))
            n = n + 1
         end do
         nf = 0
         do i = 1, size(faces)
            cf => faces(i)
            if(cf%itype == 0) exit
            if(wlabel(iw) == -1) then
               if(cf%itype /= partition_face_no + 1) cycle
            else
               if(cf%itype >= partition_face_no) cycle
               if(bc(cf%itype)%label /= wlabel(iw)) cycle
            end if
            nf = nf + 1   ! counting the number of faces
         end do
         allocate(vars(n,nf),f2f(nf),in(size(nodes)))
         nf = 0; in = 0
         do i = 1, size(faces)
            cf => faces(i)
            if(cf%itype == 0) exit
            if(wlabel(iw) == -1) then
               if(cf%itype /= partition_face_no + 1) cycle
            else
               if(cf%itype >= partition_face_no) cycle
               if(bc(cf%itype)%label /= wlabel(iw)) cycle
            end if
            call face_qv_tauw(cf,nodes,qv,qvw)
            !    if(bc(cf%itype)%igrp == inlet_type) qv = cf%right_cell%qv
            nf = nf + 1
            f2f(nf) = i
            in(faces(i)%f2n) = 1
            n = 0
            do j = 1, size(woutput_selection,1)
               if(.not.woutput_selection(j,iw)) cycle
               n = n + 1
               vars(n,nf) = face_value(cf,qv,qvw,j)
            end do
            if(ndim <= 2) write(2,20)cf%centp / b_lenref,vars(:,nf)
         end do
         ! output to tecplot
         if(ndim == 3) then
            nn = 0
            do i = 1, size(in)
               if(in(i) == 0) cycle
               nn = nn + 1
               in(i) = nn
            end do
            write(2,*)'zone T="w'//trim(i2s(wlabel(iw)))//'"'
            write(2,*)' n=',nn,',e=',nf, &
               ',varlocation=([4-'//trim(i2s(n+3))//']=cellcentered),'
            write(2,*)' f=feblock,et=quadrilateral'
            do i=1,size(in)
               if(in(i) /= 0) write(2,20)nodes(i)%xyz(1) / b_lenref
            end do
            do i=1,size(in)
               if(in(i) /= 0) write(2,20)nodes(i)%xyz(2) / b_lenref
            end do
            do i=1,size(in)
               if(in(i) /= 0) write(2,20)nodes(i)%xyz(3) / b_lenref
            end do
            do i = 1, n
               write(2,20)vars(i,:)
            end do
            do i=1,nf
               n = size(faces(f2f(i))%f2n)
               f2n(:n) = in(faces(f2f(i))%f2n)
               if(n < 4) f2n(4) = f2n(3)
               write(2,*)f2n
            end do
         end if   ! end dimension
         deallocate(vars,in,f2f)
         close(2)
      end do
      20  format(5e20.8)
   end subroutine wall_plot

   subroutine face_qv_tauw(cf,nodes,qv,tauw)
      type(face),pointer::cf
      type(node)::nodes(:)
      type(vector)::qv,tauw,gqv(ndim)
      real(rfp)::dc(neq),zk1,zk2,r1,r2,t1,t2,vj(ndim,2)
      integer::j
      qv%v = zero
      do j = 1, size(cf%f2n)
         qv%v = qv%v + nodes(cf%f2n(j))%qv%v
      end do
      qv%v = qv%v / real(size(cf%f2n),rfp)
      ! cal wall temperature
      !
      zk1 = lamda(cf%left_cell%qv,cf%left_cell)
      zk2 = lamda(cf%right_cell%qv,cf%right_cell)
      r1 = abs(dot_product(cf%vecn,cf%left_cell%centp - cf%centp))
      r2 = abs(dot_product(cf%vecn,cf%right_cell%centp - cf%centp))
      t1 = cf%left_cell%qv%v(ien)
      t2 = cf%right_cell%qv%v(ien)
      zk1 = zk1 / r1
      zk2 = zk2 / r2
      qv%v(ien) = (zk1 * t1 + zk2 * t2)/ (zk1 + zk2)
      !
      dc = diff_coef(qv,cf)
      dc(imb:ime) = mumix(qv)
      call face_gradient(cf,cf%right_cell,cf%left_cell,nodes(cf%f2n),gqv,vj)
      tauw = vis_t(qv,dc, gqv, cf%vecn)
   end subroutine face_qv_tauw

   function face_value(cf,qv,tauv,n)result(f)
      type(face),pointer::cf
      type(vector)::qv,tauv
      type(cell),pointer::cl,cr
      real(rfp)::dn,f,rhow
      integer::n
      cl => cf%left_cell
      cr => cf%right_cell
      dn = sum(cf%vecn * (cf%centp - cl%centp))
      select case(n)
      case(1)    ! pressure
         f = pf(qv)
      case(2)    ! temperature
         f = qv%v(ien)
      case(3)    ! density
         f = rhof(qv)
      case(4)    ! tw
         f = sqrt(sum(tauv%v(imb:ime)**2))
      case(5)    ! utw
         f = sqrt(sum(tauv%v(imb:ime)**2))
         f = sqrt(f / rhof(qv))
      case(6)    ! yplus
         f = sqrt(sum(tauv%v(imb:ime)**2))
         rhow = rhof(qv)
         f = sqrt(f / rhow)
         f = dn * f * rhow / mumix(qv)
      case(7)  ! heat flux
         !  f = -lamda(qv,cl) * (qv%v(ien) - cl%qv%v(ien)) / dn
         f = tauv%v(ien)
      case(8)  ! velocity x
         f = qv%v(imb)
      case(9)  ! velocity y
         f = qv%v(imb+1)
      case(10)  ! velocity z
         f = qv%v(ime)
      case(11)
         !  f = -one_third / absorption_coef(qv) * (cr%qv%v(irad) - cl%qv%v(irad)) / (dn * two)
         !  f = -one_third / absorption_coef(qv) * (qv%v(irad) - cl%qv%v(irad)) / dn
         f = - one_third / absorption_coef(qv) * sum(cf%vecn*cl%gradient(irad,:))
      case(12)
         f = qv%v(ieb)   ! Ex
      case(13)
         f = qv%v(ieb+1) ! Ey
      case(14)
         f = qv%v(iee)   ! Ez
      case(15)
         f = hf(qv)
      case(16)
         f = stagnation_pressure(qv,cf%centp)
      case(17)
         f = stagnation_temperature(qv,cf%centp)
      end select
      !
   end function face_value

   subroutine post_wall(faces,nodes,ip,interf)
      type(node)::nodes(:)
      type(face),pointer::faces(:),ff
      real(rfp)::un,u2,vecn(ndim),qinf,dn,cp,cf,rhow,ypone,taw,utaw,qw
      real(rfp)::uplus,vplus,kplus,mut,y0,rx,dc(neq)
      real(rfp)::t0,th,tc,zk0
      type(itf)::interf
      type(cell),pointer::cc
      type(vector)::qv,gqv(ndim),tawv
      real(rfp)::vj(ndim,2),xyz(3)
      integer::i,ij,ib,k1,k2,ip,j
      !  for nature convection problems
      t0 = 600._rfp
      th = 960._rfp
      tc = 240._rfp
      qv%v(ien) = t0
      !   zk0 = lamda(qv)
      !   zk0 = b_lenref/(zk0*(th-tc))

      if(s_ivis == 0) then
         open(8,file='gems.wall.'//i2s(ip))
         if(b_vinf < mytiny) then
            qinf = one
         else
            qinf = half * b_rhoinf * b_vinf**2
         end if
         do i = interf%nitf + 1, size(faces)
            if(faces(i)%itype == 0) exit
            if(bc(faces(i)%itype)%igrp /= wall_type) cycle
            ff => faces(i)
            vecn = -ff%vecn
            !    cc => faces(i)%right_cell
            cc => faces(i)%left_cell
            qv = cc%qv
            cp = (pf(qv) - g_pbase) / qinf
            xyz(:ndim) = ff%centp
            if(ndim == 2) xyz(3) = zero
            write(8,'(10e16.8)')xyz,-cp,pf(qv)
         end do
         close(8)
         return
      end if
      if(ip == 0) then
         open(8,file='gems.wall')
         write(8,*)'variables ='
         write(8,*)'x y z dy y^+ C_p C_f u_t T q rho'
         write(8,*)'zone t="wall",f="point"'
         close(8)
      end if
      open(8,file='gems.wall.'//i2s(ip))
      if(b_vinf < mytiny) then
         qinf = one
      else
         qinf = half * b_rhoinf * b_vinf**2
      end if

      !  do i = 1, size(faces)
      do i = interf%nitf + 1, size(faces)
         if(faces(i)%itype == 0) exit
         if(bc(faces(i)%itype)%igrp /= wall_type) cycle
         if(bc(faces(i)%itype)%itype == 0 ) cycle
         ff => faces(i)
         vecn = -ff%vecn
         cc => faces(i)%left_cell
         dn = dot_product(vecn, cc%centp - ff%centp)
         qv%v = zero
         do j = 1, size(ff%f2n)
            qv%v = qv%v + nodes(ff%f2n(j))%qv%v
         end do
         qv%v = qv%v / real(size(ff%f2n),rfp)
         cp = pf(qv)
         rhow = rhof(qv)
         dc = diff_coef(qv,faces(i))
         dc(imb:ime) = mumix(qv)

         call face_gradient(ff,ff%right_cell,ff%left_cell, &
            nodes(ff%f2n),gqv,vj)
         !
         tawv = vis_t(qv,dc, gqv, vecn)
         taw = sqrt(sum(tawv%v(imb:ime)**2))
         !
         utaw = sqrt(taw / rhow)
         ypone = dn * utaw * rhow / dc(imb)
         cf = taw / qinf
         xyz(:ndim) = ff%centp / b_lenref
         if(ndim == 2) xyz(3) = zero
         write(8,'(10e16.8)')xyz,dn,ypone,cp,cf,utaw,qv%v(ien),   &
            lamda(qv,cc) * tawv%v(ien),rhof(qv)
      end do
      close(8)
      !
   end subroutine post_wall

   subroutine user_define_output(cells,nodes,faces,interf,id)
      !*************************************************************

      type(cell),pointer::cells(:)
      type(node),pointer::nodes(:)
      type(face),pointer::faces(:)
      type(itf)::interf
      integer::i,id
      character(len=40)::fn
      if(s_iplot /= 5) return
      fn = 'sigma.dbs.'//i2s(id)
      open(10,file=fn)
      do i = 1, size(cells)
         write(10,*)e_resistivity(cells(i)%qv,cells(i))
      end do
      do i = 1, size(faces)
         if(faces(i)%itype == 0) exit
         write(10,*)e_resistivity(faces(i)%right_cell%qv,faces(i)%right_cell)
      end do
      close(10)
   end subroutine user_define_output
end module gems_output
