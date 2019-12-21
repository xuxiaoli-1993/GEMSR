MODULE dfd_react
  USE dfd_jacob
  use mpi
  implicit none
!
 contains

 subroutine chemical_reaction(cells)
!********************************************************
implicit none
!
  type(cell),pointer::cells(:)
! local variables
  type(cell),pointer::current_cell
  type(vector),pointer::qv
  real(rfp)::dm(nspe,neq)
  real(rfp)::rate(nspe),volume
  integer::i
  
 if(s_irealgas /= 3) return   ! no reaction

 do i = 1, size(cells)
!
  current_cell => cells(i)
  volume = current_cell%vol
  qv => current_cell%qv

  call reaction_rate(qv,rate,dm)
!
!  try 0eq of reaction 
!
!  call reaction_0eq(qv,rate,dm)
!
!
   current_cell%res%v(isb:ise) = current_cell%res%v(isb:ise) +  &
                                  rate(:nspm1) * volume
!
! jacobain matrix of source terms
!
   current_cell%dm%e(isb:ise,:) = current_cell%dm%e(isb:ise,:) +  &
                                  dm(:nspm1,:) * volume
 end do
 end subroutine chemical_reaction

 function write_reaction_rate(qv,i)result(w)
  type(vector)::qv
  real(rfp)::rate(nspe),w
  integer::i
  call reaction_rate(qv,rate)
  w = rate(i)
 end function write_reaction_rate

 subroutine reaction_rate(qv,rate,dm)
  type(vector),intent(in)::qv
  type(vector)::rhox
  real(rfp),intent(out),optional::dm(nspe,neq)
  real(rfp),intent(out)::rate(nspe)
  integer::i,j,k,l,m
  real(rfp)::xcon(nspe),yi(nspe),gk(nspe),dg,dh,hk(nspe)
  real(rfp)::prodf,prodb,sumf,sumb,suma,kf,kb,kc,q,qf,qb,af(3),ab(3)
  real(rfp)::rho,rhoinv,rhoy(nspm1),rhop,rhot,t,tinv,xcm
  real(rfp)::a,b,c,e,pd,td,yd,nu,epsi=1.e-25_rfp

   rate = zero
   if(present(dm)) dm  = zero
   
   rho = rhof(qv)
   t   = qv%v(ien)
   yi  = max(yif(qv),epsi)
   tinv = one / t
   xcon = rho * yi * species(:)%mwinv
!
!  Third body concentraction
   xcm = sum(xcon)   
   
   if(any(reaction(:)%kb_equil)) then
!    gk = (hif(qv) * tinv - sif(qv)) / species(:)%r  !gk = (h - Ts)/RT
    hk = hif(qv) * tinv / species(:)%r             !hk = h / RT
    gk = hk - sif(qv) / species(:)%r               !gk = (h - Ts)/RT
   end if
!
   do k=1, nrea
! forward reaction
     prodf = one
     kf    = zero
     do i = 1, reaction(k)%nf
       j=reaction(k)%locf(i)
       nu = reaction(k)%nurf(j)
       prodf = prodf *(xcon(j)**nu)
     end do
!
!  a*T^b*exp{E/T}
      af = reaction(k)%af
      kf = af(1) * (t**af(2)) * exp( af(3) * tinv )
!
! Backward reaction
     prodb = one
     kb = zero
     do i = 1, reaction(k)%nb
       j=reaction(k)%locb(i)
       nu = reaction(k)%nurb(j)
       prodb = prodb *(xcon(j)**nu)
     end do

    if(reaction(k)%kb_equil) then
! backward reaction rate calculated from equilibrium
      dg = 0; dh = 0    
      do i = 1, reaction(k)%n
        j = reaction(k)%loc(i)
        dg = dg + reaction(k)%nur(j) * gk(j)
        dh = dh + reaction(k)%nur(j) * hk(j)
      end do
      suma = sum(reaction(k)%nur)
      if(abs(suma) <= mytiny) then
       kc = exp(-dg)
      else
       kc = exp(-dg) * (g_patm / (t * g_runi)) ** suma
      end if
       kb = kf / kc
     else

!  a*T^b*exp{E/T}
       ab = reaction(k)%ab
       kb = ab(1) * (t**ab(2)) * exp( ab(3) * tinv )
     end if
!
!
       if(reaction(k)%third_body) then   ! exist third body
        kf = kf * xcm
	kb = kb * xcm
       end if

! the rate of progress
       qf  =  kf * prodf
       qb  =  kb * prodb
       q   =  qf - qb
     do i = 1, reaction(k)%n
       j = reaction(k)%loc(i)
       rate(j) = rate(j) + species(j)%mw * reaction(k)%nur(j) * q
     end do
!
!  Jacobian of reaction rate
!
    if(.not.present(dm)) cycle
     rhoinv = one / rho
     rhox = rhoxf(qv)
     rhop = rhox%v(ico)
     rhot = rhox%v(ien)
     rhoy = rhox%v(isb:ise)
!
     sumf = sum(reaction(k)%nurf)
     sumb = sum(reaction(k)%nurb)
     if(reaction(k)%third_body) then   ! exist third body
     sumf = sumf + one
     sumb = sumb + one
     end if
!
      pd = (qf * sumf - qb * sumb) * rhop * rhoinv
!
!                     alph                E
!
     a = (reaction(k)%af(2) - reaction(k)%af(3) * tinv ) * tinv 

   if(reaction(k)%kb_equil) then
     b = a / kf - (dh - suma) * tinv
   else
     b = (reaction(k)%ab(2) - reaction(k)%ab(3) * tinv ) * tinv
   end if

      a = a + sumf * rhot * rhoinv
      b = b + sumb * rhot * rhoinv
      td = qf * a - qb * b
     do i = 1, reaction(k)%n
       j = reaction(k)%loc(i)
       if(j >= nspe) cycle
       c = - species(j)%mw * reaction(k)%nur(j)
       dm(j,ico) = dm(j,ico) +  c * pd
       dm(j,ien) = dm(j,ien) +  c * td
!
      do l = 1, reaction(k)%n
        m = reaction(k)%loc(l)
!	if(yi(m) < epsi) cycle
        e = rhoy(m) * rhoinv
        a = sumf * e + reaction(k)%nurf(m) / yi(m) &
                     - reaction(k)%nurf(nspe) / yi(nspe)

        b = sumb * e + reaction(k)%nurb(m) / yi(m) &
                     - reaction(k)%nurb(nspe) / yi(nspe)
        if(reaction(k)%third_body) then   ! exist third body
         a = a + rho * species(m)%mwinv / xcm
         b = b + rho * species(m)%mwinv / xcm
        end if
!
        yd = qf * a - qb * b

        if(m >= nspe) cycle
        dm(j,m + isb - 1) = dm(j,m + isb - 1) + c * yd
      end do
     end do
   end do
!

 end subroutine reaction_rate
!
!
 subroutine reaction_0eq(qv,rate,dm)
  type(vector),intent(inout)::qv
  type(vector)::rhox
  real(rfp),intent(inout)::dm(nspe,neq)
  real(rfp),intent(inout)::rate(nspe)
  integer::i
  real(rfp)::dyi(nspe),yi(nspe)
  real(rfp)::rho,t,dtemp,dtime

   if(any(abs(rate) < 1.e5_rfp)) return
     dtime = 1.e-10_rfp
    
    do i = 1, 1000
     rho = rhof(qv)
     t   = qv%v(ien)
     yi  = yif(qv)
     dyi = rate * dtime / rho
     dtemp = - sum(dyi*hif(qv)) / sum(yi * htif(qv))
     qv%v(ien) = qv%v(ien) + dtemp
     qv%v(isb:ise) =    qv%v(isb:ise) + dyi(:nspm1)
     call reaction_rate(qv,rate)
   end do
     call reaction_rate(qv,rate,dm)
!
 end subroutine reaction_0eq




!********************************************************
!
 subroutine init_reaction(chemfile)
!
   character*50::terms(100),section(1000),reac(200)
   character*80::line
   integer::ne,ns,nr,nterms,nsect
   integer::i,j
   integer::nproc,id,ierr
!  temperal variables
!
   character*2, allocatable::elem(:)
   integer,allocatable::spec_elem(:,:)
   real(rfp),allocatable::mw_elem(:)
   character(len=*)::chemfile
   character(len=80)::thermfile,tranfile,ljfile
!
   namelist /prop_file/thermfile,tranfile,ljfile

  call mpi_comm_rank(mpi_comm_world,id,ierr)
  call mpi_comm_size(mpi_comm_world,nproc,ierr)

! if(id == 0) then !   read 

   open(2,file=chemfile)
   read(2,prop_file)
!
   nsect = 0
   do 
    read(2,'(80a)',end=10,err=10)line
! 
    call upstring(line(:4))
    
!    if(line(:4) == 'TRAN') then
!     call read_transport(2)
!     cycle
!    end if

    call convert_line_to_terms(line,terms,nterms)

    if(nterms == 0) cycle
  
    if(terms(nterms)(1:3) /= 'END') then
     do i = 1, nterms
      nsect = nsect + 1
      section(nsect) = terms(i)
     end do
    else

     do i = 1, nterms - 1
      nsect = nsect + 1
      section(nsect) = terms(i)
     end do

! Elements section
!
   if(section(1)(1:4) == 'ELEM') then
     ne = nsect - 1
     allocate(elem(ne),mw_elem(ne))
     elem(:) = section(2:nsect)
     call lookup_period_table(elem(:),mw_elem)
   if(id == 0) print *, 'Elements =',elem
   end if
!
! Species section
!
   if(section(1)(1:4) == 'SPEC') then
     ns = nsect - 1
! user specify molecular weight for this specie
!
     do i = 2, nsect
      if(section(i)(1:1) == '/') ns = ns - 1
     end do
!
     if(ns /= nspe) print *,'error in input ',ns,' /=',nspe
!
     allocate(spec_elem(ns,ne))
     species%mw = zero
!
     ns = 0
     do i = 2, nsect
      if(section(i)(1:1) /= '/') then
       ns = ns + 1
!     species(:)%name = section(2:nsect)
       species(ns)%name = section(i)
      else
       read(section(i)(2:),*)species(ns)%mw
      end if
     end do
!
     spec_elem = 0
     call find_elem_in_spec(species(:)%name,elem(:),spec_elem(:,:))
     do i = 1, ns
      if(species(i)%mw == zero) &
        species(i)%mw = sum(spec_elem(i,:) * mw_elem(:))
      species(i)%mwinv = one / species(i)%mw 
      species(i)%R = g_runi * species(i)%mwinv
     if(id == 0) print *, ' Species ',i,'=',species(i)%name, 'mw=',species(i)%mw
     end do
!
   end if
!
!    Nondimensional specific heat, enthalpy and extropy using 
!    NASA thermo new data format by Sanford Gordon & B.J. McBride,1996
!  
!    if(section(1)(1:4) == 'THER') then
!!     if(section(2)(:3) == 'ALL') call read_therm_old_nasa
!     if(section(2)(:3) == 'ALL') call read_NASA_therm(thermfile)
!    end if

!
! Reaction section
!
   if(section(1)(1:4) == 'REAC') then
      nr = 0
      nterms = 0
!
     read(section(2),*)nrea
     allocate(reaction(nrea))
!
     do j = 3, nsect
     if(index(section(j),'=') > 0) then
      nr = nr + 1
      section(nr) = section(j)
      reaction(nr)%kb_equil     = .false.
       do i = 1, 3
       read(section(i+j),*)reaction(nr)%af(i)
       end do
!       reaction(nr)%af(3) = - reaction(nr)%af(3) * 1.0e3_rfp / g_runi
!
       if((5+j > nsect).or.index(section(4+j),'=') > 0) then
       reaction(nr)%kb_equil = .true.
       cycle
       end if
 !      
       do i = 4, 6
       read(section(i+j),*)reaction(nr)%ab(i-3)
       end do
!       reaction(nr)%ab(3) = - reaction(nr)%ab(3) * 1.0e3_rfp / g_runi

!      if(nterms > 0) then
!       if(nterms /= 6) print *,'err in reaction',nterms,section(nr)
!       do i = 1, 3
!       read(terms(i),*)reaction(nr-1)%af(i)
!       end do

!       do i = 4, 6
!       read(terms(i),*)reaction(nr-1)%ab(i-3)
!       end do
!
!       nterms = 0
!      end if
!     else
!      nterms = nterms + 1
!      terms(nterms) = section(j)
!     end if
      end if
     end do
!
      if(nr /= nrea) print *,'error in reaction input',nr,nrea      
      call find_spec_in_reac(section(:nr),spec_elem(:,:))
!
!output the data of reactions
   if(id == 0) then
     do j = 1, nr
     write(*,'("Reaction ",i2,": ",a40)')j,section(j)
     write(*,'(6e14.5)')reaction(j)%af
     if(.not.reaction(j)%kb_equil) write(*,'(6e14.5)')reaction(j)%ab
     end do
   end if

   end if 

     nsect = 0
   end if
 end do
!
10  close(2)
    deallocate(mw_elem,spec_elem,elem)
!
    call read_NASA_therm(thermfile)
    call read_transport_cea(tranfile)
    call read_lennard_jones_potential(ljfile)

!  end if    ! end reading
 
!  if(nproc > 1) call broadcast_species(id)
!
 end subroutine init_reaction
 
 subroutine read_NASA_therm(thermfile)   ! Gordon & McBride  (1994)
 character(len=*)::thermfile
 character(len=80)::line
 integer::i,j,el,sl,n
 real(rfp)::a(5)
 real(rfp)::ch(7),cs(7)
 data ch/-1.0_rfp,zero,one, half,one_third,quarter,0.2_rfp/  !(-1,0,1,1/2,1/3,1/4,1/5)
 data cs/-0.5_rfp,-1.0_rfp,zero,one,half,one_third,quarter/  !(-1/2,-1,0,1,1/2,1/3,1/4)
! 
!  open(3,file='thermo.dat')
  open(3,file=thermfile)
  out: do i = 1, nspe
   el = len_trim(species(i)%name)
   rewind(3)
  do 
    read(3,'(80a)',end=10,err=10)line
    if(line(1:1) == '!') cycle  ! skip comment lines
    if(line(1:1) == ' ') cycle  ! skip comment lines
    sl = index(line(:15),' ') - 1
    if(sl /= el) cycle
    if(line(1:sl) == species(i)%name(1:el)) then
! first line
    read(3,1)n  ! number of T interval
!
    do j = 1, n
    read(3,2)species(i)%Titv(j),species(i)%Titv(j+1)
! second line
    read(3,3)a
    species(i)%cpp(:5,j) = a
! third line
    read(3,3)a
    species(i)%cpp(6:7,j) = a(:2)
! Fourth line

    species(i)%hhp(:7,j) = species(i)%cpp(:7,j) * ch
    species(i)%hhp(2,j)  = a(4)
    species(i)%hhp(8,j)  = species(i)%cpp(2,j)
!
!  H0/R=-a1*T^-1+b1++a3*T+a4*T^2/2+a5*T^3/3+a6*T^4/4+a7*T^5/5+a2*ln(T)
!  S0/R=-a1*T^-2/2-a2*T^-1+b2+a4*T^1+a5*T^2/2+a6*T^3/3+a7*T^4/4+a3*ln(T)

    species(i)%ssp(:7,j) = species(i)%cpp(:7,j) * cs
    species(i)%ssp(3,j)  = a(5)
    species(i)%ssp(8,j)  = species(i)%cpp(3,j)

    end do
    species(i)%interv = n
    cycle out
   end if
  end do
10 print *,' error! Can not find species ',species(i)%name,' in',thermfile
!  rewind(3)
 end do out
 close(3)

1 format(I2)
2 format(2F10.3)
3 FORMAT(5D16.8)

 end subroutine read_NASA_therm  

 subroutine read_lennard_jones_potential(ljfile)   ! Gordon & McBride  (1994)
 character(len=*)::ljfile
 character::name1*15,line*80
 integer::i,j,el,sl
 real(rfp)::cl,ce
! 
  open(3,file=ljfile)
  out: do i = 1, nspe
   el = len_trim(species(i)%name)
   rewind(3)
  do 
    read(3,'(80a)',end=10,err=10)line
    if(line(1:1) == '!') cycle  ! skip comment lines
    if(line(1:1) == ' ') cycle  ! skip comment lines
!
    read(line,1)name1,cl,ce
!
    sl = len_trim(name1)
    if(sl /= el) cycle
    if(name1(1:sl) == species(i)%name(1:el)) then
    species(i)%cl = cl
    species(i)%ce = ce
    cycle out
   end if
  end do
10 print *,' error! Can not find species ',species(i)%name,' in',ljfile
!  rewind(3)
 end do out
 close(3)

1 format(A15,2E19.10)
 end subroutine read_lennard_jones_potential

 subroutine read_transport_cea(tranfile)   ! Gordon & McBride  (1994)
 character(len=*)::tranfile
 character::name1*15,name2*15,v,c,comment*40,line*80
 integer::i,j,el,sl,n,iv,ic
! 
!  Gordon et al (1984)
!
!  lnEta = AlnT+B/T+C/T^2+D
!  lnLam = AlnT+B/T+C/T^2+D
!
  open(3,file=tranfile)
  out: do i = 1, nspe
   el = len_trim(species(i)%name)
   rewind(3)
  do 
    read(3,'(80a)',end=10,err=10)line
    if(line(1:1) == '!') cycle  ! skip comment lines
    if(line(1:1) == ' ') cycle  ! skip comment lines
    backspace (3)
    read(3,1)name1,name2,v,iv,c,ic,comment
    sl = len_trim(name2)
    if(sl /= 0) cycle   ! Interaction parameter skip
!    
    sl = len_trim(name1)
    if(sl /= el) cycle
    if(name1(1:sl) == species(i)%name(1:el)) then
! viscosity
    do j = 1, iv
    read(3,2)v,species(i)%vt(j),species(i)%vt(j+1),species(i)%mu(:,j)
    end do

! Conductivity
    do j = 1, ic
    read(3,2)c,species(i)%ct(j),species(i)%ct(j+1),species(i)%lamda(:,j)
    end do
!
    species(i)%iv = iv
    species(i)%ic = ic
    cycle out
   end if
  end do
10 print *,' error! Can not find species ',species(i)%name,' in',tranfile
!  rewind(3)
 end do out
 close(3)

1 format(A15,1x,A15,3x,A1,I1,A1,I1,2x,A40)
2 format(1x,A1,2F9.2,4E15.8)

 end subroutine read_transport_cea


! subroutine read_therm_old_nasa   ! Gordon & McBride  (1971)
! character*80::line
! integer::i,j,el,sl
! real::a(5)
! real(rfp)::c(5)
! data c/one, half,one_third,quarter,0.2_rfp/  !(1,1/2,1/3,1/4,1/5)
!! 
!  open(3,file='therm_old_nasa.dat')
!  out: do i = 1, nspe
!   el = len_trim(species(i)%name)
!   rewind(3)
!  do 
!    read(3,'(80a)',end=10,err=10)line
!    if(line(1:1) == '!') cycle  ! skip comment lines
!    if(line(1:1) == ' ') cycle  ! skip comment lines
!    sl = index(line(:24),' ') - 1
!    if(sl /= el) cycle
!    if(line(1:sl) == species(i)%name(1:el)) then
!    read(line,1)species(i)%tmin,species(i)%tmax,species(i)%tcom
!! first line
!    read(3,2)a
!    species(i)%cp_l(:5) = a
!! second line
!    read(3,2)a
!    species(i)%hh_l(1) = a(1)
!    species(i)%cp_h(:3) = a(3:5)
!! third line
!    read(3,2)a
!    species(i)%cp_h(4:5) = a(:2)
!    species(i)%hh_h(1)   = a(3)
!!
!    species(i)%hh_l(2:6) = species(i)%cp_l * c
!    species(i)%hh_h(2:6) = species(i)%cp_h * c
!!
!    cycle out
!   end if
!  end do
!10 print *,' error! Can not find species ',species(i)%name
!!  rewind(3)
! end do out
! close(3)

!1 format(45x,3F10.3)    
!2 FORMAT(5E15.8)

! end subroutine read_therm_old_nasa  

! subroutine read_transport(lu)
! integer::lu,el,sl
! character*80::line
! integer::i,j
! 
!out: do i = 1, nspe
!  read(lu,'(80a)') line
!  line = adjustl(line)
!  el = index(line,' ') - 1
!  call upstring(line(:el))
!  if(line(:3) == 'END') return
!!  
!    do j = 1, nspe
!      sl = len_trim(species(j)%name)
!      if(sl /= el) cycle
!      if(line(1:sl) == species(j)%name(1:el)) then
!         read(lu,*)species(j)%mu(:)
!         read(lu,*)species(j)%lamda(:)
!         read(lu,*)species(j)%cl,species(j)%ce
!	 species(j)%ce = boltzmann_c / species(j)%ce  ! Kb/Epsi
!         cycle out
!      end if
!    end do
!10 print *,' error! Can not find species ',line(:el)  
!!
! end do out
!!
!  read(lu,'(80a)') line
!  call upstring(line(:3))
!  if(line(:3) /= 'END') print *,'error in reading transport prop.'
!
! end subroutine read_transport

 subroutine lookup_period_table(elem,mw)
 character*2::elem(:)
 character*80::line
 real(rfp)::mw(:)
 integer::ne,i,j,el,natm = 100
 character*2::symbol(100)
 real(rfp)::atmwt(100)
 data symbol/'H ','D ','HE','LI','BE','B ','C ','N ','O ','F ', &
             'NE','NA','MG','AL','SI','P ','S ','CL','AR','K ','CA','SC',&
	     'TI','V ','CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS',&
             'SE','BR','KR','RB','SR','Y ','ZR','NB','MO','TC','RU','RH',&
             'PD','AG','CD','IN','SN','SB','TE','I ','XE','CS','BA','LA',&
             'CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM',&
             'YB','LU','HF','TA','W ','RE','OS','IR','PT','AU','HG','TL',&
             'PB','BI','PO','AT','RN','FR','RA','AC','TH','PA','U ','NP',&
             'PU','AM','CM','BK','CF','ES'/
!
!  ATOMIC WEIGHTS - Atomic Weights of the Elements 1995. Pure & Appl.
!     Chem. Vol.68, No.12, 1996, pp.2339-2359.
 data atmwt/               1.00794D0,2.014102D0,4.002602D0,6.941D0, &
       9.012182D0,10.811D0,12.0107D0,14.00674D0,15.9994D0,18.9984032D0,&
       20.1797D0,22.989770D0,24.305D0,26.981538D0,28.0855D0,30.973761D0,&
       32.066D0,35.4527D0,39.948D0,39.0983D0,40.078D0,44.95591D0,&
       47.867D0, 50.9415D0,51.9961D0,54.938049D0,&
       55.845D0,58.933200D0,58.6934D0,63.546D0,65.39D0,69.723D0,72.61D0,&
       74.92160D0,78.96D0,79.904D0,83.80D0,85.4678D0,87.62D0,88.90585D0,&
       91.224D0,92.90638D0,95.94D0,97.9072D0,101.07D0,102.9055D0,&
       106.42D0,&
       107.8682D0,112.411D0,114.818D0,118.710D0, 121.760D0,127.6D0,&
       126.90447D0,131.29D0,132.90545D0,137.327D0,138.9055D0,140.116D0,&
       140.90765D0,144.9127D0,145.D0,150.36D0,151.964D0,157.25D0,&
       158.92534D0,&
       162.50D0,164.93032D0,167.26D0,168.93421D0,173.04D0,174.967D0,&
       178.49D0,180.9479D0,183.84D0,186.207D0,190.23D0,192.217D0,&
       195.078D0,196.96655D0,200.59D0,204.3833D0,207.2D0,208.98038D0,&
       208.9824D0, 209.9871D0,&
       222.0176D0,223.0197D0,226.0254D0,227.0278D0,232.0381D0,&
       231.03588D0,238.0289D0,237.0482D0,244.0642D0,243.0614D0,&
       247.0703D0,247.0703D0,251.0587D0,252.083D0/

!   open(3,file='atomic.dat')
 ne = size(elem)
!out: do i = 1, ne
!  rewind(3)
!    el = len_trim(elem(i))
!    do 
!     read(3,'(80a)',end=10,err=10)line
!   if(line(1:1) == '!') cycle  ! skip comment lines
!   if(line(1:1) == ' ') cycle  ! skip comment lines
!   if(line(1:el) == elem(i)(1:el)) then
!    read(line(3:),*)mw(i)
!    cycle out
!   end if
!  end do
!10 print *,' error! Can not find element ',elem(i)  
! end do out
! close(3)
    
    ne = size(elem)
  out:  do i = 1, ne
     el = len_trim(elem(i))
     do j = 1,natm
      if(symbol(j)(1:el) == elem(i)(1:el)) then
       mw(i) = atmwt(j)
       cycle out
      end if
     end do
     print *,' error! Can not find element ',elem(i)  
   end do out

 end subroutine lookup_period_table

 subroutine find_elem_in_spec(spec,elem,se)
 character(*)::elem(:)
 character(*)::spec(:)
 character*2,wrk
 integer::se(:,:)
 integer::ne,ns,i,j,k,m,nw,ia,ib,l
  ne = size(elem)
  ns = size(spec)

  do k = 1, ns

!  lens = len_trim(spec(k))
 
  do i = 1, ne
  wrk = elem(i)
  nw = len_trim(wrk)
  j = index(spec(k),wrk(:nw))
  if(j <= 0) cycle
  j = j + 1
  ia = 0
  do l = j, len_trim(spec(k))+1
  ib = iachar(spec(k)(l:l)) - 48
  if(ib > 0 .and. ib <= 9) then
   ia = ia * 10 + ib
  else
   if(ia == 0) ia = 1
   exit
  end if
  end do
   se(k,i) = se(k,i) + ia
  end do
 end do
 
 end subroutine find_elem_in_spec

 subroutine find_spec_in_reac(reac,se)
 character(*)::reac(:)
 character*24::terms(50),gterm
 integer::i,j,k,len,nt1,nt2,l,ia,ll,ib
 integer::se(:,:),ich(size(se,2))
 logical::global = .false.
! 
out:  do k = 1, nrea
!
      reaction(k)%third_body = .false.
      reaction(k)%n = 0
      reaction(k)%nf = 0
      reaction(k)%nb = 0
      reaction(k)%nurf = 0
      reaction(k)%nurb = 0
!
  call convert_reac_to_specs(reac(k),terms,nt1,nt2)
      ich = 0   ! let's check the reaction equation
!      
middle:  do j = 1, nt2
   l = len_trim(terms(j))
   ia = 0
   do ll = 1, 100
   ib = iachar(terms(j)(1:1)) - 48
   if(ib > 0 .and. ib <= 9) then
    terms(j)(1:l-1) = terms(j)(2:l)
    l = l - 1
    ia =ia*10+ ib
   else
    if(ia == 0) ia = 1
    exit
   end if
   end do
   if(terms(j)(:l) == 'M') then
     reaction(k)%third_body = .true.
     cycle middle
   end if
!Global reaction model
!
   if(terms(j)(:1) == 'G') then
     gterm = terms(j)(3:)
     global = .true.
     cycle middle
   end if

in:   do i = 1, nspe
       len = len_trim(species(i)%name)
       if(terms(j)(:l) == species(i)%name(:len)) then
        if(j <= nt1) then  ! forward reaction
         ich = ich + ia * se(i,:)
         reaction(k)%nurf(i) = reaction(k)%nurf(i) + ia
	else
         ich = ich - ia * se(i,:)
         reaction(k)%nurb(i) = reaction(k)%nurb(i) + ia
	end if
        cycle middle
       end if
      end do in
!
     print *, 'can not find ',terms(j)(:l),' of ',reac(k)
  end do middle
  if(any(ich /= 0)) then
      print *,' error in equation-------------'
      print *,'   +++++', reac(k),' +++++'
!      stop
  end if
!
! stoichiometri coeficient nur = nurb - nurf
!

   reaction(k)%nur = reaction(k)%nurb - reaction(k)%nurf
!
!
   do i = 1, nspe
!
    if(reaction(k)%nurf(i) /= 0) then
     reaction(k)%nf = reaction(k)%nf + 1
     reaction(k)%locf(reaction(k)%nf) = i
    end if
!
    if(reaction(k)%nurb(i) /= 0) then
     reaction(k)%nb = reaction(k)%nb + 1
     reaction(k)%locb(reaction(k)%nb) = i
    end if
!    
    if(reaction(k)%nur(i) /= 0) then
     reaction(k)%n = reaction(k)%n + 1
     reaction(k)%loc(reaction(k)%n) = i
    end if
   end do
!   
!Global reaction model
    if(global) then
     read(gterm,*)(reaction(k)%nurf(reaction(k)%locf(l)),l=1,reaction(k)%nf)
     reaction(k)%nb = 0
     reaction(k)%kb_equil = .false.
    end if
!
  end do out

 end subroutine find_spec_in_reac

  subroutine convert_reac_to_specs(line,terms,neq,nterms)
   character(*)::line,terms(:)
   integer::nterms,ic,i,neq

   ic = 0
   nterms = 0
   do i = 1, len_trim(line)
    if(line(i:i) /= '+' .and. line(i:i) /= '=') then
      if(ic == 0) then   ! new terms
       nterms = nterms + 1
       terms(nterms) = ''
       ic = 1
      else
       ic = ic + 1
      end if
       terms(nterms)(ic:ic) = line(i:i)
     else
      if(line(i:i) == '=') neq = nterms
      if(ic /= 0) ic = 0
     end if
    end do
  end subroutine convert_reac_to_specs


  subroutine convert_line_to_terms(line,terms,nterms)
   character(*)::line,terms(:)
   integer::nterms,ic,i
   
   line = adjustl(line)
   ic = 0
   nterms = 0
   do i = 1, len_trim(line)
    if(line(i:i) == '!') exit   ! skip all comments after !
    if(line(i:i) /= ' ' .and. line(i:i) /= '	') then
      if(ic == 0) then   ! new terms
       nterms = nterms + 1
       terms(nterms) = ''
       ic = 1
      else
       ic = ic + 1
      end if
       terms(nterms)(ic:ic) = upcase(line(i:i))
     else
      if(ic /= 0) ic = 0
     end if
    end do
  end subroutine convert_line_to_terms
    

  function upcase(c)result(newc)
   character*1::c,newc
   integer::ic
   newc = c         ! ASCII  a=97 z=122 A=65 Z=90
   ic = iachar(c)
   if(ic >= 97 .and. ic <= 122) newc = achar(ic - 32)
  end function upcase

  subroutine upstring(c)
   character(*)::c
   integer::i
   do i = 1, len_trim(c)
   c(i:i) = upcase(c(i:i))
   end do
  end subroutine upstring

  function sort_string(string)result(is)
  character(*)::string(:)
  character*24::s,scopy(size(string))
  integer::is(size(string))
  integer::i,j,n
   n = size(string)
   is = (/(i,i=1,n)/)
   scopy = string
   do j = 2, n
    s = scopy(j)
    do i = j - 1, 1, -1  
    if(lle(scopy(i),s)) exit
    scopy(i+1) = scopy(i)
    is(i+1) = is(i)
    end do
    scopy(i+1) = s
    is(i+1) = j
   end do
   end function sort_string
    
 subroutine broadcast_species(id)
  character,allocatable::buffer(:)
  integer::n,m,i,ierr,id,ipos
! Thermo-property
   n = 24 + (6 + 4 + 31 + 24 + 24 + 40) * rfp
   call mpi_bcast(nrea,1,mpi_packed,0,mpi_comm_world,ierr)
   allocate(buffer(n))
   if(id /= 0) allocate(reaction(nrea))

  do i = 1, nspe
  if(id == 0) then
   ipos = 1
   call mpi_pack(species(i)%name,24,mpi_character,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%mw,rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%mwinv,rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%r,rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%cl,rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%ce,rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%interv,1,mpi_integer,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%titv,4*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%cpp,31*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%hhp,24*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%ssp,24*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%iv,1,mpi_integer,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%ic,1,mpi_integer,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%vt,4*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%ct,4*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%mu,12*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(species(i)%lamda,12*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
!
   call mpi_bcast(buffer,n,mpi_packed,0,mpi_comm_world,ierr)
!
  else
   call mpi_bcast(buffer,n,mpi_packed,0,mpi_comm_world,ierr)
   ipos = 1
   call mpi_unpack(buffer,n,ipos,species(i)%name,24,mpi_character,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%mw,rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%mwinv,rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%r,rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%cl,rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%ce,rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%interv,1,mpi_integer,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%titv,4*rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%cpp,31*rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%hhp,24*rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%ssp,24*rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%iv,1,mpi_integer,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%ic,1,mpi_integer,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%vt,4*rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%ct,4*rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%mu,12*rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,species(i)%lamda,12*rfp,mpi_byte,mpi_comm_world,ierr)
  end if

  end do
  deallocate(buffer)
!
!
! Reaction
!
  n = ( (1 + 2 * nspe ) * 3  + 3 * 2) * rfp + 2
  allocate(buffer(n))
!
  do i = 1, nrea
  if(id == 0) then
   ipos = 1
   call mpi_pack(reaction(i)%nf,1,mpi_integer,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(reaction(i)%locf,nspe,mpi_integer,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(reaction(i)%nurf,nspe*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(reaction(i)%af,3*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)

   call mpi_pack(reaction(i)%nb,1,mpi_integer,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(reaction(i)%locb,nspe,mpi_integer,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(reaction(i)%nurb,nspe*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(reaction(i)%ab,3*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)

   call mpi_pack(reaction(i)%n,1,mpi_integer,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(reaction(i)%loc,nspe,mpi_integer,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(reaction(i)%nur,nspe*rfp,mpi_byte,buffer,n,ipos,mpi_comm_world,ierr)
!
   call mpi_pack(reaction(i)%third_body,1,mpi_logical,buffer,n,ipos,mpi_comm_world,ierr)
   call mpi_pack(reaction(i)%kb_equil,1,mpi_logical,buffer,n,ipos,mpi_comm_world,ierr)
!
   call mpi_bcast(buffer,n,mpi_packed,0,mpi_comm_world,ierr)   
 else
   call mpi_bcast(buffer,n,mpi_packed,0,mpi_comm_world,ierr)   
   ipos = 1
!
   call mpi_unpack(buffer,n,ipos,reaction(i)%nf,1,mpi_integer,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,reaction(i)%locf,nspe,mpi_integer,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,reaction(i)%nurf,nspe*rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,reaction(i)%af,3*rfp,mpi_byte,mpi_comm_world,ierr)
!
   call mpi_unpack(buffer,n,ipos,reaction(i)%nb,1,mpi_integer,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,reaction(i)%locb,nspe,mpi_integer,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,reaction(i)%nurb,nspe*rfp,mpi_byte,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,reaction(i)%ab,3*rfp,mpi_byte,mpi_comm_world,ierr)
!
   call mpi_unpack(buffer,n,ipos,reaction(i)%n,1,mpi_integer,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,reaction(i)%loc,nspe,mpi_integer,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,reaction(i)%nur,nspe*rfp,mpi_byte,mpi_comm_world,ierr)
!
   call mpi_unpack(buffer,n,ipos,reaction(i)%third_body,1,mpi_logical,mpi_comm_world,ierr)
   call mpi_unpack(buffer,n,ipos,reaction(i)%kb_equil,1,mpi_logical,mpi_comm_world,ierr)
   end if
  end do

  deallocate(buffer)
!
 end subroutine broadcast_species

end module dfd_react
