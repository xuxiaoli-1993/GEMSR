module gems_precon
  use gems_state

contains

  subroutine precondition_cell(switch, cells, nodes, nadv) !cc,rhopp)
    implicit none
    type(node), pointer :: nodes(:)
    type(cell), pointer :: cc, cl, cr, cells(:), nc
    type(face), pointer :: cf
    integer, intent(in) :: switch
    type(vector) :: qv, qvc, qvav, qvn
    !  real(rfp),intent(out)::rhopp
    real(rfp) :: u2, v2, dl, mu, um, c2, coef, rhop, rho, ug, u_small
    integer :: i, k, nadv, nmax = 100, db_dummy
    ! without preconditioning
    if (switch <= 0) return
    !coef = min(nadv, nmax) / real(nmax, rfp)
    !
    ! switch = 0 no precondition
    !          1 max vel and viscous vel
    !          2 max vel and local pressure
    !          3 max vel
    !          4 only cell vel
    !          5 max vel and perturbation pressure
    !          6 cutoff with vnn 

    u_small = 1.e-6_rfp

    do k = 1, size(cells)
      cc => cells(k)

      qvc = cc % qv
      um = max(dot_product(qvc % v(imb:ime), qvc % v(imb:ime)), u_small) !cf%vecn) 
      if (switch == 4) goto 10

      !average among neighbors
      do i = 1, size(cc % sface)
        cf => cc % sface(i) % to_face
        nc => cc%scell(i)%to_cell

        qvn = nc % qv
        qv % v = (qvc % v + qvn % v) * half               

        u2 = dot_product(qv % v(imb:ime), qv % v(imb:ime)) !cf%vecn) 
        um = max(um, u2)

        !if(g_igrav > 0) then  ! gravity preconditioning DUBIOUS
        !  u2 = sum(g_gravity * (nc%centp-cc%centp))**2 / max(u2,1.e-10_rfp)
        !  um = max(um,u2)
        !end if
        !
        if (switch == 5) then
          !perturbation pressure preconditioning
          rho = rhof(cc % qv)
          u2 = 4.0_rfp * abs(cc % qv % v(ico) - qvn % v(ico)) / rho
          !u2 = 4.0_rfp * abs(cl % qv % v(ico) - cr % qv % v(ico)) / rhof(qv)
          ! if(u2 > um)print *,um,u2
          um = max(um, u2)
        end if

        dl = cc % vol / max(cf % area, mytiny) ! height

        ! unstedy preconditioning
        if (s_idt > 0 .and. s_vpre_uns>-1000) then 
          v2 = dl * s_dt_inv
          !v2 = dl * s_dt_inv

          if(s_vpre_uns > zero) then
             !v2 = min(v2, s_vpre_uns)
             v2 = s_vpre_uns
          end if
          v2 = v2 * v2
          um = max(um, v2)
        end if

        !if (switch > 1) cycle

        ! viscous preconditioning
        if (s_ivis /= 0) then
          !     dl = abs(dot_product(cr%centp - cl%centp,cf%vecn)) !cr%centp - cl%centp))  !cf%vecn))
          if (cf % itype /= 0.and.cf % itype >= partition_face_no) cycle
          if (cf % itype /= 0.and.bc(cf % itype) % igrp == geom_topology) cycle
          !     if(ndim == 3) then
          !      dl = sqrt(cf%area) / dl / dl  ! AR=(side length / height) / height
          !     else if(s_iaxis == 0) then
          !      dl = cf%area / dl / dl
          !     end if
          rho = rhof(qv)
          mu = mumix(qv)

          !     if(s_ivis == 2) mu = mu + mu_turb(qv)
          mu = mu / rho
          !mu = mu / rhof(qv) ! * (s_cfl/s_vnn)
          v2 = mu / dl
          !v2 = mu * dl ! mu * AR / Dy
          v2 = v2 * v2
          um = max(um, v2)
        end if
      end do

      if (switch == 2) then
        ! local pressure preconditioning
        rho = rhof(qv)
        um = max(um, abs(qv % v(ico))/rho * quarter)
        !um = max(um, abs(qv % v(ico))/rhof(qv) * quarter)
      end if

      if (switch == 6) um = max(um, s_vnn) ! cut-off by vnn number
      !um = max(um, s_vnn)

      !    if(switch == 3) then
      10 continue
      !    print *,um,rhop
      !     if(nadv < nmax) &
      !     um = (um * rhop ) ** coef / rhop
      !    end if

      rhop = rhopf(qvc)
      cc % rp = max(one /um, rhop)
    end do
    !
  end subroutine precondition_cell

  function length_for_precondition(cc, nodes) result(dl)
    type(node) :: nodes(:)
    type(cell) :: cc
    real(rfp) :: dx, dy, dl, p(size(cc % c2n)), dd(ndim), vel(ndim)
    integer :: i
    vel = abs(cc % qv % v(imb:ime))
    do i = 1, ndim
      p = nodes(cc % c2n) % xyz(i)
      dd(i) = maxval(p) - minval(p)
    end do
    dx = dot_product(dd, vel) / sqrt(sum(vel * vel))
    vel(2:3) = vel(3:2:-1)

    dy = dot_product(dd, vel) / sqrt(sum(vel * vel))
    dl = dy * dy / dx
    ! dl = maxval(dd)
  end function length_for_precondition

  subroutine precondition_face(switch, qv, cf, cl, cr, rhopp)
    implicit none
    integer, intent(in) :: switch
    type(face), intent(in) :: cf
    type(cell), intent(in) :: cl, cr
    type(vector), intent(in) :: qv

    real(rfp) :: u2, v2, um, mu, rho, dl, c2, ar, alpha
    real(rfp), intent(out) :: rhopp
    ! without preconditioning
    if (switch <= 0) return

    ! 
    if (cf % itype == 0) then
       rhopp = (cl % rp + cr % rp) * half
    else
       rhopp = cl % rp
    end if
    return

    c2 = sound_speed(qv)
    v2 = zero
    !
    !   u2 = abs(dot_product(qv%v(imb:ime),cf%vecn))
    u2 = dot_product(qv % v(imb:ime), qv % v(imb:ime)) !cf%vecn) 
    if (s_ivis /= 0) then
      dl = abs(dot_product(cr % centp - cl % centp, cf % vecn))
      rho = rhof(qv)
      mu = mumix(qv)
      if (s_ivis == 2) mu = mu + mu_turb(qv)
      mu = s_cfl / s_vnn * mu / rho
      v2 = mu / dl
      v2 = v2 * v2
      !   v2 = v2*(v2-sqrt(u2))/(one+v2*sqrt(u2)/c2-u2/c2)
      !   v2 = alpha * ( alpha - one) * u2 / (one + (alpha - one)*u2/c2)
    end if
    um = max(u2, v2, mytiny)
    !   um = um * um
    if (um > 0.8_rfp * c2) um = c2 ! mach > 0.5 turn off precondition
    rhopp = one / um
  end subroutine precondition_face

  !subroutine precondition(switch,gas,cells,faces)
  !  implicit none
  !  type(property),intent(in)::gas
  !  type(face)::faces(:)
  !  integer,intent(in)::switch
  !  type(cell),pointer::cells(:),cl,cr
  !  real(rfp)::vel(ndim),vecn(ndim)
  !  real(rfp)::mul,mur,rhol,rhor,p,t
  !  real(rfp)::u2,c2,dl,re
  !  integer::i
  ! without preconditioning
  !   if(switch <= 0) then
  !    do i = 1, size(cells)
  !     cl => cells(i)
  !     call qv2prime(gas,cl%qv,p,vel,t)
  !     rhol  =  rhof(gas,p,t)
  !     cl%rhopp = rhopf(gas,rhol,p,t)  ! no precondition
  !    end do  ! end of loop cells
  !   return
  !   end if
  !
  ! loop all cells
  !   do i = 1, size(cells)
  !   cl => cells(i)
  !   call qv2prime(gas,cl%qv,p,vel,t)
  !   rhol  =  rhof(gas,p,t)
  !   u2 = dot_product(vel,vel)
  !   dl = cl%vol ** (one / real(ndim,rfp))
  !   mul = gas%zmu / rhol
  !   re = 4.0_rfp * mul / dl
  !   cl%rhopp = max(u2, re*re)  
  !   end do  ! end of loop cells
  !!   
  !! loop all faces
  !   do i = 1, size(faces)
  !    cl => faces(i)%left_cell
  !    cr => faces(i)%right_cell
  !    vecn = faces(i)%vecn
  !    call qv2prime(gas,cl%qv,p,vel,t)
  !    rhol  =  rhof(gas,p,t)
  !    call qv2prime(gas,cr%qv,p,vel,t)
  !    rhor  =  rhof(gas,p,t)
  !    mul = gas%zmu / rhol
  !    mur = gas%zmu / rhor
  !    dl = dot_product(cl%centp - cr%centp, vecn)
  !    u2 = 4.0_rfp * abs(cl%qv%v(ico) - cr%qv%v(ico)) / (rhol + rhor)
  !    re = 4.0_rfp * max(mul, mur) / dl
  !    u2 = max(u2, re*re)
  !!    cl%rhopp = max(cl%rhopp, u2)
  !    if(faces(i)%itype == 0) &
  !!    cr%rhopp = max(cr%rhopp, u2)
  !  end do  ! end loop of face
  !!  
  !! turn off precondition if Mref > 0.5
  !
  !  do i = 1, size(cells)
  !   cl => cells(i)
  !   u2 = cl%rhopp
  !   call qv2prime(gas,cl%qv,p,vel,t)
  !   rhol  =  rhof(gas,p,t)
  !   c2 = sound_speed(gas,t)  
  !   if ( u2 > half * c2) then
  !    cl%rhopp = rhopf(gas,rhol,p,t)  ! no precondition
  !   else
  !    u2 = max(u2,1.e-5_rfp * abs(cl%qv%v(ico)/rhol),1.e-20_rfp)
  !    cl%rhopp = one / u2
  !   end if
  !  end do
  !!
  !!
  !end subroutine precondition

end module gems_precon
