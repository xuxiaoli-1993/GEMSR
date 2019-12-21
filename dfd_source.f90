MODULE dfd_source
    USE dfd_state
    use mpi
    implicit none
contains
    !
    subroutine source(cells, nodes)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        type(node), pointer :: nodes(:)
        ! local variables
        type(cell), pointer :: current_cell
        type(vector), pointer :: qv
        real(rfp) :: volume
        type(matrix), pointer :: dm
        integer :: i
        if (ndim == 2 .and. s_iaxis > 0) then
            if (neqf > 0) call axis_source(cells, nodes)
            if (neqm > 0) call axis_mhd_source(cells, nodes)
        end if
        if (neqm > 0.and.s_isource /= 6) call electric_current(cells)
        if (g_igrav == 1) call gravity_source(cells)
        if (nrad > 0) call radiation_source(cells)
        select case(s_isource)
        case(0)
            return
        case(1)
            call external_source(cells)
        case(2)
            call growth_source(cells)
        case(3)
            call propane_global_reaction_rate(cells)
        case(4)
            call liquid_vapor_source(cells)
        case(5)
            call cavitation_source(cells)
        case(6)
            call given_electric_current(cells)
        case(7)
            call blackbody_volumetric_radiation(cells)
        case(8)
            call Coriolis_force(cells)
        case(9) ! spurce term for solid melting modeling
            call solid_melting(cells)
        case default
            return
        end select
        !
    end subroutine source

    subroutine electric_current(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: cc
        type(vector), pointer :: res, qv
        type(matrix), pointer :: dm
        real(rfp) :: volume, e(3), ev(3), v(3), b(3), g(3, 3), w(3), eta, beta, jb(3)
        real(rfp) :: je, lorenz_force(3), vol, mat(3, 3), muvol, wmat(3, 3), sigma_p, sigma_t
        real(rfp) :: epsi = 0.01_rfp, rcoef
        integer :: i, it, nd1
        !
        v = zero ! no fluid flow
        !
        if (neqm == 0) return

        do i = 1, size(cells)
            cc => cells(i)
            res => cc % res
            dm => cc % dm
            qv => cc % qv
            it = cc % itype

            !if (it == 1) cycle

            !   volume = g_mu0 * cc%vol
            volume = cc % vol
            e = qv % v(ieb:iee)
            b = qv % v(ibb:ibe)
            select case(vc(it) % itype)
            case(0) ! Plasma zone
                b = bfun(qv, i)
                if (nmeq > 0) then
                    v = avto3(qv % v(imb:ime))
                else
                    v = zero
                end if
                if (s_iaux == 1) v(3) = qv % v(iaub)
                eta = e_resistivity(qv, cc)
                beta = hall_resistivity(qv, cc)
                g = cross2mat(beta * b)
                mat = cross2mat(v)
                ! turn off hall term
                g = zero
                w = eta
                call put_diag(g, eta, 1, 3)
                call direct_inv(g)
                ev = e + acrossb(v, b)
                jb = matmul(g, ev)
                !     jb = ev / eta
                cc % jv = jb
                !
                muvol = g_mu0 * volume
                res % v(ieb:iee) = res % v(ieb:iee) - jb * muvol !matmul(g,e*volume)
                !
                sigma_p = sigma_pf(qv, cc)
                sigma_t = sigma_tf(qv, cc)
                !
                dm % e(ieb:iee, ico) = dm % e(ieb:iee, ico) + sigma_p * eta * jb * muvol
                dm % e(ieb:iee, ien) = dm % e(ieb:iee, ien) + sigma_t * eta * jb * muvol
                !
                dm % e(ieb:iee, ibb:ibe) = dm % e(ieb:iee, ibb:ibe) + matmul(g, mat) * muvol
                dm % e(ieb:iee, ieb:iee) = dm % e(ieb:iee, ieb:iee) + g * muvol
                !
                if (neqf > 0 .and. neqm > 0) then
                    !
                    wmat = matmul(g, cross2mat(b)) * muvol
                    dm % e(ieb:iee, imb:ime) = dm % e(ieb:iee, imb:ime) - wmat(:,:ndim)
                    if (s_iaux == 1) &
                    dm % e(ieb:iee, iaub) = dm % e(ieb:iee, iaub) - wmat(:, 3)
                    !
                    lorenz_force = acrossb(jb, b)
                    je = sum(jb * qv % v(ieb:iee))
                    ! put limitation on it
                    !    rcoef =  max(zero, one - epsi * (qv%v(ien) / 1.e4_rfp)**4)
                    rcoef = one
                    je = je * rcoef
                    !  special treat for Vaccum
                    !    if(qv%v(isb) < 0.9_rfp) je = zero
                    res % v(imb:ime) = res % v(imb:ime) + lorenz_force(:ndim) * volume
                    if (s_iaux == 1) &
                    res % v(iaub) = res % v(iaub) + lorenz_force(3) * volume
                    res % v(ien) = res % v(ien) + je * volume
                    !  source jacobian
                    !
                    ! mometume
                    eta = eta * volume
                    !
                    dm % e(imb:ime, ico) = dm % e(imb:ime, ico) - sigma_p * eta * lorenz_force(:ndim)
                    dm % e(imb:ime, ien) = dm % e(imb:ime, ien) - sigma_t * eta * lorenz_force(:ndim)
                    if (s_iaux == 1) then
                        dm % e(iaub, ico) = dm % e(iaub, ico) - sigma_p * eta * lorenz_force(3)
                        dm % e(iaub, ien) = dm % e(iaub, ien) - sigma_t * eta * lorenz_force(3)
                    end if
                    !
                    !   energy
                    !    dm%e(ien,ico) = dm%e(ien,ico) - sigma_p * eta * je
                    !    dm%e(ien,ien) = dm%e(ien,ien) - sigma_t * eta * je
                    ! d(JxB)/d(V)
                    wmat = matmul(g, matcrossvec(cross2mat(b), b)) * volume
                    dm % e(imb:ime, imb:ime) = dm % e(imb:ime, imb:ime) + wmat(:ndim,:ndim)
                    if (s_iaux == 1) then
                        dm % e(iaub, imb:ime) = dm % e(iaub, imb:ime) + wmat(3,:ndim)
                        dm % e(imb:ime, iaub) = dm % e(imb:ime, iaub) + wmat(:ndim, 3)
                        dm % e(iaub, iaub) = dm % e(iaub, iaub) + wmat(3, 3)
                    end if
                    ! d(JxB)/d(B)
                    wmat = (matmul(g, matcrossvec(mat, b)) + cross2mat(jb)) * volume
                    dm % e(imb:ime, ibb:ibe) = dm % e(imb:ime, ibb:ibe) - wmat(:ndim,:)
                    if (s_iaux == 1) &
                    dm % e(iaub, ibb:ibe) = dm % e(iaub, ibb:ibe) - wmat(3,:)
                    ! d(JxB)/d(E)
                    wmat = matmul(g, cross2mat(B)) * volume
                    dm % e(imb:ime, ieb:iee) = dm % e(imb:ime, ieb:iee) + wmat(:ndim,:)
                    if (s_iaux == 1) &
                    dm % e(iaub, ieb:iee) = dm % e(iaub, ieb:iee) + wmat(3,:)
                    !
                    ! d(J.e)/d(B),d(E),d(V)
                    !    if(qv%v(isb) >= 0.9_rfp) then  ! special treat for vacuum
                    dm % e(ien, ibb:ibe) = dm % e(ien, ibb:ibe) - matmul(g, matmul(e, mat)) * volume * rcoef
                    dm % e(ien, ieb:iee) = dm % e(ien, ieb:iee) - (matmul(g, e) + jb) * volume * rcoef
                    !
                    w = matmul(g, matmul(e, cross2mat(b))) * volume * rcoef
                    dm % e(ien, imb:ime) = dm % e(ien, imb:ime) + w(:ndim)
                    if (s_iaux == 1) &
                    dm % e(ien, iaub) = dm % e(ien, iaub) + w(3)
                end if
                !   
                !    if(rcoef > zero) then
                !    dm%e(ien,ien) = dm%e(ien,ien) + & 
                !     je / rcoef * 4._rfp * epsi * (qv%v(ien)/1.e4_rfp)**3/1.e4_rfp
                !    end if

                !   end if
            case(1:) ! Conductor or Insulator zone
                eta = e_resistivity(qv, cc)
                w = volume * g_mu0 / eta
                cc % jv = e / eta
                if (eta > 1e5_rfp) then
                    w = zero
                    cc % jv = zero
                end if
                res % v(ieb:iee) = res % v(ieb:iee) - e * w !volume * g_mu0 / eta
                call diagadd(dm % e(ieb:iee, ieb:iee), w)
            end select
        end do
    end subroutine electric_current
    !
    subroutine given_electric_current(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: cc
        type(vector), pointer :: res, qv
        real(rfp) :: current, w, volume, tpeak, tmax, t
        integer :: i, it
        !
        if (neqm == 0) return

        do i = 1, size(cells)
            cc => cells(i)
            res => cc % res
            qv => cc % qv
            volume = cc % vol
            it = cc % itype
            select case(vc(it) % itype)
            case(0)
                current = vc(it) % vars(ndim + 1)
            case(1) ! time function
                tpeak = vc(it) % vars(1)
                tmax = vc(it) % vars(2)
                current = vc(it) % vars(ndim + 1)
                t = s_elapsed_time + one / s_dt_inv
                if (t < tpeak) then
                    current = t /tpeak * current
                else
                    current = (one - t / tmax) / (one - tpeak / tmax) * current
                end if
            case(2)
                !     current =  vc(it)%vars(ndim+1) * sin(two * pi * 1.e4_rfp * s_elapsed_time)
                current = 100. * sin(two * pi * 1.e4_rfp * s_elapsed_time)
                !     current =  vc(it)%vars(ndim+1) * sin(s_elapsed_time)
            end select
            w = volume * g_mu0 * current
            !     w = volume * current
            res % v(iee) = res % v(iee) - w
            cc % jv = zero
            cc % jv(3) = current
        end do
    end subroutine given_electric_current
    !
    subroutine solid_melting(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: cc
        type(face), pointer :: cf
        type(vector), pointer :: res, qv
        type(vector) :: qvf
        real(rfp) :: gl, kc, volume, darcy, w(ndim), un, h
        integer :: i, it, j
        !
        do i = 1, size(cells)
            cc => cells(i)
            res => cc % res
            qv => cc % qv
            volume = cc % vol
            ! Darcian damping force
            ! Kozeny-Carman equation
            gl = liquid_fraction(qv) + 1.e-20_rfp ! avoid singularity
            kc = metal % k0 * gl**3 / (one - gl)**2
            !
            darcy = -metal % mul * rhof(qv) / (kc * metal % rhol) * volume
            res % v(imb:ime) = res % v(imb:ime) + darcy * qv % v(imb:ime)
            w = -darcy
            call diagadd(cc % dm % e(imb:ime, imb:ime), w)
            ! Buoyancy force
            w = -metal % rhol * metal % gravity * metal % beta * volume
            res % v(imb:ime) = res % v(imb:ime) + w * (qv % v(ien) - metal % tl)
            !
            cc % dm % e(imb:ime, ien) = cc % dm % e(imb:ime, ien) - w
            ! Gravity force
            !  w =   rhof(qv) * metal%gravity  * volume
            !  res%v(imb:ime) = res%v(imb:ime) + w   
            !  correct enthalpy
            ! 
            do j = 1, size(cc % sface)
                cf => cc % sface(j) % to_face
                qvf = (cf % left_cell % qv + cf % right_cell % qv) * half
                un = rhof(qvf) * sum(cf % vecn * qvf % v(imb:ime)) * cf % area
                if (associated(cc, cf % right_cell)) un = -un
                gl = liquid_fraction(qvf)
!                 if (s_ilevelset > 0) then ! level set
!                     h = heavisidef(qv)
!                     res % v(ien) = res % v(ien) + un * (1 - h) * (one - gl) * ((metal % cps - metal % cpl) * (qvf % v(ien) - metal % ts) - metal % lheat)
!                     !    res%v(ien) = res%v(ien) + un * (1-h) * (one - gl) * metal%cps * qvf%v(ien)
!                     ! print *, qvf%v(ien)
!                 else
!                     res % v(ien) = res % v(ien) + un * (one - gl) * ((metal % cps - metal % cpl) * (qvf % v(ien) - metal % ts) - metal % lheat)
!                     !    res%v(ien) = res%v(ien) + un * (one - gl) * metal%cps * qvf%v(ien)
!                 end if
            end do
        end do
        ! 
    end subroutine solid_melting
    !

    subroutine axis_mhd_source(cells, nodes)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:), cc
        type(node), pointer :: nodes(:)
        ! local variables
        type(vector), pointer :: qv
        real(rfp) :: volume, rho, zmu, r, d, div, dvdr
        integer :: i, jcol
        if (s_iaxis == 0) return
        if (neqm == 0) return
        jcol = s_iaxis - 1
        !
        do i = 1, size(cells)
            qv => cells(i) % qv
            r = abs(cells(i) % centp(s_iaxis))
            volume = cells(i) % vol / r
            ! curl E  GradFi
            cells(i) % res % v(ibe) = cells(i) % res % v(ibe) - &
            qv % v(ieb) * volume
            cells(i) % res % v(ibb + jcol) = cells(i) % res % v(ibb + jcol) + &
            qv % v(ifb) * volume
            cells(i) % dm % e(ibe, ieb) = cells(i) % dm % e(ibe, ieb) + volume
            cells(i) % dm % e(ibb + jcol, ifb) = cells(i) % dm % e(ibb + jcol, ifb) - volume
            ! curl B
            cells(i) % res % v(iee) = cells(i) % res % v(iee) + &
            qv % v(ibb) * volume
            cells(i) % dm % e(iee, ibb) = cells(i) % dm % e(iee, ibb) - volume
            ! no GradFi
            !    cells(i)%res%v(ieb +jcol) = cells(i)%res%v(ieb + jcol) +   &
            !                           qv%v(ife)  *  volume
            ! no DivE=rho/epsi
            !    cells(i)%dm%e(ieb + jcol,ife) = cells(i)%dm%e(ieb + jcol, ife) - volume
        end do
    end subroutine axis_mhd_source

    subroutine external_source(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: current_cell
        type(vector), pointer :: res
        real(rfp) :: volume
        integer :: i, it, nd1
        !
        nd1 = ndim + 1
        !
        do i = 1, size(cells)
            current_cell => cells(i)
            res => current_cell % res
            it = current_cell % itype
            volume = current_cell % vol
            select case(vc(it) % itype)
            case(0) ! no source
                cycle
            case(1) ! moment source
                res % v(imb:ime) = res % v(imb:ime) + vc(it) % vars(:ndim) * volume
                res % v(ien) = res % v(ien) + &
                dot_product(vc(it) % vars(:ndim), current_cell % qv % v(imb:ime)) * volume
            case(2) ! energy source
                !    res%v(ien) = res%v(ien) +   vc(it)%vars(nd1) * volume
                if (current_cell % centp(2)/b_lenref < 1.and.current_cell % centp(1)/b_lenref > 5.) &
                res % v(ien) = res % v(ien) + vc(it) % vars(nd1) * volume
            case(3) ! both
                res % v(imb:ime) = res % v(imb:ime) + vc(it) % vars(:ndim) * volume
                res % v(ien) = res % v(ien) + &
                dot_product(vc(it) % vars(:ndim), current_cell % qv % v(imb:ime)) * volume
                res % v(ien) = res % v(ien) + vc(it) % vars(nd1) * volume
            end select
        end do
    end subroutine external_source

    subroutine blackbody_volumetric_radiation(cells)
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: current_cell
        type(vector), pointer :: res
        type(matrix), pointer :: dm
        real(rfp) :: volume, t
        integer :: i, it
        !  Simple Radiation model
        !  Assume blackbody volumetric radiation to the wall
        !  S = 4SigmaAlpha(Tw^4-T^4)   Tw=0
        !
        do i = 1, size(cells)
            current_cell => cells(i)
            res => current_cell % res
            dm => current_cell % dm
            it = current_cell % itype
            volume = current_cell % vol * 4.0_rfp * stephan_boltzmann
            t = current_cell % qv % v(ien)
            if (t < 1.e4_rfp) cycle
            res % v(ien) = res % v(ien) - t**4 * volume
            dm % e(ien, ien) = dm % e(ien, ien) + 4.0_rfp * volume * t**3
        end do
    end subroutine blackbody_volumetric_radiation

    subroutine growth_source(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: current_cell
        type(vector), pointer :: res, qv
        type(matrix), pointer :: dm
        real(rfp) :: volume, cdest, cprod, zv, pliq, dp, rhoi(nspe)
        real(rfp) :: t, yi(nspe), rholv, uinf, sten, rho, a1, a2, theta, tott, tt, pt
        integer :: i, it, nd1

        !Sublimation Pressure
        !
        ! ln(Ps/Pc)=Tt/T(a1*q+a2*q^2.7)
        ! q=(1-T/Tt), Tt=83.8058K,Pt=68.891kPa
        !
        tt = 83.8058_rfp
        pt = 68891.0_rfp
        a1 = -11.391604_rfp
        a2 = -0.39513431_rfp
        !

        cdest = g_cdest !1.e2_rfp
        cprod = g_cprod !4e2_rfp
        !
        do i = 1, size(cells)
            current_cell => cells(i)
            qv => current_cell % qv
            res => current_cell % res
            dm => current_cell % dm

            yi = yif(qv)
            rhoi = rhoif(qv)
            rholv = rhoi(1) * rhoi(2)
            uinf = sqrt(sum(qv % v(imb:ime)**2))
            t = qv % v(ien)
            if (t >= tt) cycle ! no condense
            tott = t / tt
            theta = one - tott
            pliq = pt * exp((a1 * theta + a2 * theta**2.7_rfp) / tott)
            !
            dp = (qv % v(ico) - pliq + g_pbase) ! rho of liquid
            volume = current_cell % vol
            rho = rhof(qv)
            if (dp < zero) then
                res % v(isb) = res % v(isb) + &
                rho * yi(2) * cdest * abs(dp) * volume

                dm % e(isb, ico) = dm % e(isb, ico) + rho * yi(2) * cdest * volume
                dm % e(isb, isb) = dm % e(isb, isb) + rho * abs(dp) * cdest * volume
            else
                res % v(isb) = res % v(isb) - &
                rho * yi(1) * cprod * abs(dp) * volume
                dm % e(isb, ico) = dm % e(isb, ico) + rho * yi(1) * cprod * volume
                dm % e(isb, isb) = dm % e(isb, isb) + rho * abs(dp) * cprod * volume

            end if

        end do
    end subroutine growth_source

    subroutine liquid_vapor_source(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: current_cell
        type(vector), pointer :: res, qv
        type(vector) :: rhox
        type(matrix), pointer :: dm
        real(rfp) :: volume, cdest, cprod, zv, pliq, dp, rhoi(nspe)
        real(rfp) :: t, yi(nspe), rholv, uinf, sten, rho, a1, a2, theta, dt, tt, pt
        integer :: i, it, nd1

        !Melting Temperature
        !
        tt = g_tref(1) !380.0_rfp
        !
        cdest = g_sref(1) !borrow these variables
        !
        do i = 1, size(cells)
            current_cell => cells(i)
            qv => current_cell % qv
            res => current_cell % res
            dm => current_cell % dm

            yi = yif(qv)
            !   uinf = sqrt(sum(qv%v(imb:ime)**2))
            t = qv % v(ien)
            if (t < tt) cycle ! no vaporization
            if (yi(1) < 1.e-10_rfp) cycle
            !print *,'too big',yi,qv%v(imb:ime)
            dt = t - tt
            volume = current_cell % vol
            rho = rhof(qv)
            pt = rho * yi(1) * cdest * dt * volume
            res % v(isb) = res % v(isb) - pt
            res % v(isb + 1) = res % v(isb + 1) + pt
            cycle
            dm % e(isb, ien) = dm % e(isb, ien) - rho * yi(1) * cdest * volume
            dm % e(isb + 1, ien) = dm % e(isb + 1, ien) + rho * yi(1) * cdest * volume
            dm % e(isb, isb) = dm % e(isb, isb) - rho * dt * cdest * volume
            rhox = rhoxf(qv)
            dm % e(isb,:) = dm % e(isb,:) - rhox % v * yi(1) * cdest * volume
            dm % e(isb + 1,:) = dm % e(isb + 1,:) + rhox % v * yi(1) * cdest * volume
        end do
    end subroutine liquid_vapor_source

    subroutine cavitation_source(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: current_cell
        type(vector), pointer :: res, qv
        type(vector) :: rhox
        type(matrix), pointer :: dm
        real(rfp) :: volume, cdest, cprod, zv, pliq, dp
        real(rfp) :: yi(nspe), rho
        integer :: i, it, nd1

        !
        zv = b_vinf / b_lenref !two / (b_rhoinf * b_vinf * b_lenref)
        cdest = g_sref(1) * zv !borrow these variables
        cprod = g_sref(2) * zv !
        !
        do i = 1, size(cells)
            current_cell => cells(i)
            qv => current_cell % qv
            res => current_cell % res
            dm => current_cell % dm

            pliq = saturationp(qv) ! calculation from saturation line
            yi = yif(qv)
            dp = (qv % v(ico) - pliq + g_pbase) / pliq ! 
            !   dp = pf(qv) / pliq
            volume = current_cell % vol
            rho = rhof(qv)
            rhox = rhoxf(qv)
            if (dp < zero) then ! vaporization
                !   if(dp < one) then   ! vaporization
                dp = abs(dp) !sqrt(abs(dp))
                zv = rho * yi(1) * cdest * (dp) * volume
                res % v(isb) = res % v(isb) - zv
                !   dm%e(isb,ico) = dm%e(isb,ico) + rhox%v(ico) * zv / rho - rho * yi(1) * cdest * volume / pliq
                dm % e(isb, isb) = dm % e(isb, isb) + rhox % v(isb) * zv / rho + rho * (dp * cdest) * volume
            else
                !   dp = sqrt(dp)
                zv = rho * yi(2) * cprod * dp * volume
                res % v(isb) = res % v(isb) + zv
                !   dm%e(isb,ico) = dm%e(isb,ico) - rhox%v(ico) * zv / rho - rho * yi(2) * cprod * volume / pliq
                dm % e(isb, isb) = dm % e(isb, isb) - rhox % v(isb) * zv / rho + rho * dp * cprod * volume
                !    y2 = 1 - y1
            end if
        end do
    end subroutine cavitation_source

    subroutine Coriolis_force(cells)
        !********************************************************
        implicit none
        !   Coriolis force = -rho*OmegaXV
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: cc
        type(vector), pointer :: qv
        real(rfp) :: volume, omega(3), rvel(3), vel(3), rho, rhop, rhot
        real(rfp) :: gm(3, 3), g(3)
        type(matrix), pointer :: dm
        integer :: i
        if (ndim == 2) return ! no 2D 
        !
        omega = zero
        omega(s_iaxis) = s_omega
        gm = cross2mat(omega)
        do i = 1, size(cells)
            cc => cells(i)
            volume = cc % vol
            qv => cc % qv
            dm => cc % dm
            rho = rhof(qv) * volume ! mass
            rhop = rhopf(qv) * volume
            rhot = rhotf(qv) * volume
            vel(:ndim) = qv % v(imb:ime)
            !    g  = acrossb(omega,rvel) + two * acrossb(omega,vel)   !centripetal and coriolis force
            g = acrossb(omega, vel) !centripetal and coriolis force
            cc % res % v(imb:ime) = cc % res % v(imb:ime) - &
            rho * g(:ndim)
            !    cc%res%v(ien) = cc%res%v(ien) - rho * dot_product(vel,g)
            ! moment
            dm % e(imb:ime, ico) = dm % e(imb:ime, ico) + max(rhop * g(:ndim), zero)
            dm % e(imb:ime, ien) = dm % e(imb:ime, ien) + max(rhot * g(:ndim), zero)
            dm % e(imb:ime, imb:ime) = dm % e(imb:ime, imb:ime) + max(rho * gm(:ndim,:ndim), zero)
            !    dm%e(imb:ime,ico) =  dm%e(imb:ime,ico) + rhop * g(:ndim)
            !    dm%e(imb:ime,ien) =  dm%e(imb:ime,ien) + rhot * g(:ndim)
            !    dm%e(imb:ime,imb:ime) =  dm%e(imb:ime,imb:ime) + rho * gm(:ndim,:ndim)
        end do

    end subroutine coriolis_force

    subroutine gravity_source(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: current_cell
        type(vector), pointer :: qv
        real(rfp) :: volume
        type(matrix), pointer :: dm
        integer :: i
        if (g_igrav == 0) return ! no gravity source
        !
        do i = 1, size(cells)
            ! calculate flow filed variables rho,p,u,t,h,etc.
            ! share left cell variables
            current_cell => cells(i)
            volume = current_cell % vol
            qv => current_cell % qv
            dm => current_cell % dm
            dm = dm - dj(qv) * volume ! source term
            current_cell % res = current_cell % res + &
            source_vector(current_cell) * volume
        end do
    end subroutine gravity_source

    !************************************************
    ! source term
    !*****************************************************
    function source_vector(scell)
        implicit none
        type(cell), intent(in) :: scell
        type(vector) :: source_vector, qv
        real(rfp) :: g(ndim), rho, vel(ndim)
        qv = scell % qv
        vel = qv % v(imb:ime)
        g = g_gravity
        rho = rhof(qv)
        !
        source_vector % v = zero
        source_vector % v(imb:ime) = rho * g(:)
        source_vector % v(ien) = rho * dot_product(vel, g)
        !
    end function source_vector

    !***********************************************************************
    !  calculate source jacobian
    !***********************************************************************
    function dj(qv)

        type(matrix) :: dj
        type(vector), intent(in) :: qv
        real(rfp) :: rhop, rhot, rho, vel(ndim)
        real(rfp) :: g(ndim)
        ! initial dj
        g = g_gravity
        rho = rhof(qv)
        rhop = rhopf(qv)
        rhot = rhotf(qv)
        vel = qv % v(imb:ime)
        dj % e = zero
        ! 
        ! momentum row
        dj % e(imb:ime, ico) = rhop * g(:)
        dj % e(imb:ime, ien) = rhot * g(:)
        !energy row
        dj % e(ien, ico) = rhop * dot_product(g, vel)
        dj % e(ien, imb:ime) = rho * g(:)
        dj % e(ien, ien) = rhot * dot_product(g, vel)

    end function dj

    subroutine axis_source(cells, nodes)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:), cc
        type(node), pointer :: nodes(:)
        ! local variables
        type(vector), pointer :: qv
        real(rfp) :: volume, rho, zmu, r, d, div, dvdr
        integer :: i, jcol
        if (s_iaxis == 0) return
        if (nmeq == 0) return
        jcol = ico + s_iaxis
        !
        do i = 1, size(cells)
            qv => cells(i) % qv
            r = abs(cells(i) % centp(s_iaxis))
            volume = cells(i) % vol / r
            cells(i) % res % v(jcol) = cells(i) % res % v(jcol) + &
            qv % v(ico) * volume
            cells(i) % dm % e(jcol, ico) = cells(i) % dm % e(jcol, ico) - volume
            !  -Mu*vr/r^2
            if (s_ivis > 0) then
                zmu = mumix(qv)
                if (s_ivis == 2) zmu = zmu + mu_turb(qv)
                zmu = zmu * volume
                cc => cells(i)
                call vis_add(cc, div, dvdr, nodes)
                !dvdr=cc%gradient(iaub,s_iaxis)
                cells(i) % res % v(jcol) = cells(i) % res % v(jcol) - &
                !      zmu * qv%v(jcol) / r 
                zmu * (two * qv % v(jcol) / r - two_third * div)
                cells(i) % dm % e(jcol, jcol) = cells(i) % dm % e(jcol, jcol) + two * zmu / r
            end if

            if (s_iaux == 1) then ! swirl flow
                rho = rhof(qv)
                cells(i) % res % v(jcol) = cells(i) % res % v(jcol) + &
                rho * qv % v(iaub) * qv % v(iaub) * volume ! vq*2/r
                cells(i) % res % v(iaub) = cells(i) % res % v(iaub) - &
                rho * qv % v(iaub) * qv % v(jcol) * volume ! vq*vr/r
                if (s_ivis > 0) then
                    cells(i) % res % v(iaub) = cells(i) % res % v(iaub) + &
                    zmu * (dvdr - qv % v(iaub) / r)
                    cells(i) % dm % e(iaub, iaub) = cells(i) % dm % e(iaub, iaub) + zmu / r
                end if
                cells(i) % dm % e(jcol, iaub) = cells(i) % dm % e(jcol, iaub) - two * rho * qv % v(iaub) * volume
                cells(i) % dm % e(iaub, iaub) = cells(i) % dm % e(iaub, iaub) + rho * qv % v(jcol) * volume
                cells(i) % dm % e(iaub, jcol) = cells(i) % dm % e(iaub, jcol) + rho * qv % v(iaub) * volume
            end if
        end do
    end subroutine axis_source

    subroutine vis_add(cc, div, dvdr, nodes)
        type(cell), pointer :: cc
        type(face), pointer :: cf
        type(node), pointer :: nodes(:)
        type(vector) :: qv
        real(rfp) :: div, dvdr, vecn(ndim)
        integer :: i
        div = zero; dvdr = zero;
        do i = 1, size(cc % sface)
            cf => cc % sface(i) % to_face
            vecn = cf % area * cf % vecn
            if (associated(cf % right_cell, cc)) vecn = -vecn
            qv = (nodes(cf % f2n(1)) % qv + nodes(cf % f2n(2)) % qv) * half
            div = div + sum(vecn * qv % v(imb:ime))
            ! if(naux > 0) dvdr = dvdr + vecn(s_iaxis) / cf%centp(s_iaxis) * qv%v(iaub)
            if (naux > 0) dvdr = dvdr + vecn(s_iaxis) * qv % v(iaub)
        end do
        div = div / cc % vol
        dvdr = dvdr / cc % vol - cc % qv % v(iaub) / cc % centp(s_iaxis)
    end subroutine vis_add

    subroutine radiation_source1(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: current_cell
        type(vector), pointer :: qv
        real(rfp) :: volume, t, s, dsdi
        type(matrix), pointer :: dm
        integer :: i
        if (nrad == 0) return ! no radiation source
        !
        do i = 1, size(cells)
            current_cell => cells(i)
            qv => current_cell % qv
            dm => current_cell % dm
            volume = current_cell % vol * absorption_coef(qv)
            t = 13889. - (13889. - 2778.) * current_cell % centp(2) / b_lenref / 0.1
            s = (qv % v(irad) - 4._rfp * stephan_boltzmann * t**4) * volume
            dm % e(irad, irad) = dm % e(irad, irad) - volume ! source term
            current_cell % res % v(irad) = current_cell % res % v(irad) + s

        end do
    end subroutine radiation_source1

    subroutine radiation_source(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: current_cell
        type(vector), pointer :: qv
        real(rfp) :: volume, t, s, dsdi
        type(matrix), pointer :: dm
        integer :: i
        if (nrad == 0) return ! no radiation source
        !
        do i = 1, size(cells)
            current_cell => cells(i)
            qv => current_cell % qv
            dm => current_cell % dm
            volume = current_cell % vol * absorption_coef(qv)
            t = qv % v(ien)
            dsdi = 16._rfp * stephan_boltzmann * t**3 * volume ! source term
            s = (qv % v(irad) - 4._rfp * stephan_boltzmann * t**4) * volume
            !
            dm % e(irad, irad) = dm % e(irad, irad) - volume ! source term
            dm % e(irad, ien) = dm % e(irad, ien) + dsdi
            current_cell % res % v(irad) = current_cell % res % v(irad) + s
            ! energy equation
            !    current_cell%res%v(ien) = current_cell%res%v(ien) +   s
            !    dm%e(ien,irad) = dm%e(ien,irad) -  volume ! source term
            !    dm%e(ien,ien) = dm%e(ien,ien) + dsdi

        end do
    end subroutine radiation_source

    subroutine propane_global_reaction_rate(cells)
        type(cell), pointer :: cells(:)
        type(vector), pointer :: qv, res
        type(vector) :: rhox
        type(matrix), pointer :: dm
        real(rfp) :: rate(nspe)
        integer :: i, j, k, l
        real(rfp) :: xcon(nspe), yi(nspe), gk(nspe), dg
        real(rfp) :: prodf, prodb, sumf, sumb, suma, kf, kb, kc, q, qf, qb, af(3), ab(3)
        real(rfp) :: rho, rhoinv, rhoy(nspm1), rhop, rhot, t, tinv, xc
        real(rfp) :: a, b, c, e, pd, td, yd, m, n, volume
        !
        !C3H8 O2 H2O CO2 N2
        !C3H8+5O2 -> 4H2O + 3CO2
        !
        rate = zero
        do i = 1, size(cells)
            qv => cells(i) % qv
            res => cells(i) % res
            dm => cells(i) % dm
            volume = cells(i) % vol

            rho = rhof(qv)
            t = qv % v(ien)
            yi = yif(qv)
            tinv = one / t
            xcon = rho * yi / g_mwi
            !
            !  a*T^b*exp{E/T}
            af(1) = 8.6e11_rfp
            af(2) = zero
            af(3) = -15098.0_rfp
            kf = af(1) * exp(af(3) * tinv)
            !
            m = 0.1_rfp; n = 1.65_rfp
            !
            ! the rate of progress
            qf = kf * xcon(1) **m * xcon(2)**n
            rate(1) = -qf * g_mwi(1)
            rate(2) = -5.0_rfp * qf * g_mwi(2)
            rate(3) = 4.0_rfp * qf * g_mwi(3)
            rate(4) = 3.0_rfp * qf * g_mwi(4)
            rate = rate * volume
            do j = 1, nspm1
                res % v(isb + j - 1) = res % v(isb + j - 1) + rate(j)
                !
                !  Jacobian of reaction rate only count temperature
                !
                dm % e(isb + j - 1, ien) = dm % e(isb + j - 1, ien) + rate(j) * af(3) * tinv * tinv
            end do

        end do
    end subroutine propane_global_reaction_rate

    subroutine mhd_setup(cells)
        type(cell), pointer :: cells(:), cc
        logical :: iexist
        integer :: i
        !
        iexist = .true.
        do while (iexist)
            inquire(file = 'from_mhd', exist = iexist)
        end do
        open(2, file = 'from_mhd')
        do i = 1, size(cells)
            if (vc(cells(i) % itype) % itype == 0) then ! plasma zone
                !  write(2)cells(i)%lorenz_force(:ndim),cells(i)%je
            end if
        end do
        close(2)
        do while (.not.iexist)
            inquire(file = 'to_mhd', exist = iexist)
        end do
        open(3, file = 'to_mhd')
        do i = 1, size(cells)
            if (vc(cells(i) % itype) % itype == 0) then ! plasma zone
                !  read(3)cells(i)%vel(:ndim)
            end if
        end do
        close(3)
        call unlink('to_mhd')

    end subroutine mhd_setup

end MODULE dfd_source
