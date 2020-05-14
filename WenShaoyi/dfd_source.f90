MODULE dfd_source
    USE dfd_state
    USE dfd_data
    USE dfd_library
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
        if (s_ilevelset > 0) call levelset_source(cells)
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
            gl = liquid_fraction(qv) + 1.e - 20_rfp ! avoid singularity
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
                if (s_ilevelset > 0) then ! level set
                    h = heavisidef(qv)
                    res % v(ien) = res % v(ien) + un * (1 - h) * (one - gl) * ((metal % cps - metal % cpl) * (qvf % v(ien) - metal % ts) - metal % lheat)
                    !    res%v(ien) = res%v(ien) + un * (1-h) * (one - gl) * metal%cps * qvf%v(ien)
                    ! print *, qvf%v(ien)
                else
                    res % v(ien) = res % v(ien) + un * (one - gl) * ((metal % cps - metal % cpl) * (qvf % v(ien) - metal % ts) - metal % lheat)
                    !    res%v(ien) = res%v(ien) + un * (one - gl) * metal%cps * qvf%v(ien)
                end if
            end do
        end do
        ! 
    end subroutine solid_melting
    !
    subroutine levelset_source(cells)
        !********************************************************
        implicit none
        !
        type(cell), pointer :: cells(:)
        ! local variables
        type(cell), pointer :: cc
        type(face), pointer :: cf
        type(vector), pointer :: res, qv
        type(vector) :: qvf
        real(rfp) :: r, us(ndim), x0(ndim), x_s, x_e, fp, volume, vp, ls_co_v, ls_e, rho, h, rhoi(nspe), hi(nspe), laser, ah, tamb, emi, laser_rad, laser_p
        real(rfp) :: vecn(ndim), kc, delta
        integer :: i, it, j

        real(rfp) :: layertime, passtime, time_thispass
        integer :: layerno, passno, id, ierr
        !  return  ! not ready yet
        !

        call mpi_comm_rank(mpi_comm_world, id, ierr)

        ah = 10.
        tamb = 300.
        emi = 0.3

        !print *,'fp',powderflow1%A / (2*laser1%rb/laser1%vx)

        passtime = laser1 % pass_length / laser1 % vx + laser1 % time_between_passes
        layertime = passtime * laser1 % passes_per_layer

        layerno = s_elapsed_time/layertime + 1
        passno = (s_elapsed_time - (layerno - 1) * layertime)/passtime + 1
        time_thispass = s_elapsed_time - (layerno - 1) * layertime - (passno - 1) * passtime

        if (time_thispass > laser1 % pass_length / laser1 % vx) then
            laser1 % x = 1.e6
            laser1 % y = 1.e6
            laser1 % z = 1.e6

            if (id == 0) then
                print *, 'laser position--->out of plane!!!'
                print *, 'layerno', layerno, 'passno', passno
            end if
        else
            if (laser1 % switch > 0.5 .and. mod(layerno, 2) == 0) then
            if (id == 0) then
                print *, 'switching!!!'
            end if

            if (laser1 % zigzag > 0.5 .and. mod(passno, 2) == 0) then
                laser1 % x = laser1 % x0 + (passno - 1) * laser1 % patch_dist
                laser1 % y = 0 + (layerno - 1) * laser1 % layer_thickness
                laser1 % z = laser1 % z0 + laser1 % pass_length - time_thispass * laser1 % vx

                if (id == 0) then
                    print *, 'laser position--->backwards'
                    print *, laser1 % x, laser1 % y, laser1 % z
                    print *, 'layerno', layerno, 'passno', passno
                end if
            else
                laser1 % x = laser1 % x0 + (passno - 1) * laser1 % patch_dist
                laser1 % y = 0 + (layerno - 1) * laser1 % layer_thickness
                laser1 % z = laser1 % z0 + time_thispass * laser1 % vx

                if (id == 0) then
                print *, 'laser position--->normal'
                print *, laser1 % x, laser1 % y, laser1 % z
                print *, 'layerno', layerno, 'passno', passno
                end if
            end if
        else
            if (laser1 % zigzag > 0.5 .and. mod(passno, 2) == 0) then
            laser1 % x = laser1 % x0 + laser1 % pass_length - time_thispass * laser1 % vx
            laser1 % y = 0 + (layerno - 1) * laser1 % layer_thickness
            laser1 % z = laser1 % z0 + (passno - 1) * laser1 % patch_dist

            if (id == 0) then
                print *, 'laser position--->backwards'
                print *, laser1 % x, laser1 % y, laser1 % z
                print *, 'layerno', layerno, 'passno', passno
            end if
        else
            laser1 % x = laser1 % x0 + time_thispass * laser1 % vx
            laser1 % y = 0 + (layerno - 1) * laser1 % layer_thickness
            laser1 % z = laser1 % z0 + (passno - 1) * laser1 % patch_dist

            if (id == 0) then
            print *, 'laser position--->normal'
            print *, laser1 % x, laser1 % y, laser1 % z
            print *, 'layerno', layerno, 'passno', passno
            end if
            end if
            end if
        end if



        do i = 1, size(cells)
            cc => cells(i)
            res => cc % res
            qv => cc % qv
            volume = cc % vol
            ! Powder addition speed fp

            r = sqrt((cc % centp(1) - laser1 % x)**2 + (cc % centp(3) - laser1 % z)**2)
            if (r > powderflow1 % rp) then
                fp = 0.

            else
                !flow rate/density / cross-section area (pay attention to the unit)
                !    fp = powderflow1%A * sqrt(1-r/powderflow1%rp) !0.01 !m/s 
                fp = powderflow1 % A / (2 * laser1 % rb/laser1 % vx) * (1 - r/powderflow1 % rp)**0.5

            end if

            cc % Vinterface = fp
            if (time_thispass > laser1 % pass_length / laser1 % vx) cc % Vinterface = 0


            !if (cc%centp(1)>0.0013 .and. cc%centp(1)<0.0015 .and. cc%centp(3)>0.0003 .and. cc%centp(3)<0.0004 .and. cc%centp(2)>-30e-6 .and. cc%centp(2)<0) then
            !if (cc%vinterface>0 .and. cc%centp(2)>-30e-6 .and. cc%centp(2)<0) then
            !    print *,'pt in',id
            !    print *,cc%centp(:)
            !    print *,fp,r
            !    print *,cc%Vinterface

            !end if

            laser = laser1 % A * laser1 % p / pi /(laser1 % rb)**2 * metal % alpha * exp(-laser1 % A * (r/ laser1 % rb)**2) - powderflow1 % h * (qv % v(ien) - tamb) - stephan_boltzmann * metal%emi * (qv%v(ien)**4-tamb**4)



            !  print *, sqrt(sum(cc%gradient(ils,:)**2))
            rho = rhof(qv)
            rhoi = rhoif(qv)
            hi = hif(qv)
            h = hf(qv)
            ls_co_v = (rhoi(1) - dummy % rhol) * fp * DiracDeltaf(qv)
            ls_e = fp * DiracDeltaf(qv)
            res % v(ico) = res % v(ico) + ls_co_v * volume
            res % v(imb:ime) = res % v(imb:ime) + ls_co_v * qv % v(imb:ime) * volume
            res % v(ien) = res % v(ien) + ls_e * (rho * (hi(1) - dummy % cpl * qv % v(ien)) + h * (rhoi(1) - dummy % rhol)) * volume + laser * DiracDeltaf(qv) * volume 
            res % v(ils) = res % v(ils) - fp * rhof(qv) * volume + ls_co_v * qv % v(ils) * volume

            ! calculate curvature
            kc = cc % curvature
            vecn = cc % gradient(ils,:)
            vecn = vecn / sqrt(sum(vecn**2)) ! normilized

            !
            delta = DiracDeltaf(qv) * volume
            ! capillary and managonia force

            ! res%v(imb:ime) = res%v(imb:ime) - (metal%sigma * vecn * kc) * delta 
            res % v(imb:ime) = res % v(imb:ime) - (-(cc % gradient(ien,:) - sum(cc % gradient(ien,:) * vecn) * vecn) * metal % sigmat/(one)) * delta
            ! wenda replaced (one-heavisidef(qv)+1.e-15_rfp) with (one)

            !  print *,'I am here'
            !
        end do
        ! 
    end subroutine levelset_source

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
            if (yi(1) < 1.e - 10_rfp) cycle
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


    subroutine levelset_calculation(cells, nadv)
        implicit none
        !
        type(cell), pointer :: cells(:), cc, gc
        integer :: nadv

        integer :: i, j
        real(rfp) :: dt, oldls, vvv, ctls, ctcell

        !print *,'calculating ls'

        dt = one / s_dt_inv

        do i = 1, size(cells)
            cc => cells(i)
            !oldls = cc%qv%v(ils)    !
            oldls = cells(i) % qvn(gear_in(1)) % v(ils)
            vvv = cc % Vinterface!2

            ctls = 0
            ctcell = 0
            do j = 1, size(cc % scell)
                gc => cc % scell(j) % to_cell
                !if (abs(gc%centp(2)-cc%centp(2))>1e-6) then
                ctls = ctls + gc % qvn(gear_in(1)) % v(ils)
                ctcell = ctcell + 1.
                !end 
            end do
            if (ctcell > 0) then
                oldls = ctls/ctcell
            else
                print *, 'error: ctcell=0'
                stop
            end if

            cc % qv % v(ils) = cells(i) % qvn(gear_in(1)) % v(ils)
            if (vvv > 1.e - 6 .or. vvv < -1.e - 6) then
                cc % qv % v(ils) = cells(i) % qvn(gear_in(1)) % v(ils) - vvv * dt
                if (cc % qv % v(ien) > metal % tl) cc % qv % v(ils) = oldls - vvv * dt
            end if
        end do

    end subroutine levelset_calculation

    subroutine calculate_vinterface(cells, nadv, id)
        type(cell), pointer :: cells(:), cc, gc
        integer :: nadv, id

        integer :: i, j, flag
        real(rfp) :: vecn1(ndim), vecn2(ndim), vtemp, totalvec(ndim), totalnb, localvec(ndim), ls1, ls2, relax

        relax = 0.25

        do i = 1, size(cells)
            cc => cells(i)
            flag = 1

            totalvec(:) = 0
            totalnb = 0
            do j = 1, size(cc % scell)
                gc => cc % scell(j) % to_cell
                if (cc % qvn(gear_in(1)) % v(ils) < 0 .and. gc % qvn(gear_in(1)) % v(ils) > 0) then
                    !                flag=1
                    totalnb = totalnb + 1.0
                    vecn1 = 0.5 * cc % qv % v(imb:ime) !*0.25+ 0.25*cc%qvn(gear_in(1))%v(imb:ime) + 0.5*cc%qvn(gear_in(2))%v(imb:ime)
                    vecn2 = 0.5 * gc % qv % v(imb:ime) !*0.25 + 0.25*gc%qvn(gear_in(1))%v(imb:ime) + 0.5*gc%qvn(gear_in(2))%v(imb:ime)
                    ls1 = cc % qvn(gear_in(1)) % v(ils)
                    ls2 = gc % qvn(gear_in(1)) % v(ils)
                    localvec(:) = vecn1(:) !+ (0-ls1)/(ls2-ls1)*(vecn2-vecn1)
                    totalvec(:) = totalvec(:) + localvec(:)
                end if

            end do

            if (totalnb > 0) then
                vecn1 = cc % gradient(ils,:)
                vecn1 = vecn1 / sqrt(sum(vecn1**2))
                vecn2 = totalvec(:)/totalnb
                vtemp = dot_product(vecn1, vecn2)
                cc % Vinterface2 = cc % vinterface + vtemp
            else
                cc % Vinterface2 = 0
            end if

        end do


    end subroutine calculate_vinterface

    subroutine Vint_extension(cells, faces, id, nadv)
        implicit none
        type(face), pointer :: faces(:), cf
        type(cell), pointer :: cells(:), cc
        integer :: id, nadv

        integer :: i, j, k, ii, jj, kk, ierr, flag, ctpt
        real(rfp) :: dist, dist2, total1, ct, ct2, cosine, sine, length, dl
        real(rfp) :: lsgrad(3), direct2(3)

        character(len = 40) :: fn

        var1_ext = 0
        var2_ext = 0
        var3_ext = 0
        do i = 1, isize_ext
            do j = 1, jsize_ext
                do k = 1, ksize_ext
                    if (map_ext(i, j, k) > 0) then
                        var1_ext(i, j, k) = cells(map_ext(i, j, k)) % qv % v(ils)
                        var2_ext(i, j, k) = cells(map_ext(i, j, k)) % vinterface2
                        !print *,'readls',var1_Ext(i,j,k),yc_ext(j)
                    end if
                end do
            end do
        end do
        call mpi_allreduce(var2_ext, var3_ext, isize_ext * jsize_ext * ksize_ext, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
        call mpi_allreduce(var1_ext, var2_ext, isize_ext * jsize_ext * ksize_ext, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)

        do i = 1, isize_ext
            do j = 1, jsize_ext
                do k = 1, ksize_ext
                    ii = i
                    jj = j
                    kk = k
                    if (i == 1) ii = i + 1
                    if (i == isize_ext) ii = i - 1
                    if (j == 1) jj = j + 1
                    if (j == jsize_ext) jj = j - 1
                    if (k == 1) kk = k + 1
                    if (k == ksize_ext) kk = k - 1

                    var2_ext(i, j, k) = var2_ext(ii, jj, kk)
                    var3_ext(i, j, k) = var3_ext(ii, jj, kk)

                end do
            end do
        end do


        !var1=scan
        !var2=levelset
        !var3=vinterface

        if (id == 0) then
            dl = dxc_ext
            flag = 1
            dist = 0

            var1_ext = 0

            ctpt = 0
            do i = 1, isize_ext
                do j = 1, jsize_ext
                    do k = 1, ksize_ext
                        if (var2_ext(i, j, k) < 0) then
                            flag = 0
                            do ii = i - 1, i + 1
                                do jj = j - 1, j + 1
                                    do kk = k - 1, k + 1
                                        dist = abs(ii - i) + abs(jj - j) + abs(kk - k)
                                        if (ii < 1 .or. ii > isize_ext .or. jj < 1 .or. jj > jsize_ext .or. kk < 1 .or. kk > ksize_ext) then
                                            ! do nothing
                                        else if (dist < 1.1 .and. var2_ext(ii, jj, kk) >= 0) then
                                            flag = 1
                                            ! print *,i,j,k,'vs',ii,jj,kk
                                            ! print *,var2_ext(i,j,k),'vs',yc_ext(j)
                                            ! print *,var2_ext(ii,jj,kk),'vs',yc_ext(jj)
                                        end if
                                    end do
                                end do
                            end do
                            if (flag == 1) then
                                var1_ext(i, j, k) = 10
                                ctpt = ctpt + 1
                                !print *,'pt on kw',i,j,k,ctpt,var1_ext(i,j,k),var2_ext(i,j,k)
                            end if
                        end if

                    end do
                end do
            end do

            var4_ext = 0
            do i = 1, isize_ext
                do j = 1, jsize_ext
                    do k = 1, ksize_ext
                        if (var1_Ext(i, j, k) > 5) then
                            ct = 0
                            ct2 = 0

                            do ii = i - 2, i + 2
                                do jj = j - 2, j + 2
                                    do kk = k - 2, k + 2
                                        dist = abs(ii - i) + abs(jj - j) + abs(kk - k)
                                        if (ii < 1 .or. ii > isize_ext .or. jj < 1 .or. jj > jsize_ext .or. kk < 1 .or. kk > ksize_ext) then
                                            ! do nothing
                                        else if (dist > 0 .and. var1_ext(ii, jj, kk) > 5) then
                                            ct = ct + 1./dist
                                            ct2 = ct2 + var3_ext(ii, jj, kk)/dist
                                        end if
                                    end do
                                end do
                            end do

                            if (ct <= 0) then
                                print *, 'error: ct==0 in averaging var3'
                                stop
                            end if

                            var4_Ext(i, j, k) = ct2/ct

                        end if
                    end do
                end do
            end do

            do i = 1, isize_ext
                do j = 1, jsize_ext
                    do k = 1, ksize_ext
                        var3_Ext(i, j, k) = var4_Ext(i, j, k) !0.25*var3_ext(i,j,k)+0.75*
                    end do
                end do
            end do

            dist = 0
            do while (flag == 1)

                if (dist <= 5) dist = dist + 0.5
                if (dist > 5) dist = dist + 0.5

                if (dist > 100) then
                    print *, 'too large dist in extension'
                    stop
                end if

                flag = 0
                if (dist < 5.) flag = 1

                ctpt = 0

                !print *,'dist show',dist,dist*dl
                do i = 1, isize_ext
                    do j = 1, jsize_ext
                        do k = 1, ksize_ext
                            if (abs(var2_ext(i, j, k)) < (dist * dl) .and. var1_ext(i, j, k) == 0) then
                                !print *,i,j,k,var2_ext(i,j,k),var1_ext(i,j,k)

                                flag = 1

                                ct = 0
                                ct2 = 0
                                total1 = 0

                                lsgrad = (/0, 1, 0/)
                                if (i > 1 .and. i < isize_ext) lsgrad(1) = (var2_ext(i + 1, j, k) - var2_ext(i - 1, j, k))/(xc_ext(i + 1) - xc_ext(i - 1))
                                if (j > 1 .and. j < jsize_ext) lsgrad(2) = (var2_ext(i, j + 1, k) - var2_ext(i, j - 1, k))/(yc_ext(j + 1) - yc_ext(j - 1))
                                if (k > 1 .and. k < ksize_ext) lsgrad(3) = (var2_ext(i, j, k + 1) - var2_ext(i, j, k - 1))/(zc_ext(k + 1) - zc_ext(k - 1))

                                length = sqrt(sum(lsgrad**2))
                                lsgrad = lsgrad/length

                                !first find nb pts with their direct2 aligning best with lsgrad
                                do ii = i - 1, i + 1
                                    do jj = j - 1, j + 1
                                        do kk = k - 1, k + 1
                                            if (ii < 1 .or. ii > isize_ext .or. jj < 1 .or. jj > jsize_ext .or. kk < 1 .or. kk > ksize_ext) then
                                                ! do nothing
                                            else if (var1_ext(ii, jj, kk) > 9) then
                                                direct2(1) = ii - i
                                                direct2(2) = jj - j
                                                direct2(3) = kk - k
                                                length = sqrt(sum(direct2**2))
                                                direct2 = direct2/length
                                                cosine = abs(dot_product(lsgrad, direct2))
                                                sine = sqrt(1 - cosine * cosine)
                                                dist2 = length * sine
                                                dist2 = max(dist2, 2e - 6)
                                                ct = ct + 1/dist2
                                                ct2 = ct2 + 1.
                                                total1 = total1 + var3_ext(ii, jj, kk)/dist2
                                            end if
                                        end do
                                    end do
                                end do
                                if (ct > 0 .and. ct2 > 3) then
                                    var3_ext(i, j, k) = total1/ct
                                    var1_ext(i, j, k) = 5
                                    ctpt = ctpt + 1
                                    !write (43,*) i,j,k,ctpt
                                else
                                    !print *,i,j,k,'skip',dist                        
                                end if

                            end if
                        end do
                    end do
                end do

                !       print *,'ct of dist',dist,ctpt

                do i = 1, isize_ext
                    do j = 1, jsize_ext
                        do k = 1, ksize_ext
                            if (var1_ext(i, j, k) > 1) then
                                !print *,i,j,k                
                                var1_ext(i, j, k) = 10
                            end if
                        end do
                    end do
                end do

            end do

            !    print *,'finish extenstion at',dist

            !close(43)

        else
            var1_ext = 0
            var2_ext = 0
            var3_ext = 0
        end if

        call mpi_allreduce(var3_ext, var2_ext, isize_ext * jsize_ext * ksize_ext, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)

        !var2=vinterface

        do i = 1, isize_ext
            do j = 1, jsize_ext
                do k = 1, ksize_ext

                    if (map_ext(i, j, k) > 0) then
                        cells(map_ext(i, j, k)) % vinterface2 = var2_ext(i, j, k)
                    end if
                end do
            end do
        end do

    end subroutine vint_extension

end MODULE dfd_source
