MODULE dfd_state
    USE dfd_type
    USE dfd_data
    USE dfd_library
    USE RPTBDB_setup, ONLY: Qvar, QLocate, p2FttRootCell
    implicit none
contains
    ! 
    !***********************************************************************
    ! calculate flow properties (update every time step)
    !***********************************************************************
    !
    ! Density-------

    function rhof(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: rhof, r, rhoi(nspe), yi(nspe), h
        rhoi = rhoif(qv)
        yi = yif(qv)
        rhof = one / sum(yi/rhoi)
        if (s_ilevelset > 0) then ! level set
            h = heavisidef(qv)
            rhof = (one - h) * rhof + h * dummy % rhol ! 
        end if
        !
    end function rhof

    function rhopf(qv)
        type(vector) :: qv
        ! real(rfp)::rhopf
        real(rfp) :: rhopf, rhopi(nspe), yi(nspe), rhoi(nspe), rho, wrho(nspe), h
        rhopi = rhopif(qv)
        yi = yif(qv)
        rhoi = rhoif(qv)
        rho = one / sum(yi / rhoi)
        wrho = rho / rhoi
        wrho = wrho * wrho * yi
        !   rhopf  =  rho * rho * sum(yi * rhopi /(rhoi* rhoi))
        rhopf = sum(wrho * rhopi)
        if (s_ilevelset > 0) then ! level set
            h = heavisidef(qv)
            rhopf = (one - h) * rhopf + h * dummy % beta ! beta is rhop for dummy
        end if
    end function rhopf


    function rhotf(qv)
        type(vector) :: qv
        real(rfp) :: rhotf, rhoti(nspe), yi(nspe), rhoi(nspe), rho, wrho(nspe), h
        rhoti = rhotif(qv)
        yi = yif(qv)
        rhoi = rhoif(qv)
        rho = one / sum(yi / rhoi)
        !   rhotf  =  rho * rho * sum(yi * rhoti /(rhoi* rhoi))
        wrho = rho / rhoi
        wrho = wrho * wrho * yi
        rhotf = sum(wrho * rhoti)
        if (s_ilevelset > 0) then ! level set
            h = heavisidef(qv)
            rhotf = (one - h) * rhotf + h * dummy % alpha ! alpha is rhot for dummy
        end if
    end function rhotf

    function rhoxf(qv)
        type(vector) :: qv, rhoxf
        real(rfp) :: rho, rhop, coef, rhot, h
        real(rfp) :: rhoti(nspe), rhopi(nspe), yi(nspe), rhoi(nspe), wrho(nspe)
        yi = yif(qv)
        rhoi = rhoif(qv)
        rho = one / sum(yi / rhoi)
        !   rhot  =  rho * rho * sum(yi * rhoti /(rhoi* rhoi))
        !   rhop  =  rho * rho * sum(yi * rhopi /(rhoi* rhoi))
        wrho = rho / rhoi
        wrho = wrho * wrho * yi
        rhoxf % v = zero
        if (nceq > 0) then
            !    rhopi = rhopif(qv)
            !    rhop = sum(wrho * rhopi)
            rhoxf % v(ico) = rhopf(qv)
        end if
        if (neeq > 0) then
            !     rhoti = rhotif(qv)
            !     rhot = sum(wrho * rhoti)
            rhoxf % v(ien) = rhotf(qv)
        end if
        !
        if (nspe > 1) then ! multi-specis
            rhoxf % v(isb:ise) = -rho * rho * (one / rhoi(:nspm1) - one / rhoi(nspe))
        end if
        !
        if (s_ilevelset > 0) then ! level set
            h = diracdeltaf(qv)
            rhoxf % v(ils) = h * (dummy % rhol - rho)
        end if
        !
        ! turbulence
        !
        return
        !  if(s_ivis == 2) then
        !   coef = one / ( one + two_third * qv%v(ike) * rhop)
        !   rhoxf%v(ike) = -two_third * rho * rhop
        !   rhoxf%v(:) = rhoxf%v(:) * coef  
        !  end if

    end function rhoxf
    !
    ! Entrapy ---------
    !

    function hf(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: hf, hi(nspe), yi(nspe), h
        hi = hif(qv)
        yi = yif(qv)
        hf = sum(hi * yi)
        if (s_ilevelset > 0) then ! level set
            h = heavisidef(qv)
            hf = (one - h) * hf + h * dummy % cpl * qv % v(ien)
        end if
    end function hf

    function hpf(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: hpf, hpi(nspe), yi(nspe), h
        hpi = hpif(qv)
        yi = yif(qv)
        hpf = sum(hpi * yi)

        if (s_ilevelset > 0) then ! level set
            h = heavisidef(qv)
            hpf = (one - h) * hpf
        end if
    end function hpf

    function htf(qv)
        type(vector) :: qv
        real(rfp) :: htf, hti(nspe), yi(nspe), h
        hti = htif(qv)
        yi = yif(qv)
        htf = sum(hti * yi)
        if (s_ilevelset > 0) then ! level set
            h = heavisidef(qv)
            htf = (one - h) * htf + h * dummy % cpl
        end if
    end function htf

    function cvf(qv)
        type(vector) :: qv
        real(rfp) :: cvf

        ! Cv = Cp - 1/rho * Dp/DT=Cp + rhoT/rhoP * (1/rho)
        cvf = htf(qv) + rhotf(qv) / (rhopf(qv) * rhof(qv))
        !
    end function cvf

    function hxf(qv)
        type(vector) :: qv, hxf
        real(rfp) :: hi(nspe), h, yi(nspe)
        hxf % v = zero
        if (nceq > 0) hxf % v(ico) = hpf(qv)
        if (nmeq > 0) hxf % v(imb:ime) = qv % v(imb:ime)
        if (neeq > 0) hxf % v(ien) = htf(qv)

        !  if(s_ivis == 2) hxf%v(ike) = five_third
        !
        if (nspe > 1) then ! multi-specis
            hi = hif(qv)
            hxf % v(isb:ise) = hi(:nspm1) - hi(nspe)
        end if
        !
        if (s_ilevelset > 0) then ! level set
            h = diracdeltaf(qv)
            hi = hif(qv)
            yi = yif(qv)
            hxf % v(ils) = h * (dummy % cpl * qv % v(ien) - sum(hi * yi))
        end if
    end function hxf

    !
    ! Entropy
    !
    function entropyf(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: entropyf, si(nspe), yi(nspe), h
        integer :: i, j, n
        si = sif(qv)
        yi = yif(qv)
        entropyf = sum(si * yi)
        if (s_ilevelset > 0) then ! level set
            h = heavisidef(qv)
            si = (one - h) * si ! dummy entropy is zero
        end if
    end function entropyf

    function sif(qv) ! entropy of each species
        type(vector), intent(in) :: qv
        real(rfp) :: p, t, t_inv, tln, sif(nspe), yi(nspe), ri(nspe)
        integer :: i, j, n
        type(Qvar) :: pro


        select case(s_irealgas)
        case(0) ! perfect gas
            p = abs(qv % v(ico) + g_pbase) / g_patm
            t = abs(qv % v(ien)) / g_troom
            ri = g_runi / g_mwi
            sif = htif(qv) * log(t) - Ri * log(p)
        case(1)
            do i = 1, nspe
                sif(i) = inter_prop(fluid(i) % s, qv, i)
            end do
        case(2, 3)
            p = abs(qv % v(ico) + g_pbase) / g_patm
            t = qv % v(ien)
            do i = 1, nspe
                n = species(i) % interv
                t = min(max(species(i) % Titv(1), t), species(i) % Titv(n + 1))
                t_inv = one / (t * t)
                tln = log(t)
                do j = 1, n
                    if (t <= species(i) % Titv(j + 1)) then
                        sif(i) = polynomial(t, species(i) % ssp(:7, j)) * t_inv
                        sif(i) = species(i) % R * (sif(i) + species(i) % ssp(8, j) * tln - log(p))
                        exit
                    end if
                end do
            end do

        case(7)
            do i = 1, nspe
                pro = inter_RPTBDB(qv, i)
                sif(i) = pro % vars(3) * 1.0e3
            end do

        case(10) ! incompressible
            p = abs(qv % v(ico) + g_pbase)
            t = abs(qv % v(ien))
            sif = htif(qv) * log(t) !- p / rhof(qv) / t
        case(100)
            sif = zero
        end select

    end function sif

    function pfromrhot(rho, t, qv)result(p)
        type(vector), intent(in) :: qv
        real(rfp) :: p, t, rho, mwinv
        ! only for perfect gas
        !
        if (s_irealgas /= 0) return
        !
        mwinv = sum(yif(qv)/g_mwi)

        p = rho * g_runi * t * mwinv

    end function pfromrhot

    function tfromrhop(rho, p, qv)result(t)
        type(vector), intent(in) :: qv
        real(rfp) :: p, t, rho, mwinv
        ! only for perfect gas
        !
        if (s_irealgas /= 0) return
        !
        mwinv = sum(yif(qv)/g_mwi)

        t = p / (rho * g_runi * mwinv)

    end function tfromrhop


    function rhoif(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: rhoif(nspe), p, t, tv, t_inv, tln
        integer :: i, j, n
        type(Qvar) :: pro

        select case(s_irealgas)
        case(0) ! perfect gas
            p = pf(qv)
            t = qv % v(ien)
            !    rhoif   = g_zi *  p * g_mwi / (g_runi * T) + g_rho0
            do i = 1, nspe
                select case(g_itype(i))
                case(0)
                    rhoif(i) = g_zi(i) * p * g_mwi(i) / (g_runi * T) + g_rho0(i)
                case(10) ! user defined
                    rhoif(i) = g_rho0(i)
                end select
            end do
        case(2, 3) ! NASA property polynormial
            p = pf(qv)
            t = qv % v(ien)
            do i = 1, nspe
                rhoif(i) = p * species(i) % mw / (g_runi * T)
            end do
        case(1)
            do i = 1, nspe
                !    rhoif(i) = inter_prop(fluid(i)%rho,qv,i)
                rhoif(i) = rhointer_prop(fluid(i) % rho, qv, i)
            end do

            !   do i = 1, nspe
            !    if(rhoif(i) <= zero) then
            !     print *, 'density become negative',i,rhoif(i),qv%v(ico),qv%v(ien)
            !     rhoif(i) = mytiny  !one 
            !    end if
            !   end do
        case(7)
            do i = 1, nspe
                pro = inter_RPTBDB(qv, i)
                rhoif(i) = pro % vars(1)
            enddo

        case(10, 100)
            rhoif = g_rho0
        case(110) ! metal model
            t = liquid_fraction(qv)
            tv = t/metal % rhol/((one - t)/metal % rhos + t/metal % rhol)
            rhoif = tv * metal % rhol + (one - tv) * metal % rhos
            !    rhoif = metal%rhol * t + (one - t) * metal%rhos 
        end select

    end function rhoif


    function rhopif(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: rhopif(nspe), t, p
        integer :: i
        type(Qvar) :: pro
        select case(s_irealgas)
        case(0) ! perfect gas
            p = qv % v(ico) + g_pbase
            t = qv % v(ien)
            do i = 1, nspe
                select case(g_itype(i))
                case(0)
                    rhopif(i) = g_zi(i) * g_mwi(i) / (g_runi * T)
                case(10) ! user defined
                    rhopif(i) = 1.e - 10_rfp ! given a small number
                end select
            end do
        case(2, 3)
            p = qv % v(ico) + g_pbase
            t = qv % v(ien)
            do i = 1, nspe
                rhopif(i) = species(i) % mw / (g_runi * T)
            end do
        case(1)
            do i = 1, nspe
                !    rhopif(i) = inter_prop(fluid(i)%rhop,qv,i)
                rhopif(i) = tinv_inter_prop(fluid(i) % rhop, qv, i)
            end do

        case(7)
            do i = 1, nspe
                pro = inter_RPTBDB(qv, i)
                rhopif(i) = pro % vars(4) * 1.0e - 6
            end do

        case(10)
            rhopif = zero !1.e-30_rfp   ! incompressible
        case(100)
            rhopif = zero
        case(110)
            rhopif = zero
        end select

    end function rhopif

    function rhotif(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: rhotif(nspe), t, p
        integer :: i
        type(Qvar) :: pro

        select case(s_irealgas)
        case(0) ! perfect gas
            p = qv % v(ico) + g_pbase
            t = qv % v(ien)
            do i = 1, nspe
                select case(g_itype(i))
                case(0)
                    rhotif(i) = -g_zi(i) * p * g_mwi(i) / (g_runi * T * T)
                case(10) ! user defined
                    rhotif(i) = zero
                end select
            end do
        case(2, 3)
            p = qv % v(ico) + g_pbase
            t = qv % v(ien)
            do i = 1, nspe
                rhotif(i) = -p * species(i) % mw / (g_runi * T * T)
            end do
        case(1)
            do i = 1, nspe
                !    rhotif(i) = inter_prop(fluid(i)%rhot,qv,i)
                rhotif(i) = tinv_inter_prop(fluid(i) % rhot, qv, i)
            end do
        case(7)
            do i = 1, nspe
                pro = inter_RPTBDB(qv, i)
                rhotif(i) = pro % vars(5)
            end do

        case(10, 100)
            rhotif = zero
        case(110)
            rhotif = liquid_fraction_t(qv) * (metal % rhol - metal % rhos)
        end select

    end function rhotif

    function hif(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: hif(nspe), t, t_inv, tln, cp, a1, a2, tr
        integer :: i, j, n
        type(Qvar) :: pro

        select case(s_irealgas)
        case(0) ! perfect gas
            !   t    = qv%v(ien) - g_htref
            do i = 1, nspe
                select case(g_itype(i))
                case(0)
                    hif(i) = g_hrefi(i) + g_cpi(i) * (qv % v(ien) - g_htref(i)) ! linear function of hi
                case(10) ! sublimation mixture
                    a1 = -11.391604_rfp; a2 = -0.39513431; tr = 83.8058
                    t = min(one, qv % v(ien) / tr)
                    hif(i) = g_hrefi(i) + g_cpi(i) * (qv % v(ien) - g_htref(i)) + &
                    a1 * g_runi / g_mwi(i) * tr + &
                    a2 * g_runi/g_mwi(i) * tr * ((one - t)**2.7_rfp + 2.7_rfp * t * (one - t)**1.7_rfp)
                end select
            end do
        case(1)
            do i = 1, nspe
                hif(i) = inter_prop(fluid(i) % h, qv, i)
            end do
        case(2, 3)
            t = qv % v(ien)
            do i = 1, nspe
                n = species(i) % interv
                t = min(max(species(i) % Titv(1), t), species(i) % Titv(n + 1))
                t_inv = one / t
                tln = log(t)
                do j = 1, n
                    if (t <= species(i) % Titv(j + 1)) then
                        hif(i) = polynomial(t, species(i) % hhp(:7, j)) * t_inv
                        hif(i) = species(i) % R * (hif(i) + species(i) % hhp(8, j) * tln)
                        exit
                    end if
                end do
                !
                ! extrapolation
                if (qv % v(ien) /= t) then
                    cp = species(i) % R * polynomial(t, species(i) % cpp(:, j)) * t_inv / t
                    hif(i) = hif(i) + cp * (qv % v(ien) - t)
                end if
            end do

        case(7)
            do i = 1, nspe
                pro = inter_RPTBDB(qv, i)
                hif(i) = pro % vars(2) * 1.0e3
            end do

        case(10) ! imcompressible
            hif = g_hrefi + g_cpi * (qv % v(ien) - g_htref) + pf(qv) / rhoif(qv) ! incompressible
        case(100)
            hif = g_hrefi + g_cpi * (qv % v(ien) - g_htref) ! solid
        case(110)
            a1 = liquid_fraction(qv)
            cp = metal % cpl * a1 + metal % cps * (one - a1)
            hif = ((metal % cps - metal % cpl) * metal % ts + metal % lheat) * a1 + cp * qv % v(ien)
            !   hif = metal%lheat * a1 + cp * qv%v(ien)
        end select

    end function hif

    function hpif(qv)
        type(vector) :: qv
        real(rfp) :: hpif(nspe)
        integer :: i
        type(Qvar) :: pro

        select case(s_irealgas)
        case(0, 2, 3) ! perfect gas
            hpif = zero
        case(1)
            do i = 1, nspe
                hpif(i) = inter_prop(fluid(i) % hp, qv, i)
            end do
        case(7)
            do i = 1, nspe
                pro = inter_RPTBDB(qv, i)
                hpif(i) = pro % vars(6) * 1.0e - 3
            end do

        case(10)
            hpif = one / rhoif(qv)
        case(100)
            hpif = zero
        case(110)
            hpif = zero
        end select

    end function hpif

    function htif(qv)
        type(vector) :: qv
        real(rfp) :: htif(nspe), t, t2_inv, a1, a2, tr
        integer :: i, j, n
        type(Qvar) :: pro

        !
        select case(s_irealgas)
        case(0, 10, 100) ! perfect gas
            do i = 1, nspe
                select case(g_itype(i))
                case(0)
                    htif(i) = g_cpi(i)
                case(10) ! sublimation mixture
                    a1 = -11.391604_rfp; a2 = -0.39513431; tr = 83.8058
                    t = min(one, qv % v(ien) / tr)
                    htif(i) = g_cpi(i) - &
                    a2 * g_runi / g_mwi(i) * 4.59_rfp * t * (one - t)**0.7_rfp
                end select
            end do
        case(110)
            a1 = liquid_fraction(qv)
            a2 = liquid_fraction_t(qv)
            tr = metal % cpl * a1 + metal % cps * (one - a1)
            htif = ((metal % cps - metal % cpl) * metal % ts + metal % lheat) * a2 + tr + (metal % cpl - metal % cps) * qv % v(ien) * a2
            !   htif = metal%lheat * a2 + tr  + (metal%cpl - metal%cps) * qv%v(ien) * a2
        case(1)
            do i = 1, nspe
                htif(i) = inter_prop(fluid(i) % ht, qv, i)
            end do
        case(2, 3)
            t = qv % v(ien)
            do i = 1, nspe
                n = species(i) % interv
                t = min(max(species(i) % Titv(1), t), species(i) % Titv(n + 1))
                t2_inv = one / (t * t)
                do j = 1, n
                    if (t <= species(i) % Titv(j + 1)) then
                        htif(i) = species(i) % R * polynomial(t, species(i) % cpp(:, j)) * t2_inv
                        exit
                    end if
                end do
            end do
        case(7)
            do i = 1, nspe
                pro = inter_RPTBDB(qv, i)
                htif(i) = pro % vars(7) * 1.0e3
            end do

        end select
    end function htif

    function pf(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: pf
        pf = qv % v(ico) + g_pbase
        if (pf <= zero) then
            !   print *,'negative pressure',qv%v
            pf = 1.e - 10_rfp
        end if
        ! if(s_ivis == 2 ) pf = pf - two_third * qv%v(ike) * rhof(qv)
    end function pf

    function h0f(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: h0f, hi(nspe), yi(nspe)
        h0f = hf(qv)
        if (nmeq > 0) h0f = h0f + half * dot_product(qv % v(imb:ime), qv % v(imb:ime))
        !  if(s_ivis == 2) &
        !  h0f = h0f + qv%v(ike)  
    end function h0f

    subroutine qvfromhs(h, s, qv)
        type(vector) :: qv
        real(rfp) :: h, s, p, t, hnew, snew, yi(nspe), cpmix, rmix
        real(rfp) :: sp, st, hp, ht, rho, a(2, 2), b(2), e, dp, dt, dh, ds
        real(rfp) :: alpha1 = 1.e - 4_rfp, alpha2 = 1.e - 2_rfp
        integer :: i, n = 10
        !
        select case(s_irealgas)
        case(0) ! perfect gas
            cpmix = htf(qv)
            t = (h - hf(qv)) / cpmix + qv % v(ien) ! cp is constant
            rmix = cpmix - cvf(qv)
            p = exp((cpmix * log(t/g_troom) - s)/rmix) * g_patm
            qv % v(ico) = p - g_pbase
            qv % v(ien) = t
        case default
            !
            !Use newton iteration method to solve ds=dh=0
            do i = 1, n
                hnew = hf(qv)
                !   snew = entropyf(qv)
                !   ds = s - snew 
                dh = h - hnew
                e = abs(dh/h)
                if (e < 1.e - 6_rfp) exit
                rho = rhof(qv)
                t = qv % v(ien)
                hp = hpf(qv)
                ht = htf(qv)
                dp = rho * dh
                dt = (dh - hp * dp) / ht
                !   dp = sign(one,dp) * min(abs(dp),alpha * pf(qv))
                !   dt = sign(one,dt) * min(abs(dt),alpha * qv%v(ien))
                qv % v(ico) = qv % v(ico) + dp
                qv % v(ien) = qv % v(ien) + dt
            end do
        end select

    end subroutine qvfromhs

    function mesh_vel(r)result(om)
        real(rfp) :: r(ndim), om(ndim), om1(3)

        om1 = zero; om = s_body_vel
        if (abs(s_omega) < mytiny) return
        om1(s_iaxis) = s_omega
        om1 = acrossb(om1, avto3(r))
        om = om + om1(:ndim)
    end function mesh_vel

    function stagnation_pressure(qv, r)result(p0)
        type(vector) :: qv, qvs
        real(rfp) :: p0, s, h0, mach, r(ndim), omega(3)
        !    mach = machf(qv)
        qvs = qv
        h0 = h0f(qvs)
        s = entropyf(qvs)
        call qvfromhs(h0, s, qvs)
        p0 = qvs % v(ico) + g_pbase
    end function stagnation_pressure

    function stagnation_temperature(qv, r)result(t0)
        type(vector) :: qv, qvs
        real(rfp) :: t0, s, h0, mach, r(ndim), omega(3)
        !    mach = machf(qv)
        qvs = qv
        s = entropyf(qvs)
        h0 = h0f(qvs)
        call qvfromhs(h0, s, qvs)
        t0 = qvs % v(ien)
    end function stagnation_temperature

    function pitot_pressure(qv)result(pt)
        type(vector) :: qv, qvs
        real(rfp) :: rhoold, uold, told, h0, s, mach, pt
        real(rfp) :: rhop, rhot, rho, u, u2, p, t, a(3, 3), b(3)
        real(rfp) :: gamma, gamma1, gOgamma1
        integer :: i, n = 20
        mach = machf(qv)
        !
        p = pf(qv)
        gamma = htf(qv) / cvf(qv)
        gamma1 = one / (gamma - one)
        gOgamma1 = gamma * gamma1
        if (mach < one) then
            pt = p * (one + (gamma - one) * half * mach**2)**gOgamma1
        else
            pt = p * (half * (gamma + one) * mach**2)**gOgamma1
            pt = pt /(two * gamma/(gamma + one) * mach**2 - (gamma - one)/(gamma + one))**gamma1
        end if
        return
        !
        if (mach <= one) then !subsonic
            s = entropyf(qv)
            h0 = h0f(qv)
            qvs = qv
            call qvfromhs(h0, s, qvs)
        else
            !Use newton iteration method to calculate normal shock
            qvs % v(ico) = pf(qv) * 2._rfp - g_pbase
            qvs % v(imb) = qv % v(imb) * 0.5_rfp
            qvs % v(ien) = qv % v(ien) * 0.25_rfp
            rhoold = rhof(qv)
            uold = qv % v(imb)
            told = qv % v(ien)
            do i = 1, n
                rho = rhof(qvs)
                p = qvs % v(ico)
                u = qvs % v(imb); u2 = u * u
                t = qvs % v(ien)
                b(1) = rho * u - rhoold * uold
                b(2) = rho * u2 + p - rhoold * uold**2 - qv % v(ico)
                b(3) = hf(qvs) + half * u2 - hf(qv) - half * uold**2
                rhop = rhopf(qvs); rhot = rhotf(qvs)
                a(1, 1) = rhop * u
                a(1, 2) = rho
                a(1, 3) = rhot * u
                a(2, 1) = rhop * u2 + one
                a(2, 2) = two * rho * u
                a(2, 3) = rhot * u2
                a(3, 1) = hpf(qvs)
                a(3, 2) = u
                a(3, 3) = htf(qvs)
                call direct_inv(a)
                b = matmul(a, -b)
                qvs % v(ico) = qvs % v(ico) + b(1)
                qvs % v(imb) = qvs % v(imb) + b(2)
                qvs % v(ien) = qvs % v(ien) + b(3)
                if (sum(abs(b)) < 1.e - 10_rfp) exit
            end do
            mach = machf(qvs)
            if (mach > one) print *, 'soluton failed', mach
            s = entropyf(qvs)
            h0 = hf(qvs) + half * qvs % v(imb)**2
            call qvfromhs(h0, s, qvs)
        end if
        pt = qvs % v(ico) + g_pbase
    end function pitot_pressure

    subroutine qvfromh0(rho, h0, qv, qvl, qvr)
        type(vector) :: qv, qvl, qvr
        real(rfp) :: rho, h0, rhoh, hh, h, dh, drho, d
        real(rfp) :: rhop, rhot, hp, ht, pmin, pmax, tmin, tmax
        real(rfp) :: d1, d2, d0, p1, p2, p0, cp, t, p
        integer :: i, n = 1
        !  h = h0 - half * sum(qv%v(imb:ime) **2)
        !  cp = (htf(qvl) + htf(qvr)) * half
        !  t = (h - hf(qvl) + cp * qvl%v(ien)) / cp
        !  p = t * sqrt((pf(qvl)/qvl%v(ien)) * (pf(qvr) / qvr%v(ien)))
        !  qv%v(ien) = t
        !  qv%v(ico) = p - g_pbase
        !  return
        !
        ! if(s_irealgas ==10.or.s_irealgas==100) return
        !
        pmin = min(qvl % v(ico), qvr % v(ico))
        pmax = max(qvl % v(ico), qvr % v(ico))
        tmin = min(qvl % v(ien), qvr % v(ien))
        tmax = min(qvl % v(ien), qvr % v(ien))

        h = h0 - half * sum(qv % v(imb:ime) **2)

        dh = h - hf(qv)
        d1 = dh / htf(qv)
        qv % v(ien) = d1 + qv % v(ien)
        drho = rho - rhof(qv)
        d2 = drho / rhopf(qv)
        qv % v(ico) = d2 + qv % v(ico)
        return


        do i = 1, n
            if (abs(dh/h) < 1.e - 5_rfp) exit
            if (abs(drho/rho) < 1.e - 5_rfp) exit
            ht = htf(qv)
            hp = hpf(qv)
            rhop = rhopf(qv)
            rhot = rhotf(qv)
            d = rhop * ht - rhot * hp
            d1 = (ht * dh - rhot * drho) / d
            d2 = (-hp * dh + rhop * drho) / d
            qv % v(ien) = d1 + qv % v(ien)
            qv % v(ico) = d2 + qv % v(ico)
        end do
        if (abs(dh) > one .or. abs(d0) > 0.01) print *, h - hf(qv), rho - rhof(qv)
        return
        !
        !
        !  find root for pressure
        d1 = rhof(qvl) - rho
        d2 = rhof(qvr) - rho
        p1 = qvl % v(ico)
        p2 = qvr % v(ico)
        p0 = qv % v(ico)
        !
        do i = 1, n
            d0 = rhof(qv) - rho
            if (abs(d0 / rho) < 1.e - 8_rfp) exit
            if (d1 * d0 <= zero) then
                p2 = p0
                d2 = d0
            else
                p1 = p0
                d1 = d0
            end if
            p0 = half * (p1 + p2)
            qv % v(ico) = p0
        end do
        !
        !  find root for temperature
        !
        d1 = hf(qvl) - h
        d2 = hf(qvr) - h
        p1 = qvl % v(ien)
        p2 = qvr % v(ien)
        p0 = qv % v(ien)
        !
        do i = 1, n
            d0 = hf(qv) - h
            if (abs(d0 / h) < 1.e - 8_rfp) exit
            if (d1 * d0 <= zero) then
                p2 = p0
                d2 = d0
            else
                p1 = p0
                d1 = d0
            end if
            p0 = half * (p1 + p2)
            qv % v(ien) = p0
        end do
        !
        return



        qv % v(ico) = half * exp(log(pmin + g_pbase) + log(pmax + g_pbase)) - g_pbase
        qv % v(ien) = half * exp(log(tmin) + log(tmax))
        !

    end subroutine qvfromh0

    !
    ! Transport properties
    !
    function mumix(qv)result(mu)
        type(vector), intent(in) :: qv
        real(rfp) :: mu, rho, e, f, mij, mui(nspe)
        real(rfp) :: xcon(nspe), molmu(nspe), h
        integer :: i, j
        type(Qvar) :: pro

        select case(s_irealgas)
        case(0, 10)
            mui = mu_suth(qv)
            mu = sum(yif(qv) * mui)
        case(1)
            do i = 1, nspe
                mui(i) = inter_prop(fluid(i) % mu, qv, i)
            end do
            mu = sum(yif(qv) * mui)
        case(2, 3)
            !   if(g_zmu < zero) then
            !    mu = abs(g_zmu)
            !    return
            !   end if
            !
            !the viscosity of mixture is determined from Wilk(1950)'s mixing rule and
            !modifiel by Bird et al.(1960)
            !
            !   Mumix=SUM_i{[X_i]Mu_i/SUM_j([X_j]F_ij)}
            !   F_ij = 1/8^1/2 (1+(Mu_i/Mu_j)^1/2 (W_j/W_i)^1/4)^2 / (1+W_i/W_j)^1/2
            ! or
            !   Mumix=SUM_i{[X_i]/SUM_j([X_j]E_ij)}
            !   E_ij = 1/8^1/2 ((1/Mu_i)^1/2+(1/Mu_j)^1/2 (MW_j/MW_i)^1/4)^2 / (1+M_i/M_j)^1/2

            !  rho = rhof(qv)
            !  yi = yif(qv)
            !  xcon(:) = rho * yif(qv) * species(:)%mwinv
            xcon(:) = yif(qv) * species(:) % mwinv ! concentration / rho
            molmu(:) = sqrt(one / molecular_viscosity(qv))
            mu = zero
            do i = 1, nspe
                e = zero
                do j = 1, nspe
                    mij = species(i) % mw * species(j) % mwinv
                    f = molmu(i) + molmu(j) / mij**0.25_rfp
                    f = f * f / sqrt(one + mij)
                    e = e + xcon(j) * f
                end do
                mu = mu + xcon(i) / e
            end do
            mu = two_sqrt2 * mu
        case(7)
            do i = 1, nspe
                pro = inter_RPTBDB(qv, i)
                mui(i) = pro % vars(8) * 1.e - 6
            end do
            mu = sum(yif(qv) * mui)
        case(100)
            mu = zero
        case(110)
            mu = metal % mul * rhof(qv) / metal % rhol
        end select
        if (s_ilevelset > 0) then ! level set
            h = heavisidef(qv)
            mu = (one - h) * mu + h * dummy % mul ! 
            if (h < 0.99) then
                mu = metal % mul * rhof(qv) / metal % rhol
            else
                mu = dummy % mul
            end if
        end if
    end function mumix

    function mu_suth(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: mu_suth(nspe), t1, s1
        integer :: i
        do i = 1, nspe
            if (g_zmu(i) < zero) then
                mu_suth(i) = -g_zmu(i)
            else
                !sutherland's law
                t1 = qv % v(ien) / g_tref(i) !  T/T0
                s1 = g_sref(i) !  S/T0
                mu_suth(i) = g_zmu(i) * t1 * sqrt(t1) * (one + s1) / (t1 + s1)
            end if
        end do
    end function mu_suth

    function lamda(qv, cc)
        type(cell) :: cc
        type(vector), intent(in) :: qv
        real(rfp) :: lamda, a, av, yi(nspe), lamdai(nspe)
        real(rfp) :: mole_f(nspe), mollamda(nspe), h
        integer :: i
        type(Qvar) :: pro
        if (cc % itype == partition_face_no) then ! special for interface cell
            lamda = cc % rp
            return
        end if

        select case(s_irealgas)
        case(0, 10)
            lamdai = mu_suth(qv) * htif(qv) / g_pr
            lamda = sum(yif(qv) * lamdai)
        case(1)
            do i = 1, nspe
                lamdai(i) = inter_prop(fluid(i) % lamd, qv, i)
            end do
            lamda = sum(yif(qv) * lamdai)
        case(2, 3)
            !
            !the mixture-averaged thermal conductivity is determined by Mathur et al. (1967)!
            !   Lamda=1/2 (SUM_i{X_iLamda_i} + 1/SUM_i{X_i/Lamda_i})
            yi = yif(qv)
            a = one / sum(yi(:) * species(:) % mwinv)
            mole_f(:) = yi(:) * species(:) % mwinv * a
            mollamda(:) = molecular_conductivity(qv)
            lamda = half * (sum(mole_f * mollamda) + one / sum(mole_f / mollamda))
            !
        case(7)
            do i = 1, nspe
                pro = inter_RPTBDB(qv, i)
                lamdai(i) = pro % vars(9) * 1.e - 3
            end do
            lamda = sum(yif(qv) * lamdai)
        case(100)
            lamda = g_pr(1)
        case(110)
            a = liquid_fraction(qv)
            av = a/metal % rhol/((one - a)/metal % rhos + a/metal % rhol)
            lamda = one / (av / metal % kl + (one - av)/metal % ks)
            !    lamda = one / (a / metal%kl + (one - a)/metal%ks)
        end select

        if (s_ilevelset > 0) then ! level set
            h = heavisidef(qv)
            lamda = (one - h) * lamda + h * dummy % kl ! 
        end if

    end function lamda

    function molecular_viscosity(qv)result(mu)
        type(vector), intent(in) :: qv
        real(rfp) :: mu(nspe), t, t_inv, tln, a(4)
        integer :: i, j, n
        t = qv % v(ien)
        ! 
        do i = 1, nspe
            n = species(i) % iv
            t = min(max(species(i) % vt(1), t), species(i) % vt(n + 1))
            t_inv = one / t
            tln = log(t)
            do j = 1, n
                if (t <= species(i) % vt(j + 1)) then
                    a = species(i) % mu(:, j)
                    mu(i) = exp(a(1) * tln + t_inv * (a(2) + a(3) * t_inv) + a(4)) * 1.e - 7_rfp
                    exit
                end if
            end do
        end do
        ! 
    end function molecular_viscosity

    function mole_fraction(qv)result(x)
        type(vector), intent(in) :: qv
        real(rfp) :: x(nspe), yi(nspe), a
        yi = yif(qv)
        a = one / sum(yi(:) * species(:) % mwinv)
        x(:) = yi(:) * species(:) % mwinv * a
    end function mole_fraction

    function molecular_conductivity(qv)result(lamda)
        type(vector), intent(in) :: qv
        real(rfp) :: lamda(nspe), t, t_inv, tln, a(4)
        integer :: i, j, n
        t = qv % v(ien)
        ! 
        do i = 1, nspe
            n = species(i) % ic
            t = min(max(species(i) % ct(1), t), species(i) % ct(n + 1))
            t_inv = one / t
            tln = log(t)
            do j = 1, n
                if (t <= species(i) % ct(j + 1)) then
                    a = species(i) % lamda(:, j)
                    lamda(i) = exp(a(1) * tln + t_inv * (a(2) + a(3) * t_inv) + a(4)) * 1.e - 4_rfp
                    exit
                end if
            end do
        end do
        ! 
    end function molecular_conductivity

    function mass_diffusivity(qv)result(d)
        type(vector), intent(in) :: qv
        real(rfp) :: d(nspe)
        real(rfp) :: yi(nspe), mole_fraction(nspe), mu
        real(rfp) :: a_mole, rho, t, p, length, epsT, epsT2, omegad, dij, sumd
        integer :: i, j
        !  if(s_irec == 0) then
        !   d = zero
        !   return
        !  end if
        select case(s_irealgas)
        case(0, 1, 10)
            mu = mumix(qv)
            d = mu / g_sch
            !   d = zero
        case(7)
            mu = mumix(qv)
            d = mu / g_sch
        case(2, 3)
            !            0.0266T^3/2
            !   D_kl = ------------------------------
            !            pW_kl^1/2Epsi_kl^2Omega_D
            !
            yi = yif(qv)
            t = qv % v(ien)
            p = pf(qv)
            rho = rhof(qv)
            a_mole = one / sum(yi(:) * species(:) % mwinv)
            mole_fraction(:) = yi(:) * species(:) % mwinv * a_mole + 1.e - 12_rfp
            ! 
            do i = 1, nspe
                !   t = min(max(species(i)%tmin,t),species(i)%tmax)
                sumd = zero !mytiny
                do j = 1, nspe
                    if (i == j) cycle
                    ! Lennard-Jones characteristic length
                    length = half * (species(i) % cl + species(j) % cl)
                    ! the reduced temperature
                    epsT = t / sqrt(species(i) % ce * species(j) % ce)
                    epsT2 = one / (half + epsT)
                    !  
                    omegad = epsT**(-0.145_rfp) + epsT2 * epsT2
                    !  
                    a_mole = two / (species(i) % mwinv + species(j) % mwinv)
                    !   dij = (0.0266_rfp * t * sqrt(t/a_mole))/(p * length * length * omegad)
                    dij = (0.0266_rfp * t**(3/2))/(p * sqrt(a_mole) * length * length * omegad)
                    sumd = sumd + mole_fraction(j) / dij
                end do
                if (sumd < mytiny) then
                    d(i) = zero
                else
                    d(i) = rho * abs(one - mole_fraction(i)) / sumd
                end if
            end do
        case(110)
            d = metal % dl
        end select
    end function mass_diffusivity

    function absorption_coef(qv)result(alpha)
        type(vector) :: qv
        ! type(cell),optional::cc
        real(rfp) :: alpha
        select case(s_irealgas)
        case(0)
            alpha = g_absorptivity
            !   alpha = 100. * exp(-cc%centp(2) / b_lenref /0.1)
        case(1)
            alpha = inter_prop(fluid(1) % alpha, qv, 1)
        end select
        !
    end function absorption_coef

    function mu_turb(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: mu_turb
        !  mu_turb = min(rhof(qv) * abs(qv%v(ike) / qv%v(iom)),0.001_rfp)
        mu_turb = rhof(qv) * abs(qv % v(ike) / qv % v(iom))
    end function mu_turb


    function sound_speed(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: rhop, rhot, hp, ht, rho, dpp
        real(rfp) :: sound_speed
        rho = rhof(qv)
        rhop = rhopf(qv)
        rhot = rhotf(qv)
        ht = htf(qv)
        hp = hpf(qv)
        dpp = rhop + rhot * (one - rho * hp) / (rho * ht)
        if (dpp <= mytiny) then
            sound_speed = 1.e8_rfp
        else
            sound_speed = abs(one / dpp)
        end if

    end function sound_speed

    function machf(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: vel(ndim), machf
        vel = qv % v(imb:ime)
        machf = sqrt(sum(vel * vel)/sound_speed(qv))
    end function machf

    function yif(qv)
        type(vector), intent(in) :: qv
        real(rfp) :: yif(nspe)
        if (nspe == 1) then
            yif(1) = one
            return
        end if
        yif(:nspm1) = min(max(zero, qv % v(isb:ise)), one)
        !  where(yif < 1.e-20_rfp) yif = zero
        yif(nspe) = max(zero, one - sum(yif(:nspm1)))
    end function yif

    function volume_fraction(qv, ispec)result(vf)
        type(vector), intent(in) :: qv
        real(rfp) :: vf, yi(nspe), rhoi(nspe)
        integer :: ispec
        yi = yif(qv)
        rhoi = rhoif(qv)
        vf = rhof(qv) * yi(ispec) / rhoi(ispec)
    end function volume_fraction

    function concentration(qv, ispec)result(vf)
        type(vector), intent(in) :: qv
        real(rfp) :: vf, yi(nspe)
        integer :: ispec
        yi = yif(qv)
        vf = rhof(qv) * yi(ispec) / g_mwi(ispec)
    end function concentration

    !function volume_fraction(qv)result(vf)
    !  type(vector),intent(in)::qv
    !  real(rfp)::vf,y(nspe),r(nspe)
    !   y = yif(qv)
    !   r = rhoif(qv)
    !  vf = rhof(qv) * y(1) / r(1)
    !end function volume_fraction

    function tinv_inter_prop(f, qv, isp)result(fi)
        type(vector) :: qv
        real(rfp), intent(in) :: f(:,:,:)
        real(rfp) :: fp, ft, fc, p, t, fi, t1, t2
        integer, intent(in) :: isp
        integer :: i, j, k, np, nt, nc
        !  p = max(log(pf(qv)),fluid(isp)%p0)
        p = log(pf(qv))
        t = qv % v(ien)
        i = location(p, fluid(isp) % p0, fluid(isp) % dp, fluid(isp) % np, fp)
        j = location(t, fluid(isp) % t0, fluid(isp) % dt, fluid(isp) % nt, ft)
        ! only mixture allowed
        k = 1
        fc = 0.0
        t = one / t
        t1 = one / (fluid(isp) % t0 + real(j - 1, rfp) * fluid(isp) % dt)
        t2 = one / (fluid(isp) % t0 + real(j, rfp) * fluid(isp) % dt)
        ft = (t - t1) / (t2 - t1)
        !
        !
        fi = f(i, j, k) * (one - fp) * (one - ft) * (one - fc) + &
        f(i + 1, j, k) * fp * (one - ft) * (one - fc) + &
        f(i, j + 1, k) * (one - fp) * ft * (one - fc) + &
        !       f(i  ,j  ,k+1) * ( one - fp) * ( one - ft) *         fc + &
        f(i + 1, j + 1, k) * fp * (ft) * (one - fc) !+ &
        !       f(i  ,j+1,k+1) * ( one - fp) * (       ft) * (       fc) + &
        !       f(i+1,j  ,k+1) * (       fp) * ( one - ft) * (       fc) + &
        !       f(i+1,j+1,k+1) * (       fp) * (       ft) * (       fc) 
        !
    end function tinv_inter_prop

    function saturationp(qv)result(p)
        type(vector) :: qv
        real(rfp) :: t, p, ft
        integer :: j
        t = qv % v(ien)
        j = location(t, fluid(1) % t0, fluid(1) % dt, fluid(1) % nt, ft)
        p = fluid(1) % satp(j) * (one - ft) + fluid(1) % satp(j + 1) * ft
    end function saturationp

    function inter_prop(f, qv, isp)result(fi)
        type(vector) :: qv
        real(rfp), intent(in) :: f(:,:,:)
        real(rfp) :: fp, ft, fc, p, t, fi
        integer, intent(in) :: isp
        integer :: i, j, k, np, nt, nc
        !  p = max(log(pf(qv)),fluid(isp)%p0)
        p = log(pf(qv))
        t = qv % v(ien)
        i = location(p, fluid(isp) % p0, fluid(isp) % dp, fluid(isp) % np, fp)
        j = location(t, fluid(isp) % t0, fluid(isp) % dt, fluid(isp) % nt, ft)
        ! only mixture allowed
        k = 1
        fc = 0.0
        !
        fi = f(i, j, k) * (one - fp) * (one - ft) * (one - fc) + &
        f(i + 1, j, k) * fp * (one - ft) * (one - fc) + &
        f(i, j + 1, k) * (one - fp) * ft * (one - fc) + &
        !       f(i  ,j  ,k+1) * ( one - fp) * ( one - ft) *         fc + &
        f(i + 1, j + 1, k) * fp * (ft) * (one - fc) !+ &
        !       f(i  ,j+1,k+1) * ( one - fp) * (       ft) * (       fc) + &
        !       f(i+1,j  ,k+1) * (       fp) * ( one - ft) * (       fc) + &
        !       f(i+1,j+1,k+1) * (       fp) * (       ft) * (       fc) 
        !
    end function inter_prop

    function rhointer_prop(f, qv, isp)result(fi)
        type(vector) :: qv
        real(rfp), intent(in) :: f(:,:,:)
        real(rfp) :: fp, ft, fc, p, t, fi, t1, t2
        integer, intent(in) :: isp
        integer :: i, j, k, np, nt, nc
        !  p = max(log(pf(qv)),fluid(isp)%p0)
        p = log(pf(qv))
        t = qv % v(ien)
        i = location(p, fluid(isp) % p0, fluid(isp) % dp, fluid(isp) % np, fp)
        j = location(t, fluid(isp) % t0, fluid(isp) % dt, fluid(isp) % nt, ft)
        t = log(t)
        t1 = log(fluid(isp) % t0 + real(j - 1, rfp) * fluid(isp) % dt)
        t2 = log(fluid(isp) % t0 + real(j, rfp) * fluid(isp) % dt)
        ft = (t - t1) / (t2 - t1)
        ! only mixture allowed
        k = 1
        fc = 0.0
        !
        fi = log(f(i, j, k)) * (one - fp) * (one - ft) * (one - fc) + &
        log(f(i + 1, j, k)) * fp * (one - ft) * (one - fc) + &
        log(f(i, j + 1, k)) * (one - fp) * ft * (one - fc) + &
        !       f(i  ,j  ,k+1) * ( one - fp) * ( one - ft) *         fc + &
        log(f(i + 1, j + 1, k)) * fp * (ft) * (one - fc) !+ &
        !       f(i  ,j+1,k+1) * ( one - fp) * (       ft) * (       fc) + &
        !       f(i+1,j  ,k+1) * (       fp) * ( one - ft) * (       fc) + &
        !       f(i+1,j+1,k+1) * (       fp) * (       ft) * (       fc) 
        !
        fi = exp(fi)
        !
    end function rhointer_prop

    function location(x, x0, dx, n, fac)result(loc)
        real(rfp), intent(in) :: x0, dx, x
        real(rfp) :: fac
        integer :: loc, n
        fac = (x - x0) / dx + one
        loc = int(fac)
        fac = fac - loc
        if (loc < 1) then
            fac = loc + fac - one
            loc = 1
        else if (loc >= n) then
            fac = fac + loc - real(n + 1, rfp)
            loc = n - 1
        end if
    end function location

    function location1(x, x0, n)result(loc)
        real(rfp), intent(in) :: x0(:), x
        integer :: loc, n, i
        if (n /= size(x0)) then
            print *, 'error in x0 array', n, size(x0)
            return
        end if
        !
        loc = n - 1
        do i = 1, n
            if (x <= x0(i)) then
                loc = i - 1
                exit
            end if
        end do
        loc = max(1, loc)
    end function location1

    subroutine divf(cells, faces, nodes)
        !********************************************************
        implicit none
        type(node), pointer :: nodes(:)
        type(cell), pointer :: cells(:)
        type(face), pointer :: faces(:), fs
        ! local variables
        type(vector) :: qv
        real(rfp) :: u, ru, ra, vref, vecn(ndim)
        ! variables for computation
        integer :: i, k
        !
        cells % srf = zero
        cells % rp = zero

        do i = 1, size(faces)
            fs => faces(i)
            vecn = fs % vecn * fs % area
            qv % v = zero
            do k = 1, neq
                qv % v(k) = sum(nodes(fs % f2n) % qv % v(k))
            end do
            qv % v = qv % v / real(size(fs % f2n), rfp)
            u = sum(qv % v(imb:imb + ndim - 1) * vecn)
            ru = u * rhof(qv)
            fs % left_cell % srf = fs % left_cell % srf - u
            fs % left_cell % rp = fs % left_cell % rp - ru

            fs % right_cell % srf = fs % right_cell % srf + u
            fs % right_cell % rp = fs % right_cell % rp + ru
        end do

        cells % srf = cells % srf / cells % vol
        cells % rp = cells % rp / cells % vol
        qv % v(ien) = 600._rfp
        qv % v(ico) = zero
        ra = g_pr(1) * abs(g_gravity(2)) * rhof(qv)**2 * 720._rfp * b_lenref**3/600._rfp/mumix(qv)**2
        vref = mumix(qv) * sqrt(ra)/rhof(qv)/b_lenref
        print *, 'Ra=', ra, ' vref=', vref
        print *, 'divuMax=', maxval(cells % srf) * b_lenref/vref
        print *, 'divRuMax=', maxval(cells % rp) * b_lenref/vref
        print *, 'divuMin=', minval(cells % srf) * b_lenref/vref
        print *, 'divRuMin=', minval(cells % rp) * b_lenref/vref
        print *, 'divuav=', sum(abs(cells % srf))/real(size(cells), rfp) * b_lenref/vref
        print *, 'divRuav=', sum(abs(cells % rp))/real(size(cells), rfp) * b_lenref/vref
        print *, 'Pmax/p0=', (maxval(cells % qv % v(ico)) + g_pbase)/g_pbase
        print *, 'Pmin/p0=', (minval(cells % qv % v(ico)) + g_pbase)/g_pbase
    end subroutine divf

    function e_resistivity(qv, cc)result(e)
        type(cell) :: cc
        type(vector), intent(in) :: qv
        real(rfp) :: e, sigma, t, sig(nspe)
        integer :: it, i
        it = cc % itype
        if (it == partition_face_no.and.neqm > 0) then
            e = cc % er
            return
        end if
        select case(s_irealgas)
        case(0, 2, 3, 10) ! perfect gas
            !     sigma = s_vnn   !vc(it)%vars(1)
            !      sigma = sum(g_sigma*yif(qv))
            sigma = one / sum(yif(qv)/g_sigma)
            !    if(qv%v(isb) > 0.9_rfp) then
            t = max(min(qv % v(ien), 1.2e5_rfp), 1000._rfp)
            !      sigma = (t - 4.e3_rfp) / (1.2e5_rfp-4.e3_rfp) * (1.e4_rfp-200._rfp) + 200._rfp
            !      sigma = max(sigma, g_sigma(1))
            !     end if
            sigma = exp(12.648_rfp - 54550._rfp/t)
        case(1)
            do i = 1, nspe
                sig(i) = rhointer_prop(fluid(i) % sigma, qv, i)
            end do
            sigma = sum(yif(qv) * sig)
            sigma = max(1.e - 3_rfp, sigma)
        case(100)
            sigma = g_sigma(1)
            sigma = vc(it) % vars(1) + mytiny
            !    e = min(one,max(zero,cc%centp(1)))
            !    sigma = 1.e2*(exp(3.*e)-one) + 1.e-5
            !    e = min(one,max(zero,cc%centp(2)/0.1))
            !    sigma = 1.e3*(exp(3.0) - exp(3.*e)) + 1.e-5
            if (sigma < zero) sigma = one / cc % er
        end select
        e = one / sigma ! 1/(sigam)
    end function e_resistivity

    function sigma_pf(qv, cc)result(e)
        type(cell) :: cc
        type(vector), intent(in) :: qv
        real(rfp) :: e, sigma, t
        integer :: it
        it = cc % itype
        !    if(vc(it)%itype == 0) then
        e = zero
        !    end if
    end function sigma_pf

    function sigma_tf(qv, cc)result(e)
        type(cell) :: cc
        type(vector), intent(in) :: qv
        real(rfp) :: e, sigma, t
        integer :: it
        it = cc % itype
        e = zero
        return
        !
        if (vc(it) % itype == 0) then
            if (vc(it) % vars(1) < 0) then
                e = zero
            else
                t = max(qv % v(ien), 300._rfp)
                sigma = exp(12.645 - 54550./t)
                e = 54550.0 / t**2 * sigma
            end if
        end if
    end function sigma_tf

    function hall_resistivity(qv, cc)result(e)
        type(cell) :: cc
        type(vector), intent(in) :: qv
        real(rfp) :: e, ne, ei(nspe)
        integer :: it, i
        it = cc % itype
        select case(s_irealgas)
        case(0, 2, 3, 10) ! perfect gas
            ne = sum(yif(qv) * g_ne)
        case(1)
            do i = 1, nspe
                ei(i) = rhointer_prop(fluid(i) % ne, qv, i)
            end do
            ne = sum(yif(qv) * ei)
        case(100)
            ne = 1.e30_rfp !vc(it)%vars(2)
        end select
        !    ne = 1.e20_rfp
        e = one / (ne * g_elec) ! 1/(ne*ec)
    end function hall_resistivity


    function bthetaf(r1, r2, current, r) result(btheta)
        real(rfp), intent(in) :: r1, r2, current, r
        real(rfp) :: btheta
        real(rfp) :: ci
        ci = g_mu0 * current * half / pi
        if (r < r1) then
            btheta = zero
        else if (r < r2) then
            btheta = (r * r - r1 * r1)/(r2 * r2 - r1 * r1) * ci / r
        else
            btheta = ci / r
        end if
    end function bthetaf

    function bfun(qv, i)
        type(vector), intent(in) :: qv
        real(rfp) :: bfun(3)
        integer :: i
        select case(s_ibfield)
        case(1)
            bfun = qv % v(ibb:ibe) + g_b0
        case(2)
            bfun = qv % v(ibb:ibe)
            bfun(3) = bfun(3) + e_b0(i) ! external magnetic field in z-direction
        case default
            bfun = qv % v(ibb:ibe)
        end select
    end function bfun

    function inter_RPTBDB(qv, irgas)result(pro)
        type(Qvar) :: pro
        type(vector) :: qv
        real(rfp) :: p, t, epsi = 1.e - 10_rfp
        integer :: ierr, irgas
        type(p2FttRootCell) :: RBB

        ierr = 0
        if (abs(s_origP(irgas) - qv % v(ico)) < epsi.and.abs(s_origT(irgas) - qv % v(ien)) < epsi) then
            pro = s_pro(irgas)
            return ! do not need re-interpolate
        else
            s_origP(irgas) = qv % v(ico); s_origT(irgas) = qv % v(ien)
        end if
        p = max(min(1., log10(pf(qv)) - 6.0), 0.778)
        t = max(min(400., qv % v(ien)), 60.)
        !  p = log10(pf(qv)) - 6.0_rfp
        !  t = qv%v(ien)
        pro = QLocate(RB(irgas) % p, t, p, ierr)
        s_pro(irgas) = pro

    end function inter_RPTBDB

    function liquid_fraction(qv)result(g)
        real(rfp) :: g
        type(vector) :: qv
        ! linear liquid fraction model for steel
        g = max(zero, min(one, (qv % v(ien) - metal % ts) / (metal % tl - metal % ts)))
    end function liquid_fraction

    function liquid_fraction_t(qv)result(g)
        real(rfp) :: g, t
        type(vector) :: qv
        ! linear liquid fraction model for steel
        t = qv % v(ien)
        if (t < metal % ts.or.t >= metal % tl) then
            g = zero
        else
            g = one / (metal % tl - metal % ts)
        end if
    end function liquid_fraction_t

    function heavisidef(qv)result(g)
        ! Heaviside function
        real(rfp) :: g, fi, epsi
        type(vector) :: qv
        if (s_ilevelset == 0) return ! no level set
        epsi = dummy % k0 ! stored critiria
        fi = qv % v(ils) / epsi
        if (fi < -one) then
            g = zero
        else if (fi > one) then
            g = one
        else
            g = half * (one + fi + sin(pi * fi) / pi)
        end if
    end function heavisidef

    function DiracDeltaf(qv)result(g)
        ! Dirac Delta function
        real(rfp) :: g, fi, epsi
        type(vector) :: qv
        if (s_ilevelset == 0) return ! no level set
        epsi = one / dummy % k0 ! stored critiria
        fi = qv % v(ils) * epsi
        if (fi < -one) then
            g = zero
        else if (fi > one) then
            g = zero
        else
            g = half * epsi * (one + cos(pi * fi))
        end if
    end function diracdeltaf

end MODULE dfd_state
