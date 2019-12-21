MODULE dfd_data
    use dfd_type
    USE RPTBDB_setup, ONLY: p2FttRootCell, Qvar
    implicit none
    ! boundary value
    type(bc_type), pointer :: bc(:), vc(:) ! vc is volume conditions
    real(rfp) :: b_lenref, b_re, b_p0, b_pinf, b_rhoinf, b_vinf, b_nblade

    !   real(kind=rfp)::b_p0,b_t0,b_rhoinf,b_pback,b_pinf,b_velinf(ndim),b_tinf
    !   real(kind=rfp)::b_velref,b_lenref,b_re,b_u2,b_mf,b_wallvel(ndim),b_wt,b_t1,b_p1
    !   integer::b_p0t0
    ! Multi-physical zone
    type(phy), pointer :: pz(:)
    ! solution_control
    integer :: s_nstart, s_nend, s_nstep !--- min and max. number of iterations
    integer :: s_init !--- init == 0 => initial from constant
    !    init /= 0 => initial from out file
    integer :: s_ipre !--- ipre  == 0 => without preconditioning                                       !    ipre  /= 0 => with preconditioning
    integer :: s_ivis !--- ivis == 0 =>  inviscous flow(Euler eq)
    !    itime = 1 => viscous flow(N-S eq)
    !    itime = 2 => viscous flow(Turbulence)
    integer :: s_visb !  ==0 automatical BC for k,omega, ==1 specific
    integer :: s_irec ! =1  least square  2 = Ding's method
    integer :: s_ischeme ! =0 Roe's 1 HLL 2 HLLC 3 Godunov's
    ! =-1 Roe's with algebraic average
    integer :: s_irealgas ! indicator of real gas
    integer :: s_ialg !--- ialg == 0 first order
    !----     == 1 or -1 second order negative number with limiter
    !----     == 2 or -2 third order negative number with limiter
    integer :: s_imeth !--- imeth == 0 explicit
    !----      == 1 block point gauess seideil
    !----      == 2 block DILU
    !----      == 3 block GMRES with D-ILU Preconditioning
    integer :: s_nrec
    integer :: s_iscr !--- iscr == 0 => no screen output
    !    iscr == 1 => write convergence to screen
    integer :: s_ifile !--- ifile == 0 => no file output
    !    ifile == 1 => write convergence to file
    integer :: s_iplot !--- 0 no plot 1 tecplot 2 fieldview
    integer :: s_nplot !--- each npplot iterations to plot
    integer :: s_isub !--- itrations of sub SOR
    integer :: s_ifmt !--- grid type
    integer :: s_iaxis !--- axis symmetry = 1 r->x r -> y
    integer :: s_iperiod !--- periodic boundary conditions
    integer :: s_istop, s_ni, s_nj, s_nk
    integer :: s_isource ! indicator for volume condition > 0 has volume source
    integer :: s_ilevelset ! indicator for volume condition > 0 has volume source
    integer :: s_iaux ! indicator for auxiliary equations
    integer :: s_ibfield ! indicator for external magnetic field
    real(kind = rfp) :: s_cfl, s_cflmin, s_cflmax, s_vnn, s_damping !--- cfl number
    real(kind = rfp) :: s_cflm, s_cflmm
    real(kind = rfp) :: s_errm, s_alpha = zero !--- max. allowable error
    real(kind = rfp) :: s_cspeed = 2.99792458e8_rfp !--- light speed
    real(kind = rfp) :: s_bspeed !--- vitural magnitic speed
    real(kind = rfp) :: s_espeed !--- vitural electric speed
    real(kind = rfp) :: s_fspeed !--- vitural diffusion speed
    real(kind = rfp) :: s_pcenter(ndim) !--- center of rotate axis
    real(kind = rfp) :: s_omega, s_omega0, s_damping_o !--- rotate speed (1/s)
    real(kind = rfp) :: s_body_vel(ndim) !--- rigid body moving velocity (m/s)
    ! 
    real(kind = rfp) :: s_origP(nspe) = zero, s_origT(nspe) = zero
    type(Qvar) :: s_pro(nspe)

    type(vector) :: s_err
    !
    !  Unsteady parameters
    !
    integer :: s_idt !--- itime == 0 => steady case
    !    itime /= 0 => unsteady case
    integer :: s_nt !--- number of iterations per DT cycle	
    integer :: s_imonit, npos !--- output redults in every time step
    integer, pointer :: s_ipos(:) ! max monitoring points
    real(kind = rfp) :: s_dt_inv !--- inversion of real time step(deltT)
    real(kind = rfp) :: s_elapsed_time
    real(rfp), pointer :: gear_c(:) !Coeffient for Gear Methods
    integer, pointer :: gear_in(:) !previous value location in qvn
    !     
    !** Perfect gas property data structure
    !
    real(kind = rfp) :: g_runi = 8314.9 !--- universal gas constant
    real(kind = rfp) :: g_patm = 1.013e5_rfp !--- 1atm = 101.3Kpa
    real(kind = rfp) :: g_troom = 278.0 ! Kelvin
    !
    real(kind = rfp) :: g_mwi(nspe) !--- molecular weight
    real(kind = rfp) :: g_zi(nspe) !--- molecular weight
    real(kind = rfp) :: g_rho0(nspe) !--- molecular weight
    real(kind = rfp) :: g_r !--- gas constant
    real(kind = rfp) :: g_gm !--- the ratio of specific heats
    real(kind = rfp) :: g_cpi(nspe) !--- the specific heat at constant pressure
    real(kind = rfp) :: g_cvi(nspe) !--- the specific heat at constant volume
    real(kind = rfp) :: g_hrefi(nspe), g_htref(nspe) !--- ref h and ref t

    real(kind = rfp) :: g_sigma(nspe) !--- the kinematic viscosity
    real(kind = rfp) :: g_ne(nspe) !--- the kinematic viscosity
    real(kind = rfp) :: g_zmu(nspe) !--- the kinematic viscosity
    real(kind = rfp) :: g_tref(nspe), g_sref(nspe) !--- used for Sutherland's law
    real(kind = rfp) :: g_pr(nspe) !--- the Prandtl number(0.72)
    real(kind = rfp) :: g_prt = 0.9_rfp !--- Turbulent the Prandtl number for turbulence
    real(kind = rfp) :: g_sch = 0.62_rfp !--- the Schmidt number
    real(kind = rfp) :: g_scht = 1._rfp !--- the Turbulent Schmidt number

    real(kind = rfp) :: g_pbase, g_rhobase !--- base pressure and density
    real(kind = rfp) :: g_gravity(ndim) !--- gravity
    real(kind = rfp) :: g_absorptivity !--- absorptance or absorptivity
    real(kind = rfp) :: g_cdest, g_cprod ! coefficients of rate of phase change
    integer :: g_igrav ! > 0 count gravity
    integer :: g_itype(nspe)
    !
    ! Reaction coefficients
    !
    integer :: nrea
    type(therm) :: species(nspe)
    type(react), pointer :: reaction(:)
    !
    ! real fluid
    type(gas) :: fluid(nspe)
    ! solid
    type(solid) :: metal
    type(laser) :: laser1
    type(powderflow) :: powderflow1

    !
    !Turbelence model(k-omega two equation model)
    !
    real(kind = rfp) :: t_alpha = 13.0_rfp / 25.0_rfp
    !  real(kind=rfp)::t_beta0  = 9.0_rfp / 125.0_rfp
    real(kind = rfp) :: t_beta0 = 0.0708_rfp ! Wilcox(2006)
    real(kind = rfp) :: t_betas0 = 9.0_rfp / 100.0_rfp
    real(kind = rfp) :: t_betas3 = 729.0_rfp / 1.0e6_rfp
    real(kind = rfp) :: t_psi = 0.5_rfp
    !  real(kind=rfp)::t_psis  = 0.5_rfp  
    real(kind = rfp) :: t_psis = 0.6_rfp ! Wilcox(2006)
    real(kind = rfp) :: t_psido = one/ 8._rfp ! Wilcox(2006)
    real(kind = rfp) :: t_keinf, t_omegainf, t_sr
    real(kind = rfp) :: t_cdes = -one !  > 0 indicate to DES method
    !
    !!Turbulence model(Baldwin-Lomax)
    ! 
    real(rfp), pointer :: zmut_bl(:)
    integer, pointer :: p2f_bl(:)
    type(wallbl), pointer :: wallbc(:)
    !
    real(rfp) :: stephan_boltzmann = 5.672e-8_rfp ! Wm^-2K^-4
    ! Output selection
    logical :: output_selection(100)
    character(len = 40) :: vname(100), wname(100)
    integer :: nwplot, ivplot
    integer, pointer :: wlabel(:)
    logical, pointer :: woutput_selection(:,:)
    !MHD
    real(kind = rfp) :: g_mu0 = 0.1256637e-5_rfp ! Permeability of free space (N/Amp^2)
    real(kind = rfp) :: g_epsi0 = 8.85418782e-12_rfp ! Permittivity of free space (S^4A^2/(m^-3kg^-1))
    real(kind = rfp) :: g_mu0inv = 7.9577472e5_rfp ! Permeability of vacuum (N/Amp^2
    real(kind = rfp) :: g_elec = 1.60217733e-19_rfp ! Electron charge  coul=Amp Second
    real(kind = rfp) :: g_b0(3) ! external b field
    real(kind = rfp), pointer :: e_b0(:) ! external b field

    ! schedule
    integer :: s_ischedule
    type(job) :: schedule
    !
    type(p2FttRootCell) :: RB(nspe)
    !
    ! swirl flow & interface
    integer :: ispecialface = 0
    !  integer::iswirl=0
    !  integer,pointer::facelist(:)
    !
    !  Maxwell solver
    integer :: m_input = 0
    ! Units
    real(rfp) :: psi2pa = 6894.76_rfp
    real(rfp) :: kg2lb = 2.20462262_rfp

    ! Wenda's single core treatment
    integer :: isize_ext, jsize_ext, ksize_ext
    real(rfp), allocatable :: xc_ext(:), yc_ext(:), zc_ext(:)
    real(rfp) :: dxc_ext, dyc_ext, dzc_ext

    double precision, allocatable :: var1_ext(:,:,:), var2_ext(:,:,:), var3_ext(:,:,:), var4_ext(:,:,:)
    integer, allocatable :: map_ext(:,:,:)

    !level-set properties
    real(kind = rfp) :: ls_thickness, material, intf_treatment, cfl_adv, cfl_RI, init, adv_freq, RI_freq, RI_quit, &
        adv_nRK, RI_nRK
    real(kind = rfp) :: ls_cfl_adv
    real(kind = rfp) :: ls_cfl_RI
    integer :: ls_material
    integer :: ls_intf_treatment
    integer :: ls_init
    integer :: ls_adv_freq
    integer :: ls_RI_quit
    integer :: ls_RI_freq
    integer :: ls_adv_nRK
    integer :: ls_RI_nRK
    real(kind = rfp) :: ls_xi

end MODULE dfd_data
