module gems_type
  use gems_constant
  !** define a vector data structure
  !
  type phy
    integer :: neq, imeth, ndim
    integer, pointer :: loc(:)
  end type phy

  type vector
    real(kind = rfp), dimension(neq) :: v
  end type vector
  !
  !** define a matrix data structure
  !
  type matrix
    real(kind = rfp), dimension(neq, neq) :: e
  end type matrix
  !
  type neighbour
    !   type(cell),pointer::to_cell
    type(face), pointer :: to_face
  end type neighbour

  type neighbour_cell
    type(cell), pointer :: to_cell => NULL()
  end type neighbour_cell

  type chain
    type(chain), pointer :: up_chain, next_chain
    type(cell), pointer :: cc
    type(face), pointer :: cf
    !   type(matrix)::cm
    !   type(vector)::res
  end type chain

  type raw_chain
    type(chain), pointer :: bchain(:)
    integer, pointer :: nchain(:), icell(:)
    integer :: npgs
  end type raw_chain

  type node
    real(kind = rfp) :: xyz(ndim) ! coordination
    !   real(kind=rfp)::vol         ! volume of surronding the node
    type(vector) :: qv ! the solution at each node
    !   type(vector)::gqv(ndim)     ! gradient qv

    !for particle tracking
    real(rfp), pointer :: precoil
  end type node

  type bnode
    integer :: n, itype
    type(node), pointer :: tn
    type(neighbour_cell), pointer :: scell(:)
  end type bnode

  type cell
    integer :: itype ! index of color
    type(vector) :: qv, dqv ! q and dq variables
    type(vector) :: res ! Residual of control volume
    type(matrix) :: dm !,lm,um          ! Diagonal matrix
    real(kind = rfp), pointer :: gradient(:,:) ! gradient for interpolation
    real(kind = rfp) :: srf, srfm, rp, er ! max spectral radius flux (c*
    real(kind = rfp) :: vol ! volume of control voflume
    real(kind = rfp) :: centp(ndim) !centeriod of control volume
    !   real(kind=rfp)::rhopp       !Preconditioning dr/dp'
    real(kind = rfp), pointer :: weight(:) !Weights for viscous term
    !   real(rfp)::vj(ndim)  ! electrical current
    !
    integer, pointer :: c2n(:) ! point to nodes
    type(neighbour_cell), pointer :: scell(:) ! surrounding neighbour cells
    type(neighbour), pointer :: sface(:) ! surrounding faces
    real(kind = rfp), pointer :: svd(:,:)
    ! singular value decomposition for least square
    type(vector), pointer :: qvn(:) ! previous qv variable
    real(rfp), pointer :: zmut ! used for Baldwin-Lomax
    real(rfp), pointer :: jv(:) ! electrical current
  end type cell

  type face
    integer :: itype !itype = 0  interior face
    !      = 1  inflow or outflow boundary face
    !      = 2  inviscous boundary face
    !      = 3  viscous wall face
    !      = 5  inlet
    !      = 6  outlet
    type(cell), pointer :: left_cell ! pointer to left cell 
    type(cell), pointer :: right_cell ! pointer to right cell
    ! smaller index of two cells put in left cell
    real(kind = rfp) :: area ! area of this face
    type(matrix) :: ajl, ajr !left and right Jacobian !aav average A
    real(kind = rfp) :: vecn(ndim) ! unit inward normal vector to the left cell
    real(kind = rfp) :: centp(ndim), vl, vr ! volume of left and right cell
    integer, pointer :: f2n(:) ! pointer to nodes
    real(rfp), pointer :: avec(:,:) ! normal vector on the nodes, calculated in face_reconstruction (disc)
  end type face

  type gas
    integer :: np, nt, nc
    real(rfp) :: p0, t0, c0, dp, dt, dc
    real(rfp), pointer :: rho(:,:,:), rhop(:,:,:), rhot(:,:,:), rhoc(:,:,:)
    real(rfp), pointer :: h(:,:,:), hp(:,:,:), ht(:,:,:), hc(:,:,:)
    real(rfp), pointer :: s(:,:,:), mu(:,:,:), lamd(:,:,:), sigma(:,:,:), ne(:,:,:)
    real(rfp), pointer :: satp(:), alpha(:,:,:)
  end type gas

  type solid
    integer :: np, nt, nc
    ! specific heat capacity, thermal conductivity, latent heat
    real(rfp) :: cpl, cps, kl, ks, lheat_melting, lheat_evaporation
    ! density, liquid & solid temperature
    real(rfp) :: rhol, rhos, ts, tl
    ! Liquid viscosity,thermal expansion coefficient, surface tension co., 
    ! Darcy co., Mass diffusivity
    real(rfp) :: mul, beta, sigma, sigmaT, k0, dl, alpha, emi
    real(rfp) :: gravity(ndim)

    ! real and complex refractive index
    real(rfp) :: n, k
    ! Molecular weight
    real(rfp) :: Mw
    ! Reference pressure and temperature for Clayperon-Clasius Equation
    real(rfp) :: pb, Tb
    ! switch for source terms
    integer :: Darcy_on, Thermal_Buoyancy_on, Gravity_on, Enthalpy_Correction_on, Evaporation_on
  end type solid

  type bc_type
    integer :: label, itype, igrp, n
    real(rfp), pointer :: vars(:)
    integer, pointer :: facelist(:)
    real(rfp) :: area
  end type bc_type

  type wallface
    real(kind = rfp) :: vecn(ndim) ! unit inward normal vector to the left cell
    real(kind = rfp) :: centp(ndim)
    real(rfp) :: utawonu, ymax, fmax, umax
  end type wallface

  type wallbl
    integer :: n
    type(wallface), pointer :: faces(:)
  end type wallbl

  type itf ! interface between partiti ns
    integer :: nitf
    integer :: nn, nns !  receive
    integer, pointer :: np(:), nnp(:)
    type(cell), pointer :: pcell(:)

    integer :: nc, ncs
    integer, pointer :: cp(:), ncp(:) ! send
    type(neighbour_cell), pointer :: scell(:)
  end type itf

  type therm
    character*24 :: name
    real(kind = rfp) :: mw, mwinv !--- molecular weight of species
    real(kind = rfp) :: R !--- gas constant
    ! Lennard-Jones potential
    !
    real(kind = rfp) :: cl ! Lennard-Jones characteristic length
    real(kind = rfp) :: ce ! Lennard-Jones characteristic energy divied by
    ! Boltzmann constant epsi / kb
    !
    !--- NASA 4th-order the specific heat at constant pressure
    ! Caution
    ! all coefficient used here are specific properties not the molar units
    ! cp = Cp/MW
    !
    !    real(kind=rfp)::tmin,tcom,tmax
    !    real(kind=rfp)    ::cp_l(npoly) ! lower temperature
    !    real(kind=rfp)    ::cp_h(npoly) ! high temperature
    !--- 4th-order the enthalpy  hh = cp/n
    !
    !    real(kind=rfp)    ::hh_l(npoly + 1) ! lower temperature
    !    real(kind=rfp)    ::hh_h(npoly + 1) ! high temperature

    ! The NASA thermo data file format by
    !Sanford Gordon and Bonnie J. McBride,1994
    !
    integer :: interv
    real(rfp) :: Titv(4) ! Max number of T intervals max is three
    real(rfp) :: cpp(7, 3), hhp(8, 3), ssp(8, 3)
    !
    !  Cp/R =a1*T^-2+a2*T^-1+a3+a4*T+a5*T^2+a6*T^3+a7*T^4
    !  H0/RT=-a1*T^-2+a2*T^-1*ln(T)+a3+a4*T/2+a5*T^2/3+a6*T^3/4+a7*T^4/5+b1/T
    !  S0/R =-a1*T^-2/2+a2*T^-1 +a3*ln(T)+a4*T+a5*T^2/2+a6*T^3/3+a7*T^4/4+b2
    !
    !--- 4th-order the viscosity and thermal conductivity

    integer :: iv, ic
    real(rfp) :: vt(4), ct(4) ! Max number of T intervals is three
    real(kind = rfp) :: mu(4, 3), lamda(4, 3)
    !  LnMU    = AlnT+B/T+C/T^2+D
    !  LnLamda = AlnT+B/T+C/T^2+D
    !
  end type therm

  type react
    ! forward reaction
    integer :: nf, locf(nspe)
    real(kind = rfp) :: nurf(nspe) ! stoichiometri coeficient of reactant
    real(kind = rfp) :: af(3) !,alphf,ef  !Kf=AT^alphexp{-E/RT}
    ! Backward reaction
    integer :: nb, locb(nspe)
    real(kind = rfp) :: nurb(nspe) ! stoichiometri coeficient of product
    real(kind = rfp) :: ab(3) !,alphb,eb  !Kb=AT^alphexp{-E/RT}
    !
    integer :: n, loc(nspe)
    real(kind = rfp) :: nur(nspe) ! stoichiometri coeficient nur = nurb - nurf
    logical :: third_body, kb_equil
  end type react

  type job
    integer :: n
    integer, pointer :: steps(:), ischeme(:), ialg(:), iprecon(:)
    real(rfp), pointer :: cfl(:)
  end type job

end module gems_type
