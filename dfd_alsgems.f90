!**********************************************
!Akward Level-Set GEMS:
!Xuxiao Li
!02/12/2018
!treat level-set equations seperatally 

MODULE dfd_constant
!
!** set some constant
!
 implicit none
! set calculation precision
!  4  -- single precision
!  8  -- double precision
! 12 -- 12 bytes precision (byte = n * bite)
! 16 -- 16 bytes precision 

   integer, parameter :: rfp  = selected_real_kind(8)

   integer, parameter :: ndim = 3   ! - number of dimension
   integer, parameter :: nceq = 1   ! - number of continuous eq
   integer, parameter :: nmeq = 3   ! - number of momentum eq
   integer, parameter :: neeq = 1   ! - number of energy eq
   integer, parameter :: nteq = 0   ! - number of turblent eq
   integer, parameter :: naux = 0   ! - number of auxiliary eq
   integer, parameter :: nrad = 0   ! - number of radiation eq
   integer, parameter :: nspe = 1   ! - number of specis
   integer, parameter :: nspm1 = nspe - 1   ! - number of specis - 1
!   integer, parameter :: nrea = 1   ! - number of reactions
   integer, parameter :: nmag = 0
   integer, parameter :: nele = 0
   integer, parameter :: nfi  = 0
!
   integer, parameter :: neqf  = nceq + nmeq + neeq + nteq + naux + nrad + nspm1
   integer, parameter :: neqm  = nmag + nele + nfi
   integer, parameter :: neq  =  neqf + neqm
                                ! - number of equations to be solved 
   integer, parameter :: neq4   = neq + 5
   integer, parameter :: neq3d  = neq + 4 + ndim
!
   real(rfp), parameter :: pi = 3.141592654
!
   integer, parameter :: ico = 1         ! - the number of continute eq
   integer, parameter :: imb = nceq + 1   ! - the begin number of momunt eq
   integer, parameter :: ime = imb + ndim  - 1   
                                         ! - the end number of momunt eq
   integer, parameter :: ien = nceq + nmeq + 1   ! - the begin number of energy eq
   integer, parameter :: ike = nceq + nmeq + neeq + 1   ! - the begin number of k eq
   integer, parameter :: iom = ike + 1   ! - the begin number of omega eq
   integer, parameter :: isb = nceq + nmeq + neeq + nteq + 1   ! - the begin number of multispecies
   integer, parameter :: ise = isb + nspe - 2 ! - the end number of multispecies
   integer, parameter :: iaub = nceq + nmeq + neeq + nteq + nspm1 +  1   ! - the begin number of auxiliary equation
 
   integer, parameter :: iaue = iaub + naux - 1   ! - the end number of auxiliary equation
   integer, parameter :: irad = nceq + nmeq + neeq + nteq + naux + nspm1 + 1   ! - the begin number of auxiliary equation
!
   integer, parameter :: ibb = irad + nrad         ! the begin number of b equ.
   integer, parameter :: ibe = ibb + 2   ! the end number of b equ
   integer, parameter :: ieb = ibe + 1   ! the begin number of E equ.
   integer, parameter :: iee = ieb + 2   ! the end number of E equ.
   integer, parameter :: ifb = iee + 1   ! the number of Fi equ.
   integer, parameter :: ife = ifb + 1   ! the number of Fi equ.

  !interface capture
   integer, parameter :: usels = 1
!   
   integer, parameter :: inlet_type = 1  !inflow or outflow
   integer, parameter :: outlet_type  = 2
   integer, parameter :: farfield_type = 3
   integer, parameter :: wall_type  = 4
   integer, parameter :: geom_topology = 5
   integer, parameter :: mhd_type   = 6
   integer, parameter :: partition_face_no = 19621011
!
   real(rfp), parameter    :: zero      = 0.0_rfp
   real(rfp), parameter    :: one       = 1.0_rfp
   real(rfp), parameter    :: two       = 2.0_rfp
   real(rfp), parameter    :: three     = 3.0_rfp
   real(rfp), parameter    :: half      = 0.5_rfp
   real(rfp), parameter    :: quarter   = 0.25_rfp
   real(rfp), parameter    :: one_third = one / three
   real(rfp), parameter    :: two_third = two / three
   real(rfp), parameter    :: four_third = two * two_third
   real(rfp), parameter    :: five_third = 5.0_rfp / 3.0_rfp
   real(rfp), parameter    :: three_second = three / two    
   real(rfp), parameter    :: mytiny = 1.e-64_rfp  !tiny(one)
   real(rfp), parameter    :: myhuge = 1.e50_rfp      !huge(one)
  
!
! for reaction
!
  integer, parameter :: npoly = 5   ! 4th-order of NASA polynormial
  real(rfp), parameter :: boltzmann_c = 1.3806503e-23_rfp
  real(rfp), parameter :: sqrt2 = 1.414213562_rfp
  real(rfp), parameter    :: two_sqrt2 = two * sqrt2
!
!   methods
!
  integer, parameter :: mfluid = 0
  integer, parameter :: mturb  = 1
  integer, parameter :: msolid = 2
  integer, parameter :: mmaxwell = 3
  integer, parameter :: mfluid_maxwell = 4
  integer, parameter :: msolid_maxwell = 5
 

END MODULE dfd_constant
