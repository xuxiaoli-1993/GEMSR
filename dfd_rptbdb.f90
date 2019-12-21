!///////////////////////////////////////////////////////////////////////
!/                                                                     /
!/                (NIST) RefProp -- Tree-Based DataBase                /
!/                           (RPTBDB)                                  /
!/                                                                     / 
!/                          version 1.0                                /
!/                                                                     / 
!/                                                                     /
!/ SUMMARY: This program a tree data structure to store and retrieve   /
!/          the NIST RefProp for CFD calcutation                       /
!/ RELEASE: In House Code                                              /
!/ USAGE  : You may freely copy all or part of the code for research   /
!/          and teaching purposes. However, the author does not        /
!/          guarantee results from the code or any part of it. The     / 
!/          code or any part of it may not be sold or used for a       /
!/          commercial purpose without our consent.                    /
!/ AUTHOR : Guoping Xia, Ph. D                                         /
!/ ADVISOR: Charles L. Merkle, Professor                               /
!/ E-MAILS: xiag@purdue.edu                                            / 
!/          merkle@purdue.edu                                          /
!/ ORG    : Purdue University, West Lafayette                          /
!/ DATE OF ORIGINAL WORK:   Aug., 2004                                 /
!/                                                                     /
!***********************  DO NOT REMOVE THIS BANNER  ******************/
!

MODULE RPTBDB_constant
    !
    implicit none
    !
    integer, parameter :: rfp = selected_real_kind(8)
    integer, parameter :: ndim = 2, ndim1 = ndim + 1
    integer, parameter :: FTT_CELLS = 4
    integer, parameter :: FTT_NEIGHBORS = 4
    integer, parameter :: FTT_VERTICES = 4
    integer, parameter :: FTT_MAX_NEW_VERTICES = 2 * ndim**(ndim - 1)
    integer, parameter :: FTT_BOUNDARY = 1, FTT_FINE_FINE = 2, FTT_FINE_COARSE = 3
    integer, parameter :: FTT_RIGHT = 1, FTT_LEFT = 2, FTT_TOP = 3, FTT_BOTTOM = 4
    integer, parameter :: FTT_FRONT = 5, FTT_BACK = 6
    integer, parameter :: FTT_PRE_ORDER = 1, FTT_POST_ORDER = 2
    integer, parameter :: nvar = 9, ni = 3, nj = 3
    integer, parameter :: num_of_bnd_cells = 4
    integer, parameter :: ncmax = 20
    logical, parameter :: FALSE = .false., TRUE = .true.
    real(rfp), parameter :: one = 1.0_rfp, zero = 0.0_rfp, half = 0.5_rfp
    real(rfp), parameter :: tiny = 1.e-15_rfp
    real(rfp), parameter :: pi = 3.1415926

    integer :: tri_cell_vertices(4, 64) = &
    reshape((/ 1, 1, 2, 3, &
    2, 1, 3, 4, &
    3, 1, 2, 5, &
    4, 2, 3, 5, &
    5, 3, 4, 5, &
    6, 1, 6, 4, &
    7, 2, 3, 6, &
    8, 3, 4, 6, &
    9, 1, 2, 7, &
    10, 1, 7, 4, &
    11, 3, 4, 7, &
    12, 1, 2, 8, &
    13, 1, 8, 4, &
    14, 2, 3, 8, &
    15, 1, 6, 5, &
    16, 2, 3, 6, &
    17, 3, 4, 5, &
    18, 3, 5, 6, &
    19, 1, 2, 7, &
    20, 1, 7, 5, &
    21, 3, 4, 5, &
    22, 3, 5, 7, &
    23, 1, 2, 5, &
    24, 2, 3, 8, &
    25, 2, 8, 5, &
    26, 4, 5, 8, &
    27, 1, 6, 4, &
    28, 2, 7, 6, &
    29, 3, 4, 7, &
    30, 4, 6, 7, &
    31, 1, 6, 8, &
    32, 1, 8, 4, &
    33, 2, 3, 6, &
    34, 3, 8, 6, &
    35, 1, 2, 7, &
    36, 1, 8, 4, &
    37, 1, 7, 8, &
    38, 3, 8, 7, &
    39, 1, 6, 5, &
    40, 2, 7, 6, &
    41, 3, 4, 5, &
    42, 3, 5, 7, &
    43, 5, 6, 7, &
    44, 1, 6, 5, &
    45, 2, 3, 6, &
    46, 3, 8, 6, &
    47, 4, 5, 8, &
    48, 5, 6, 8, &
    49, 1, 2, 7, &
    50, 1, 7, 5, &
    51, 3, 8, 7, &
    52, 4, 5, 8, &
    53, 5, 7, 8, &
    54, 1, 6, 8, &
    55, 1, 8, 4, &
    56, 2, 7, 6, &
    57, 3, 8, 7, &
    58, 6, 7, 8, &
    59, 1, 6, 5, &
    60, 2, 7, 6, &
    61, 3, 8, 7, &
    62, 4, 5, 8, &
    63, 5, 6, 8, &
    64, 6, 7, 8 /), (/4, 64/))

    type FttVector
        real(rfp) :: xyz(ndim)
    end type FttVector

    type Qvar
        real(rfp) :: vars(nvar)
    end type Qvar

    type vars
        !  type(FttVector):: pos(FTT_CELLS)
        type(Qvar) :: q(FTT_CELLS)
    end type vars

    type RPTBDB_node
        integer :: id
        type(FttVector) :: pos, org_pos
        type(Qvar) :: q
    end type RPTBDB_node

    type p2node
        type(RPTBDB_node), pointer :: p
    end type p2node

    type cell
        type(p2node) :: c2n(FTT_MAX_NEW_VERTICES + FTT_VERTICES)
    end type cell

    TYPE ACM_COF_C2
        REAL :: Q00(3), Q01(3), Q02(3), Q03(3), Q04(3), Q05(3), &
        Q06(3), Q07(3), Q08(3), Q09(3), &
        Q12(2), Q13(2), Q14(2), Q15(2), Q16(2), Q17(2), Q18(2), &
        Q23(2), Q24(2), Q25(2), Q26(2), Q27(2), &
        Q34(2), Q35(2), Q36(2), &
        Q45(2), Q11, Q22, Q33, Q44, &
        DX(3), DY(3), DX2(3), DY2(3), SL2(3), &
        AP, BP, CP, DP
    END TYPE ACM_COF_C2

    TYPE ACM_COF_C1
        REAL :: P0(3), P1(3), P2(3), P3(3), P4(3), P5(3), &
        P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
        DX(3), DY(3), DX2(3), DY2(3), SL2(3), &
        AP, BP, CP, DP
    END TYPE ACM_COF_C1

    type TriCell
        type(p2node) :: c2n(3)
        type(ACM_COF_C1), pointer :: acf1, acf2
        type(ACM_COF_C2), pointer :: acf
        integer :: icheck, icheck1, icheck2
    end type TriCell

    type p2TriCell
        type(TriCell), pointer :: p
    end type p2TriCell

    type FttCellNeighborsEle
        type(FttCell), pointer :: p
    end type FttCellNeighborsEle

    type FttCellNeighbors
        ! right, left, top, bottom, front, back 
        type(FttCellNeighborsEle) :: e(FTT_NEIGHBORS)
    end type FttCellNeighbors

    type FttCell
        integer :: flags
        integer :: level
        integer :: id

        type(p2node) :: c2n(FTT_MAX_NEW_VERTICES + FTT_VERTICES)
        type(FttCellNeighbors), pointer :: neighbors

        type(FttOct), pointer :: parent, children
        type(FttVector) :: center

        integer :: TriGeom, NumTriCells
        type(p2TriCell), pointer :: tri_cells(:)

    end type FttCell

    type p2FttCell
        type(FttCell), pointer :: p
    end type p2FttCell

    type FttCellChildren
        type(FttCell) :: c(FTT_CELLS)
    end type FttCellChildren

    type FttOct
        integer :: level
        type(FttCell), pointer :: parent

        type(p2FttCell) :: cell(FTT_CELLS)
    end type FttOct

    type SplineSeg
        real(rfp), pointer :: x(:), y(:), cm(:)
    end type SplineSeg

    type FttRootCell

        type(FttCell), pointer :: p
        integer :: nnodes, ncells, nFttCells, nTriCells
        integer, pointer :: children_index(:,:)
        integer :: inew, icst, ix, iy, icmb, if3vars(3)
        real(rfp) :: MolWgt
        type(p2FttCell), pointer :: cells(:)
        type(RPTBDB_node), pointer :: nodes(:)
        type(SplineSeg) :: BndCur(1:4)
        character(len = 17), pointer :: var_nam(:)
        type(FttCell), pointer :: bnd_cell(:)
        integer :: i_con_sel

        real(rfp), pointer :: RPTBDB_WK(:,:), RPTBDB_GIBBS(:), RPTBDB_WK1(:,:), RPTBDB_WK2(:,:)

    end type FttRootCell

    type p2FttRootCell
        type(FttRootCell), pointer :: p
    end type p2FttRootCell


END MODULE RPTBDB_constant

!      ALGORITHM 684, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 16, NO. 3, PP. 253-257.

MODULE ACM_684_C1
    USE RPTBDB_constant

CONTAINS
    !           
    SUBROUTINE C1PNT(ND, XyD, ZD, XP, YP, WK, &
        ID, ZP, ICHECK, &
        P0, P1, P2, P3, P4, P5, &
        P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
        DX, DY, DX2, DY2, SL2, &
        AP, BP, CP, DP)


        REAL :: P0(3), P1(3), P2(3), P3(3), P4(3), P5(3), &
        P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
        DX(3), DY(3), DX2(3), DY2(3), SL2(3), &
        AP, BP, CP, DP

        DIMENSION ZP(3), ID(3)
        REAL(RFP) :: WK(ND, 5)
        DIMENSION PDL(5, 3), X(3), Y(3), Z(3), XYD(3, 2), ZD(3)


        !
        !***********************************************************
        !
        !     INTERPOLATE Z-VALUE OF A POINT IN A C1-REPRESENTATION
        !     -----------------------------------------------------
        !                   ON A SET OF TRIANGLES
        !                   ---------------------
        !
        !!
        !     AUTHOR:
        !     A.PREUSSER
        !     FRITZ-HABER-INSTITUT DER MPG
        !     FARADAYWEG 4-6
        !     D-1000 BERLIN 33
        !
        !
        !
        !
        ! GIVEN A TRIANGULATION OF A SET OF POINTS IN THE PLANE,
        ! IN A DATA STRUCTURE AS DEVELOPED BY R. RENKA FOR ACM ALG. 624,
        ! AND VALUES AND PARTIAL DERIVATIVES AT THESE POINTS,
        ! THIS ROUTINE COMPUTES A VALUE ZP AT XP,YP.
        !
        ! THE PREREQUISITE VALUES CAN BE COMPUTED BY SUBROUTINE C1INI.
        !
        ! INPUT PARAMETERS
        ! ----------------
        !       ND -    NUMBER OF NODES IN THE MESH.
        !               ND .GE. 3.
        !               THE VALUE OF ND MUST BE EXACTLY
        !               THE SAME AS SPECIFIED FOR THE CALL
        !               OF SUBROUTINE C2INI.
        !
        !       XD,YD - VECTORS OF COORDINATES OF THE
        !               NODES IN THE MESH.
        !
        !          ZD - VECTOR OF DATA VALUES AT THE
        !               NODES.
        !
        !       XP,YP - COORDINATES OF A POINT AT WHICH
        !               FUNCTION VALUE ZP IS REQUIRED
        !
        !          PD - ND BY 5 ARRAY CONTAINING
        !               PARTIAL DERIVATIVES AT THE
        !               NODES.
        !               FIRST INDEX INDICATES POINT,
        !               SECOND INDEX DERIVATIVE:
        !               COLUMN 1: ZX
        !                      2: ZY
        !                      3: ZXX
        !                      4: ZXY
        !                      5: ZYY
        !
        !        IREN - INTEGER ARRAY OF LENGTH 7*ND.
        !               CAN BE COMPUTED BY A CALL TO SUBROUTINE
        !               C2INI. IT CONTAINS THE DATA FOR
        !               THE DEFINITION OF THE TRIANGLES.
        !               (ARRAYS IADJ AND IEND OF ALG. 624)
        !
        !         IST - INDEX OF THE STARTING NODE IN
        !               THE SEARCH FOR A TRIANGLE CON-
        !               TAINING (XP,YP).  1 .LE. IST
        !               .LE. ND.  THE OUTPUT VALUE OF
        !               IST FROM A PREVIOUS CALL MAY
        !               BE A GOOD CHOICE.
        !
        !     IADJ AND IEND MAY BE CREATED BY TRMESH (ACM ALG 624)
        !
        !   INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
        !   ROUTINE.
        !
        ! OUTPUT PARAMETERS
        ! -----------------
        !
        !         IST - INDEX OF ONE OF THE VERTICES OF
        !               THE TRIANGLE CONTAINING (XP,YP),
        !
        !          ZP - VALUE TO BE COMPUTED AT (XP,YP)
        !               =0, IF XP,YP IS OUTSIDE THE TRIANGULATION
        !
        !       IOUTS - =0, IF (XP,YP) INSIDE THE TRIANGULATION
        !                   AND ZP CONTAINS THE INTERPOLATED VALUE.
        !               =1, IF (XP,YP) IS OUTSIDE THE TRIANGULATION
        !                   AND ZP=0.
        !
        !     REMARK
        !     ------
        !     THIS ROUTINE SAVES THE INDICES OF THE THREE VERTICES OF
        !     THE TRIANGLE IN WHICH THE POINT XP,YP IS FOUND.
        !     IF, IN THE NEXT CALL, THE SAME INDICES ARE FOUND,
        !     IT ASSUMES THAT /C1PCO/ AND /C1CCO/ HAVE NOT BEEN
        !     ALTERED BY THE USER. THEN THE COEFFICIENTS FOR
        !     THE TRIANGLE ARE NOT CALCULATED BUT TAKEN FROM /C1PCO/
        !     AND /C1CCO/.
        !***********************************************************
        !
        SAVE I1OLD, I2OLD, I3OLD, X, Y
        DATA I1OLD, I2OLD, I3OLD /0, 0, 0/
        !
        IST = ID(1)
        !
        !     LOAD PARTIAL DERIVATIVES
        !
        DO 100 J = 1, 5
            PDL(J, 1) = WK(ID(1), J)
            PDL(J, 2) = WK(ID(2), J)
            PDL(J, 3) = WK(ID(3), J)
            100 CONTINUE
            !
            !     LOAD COORDINATES
            !

            X(1) = XYD(1, 1)
            Y(1) = XYD(1, 2)
            X(2) = XYD(2, 1)
            Y(2) = XYD(2, 2)
            X(3) = XYD(3, 1)
            Y(3) = XYD(3, 2)
            Z(1) = ZD(1)
            Z(2) = ZD(2)
            Z(3) = ZD(3)
            !
            IF (ICHECK == 0) THEN
                !     COMPUTE COEFFICIENTS ALONG SIDES
                CALL C1SIDE(X, Y, Z, PDL, 2, &
                P0, P1, P2, P3, P4, P5, &
                P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
                DX, DY, DX2, DY2, SL2, &
                AP, BP, CP, DP)
                !
                !     COMPUTE COEFFICIENTS FOR INSIDE OF TRIANGLE
                !     AND CONSTANTS FOR TRANSF. CARTESIAN TO TRIANGULAR
                CALL C1INSD(PDL, &
                P0, P1, P2, P3, P4, P5, &
                P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
                DX, DY, DX2, DY2, SL2, &
                AP, BP, CP, DP)
                !
                ICHECK = 1
            ENDIF
            !
            !     EVALUATE POLYNOMIAL AT XP,YP

            XR = XP - X(3)
            YR = YP - Y(3)
            CALL C1HORN(XR, YR, ZP, &
            P0, P1, P2, P3, P4, P5, &
            P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
            DX, DY, DX2, DY2, SL2, &
            AP, BP, CP, DP)
            !
            RETURN
        END SUBROUTINE C1PNT


        SUBROUTINE C1SIDE(X, Y, Z, PDL, NSIDES, &
            P0, P1, P2, P3, P4, P5, &
            P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
            DX, DY, DX2, DY2, SL2, &
            AP, BP, CP, DP)
            !
            REAL :: P0(3), P1(3), P2(3), P3(3), P4(3), P5(3), &
            P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
            DX(3), DY(3), DX2(3), DY2(3), SL2(3), &
            AP, BP, CP, DP
            !
            !
            !     C1-GRID
            !
            !     COMPUTATION OF COEFFICIENTS FOR POLYNOMIALS ALONG SIDES
            !
            !     AUTHOR      : A. PREUSSER, 1989
            !                   FRITZ-HABER-INSTITUT
            !                   DER MAX-PLANCK-GESELLSCHAFT
            !                   FARADAYWEG 4-6
            !                   D-1000 BERLIN 33
            !
            !
            !     X        X-COORDINATES OF VERTICES
            !     Y        Y-COORDINATES OF VERTICES
            !     PDL(J,K) PARTIAL DERIVATIVES AT VERTEX K
            !              (ZX,ZY,ZXX,ZXY,ZYY)
            !     NSIDES   NUMBER OF POLYNOMIALS TO BE DETERMINED
            !              2, FOR SIDE 1 AND 2
            !              3, FOR SIDE 1, 2, AND 3
            !
            !                                  2
            !                                  *
            !!                               /  .
            !                              /     .
            !                            /        .
            !                          /           .
            !                        /              .
            !             SIDE(1)  /                 .
            !                    /                    . SIDE(3)
            !                  /                       .
            !                /                          .
            !              /                             .
            !            /                                .
            !          /                                   .
            !        /                                      .
            !    3 *-----------------------------------------* 1
            !                      SIDE(2)
            !
            !
            !
            !     P0...P5     COEFFICIENTS OF POLYNOMIALS ALONG THE SIDES
            !     P11...P41   COEFFICIENTS FOR BIVARIATE POLYNOMIAL INSIDE TRIANGLE
            !
            COMMON /C1DCO/ ZS(3, 3), ZSS(3, 3)
            !
            !
            !     OUTPUT OF THIS ROUTINE
            !     ----------------------
            !     DX,DY     COORDINATE DIFFERENCES OF ENDPOINTS OF SIDES
            !     DX2,DY2   DX**2,DY**2
            !     SL2       SIDE LENGTH **2
            !     ZS        PARTIAL DERIVATIVES WITH RESPECT TO TRIANGULAR
            !               COORDINATES L1,L2,L3 (CORRESPOND TO DIRECTIONS
            !               OF SIDES 2,1,3 RESPECTIVELY)
            !               FIRST INDEX: VARIABLE
            !               SECOND INDEX: POINT
            !     ZSS       SECOND PARTIAL DERIVATIVES
            !               SIMILAR TO ZS
            !     P0...P5   COEFFICIENTS FOR POLYNOMIALS ALONG THE THREE
            !               SIDES
            !
            !
            DIMENSION PDL(5, 3), X(3), Y(3), Z(3)
            !
            !     COMPUTE COORDINATE DIFFERENCES FOR SIDES
            DO 60 I = 1, 3
                NP1 = 3
                NP2 = 2
                IF (I .EQ. 3) NP1 = 1
                IF (I .EQ. 2) NP2 = 1
                DX(I) = X(NP2) - X(NP1)
                DY(I) = Y(NP2) - Y(NP1)
                DX2(I) = DX(I) * DX(I)
                DY2(I) = DY(I) * DY(I)
                SL2(I) = DX2(I) + DY2(I)
                60 CONTINUE
                !
                !     CONVERTS THE PARTIAL DERIVATIVES AT THE VERTICES OF THE
                !     TRIANGLE FROM THE CARTESIAN X-Y-SYSTEM TO THE TRIANGULAR
                !     L1-L2-SYSTEM USING THE CHAIN RULE
                DXDY1 = 2.0 * DX(1) * DY(1)
                DXDY2 = 2.0 * DX(2) * DY(2)
                DXDY3 = 2.0 * DX(3) * DY(3)
                DO 260 K = 1, 3
                    ZS(1, K) = DX(2) * PDL(1, K) + DY(2) * PDL(2, K)
                    ZS(2, K) = DX(1) * PDL(1, K) + DY(1) * PDL(2, K)
                    ZSS(1, K) = DX2(2) * PDL(3, K) + DXDY2 * PDL(4, K) + DY2(2) * PDL(5, K)
                    ZSS(2, K) = DX2(1) * PDL(3, K) + DXDY1 * PDL(4, K) + DY2(1) * PDL(5, K)
                    260 CONTINUE
                    !
                    !     OPTIONALLY, COMPUTE PARTIAL DERIVATIVES FOR THIRD
                    !     VARIABLE L3 (AT POINTS 1 AND 2)
                    IF (NSIDES .LT. 3) GOTO 280
                    DO 270 K = 1, 2
                        ZS(3, K) = DX(3) * PDL(1, K) + DY(3) * PDL(2, K)
                        ZSS(3, K) = DX2(3) * PDL(3, K) + DXDY3 * PDL(4, K) + DY2(3) * PDL(5, K)
                        270 CONTINUE
                        280 CONTINUE
                        !
                        !     CALCULATES THE COEFFICIENTS OF THE POLYNOMIALS ALONG
                        !     THE TWO (OR THREE) SIDES OF THE TRIANGLE
                        DO 300 I = 1, NSIDES
                            NP1 = 3
                            NP2 = 2
                            IF (I .EQ. 3) NP1 = 1
                            IF (I .EQ. 2) NP2 = 1
                            J = NP2
                            IF (I .EQ. 3) J = 3
                            P0(I) = Z(NP1)
                            P1(I) = ZS(J, NP1)
                            P2(I) = 0.5 * ZSS(J, NP1)
                            H1 = Z(NP2) - P0(I) - P1(I) - P2(I)
                            H2 = ZS(J, NP2) - P1(I) - ZSS(J, NP1)
                            H3 = ZSS(J, NP2) - ZSS(J, NP1)
                            P3(I) = 10.0 * H1 - 4.0 * H2 + 0.5 * H3
                            P4(I) = -15.0 * H1 + 7.0 * H2 - H3
                            P5(I) = 6.0 * H1 - 3.0 * H2 + 0.5 * H3
                            300 CONTINUE
                            !
                            RETURN
                        END SUBROUTINE C1SIDE

                        SUBROUTINE C1INSD(PDL, &
                            P0, P1, P2, P3, P4, P5, &
                            P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
                            DX, DY, DX2, DY2, SL2, &
                            AP, BP, CP, DP)
                            !
                            REAL :: P0(3), P1(3), P2(3), P3(3), P4(3), P5(3), &
                            P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
                            DX(3), DY(3), DX2(3), DY2(3), SL2(3), &
                            AP, BP, CP, DP
                            !
                            !     C1-GRID
                            !
                            !     COMPUTATION OF POLYNOMIAL COEFFICIENTS P11...P41 FOR
                            !     BIVARIATE POLYNOMIAL INSIDE TRIANGLE
                            !     AND
                            !     CONSTANTS FOR TRANSFORMATION OF
                            !     CARTESIAN TO TRIANGULAR COORDINATES
                            !
                            !     AUTHOR      : A. PREUSSER, 1989
                            !                   FRITZ-HABER-INSTITUT
                            !                   DER MAX-PLANCK-GESELLSCHAFT
                            !                   FARADAYWEG 4-6
                            !                   D-1000 BERLIN 33
                            !
                            DIMENSION PDL(5, 3), Z12(3)
                            !
                            !     PDL     PARTIAL DERIVATIVES AT VERTICES 1,2,3
                            !             FIRST INDEX: 1, ZX
                            !                        : 2, ZY
                            !                        : 3, ZXX
                            !                        : 4, ZXY
                            !                        : 5, ZYY
                            !
                            !     Z12     MIXED PARTIAL DERIVATIVES WITH RESPECT TO L1,L2
                            !             AT THE THREE VERTICES 1,2,3
                            !
                            !     P0...P5     COEFFICIENTS OF POLYNOMIALS ALONG THE SIDES
                            !     P11...P41   COEFFICIENTS FOR BIVARIATE POLYNOMIAL INSIDE TRIANGLE
                            !
                            !
                            COMMON /C1DCO/ ZS(3, 3), ZSS(3, 3)
                            !
                            !
                            !     INPUT TO THIS ROUTINE
                            !     ---------------------
                            !     VARIABLES FROM /C1CCO/ AND /C1DCO/ COMPUTED IN C1SIDE
                            !
                            !
                            !     OUTPUT OF THIS ROUTINE
                            !     ----------------------
                            !     P11...P41 POLYNOMIAL COEFFICIENTS FOR INSIDE REPRESENTATION
                            !     AP,BP,CP,DP    CONST. FOR TRANSF. CART. TO TRIANG. SYSTEM
                            !
                            !
                            !--------------------------------------------------------------------
                            !
                            !     COMPUTE CONST. FOR TRANSF. CART.-TRIANG. COORD. SYSTEM
                            AD = DX(2) * DY(1)
                            BC = DX(1) * DY(2)
                            DLT = AD - BC
                            AP = DY(1)/DLT
                            BP = -DX(1)/DLT
                            CP = -DY(2)/DLT
                            DP = DX(2)/DLT
                            !
                            !     COMPUTE MIXED PARTIALS
                            ADBC = AD + BC
                            AB = DX(2) * DX(1)
                            CD = DY(1) * DY(2)
                            DO 100 K = 1, 3
                                100 Z12(K) = AB * PDL(3, K) + ADBC * PDL(4, K) + CD * PDL(5, K)
                                !
                                !     COMPUTE COEFFICIENTS FOR REPRESENTATION INSIDE
                                !             TRIANGLE
                                ABCD5 = 5. * (AB + CD)
                                P14 = P5(1) * ABCD5/SL2(1)
                                P41 = P5(2) * ABCD5/SL2(2)
                                P11 = Z12(3)
                                H1 = ZS(2, 1) - P1(1) - P11 - P41
                                H2 = Z12(1) - P11 - 4.0 * P41
                                P21 = 3.0 * H1 - H2
                                P31 = -2.0 * H1 + H2
                                H1 = ZS(1, 2) - P1(2) - P11 - P14
                                H2 = Z12(2) - P11 - 4.0 * P14
                                P12 = 3.0 * H1 - H2
                                P13 = -2.0 * H1 + H2
                                H1 = 0.5 * ZSS(2, 1) - P2(1) - P12
                                H2 = 0.5 * ZSS(1, 2) - P2(2) - P21
                                E1 = 2.5 * (SL2(2) - SL2(1))/SL2(3) + 2.5
                                !     EQUIVALENT TO
                                !     E1= -5.*(DX(2)*DX(3)+DY(2)*DY(3))/SL2(3)
                                G1 = 3. - E1
                                G2 = -2. + E1
                                G3 = E1 * (P5(1) - P5(2) + P41 - P14) &
                                +P14 - 4. * P41 + 5. * P5(2)
                                P22 = G1 * H1 + G2 * H2 + G3
                                P32 = H1 - P22
                                P23 = H2 - P22
                                RETURN
                            END SUBROUTINE C1INSD

                            SUBROUTINE C1HORN(X, Y, Z, &
                                P0, P1, P2, P3, P4, P5, &
                                P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
                                DX, DY, DX2, DY2, SL2, &
                                AP, BP, CP, DP)
                                !
                                REAL :: P0(3), P1(3), P2(3), P3(3), P4(3), P5(3), &
                                P11, P12, P13, P14, P21, P22, P23, P31, P32, P41, &
                                DX(3), DY(3), DX2(3), DY2(3), SL2(3), &
                                AP, BP, CP, DP, Z(3)
                                !
                                !     C1-GRID
                                !
                                !     EVALUATION OF BIVARIATE QUINTIC POLYNOMIAL
                                !
                                !     AUTHOR      : A. PREUSSER, 1989
                                !                   FRITZ-HABER-INSTITUT
                                !                   DER MAX-PLANCK-GESELLSCHAFT
                                !                   FARADAYWEG 4-6
                                !                   D-1000 BERLIN 33
                                !
                                !
                                !    INPUT PARAMETERS
                                !    ----------------
                                !    X,Y      CARTESIAN COORDINATES RELATIVE TO VERTEX (3)
                                !             OF TRIANGLE
                                !
                                !    OUTPUT PARAMETER
                                !    ----------------
                                !    Z        VALUE OF POLYNOMIAL AT X,Y
                                !
                                !
                                !     INPUT TO THIS ROUTINE
                                !     ---------------------
                                !     ALL VARIABLES FROM /C1PCO/ (POLYNOMIAL COEFFICIENTS)
                                !     AP,BP,CP,DP    CONSTANTS DERIVED FROM COORDINATE DIFFERENCES
                                !                    ON SIDES
                                !
                                ! ------------------------------------------------------------------
                                !
                                !     TRANSFORMATION FROM CARTESIAN TO TRIANGULAR COORDINATES
                                !                                                (U=L1, V=L2)
                                U = AP * X + BP * Y
                                V = CP * X + DP * Y
                                !
                                if (u > 1.0001 .or. u < -0.0001 .or. &
                                    v > 1.0001 .or. v < -0.0001) then
                                    print *, ' out of range '
                                    !pause
                                endif
                                !
                                !     EVALUATION BY HORNER'S SCHEME (BIVARIATE)
                                H0 = P0(1) + V * (P1(1) + V * (P2(1) + V * (P3(1) + V * (P4(1) + V * P5(1)))))
                                H1 = P1(2) + V * (P11 + V * (P12 + V * (P13 + V * P14)))
                                H2 = P2(2) + V * (P21 + V * (P22 + V * P23))
                                H3 = P3(2) + V * (P31 + V * P32)
                                H4 = P4(2) + V * P41
                                H5 = P5(2)
                                Z(1) = H0 + U * (H1 + U * (H2 + U * (H3 + U * (H4 + U * H5))))
                                !     EVALUATION OF ZU      
                                ZU = H1 + U * (2 * H2 + U * (3 * H3 + U * (4 * H4 &
                                +U * 5 * H5)))
                                !
                                !     NOW EXCHANGE TO EVALUATE ZV
                                H0 = P1(1) + V * (2 * P2(1) + V * (3 * P3(1) + V * (4 * P4(1) + V * 5 * P5(1))))
                                H1 = P11 + V * (2 * P12 + V * (3 * P13 + V * 4 * P14))
                                H2 = P21 + V * (2 * P22 + V * 3 * P23)
                                H3 = P31 + V * 2 * P32
                                H4 = P41
                                H5 = 0
                                ZV = H0 + U * (H1 + U * (H2 + U * (H3 + U * (H4 + U * H5))))

                                !     CONVERT DERIVATIVES OF U & V TO DERIVATIVES OF X & Y

                                Z(2) = AP * ZU + CP * ZV
                                Z(3) = BP * ZU + DP * ZV

                                RETURN
                            END SUBROUTINE C1HORN

                        END MODULE ACM_684_C1

                        !      ALGORITHM 684, COLLECTED ALGORITHMS FROM ACM.
                        !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
                        !      VOL. 16, NO. 3, PP. 253-257.

                        MODULE ACM_684_C2
                            USE RPTBDB_constant

                        CONTAINS

                            SUBROUTINE C2PNT(ND, XYD, ZD, XP, YP, WK, &
                                ID, ZP, ICHECK, &
                                Q00, Q01, Q02, Q03, Q04, Q05, &
                                Q06, Q07, Q08, Q09, &
                                Q12, Q13, Q14, Q15, Q16, Q17, Q18, &
                                Q23, Q24, Q25, Q26, Q27, &
                                Q34, Q35, Q36, &
                                Q45, Q11, Q22, Q33, Q44, &
                                DX, DY, DX2, DY2, SL2, &
                                AP, BP, CP, DP)

                                REAL :: Q00(3), Q01(3), Q02(3), Q03(3), Q04(3), Q05(3), &
                                Q06(3), Q07(3), Q08(3), Q09(3), &
                                Q12(2), Q13(2), Q14(2), Q15(2), Q16(2), Q17(2), Q18(2), &
                                Q23(2), Q24(2), Q25(2), Q26(2), Q27(2), &
                                Q34(2), Q35(2), Q36(2), &
                                Q45(2), Q11, Q22, Q33, Q44, &
                                DX(3), DY(3), DX2(3), DY2(3), SL2(3), &
                                AP, BP, CP, DP
                                !
                                !
                                !D    DOUBLE PRECISION ZP
                                DIMENSION ZP(6), ID(3)
                                REAL(RFP) :: WK(ND, 14)
                                DIMENSION PDL(14, 3), X(3), Y(3), Z(3), XYD(3, 2), ZD(3)

                                SAVE I1OLD, I2OLD, I3OLD, X, Y
                                DATA I1OLD, I2OLD, I3OLD /0, 0, 0/

                                I1 = ID(1)
                                I2 = ID(2)
                                I3 = ID(3)
                                IST = I1

                                K1 = 1
                                K2 = 2
                                K3 = 3

                                DO 100 J = 1, 14
                                    PDL(J, K1) = WK(I1, J)
                                    PDL(J, K2) = WK(I2, J)
                                    PDL(J, K3) = WK(I3, J)
                                    100 CONTINUE

                                    X(K1) = XYD(1, 1)
                                    Y(K1) = XYD(1, 2)
                                    X(K2) = XYD(2, 1)
                                    Y(K2) = XYD(2, 2)
                                    X(K3) = XYD(3, 1)
                                    Y(K3) = XYD(3, 2)
                                    Z(K1) = ZD(K1)
                                    Z(K2) = ZD(K2)
                                    Z(K3) = ZD(K3)

                                    IF (ICHECK == 0) THEN

                                        CALL C2SIDE(X, Y, Z, PDL, 2, &
                                        Q00, Q01, Q02, Q03, Q04, Q05, &
                                        Q06, Q07, Q08, Q09, &
                                        Q12, Q13, Q14, Q15, Q16, Q17, Q18, &
                                        Q23, Q24, Q25, Q26, Q27, &
                                        Q34, Q35, Q36, &
                                        Q45, Q11, Q22, Q33, Q44, &
                                        DX, DY, DX2, DY2, SL2, &
                                        AP, BP, CP, DP)

                                        CALL C2INSD(PDL, &
                                        Q00, Q01, Q02, Q03, Q04, Q05, &
                                        Q06, Q07, Q08, Q09, &
                                        Q12, Q13, Q14, Q15, Q16, Q17, Q18, &
                                        Q23, Q24, Q25, Q26, Q27, &
                                        Q34, Q35, Q36, &
                                        Q45, Q11, Q22, Q33, Q44, &
                                        DX, DY, DX2, DY2, SL2, &
                                        AP, BP, CP, DP)

                                    ENDIF

                                    XR = XP - X(3)
                                    YR = YP - Y(3)
                                    CALL C2HORN(XR, YR, ZP, &
                                    Q00, Q01, Q02, Q03, Q04, Q05, &
                                    Q06, Q07, Q08, Q09, &
                                    Q12, Q13, Q14, Q15, Q16, Q17, Q18, &
                                    Q23, Q24, Q25, Q26, Q27, &
                                    Q34, Q35, Q36, &
                                    Q45, Q11, Q22, Q33, Q44, &
                                    DX, DY, DX2, DY2, SL2, &
                                    AP, BP, CP, DP)
                                    RETURN
                                END SUBROUTINE C2PNT

                                SUBROUTINE C2SIDE(X, Y, Z, PDL, NSIDES, &
                                    Q00, Q01, Q02, Q03, Q04, Q05, &
                                    Q06, Q07, Q08, Q09, &
                                    Q12, Q13, Q14, Q15, Q16, Q17, Q18, &
                                    Q23, Q24, Q25, Q26, Q27, &
                                    Q34, Q35, Q36, &
                                    Q45, Q11, Q22, Q33, Q44, &
                                    DX, DY, DX2, DY2, SL2, &
                                    AP, BP, CP, DP)

                                    REAL :: Q00(3), Q01(3), Q02(3), Q03(3), Q04(3), Q05(3), &
                                    Q06(3), Q07(3), Q08(3), Q09(3), &
                                    Q12(2), Q13(2), Q14(2), Q15(2), Q16(2), Q17(2), Q18(2), &
                                    Q23(2), Q24(2), Q25(2), Q26(2), Q27(2), &
                                    Q34(2), Q35(2), Q36(2), &
                                    Q45(2), Q11, Q22, Q33, Q44, &
                                    DX(3), DY(3), DX2(3), DY2(3), SL2(3), &
                                    AP, BP, CP, DP


                                    !*D    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
                                    REAL X, Y, Z, PDL

                                    COMMON /C2DCO/ Z10(3, 3), Z20(3, 3), Z30(3, 3), Z40(3, 3), &
                                    Z21(2, 3), Z31(2, 3), Z11(3), Z22(3)

                                    DIMENSION PDL(14, 3), X(3), Y(3), Z(3)
                                    DO 60 I = 1, 3
                                        NP1 = 3
                                        NP2 = 2
                                        IF (I .EQ. 3) NP1 = 1
                                        IF (I .EQ. 2) NP2 = 1
                                        XNP1 = X(NP1)
                                        YNP1 = Y(NP1)
                                        DX(I) = X(NP2) - XNP1
                                        DY(I) = Y(NP2) - YNP1
                                        DX2(I) = DX(I) * DX(I)
                                        DY2(I) = DY(I) * DY(I)
                                        SL2(I) = DX2(I) + DY2(I)
                                        60 CONTINUE

                                        DO 260 I = 1, NSIDES
                                            J = 3 - MOD(I, 3)
                                            E2XY = 2. * DX(J) * DY(J)
                                            DX3 = DX2(J) * DX(J)
                                            E3X2Y = 3. * DX2(J) * DY(J)
                                            E3Y2X = 3. * DY2(J) * DX(J)
                                            DY3 = DY2(J) * DY(J)
                                            DX4 = DX3 * DX(J)
                                            DY4 = DY3 * DY(J)
                                            D4X3Y = 4. * DX3 * DY(J)
                                            D4Y3X = 4. * DY3 * DX(J)
                                            D6X2Y2 = 6. * DX2(J) * DY2(J)

                                            NK = 3
                                            IF (I .GT. 2) NK = 2
                                            DO 230 K = 1, NK
                                                Z10(I, K) = DX(J) * PDL(1, K) + DY(J) * PDL(2, K)
                                                Z20(I, K) = DX2(J) * PDL(3, K) + E2XY * PDL(4, K) + DY2(J) * PDL(5, K)
                                                Z30(I, K) = DX3 * PDL(6, K) + E3X2Y * PDL(7, K) + &
                                                DY3 * PDL(9, K) + E3Y2X * PDL(8, K)
                                                Z40(I, K) = DX4 * PDL(10, K) + D4X3Y * PDL(11, K) + &
                                                DY4 * PDL(14, K) + D4Y3X * PDL(13, K) + &
                                                D6X2Y2 * PDL(12, K)
                                                230 CONTINUE

                                                260 CONTINUE

                                                DO 300 I = 1, NSIDES
                                                    J = 3 - MOD(I, 3)

                                                    NP1 = 3
                                                    NP2 = 2
                                                    IF (I .EQ. 3) NP1 = 1
                                                    IF (I .EQ. 2) NP2 = 1

                                                    Q00(I) = Z(NP1)
                                                    Q01(I) = Z10(J, NP1)
                                                    Q02(I) = Z20(J, NP1) * 0.5
                                                    Q03(I) = Z30(J, NP1)/6.
                                                    Q04(I) = Z40(J, NP1)/24.
                                                    Z1MZ2 = Q00(I) - Z(NP2)
                                                    HQ4 = Q04(I) * 5.
                                                    HQ42 = HQ4 + HQ4
                                                    Q05(I) = Z40(J, NP2)/24. - Z30(J, NP2) + 10.5 * Z20(J, NP2) &
                                                    -56. * Z10(J, NP2) - HQ4 &
                                                    -15. * Q03(I) - 35. * Q02(I) - 70. * Q01(I) - 126. * Z1MZ2
                                                    Q06(I) = (-Z40(J, NP2) + 23. * Z30(J, NP2))/6. - 38.5 * Z20(J, NP2) &
                                                    +196. * Z10(J, NP2) + HQ42 &
                                                    +40. * Q03(I) + 105. * Q02(I) + 224. * Q01(I) + 420. * Z1MZ2
                                                    Q07(I) = 0.25 * Z40(J, NP2) - 5.5 * Z30(J, NP2) + 53. * Z20(J, NP2) &
                                                    -260. * Z10(J, NP2) - HQ42 &
                                                    -45. * Q03(I) - 126. * Q02(I) - 280. * Q01(I) - 540 * Z1MZ2
                                                    Q08(I) = -Z40(J, NP2)/6. + 3.5 * Z30(J, NP2) - 32.5 * Z20(J, NP2) &
                                                    +155. * Z10(J, NP2) + HQ4 &
                                                    +24. * Q03(I) + 70. * Q02(I) + 160. * Q01(I) + 315. * Z1MZ2
                                                    Q09(I) = (Z40(J, NP2) - 20 * Z30(J, NP2))/24. + 7.5 * Z20(J, NP2) &
                                                    -35. * Z10(J, NP2) - Q04(I) &
                                                    -5. * Q03(I) - 15. * Q02(I) - 35. * Q01(I) - 70. * Z1MZ2

                                                    300 CONTINUE

                                                    RETURN
                                                END SUBROUTINE C2SIDE

                                                SUBROUTINE C2INSD(PDL, &
                                                    Q00, Q01, Q02, Q03, Q04, Q05, &
                                                    Q06, Q07, Q08, Q09, &
                                                    Q12, Q13, Q14, Q15, Q16, Q17, Q18, &
                                                    Q23, Q24, Q25, Q26, Q27, &
                                                    Q34, Q35, Q36, &
                                                    Q45, Q11, Q22, Q33, Q44, &
                                                    DX, DY, DX2, DY2, SL2, &
                                                    AP, BP, CP, DP)

                                                    REAL :: Q00(3), Q01(3), Q02(3), Q03(3), Q04(3), Q05(3), &
                                                    Q06(3), Q07(3), Q08(3), Q09(3), &
                                                    Q12(2), Q13(2), Q14(2), Q15(2), Q16(2), Q17(2), Q18(2), &
                                                    Q23(2), Q24(2), Q25(2), Q26(2), Q27(2), &
                                                    Q34(2), Q35(2), Q36(2), &
                                                    Q45(2), Q11, Q22, Q33, Q44, &
                                                    DX(3), DY(3), DX2(3), DY2(3), SL2(3), &
                                                    AP, BP, CP, DP

                                                    !*D    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
                                                    REAL PDL
                                                    DIMENSION PDL(14, 3)

                                                    COMMON /C2DCO/ Z10(3, 3), Z20(3, 3), Z30(3, 3), Z40(3, 3), &
                                                    Z21(2, 3), Z31(2, 3), Z11(3), Z22(3)
                                                    !
                                                    PARAMETER (ONE = 1., TWO = 2., FIVE = 5.)

                                                    DIMENSION XJ2XI(2), YJ2YI(2)

                                                    ADBC = AD + BC
                                                    AD = DX(2) * DY(1)
                                                    BC = DX(1) * DY(2)
                                                    DLT = AD - BC
                                                    AP = DY(1)/DLT
                                                    BP = -DX(1)/DLT
                                                    CP = -DY(2)/DLT
                                                    DP = DX(2)/DLT

                                                    AB = DX(2) * DX(1)
                                                    CD = DY(1) * DY(2)
                                                    ABCD = AB + CD

                                                    DO 100 I = 1, 2
                                                        J = 3 - I
                                                        XJ2XI(I) = DX2(J) * DX(I)
                                                        YJ2YI(I) = DY2(J) * DY(I)
                                                        XJ2YI = DX2(J) * DY(I)
                                                        YJ2XI = DY2(J) * DX(I)
                                                        E2XXY = AB * DY(J) * 2. + XJ2YI
                                                        E2YYX = CD * DX(J) * 2. + YJ2XI
                                                        XJ3XI = XJ2XI(I) * DX(J)
                                                        YJ3YI = YJ2YI(I) * DY(J)
                                                        E3X2XY = 3. * XJ2XI(I) * DY(J) + XJ2YI * DX(J)
                                                        E3Y2YX = 3. * YJ2YI(I) * DX(J) + YJ2XI * DY(J)
                                                        E3J2 = 3. * (XJ2YI * DY(J) + YJ2XI * DX(J))
                                                        DO 40 K = 1, 3
                                                            Z21(I, K) = XJ2XI(I) * PDL(6, K) + YJ2YI(I) * PDL(9, K) + &
                                                            E2XXY * PDL(7, K) + E2YYX * PDL(8, K)
                                                            Z31(I, K) = XJ3XI * PDL(10, K) + YJ3YI * PDL(14, K) + &
                                                            E3X2XY * PDL(11, K) + E3Y2YX * PDL(13, K) + &
                                                            E3J2 * PDL(12, K)
                                                            40 CONTINUE
                                                            100 CONTINUE

                                                            DX1X22 = DX2(1) * DX2(2)
                                                            DY1Y22 = DY2(1) * DY2(2)
                                                            E2IJX = 2. * (XJ2XI(1) * DY(1) + XJ2XI(2) * DY(2))
                                                            E2IJY = 2. * (YJ2YI(1) * DX(1) + YJ2YI(2) * DX(2))
                                                            E12 = DX2(2) * DY2(1) + DX2(1) * DY2(2) + 4. * AB * CD
                                                            DO 200 K = 1, 3
                                                                Z11(K) = AB * PDL(3, K) + ADBC * PDL(4, K) + CD * PDL(5, K)
                                                                Z22(K) = DX1X22 * PDL(10, K) + DY1Y22 * PDL(14, K) + &
                                                                E2IJX * PDL(11, K) + E2IJY * PDL(13, K) + &
                                                                E12 * PDL(12, K)
                                                                200 CONTINUE

                                                                Q11 = Z11(3)
                                                                Q22 = Z22(3) * 0.25
                                                                T6 = ONE/6.
                                                                DO 500 I = 1, 2
                                                                    J = 3 - I

                                                                    B = SL2(I)
                                                                    D = ABCD
                                                                    BDNO = (ABS(B) + ABS(D)) * 0.5
                                                                    B = B/BDNO
                                                                    D = D/BDNO
                                                                    DB = D/B
                                                                    DB2 = DB * DB

                                                                    Q18(I) = 9. * Q09(I) * DB
                                                                    Q12(I) = Z21(J, 3) * 0.5
                                                                    Q13(I) = Z31(J, 3) * T6

                                                                    Q184 = Q18(I) * 4.
                                                                    Q134 = Q13(I) * 4.
                                                                    Z31JJ6 = Z31(J, J) * T6
                                                                    Z31JJ2 = Z31(J, J) * 0.5
                                                                    ZTDIF = Z10(I, 3) - Z10(I, J)
                                                                    Q14(I) = -Z31JJ6 + 2.5 * Z21(J, J) - 15. * Z11(J) + Q18(I) &
                                                                    -Q134 - 10. * Q12(I) - 20. * Q11 - 35. * ZTDIF

                                                                    Q15(I) = Z31JJ2 - 7. * Z21(J, J) + 39. * Z11(J) - Q184 &
                                                                    +6. * Q13(I) + 20. * Q12(I) + 45. * Q11 + 84. * ZTDIF

                                                                    Q16(I) = -Z31JJ2 + 6.5 * Z21(J, J) - 34. * Z11(J) + 6. * Q18(I) &
                                                                    -Q134 - 15. * Q12(I) - 36. * Q11 - 70. * ZTDIF

                                                                    Q17(I) = Z31JJ6 - 2. * Z21(J, J) + 10. * Z11(J) - Q184 &
                                                                    +Q13(I) + 4. * Q12(I) + 10. * Q11 + 20. * ZTDIF

                                                                    Q26(I) = 7. * Q17(I) * DB - 28. * Q08(I) * DB2
                                                                    Q27(I) = 8. * Q18(I) * DB - 36. * Q09(I) * DB2

                                                                    Z20JM3 = Z20(I, J) - Z20(I, 3)
                                                                    Z224 = Z22(J) * 0.25
                                                                    Q263 = Q26(I) * 3.
                                                                    Q223 = Q22 * 3.
                                                                    Q25(I) = -6. * Q27(I) - Q263 - Q22 &
                                                                    -1.5 * Z21(I, 3) + Z224 &
                                                                    -1.5 * Z21(I, J) + 3. * Z20JM3

                                                                    Q24(I) = 8. * Q27(I) + Q263 + Q223 &
                                                                    +4. * Z21(I, 3) - 0.5 * Z22(J) &
                                                                    +3.5 * Z21(I, J) - 7.5 * Z20JM3

                                                                    Q23(I) = -3. * Q27(I) - Q26(I) - Q223 &
                                                                    -3. * Z21(I, 3) + Z224 &
                                                                    -2. * Z21(I, J) + 5. * Z20JM3
                                                                    500 CONTINUE

                                                                    B = DX(1) * DX(3) + DY(1) * DY(3)
                                                                    D = -DX(2) * DX(3) - DY(2) * DY(3)
                                                                    BDNO = (ABS(B) + ABS(D)) * 0.5
                                                                    B = B/BDNO
                                                                    D = D/BDNO
                                                                    BB = B * B
                                                                    DD = D * D
                                                                    BBPDD = BB + DD
                                                                    BPD = B + D
                                                                    BPD2 = BPD * BPD
                                                                    BMD = B - D
                                                                    BMD2 = BMD * BMD
                                                                    BD = B * D

                                                                    T3 = ONE/3.
                                                                    T23 = TWO/3.

                                                                    BB9 = BB * 9.
                                                                    BB11 = BB * 11.

                                                                    DD2 = DD * 2.
                                                                    DD18 = DD * 18.

                                                                    BD3 = BD * 3.
                                                                    BD9 = BD * 9.
                                                                    BD15 = BD * 15.

                                                                    C1 = DD * 32. - BD15 + BB * 25.
                                                                    C2 = DD * 10. + BD15 + BB * 17.
                                                                    C3 = DD * 12. - BD3 + BB * 5.
                                                                    C4 = DD2 + BD3 + BB9
                                                                    C5 = DD * 25. - BD9 + BB * 18.
                                                                    C6 = DD * 6. + BD15 + BB * 13.
                                                                    C7 = DD18 - BD3 + BB11
                                                                    C8 = DD2 + BD15 + BB9
                                                                    C9 = DD18 + BD9 + BB11
                                                                    C0 = DD * 8. + BD15 + BB * 15.

                                                                    PART1 = (Q03(1) * C1 + Q03(2) * C2 + &
                                                                    Q04(1) * C3 + Q04(2) * C4 + &
                                                                    Q13(1) * C5 + Q13(2) * C6 + &
                                                                    Q14(1) * C3 + Q14(2) * C4 + &
                                                                    Q23(1) * C7 + Q23(2) * C8 + &
                                                                    Q24(1) * C9 + Q24(2) * C0) * 3.
                                                                    PART2 = Z31(2, 1)*(3.5 * BBPDD - BD3) + &
                                                                    Z31(1, 2) * BBPDD * 2. + &
                                                                    Z40(1, 2)*(-DD * 0.25 - BD * 0.375 - BB * 1.125) + &
                                                                    Z40(2, 1)*(-DD * 1.5 + BD * 0.375 - BB * 0.625) - &
                                                                    (Z30(1, 2) * C2 + Z30(2, 1) * C1) * 0.5
                                                                    PART3 = (14. * Q27(2) + 27. * Q26(2) + 30. * Q25(2) &
                                                                    -T3 * Q18(2) - 84. * Q27(1) - 15. * Q26(1) &
                                                                    +30. * Q25(1) + (226. + T23) * Q18(1) + 63. * Q17(1) &
                                                                    -399. * Q09(1) - 84. * Q08(1)) * DD
                                                                    PART4 = (151. * Q27(2) + 96. * Q26(2) + 60. * Q25(2) &
                                                                    -(81. + T23) * Q18(2) - 21. * Q17(2) - 3. * Q09(2) &
                                                                    +129. * Q27(1) + 96. * Q26(1) + 60. * Q25(1) &
                                                                    -(53. + T23) * Q18(1) - 21. * Q17(1) - 39. * Q09(1)) * BD
                                                                    PART5 = (-91. * Q27(2) - 15. * Q26(2) + 30. * Q25(2) &
                                                                    +(230. + T23) * Q18(2) + 63. * Q17(2) - 399. * Q09(2) &
                                                                    -84. * Q08(2) + 21. * Q27(1) + 27. * Q26(1) &
                                                                    +30. * Q25(1) - (4. + T3) * Q18(1)) * BB
                                                                    Q45(1) = -(PART1 + PART2 + PART3 + PART4 + PART5)/BPD2

                                                                    Q45(2) = Q45(1) + &
                                                                    ((-36. * (Q09(1) - Q09(2)) + 20. * (Q18(1) - Q18(2)) &
                                                                    -8. * (Q27(1) - Q27(2))) * BD &
                                                                    +(45. * (Q03(1) - Q03(2)) + 9. * (Q04(1) - Q04(2)) &
                                                                    +36. * (Q13(1) - Q13(2)) + 9. * (Q14(1) - Q14(2)) &
                                                                    +27. * (Q23(1) - Q23(2)) + 9. * (Q24(1) - Q24(2)) &
                                                                    +7. * (Q27(1) - Q27(2)) - 4. * (Q18(1) - Q18(2)) &
                                                                    +0.375 * (Z40(1, 2) - Z40(2, 1)) + 7.5 * (Z30(1, 2) - Z30(2, 1)) &
                                                                    -1.5 * (Z31(1, 2) - Z31(2, 1))) * BMD2)/BPD2
                                                                    DM2B = D - B * 2

                                                                    PART1 = ((Q03(2) - Q03(1)) * 5. + (Q13(2) - Q13(1)) * 4. + &
                                                                    (Q23(2) - Q23(1)) * 3. + Q24(2) - Q24(1) + &
                                                                    Q14(2) - Q14(1) + Q04(2) - Q04(1)) * DM2B * 24.
                                                                    PART2 = ((Z30(2, 1) - Z30(1, 2)) * 20. + (Z31(1, 2) - Z31(2, 1)) * 4. + &
                                                                    Z40(2, 1) - Z40(1, 2)) * DM2B
                                                                    PART3 = (2. * Q27(2) - Q18(2) - 7. * Q27(1) + 8. * Q18(1) - 9. * Q09(1)) * D * 8.
                                                                    PART4 = (2. * Q27(1) - Q18(1) - 7. * Q27(2) + 8. * Q18(2) - 9. * Q09(2)) * B * 8.
                                                                    PART5 = (Q45(2) * 8. - Q45(1) * 16.) * BPD
                                                                    Q36(1) = -(PART1 + PART2 + PART3 + PART4 + PART5)/(24. * BPD)

                                                                    T24 = ONE/24.
                                                                    T56 = FIVE/6.
                                                                    Q36(2) = -Z40(2, 1) * T24 + Z31(2, 1) * T6 &
                                                                    -Z30(2, 1) * T56 + Q45(2) - Q24(2) &
                                                                    -3. * Q23(2) - Q14(2) - 4. * Q13(2) - Q04(2) &
                                                                    -5. * Q03(2) + Z40(1, 2) * T24 - Z31(1, 2) * T6 &
                                                                    +Z30(1, 2) * T56 - Q45(1) + Q36(1) &
                                                                    +Q24(1) + 3. * Q23(1) + Q14(1) + 4. * Q13(1) &
                                                                    +Q04(1) + 5. * Q03(1)

                                                                    PART1 = (5. * Z40(2, 1) - 14. * Z31(2, 1) + 7. * Z40(1, 2))/36.
                                                                    PART2 = 10. * (-Q26(2) - Q25(2) - Q26(1) - Q25(1)) - &
                                                                    8. * Q23(2) - T23 * Q36(1) - 6. * Q23(1)
                                                                    PART3 = (8. * Z30(2, 1) + 5. * Q45(2) - 35. * Q27(2) &
                                                                    -48. * (Q24(2) + Q13(2) + Q03(1)) &
                                                                    -21. * (Q14(2) + Q04(2)) - 60. * Q03(2) &
                                                                    -2. * Z31(1, 2) + 10. * Z30(1, 2) &
                                                                    -Q45(1) - 35. * Q27(1) - 42. * Q24(1) &
                                                                    -15. * Q14(1) - 37.5 * Q13(1) - 15. * Q04(1))/4.5
                                                                    Q35(2) = -PART1 - PART2 - PART3

                                                                    PART1 = 2. * (Q36(2) - Q36(1) + Q23(2) - Q23(1))
                                                                    PART2 = 3. * (Q13(2) - Q13(1))
                                                                    PART3 = 4. * (Q03(2) - Q03(1))
                                                                    PART4 = T23 * (Z30(2, 1) - Z30(1, 2))
                                                                    PART5 = T6 * (Z31(1, 2) - Z31(2, 1))
                                                                    Q35(1) = Q35(2) + &
                                                                    PART1 + PART2 + PART3 + PART4 + PART5

                                                                    Q33 = -2. * Q23(2) - 3. * Q13(2) - 4. * Q03(2) &
                                                                    -Z31(1, 2) * T6 + T23 * Z30(1, 2) + 2. * Q36(1) + Q35(1)

                                                                    Q34(2) = -0.75 * Q33 + Z31(2, 1) * T24 &
                                                                    -1.5 * Q36(2) - 1.25 * Q35(2) - 0.5 * Q23(1) - 0.25 * Q13(1)

                                                                    Q34(1) = -0.75 * Q33 - 0.5 * Q23(2) &
                                                                    -0.25 * Q13(2) + Z31(1, 2) * T24 - 1.5 * Q36(1) - 1.25 * Q35(1)

                                                                    Q44 = +Z40(2, 1) * T24 - Q45(2) &
                                                                    -Q23(2) - 2. * Q13(2) - 3. * Q03(2) - Z31(1, 2) * T6 &
                                                                    +0.5 * Z30(1, 2) + 3. * Q36(1) + 2. * Q35(1) &
                                                                    -Q24(1) - Q14(1) - Q04(1)

                                                                    RETURN
                                                                END SUBROUTINE C2INSD

                                                                SUBROUTINE C2HORN(X, Y, Z, &
                                                                    Q00, Q01, Q02, Q03, Q04, Q05, &
                                                                    Q06, Q07, Q08, Q09, &
                                                                    Q12, Q13, Q14, Q15, Q16, Q17, Q18, &
                                                                    Q23, Q24, Q25, Q26, Q27, &
                                                                    Q34, Q35, Q36, &
                                                                    Q45, Q11, Q22, Q33, Q44, &
                                                                    DX, DY, DX2, DY2, SL2, &
                                                                    AP, BP, CP, DP)

                                                                    REAL :: Q00(3), Q01(3), Q02(3), Q03(3), Q04(3), Q05(3), &
                                                                    Q06(3), Q07(3), Q08(3), Q09(3), &
                                                                    Q12(2), Q13(2), Q14(2), Q15(2), Q16(2), Q17(2), Q18(2), &
                                                                    Q23(2), Q24(2), Q25(2), Q26(2), Q27(2), &
                                                                    Q34(2), Q35(2), Q36(2), &
                                                                    Q45(2), Q11, Q22, Q33, Q44, &
                                                                    DX(3), DY(3), DX2(3), DY2(3), SL2(3), &
                                                                    AP, BP, CP, DP, Z(6)

                                                                    !*D    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
                                                                    REAL X, Y
                                                                    U = AP * X + BP * Y
                                                                    V = CP * X + DP * Y

                                                                    if (u > 1 .or. u < 0 .or. v > 1 .or. v < 0) then
                                                                        print *, ' out of range '
                                                                        !pause
                                                                    endif

                                                                    H0 = Q00(2) + V * (Q01(1) + V * (Q02(1) + V * (Q03(1) + V * (Q04(1) + V * (Q05(1) &
                                                                    +V * (Q06(1) + V * (Q07(1) + V * (Q08(1) + V * Q09(1)))))))))
                                                                    H1 = Q01(2) + V * (Q11 + V * (Q12(1) + V * (Q13(1) + V * (Q14(1) &
                                                                    +V * (Q15(1) + V * (Q16(1) + V * (Q17(1) + V * Q18(1))))))))
                                                                    H2 = Q02(2) + V * (Q12(2) + V * (Q22 + V * (Q23(1) + V * (Q24(1) &
                                                                    +V * (Q25(1) + V * (Q26(1) + V * Q27(1)))))))
                                                                    H3 = Q03(2) + V * (Q13(2) + V * (Q23(2) + V * (Q33 + V * (Q34(1) &
                                                                    +V * (Q35(1) + V * Q36(1))))))
                                                                    H4 = Q04(2) + V * (Q14(2) + V * (Q24(2) + V * (Q34(2) + V * (Q44 &
                                                                    +V * Q45(1)))))
                                                                    H5 = Q05(2) + V * (Q15(2) + V * (Q25(2) + V * (Q35(2) + V * Q45(2))))
                                                                    H6 = Q06(2) + V * (Q16(2) + V * (Q26(2) + V * Q36(2)))
                                                                    H7 = Q07(2) + V * (Q17(2) + V * Q27(2))
                                                                    H8 = Q08(2) + V * Q18(2)
                                                                    H9 = Q09(2)
                                                                    Z(1) = H0 + U * (H1 + U * (H2 + U * (H3 + U * (H4 &
                                                                    +U * (H5 + U * (H6 + U * (H7 + U * (H8 + U * H9))))))))
                                                                    ZU = H1 + U * (2 * H2 + U * (3 * H3 + U * (4 * H4 &
                                                                    +U * (5 * H5 + U * (6 * H6 + U * (7 * H7 + U * (8 * H8 + 9 * U * H9)))))))

                                                                    ZUU = 2 * 1 * H2 + U * (3 * 2 * H3 + U * (4 * 3 * H4 + U * (5 * 4 * H5 + U * (6 * 5 * H6 &
                                                                    +U * (7 * 6 * H7 + U * (8 * 7 * H8 + 9 * 8 * U * H9))))))

                                                                    H0 = Q01(1) + V * (2 * Q02(1) + V * (3 * Q03(1) + V * (4 * Q04(1) &
                                                                    +V * (5 * Q05(1) + V * (6 * Q06(1) + V * (7 * Q07(1) &
                                                                    +V * (8 * Q08(1) + V * 9 * Q09(1))))))))
                                                                    H1 = Q11 + V * (2 * Q12(1) + V * (3 * Q13(1) + V * (4 * Q14(1) + V * (5 * Q15(1) &
                                                                    +V * (6 * Q16(1) + V * (7 * Q17(1) + V * 8 * Q18(1)))))))
                                                                    H2 = Q12(2) + V * (2 * Q22 + V * (3 * Q23(1) + V * (4 * Q24(1) &
                                                                    +V * (5 * Q25(1) + V * (6 * Q26(1) + V * 7 * Q27(1))))))
                                                                    H3 = Q13(2) + V * (2 * Q23(2) + V * (3 * Q33 + V * (4 * Q34(1) &
                                                                    +V * (5 * Q35(1) + V * 6 * Q36(1)))))
                                                                    H4 = Q14(2) + V * (2 * Q24(2) + V * (3 * Q34(2) + V * (4 * Q44 &
                                                                    +V * 5 * Q45(1))))
                                                                    H5 = Q15(2) + V * (2 * Q25(2) + V * (3 * Q35(2) + V * 4 * Q45(2)))
                                                                    H6 = Q16(2) + V * (2 * Q26(2) + V * 3 * Q36(2))
                                                                    H7 = Q17(2) + V * 2 * Q27(2)
                                                                    H8 = Q18(2)
                                                                    H9 = 0


                                                                    ZV = H0 + U * (H1 + U * (H2 + U * (H3 + U * (H4 &
                                                                    +U * (H5 + U * (H6 + U * (H7 + U * (H8 + U * H9))))))))
                                                                    H0 = 2 * 1 * Q02(1) + V * (3 * 2 * Q03(1) + V * (4 * 3 * Q04(1) &
                                                                    +V * (5 * 4 * Q05(1) + V * (6 * 5 * Q06(1) + V * (7 * 6 * Q07(1) &
                                                                    +V * (8 * 7 * Q08(1) + V * 9 * 8 * Q09(1)))))))
                                                                    H1 = 2 * 1 * Q12(1) + V * (3 * 2 * Q13(1) + V * (4 * 3 * Q14(1) + V * (5 * 4 * Q15(1) &
                                                                    +V * (6 * 5 * Q16(1) + V * (7 * 6 * Q17(1) + V * 8 * 7 * Q18(1))))))
                                                                    H2 = 2 * 1 * Q22 + V * (3 * 2 * Q23(1) + V * (4 * 3 * Q24(1) &
                                                                    +V * (5 * 4 * Q25(1) + V * (6 * 5 * Q26(1) + V * 7 * 6 * Q27(1)))))
                                                                    H3 = 2 * 1 * Q23(2) + V * (3 * 2 * Q33 + V * (4 * 3 * Q34(1) &
                                                                    +V * (5 * 4 * Q35(1) + V * 6 * 5 * Q36(1))))
                                                                    H4 = 2 * 1 * Q24(2) + V * (3 * 2 * Q34(2) + V * (4 * 3 * Q44 &
                                                                    +V * 5 * 4 * Q45(1)))
                                                                    H5 = 2 * 1 * Q25(2) + V * (3 * 2 * Q35(2) + V * 4 * 3 * Q45(2))
                                                                    H6 = 2 * 1 * Q26(2) + V * 3 * 2 * Q36(2)
                                                                    H7 = 2 * 1 * Q27(2)
                                                                    H8 = 0
                                                                    H9 = 0


                                                                    ZVV = H0 + U * (H1 + U * (H2 + U * (H3 + U * (H4 &
                                                                    +U * (H5 + U * (H6 + U * (H7 + U * (H8 + U * H9))))))))

                                                                    H0 = 0
                                                                    H1 = Q11 + V * (2 * Q12(1) + V * (3 * Q13(1) + V * (4 * Q14(1) &
                                                                    +V * (5 * Q15(1) + V * (6 * Q16(1) + V * (7 * Q17(1) &
                                                                    +V * 8 * Q18(1)))))))
                                                                    H2 = Q12(2) + V * (2 * Q22 + V * (3 * Q23(1) + V * (4 * Q24(1) &
                                                                    +V * (5 * Q25(1) + V * (6 * Q26(1) + V * 7 * Q27(1))))))
                                                                    H3 = Q13(2) + V * (2 * Q23(2) + V * (3 * Q33 + V * (4 * Q34(1) &
                                                                    +V * (5 * Q35(1) + V * 6 * Q36(1)))))
                                                                    H4 = Q14(2) + V * (2 * Q24(2) + V * (3 * Q34(2) + V * (4 * Q44 &
                                                                    +V * 5 * Q45(1))))
                                                                    H5 = Q15(2) + V * (2 * Q25(2) + V * (3 * Q35(2) + V * 4 * Q45(2)))
                                                                    H6 = Q16(2) + V * (2 * Q26(2) + V * 3 * Q36(2))
                                                                    H7 = Q17(2) + V * 2 * Q27(2)
                                                                    H8 = Q18(2)
                                                                    H9 = 0

                                                                    ZUV = H1 + U * (2 * H2 + U * (3 * H3 + U * (4 * H4 &
                                                                    +U * (5 * H5 + U * (6 * H6 + U * (7 * H7 + U * (8 * H8 + 9 * U * H9)))))))


                                                                    Z(2) = AP * ZU + CP * ZV
                                                                    Z(3) = BP * ZU + DP * ZV

                                                                    Z(4) = AP * (AP * ZUU + CP * ZUV) + CP * (AP * ZUV + CP * ZVV)
                                                                    Z(5) = BP * (BP * ZUU + DP * ZUV) + DP * (BP * ZUV + DP * ZVV)
                                                                    Z(6) = AP * (BP * ZUU + DP * ZUV) + CP * (BP * ZUV + DP * ZVV)

                                                                    RETURN
                                                                END SUBROUTINE C2HORN

                                                            END MODULE ACM_684_C2


                                                            MODULE RPTBDB_common_data
                                                                USE RPTBDB_constant

                                                                integer :: iTriCellOutput = 1
                                                                integer :: RPTBDB_Num_Tri_Cell
                                                                type(p2TriCell), pointer :: RPTBDB_Tri_Cells(:)
                                                                integer, pointer :: RPTBDB_N2N(:), RPTBDB_NEND(:)

                                                            END MODULE RPTBDB_common_data


                                                            MODULE RPTBDB_data_structure
                                                                !
                                                                USE RPTBDB_common_data
                                                                USE ACM_684_C1
                                                                USE ACM_684_C2
                                                                implicit none

                                                            contains

                                                                logical function FTT_CELL_IS_ROOT(cell)
                                                                    Type(FttCell), pointer :: cell

                                                                    if (associated(cell % parent)) then
                                                                        FTT_CELL_IS_ROOT = FALSE
                                                                    else
                                                                        FTT_CELL_IS_ROOT = TRUE
                                                                    endif

                                                                end function FTT_CELL_IS_ROOT

                                                                logical function FTT_CELL_IS_LEAF(cell)
                                                                    Type(FttCell), pointer :: cell

                                                                    if (associated(cell % children)) then
                                                                        FTT_CELL_IS_LEAF = FALSE
                                                                    else
                                                                        FTT_CELL_IS_LEAF = TRUE
                                                                    endif

                                                                end function FTT_CELL_IS_LEAF

                                                                function ftt_cell_level(cell)
                                                                    type(FttCell), pointer :: cell
                                                                    integer :: ftt_cell_level

                                                                    ftt_cell_level = cell % level

                                                                end function ftt_cell_level

                                                                function ftt_cell_id(cell)
                                                                    type(FttCell), pointer :: cell
                                                                    integer :: ftt_cell_id

                                                                    ftt_cell_id = cell % flags

                                                                end function ftt_cell_id

                                                                recursive subroutine cell_traverse_leafs(cell0, max_depth, f)

                                                                    INTERFACE
                                                                        SUBROUTINE f(cell0)
                                                                            USE RPTBDB_constant
                                                                            USE RPTBDB_common_data
                                                                            TYPE(FttCell), pointer :: cell0
                                                                        END SUBROUTINE f
                                                                    END INTERFACE

                                                                    type(FttCell), pointer :: cell0, c
                                                                    integer :: max_depth, n
                                                                    type(FttOct), pointer :: children

                                                                    if (max_depth >= 0 .and. ftt_cell_level(cell0) > max_depth) return

                                                                    if (FTT_CELL_IS_LEAF(cell0)) then
                                                                        call f(cell0)
                                                                    else
                                                                        do n = 1, FTT_CELLS
                                                                            c => cell0 % children % cell(n) % p
                                                                            !     if ( .not. FTT_CELL_IS_DESTROYED (c))
                                                                            call cell_traverse_leafs(c, max_depth, f)
                                                                        enddo
                                                                    endif

                                                                end subroutine cell_traverse_leafs

                                                                recursive function ftt_cell_locate(root, target) result(result_p)

                                                                    type(FttCell), pointer :: root, c, located, result_p
                                                                    type(FttVector) :: target, center
                                                                    integer :: i, n
                                                                    type(FttOct), pointer :: children

                                                                    if (.not.associated(root % children)) then
                                                                        result_p => root
                                                                        return
                                                                    endif

                                                                    children => root % children

                                                                    if (target % xyz(1) <= root % center % xyz(1)) then
                                                                        if (target % xyz(2) <= root % center % xyz(2)) then
                                                                            n = 1
                                                                        else
                                                                            n = 4
                                                                        endif
                                                                    else
                                                                        if (target % xyz(2) <= root % center % xyz(2)) then
                                                                        n = 2
                                                                    else
                                                                        n = 3
                                                                        endif
                                                                    endif

                                                                    c => children % cell(n) % p
                                                                    located => ftt_cell_locate(c, target)
                                                                    if (associated(located)) then
                                                                        result_p => located
                                                                        return
                                                                    endif

                                                                    print *, 'Wrong, too far in ftt_cell_locate'

                                                                end function ftt_cell_locate

                                                                subroutine TriCellNumbering(cell)
                                                                    type(FttCell), pointer :: cell
                                                                    integer :: ib, ie, n, i, j

                                                                    cell % NumTriCells = 2
                                                                    cell % TriGeom = 0
                                                                    do j = 1, FTT_VERTICES
                                                                        if (associated(cell % c2n(FTT_VERTICES + j) % p)) then
                                                                            cell % NumTriCells = cell % NumTriCells + 1
                                                                            cell % TriGeom = cell % TriGeom + (j + FTT_VERTICES) * 10**(FTT_VERTICES - j)
                                                                        endif
                                                                    enddo

                                                                    allocate(cell % tri_cells(cell % NumTriCells))
                                                                    select case (cell % TriGeom)

                                                                    case(0)
                                                                        ib = 1; ie = 2
                                                                    case(5000)
                                                                        ib = 3; ie = 5
                                                                    case(600)
                                                                        ib = 6; ie = 8
                                                                    case(70)
                                                                        ib = 9; ie = 11
                                                                    case(8)
                                                                        ib = 12; ie = 14
                                                                    case(5600)
                                                                        ib = 15; ie = 18
                                                                    case(5070)
                                                                        ib = 19; ie = 22
                                                                    case(5008)
                                                                        ib = 23; ie = 26
                                                                    case(670)
                                                                        ib = 27; ie = 30
                                                                    case(608)
                                                                        ib = 31; ie = 34
                                                                    case(78)
                                                                        ib = 35; ie = 38
                                                                    case(5670)
                                                                        ib = 39; ie = 43
                                                                    case(5608)
                                                                        ib = 44; ie = 48
                                                                    case(5078)
                                                                        ib = 49; ie = 53
                                                                    case(678)
                                                                        ib = 54; ie = 58
                                                                    case(5678)
                                                                        ib = 59; ie = 64
                                                                    case default
                                                                        print *, 'Wrong in Triangulation of FttCells'
                                                                        stop
                                                                    end select

                                                                    if ((ie - ib + 1) /= cell % NumTriCells) then
                                                                        print*, 'Wrong in TriCellNumbering', ib, ie, cell % numTriCells
                                                                        stop
                                                                    endif

                                                                    do i = 1, cell % numTriCells
                                                                        RPTBDB_Num_Tri_Cell = RPTBDB_Num_Tri_Cell + 1
                                                                        allocate(cell % tri_cells(i) % p)
                                                                        cell % tri_cells(i) % p % icheck = 0; cell % tri_cells(i) % p % icheck1 = 0; cell % tri_cells(i) % p % icheck2 = 0
                                                                        RPTBDB_Tri_Cells(RPTBDB_Num_Tri_Cell) % p => cell % tri_cells(i) % p
                                                                        do j = 1, 3
                                                                            cell % tri_cells(i) % p % c2n(j) % p => cell % c2n(tri_cell_vertices(j + 1, ib + i - 1)) % p
                                                                        enddo
                                                                    enddo

                                                                end subroutine TriCellNumbering

                                                                function i2s(i)
                                                                    character(len = 4) :: i2s
                                                                    integer, intent(in) :: i
                                                                    i2s = ''
                                                                    select case(i)
                                                                    case(0:9)
                                                                        write(i2s, '(i1)') i
                                                                    case(10:99)
                                                                        write(i2s, '(i2)') i
                                                                    case default
                                                                        write(i2s, '(4h****)')
                                                                    end select
                                                                end function i2s

                                                                function LinearInterp(sg, xv)

                                                                    type(SplineSeg) :: sg
                                                                    real(rfp) :: LinearInterp

                                                                    integer :: n, i, m, j, k, count, sele
                                                                    real(rfp) :: xv, fact, e, ff, g, h

                                                                    n = size(sg % x)

                                                                    if (xv < sg % x(1) .or. xv > sg % x(n)) then
                                                                        print *, 'Interpolation out of range in CubSplInterp'
                                                                        stop
                                                                    endif

                                                                    count = aint((xv - sg % X(1))/(sg % X(n) - sg % X(1))*(n - 1))
                                                                    count = min(count + 1, n - 1)

                                                                    if (xv < sg % x(count)) then
                                                                        do while (xv < sg % x(count))
                                                                            count = count - 1
                                                                        enddo
                                                                    endif

                                                                    if (xv > sg % x(count + 1)) then
                                                                        do while (xv > sg % x(count + 1))
                                                                            count = count + 1
                                                                        enddo
                                                                    endif

                                                                    fact = sg % X(count + 1) - sg % X(count)

                                                                    LinearInterp = (sg % X(count + 1) - xv) * sg % Y(count) + (xv - sg % X(count)) * sg % Y(count + 1)
                                                                    LinearInterp = LinearInterp / fact

                                                                end function LinearInterp

                                                                function OrgPos(RPTBDB_root, pos)

                                                                    type(FttRootCell), pointer :: RPTBDB_root
                                                                    TYPE(FttVector) :: pos
                                                                    Type(FttVector) :: OrgPos
                                                                    Real(rfp) :: a, x0, y0, xu, xl, yu, yl
                                                                    integer :: m, n

                                                                    m = size(RPTBDB_root % BndCur(1) % x)
                                                                    n = aint(pos % xyz(2) * (m - 1)) + 1
                                                                    a = (n - pos % xyz(2)*(m - 1)) * RPTBDB_root % BndCur(1) % x(n)
                                                                    if (n < m) a = a + (pos % xyz(2)*(m - 1) - n + 1) * RPTBDB_root % BndCur(1) % x(n + 1)
                                                                    !   xl = CubSplInterp(BndCur(1),a)
                                                                    xl = LinearInterp(RPTBDB_root % BndCur(1), a)

                                                                    m = size(RPTBDB_root % BndCur(2) % x)
                                                                    n = aint(pos % xyz(1) * (m - 1)) + 1
                                                                    a = (n - pos % xyz(1)*(m - 1)) * RPTBDB_root % BndCur(2) % x(n)
                                                                    if (n < m) a = a + (pos % xyz(1)*(m - 1) - n + 1) * RPTBDB_root % BndCur(2) % x(n + 1)
                                                                    !   yl = CubSplInterp(BndCur(2),a)
                                                                    yl = LinearInterp(RPTBDB_root % BndCur(2), a)

                                                                    m = size(RPTBDB_root % BndCur(3) % x)
                                                                    n = aint(pos % xyz(2) * (m - 1)) + 1
                                                                    a = (n - pos % xyz(2)*(m - 1)) * RPTBDB_root % BndCur(3) % x(n)
                                                                    if (n < m) a = a + (pos % xyz(2)*(m - 1) - n + 1) * RPTBDB_root % BndCur(3) % x(n + 1)
                                                                    !   xu = CubSplInterp(BndCur(3),a)    
                                                                    xu = LinearInterp(RPTBDB_root % BndCur(3), a)

                                                                    m = size(RPTBDB_root % BndCur(4) % x)
                                                                    n = aint(pos % xyz(1) * (m - 1)) + 1
                                                                    a = (n - pos % xyz(1)*(m - 1)) * RPTBDB_root % BndCur(4) % x(n)
                                                                    if (n < m) a = a + (pos % xyz(1)*(m - 1) - n + 1) * RPTBDB_root % BndCur(4) % x(n + 1)
                                                                    !   yu = CubSplInterp(BndCur(4),a)
                                                                    yu = LinearInterp(RPTBDB_root % BndCur(4), a)

                                                                    if (RPTBDB_root % icst == 1) then
                                                                        OrgPos % xyz(1) = xl + pos % xyz(1) * (xu - xl)
                                                                        OrgPos % xyz(2) = yl + pos % xyz(2) * (yu - yl)
                                                                    else
                                                                        OrgPos % xyz(2) = xl + pos % xyz(1) * (xu - xl)
                                                                        OrgPos % xyz(1) = yl + pos % xyz(2) * (yu - yl)
                                                                    endif

                                                                end function OrgPos


                                                                function CellInterp(RPTBDB_root, cell, pos, org_pos)
                                                                    type(FttRootCell), pointer :: RPTBDB_root
                                                                    type(FttCell), pointer :: cell, p
                                                                    TYPE(FttVector) :: pos
                                                                    TYPE(FttVector), optional :: org_pos
                                                                    Type(Qvar) :: CellInterp, cof_x, cof_y, q_n1, q_n2, q_n3
                                                                    integer :: i, n, n1, n2, n3, icell, id(3)
                                                                    real(rfp) :: xn, yn
                                                                    real :: xyzd(3, 2), zp(6), GD(3), x_org, y_org, rho(3), h(3), zp1(3), zp2(3)
                                                                    real(rfp) :: cof_z, xyz_n1(ndim), xyz_n2(ndim), xyz_n3(ndim)
                                                                    type(TriCell), pointer :: c
                                                                    real(rfp) :: lg10e, T, Pre, G, Gp, Gt, Gpp, Gtt, Gpt

                                                                    lg10e = log10(exp(1.))


                                                                    xn = (pos % xyz(1) - cell % c2n(1) % p % pos % xyz(1)) &
                                                                    / (cell % c2n(2) % p % pos % xyz(1) - cell % c2n(1) % p % pos % xyz(1))
                                                                    yn = (pos % xyz(2) - cell % c2n(1) % p % pos % xyz(2)) &
                                                                    / (cell % c2n(4) % p % pos % xyz(2) - cell % c2n(1) % p % pos % xyz(2))

                                                                    select case (cell % TriGeom)

                                                                    case(0)
                                                                        if (yn <= xn) then
                                                                            icell = 1
                                                                        else
                                                                            icell = 2
                                                                        endif
                                                                    case(5000)
                                                                        if ((xn + 2 * yn) <= 1) then
                                                                            icell = 1
                                                                        else if ((yn - 0.5) <= 0.5 * xn) then
                                                                            icell = 2
                                                                        else
                                                                            icell = 3
                                                                        endif
                                                                    case(600)
                                                                        if ((2 * xn + yn) <= 1) then
                                                                            icell = 1
                                                                        else if (yn <= 2 * (xn - 0.5)) then
                                                                            icell = 2
                                                                        else
                                                                            icell = 3
                                                                        endif
                                                                    case(70)
                                                                        if (yn <= 0.5 * xn) then
                                                                            icell = 1
                                                                        else if ((xn + 2 * (yn - 0.5)) <= 1) then
                                                                            icell = 2
                                                                        else
                                                                            icell = 3
                                                                        endif
                                                                    case(8)
                                                                        if (yn >= 2 * xn) then
                                                                            icell = 2
                                                                        else if ((2 * (xn - 0.5) + yn) <= 1)then
                                                                            icell = 1
                                                                        else
                                                                            icell = 3
                                                                        endif
                                                                    case(5600)
                                                                        if ((xn + yn) <= 0.5) then
                                                                            icell = 1
                                                                        else if (yn <= 2 * (xn - 0.5)) then
                                                                            icell = 2
                                                                        else if ((yn - 0.5) >= 0.5 * xn) then
                                                                            icell = 3
                                                                        else
                                                                            icell = 4
                                                                        endif
                                                                    case(5070)
                                                                        if (yn <= 0.5 * xn) then
                                                                            icell = 1
                                                                        else if (yn <= 0.5) then
                                                                            icell = 2
                                                                        else if ((yn - 0.5) <= 0.5 * xn) then
                                                                            icell = 4
                                                                        else
                                                                            icell = 3
                                                                        endif
                                                                    case(5008)
                                                                        if ((xn + 2 * yn) <= 1) then
                                                                            icell = 1
                                                                        else if ((2 * (xn - 0.5) + yn) >= 1)then
                                                                            icell = 2
                                                                        else if ((yn - 0.5) >= xn) then
                                                                            icell = 4
                                                                        else
                                                                            icell = 3
                                                                        endif
                                                                    case(670)
                                                                        if ((2 * xn + yn) <= 1) then
                                                                            icell = 1
                                                                        else if (yn <= (xn - 0.5)) then
                                                                            icell = 2
                                                                        else if (xn + 2 * (yn - 0.5) >= 1) then
                                                                            icell = 3
                                                                        else
                                                                            icell = 4
                                                                        endif
                                                                    case(608)
                                                                        if (yn >= 2 * xn) then
                                                                            icell = 2
                                                                        else if (xn <= 0.5) then
                                                                            icell = 1
                                                                        else if (yn <= 2 * (xn - 0.5)) then
                                                                            icell = 3
                                                                        else
                                                                            icell = 4
                                                                        endif
                                                                    case(78)
                                                                        if (yn <= 0.5 * xn) then
                                                                            icell = 1
                                                                        else if (yn >= 2 * xn) then
                                                                            icell = 2
                                                                        else if ((xn - 0.5 + yn - 0.5) >= 0.5) then
                                                                            icell = 4
                                                                        else
                                                                            icell = 3
                                                                        endif
                                                                    case(5670)
                                                                        if ((xn + yn) <= 0.5) then
                                                                            icell = 1
                                                                        else if (yn <= (xn - 0.5)) then
                                                                            icell = 2
                                                                        else if (yn <= 0.5) then
                                                                            icell = 5
                                                                        else if ((yn - 0.5) >= 0.5 * xn) then
                                                                            icell = 3
                                                                        else
                                                                            icell = 4
                                                                        endif
                                                                    case(5608)
                                                                        if ((xn + yn) <= 0.5) then
                                                                            icell = 1
                                                                        else if (yn <= 2 * (xn - 0.5)) then
                                                                            icell = 2
                                                                        else if (xn >= 0.5) then
                                                                            icell = 3
                                                                        else if ((yn - 0.5) >= xn) then
                                                                            icell = 4
                                                                        else
                                                                            icell = 5
                                                                        endif
                                                                    case(5078)
                                                                        if (yn <= 0.5 * xn) then
                                                                            icell = 1
                                                                        else if (yn <= 0.5) then
                                                                            icell = 2
                                                                        else if ((xn + yn) >= 1.5) then
                                                                            icell = 3
                                                                        else if ((yn - 0.5) >= xn) then
                                                                            icell = 4
                                                                        else
                                                                            icell = 5
                                                                        endif
                                                                    case(678)
                                                                        if (yn >= 2 * xn) then
                                                                            icell = 2
                                                                        else if (xn <= 0.5) then
                                                                            icell = 1
                                                                        else if (yn <= (xn - 0.5)) then
                                                                            icell = 3
                                                                        else if ((xn + yn) >= 1.5) then
                                                                            icell = 4
                                                                        else
                                                                            icell = 5
                                                                        endif
                                                                    case(5678)
                                                                        if ((xn + yn) <= 0.5) then
                                                                            icell = 1
                                                                        else if (yn <= (xn - 0.5)) then
                                                                            icell = 2
                                                                        else if ((xn + yn) >= 1.5) then
                                                                            icell = 3
                                                                        else if ((yn - 0.5) >= xn) then
                                                                            icell = 4
                                                                        else if (xn <= 0.5) then
                                                                            icell = 5
                                                                        else
                                                                            icell = 6
                                                                        endif
                                                                    case default
                                                                        print *, 'Wrong in Selection of Triangles inside FttCell'
                                                                        stop
                                                                    end select

                                                                    x_org = org_pos % xyz(1); y_org = org_pos % xyz(2)
                                                                    c => cell % tri_cells(icell) % p

                                                                    do i = 1, 3
                                                                        id(i) = c % c2n(i) % p % id
                                                                        xyzd(i,:) = c % c2n(i) % p % org_pos % xyz(:)
                                                                        rho(i) = RPTBDB_root % nodes(id(i)) % q % vars(1)
                                                                        h(i) = RPTBDB_root % nodes(id(i)) % q % vars(2)
                                                                    enddo

                                                                    ! Linear
                                                                    !
                                                                    xyz_n1(:) = RPTBDB_root % nodes(id(1)) % pos % xyz(:)
                                                                    xyz_n2(:) = RPTBDB_root % nodes(id(2)) % pos % xyz(:)
                                                                    xyz_n3(:) = RPTBDB_root % nodes(id(3)) % pos % xyz(:)

                                                                    q_n1 = RPTBDB_root % nodes(id(1)) % q
                                                                    q_n2 = RPTBDB_root % nodes(id(2)) % q
                                                                    q_n3 = RPTBDB_root % nodes(id(3)) % q

                                                                    cof_z = (xyz_n2(1) - xyz_n1(1)) * (xyz_n3(2) - xyz_n1(2)) &
                                                                    -(xyz_n3(1) - xyz_n1(1)) * (xyz_n2(2) - xyz_n1(2))
                                                                    !

                                                                    cof_y % vars(:) = (q_n2 % vars(:) - q_n1 % vars(:)) * (xyz_n3(1) - xyz_n1(1)) &
                                                                    -(q_n3 % vars(:) - q_n1 % vars(:)) * (xyz_n2(1) - xyz_n1(1))
                                                                    !
                                                                    cof_x % vars(:) = (q_n3 % vars(:) - q_n1 % vars(:)) * (xyz_n2(2) - xyz_n1(2)) &
                                                                    -(q_n2 % vars(:) - q_n1 % vars(:)) * (xyz_n3(2) - xyz_n1(2))
                                                                    !     
                                                                    CellInterp % vars(:) = q_n1 % vars(:) - ((pos % xyz(2) - xyz_n1(2)) * cof_y % vars(:) &
                                                                    +(pos % xyz(1) - xyz_n1(1)) * cof_x % vars(:)) / cof_z
                                                                    ! 
                                                                    !  return
                                                                    !  
                                                                    select case(RPTBDB_root % i_con_sel)
                                                                    case(0)
                                                                        return
                                                                    case(1)
                                                                        !   
                                                                        if (c % icheck1 == 0) allocate(c % acf1, c % acf2)

                                                                        CALL C1PNT(RPTBDB_root % nnodes, XYZD, rho, &
                                                                        x_org, y_org, RPTBDB_root % RPTBDB_WK1, &
                                                                        ID, ZP1, c % ICHECK1, &
                                                                        c % acf1 % P0, c % acf1 % P1, c % acf1 % P2, c % acf1 % P3, c % acf1 % P4, c % acf1 % P5, &
                                                                        c % acf1 % P11, c % acf1 % P12, c % acf1 % P13, c % acf1 % P14, c % acf1 % P21, &
                                                                        c % acf1 % P22, c % acf1 % P23, c % acf1 % P31, c % acf1 % P32, c % acf1 % P41, &
                                                                        c % acf1 % DX, c % acf1 % DY, c % acf1 % DX2, c % acf1 % DY2, c % acf1 % SL2, &
                                                                        c % acf1 % AP, c % acf1 % BP, c % acf1 % CP, c % acf1 % DP)

                                                                        CALL C1PNT(RPTBDB_root % nnodes, XYZD, h, &
                                                                        x_org, y_org, RPTBDB_root % RPTBDB_WK2, &
                                                                        ID, ZP2, c % ICHECK2, &
                                                                        c % acf2 % P0, c % acf2 % P1, c % acf2 % P2, c % acf2 % P3, c % acf2 % P4, c % acf2 % P5, &
                                                                        c % acf2 % P11, c % acf2 % P12, c % acf2 % P13, c % acf2 % P14, c % acf2 % P21, &
                                                                        c % acf2 % P22, c % acf2 % P23, c % acf2 % P31, c % acf2 % P32, c % acf2 % P41, &
                                                                        c % acf2 % DX, c % acf2 % DY, c % acf2 % DX2, c % acf2 % DY2, c % acf2 % SL2, &
                                                                        c % acf2 % AP, c % acf2 % BP, c % acf2 % CP, c % acf2 % DP)
                                                                        !
                                                                        T = x_org
                                                                        Pre = 10**y_org * 1000

                                                                        CellInterp % vars(1) = zp1(1)
                                                                        CellInterp % vars(2) = zp2(1)

                                                                        CellInterp % vars(4) = zp1(3) * lg10e/Pre * 1000
                                                                        CellInterp % vars(5) = zp1(2)
                                                                        CellInterp % vars(6) = zp2(3) * lg10e/Pre * 1000
                                                                        CellInterp % vars(7) = zp2(2)

                                                                    case(2)

                                                                        if (c % icheck == 0) allocate(c % acf)

                                                                        do i = 1, 3
                                                                            GD(i) = RPTBDB_root % RPTBDB_GIBBS(id(i))
                                                                        enddo

                                                                        CALL C2PNT(RPTBDB_root % nnodes, XYZD, GD, &
                                                                        x_org, y_org, RPTBDB_root % RPTBDB_WK, &
                                                                        ID, ZP, c % ICHECK, &
                                                                        c % acf % Q00, c % acf % Q01, c % acf % Q02, c % acf % Q03, c % acf % Q04, c % acf % Q05, &
                                                                        c % acf % Q06, c % acf % Q07, c % acf % Q08, c % acf % Q09, &
                                                                        c % acf % Q12, c % acf % Q13, c % acf % Q14, c % acf % Q15, c % acf % Q16, c % acf % Q17, c % acf % Q18, &
                                                                        c % acf % Q23, c % acf % Q24, c % acf % Q25, c % acf % Q26, c % acf % Q27, &
                                                                        c % acf % Q34, c % acf % Q35, c % acf % Q36, &
                                                                        c % acf % Q45, c % acf % Q11, c % acf % Q22, c % acf % Q33, c % acf % Q44, &
                                                                        c % acf % DX, c % acf % DY, c % acf % DX2, c % acf % DY2, c % acf % SL2, &
                                                                        c % acf % AP, c % acf % BP, c % acf % CP, c % acf % DP)
                                                                        !
                                                                        T = x_org
                                                                        Pre = 10**y_org * 1000
                                                                        G = ZP(1)
                                                                        Gt = ZP(2)
                                                                        Gp = ZP(3) * lg10e/Pre
                                                                        Gtt = ZP(4)
                                                                        Gpp = ZP(5) * (lg10e/Pre)**2 + ZP(3) * (-lg10e/Pre**2)
                                                                        Gpt = ZP(6) * lg10e/Pre

                                                                        CellInterp % vars(1) = one / Gp
                                                                        CellInterp % vars(2) = G - T * Gt
                                                                        CellInterp % vars(3) = -Gt
                                                                        CellInterp % vars(4) = -Gpp/Gp**2 * 1000
                                                                        CellInterp % vars(5) = -Gpt/Gp**2
                                                                        CellInterp % vars(6) = (Gp - T * Gpt) * 1000
                                                                        CellInterp % vars(7) = -T * Gtt

                                                                    end select

                                                                end function CellInterp

                                                                function CellInterpPartiallyLinear(cell, pos)
                                                                    type(FttCell), pointer :: cell
                                                                    TYPE(FttVector) :: pos
                                                                    Type(Qvar) :: CellInterpPartiallyLinear, cof_x, cof_y
                                                                    integer :: n1, n2, n3
                                                                    real(rfp) :: cof_z

                                                                    if (SUM(pos % xyz(:)) < SUM(cell % c2n(2) % p % pos % xyz(:))) then
                                                                        n1 = 1; n2 = 2; n3 = 4
                                                                    else
                                                                        n1 = 2; n2 = 3; n3 = 4
                                                                    endif

                                                                    cof_z = (cell % c2n(n2) % p % pos % xyz(1) - cell % c2n(n1) % p % pos % xyz(1)) * &
                                                                    (cell % c2n(n3) % p % pos % xyz(2) - cell % c2n(n1) % p % pos % xyz(2)) &
                                                                    -(cell % c2n(n3) % p % pos % xyz(1) - cell % c2n(n1) % p % pos % xyz(1)) * &
                                                                    (cell % c2n(n2) % p % pos % xyz(2) - cell % c2n(n1) % p % pos % xyz(2))
                                                                    !
                                                                    cof_y % vars(:) = (cell % c2n(n2) % p % q % vars(:) - cell % c2n(n1) % p % q % vars(:)) * &
                                                                    (cell % c2n(n3) % p % pos % xyz(1) - cell % c2n(n1) % p % pos % xyz(1)) &
                                                                    -(cell % c2n(n3) % p % q % vars(:) - cell % c2n(n1) % p % q % vars(:)) * &
                                                                    (cell % c2n(n2) % p % pos % xyz(1) - cell % c2n(n1) % p % pos % xyz(1))
                                                                    !
                                                                    cof_x % vars(:) = (cell % c2n(n3) % p % q % vars(:) - cell % c2n(n1) % p % q % vars(:)) * &
                                                                    (cell % c2n(n2) % p % pos % xyz(2) - cell % c2n(n1) % p % pos % xyz(2)) &
                                                                    -(cell % c2n(n2) % p % q % vars(:) - cell % c2n(n1) % p % q % vars(:)) * &
                                                                    (cell % c2n(n3) % p % pos % xyz(2) - cell % c2n(n1) % p % pos % xyz(2))
                                                                    !     
                                                                    CellInterpPartiallyLinear % vars(:) = cell % c2n(n1) % p % q % vars(:) &
                                                                    -((pos % xyz(2) - cell % c2n(n1) % p % pos % xyz(2)) * cof_y % vars(:) &
                                                                    +(pos % xyz(1) - cell % c2n(n1) % p % pos % xyz(1)) * cof_x % vars(:)) / cof_z

                                                                end function CellInterpPartiallyLinear

                                                            end MODULE RPTBDB_data_structure

                                                            MODULE RPTBDB_spline
                                                                !
                                                                USE RPTBDB_data_structure
                                                                implicit none

                                                            contains

                                                                function CubSplItp(sg, xv) result(yv_out)

                                                                    type(SplineSeg) :: sg
                                                                    real(rfp) :: yv_out

                                                                    integer :: n, i, m, j, k, count, sele
                                                                    real(rfp) :: xv, fact, e, ff, g, h

                                                                    n = size(sg % x)

                                                                    !      if ( xv < sg%x(1) .or. xv > sg%x(n) ) then
                                                                    !        print *, 'Interpolation out of range '
                                                                    !        stop
                                                                    !      endif

                                                                    count = aint((xv - sg % X(1))/(sg % X(n) - sg % X(1))*(n - 1))
                                                                    count = min(count + 1, n - 1)

                                                                    if (xv < sg % x(count)) then
                                                                        do while (xv < sg % x(count))
                                                                            count = count - 1
                                                                        enddo
                                                                    endif

                                                                    if (xv > sg % x(count + 1)) then
                                                                        do while (xv > sg % x(count + 1))
                                                                            count = count + 1
                                                                        enddo
                                                                    endif

                                                                    fact = sg % X(count + 1) - sg % X(count)
                                                                    e = sg % CM(count - 1) / (6. * fact)
                                                                    ff = sg % CM(count) / (6. * fact)
                                                                    g = (sg % Y(count)/fact) - (sg % CM(count - 1) * fact/6.)
                                                                    h = (sg % Y(count + 1)/fact) - (sg % CM(count) * fact/6.)
                                                                    yv_out = e * (sg % X(count + 1) - xv) **3 + ff * (xv - sg % X(count)) **3 + &
                                                                    g * (sg % X(count + 1) - xv) + h * (xv - sg % X(count))

                                                                end function CubSplItp

                                                            end module RPTBDB_spline

                                                            Module RPTBDB_setup
                                                                !
                                                                USE RPTBDB_common_data
                                                                USE RPTBDB_spline
                                                                implicit none

                                                            contains

                                                                function ReadDatabase(dbf) result(RPTBDB_root)

                                                                    type(FttRootCell), pointer :: RPTBDB_root, newRPTBDB_root
                                                                    type(FttCell), pointer :: p, newFttCell
                                                                    integer :: i, j, k, n, idbs, izero, id, ip, ine, ic(FTT_CELLS), iv(FTT_VERTICES + FTT_MAX_NEW_VERTICES), num_bc
                                                                    type(p2FttCell), allocatable :: FttCells(:)
                                                                    character(len = 40), intent(in) :: dbf
                                                                    logical :: yes
                                                                    real(rfp), allocatable :: xn(:), yn(:)

                                                                    nullify(RPTBDB_root)
                                                                    allocate(newRPTBDB_root)
                                                                    RPTBDB_root => newRPTBDB_root

                                                                    allocate(RPTBDB_root % var_nam(10))

                                                                    !  RPTBDB_root%var_nam = (/'T', '1/T', '-1/T', 'P', 'log(P)', 'rho', 'log(rho)', '1/rho', 'h', 's' /)

                                                                    !   RPTBDB_root%var_nam = (/'T [K]            ', '1/T [K]          ', '-1/T [K]         ', 'P [MPa]          ', 'log(P) [MPa]     ' ,  &
                                                                    !                           'rho [kg/m^3]     ', 'log(rho) [kg/m^3]', '1/rho [kg/m^3]   ', 'h [KJ/Kg]        ', 's [KJ/Kg-K]      ' /)

                                                                    !   write(*,*) '/////    Welcome to RPTBDB, Tree-based database for NIST RefProp    \\\\\'
                                                                    !   write(*,*) '/-----------------------------------------------------------------------\'
                                                                    !   write(*,*) '/  RPTBDB provides the database and interpolation module for the        \'
                                                                    !   write(*,*) '/ propertiesof fluid within the region you specified.                   \'
                                                                    !   write(*,*) '/////---------------------------------------------------------------\\\\\'
                                                                    !   write(*,*)

                                                                    izero = 0
                                                                    idbs = 2
                                                                    allocate(RPTBDB_root % bnd_cell(num_of_bnd_cells))
                                                                    num_bc = 0

                                                                    inquire(file = dbf, exist = yes)
                                                                    if (.not.yes) then
                                                                        write(*, *) 'Sorry, could not locate the database file', dbf
                                                                        stop
                                                                    end if

                                                                    open (idbs, file = dbf, status = 'old')

                                                                    read(idbs, *) RPTBDB_root % ix, RPTBDB_root % iy, RPTBDB_root % icst, &
                                                                    RPTBDB_root % if3vars(1:3), RPTBDB_root % i_con_sel
                                                                    read(idbs, *) RPTBDB_root % MolWgt

                                                                    do i = 1, 4
                                                                        n = size(RPTBDB_root % BndCur(i) % X)
                                                                        read(idbs, *) n
                                                                        allocate(RPTBDB_root % BndCur(i) % X(1:n), RPTBDB_root % BndCur(i) % Y(1:n), RPTBDB_root % BndCur(i) % CM(0:n))
                                                                        read(idbs, *) RPTBDB_root % BndCur(i) % X(1:n)
                                                                        read(idbs, *) RPTBDB_root % BndCur(i) % Y(1:n)
                                                                        read(idbs, *) RPTBDB_root % BndCur(i) % CM(0:n)
                                                                    enddo

                                                                    read(idbs, *) RPTBDB_root % nnodes
                                                                    allocate(RPTBDB_root % nodes(RPTBDB_root % nnodes))
                                                                    allocate(RPTBDB_N2N(10 * RPTBDB_root % nnodes), RPTBDB_NEND(RPTBDB_root % nnodes))

                                                                    do i = 1, RPTBDB_root % nnodes
                                                                        read(idbs, *) RPTBDB_root % nodes(i) % pos % xyz(1:ndim)
                                                                        read(idbs, *) RPTBDB_root % nodes(i) % q % vars(1:nvar)
                                                                        RPTBDB_root % nodes(i) % org_pos = OrgPos(RPTBDB_root, RPTBDB_root % nodes(i) % pos)
                                                                        RPTBDB_root % nodes(i) % id = i
                                                                    enddo

                                                                    if (RPTBDB_root % i_con_sel == 1) then
                                                                        allocate(RPTBDB_root % RPTBDB_WK1(RPTBDB_root % nnodes, 5), &
                                                                        RPTBDB_root % RPTBDB_WK2(RPTBDB_root % nnodes, 5))
                                                                        do i = 1, RPTBDB_root % nnodes
                                                                            read(idbs, *) RPTBDB_root % RPTBDB_WK1(i, 1:5)
                                                                            read(idbs, *) RPTBDB_root % RPTBDB_WK2(i, 1:5)
                                                                        enddo
                                                                    endif

                                                                    if (RPTBDB_root % i_con_sel == 2) then
                                                                        allocate(RPTBDB_root % RPTBDB_WK(RPTBDB_root % nnodes, 14), &
                                                                        RPTBDB_root % RPTBDB_GIBBS(RPTBDB_root % nnodes))
                                                                        do i = 1, RPTBDB_root % nnodes
                                                                            read(idbs, *) RPTBDB_root % RPTBDB_WK(i, 1:14)
                                                                        enddo
                                                                    endif

                                                                    read(idbs, *) RPTBDB_root % nFttCells
                                                                    allocate(FttCells(RPTBDB_root % nFttCells))
                                                                    do i = 1, RPTBDB_root % nFttCells
                                                                        nullify(FttCells(i) % p)
                                                                        allocate(newFttCell)
                                                                        nullify(newFttCell % children)
                                                                        FttCells(i) % p => newFttCell
                                                                        FttCells(i) % p % id = i
                                                                    enddo

                                                                    RPTBDB_root % p => FttCells(1) % p

                                                                    do i = 1, RPTBDB_root % nFttCells

                                                                        read(idbs, *) id
                                                                        p => FttCells(id) % p
                                                                        read(idbs, *) p % flags, p % level
                                                                        read(idbs, *) p % center % xyz
                                                                        read(idbs, *) ip
                                                                        if (ip /= izero) then
                                                                            if (.not.associated(FttCells(ip) % p % children)) allocate(FttCells(ip) % p % children)
                                                                            FttCells(ip) % p % children % parent => FttCells(ip) % p
                                                                            FttCells(ip) % p % children % level = p % level - 1
                                                                            p % parent => FttCells(ip) % p % children
                                                                        endif

                                                                        read(idbs, *) (ic(j), j = 1, FTT_CELLS)
                                                                        if (ic(1) /= izero) then
                                                                            if (.not.associated(p % children)) allocate(p % children)
                                                                            do j = 1, FTT_CELLS
                                                                                p % children % cell(j) % p => FttCells(ic(j)) % p
                                                                            enddo
                                                                        endif

                                                                        read(idbs, *) iv(1:FTT_VERTICES + FTT_MAX_NEW_VERTICES)
                                                                        do j = 1, FTT_VERTICES + FTT_MAX_NEW_VERTICES
                                                                            if (iv(j) > 0) then
                                                                                p % c2n(j) % p => RPTBDB_root % nodes(iv(j))
                                                                            else
                                                                                nullify(p % c2n(j) % p)
                                                                            endif
                                                                        enddo

                                                                        allocate(p % neighbors)
                                                                        do j = 1, FTT_NEIGHBORS
                                                                            read(idbs, *) ine
                                                                            if (ine /= izero) then
                                                                                if (ine > 0) then
                                                                                    p % neighbors % e(j) % p => FttCells(ine) % p
                                                                                else
                                                                                    do k = 1, num_bc
                                                                                        if (ine == RPTBDB_root % bnd_cell(k) % id) exit
                                                                                    enddo
                                                                                    if (k > num_bc) then
                                                                                    num_bc = num_bc + 1
                                                                                    RPTBDB_root % bnd_cell(num_bc) % id = ine
                                                                                    p % neighbors % e(j) % p => RPTBDB_root % bnd_cell(num_bc)
                                                                                    nullify(RPTBDB_root % bnd_cell(num_bc) % children, RPTBDB_root % bnd_cell(num_bc) % parent, &
                                                                                    RPTBDB_root % bnd_cell(num_bc) % neighbors)
                                                                                else
                                                                                    p % neighbors % e(j) % p => RPTBDB_root % bnd_cell(k)
                                                                                    endif
                                                                                endif
                                                                            endif
                                                                        enddo

                                                                    enddo

                                                                    allocate(RPTBDB_Tri_Cells(4 * RPTBDB_root % nFttCells))

                                                                    RPTBDB_Num_Tri_Cell = 0
                                                                    call cell_traverse_leafs(RPTBDB_root % p, -1, TriCellNumbering)
                                                                    RPTBDB_root % nTriCells = RPTBDB_Num_Tri_Cell

                                                                    call NodeConnectivity(RPTBDB_root)

                                                                    allocate(xn(RPTBDB_root % nnodes), yn(RPTBDB_root % nnodes))
                                                                    do i = 1, RPTBDB_root % nnodes
                                                                        xn(i) = RPTBDB_root % nodes(i) % org_pos % xyz(1)
                                                                        yn(i) = RPTBDB_root % nodes(i) % org_pos % xyz(2)
                                                                        !
                                                                        !   Here Gibb's function is for log(p) and T only
                                                                        if (RPTBDB_root % i_con_sel == 2) &
                                                                        RPTBDB_root % RPTBDB_GIBBS(i) = RPTBDB_root % nodes(i) % q % vars(2) - &
                                                                        RPTBDB_root % nodes(i) % org_pos % xyz(1) * RPTBDB_root % nodes(i) % q % vars(3)
                                                                        !
                                                                    enddo

                                                                    call OutputTriCellsToTecplot(RPTBDB_root)
                                                                    iTriCellOutput = iTriCellOutput + 1

                                                                    close(idbs)
                                                                    deallocate(FttCells, xn, yn)

                                                                end function ReadDatabase
                                                                !

                                                                subroutine OutputTriCellsToTecplot(RPTBDB_root)

                                                                    type(FttRootCell), pointer :: RPTBDB_root
                                                                    integer :: i, j, id, n
                                                                    type(FttVector), allocatable :: pos1(:)
                                                                    type(FttVector) :: center_pos1
                                                                    integer, allocatable :: num_cell_at_vertex(:)
                                                                    type(Qvar), allocatable :: qnode(:)
                                                                    real(rfp) :: lg10e, T, Pre, G, Gp, Gt, Gpp, Gtt, Gpt
                                                                    real(rfp), allocatable :: rho1(:), h1(:), s1(:), rhop1(:), rhot1(:), hp1(:), ht1(:)

                                                                    n = RPTBDB_root % nnodes
                                                                    allocate(rho1(n), h1(n), s1(n), rhop1(n), rhot1(n), hp1(n), ht1(n))
                                                                    lg10e = log10(exp(1.))

                                                                    open(2, file = 'RPTBDB_TriCells.plt' // i2s(iTriCellOutput))

                                                                    write(2, *) 'variables = '

                                                                    write(2, *) '"', RPTBDB_root % var_nam(RPTBDB_root % ix), '"'
                                                                    write(2, *) '"', RPTBDB_root % var_nam(RPTBDB_root % iy), '"'
                                                                    write(2, *) ' rho, h, s, rhop, rhot, hp, ht,'

                                                                    if (RPTBDB_root % i_con_sel == 1) then
                                                                        write(2, *) 'rhop0, rhot0, hp0, ht0'
                                                                        write(2, *) 'wk11, wk12, wk13, wk14, wk15, wk21, wk22, wk23, wk24 wk25'
                                                                    else if (RPTBDB_root % i_con_sel == 2) then
                                                                        write(2, *) ' rho0, h0, s0, rhop0, rhot0, hp0, ht0'
                                                                        write(2, *) 'wk1, wk2, wk3, wk4, wk5, wk6, wk7, wk8, wk9, wk10, wk11, wk12, wk13, wk14'
                                                                    endif

                                                                    write(2, *) 'zone n=', RPTBDB_root % nnodes, ',e=', RPTBDB_root % nTriCells, ', f=feblock,et=triangle'

                                                                    do i = 1, ndim
                                                                        do j = 1, RPTBDB_root % nnodes
                                                                            write(2, 20) RPTBDB_root % nodes(j) % org_pos % xyz(i)
                                                                        enddo
                                                                    enddo

                                                                    do i = 1, nvar - 2
                                                                        do j = 1, RPTBDB_root % nnodes
                                                                            write(2, 20) RPTBDB_root % nodes(j) % q % vars(i)
                                                                        enddo
                                                                    enddo

                                                                    if (RPTBDB_root % i_con_sel == 2) then

                                                                        do j = 1, RPTBDB_root % nnodes
                                                                            T = RPTBDB_root % nodes(j) % org_pos % xyz(1)
                                                                            Pre = 10**RPTBDB_root % nodes(j) % org_pos % xyz(2) * 1000
                                                                            G = RPTBDB_root % RPTBDB_GIBBS(j)
                                                                            Gt = RPTBDB_root % RPTBDB_WK(j, 1)
                                                                            Gp = RPTBDB_root % RPTBDB_WK(j, 2) * lg10e/Pre
                                                                            Gtt = RPTBDB_root % RPTBDB_WK(j, 3)
                                                                            Gpp = RPTBDB_root % RPTBDB_WK(j, 5) * (lg10e/Pre)**2 + RPTBDB_root % RPTBDB_WK(j, 2) * (-lg10e/Pre**2)
                                                                            Gpt = RPTBDB_root % RPTBDB_WK(j, 4) * lg10e/Pre

                                                                            rho1(j) = one / Gp
                                                                            h1(j) = G - T * Gt
                                                                            s1(j) = -Gt
                                                                            rhop1(j) = -Gpp/Gp**2 * 1000
                                                                            rhot1(j) = -Gpt/Gp**2
                                                                            hp1(j) = (Gp - T * Gpt) * 1000
                                                                            ht1(j) = -T * Gtt
                                                                        enddo

                                                                        write(2, 20) rho1
                                                                        write(2, 20) h1
                                                                        write(2, 20) s1
                                                                        write(2, 20) rhop1
                                                                        write(2, 20) rhot1
                                                                        write(2, 20) hp1
                                                                        write(2, 20) ht1
                                                                        do i = 1, 14
                                                                            write(2, 20) RPTBDB_root % RPTBDB_WK(:, i)
                                                                        enddo

                                                                    else if (RPTBDB_root % i_con_sel == 1) then

                                                                        do j = 1, RPTBDB_root % nnodes
                                                                            T = RPTBDB_root % nodes(j) % org_pos % xyz(1)
                                                                            Pre = 10**RPTBDB_root % nodes(j) % org_pos % xyz(2) * 1000

                                                                            rhop1(j) = RPTBDB_root % RPTBDB_WK1(j, 2) * lg10e/Pre * 1000
                                                                            rhot1(j) = RPTBDB_root % RPTBDB_WK1(j, 1)
                                                                            hp1(j) = RPTBDB_root % RPTBDB_WK2(j, 2) * lg10e/Pre * 1000
                                                                            ht1(j) = RPTBDB_root % RPTBDB_WK2(j, 1)
                                                                        enddo

                                                                        write(2, 20) rhop1
                                                                        write(2, 20) rhot1
                                                                        write(2, 20) hp1
                                                                        write(2, 20) ht1

                                                                        do i = 1, 5
                                                                            write(2, 20) RPTBDB_root % RPTBDB_WK1(:, i)
                                                                        enddo

                                                                        do i = 1, 5
                                                                            write(2, 20) RPTBDB_root % RPTBDB_WK2(:, i)
                                                                        enddo

                                                                    endif

                                                                    do i = 1, RPTBDB_root % nTriCells
                                                                        do j = 1, 3
                                                                            write(2, 10) RPTBDB_Tri_Cells(i) % p % c2n(j) % p % id
                                                                        enddo
                                                                    enddo

                                                                    10 format(10I8)
                                                                    20 format(5e20.8)

                                                                    close(2)
                                                                    deallocate(rho1, h1, s1, rhop1, rhot1, hp1, ht1)
                                                                    RPTBDB_num_Tri_Cell = 0
                                                                    deallocate(RPTBDB_NEND, RPTBDB_N2N, RPTBDB_Tri_Cells)

                                                                end subroutine OutputTriCellsToTecplot


                                                                subroutine NodeConnectivity(RPTBDB_root)
                                                                    type(FttRootCell), pointer :: RPTBDB_root
                                                                    type(RPTBDB_node), pointer :: p0, p1
                                                                    type(FttVector) :: pos0, pos1, posb
                                                                    integer :: i, j, k, ib, jcheck(14), jmax, iend
                                                                    integer, allocatable :: n2n(:,:), nsrn(:)
                                                                    real(rfp) :: alpha0, alpha1, alpha_max

                                                                    !    if ( RPTBDB_root%i_con_sel < 2 )  return

                                                                    allocate(n2n(RPTBDB_root % nnodes, 12), nsrn(RPTBDB_root % nnodes))

                                                                    nsrn(:) = 0
                                                                    n2n(:,:) = 0

                                                                    do i = 1, RPTBDB_Num_Tri_Cell
                                                                        do j = 1, 3
                                                                            k = j + 1
                                                                            if (k > 3) k = k - 3
                                                                            call AddNeighNode(RPTBDB_root, RPTBDB_Tri_Cells(i) % p % c2n(j) % p, &
                                                                            RPTBDB_Tri_Cells(i) % p % c2n(k) % p, n2n, nsrn)
                                                                            k = j + 2
                                                                            if (k > 3) k = k - 3
                                                                            call AddNeighNode(RPTBDB_root, RPTBDB_Tri_Cells(i) % p % c2n(j) % p, &
                                                                            RPTBDB_Tri_Cells(i) % p % c2n(k) % p, n2n, nsrn)
                                                                        enddo
                                                                    enddo

                                                                    do i = 1, RPTBDB_root % nnodes
                                                                        ib = 0
                                                                        p0 => RPTBDB_root % nodes(i)
                                                                        pos0 = p0 % pos
                                                                        posb = pos0
                                                                        do j = 1, ndim
                                                                            if (pos0 % xyz(j) > 0.99999) then
                                                                                posb % xyz(j) = posb % xyz(j) + 0.2
                                                                                ib = ib + 1
                                                                            endif
                                                                            if (pos0 % xyz(j) < 0.00001) then
                                                                                posb % xyz(j) = posb % xyz(j) - 0.2
                                                                                ib = ib + 1
                                                                            endif
                                                                        enddo

                                                                        if (ib > 0) then
                                                                            nsrn(i) = nsrn(i) + 1
                                                                            n2n(i, nsrn(i)) = -1
                                                                        else
                                                                            posb = RPTBDB_root % nodes(n2n(i, nsrn(i))) % pos
                                                                        endif

                                                                        iend = nsrn(i)
                                                                        if (i > 1) iend = iend + RPTBDB_NEND(i - 1)
                                                                        if (iend > 10 * RPTBDB_root % nnodes) then
                                                                            print *, 'Not enough space allocated for N2N', iend, 10 * RPTBDB_root % nnodes
                                                                            stop
                                                                        end if
                                                                        RPTBDB_NEND(i) = iend
                                                                        RPTBDB_N2N(iend) = max(n2n(i, nsrn(i)), 0)

                                                                        jcheck(:) = 0
                                                                        jmax = 1
                                                                        alpha0 = atan2(posb % xyz(2) - pos0 % xyz(2), posb % xyz(1) - pos0 % xyz(1))
                                                                        do k = 1, nsrn(i) - 1
                                                                            alpha_max = -5 * pi
                                                                            do j = 1, nsrn(i) - 1
                                                                                if (jcheck(j) /= 1) then
                                                                                    pos1 = RPTBDB_root % nodes(n2n(i, j)) % pos
                                                                                    alpha1 = atan2(pos1 % xyz(2) - pos0 % xyz(2), pos1 % xyz(1) - pos0 % xyz(1))
                                                                                    if (alpha1 < alpha0) alpha1 = alpha1 + pi * 2
                                                                                    if (alpha1 > alpha_max) then
                                                                                        jmax = j
                                                                                        alpha_max = alpha1
                                                                                    endif
                                                                                endif
                                                                            enddo
                                                                            jcheck(jmax) = 1
                                                                            RPTBDB_N2N(iend - k) = n2n(i, jmax)
                                                                        enddo

                                                                    enddo

                                                                    !    write(*,*) 'Triangulation with ', iend, ' total node connectivity, ', &
                                                                    !                10*RPTBDB_root%nnodes, ' spaces allocated'

                                                                    deallocate(n2n, nsrn)

                                                                end subroutine NodeConnectivity

                                                                subroutine AddNeighNode(RPTBDB_root, p0, p1, n2n, nsrn)
                                                                    type(FttRootCell), pointer :: RPTBDB_root
                                                                    integer :: n2n(:,:), nsrn(:)
                                                                    integer :: n1, n0, j
                                                                    type(RPTBDB_node), pointer :: p0, p1

                                                                    n0 = p0 % id
                                                                    n1 = p1 % id

                                                                    do j = 1, nsrn(n0)
                                                                        if (n2n(n0, j) == n1) return
                                                                    enddo

                                                                    nsrn(n0) = nsrn(n0) + 1
                                                                    n2n(n0, nsrn(n0)) = n1

                                                                end subroutine AddNeighNode

                                                                function QLocate(p, x, y, ierr)

                                                                    type(Qvar) :: QLocate
                                                                    type(FttRootCell), pointer :: p
                                                                    type(FttCell), pointer :: p1, found
                                                                    type(FttVector) :: targetp, org_pos
                                                                    integer :: i, n, m, ierr
                                                                    real(rfp) :: x, y, x0, y0, xp, yp, xu, xl, xe, ye

                                                                    ierr = 0

                                                                    if (p % icst == 1) then
                                                                        x0 = x
                                                                        y0 = y
                                                                    else
                                                                        x0 = y
                                                                        y0 = x
                                                                    endif

                                                                    ye = y
                                                                    if (y0 < p % BndCur(2) % Y(1) .or. y0 > p % BndCur(4) % Y(1)) then

                                                                        ye = max(p % BndCur(4) % Y(1), min(y0, p % BndCur(2) % Y(1)))
                                                                        ierr = 1
                                                                        !    return

                                                                    endif

                                                                    xl = CubSplItp(p % BndCur(1), ye)

                                                                    xu = CubSplItp(p % BndCur(3), ye)

                                                                    if ((p % icst == 1 .and. (x0 < xl .or. x0 > xu)) .or. &
                                                                        (p % icst == 2 .and. (x0 > xl .or. x0 < xu))) then
                                                                        ierr = 1
                                                                        !    return
                                                                    endif

                                                                    xp = (x0 - xl)/(xu - xl)
                                                                    yp = (y0 - p % BndCur(2) % Y(1))/(p % BndCur(4) % Y(1) - p % BndCur(2) % Y(1))

                                                                    targetp % xyz(1) = min(max(xp, 0.), 1.)
                                                                    targetp % xyz(2) = min(max(yp, 0.), 1.)
                                                                    !   org_pos = OrgPos(p,targetp)
                                                                    org_pos % xyz(1) = x0
                                                                    org_pos % xyz(2) = y0

                                                                    p1 => p % p
                                                                    found => ftt_cell_locate(p1, targetp)

                                                                    targetp % xyz(1) = xp
                                                                    targetp % xyz(2) = yp
                                                                    QLocate = CellInterp(p, found, targetp, org_pos)
                                                                    !   QLocate = CellInterpPartiallyLinear(found,targetp)

                                                                end function QLocate


                                                            END module RPTBDB_setup
