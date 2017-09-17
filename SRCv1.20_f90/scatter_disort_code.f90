! this only differs by "SUBROUTINE GETMOM" from Dave Turners LBLDIS
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:44
 
! which is in scatter_disort_misc.f
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! the same as the f77 version picked up from Warren Wiscombe ftp site on
! climate.gsfc.nasa.gov in August 20000, except for one or two modifications
! indicated by cccccccccccc sergio
! and every REAL variable --> DOUBLE PRECISION AFTER compilation
! this means that R1MACH has been replaced by D1MACH

! functions/subroutines isamx, sasum, sscal, sswap,saxpy, sdot have been
! renamed to isamx_dis, sasum_dis, sscal_dis, sswap_dis,saxpy_dis, sdot_dis
! so that there is nbo conflict with the BLAS library

! rewrote plkavg so that it does not do any integrating, but just computes
! Plank radiance (if you use the integration feature, divide out by
! delta(wavenumber)
! ie DISORT plkavg ---> plkavg_orig
!           plkavg ---> ttorad(w,T)

! also rewrote write( * ----> write( kStdWarn
! this meant that a few fcns/subroutine which needed include 'scatterparam.f90'
! have to be modified so that maxmom, maxulv are not sent in as arguments
! and reused inside of scatterparam.f90

! by Frank Evans sugggestion, we turn DELTAM = .FALSE.
! First, you don't need to un-delta scale.  Just run DISORT with its
! delta scaling flag off, but using the delta scaled sscatmie output.
! This would match rtspec better.


! RCS version control information:
! $Header: DISORTsp.f,v 2.1 2000/04/04 18:21:55 laszlo Exp $
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      SUBROUTINE DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,
!     &                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
!     &                   UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
!     &                   FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS,
!     &                   PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY,
!     &                   MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, RFLDN,
!     &                   FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )

! since they are in  include '../INCLUDE/scatterparam.f90'
! got rid of MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM
!      INTEGER   MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI,
!     &          MXSQT
!      PARAMETER ( MXCLY = maxcly, MXULV = maxcly + 1, MXCMU = maxcmu,
!     &            MXUMU = maxumu, MXPHI = maxphi,
!     &            MI = MXCMU / 2, MI9M2 = 9*MI - 2,
!     &            NNLYRI = MXCMU*MXCLY, MXSQT = maxsqt )

SUBROUTINE DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,  &
    WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,  &
    UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,  &
FISOT, LAMBER, ALBEDO  BTEMP, TTEMP, TEMIS,  &
      PLANK, ONLYFL, ACCUR, PRNT, HEADER, RFLDIR, RFLDN,  &
      FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )
  
  
! *******************************************************************
!       Plane-parallel discrete ordinates radiative transfer program
!             ( see DISORT.doc for complete documentation )
! *******************************************************************
  
! +------------------------------------------------------------------+
!  Calling Tree (omitting calls to ERRMSG):
!  (routines in parentheses are not in this file)
  
!  DISORT-+-(D1MACH)
!         +-SLFTST-+-(TSTBAD)
!         +-ZEROIT
!         +-CHEKIN-+-(WRTBAD)
!         |        +-(WRTDIM)
!         |        +-DREF
!         +-ZEROAL
!         +-SETDIS-+-QGAUSN-+-(D1MACH)
!         +-PRTINP
!         +-ALBTRN-+-LEPOLY
!         |        +-ZEROIT
!         |        +-SOLEIG-+-ASYMTX-+-(D1MACH)
!         |        +-TERPEV
!         |        +-SETMTX-+-ZEROIT
!         |        +-(DGBCO)
!         |        +-SOLVE1-+-ZEROIT
!         |        |        +-(DGBSL)
!         |        +-ALTRIN
!         |        +-SPALTR
!         |        +-PRALTR
!         +-PLKAVG-+-(D1MACH)
!         +-LEPOLY
!         +-SURFAC-+-QGAUSN-+-(D1MACH)
!         |        +-BDREF
!         |        +-ZEROIT
!         +-SOLEIG-+-ASYMTX-+-(D1MACH)
!         +-UPBEAM-+-(DGECO)
!         |        +-(DGESL)
!         +-UPISOT-+-(DGECO)
!         |        +-(DGESL)
!         +-TERPEV
!         +-TERPSO
!         +-SETMTX-+-ZEROIT
!         +-SOLVE0-+-ZEROIT
!         |        +-(DGBCO)
!         |        +-(DGBSL)
!         +-FLUXES--ZEROIT
!         +-ZEROIT
!         +-USRINT
!         +-CMPINT
!         +-PRAVIN
!         +-ZEROIT
!         +-RATIO--(D1MACH)
!         +-INTCOR-+-SINSCA
!         |        +-SECSCA-+-XIFUNC
!         +-PRTINT
  
! *** Intrinsic Functions used in DISORT package which take
!     non-negligible amount of time:
  
!    EXP :  Called by- ALBTRN, ALTRIN, CMPINT, FLUXES, SETDIS,
!                      SETMTX, SPALTR, USRINT, PLKAVG
  
!    SQRT : Called by- ASYMTX, SOLEIG
  
! +-------------------------------------------------------------------+
  
!  Index conventions (for all DO-loops and all variable descriptions):
  
!     IU     :  for user polar angles
  
!  IQ,JQ,KQ  :  for computational polar angles ('quadrature angles')
  
!   IQ/2     :  for half the computational polar angles (just the ones
!               in either 0-90 degrees, or 90-180 degrees)
  
!     J      :  for user azimuthal angles
  
!     K,L    :  for Legendre expansion coefficients or, alternatively,
!               subscripts of associated Legendre polynomials
  
!     LU     :  for user levels
  
!     LC     :  for computational layers (each having a different
!               single-scatter albedo and/or phase function)
  
!    LEV     :  for computational levels
  
!    MAZIM   :  for azimuthal components in Fourier cosine expansion
!               of intensity and phase function
  
! +------------------------------------------------------------------+
  
!               I N T E R N A L    V A R I A B L E S
  
!   AMB(IQ/2,IQ/2)    First matrix factor in reduced eigenvalue problem
!                     of Eqs. SS(12), STWJ(8E), STWL(23f)
!                     (used only in SOLEIG)
  
!   APB(IQ/2,IQ/2)    Second matrix factor in reduced eigenvalue problem
!                     of Eqs. SS(12), STWJ(8E), STWL(23f)
!                     (used only in SOLEIG)
  
!   ARRAY(IQ,IQ)      Scratch matrix for SOLEIG, UPBEAM and UPISOT
!                     (see each subroutine for definition)
  
!   B()               Right-hand side vector of Eq. SC(5) going into
!                     SOLVE0,1;  returns as solution vector
!                     vector  L, the constants of integration
  
!   BDR(IQ/2,0:IQ/2)  Bottom-boundary bidirectional reflectivity for a
!                     given azimuthal component.  First index always
!                     refers to a computational angle.  Second index:
!                     if zero, refers to incident beam angle UMU0;
!                     if non-zero, refers to a computational angle.
  
!   BEM(IQ/2)         Bottom-boundary directional emissivity at compu-
!                     tational angles.
  
!   BPLANK            Intensity emitted from bottom boundary
  
!   CBAND()           Matrix of left-hand side of the linear system
!                     Eq. SC(5), scaled by Eq. SC(12);  in banded
!                     form required by LINPACK solution routines
  
!   CC(IQ,IQ)         C-sub-IJ in Eq. SS(5)
  
!   CMU(IQ)           Computational polar angles (Gaussian)
  
!   CWT(IQ)           Quadrature weights corresponding to CMU
  
!   CORINT            When set TRUE, correct intensities for
!                     delta-scaling effects (see Nakajima and Tanaka,
!                     1988). When FALSE, intensities are not corrected.
!                     In general, CORINT should be set true when beam
!                     source is present (FBEAM is not zero) and DELTAM
!                     is TRUE in a problem including scattering.
!                     However, execution is faster when CORINT is FALSE,
!                     and intensities outside the aureole may still be
!                     accurate enough.  When CORINT is TRUE, it is
!                     important to have a sufficiently high order of
!                     Legendre approximation of the phase function. This
!                     is because the intensities are corrected by
!                     calculating the single-scattered radiation, for
!                     which an adequate representation of the phase
!                     function is crucial.  In case of a low order
!                     Legendre approximation of an otherwise highly
!                     anisotropic phase function, the intensities might
!                     actually be more accurate when CORINT is FALSE.
!                     When only fluxes are calculated (ONLYFL is TRUE),
!                     or there is no beam source (FBEAM=0.0), or there
!                     is no scattering (SSALB=0.0 for all layers) CORINT
!                     is set FALSE by the code.
  
!   DELM0             Kronecker delta, delta-sub-M0, where M = MAZIM
!                     is the number of the Fourier component in the
!                     azimuth cosine expansion
  
!   DELTAM            TRUE,  use delta-M method ( see Wiscombe, 1977 );
!                     FALSE, do not use delta-M method. In general, for
!                     a given number of streams, intensities and
!                     fluxes will be more accurate for phase functions
!                     with a large forward peak if DELTAM is set true.
!                     Intensities close to the forward scattering
!                     direction are often less accurate, however, when
!                     the delta-M method is applied. The intensity
!                     correction of Nakajima and Tanaka is used to
!                     improve the accuracy of the intensities.
  
!   DITHER            Small quantity subtracted from single-scattering
!                     albedos of unity, in order to avoid using special
!                     case formulas;  prevents an eigenvalue of exactly
!                     zero from occurring, which would cause an
!                     immediate overflow
  
!   DTAUCP(LC)        Computational-layer optical depths (delta-M-scaled
!                     if DELTAM = TRUE, otherwise equal to DTAUC)
  
!   EMU(IU)           Bottom-boundary directional emissivity at user
!                     angles.
  
!   EVAL(IQ)          Temporary storage for eigenvalues of Eq. SS(12)
  
!   EVECC(IQ,IQ)      Complete eigenvectors of SS(7) on return from
!                     SOLEIG; stored permanently in  GC
  
!   EXPBEA(LC)        Transmission of direct beam in delta-M optical
!                     depth coordinates
  
!   FLYR(LC)          Separated fraction in delta-M method
  
!   GL(K,LC)          Phase function Legendre polynomial expansion
!                     coefficients, calculated from PMOM by
!                     including single-scattering albedo, factor
!                     2K+1, and (if DELTAM=TRUE) the delta-M
!                     scaling
  
!   GC(IQ,IQ,LC)      Eigenvectors at polar quadrature angles,
!                     g  in Eq. SC(1)
  
!   GU(IU,IQ,LC)      Eigenvectors interpolated to user polar angles
!                     ( g  in Eqs. SC(3) and S1(8-9), i.e.
!                       G without the L factor )
  
!   IPVT(LC*IQ)       Integer vector of pivot indices for LINPACK
!                     routines
  
!   KK(IQ,LC)         Eigenvalues of coeff. matrix in Eq. SS(7)
  
!   KCONV             Counter in azimuth convergence test
  
!   LAYRU(LU)         Computational layer in which user output level
!                     UTAU(LU) is located
  
!   LL(IQ,LC)         Constants of integration L in Eq. SC(1),
!                     obtained by solving scaled version of Eq. SC(5)
  
!   LYRCUT            TRUE, radiation is assumed zero below layer
!                     NCUT because of almost complete absorption
  
!   NAZ               Number of azimuthal components considered
  
!   NCUT              Computational layer number in which absorption
!                     optical depth first exceeds ABSCUT
  
!   OPRIM(LC)         Single scattering albedo after delta-M scaling
  
!   PASS1             TRUE on first entry, FALSE thereafter
  
!   PKAG(0:LC)        Integrated Planck function for internal emission
  
!   PRNTU0(L)         logical flag to trigger printing of azimuthally-
!                     averaged intensities:
!                       L    quantities printed
!                      --    ------------------
!                       1    azimuthally-averaged intensities at user
!                               levels and computational polar angles
!                       2    azimuthally-averaged intensities at user
!                               levels and user polar angles
  
!   PSI0(IQ)          Sum just after square bracket in  Eq. SD(9)
  
!   PSI1(IQ)          Sum in  Eq. STWL(31d)
  
!   RMU(IU,0:IQ)      Bottom-boundary bidirectional reflectivity for a
!                     given azimuthal component.  First index always
!                     refers to a user angle.  Second index:
!                     if zero, refers to incident beam angle UMU0;
!                     if non-zero, refers to a computational angle.
  
!   SQT(k)            Square root of k (used only in LEPOLY for
!                     computing associated Legendre polynomials)
  
!   TAUC(0:LC)        Cumulative optical depth (un-delta-M-scaled)
  
!   TAUCPR(0:LC)      Cumulative optical depth (delta-M-scaled if
!                     DELTAM = TRUE, otherwise equal to TAUC)
  
!   TPLANK            Intensity emitted from top boundary
  
!   UUM(IU,LU)        Expansion coefficients when the intensity
!                     (u-super-M) is expanded in Fourier cosine series
!                     in azimuth angle
  
!   U0C(IQ,LU)        Azimuthally-averaged intensity at quadrature
!                     angle
  
!   U0U(IU,LU)        If ONLYFL = FALSE, azimuthally-averaged intensity
!                     at user angles and user levels
  
!                     If ONLYFL = TRUE and MAXUMU.GE.NSTR,
!                     azimuthally-averaged intensity at computational
!                     (Gaussian quadrature) angles and user levels;
!                     the corresponding quadrature angle cosines are
!                     returned in UMU.  If MAXUMU.LT.NSTR, U0U will be
!                     zeroed, and UMU, NUMU will not be set.
  
!   UTAUPR(LU)        Optical depths of user output levels in delta-M
!                     coordinates;  equal to  UTAU(LU) if no delta-M
  
!   WK()              scratch array
  
!   XR0(LC)           X-sub-zero in expansion of thermal source func-
!                     tion preceding Eq. SS(14)(has no mu-dependence);
!                     b-sub-zero in Eq. STWL(24d)
  
!   XR1(LC)           X-sub-one in expansion of thermal source func-
!                     tion; see  Eqs. SS(14-16); b-sub-one in STWL(24d)
  
!   YLM0(L)           Normalized associated Legendre polynomial
!                     of subscript L at the beam angle (not saved
!                     as function of superscipt M)
  
!   YLMC(L,IQ)        Normalized associated Legendre polynomial
!                     of subscript L at the computational angles
!                     (not saved as function of superscipt M)
  
!   YLMU(L,IU)        Normalized associated Legendre polynomial
!                     of subscript L at the user angles
!                     (not saved as function of superscipt M)
  
!   Z()               scratch array used in SOLVE0, ALBTRN to solve
!                     a linear system for the constants of integration
  
!   Z0(IQ)            Solution vectors Z-sub-zero of Eq. SS(16)
  
!   Z0U(IU,LC)        Z-sub-zero in Eq. SS(16) interpolated to user
!                     angles from an equation derived from SS(16)
  
!   Z1(IQ)            Solution vectors Z-sub-one  of Eq. SS(16)
  
!   Z1U(IU,LC)        Z-sub-one in Eq. SS(16) interpolated to user
!                     angles from an equation derived from SS(16)
  
!   ZBEAM(IU,LC)      Particular solution for beam source
  
!   ZJ(IQ)            Right-hand side vector  X-sub-zero in
!                     Eq. SS(19), also the solution vector
!                     Z-sub-zero after solving that system
  
!   ZZ(IQ,LC)         Permanent storage for the beam source vectors ZJ
  
!   ZPLK0(IQ,LC)      Permanent storage for the thermal source
!                     vectors  Z0  obtained by solving  Eq. SS(16)
  
!   ZPLK1(IQ,LC)      Permanent storage for the thermal source
!                     vectors  Z1  obtained by solving  Eq. SS(16)
  
! +-------------------------------------------------------------------+
  
!  LOCAL SYMBOLIC DIMENSIONS (have big effect on storage requirements):
  
!       MXCLY  = Max no. of computational layers
!       MXULV  = Max no. of output levels
!       MXCMU  = Max no. of computation polar angles
!       MXUMU  = Max no. of output polar angles
!       MXPHI  = Max no. of output azimuthal angles
!       MXSQT  = Max no. of square roots of integers (for LEPOLY)
! +-------------------------------------------------------------------+
  
!     .. Parameters ..
  
!ccccccccc sergio modified this a little
  
  INTEGER, INTENT(IN)                      :: NLYR
  REAL, INTENT(IN)                         :: DTAUC( MAXCLY )
  REAL, INTENT(IN OUT)                     :: SSALB( MAXCLY )
  INTEGER, INTENT(IN OUT)                  :: NMOM
  REAL, INTENT(IN OUT)                     :: PMOM( 0:MAXMOM, MAXCLY )
  REAL, INTENT(IN)                         :: TEMPER( 0:MAXCLY )
  REAL, INTENT(IN)                         :: WVNMLO
  REAL, INTENT(IN)                         :: WVNMHI
  LOGICAL, INTENT(IN OUT)                  :: USRTAU
  INTEGER, INTENT(IN)                      :: NTAU
  REAL, INTENT(IN OUT)                     :: UTAU( MAXULV )
  INTEGER, INTENT(IN)                      :: NSTR
  LOGICAL, INTENT(IN)                      :: USRANG
  INTEGER, INTENT(IN)                      :: NUMU
  REAL, INTENT(IN OUT)                     :: UMU( MAXUMU )
  INTEGER, INTENT(IN)                      :: NPHI
  REAL, INTENT(IN)                         :: PHI( MAXPHI )
  INTEGER, INTENT(IN)                      :: IBCND
  REAL, INTENT(IN)                         :: FBEAM
  REAL, INTENT(IN)                         :: UMU0
  REAL, INTENT(IN)                         :: PHI0
  REAL, INTENT(IN OUT)                     :: FISOT
  LOGICAL, INTENT(IN OUT)                  :: LAMBER
  NO TYPE, INTENT(IN OUT)                  :: ALBEDO  BT
    REAL, INTENT(IN)                         :: TTEMP
    REAL, INTENT(IN)                         :: TEMIS
    LOGICAL, INTENT(IN)                      :: PLANK
    LOGICAL, INTENT(IN)                      :: ONLYFL
    REAL, INTENT(IN)                         :: ACCUR
    LOGICAL, INTENT(IN)                      :: PRNT( 5 )
    CHARACTER (LEN=127), INTENT(IN OUT)      :: HEADER
    REAL, INTENT(IN OUT)                     :: RFLDIR( MAXULV )
    REAL, INTENT(IN OUT)                     :: RFLDN( MAXULV )
    REAL, INTENT(IN OUT)                     :: FLUP( MAXULV )
    REAL, INTENT(IN OUT)                     :: DFDT( MAXULV )
    REAL, INTENT(IN OUT)                     :: UAVG( MAXULV )
    REAL, INTENT(OUT)                        :: UU( MAXUMU, MAXULV, MAXPHI
    REAL, INTENT(IN OUT)                     :: ALBMED( MAXUMU )
    REAL, INTENT(IN OUT)                     :: TRNMED( MAXUMU )
    IMPLICIT NONE
    
    INCLUDE '../INCLUDE/scatterparam.f90'
    
    INTEGER ::
    INTEGER, PARAMETER :: MXCLY = maxcly
    INTEGER, PARAMETER :: MXULV = maxcly + 1
    INTEGER, PARAMETER :: MXCMU = maxcmu
    INTEGER, PARAMETER :: MXUMU = maxumu
    INTEGER, PARAMETER :: MXPHI = maxphi
    INTEGER, PARAMETER :: MI = MXCMU / 2
    INTEGER, PARAMETER :: MI9M2 = 9*MI - 2
    INTEGER, PARAMETER :: NNLYRI = MXCMU*MXCLY
    INTEGER, PARAMETER :: MXSQT = maxsqt
!     ..
!     .. Scalar Arguments ..
    
! got rid of MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM by using scatterparam.f90
    
    
    
!cccccccc sergio
!ccccccc      INTEGER   IBCND, MAXCLY, MAXMOM, MAXPHI, MAXULV, MAXUMU, NLYR,
!ccccccc    &          NMOM, NPHI, NSTR, NTAU, NUMU
    
    REAL :: ALBEDO  BTEMP,  &
          
!     ..
!     .. Array Arguments ..
      
      
      
!     ..
!     .. Local Scalars ..
      
      LOGICAL :: COMPAR, CORINT, DELTAM, LYRCUT, PASS1
      INTEGER :: IQ, IU, J, KCONV, L, LC, LEV, LU, MAZIM, NAZ, NCOL,  &
          NCOS, NCUT, NN, NS
      REAL :: ANGCOS, AZERR, AZTERM, BPLANK, COSPHI, DELM0,  &
          DITHER, DUM, PI, RPD, SGN, TPLANK
!     ..
!     .. Local Arrays ..
      
      LOGICAL :: PRNTU0( 2 )
      INTEGER :: IPVT( NNLYRI ), LAYRU( MXULV )
      
      REAL :: AMB( MI, MI ), APB( MI, MI ), ARRAY( MXCMU, MXCMU ),  &
          B( NNLYRI ), BDR( MI, 0:MI ), BEM( MI ),  &
          CBAND( MI9M2, NNLYRI ), CC( MXCMU, MXCMU ),  &
          CMU( MXCMU ), CWT( MXCMU ), DTAUCP( MXCLY ),  &
          EMU( MXUMU ), EVAL( MI ), EVECC( MXCMU, MXCMU ),  &
          EXPBEA( 0:MXCLY ), FLDIR( MXULV ), FLDN( MXULV ),  &
          FLYR( MXCLY ), GC( MXCMU, MXCMU, MXCLY ),  &
          GL( 0:MXCMU, MXCLY ), GU( MXUMU, MXCMU, MXCLY ),  &
          KK( MXCMU, MXCLY ), LL( MXCMU, MXCLY ), OPRIM( MXCLY ),  &
          PHASA( MXCLY ), PHAST( MXCLY ), PHASM( MXCLY ),  &
          PHIRAD( MXPHI ), PKAG( 0:MXCLY ), PSI0( MXCMU ),  &
          PSI1( MXCMU ), RMU( MXUMU, 0:MI ), SQT( MXSQT ),  &
          TAUC( 0:MXCLY ), TAUCPR( 0:MXCLY ), U0C( MXCMU, MXULV ),  &
          U0U( MXUMU, MXULV ), UTAUPR( MXULV ),  &
          UUM( MXUMU, MXULV ), WK( MXCMU ), XR0( MXCLY ),  &
          XR1( MXCLY ), YLM0( 0:MXCMU ), YLMC( 0:MXCMU, MXCMU ),  &
          YLMU( 0:MXCMU, MXUMU ), Z( NNLYRI ), Z0( MXCMU ),  &
          Z0U( MXUMU, MXCLY ), Z1( MXCMU ), Z1U( MXUMU, MXCLY ),  &
          ZBEAM( MXUMU, MXCLY )
      REAL :: ZJ( MXCMU ), ZPLK0( MXCMU, MXCLY ),  &
          ZPLK1( MXCMU, MXCLY ), ZZ( MXCMU, MXCLY )
      
      REAL :: AAD( MI, MI ), EVALD( MI ), EVECCD( MI, MI ), WKD( MXCMU )
!     ..
!     .. External Functions ..
      
      REAL :: PLKAVG, D1MACH, RATIO
      EXTERNAL  PLKAVG, D1MACH, RATIO
!     ..
!     .. External Subroutines ..
      
      EXTERNAL  ALBTRN, CHEKIN, CMPINT, FLUXES, INTCOR, LEPOLY, PRAVIN,  &
          PRTINP, PRTINT, SETDIS, SETMTX, SLFTST, SOLEIG, SOLVE0,  &
          SURFAC, TERPEV, TERPSO, UPBEAM, UPISOT, USRINT, ZEROAL, ZEROIT
!     ..
!     .. Intrinsic Functions ..
      
      INTRINSIC ABS, ASIN, COS, FLOAT, LEN, MAX, SQRT
!     ..
      SAVE      DITHER, PASS1, PI, RPD, SQT
      DATA      PASS1 / .TRUE. /, PRNTU0 / 2*.FALSE. /
      
      DELTAM = .FALSE.       !!!!!*** to agree better with RTSPEC *** !!!!!
      DELTAM = .TRUE.        !!!!!**** originally set by Stamnes **** !!!!!!
      CORINT = .TRUE.
      
      IF( PASS1 ) THEN
        
        PI     = 2.*ASIN( 1.0 )
        DITHER = 10.*D1MACH( 4 )
        
!                            ** Must dither more on high (14-digit)
!                            ** precision machine
        
        IF( DITHER < 1.E-10 ) DITHER = 10.*DITHER
        
        RPD  = PI / 180.0
        
        DO  NS = 1, MXSQT
          SQT( NS ) = SQRT( FLOAT( NS ) )
        END DO
!                            ** Set input values for self-test.
!                            ** Be sure SLFTST sets all print flags off.
        COMPAR = .FALSE.
        
        CALL SLFTST( CORINT, ACCUR, ALBEDO  BTEMP, DELTAM, DTAUC( 1 ),  &
              FBEAM, FISOT, IBCND, LAMBER, NLYR, PLANK, NPHI,  &
              NUMU, NSTR, NTAU, ONLYFL, PHI( 1 ), PHI0, NMOM,  &
              PMOM( 0,1 ), PRNT, PRNTU0, SSALB( 1 ), TEMIS,  &
              TEMPER( 0 ), TTEMP, UMU( 1 ), USRANG, USRTAU,  &
              UTAU( 1 ), UMU0, WVNMHI, WVNMLO, COMPAR, DUM, DUM, DUM, DUM )
        END IF
        
        
        20 CONTINUE
        
!sergio comment out this header as it puts too much stuff into warning.msg
!      IF( .NOT.PASS1 .AND. LEN( HEADER ).NE.0 )
!     &    WRITE( KSTDWARN,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' )
!     &    ' DISORT: '//HEADER
        
!                                  ** Calculate cumulative optical depth
!                                  ** and dither single-scatter albedo
!                                  ** to improve numerical behavior of
!                                  ** eigenvalue/vector computation
        CALL ZEROIT( TAUC, MXCLY + 1 )
        
        DO  LC = 1, NLYR
          IF( SSALB( LC ) == 1.0 ) SSALB( LC ) = 1.0 - DITHER
          TAUC( LC ) = TAUC( LC - 1 ) + DTAUC( LC )
        END DO
!                                ** Check input dimensions and variables
        
        CALL CHEKIN( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,  &
            WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,  &
            NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,  &
        PHI0, FISOT, LAMBER, ALBEDO  BTEMP, TTEMP,  &
              TEMIS, PLANK, ONLYFL, DELTAM, CORINT, ACCUR,  &
              TAUC, MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM,  &
              MXCLY, MXULV, MXUMU, MXCMU, MXPHI, MXSQT )
          
!                                 ** Zero internal and output arrays
          
          CALL  ZEROAL( MXCLY, EXPBEA(1), FLYR, OPRIM, PHASA, PHAST, PHASM,  &
              TAUCPR(1), XR0, XR1, MXCMU, CMU, CWT, PSI0, PSI1, WK, Z0, Z1, ZJ,  &
              MXCMU+1, YLM0, MXCMU**2, ARRAY, CC, EVECC,  &
              (MXCMU+1)*MXCLY, GL, (MXCMU+1)*MXCMU, YLMC,  &
              (MXCMU+1)*MXUMU, YLMU, MXCMU*MXCLY, KK, LL, ZZ, ZPLK0, ZPLK1,  &
              MXCMU**2*MXCLY, GC, MXULV, LAYRU, UTAUPR,  &
              MXUMU*MXCMU*MXCLY, GU, MXUMU*MXCLY, Z0U, Z1U, ZBEAM,  &
              MI, EVAL, MI**2, AMB, APB,  &
              NNLYRI, IPVT, Z, MAXULV, RFLDIR, RFLDN, FLUP, UAVG, DFDT,  &
              MAXUMU, ALBMED, TRNMED, MXUMU*MXULV, U0U,  &
              MAXUMU*MAXULV*MAXPHI, UU )
          
!                                 ** Perform various setup operations
          
          CALL SETDIS( CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM, FLYR,  &
              GL, IBCND, LAYRU, LYRCUT, MAXMOM, MAXUMU, MXCMU,  &
              NCUT, NLYR, NTAU, NN, NSTR, PLANK, NUMU, ONLYFL,  &
              CORINT, OPRIM, PMOM, SSALB, TAUC, TAUCPR, UTAU,  &
              UTAUPR, UMU, UMU0, USRTAU, USRANG )
          
!                                 ** Print input information
          IF( PRNT( 1 ) )
!sergio got rid of maxmom
!     &    CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER,
!     &                 WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,
!     &                 NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT,
!     &                 LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM,
!     &                 PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT,
!     &                 OPRIM, TAUC, TAUCPR, MAXMOM, PRNT( 5 ) )  &
          CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER,  &
              WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,  &
              NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT,  &
          LAMBER, ALBEDO  BTEMP, TTEMP, TEMIS, DELTAM,  &
                PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT,  &
                OPRIM, TAUC, TAUCPR, PRNT( 5 ) )
            
!                              ** Handle special case for getting albedo
!                              ** and transmissivity of medium for many
!                              ** beam angles at once
            IF( IBCND == 1 ) THEN
              
              CALL ALBTRN( ALBEDO  AMB, APB, ARRAY, B, BDR, CBAND, CC, CMU,  &
                    CWT, DTAUCP, EVAL, EVECC, GL, GC, GU, IPVT, KK,  &
                    LL, NLYR, NN, NSTR, NUMU, PRNT, TAUCPR, UMU, U0U,  &
                    WK, YLMC, YLMU, Z, WKD, MI, MI9M2, MAXUMU, MXCMU,  &
                    MXUMU, NNLYRI, SQT, ALBMED, TRNMED )
                RETURN
                
              END IF
!                                   ** Calculate Planck functions
              IF( .NOT.PLANK ) THEN
                
                BPLANK = 0.0
                TPLANK = 0.0
                CALL ZEROIT( PKAG,  MXCLY + 1 )
                
              ELSE
                
                TPLANK = TEMIS*PLKAVG( WVNMLO, WVNMHI, TTEMP )
                BPLANK =       PLKAVG( WVNMLO, WVNMHI, BTEMP )
                
                DO  LEV = 0, NLYR
                  PKAG( LEV ) = PLKAVG( WVNMLO, WVNMHI, TEMPER( LEV ) )
                END DO
                
              END IF
              
              
! ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
!           (EQ STWJ 5, STWL 6)
              
              KCONV  = 0
              NAZ    = NSTR - 1
!                                    ** Azimuth-independent case
              
              IF( FBEAM == 0.0 .OR. ABS(1.-UMU0) < 1.E-5 .OR. ONLYFL .OR.  &
                  ( NUMU == 1 .AND. ABS(1.-UMU(1)) < 1.E-5 ) .OR.  &
                  ( NUMU == 1 .AND. ABS(1.+UMU(1)) < 1.E-5 ) .OR.  &
                  ( NUMU == 2 .AND. ABS(1.+UMU(1)) < 1.E-5 .AND.  &
                  ABS(1.-UMU(2)) < 1.E-5 ) ) NAZ = 0
              
              
              DO  MAZIM = 0, NAZ
                
                IF( MAZIM == 0 ) DELM0  = 1.0
                IF( MAZIM > 0 ) DELM0  = 0.0
                
!                             ** Get normalized associated Legendre
!                             ** polynomials for
!                             ** (a) incident beam angle cosine
!                             ** (b) computational and user polar angle
!                             **     cosines
                IF( FBEAM > 0.0 ) THEN
                  
                  NCOS   = 1
                  ANGCOS = -UMU0
                  
                  CALL LEPOLY( NCOS, MAZIM, MXCMU, NSTR-1, ANGCOS, SQT, YLM0 )
                  
                END IF
                
                
                IF( .NOT.ONLYFL .AND. USRANG )  &
                    CALL LEPOLY( NUMU, MAZIM, MXCMU, NSTR-1, UMU, SQT, YLMU )
                
                CALL LEPOLY( NN, MAZIM, MXCMU, NSTR-1, CMU, SQT, YLMC )
                
!                       ** Get normalized associated Legendre polys.
!                       ** with negative arguments from those with
!                       ** positive arguments; Dave/Armstrong Eq. (15),
!                       ** STWL(59)
                SGN  = -1.0
                
                DO  L = MAZIM, NSTR - 1
                  
                  SGN  = -SGN
                  
                  DO  IQ = NN + 1, NSTR
                    YLMC( L, IQ ) = SGN*YLMC( L, IQ - NN )
                  END DO
                  
                END DO
!                                 ** Specify users bottom reflectivity
!                                 ** and emissivity properties
                IF( .NOT.LYRCUT )  &
                CALL SURFAC( ALBEDO  DELM0, CMU, FBEAM, LAMBER, MI, MAZIM,  &
                      MXUMU, NN, NUMU, ONLYFL, PI, UMU, UMU0,  &
                      USRANG, WVNMLO, WVNMHI, BDR, EMU, BEM, RMU )
                  
                  
! ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============
                  
                  DO  LC = 1, NCUT
                    
!                      ** Solve eigenfunction problem in Eq. STWJ(8B),
!                      ** STWL(23f); return eigenvalues and eigenvectors
                    
                    CALL SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL( 0,LC ), MI,  &
                        MAZIM, MXCMU, NN, NSTR, YLMC, CC, EVECC, EVAL,  &
                        KK( 1,LC ), GC( 1,1,LC ), WKD )
                    
!                                  ** Calculate particular solutions of
!                                  ** Eq. SS(18), STWL(24a) for incident
!                                  ** beam source
                    IF( FBEAM > 0.0 )  &
                        CALL UPBEAM( ARRAY, CC, CMU, DELM0, FBEAM, GL( 0,LC ),  &
                        IPVT, MAZIM, MXCMU, NN, NSTR, PI, UMU0, WK,  &
                        YLM0, YLMC, ZJ, ZZ( 1,LC ) )
                    
!                              ** Calculate particular solutions of Eq.
!                              ** SS(15), STWL(25) for thermal emission
!                              ** source
                    
                    IF( PLANK .AND. MAZIM == 0 ) THEN
                      
                      XR1( LC ) = 0.0
                      
                      IF( DTAUCP( LC ) > 0.0 ) XR1( LC ) =  &
                          ( PKAG( LC ) - PKAG( LC-1 ) ) / DTAUCP( LC )
                      
                      XR0( LC ) = PKAG( LC-1 ) - XR1( LC )*TAUCPR( LC-1 )
                      
                      CALL UPISOT( ARRAY, CC, CMU, IPVT, MXCMU, NN, NSTR,  &
                          OPRIM( LC ), WK, XR0( LC ), XR1( LC ),  &
                          Z0, Z1, ZPLK0( 1,LC ), ZPLK1( 1,LC ) )
                    END IF
                    
                    
                    IF( .NOT.ONLYFL .AND. USRANG ) THEN
                      
!                                            ** Interpolate eigenvectors
!                                            ** to user angles
                      
                      CALL TERPEV( CWT, EVECC, GL( 0,LC ), GU( 1,1,LC ), MAZIM,  &
                          MXCMU, MXUMU, NN, NSTR, NUMU, WK, YLMC, YLMU )
!                                            ** Interpolate source terms
!                                            ** to user angles
                      
                      CALL TERPSO( CWT, DELM0, FBEAM, GL( 0,LC ), MAZIM, MXCMU,  &
                          PLANK, NUMU, NSTR, OPRIM( LC ), PI, YLM0,  &
                          YLMC, YLMU, PSI0, PSI1, XR0( LC ),  &
                          XR1( LC ), Z0, Z1, ZJ, ZBEAM( 1,LC ),  &
                          Z0U( 1,LC ), Z1U( 1,LC ) )
                      
                    END IF
                    
                  END DO
                  
! ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============
                  
                  
!                      ** Set coefficient matrix of equations combining
!                      ** boundary and layer interface conditions
                  
                  CALL SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,  &
                      LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT,  &
                      NNLYRI, NN, NSTR, TAUCPR, WK )
                  
!                      ** Solve for constants of integration in homo-
!                      ** geneous solution (general boundary conditions)
                  
                  CALL SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,  &
                      FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM, MI,  &
                      MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, NNLYRI, PI,  &
                      TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )
                  
!                                  ** Compute upward and downward fluxes
                  
                  IF( MAZIM == 0 )
!sergio comment out maxulv
!     &       CALL FLUXES( CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,
!     &                    MAXULV, MXCMU, MXULV, NCUT, NN, NSTR, NTAU,
!     &                    PI, PRNT, PRNTU0( 1 ), SSALB, TAUCPR, UMU0,
!     &                    UTAU, UTAUPR, XR0, XR1, ZZ, ZPLK0, ZPLK1,
!     &                    DFDT, FLUP, FLDN, FLDIR, RFLDIR, RFLDN, UAVG,
!     &                    U0C )
                   CALL FLUXES( CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,  &
                      MXCMU, MXULV, NCUT, NN, NSTR, NTAU,  &
                      PI, PRNT, PRNTU0( 1 ), SSALB, TAUCPR, UMU0,  &
                      UTAU, UTAUPR, XR0, XR1, ZZ, ZPLK0, ZPLK1,  &
                      DFDT, FLUP, FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C )
                  
                  IF( ONLYFL ) THEN
                    
                    IF( MAXUMU >= NSTR ) THEN
!                                     ** Save azimuthal-avg intensities
!                                     ** at quadrature angles
                      DO  LU = 1, NTAU
                        
                        DO  IQ = 1, NSTR
                          U0U( IQ, LU ) = U0C( IQ, LU )
                        END DO
                        
                      END DO
                      
                    END IF
                    
                    GO TO 190
                    
                  END IF
                  
                  
                  CALL ZEROIT( UUM, MXUMU*MXULV )
                  
                  IF( USRANG ) THEN
!                                     ** Compute azimuthal intensity
!                                     ** components at user angles
                    
                    CALL USRINT( BPLANK, CMU, CWT, DELM0, DTAUCP, EMU, EXPBEA,  &
                        FBEAM, FISOT, GC, GU, KK, LAMBER, LAYRU, LL,  &
                        LYRCUT, MAZIM, MXCMU, MXULV, MXUMU, NCUT, NLYR,  &
                        NN, NSTR, PLANK, NUMU, NTAU, PI, RMU, TAUCPR,  &
                        TPLANK, UMU, UMU0, UTAUPR, WK, ZBEAM, Z0U, Z1U,  &
                        ZZ, ZPLK0, ZPLK1, UUM )
                    
                  ELSE
!                                     ** Compute azimuthal intensity
!                                     ** components at quadrature angles
                    
                    CALL CMPINT( FBEAM, GC, KK, LAYRU, LL, LYRCUT, MAZIM, MXCMU,  &
                        MXULV, MXUMU, NCUT, NN, NSTR, PLANK, NTAU,  &
                        TAUCPR, UMU0, UTAUPR, ZZ, ZPLK0, ZPLK1, UUM )
                    
                  END IF
                  
                  
                  IF( MAZIM == 0 ) THEN
!                               ** Save azimuthally averaged intensities
                    
                    DO  LU = 1, NTAU
                      
                      DO  IU = 1, NUMU
                        U0U( IU, LU ) = UUM( IU, LU )
                        
                        DO  J = 1, NPHI
                          UU( IU, LU, J ) = UUM( IU, LU )
                        END DO
                        
                      END DO
                      
                    END DO
!                              ** Print azimuthally averaged intensities
!                              ** at user angles
                    
                    IF( PRNTU0( 2 ) )  &
                        CALL PRAVIN( UMU, NUMU, MXUMU, UTAU, NTAU, U0U )
                    
                    IF( NAZ > 0 ) THEN
                      
                      CALL ZEROIT( PHIRAD, MXPHI )
                      DO  J = 1, NPHI
                        PHIRAD( J ) = RPD*( PHI( J ) - PHI0 )
                      END DO
                      
                    END IF
                    
                    
                  ELSE
!                                ** Increment intensity by current
!                                ** azimuthal component (Fourier
!                                ** cosine series);  Eq SD(2), STWL(6)
                    AZERR  = 0.0
                    
                    DO  J = 1, NPHI
                      
                      COSPHI = COS( MAZIM*PHIRAD( J ) )
                      
                      DO  LU = 1, NTAU
                        
                        DO  IU = 1, NUMU
                          AZTERM = UUM( IU, LU )*COSPHI
                          UU( IU, LU, J ) = UU( IU, LU, J ) + AZTERM
                          AZERR  = MAX( AZERR,  &
                              RATIO( ABS(AZTERM), ABS(UU(IU,LU,J)) ) )
                        END DO
                        
                      END DO
                      
                    END DO
                    
                    IF( AZERR <= ACCUR ) KCONV  = KCONV + 1
                    
                    IF( KCONV >= 2 ) GO TO 190
                    
                  END IF
                  
                END DO
                
! ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============
                
                
                190 CONTINUE
                
!                                    ** Apply Nakajima/Tanaka intensity
!                                    ** corrections
                
                IF( CORINT )  &
                    CALL INTCOR( DITHER, FBEAM, FLYR, LAYRU, LYRCUT, MAXMOM,  &
                    MAXULV, MAXUMU, NMOM, NCUT, NPHI, NSTR, NTAU,  &
                    NUMU, OPRIM, PHASA, PHAST, PHASM, PHIRAD, PI,  &
                    RPD, PMOM, SSALB, DTAUC, TAUC, TAUCPR, UMU,  &
                    UMU0, UTAU, UTAUPR, UU )
                
!                                          ** Print intensities
                
                IF( PRNT( 3 ) .AND. .NOT.ONLYFL )
!sergio comment out maxulv,maxumu
!     &    CALL PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
!     &                 MAXUMU )  &
                CALL PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI)
                
                
                IF( PASS1 ) THEN
!                                    ** Compare test case results with
!                                    ** correct answers and abort if bad
                  COMPAR = .TRUE.
                  
                  CALL SLFTST( CORINT, ACCUR, ALBEDO  BTEMP, DELTAM, DTAUC( 1 ),  &
                        FBEAM, FISOT, IBCND, LAMBER, NLYR, PLANK, NPHI,  &
                        NUMU, NSTR, NTAU, ONLYFL, PHI( 1 ), PHI0, NMOM,  &
                        PMOM( 0,1 ), PRNT, PRNTU0, SSALB( 1 ), TEMIS,  &
                        TEMPER( 0 ), TTEMP, UMU( 1 ), USRANG, USRTAU,  &
                        UTAU( 1 ), UMU0, WVNMHI, WVNMLO, COMPAR,  &
                        FLUP( 1 ), RFLDIR( 1 ), RFLDN( 1 ), UU( 1,1,1 ) )
                    
                    PASS1  = .FALSE.
                    GO TO 20
                    
                  END IF
                  
                  
                  RETURN
                END SUBROUTINE DISORT
                
!************************************************************************
                
                SUBROUTINE ASYMTX( AA, EVEC, EVAL, M, IA, IEVEC, IER, WK )
                
!    ======  S I N G L E    P R E C I S I O N    V E R S I O N  ======
                
!       Solves eigenfunction problem for real asymmetric matrix
!       for which it is known a priori that the eigenvalues are real.
                
!       This is an adaptation of a subroutine EIGRF in the IMSL
!       library to use real instead of complex arithmetic, accounting
!       for the known fact that the eigenvalues and eigenvectors in
!       the discrete ordinate solution are real.  Other changes include
!       putting all the called subroutines in-line, deleting the
!       performance index calculation, updating many DO-loops
!       to Fortran77, and in calculating the machine precision
!       TOL instead of specifying it in a data statement.
                
!       EIGRF is based primarily on EISPACK routines.  The matrix is
!       first balanced using the Parlett-Reinsch algorithm.  Then
!       the Martin-Wilkinson algorithm is applied.
                
!       There is a statement 'J  = WK( I )' that converts a REAL
!       variable to an integer variable, that seems dangerous
!       to us in principle, but seems to work fine in practice.
                
!       References:
!          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
!             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
!             Sources and Development of Mathematical Software,
!             Prentice-Hall, Englewood Cliffs, NJ
!         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
!             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
!         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
!             Clarendon Press, Oxford
                
                
!   I N P U T    V A R I A B L E S:
                
!       AA    :  input asymmetric matrix, destroyed after solved
                
!        M    :  order of  AA
                
!       IA    :  first dimension of  AA
                
!    IEVEC    :  first dimension of  EVEC
                
                
!   O U T P U T    V A R I A B L E S:
                
!       EVEC  :  (unnormalized) eigenvectors of  AA
!                   ( column J corresponds to EVAL(J) )
                
!       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M )
                
!       IER   :  if .NE. 0, signals that EVAL(IER) failed to converge;
!                   in that case eigenvalues IER+1,IER+2,...,M  are
!                   correct but eigenvalues 1,...,IER are set to zero.
                
                
!   S C R A T C H   V A R I A B L E S:
                
!       WK    :  work area ( dimension at least 2*M )
                
!   Called by- SOLEIG
!   Calls- D1MACH, ERRMSG
! +-------------------------------------------------------------------+
                
!     .. Scalar Arguments ..
                
                
                REAL, INTENT(IN OUT)                     :: AA( IA, M )
                REAL, INTENT(OUT)                        :: EVEC( IEVEC, M )
                REAL, INTENT(OUT)                        :: EVAL( M )
                INTEGER, INTENT(IN OUT)                  :: M
                INTEGER, INTENT(IN)                      :: IA
                INTEGER, INTENT(IN)                      :: IEVEC
                INTEGER, INTENT(OUT)                     :: IER
                REAL, INTENT(OUT)                        :: WK( * )
                
!     ..
!     .. Array Arguments ..
                
                
!     ..
!     .. Local Scalars ..
                
                LOGICAL :: NOCONV, NOTLAS
                INTEGER :: I, II, IN, J, K, KA, KKK, L, LB, LLL, N, N1, N2
                REAL :: C1, C2, C3, C4, C5, C6, COL, DISCRI, F, G, H, ONE, P, Q,  &
                    R, REPL, RNORM, ROW, S, SCALE, SGN, T, TOL, UU, VV, W,  &
                    X, Y, Z, ZERO
!     ..
!     .. External Functions ..
                
                REAL :: D1MACH
                EXTERNAL  D1MACH
!     ..
!     .. External Subroutines ..
                
                EXTERNAL  ERRMSG
!     ..
!     .. Intrinsic Functions ..
                
                INTRINSIC ABS, MIN, SIGN, SQRT
!     ..
                DATA      C1 / 0.4375 / , C2 / 0.5 / , C3 / 0.75 / , C4 / 0.95 / ,  &
                    C5 / 16.0 / , C6 / 256.0 / , ZERO / 0.0 / , ONE / 1.0 /
                
                
                IER  = 0
                TOL  = D1MACH( 4 )
                
                IF( M < 1 .OR. IA < M .OR. IEVEC < M )  &
                    CALL ERRMSG( 'ASYMTX--bad input variable(s)', .TRUE. )
                
                
!                           ** Handle 1x1 and 2x2 special cases
                IF( M == 1 ) THEN
                  
                  EVAL( 1 )   = AA( 1,1 )
                  EVEC( 1,1 ) = ONE
                  RETURN
                  
                ELSE IF( M == 2 ) THEN
                  
                  DISCRI = ( AA( 1,1 ) - AA( 2,2 ) )**2 + 4.*AA( 1,2 )*AA( 2,1 )
                  
                  IF( DISCRI < ZERO )  &
                      CALL ERRMSG( 'ASYMTX--complex evals in 2x2 case',.TRUE. )
                  
                  SGN  = ONE
                  
                  IF( AA( 1,1 ) < AA( 2,2 ) ) SGN  = - ONE
                  
                  EVAL( 1 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) + SGN*SQRT( DISCRI ) )
                  EVAL( 2 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) - SGN*SQRT( DISCRI ) )
                  EVEC( 1,1 ) = ONE
                  EVEC( 2,2 ) = ONE
                  
                  IF( AA( 1,1 ) == AA( 2,2 ) .AND.  &
                        ( AA( 2,1 ) == ZERO .OR. AA( 1,2 ) == ZERO ) ) THEN
                    
                    RNORM = ABS( AA( 1,1 ) ) + ABS( AA( 1,2 ) ) +  &
                        ABS( AA( 2,1 ) ) + ABS( AA( 2,2 ) )
                    W     = TOL * RNORM
                    EVEC( 2,1 ) =   AA( 2,1 ) / W
                    EVEC( 1,2 ) = - AA( 1,2 ) / W
                    
                  ELSE
                    
                    EVEC( 2,1 ) = AA( 2,1 ) / ( EVAL( 1 ) - AA( 2,2 ) )
                    EVEC( 1,2 ) = AA( 1,2 ) / ( EVAL( 2 ) - AA( 1,1 ) )
                    
                  END IF
                  
                  RETURN
                  
                END IF
                
!                                        ** Initialize output variables
                IER  = 0
                
                DO  I = 1, M
                  
                  EVAL( I ) = ZERO
                  
                  DO  J = 1, M
                    EVEC( I, J ) = ZERO
                  END DO
                  
                  EVEC( I, I ) = ONE
                  
                END DO
!                  ** Balance the input matrix and reduce its norm by
!                  ** diagonal similarity transformation stored in WK;
!                  ** then search for rows isolating an eigenvalue
!                  ** and push them down
                RNORM  = ZERO
                L  = 1
                K  = M
                
                30 CONTINUE
                KKK  = K
                
                DO  J = KKK, 1, -1
                  
                  ROW  = ZERO
                  
                  DO  I = 1, K
                    IF( I /= J ) ROW  = ROW + ABS( AA( J,I ) )
                  END DO
                  
                  IF( ROW == ZERO ) THEN
                    
                    WK( K ) = J
                    
                    IF( J /= K ) THEN
                      
                      DO  I = 1, K
                        REPL   = AA( I, J )
                        AA( I, J ) = AA( I, K )
                        AA( I, K ) = REPL
                      END DO
                      
                      DO  I = L, M
                        REPL   = AA( J, I )
                        AA( J, I ) = AA( K, I )
                        AA( K, I ) = REPL
                      END DO
                      
                    END IF
                    
                    K  = K - 1
                    GO TO  30
                    
                  END IF
                  
                END DO
!                                  ** Search for columns isolating an
!                                  ** eigenvalue and push them left
                80 CONTINUE
                LLL  = L
                
                DO  J = LLL, K
                  
                  COL  = ZERO
                  
                  DO  I = L, K
                    IF( I /= J ) COL  = COL + ABS( AA( I,J ) )
                  END DO
                  
                  IF( COL == ZERO ) THEN
                    
                    WK( L ) = J
                    
                    IF( J /= L ) THEN
                      
                      DO  I = 1, K
                        REPL       = AA( I, J )
                        AA( I, J ) = AA( I, L )
                        AA( I, L ) = REPL
                      END DO
                      
                      DO  I = L, M
                        REPL       = AA( J, I )
                        AA( J, I ) = AA( L, I )
                        AA( L, I ) = REPL
                      END DO
                      
                    END IF
                    
                    L  = L + 1
                    GO TO  80
                    
                  END IF
                  
                END DO
                
!                           ** Balance the submatrix in rows L through K
                DO  I = L, K
                  WK( I ) = ONE
                END DO
                
                140 CONTINUE
                NOCONV = .FALSE.
                
                DO  I = L, K
                  
                  COL  = ZERO
                  ROW  = ZERO
                  
                  DO  J = L, K
                    
                    IF( J /= I ) THEN
                      COL  = COL + ABS( AA( J,I ) )
                      ROW  = ROW + ABS( AA( I,J ) )
                    END IF
                    
                  END DO
                  
                  F  = ONE
                  G  = ROW / C5
                  H  = COL + ROW
                  
                  160    CONTINUE
                  IF( COL < G ) THEN
                    
                    F    = F*C5
                    COL  = COL*C6
                    GO TO  160
                    
                  END IF
                  
                  G  = ROW*C5
                  
                  170    CONTINUE
                  IF( COL >= G ) THEN
                    
                    F    = F / C5
                    COL  = COL / C6
                    GO TO  170
                    
                  END IF
!                                                ** Now balance
                  IF( ( COL + ROW ) / F < C4*H ) THEN
                    
                    WK( I ) = WK( I )*F
                    NOCONV = .TRUE.
                    
                    DO  J = L, M
                      AA( I, J ) = AA( I, J ) / F
                    END DO
                    
                    DO  J = 1, K
                      AA( J, I ) = AA( J, I )*F
                    END DO
                    
                  END IF
                  
                END DO
                
                
                IF( NOCONV ) GO TO  140
!                                   ** Is A already in Hessenberg form?
                
                IF( K-1 < L+1 ) GO TO  350
                
!                                   ** Transfer A to a Hessenberg form
                DO  N = L + 1, K - 1
                  
                  H  = ZERO
                  WK( N + M ) = ZERO
                  SCALE  = ZERO
!                                                 ** Scale column
                  DO  I = N, K
                    SCALE  = SCALE + ABS( AA( I,N - 1 ) )
                  END DO
                  
                  IF( SCALE /= ZERO ) THEN
                    
                    DO  I = K, N, -1
                      WK( I + M ) = AA( I, N - 1 ) / SCALE
                      H  = H + WK( I + M )**2
                    END DO
                    
                    G  = -SIGN( SQRT( H ), WK( N + M ) )
                    H  = H - WK( N + M )*G
                    WK( N + M ) = WK( N + M ) - G
                    
!                                            ** Form (I-(U*UT)/H)*A
                    DO  J = N, M
                      
                      F  = ZERO
                      DO  I = K, N, -1
                        F  = F + WK( I + M )*AA( I, J )
                      END DO
                      
                      DO  I = N, K
                        AA( I, J ) = AA( I, J ) - WK( I + M )*F / H
                      END DO
                      
                    END DO
                    
!                                    ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
                    DO  I = 1, K
                      
                      F  = ZERO
                      DO  J = K, N, -1
                        F  = F + WK( J + M )*AA( I, J )
                      END DO
                      
                      DO  J = N, K
                        AA( I, J ) = AA( I, J ) - WK( J + M )*F / H
                      END DO
                      
                    END DO
                    
                    WK( N + M ) = SCALE*WK( N + M )
                    AA( N, N - 1 ) = SCALE*G
                    
                  END IF
                  
                END DO
                
                
                DO  N = K - 2, L, -1
                  
                  N1   = N + 1
                  N2   = N + 2
                  F  = AA( N + 1, N )
                  
                  IF( F /= ZERO ) THEN
                    
                    F  = F*WK( N + 1 + M )
                    
                    DO  I = N + 2, K
                      WK( I + M ) = AA( I, N )
                    END DO
                    
                    IF( N + 1 <= K ) THEN
                      
                      DO  J = 1, M
                        
                        G  = ZERO
                        DO  I = N + 1, K
                          G  = G + WK( I + M )*EVEC( I, J )
                        END DO
                        
                        G  = G / F
                        
                        DO  I = N + 1, K
                          EVEC( I, J ) = EVEC( I, J ) + G*WK( I + M )
                        END DO
                        
                      END DO
                      
                    END IF
                    
                  END IF
                  
                END DO
                
                
                350 CONTINUE
                N  = 1
                
                DO  I = 1, M
                  
                  DO  J = N, M
                    RNORM  = RNORM + ABS( AA( I,J ) )
                  END DO
                  
                  N  = I
                  IF( I < L .OR. I > K ) EVAL( I ) = AA( I, I )
                  
                END DO
                
                N  = K
                T  = ZERO
                
!                                      ** Search for next eigenvalues
                380 CONTINUE
                IF( N < L ) GO TO  530
                
                IN   = 0
                N1   = N - 1
                N2   = N - 2
!                          ** Look for single small sub-diagonal element
                390 CONTINUE
                
                DO  I = L, N
                  LB   = N + L - I
                  
                  IF( LB == L ) GO TO  410
                  
                  S  = ABS( AA( LB - 1,LB - 1 ) ) + ABS( AA( LB,LB ) )
                  IF( S == ZERO ) S  = RNORM
                  
                  IF( ABS( AA( LB,LB - 1 ) ) <= TOL*S ) GO TO  410
                  
                END DO
                
                
                410 CONTINUE
                X  = AA( N, N )
                
                IF( LB == N ) THEN
!                                        ** One eigenvalue found
                  AA( N, N ) = X + T
                  EVAL( N ) = AA( N, N )
                  N  = N1
                  GO TO  380
                  
                END IF
                
                
                Y  = AA( N1, N1 )
                W  = AA( N, N1 )*AA( N1, N )
                
                IF( LB == N1 ) THEN
!                                        ** Two eigenvalues found
                  P  = ( Y - X )*C2
                  Q  = P**2 + W
                  Z  = SQRT( ABS( Q ) )
                  AA( N, N ) = X + T
                  X  = AA( N, N )
                  AA( N1, N1 ) = Y + T
!                                        ** Real pair
                  Z  = P + SIGN( Z, P )
                  EVAL( N1 ) = X + Z
                  EVAL( N ) = EVAL( N1 )
                  
                  IF( Z /= ZERO ) EVAL( N ) = X - W / Z
                  
                  X  = AA( N, N1 )
!                                  ** Employ scale factor in case
!                                  ** X and Z are very small
                  R  = SQRT( X*X + Z*Z )
                  P  = X / R
                  Q  = Z / R
!                                             ** Row modification
                  
                  DO  J = N1, M
                    Z  = AA( N1, J )
                    AA( N1, J ) = Q*Z + P*AA( N, J )
                    AA( N, J ) = Q*AA( N, J ) - P*Z
                  END DO
!                                             ** Column modification
                  
                  DO  I = 1, N
                    Z  = AA( I, N1 )
                    AA( I, N1 ) = Q*Z + P*AA( I, N )
                    AA( I, N ) = Q*AA( I, N ) - P*Z
                  END DO
!                                          ** Accumulate transformations
                  
                  DO  I = L, K
                    Z  = EVEC( I, N1 )
                    EVEC( I, N1 ) = Q*Z + P*EVEC( I, N )
                    EVEC( I, N ) = Q*EVEC( I, N ) - P*Z
                  END DO
                  
                  N  = N2
                  GO TO  380
                  
                END IF
                
                
                IF( IN == 30 ) THEN
                  
!                    ** No convergence after 30 iterations; set error
!                    ** indicator to the index of the current eigenvalue
                  IER  = N
                  RETURN
                  
                END IF
!                                                  ** Form shift
                
                IF( IN == 10 .OR. IN == 20 ) THEN
                  
                  T  = T + X
                  
                  DO  I = L, N
                    AA( I, I ) = AA( I, I ) - X
                  END DO
                  
                  S  = ABS( AA( N,N1 ) ) + ABS( AA( N1,N2 ) )
                  X  = C3*S
                  Y  = X
                  W  = -C1*S**2
                  
                END IF
                
                
                IN   = IN + 1
                
!                ** Look for two consecutive small sub-diagonal elements
                
                DO  J = LB, N2
                  I  = N2 + LB - J
                  Z  = AA( I, I )
                  R  = X - Z
                  S  = Y - Z
                  P  = ( R*S - W ) / AA( I + 1, I ) + AA( I, I + 1 )
                  Q  = AA( I + 1, I + 1 ) - Z - R - S
                  R  = AA( I + 2, I + 1 )
                  S  = ABS( P ) + ABS( Q ) + ABS( R )
                  P  = P / S
                  Q  = Q / S
                  R  = R / S
                  
                  IF( I == LB ) GO TO  470
                  
                  UU   = ABS( AA( I,I - 1 ) )*( ABS( Q ) + ABS( R ) )
                  VV   = ABS( P )*( ABS( AA( I - 1,I - 1 ) ) + ABS( Z ) +  &
                      ABS( AA( I + 1,I + 1 ) ) )
                  
                  IF( UU <= TOL*VV ) GO TO  470
                  
                END DO
                
                470 CONTINUE
                
                AA( I + 2, I ) = ZERO
                
                DO  J = I + 3, N
                  AA( J, J - 2 ) = ZERO
                  AA( J, J - 3 ) = ZERO
                END DO
                
!             ** Double QR step involving rows K to N and columns M to N
                
                DO  KA = I, N1
                  
                  NOTLAS = KA /= N1
                  
                  IF( KA == I ) THEN
                    
                    S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
                    
                    IF( LB /= I ) AA( KA, KA - 1 ) = - AA( KA, KA - 1 )
                    
                  ELSE
                    
                    P  = AA( KA, KA - 1 )
                    Q  = AA( KA + 1, KA - 1 )
                    R  = ZERO
                    
                    IF( NOTLAS ) R  = AA( KA + 2, KA - 1 )
                    
                    X  = ABS( P ) + ABS( Q ) + ABS( R )
                    
                    IF( X == ZERO ) CYCLE
                    
                    P  = P / X
                    Q  = Q / X
                    R  = R / X
                    S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
                    AA( KA, KA - 1 ) = -S*X
                    
                  END IF
                  
                  P  = P + S
                  X  = P / S
                  Y  = Q / S
                  Z  = R / S
                  Q  = Q / P
                  R  = R / P
!                                              ** Row modification
                  
                  DO  J = KA, M
                    P  = AA( KA, J ) + Q*AA( KA + 1, J )
                    
                    IF( NOTLAS ) THEN
                      
                      P  = P + R*AA( KA + 2, J )
                      AA( KA + 2, J ) = AA( KA + 2, J ) - P*Z
                      
                    END IF
                    
                    AA( KA + 1, J ) = AA( KA + 1, J ) - P*Y
                    AA( KA, J ) = AA( KA, J ) - P*X
                    
                  END DO
!                                                 ** Column modification
                  
                  DO  II = 1, MIN( N, KA + 3 )
                    P  = X*AA( II, KA ) + Y*AA( II, KA + 1 )
                    
                    IF( NOTLAS ) THEN
                      
                      P  = P + Z*AA( II, KA + 2 )
                      AA( II, KA + 2 ) = AA( II, KA + 2 ) - P*R
                      
                    END IF
                    
                    AA( II, KA + 1 ) = AA( II, KA + 1 ) - P*Q
                    AA( II, KA ) = AA( II, KA ) - P
                    
                  END DO
!                                          ** Accumulate transformations
                  
                  DO  II = L, K
                    
                    P  = X*EVEC( II, KA ) + Y*EVEC( II, KA + 1 )
                    
                    IF( NOTLAS ) THEN
                      
                      P  = P + Z*EVEC( II, KA + 2 )
                      EVEC( II, KA + 2 ) = EVEC( II, KA + 2 ) - P*R
                      
                    END IF
                    
                    EVEC( II, KA + 1 ) = EVEC( II, KA + 1 ) - P*Q
                    EVEC( II, KA ) = EVEC( II, KA ) - P
                    
                  END DO
                  
                END DO
                
                GO TO  390
                
!                     ** All evals found, now backsubstitute real vector
                530 CONTINUE
                
                IF( RNORM /= ZERO ) THEN
                  
                  DO  N = M, 1, -1
                    
                    N2   = N
                    AA( N, N ) = ONE
                    
                    DO  I = N - 1, 1, -1
                      
                      W  = AA( I, I ) - EVAL( N )
                      
                      IF( W == ZERO ) W  = TOL*RNORM
                      
                      R  = AA( I, N )
                      DO  J = N2, N - 1
                        R  = R + AA( I, J )*AA( J, N )
                      END DO
                      
                      AA( I, N ) = -R / W
                      N2   = I
                      
                    END DO
                    
                  END DO
!                      ** End backsubstitution vectors of isolated evals
                  
                  DO  I = 1, M
                    
                    IF( I < L .OR. I > K ) THEN
                      
                      DO  J = I, M
                        EVEC( I, J ) = AA( I, J )
                      END DO
                      
                    END IF
                    
                  END DO
!                                   ** Multiply by transformation matrix
                  
                  IF( K /= 0 ) THEN
                    
                    DO  J = M, L, -1
                      
                      DO  I = L, K
                        
                        Z  = ZERO
                        DO  N = L, MIN( J, K )
                          Z  = Z + EVEC( I, N )*AA( N, J )
                        END DO
                        
                        EVEC( I, J ) = Z
                        
                      END DO
                      
                    END DO
                    
                  END IF
                  
                END IF
                
                
                DO  I = L, K
                  
                  DO  J = 1, M
                    EVEC( I, J ) = EVEC( I, J )*WK( I )
                  END DO
                  
                END DO
!                           ** Interchange rows if permutations occurred
                
                DO  I = L - 1, 1, -1
                  
                  J  = WK( I )
                  
                  IF( I /= J ) THEN
                    
                    DO  N = 1, M
                      REPL   = EVEC( I, N )
                      EVEC( I, N ) = EVEC( J, N )
                      EVEC( J, N ) = REPL
                    END DO
                    
                  END IF
                  
                END DO
                
                
                DO  I = K + 1, M
                  
                  J  = WK( I )
                  
                  IF( I /= J ) THEN
                    
                    DO  N = 1, M
                      REPL   = EVEC( I, N )
                      EVEC( I, N ) = EVEC( J, N )
                      EVEC( J, N ) = REPL
                    END DO
                    
                  END IF
                  
                END DO
                
                
                RETURN
              END SUBROUTINE ASYMTX
              
!************************************************************************
              
              SUBROUTINE CMPINT( FBEAM, GC, KK, LAYRU, LL, LYRCUT, MAZIM, MXCMU,  &
                  MXULV, MXUMU, NCUT, NN, NSTR, PLANK, NTAU,  &
                  TAUCPR, UMU0, UTAUPR, ZZ, ZPLK0, ZPLK1, UUM )
              
!          Calculates the Fourier intensity components at the quadrature
!          angles for azimuthal expansion terms (MAZIM) in Eq. SD(2),
!          STWL(6)
              
              
!    I N P U T    V A R I A B L E S:
              
!       KK      :  Eigenvalues of coeff. matrix in Eq. SS(7), STWL(23b)
              
!       GC      :  Eigenvectors at polar quadrature angles in Eq. SC(1)
              
!       LL      :  Constants of integration in Eq. SC(1), obtained
!                  by solving scaled version of Eq. SC(5);
!                  exponential term of Eq. SC(12) not included
              
!       LYRCUT  :  Logical flag for truncation of computational layer
              
!       MAZIM   :  Order of azimuthal component
              
!       NCUT    :  Number of computational layer where absorption
!                  optical depth exceeds ABSCUT
              
!       NN      :  Order of double-Gauss quadrature (NSTR/2)
              
!       TAUCPR  :  Cumulative optical depth (delta-M-scaled)
              
!       UTAUPR  :  Optical depths of user output levels in delta-M
!                  coordinates;  equal to UTAU if no delta-M
              
!       ZZ      :  Beam source vectors in Eq. SS(19), STWL(24b)
              
!       ZPLK0   :  Thermal source vectors Z0, by solving Eq. SS(16),
!                  Y-sub-zero in STWL(26ab)
              
!       ZPLK1   :  Thermal source vectors Z1, by solving Eq. SS(16),
!                  Y-sub-one in STWL(26ab)
              
!       (Remainder are 'DISORT' input variables)
              
              
!    O U T P U T   V A R I A B L E S:
              
!       UUM     :  Fourier components of the intensity in Eq. SD(12)
!                    (at polar quadrature angles)
              
              
!    I N T E R N A L   V A R I A B L E S:
              
!       FACT    :  EXP( - UTAUPR / UMU0 )
!       ZINT    :  intensity of M=0 case, in Eq. SC(1)
              
!   Called by- DISORT
! +--------------------------------------------------------------------
              
!     .. Scalar Arguments ..
              
              
              REAL, INTENT(IN)                         :: FBEAM
              REAL, INTENT(IN)                         :: GC( MXCMU, MXCMU, * )
              REAL, INTENT(IN OUT)                     :: KK( MXCMU, * )
              INTEGER, INTENT(IN)                      :: LAYRU( * )
              REAL, INTENT(IN)                         :: LL( MXCMU, * )
              LOGICAL, INTENT(IN)                      :: LYRCUT
              INTEGER, INTENT(IN)                      :: MAZIM
              INTEGER, INTENT(IN OUT)                  :: MXCMU
              INTEGER, INTENT(IN OUT)                  :: MXULV
              INTEGER, INTENT(IN OUT)                  :: MXUMU
              INTEGER, INTENT(IN)                      :: NCUT
              INTEGER, INTENT(IN)                      :: NN
              INTEGER, INTENT(IN)                      :: NSTR
              LOGICAL, INTENT(IN)                      :: PLANK
              INTEGER, INTENT(IN)                      :: NTAU
              REAL, INTENT(IN OUT)                     :: TAUCPR( 0:* )
              REAL, INTENT(IN OUT)                     :: UMU0
              REAL, INTENT(IN OUT)                     :: UTAUPR( MXULV )
              REAL, INTENT(IN OUT)                     :: ZZ( MXCMU, * )
              REAL, INTENT(IN)                         :: ZPLK0( MXCMU, * )
              REAL, INTENT(IN OUT)                     :: ZPLK1( MXCMU, * )
              REAL, INTENT(OUT)                        :: UUM( MXUMU, MXULV )
              
              
              
!     ..
!     .. Array Arguments ..
              
              
              
!     ..
!     .. Local Scalars ..
              
              INTEGER :: IQ, JQ, LU, LYU
              REAL :: ZINT
!     ..
!     .. Intrinsic Functions ..
              
              INTRINSIC EXP
!     ..
              
!                                       ** Loop over user levels
              DO  LU = 1, NTAU
                
                LYU  = LAYRU( LU )
                
                IF( LYRCUT .AND. LYU > NCUT ) CYCLE
                
                DO  IQ = 1, NSTR
                  
                  ZINT = 0.0
                  
                  DO  JQ = 1, NN
                    ZINT = ZINT + GC( IQ, JQ, LYU ) * LL( JQ, LYU ) *  &
                        EXP( -KK( JQ,LYU )* ( UTAUPR( LU ) - TAUCPR( LYU ) ) )
                  END DO
                  
                  DO  JQ = NN + 1, NSTR
                    ZINT = ZINT + GC( IQ, JQ, LYU ) * LL( JQ, LYU ) *  &
                        EXP( -KK( JQ,LYU )* ( UTAUPR( LU ) - TAUCPR( LYU-1 ) ) )
                  END DO
                  
                  UUM( IQ, LU ) = ZINT
                  
                  IF( FBEAM > 0.0 ) UUM( IQ, LU ) = ZINT +  &
                      ZZ( IQ, LYU )*EXP( -UTAUPR( LU )/UMU0 )
                  
                  IF( PLANK .AND. MAZIM == 0 )  &
                      UUM( IQ, LU ) = UUM( IQ, LU ) + ZPLK0( IQ,LYU ) +  &
                      ZPLK1( IQ,LYU ) * UTAUPR( LU )
                END DO
                
              END DO
              
              
              RETURN
            END SUBROUTINE CMPINT
            
!************************************************************************
            
            SUBROUTINE FLUXES( CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,  &
                MXCMU, MXULV, NCUT, NN, NSTR, NTAU,  &
                PI, PRNT, PRNTU0, SSALB, TAUCPR, UMU0, UTAU,  &
                UTAUPR, XR0, XR1, ZZ, ZPLK0, ZPLK1, DFDT,  &
                FLUP, FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C )
            
            
            REAL, INTENT(IN)                         :: CMU( MXCMU )
            REAL, INTENT(IN)                         :: CWT( MXCMU )
            REAL, INTENT(IN)                         :: FBEAM
            REAL, INTENT(IN)                         :: GC( MXCMU, MXCMU, * )
            REAL, INTENT(IN OUT)                     :: KK( MXCMU, * )
            INTEGER, INTENT(IN)                      :: LAYRU( MXULV )
            REAL, INTENT(IN)                         :: LL( MXCMU, * )
            LOGICAL, INTENT(IN)                      :: LYRCUT
            INTEGER, INTENT(IN OUT)                  :: MXCMU
            INTEGER, INTENT(IN OUT)                  :: MXULV
            INTEGER, INTENT(IN)                      :: NCUT
            INTEGER, INTENT(IN)                      :: NN
            INTEGER, INTENT(IN)                      :: NSTR
            INTEGER, INTENT(IN)                      :: NTAU
            REAL, INTENT(IN)                         :: PI
            LOGICAL, INTENT(IN)                      :: PRNT( * )
            LOGICAL, INTENT(IN)                      :: PRNTU0
            REAL, INTENT(IN)                         :: SSALB( * )
            REAL, INTENT(IN OUT)                     :: TAUCPR( 0:* )
            REAL, INTENT(IN)                         :: UMU0
            REAL, INTENT(IN)                         :: UTAU( MAXULV )
            REAL, INTENT(IN)                         :: UTAUPR( MXULV )
            REAL, INTENT(IN)                         :: XR0( * )
            REAL, INTENT(IN)                         :: XR1( * )
            REAL, INTENT(IN)                         :: ZZ( MXCMU, * )
            REAL, INTENT(IN)                         :: ZPLK0( MXCMU, * )
            REAL, INTENT(IN OUT)                     :: ZPLK1( MXCMU, * )
            REAL, INTENT(OUT)                        :: DFDT( MAXULV )
            REAL, INTENT(OUT)                        :: FLUP( MAXULV )
            REAL, INTENT(OUT)                        :: FLDN( MXULV )
            REAL, INTENT(OUT)                        :: FLDIR( MXULV )
            REAL, INTENT(OUT)                        :: RFLDIR( MAXULV )
            REAL, INTENT(OUT)                        :: RFLDN( MAXULV )
            REAL, INTENT(OUT)                        :: UAVG( MAXULV )
            REAL, INTENT(OUT)                        :: U0C( MXCMU, MXULV )
            IMPLICIT NONE
            
            INCLUDE '../INCLUDE/scatterparam.f90'
            
!       Calculates the radiative fluxes, mean intensity, and flux
!       derivative with respect to optical depth from the m=0 intensity
!       components (the azimuthally-averaged intensity)
            
            
!    I N P U T     V A R I A B L E S:
            
!       CMU      :  Abscissae for Gauss quadrature over angle cosine
            
!       CWT      :  Weights for Gauss quadrature over angle cosine
            
!       GC       :  Eigenvectors at polar quadrature angles, SC(1)
            
!       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7), STWL(23b)
            
!       LAYRU    :  Layer number of user level UTAU
            
!       LL       :  Constants of integration in Eq. SC(1), obtained
!                   by solving scaled version of Eq. SC(5);
!                   exponential term of Eq. SC(12) not included
            
!       LYRCUT   :  Logical flag for truncation of comput. layer
            
!       NN       :  Order of double-Gauss quadrature (NSTR/2)
            
!       NCUT     :  Number of computational layer where absorption
!                   optical depth exceeds ABSCUT
            
!       PRNTU0   :  TRUE, print azimuthally-averaged intensity at
!                   quadrature angles
            
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
            
!       UTAUPR   :  Optical depths of user output levels in delta-M
!                   coordinates;  equal to UTAU if no delta-M
            
!       XR0      :  Expansion of thermal source function in Eq. SS(14),
!                   STWL(24c)
            
!       XR1      :  Expansion of thermal source function Eq. SS(16),
!                   STWL(24c)
            
!       ZZ       :  Beam source vectors in Eq. SS(19), STWL(24b)
            
!       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16),
!                   Y0 in STWL(26b)
            
!       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16),
!                   Y1 in STWL(26a)
            
!       (remainder are DISORT input variables)
            
            
!    O U T P U T     V A R I A B L E S:
            
!       U0C      :  Azimuthally averaged intensities
!                   ( at polar quadrature angles )
            
!       (RFLDIR, RFLDN, FLUP, DFDT, UAVG are DISORT output variables)
            
            
!    I N T E R N A L       V A R I A B L E S:
            
!       DIRINT   :  Direct intensity attenuated
!       FDNTOT   :  Total downward flux (direct + diffuse)
!       FLDIR    :  Direct-beam flux (delta-M scaled)
!       FLDN     :  Diffuse down-flux (delta-M scaled)
!       FNET     :  Net flux (total-down - diffuse-up)
!       FACT     :  EXP( - UTAUPR / UMU0 )
!       PLSORC   :  Planck source function (thermal)
!       ZINT     :  Intensity of m = 0 case, in Eq. SC(1)
            
!   Called by- DISORT
!   Calls- ZEROIT
! +-------------------------------------------------------------------+
            
!     .. Scalar Arguments ..
            
            
!sergio comment out maxulv
!      INTEGER   MAXULV, MXCMU, MXULV, NCUT, NN, NSTR, NTAU
            
            
!     ..
!     .. Array Arguments ..
            
            
            
            
!     ..
!     .. Local Scalars ..
            
            INTEGER :: IQ, JQ, LU, LYU
            REAL :: ANG1, ANG2, DIRINT, FACT, FDNTOT, FNET, PLSORC, ZINT
!     ..
!     .. External Subroutines ..
            
            EXTERNAL  ZEROIT
!     ..
!     .. Intrinsic Functions ..
            
            INTRINSIC EXP
!     ..
            
            
            IF( PRNT( 2 ) ) WRITE ( KSTDWARN, '(//,21X,A,/,2A,/,2A,/)' )  &
                '<----------------------- FLUXES ----------------------->',  &
                '   Optical  Compu    Downward    Downward    Downward     ',  &
                ' Upward                    Mean      Planck   d(Net Flux)',  &
                '     Depth  Layer      Direct     Diffuse       Total     ',  &
                'Diffuse         Net   Intensity      Source   / d(Op Dep)'
            
!                                        ** Zero DISORT output arrays
            CALL ZEROIT( U0C, MXULV*MXCMU )
            CALL ZEROIT( FLDIR, MXULV )
            CALL ZEROIT( FLDN, MXULV )
            
!                                        ** Loop over user levels
            DO  LU = 1, NTAU
              
              LYU  = LAYRU( LU )
              
              IF( LYRCUT .AND. LYU > NCUT ) THEN
!                                                ** No radiation reaches
!                                                ** this level
                FDNTOT = 0.0
                FNET   = 0.0
                PLSORC = 0.0
                GO TO 70
                
              END IF
              
              
              IF( FBEAM > 0.0 ) THEN
                
                FACT         = EXP( -UTAUPR( LU ) / UMU0 )
                DIRINT       = FBEAM*FACT
                FLDIR( LU )  = UMU0*( FBEAM*FACT )
                RFLDIR( LU ) = UMU0*FBEAM * EXP( -UTAU( LU ) / UMU0 )
                
              ELSE
                
                DIRINT       = 0.0
                FLDIR( LU )  = 0.0
                RFLDIR( LU ) = 0.0
                
              END IF
              
              
              DO  IQ = 1, NN
                
                ZINT = 0.0
                
                DO  JQ = 1, NN
                  ZINT = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*  &
                      EXP( -KK( JQ,LYU )*( UTAUPR( LU ) - TAUCPR( LYU ) ) )
                END DO
                
                DO  JQ = NN + 1, NSTR
                  ZINT = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*  &
                      EXP( -KK( JQ,LYU )*( UTAUPR( LU ) - TAUCPR( LYU-1 ) ) )
                END DO
                
                U0C( IQ, LU ) = ZINT
                
                IF( FBEAM > 0.0 ) U0C( IQ, LU ) = ZINT + ZZ( IQ, LYU )*FACT
                
                U0C( IQ, LU ) = U0C( IQ, LU ) + ZPLK0( IQ,LYU ) +  &
                    ZPLK1( IQ,LYU )*UTAUPR( LU )
                UAVG( LU ) = UAVG( LU ) + CWT( NN + 1 - IQ )*U0C( IQ, LU )
                FLDN( LU ) = FLDN( LU ) + CWT( NN + 1 - IQ )*  &
                    CMU( NN + 1 - IQ )*U0C( IQ, LU )
              END DO
              
              
              DO  IQ = NN + 1, NSTR
                
                ZINT = 0.0
                
                DO  JQ = 1, NN
                  ZINT = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*  &
                      EXP( -KK( JQ,LYU )*( UTAUPR( LU ) - TAUCPR( LYU ) ) )
                END DO
                
                DO  JQ = NN + 1, NSTR
                  ZINT = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*  &
                      EXP( -KK( JQ,LYU )*( UTAUPR( LU ) - TAUCPR( LYU-1 ) ) )
                END DO
                
                U0C( IQ, LU ) = ZINT
                
                IF( FBEAM > 0.0 ) U0C( IQ, LU ) = ZINT + ZZ( IQ, LYU )*FACT
                
                U0C( IQ, LU ) = U0C( IQ, LU ) + ZPLK0( IQ,LYU ) +  &
                    ZPLK1( IQ,LYU )*UTAUPR( LU )
                UAVG( LU ) = UAVG( LU ) + CWT( IQ - NN )*U0C( IQ, LU )
                FLUP( LU ) = FLUP( LU ) + CWT( IQ - NN )*CMU( IQ - NN )*  &
                    U0C( IQ, LU )
              END DO
              
              
              FLUP( LU )  = 2.*PI*FLUP( LU )
              FLDN( LU )  = 2.*PI*FLDN( LU )
              FDNTOT      = FLDN( LU ) + FLDIR( LU )
              FNET        = FDNTOT - FLUP( LU )
              RFLDN( LU ) = FDNTOT - RFLDIR( LU )
              UAVG( LU )  = ( 2.*PI*UAVG( LU ) + DIRINT ) / ( 4.*PI )
              PLSORC      = XR0( LYU ) + XR1( LYU )*UTAUPR( LU )
              DFDT( LU )  = ( 1. - SSALB( LYU ) ) * 4.*PI *  &
                  ( UAVG( LU ) - PLSORC )
              
              70    CONTINUE
              IF( PRNT( 2 ) ) WRITE ( KSTDWARN, '(F10.4,I7,1P,7E12.3,E14.3)' )  &
                  UTAU( LU ), LYU, RFLDIR( LU ), RFLDN( LU ), FDNTOT,  &
                  FLUP( LU ), FNET, UAVG( LU ), PLSORC, DFDT( LU )
              
            END DO
            
            
            IF( PRNTU0 ) THEN
              
              WRITE ( KSTDWARN, '(//,2A)' ) ' ******** AZIMUTHALLY AVERAGED ',  &
                  'INTENSITIES ( at polar quadrature angles ) *******'
              
              DO  LU = 1, NTAU
                
                WRITE ( KSTDWARN, '(/,A,F10.4,//,2A)' )  &
                    ' Optical depth =', UTAU( LU ),  &
                    '     Angle (deg)   cos(Angle)     Intensity',  &
                    '     Angle (deg)   cos(Angle)     Intensity'
                
                DO  IQ = 1, NN
                  ANG1 = ( 180./PI )*ACOS( CMU( 2 *NN-IQ+1 ) )
                  ANG2 = ( 180./PI )*ACOS( CMU( IQ ) )
                  WRITE ( KSTDWARN, '(2(0P,F16.4,F13.5,1P,E14.3))' )  &
                      ANG1, CMU(2*NN-IQ+1), U0C(IQ,LU),  &
                      ANG2, CMU(IQ),        U0C(IQ+NN,LU)
                END DO
                
              END DO
              
            END IF
            
            
            RETURN
          END SUBROUTINE FLUXES
          
          SUBROUTINE INTCOR( DITHER, FBEAM, FLYR, LAYRU, LYRCUT, MAXMOM,  &
              MAXULV, MAXUMU, NMOM, NCUT, NPHI, NSTR, NTAU,  &
              NUMU, OPRIM, PHASA, PHAST, PHASM, PHIRAD, PI,  &
              RPD, PMOM, SSALB, DTAUC, TAUC, TAUCPR, UMU, UMU0, UTAU, UTAUPR, UU )
          
!       Corrects intensity field by using Nakajima-Tanaka algorithm
!       (1988). For more details, see Section 3.6 of STWL NASA report.
          
!                I N P U T   V A R I A B L E S
          
!       DITHER  10 times machine precision
          
!       DTAUC   computational-layer optical depths
          
!       FBEAM   incident beam radiation at top
          
!       FLYR    separated fraction in delta-M method
          
!       LAYRU   index of UTAU in multi-layered system
          
!       LYRCUT  logical flag for truncation of computational layer
          
!       NMOM    number of phase function Legendre coefficients supplied
          
!       NCUT    total number of computational layers considered
          
!       NPHI    number of user azimuthal angles
          
!       NSTR    number of polar quadrature angles
          
!       NTAU    number of user-defined optical depths
          
!       NUMU    number of user polar angles
          
!       OPRIM   delta-M-scaled single-scatter albedo
          
!       PHIRAD  azimuthal angles in radians
          
!       PMOM    phase function Legendre coefficients (K, LC)
!                   K = 0 to NMOM, LC = 1 to NLYR with PMOM(0,LC)=1
          
!       RPD     PI/180
          
!       SSALB   single scattering albedo at computational layers
          
!       TAUC    optical thickness at computational levels
          
!       TAUCPR  delta-M-scaled optical thickness
          
!       UMU     cosine of emergent angle
          
!       UMU0    cosine of incident zenith angle
          
!       UTAU    user defined optical depths
          
!       UTAUPR  delta-M-scaled version of UTAU
          
!                O U T P U T   V A R I A B L E S
          
!       UU      corrected intensity field; UU(IU,LU,J)
!                         IU=1,NUMU; LU=1,NTAU; J=1,NPHI
          
!                I N T E R N A L   V A R I A B L E S
          
!       CTHETA  cosine of scattering angle
!       DTHETA  angle (degrees) to define aureole region as
!                    direction of beam source +/- DTHETA
!       PHASA   actual (exact) phase function
!       PHASM   delta-M-scaled phase function
!       PHAST   phase function used in TMS correction; actual phase
!                    function divided by (1-FLYR*SSALB)
!       PL      ordinary Legendre polynomial of degree l, P-sub-l
!       PLM1    ordinary Legendre polynomial of degree l-1, P-sub-(l-1)
!       PLM2    ordinary Legendre polynomial of degree l-2, P-sub-(l-2)
!       THETA0  incident zenith angle (degrees)
!       THETAP  emergent angle (degrees)
!       USSNDM  single-scattered intensity computed by using exact
!                   phase function and scaled optical depth
!                   (first term in STWL(68a))
!       USSP    single-scattered intensity from delta-M method
!                   (second term in STWL(68a))
!       DUIMS   intensity correction term from IMS method
!                   (delta-I-sub-IMS in STWL(A.19))
          
!   Called by- DISORT
!   Calls- SINSCA, SECSCA
          
! +-------------------------------------------------------------------+
          
!     .. Scalar Arguments ..
          
          
          REAL, INTENT(IN)                         :: DITHER
          REAL, INTENT(IN OUT)                     :: FBEAM
          REAL, INTENT(IN)                         :: FLYR( * )
          INTEGER, INTENT(IN)                      :: LAYRU( * )
          LOGICAL, INTENT(IN OUT)                  :: LYRCUT
          INTEGER, INTENT(IN)                      :: MAXMOM
          INTEGER, INTENT(IN OUT)                  :: MAXULV
          INTEGER, INTENT(IN OUT)                  :: MAXUMU
          INTEGER, INTENT(IN)                      :: NMOM
          INTEGER, INTENT(IN)                      :: NCUT
          INTEGER, INTENT(IN)                      :: NPHI
          INTEGER, INTENT(IN)                      :: NSTR
          INTEGER, INTENT(IN)                      :: NTAU
          INTEGER, INTENT(IN)                      :: NUMU
          REAL, INTENT(IN OUT)                     :: OPRIM( * )
          REAL, INTENT(OUT)                        :: PHASA( * )
          REAL, INTENT(OUT)                        :: PHAST( * )
          REAL, INTENT(OUT)                        :: PHASM( * )
          REAL, INTENT(IN OUT)                     :: PHIRAD( * )
          REAL, INTENT(IN OUT)                     :: PI
          REAL, INTENT(IN)                         :: RPD
          REAL, INTENT(IN OUT)                     :: PMOM( 0:MAXMOM, * )
          REAL, INTENT(IN)                         :: SSALB( * )
          REAL, INTENT(IN OUT)                     :: DTAUC( * )
          REAL, INTENT(IN OUT)                     :: TAUC( 0:* )
          REAL, INTENT(IN OUT)                     :: TAUCPR( 0:* )
          REAL, INTENT(IN)                         :: UMU( * )
          REAL, INTENT(IN)                         :: UMU0
          REAL, INTENT(IN)                         :: UTAU( * )
          REAL, INTENT(IN OUT)                     :: UTAUPR( * )
          REAL, INTENT(OUT)                        :: UU( MAXUMU, MAXULV, * )
          
          
          
!     ..
!     .. Array Arguments ..
          
          
          
!     ..
!     .. Local Scalars ..
          
          INTEGER :: IU, JP, K, LC, LTAU, LU
          REAL :: CTHETA, DTHETA, DUIMS, PL, PLM1, PLM2, THETA0, THETAP,  &
              USSNDM, USSP
!     ..
!     .. External Functions ..
          
          REAL :: SECSCA, SINSCA
          EXTERNAL  SECSCA, SINSCA
!     ..
!     .. Intrinsic Functions ..
          
          INTRINSIC ABS, ACOS, COS, SQRT
!     ..
          
          
          DTHETA = 10.
          
!                                ** Start loop over zenith angles
          
          DO  IU = 1, NUMU
            
            IF( UMU( IU ) < 0. ) THEN
              
!                                ** Calculate zenith angles of icident
!                                ** and emerging directions
              
              THETA0 = ACOS( -UMU0 ) / RPD
              THETAP = ACOS( UMU( IU ) ) / RPD
              
            END IF
            
!                                ** Start loop over azimuth angles
            
            DO  JP = 1, NPHI
              
!                                ** Calculate cosine of scattering
!                                ** angle, Eq. STWL(4)
              
              CTHETA = -UMU0*UMU( IU ) + SQRT( ( 1.-UMU0**2 )*  &
                  ( 1.-UMU( IU )**2 ) )*COS( PHIRAD( JP ) )
              
!                                ** Initialize phase function
              DO  LC = 1, NCUT
                
                PHASA( LC ) = 1.
                PHASM( LC ) = 1.
                
              END DO
!                                ** Initialize Legendre poly. recurrence
              PLM1 = 1.
              PLM2 = 0.
              
              DO  K = 1, NMOM
!                                ** Calculate Legendre polynomial of
!                                ** P-sub-l by upward recurrence
                
                PL   = ( ( 2 *K-1 )*CTHETA*PLM1 - ( K-1 )*PLM2 ) / K
                PLM2 = PLM1
                PLM1 = PL
!                                ** Calculate actual phase function
                DO  LC = 1, NCUT
                  
                  PHASA( LC ) = PHASA( LC ) + ( 2*K + 1 )*PL*PMOM( K, LC )
                  
                END DO
                
!                                ** Calculate delta-M transformed
!                                ** phase function
                IF( K <= NSTR - 1 ) THEN
                  
                  DO  LC = 1, NCUT
                    
                    PHASM( LC ) = PHASM( LC ) + ( 2*K + 1 ) * PL *  &
                        ( PMOM( K,LC ) - FLYR( LC ) ) / ( 1. - FLYR( LC ) )
                  END DO
                  
                END IF
                
              END DO
              
              
!                                ** Apply TMS method, Eq. STWL(68)
              DO  LC = 1, NCUT
                
                PHAST( LC ) = PHASA(LC) / ( 1. - FLYR(LC) * SSALB(LC) )
                
              END DO
              
              DO  LU = 1, NTAU
                
                IF( .NOT.LYRCUT .OR. LAYRU( LU ) < NCUT ) THEN
                  
                  USSNDM  = SINSCA( DITHER, LAYRU( LU ), NCUT, PHAST,  &
                      SSALB, TAUCPR, UMU( IU ), UMU0, UTAUPR( LU ), FBEAM, PI )
                  
                  USSP    = SINSCA( DITHER, LAYRU( LU ), NCUT, PHASM,  &
                      OPRIM, TAUCPR, UMU( IU ), UMU0, UTAUPR( LU ), FBEAM, PI )
                  
                  UU( IU, LU, JP ) = UU( IU, LU, JP ) + USSNDM - USSP
                  
                END IF
                
              END DO
              
              IF( UMU(IU) < 0. .AND. ABS( THETA0-THETAP ) <= DTHETA) THEN
                
!                                ** Emerging direction is in the aureole
!                                ** (theta0 +/- dtheta). Apply IMS
!                                ** method for correction of secondary
!                                ** scattering below top level.
                
                LTAU = 1
                
                IF( UTAU( 1 ) <= DITHER ) LTAU = 2
                
                DO  LU = LTAU, NTAU
                  
                  IF( .NOT.LYRCUT .OR. LAYRU( LU ) < NCUT ) THEN
                    
                    DUIMS = SECSCA( CTHETA, FLYR, LAYRU( LU ), MAXMOM,  &
                        NMOM, NSTR, PMOM, SSALB, DTAUC,  &
                        TAUC, UMU( IU ), UMU0, UTAU( LU ), FBEAM, PI )
                    
                    UU( IU, LU, JP ) = UU( IU, LU, JP ) - DUIMS
                    
                  END IF
                  
                END DO
                
              END IF
!                                ** End loop over azimuth angles
            END DO
            
!                                ** End loop over zenith angles
          END DO
          
          
          RETURN
        END SUBROUTINE INTCOR
!************************************************************************
        
        REAL FUNCTION  SECSCA( CTHETA, FLYR, LAYRU, MAXMOM, NMOM, NSTR,  &
            PMOM, SSALB, DTAUC, TAUC, UMU, UMU0, UTAU, FBEAM, PI )
        
!          Calculates secondary scattered intensity of EQ. STWL (A7)
        
!                I N P U T   V A R I A B L E S
        
!        CTHETA  cosine of scattering angle
        
!        DTAUC   computational-layer optical depths
        
!        FLYR    separated fraction f in Delta-M method
        
!        LAYRU   index of UTAU in multi-layered system
        
!        MAXMOM  maximum number of phase function moment coefficients
        
!        NMOM    number of phase function Legendre coefficients supplied
        
!        NSTR    number of polar quadrature angles
        
!        PMOM    phase function Legendre coefficients (K, LC)
!                K = 0 to NMOM, LC = 1 to NLYR, with PMOM(0,LC)=1
        
!        SSALB   single scattering albedo of computational layers
        
!        TAUC    cumulative optical depth at computational layers
        
!        UMU     cosine of emergent angle
        
!        UMU0    cosine of incident zenith angle
        
!        UTAU    user defined optical depth for output intensity
        
!        FBEAM   incident beam radiation at top
        
!        PI       3.1415...
        
!   LOCAL VARIABLES
        
!        PSPIKE  2*P"-P"**2, where P" is the residual phase function
!        WBAR    mean value of single scattering albedo
!        FBAR    mean value of separated fraction f
!        DTAU    layer optical depth
!        STAU    sum of layer optical depths between top of atmopshere
!                and layer LAYRU
        
!   Called by- INTCOR
!   Calls- XIFUNC
! +-------------------------------------------------------------------+
        
!     .. Scalar Arguments ..
        
        REAL, INTENT(IN)                         :: CTHETA
        REAL, INTENT(IN)                         :: FLYR( * )
        INTEGER, INTENT(IN)                      :: LAYRU
        INTEGER, INTENT(IN OUT)                  :: MAXMOM
        INTEGER, INTENT(IN)                      :: NMOM
        INTEGER, INTENT(IN)                      :: NSTR
        REAL, INTENT(IN)                         :: PMOM( 0:MAXMOM, * )
        REAL, INTENT(IN)                         :: SSALB( * )
        REAL, INTENT(IN)                         :: DTAUC( * )
        REAL, INTENT(IN)                         :: TAUC( 0:* )
        REAL, INTENT(IN OUT)                     :: UMU
        REAL, INTENT(IN)                         :: UMU0
        REAL, INTENT(IN)                         :: UTAU
        REAL, INTENT(IN)                         :: FBEAM
        REAL, INTENT(IN)                         :: PI
        
        
!     ..
!     .. Array Arguments ..
        
!     ..
!     .. Local Scalars ..
        INTEGER :: K, LYR
        REAL :: DTAU, FBAR, GBAR, PL, PLM1, PLM2, PSPIKE, STAU, UMU0P,  &
            WBAR, ZERO
!     ..
!     .. External Functions ..
        REAL :: XIFUNC
        EXTERNAL  XIFUNC
!     ..
        
        ZERO = 1E-4
        
!                          ** Calculate vertically averaged value of
!                          ** single scattering albedo and separated
!                          ** fraction f, Eq. STWL (A.15)
        
        DTAU = UTAU - TAUC( LAYRU - 1 )
        WBAR = SSALB( LAYRU ) * DTAU
        FBAR = FLYR( LAYRU ) * WBAR
        STAU = DTAU
        
        DO  LYR = 1, LAYRU - 1
          
          WBAR = WBAR + SSALB( LYR ) * DTAUC( LYR )
          FBAR = FBAR + SSALB( LYR ) * DTAUC( LYR ) * FLYR( LYR )
          STAU = STAU + DTAUC( LYR )
          
        END DO
        
        IF( WBAR <= ZERO .OR.  &
              FBAR <= ZERO .OR. STAU <= ZERO .OR.FBEAM <= ZERO ) THEN
          
          SECSCA = 0.0
          RETURN
          
        END IF
        
        FBAR  = FBAR / WBAR
        WBAR  = WBAR / STAU
        
        
!                          ** Calculate PSPIKE=(2P"-P"**2)
        PSPIKE = 1.
        GBAR   = 1.
        PLM1    = 1.
        PLM2    = 0.
!                                   ** PSPIKE for L<=2N-1
        DO  K = 1, NSTR - 1
          
          PL   = ( ( 2 *K-1 )*CTHETA*PLM1 - ( K-1 )*PLM2 ) / K
          PLM2  = PLM1
          PLM1  = PL
          
          PSPIKE = PSPIKE + ( 2.*GBAR - GBAR**2 )*( 2*K + 1 )*PL
          
        END DO
!                                   ** PSPIKE for L>2N-1
        DO  K = NSTR, NMOM
          
          PL   = ( ( 2 *K-1 )*CTHETA*PLM1 - ( K-1 )*PLM2 ) / K
          PLM2  = PLM1
          PLM1  = PL
          
          DTAU = UTAU - TAUC( LAYRU - 1 )
          
          GBAR = PMOM( K, LAYRU ) * SSALB( LAYRU ) * DTAU
          
          DO  LYR = 1, LAYRU - 1
            GBAR = GBAR + PMOM( K, LYR ) * SSALB( LYR ) * DTAUC( LYR )
          END DO
          
          IF( FBAR*WBAR*STAU <= ZERO ) THEN
            GBAR   = 0.0
          ELSE
            GBAR   = GBAR / ( FBAR*WBAR*STAU )
          END IF
          
          PSPIKE = PSPIKE + ( 2.*GBAR - GBAR**2 )*( 2*K + 1 )*PL
          
        END DO
        
        UMU0P = UMU0 / ( 1. - FBAR*WBAR )
        
!                              ** Calculate IMS correction term,
!                              ** Eq. STWL (A.13)
        
        SECSCA = FBEAM / ( 4.*PI ) * ( FBAR*WBAR )**2 / ( 1.-FBAR*WBAR ) *  &
            PSPIKE * XIFUNC( -UMU, UMU0P, UMU0P, UTAU )
        
        
        RETURN
      END FUNCTION  SECSCA
      
!************************************************************************
      
      SUBROUTINE SETDIS( CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM,  &
          FLYR, GL, IBCND, LAYRU, LYRCUT, MAXMOM, MAXUMU,  &
          MXCMU, NCUT, NLYR, NTAU, NN, NSTR, PLANK, NUMU,  &
          ONLYFL, CORINT, OPRIM, PMOM, SSALB, TAUC,  &
          TAUCPR, UTAU, UTAUPR, UMU, UMU0, USRTAU, USRANG )
      
!          Perform miscellaneous setting-up operations
      
!    INPUT :  all are DISORT input variables (see DOC file)
      
      
!    O U T P U T     V A R I A B L E S:
      
!       NTAU,UTAU   if USRTAU = FALSE (defined in DISORT.doc)
!       NUMU,UMU    if USRANG = FALSE (defined in DISORT.doc)
      
!       CMU,CWT     computational polar angles and
!                   corresponding quadrature weights
      
!       EXPBEA      transmission of direct beam
      
!       FLYR        separated fraction in delta-M method
      
!       GL          phase function Legendre coefficients multiplied
!                   by (2L+1) and single-scatter albedo
      
!       LAYRU       Computational layer in which UTAU falls
      
!       LYRCUT      flag as to whether radiation will be zeroed
!                   below layer NCUT
      
!       NCUT        computational layer where absorption
!                   optical depth first exceeds  ABSCUT
      
!       NN          NSTR / 2
      
!       OPRIM       delta-M-scaled single-scatter albedo
      
!       TAUCPR      delta-M-scaled optical depth
      
!       UTAUPR      delta-M-scaled version of  UTAU
      
!   Called by- DISORT
!   Calls- QGAUSN, ERRMSG
! ---------------------------------------------------------------------
      
!     .. Scalar Arguments ..
      
      
      REAL, INTENT(OUT)                        :: CMU( MXCMU )
      REAL, INTENT(OUT)                        :: CWT( MXCMU )
      LOGICAL, INTENT(IN)                      :: DELTAM
      REAL, INTENT(IN)                         :: DTAUC( * )
      REAL, INTENT(OUT)                        :: DTAUCP( * )
      REAL, INTENT(OUT)                        :: EXPBEA( 0:* )
      REAL, INTENT(IN)                         :: FBEAM
      REAL, INTENT(OUT)                        :: FLYR( * )
      REAL, INTENT(OUT)                        :: GL( 0:MXCMU, * )
      INTEGER, INTENT(IN)                      :: IBCND
      INTEGER, INTENT(OUT)                     :: LAYRU( * )
      LOGICAL, INTENT(OUT)                     :: LYRCUT
      INTEGER, INTENT(IN OUT)                  :: MAXMOM
      INTEGER, INTENT(IN OUT)                  :: MAXUMU
      INTEGER, INTENT(IN OUT)                  :: MXCMU
      INTEGER, INTENT(OUT)                     :: NCUT
      INTEGER, INTENT(IN OUT)                  :: NLYR
      INTEGER, INTENT(OUT)                     :: NTAU
      INTEGER, INTENT(OUT)                     :: NN
      INTEGER, INTENT(IN)                      :: NSTR
      LOGICAL, INTENT(IN OUT)                  :: PLANK
      INTEGER, INTENT(OUT)                     :: NUMU
      LOGICAL, INTENT(IN)                      :: ONLYFL
      LOGICAL, INTENT(OUT)                     :: CORINT
      REAL, INTENT(OUT)                        :: OPRIM( * )
      REAL, INTENT(OUT)                        :: PMOM( 0:MAXMOM, * )
      REAL, INTENT(IN)                         :: SSALB( * )
      REAL, INTENT(IN)                         :: TAUC( 0:* )
      REAL, INTENT(OUT)                        :: TAUCPR( 0:* )
      REAL, INTENT(OUT)                        :: UTAU( * )
      REAL, INTENT(OUT)                        :: UTAUPR( * )
      REAL, INTENT(OUT)                        :: UMU( MAXUMU )
      REAL, INTENT(IN)                         :: UMU0
      LOGICAL, INTENT(IN OUT)                  :: USRTAU
      LOGICAL, INTENT(IN)                      :: USRANG
      
      
      
!     ..
!     .. Array Arguments ..
      
      
      
!     ..
!     .. Local Scalars ..
      
      INTEGER :: IQ, IU, K, LC, LU
      REAL :: ABSCUT, ABSTAU, F, YESSCT
!     ..
!     .. External Subroutines ..
      
      EXTERNAL  ERRMSG, QGAUSN
!     ..
!     .. Intrinsic Functions ..
      
      INTRINSIC ABS, EXP
!     ..
      DATA      ABSCUT / 10. /
      
      
      IF( .NOT.USRTAU ) THEN
!                              ** Set output levels at computational
!                              ** layer boundaries
        NTAU  = NLYR + 1
        
        DO  LC = 0, NTAU - 1
          UTAU( LC + 1 ) = TAUC( LC )
        END DO
        
      END IF
!                        ** Apply delta-M scaling and move description
!                        ** of computational layers to local variables
      EXPBEA( 0 ) = 1.0
      TAUCPR( 0 ) = 0.0
      ABSTAU      = 0.0
      YESSCT      = 0.0
      
      DO  LC = 1, NLYR
        
        YESSCT = YESSCT + SSALB( LC )
        
        PMOM( 0, LC ) = 1.0
        
        IF( ABSTAU < ABSCUT ) NCUT  = LC
        
        ABSTAU = ABSTAU + ( 1. - SSALB( LC ) )*DTAUC( LC )
        
        IF( .NOT.DELTAM ) THEN
          
          OPRIM( LC )  = SSALB( LC )
          DTAUCP( LC ) = DTAUC( LC )
          TAUCPR( LC ) = TAUC( LC )
          
          DO  K = 0, NSTR - 1
            GL( K, LC ) = ( 2*K + 1 )*OPRIM( LC )*PMOM( K, LC )
          END DO
          
          F  = 0.0
          
          
        ELSE
!                                    ** Do delta-M transformation
          
          F  = PMOM( NSTR, LC )
          OPRIM( LC )  = SSALB( LC )*( 1. - F ) / ( 1. - F*SSALB(LC) )
          DTAUCP( LC ) = ( 1. - F*SSALB( LC ) )*DTAUC( LC )
          TAUCPR( LC ) = TAUCPR( LC - 1 ) + DTAUCP( LC )
          
          DO  K = 0, NSTR - 1
            GL( K, LC ) = ( 2*K + 1 )*OPRIM( LC )*  &
                ( PMOM( K,LC ) - F ) / ( 1. - F )
          END DO
          
        END IF
        
        FLYR( LC ) = F
        EXPBEA( LC ) = 0.0
        
        IF( FBEAM > 0.0 ) EXPBEA( LC ) = EXP( -TAUCPR( LC )/UMU0 )
        
      END DO
!                      ** If no thermal emission, cut off medium below
!                      ** absorption optical depth = ABSCUT ( note that
!                      ** delta-M transformation leaves absorption
!                      ** optical depth invariant ).  Not worth the
!                      ** trouble for one-layer problems, though.
      LYRCUT = .FALSE.
      
      IF( ABSTAU >= ABSCUT .AND. .NOT.PLANK .AND. IBCND /= 1 .AND.  &
          NLYR > 1 ) LYRCUT = .TRUE.
      
      IF( .NOT.LYRCUT ) NCUT = NLYR
      
!                             ** Set arrays defining location of user
!                             ** output levels within delta-M-scaled
!                             ** computational mesh
      DO  LU = 1, NTAU
        
        DO  LC = 1, NLYR
          
          IF( UTAU( LU ) >= TAUC( LC-1 ) .AND.  &
              UTAU( LU ) <= TAUC( LC ) ) GO TO 60
          
        END DO
        LC   = NLYR
        
        60    CONTINUE
        UTAUPR( LU ) = UTAU( LU )
        IF( DELTAM ) UTAUPR( LU ) = TAUCPR( LC - 1 ) +  &
            ( 1. - SSALB( LC )*FLYR( LC ) )* ( UTAU( LU ) - TAUC( LC-1 ) )
        LAYRU( LU ) = LC
        
      END DO
!                      ** Calculate computational polar angle cosines
!                      ** and associated quadrature weights for Gaussian
!                      ** quadrature on the interval (0,1) (upward)
      NN   = NSTR / 2
      
      CALL QGAUSN( NN, CMU, CWT )
!                                  ** Downward (neg) angles and weights
      DO  IQ = 1, NN
        CMU( IQ + NN ) = -CMU( IQ )
        CWT( IQ + NN ) = CWT( IQ )
      END DO
      
      
      IF( FBEAM > 0.0 ) THEN
!                               ** Compare beam angle to comput. angles
        DO  IQ = 1, NN
          
          IF( ABS( UMU0-CMU( IQ ) )/UMU0 < 1.E-4 ) CALL ERRMSG(  &
              'SETDIS--beam angle=computational angle; change NSTR', .True. )
          
        END DO
        
      END IF
      
      
      IF( .NOT.USRANG .OR. ( ONLYFL.AND.MAXUMU >= NSTR ) ) THEN
        
!                                   ** Set output polar angles to
!                                   ** computational polar angles
        NUMU = NSTR
        
        DO  IU = 1, NN
          UMU( IU ) = -CMU( NN + 1 - IU )
        END DO
        
        DO  IU = NN + 1, NSTR
          UMU( IU ) = CMU( IU - NN )
        END DO
        
      END IF
      
      
      IF( USRANG .AND. IBCND == 1 ) THEN
        
!                               ** Shift positive user angle cosines to
!                               ** upper locations and put negatives
!                               ** in lower locations
        DO  IU = 1, NUMU
          UMU( IU + NUMU ) = UMU( IU )
        END DO
        
        DO  IU = 1, NUMU
          UMU( IU ) = -UMU( 2*NUMU + 1 - IU )
        END DO
        
        NUMU = 2*NUMU
        
      END IF
      
!                               ** Turn off intensity correction when
!                               ** only fluxes are calculated, there
!                               ** is no beam source, no scattering,
!                               ** or delta-M transformation is not
!                               ** applied
      
      IF( ONLYFL .OR. FBEAM == 0.0 .OR. YESSCT == 0.0 .OR.  &
          .NOT.DELTAM )  CORINT = .FALSE.
      
      
      RETURN
    END SUBROUTINE SETDIS
    
!************************************************************************
    
    SUBROUTINE SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,  &
        LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT,  &
        NNLYRI, NN, NSTR, TAUCPR, WK )
    
!        Calculate coefficient matrix for the set of equations
!        obtained from the boundary conditions and the continuity-
!        of-intensity-at-layer-interface equations;  store in the
!        special banded-matrix format required by LINPACK routines
    
    
!    I N P U T      V A R I A B L E S:
    
!       BDR      :  surface bidirectional reflectivity
    
!       CMU,CWT     abscissae, weights for Gauss quadrature
!                   over angle cosine
    
!       DELM0    :  Kronecker delta, delta-sub-m0
    
!       GC       :  Eigenvectors at polar quadrature angles, SC(1)
    
!       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7), STWL(23b)
    
!       LYRCUT   :  Logical flag for truncation of computational layers
    
!       NN       :  Number of streams in a hemisphere (NSTR/2)
    
!       NCUT     :  Total number of computational layers considered
    
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
    
!       (remainder are DISORT input variables)
    
    
!   O U T P U T     V A R I A B L E S:
    
!       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
!                   scaled by Eq. SC(12); in banded form required
!                   by LINPACK solution routines
    
!       NCOL     :  Number of columns in CBAND
    
    
!   I N T E R N A L    V A R I A B L E S:
    
!       IROW     :  Points to row in CBAND
!       JCOL     :  Points to position in layer block
!       LDA      :  Row dimension of CBAND
!       NCD      :  Number of diagonals below or above main diagonal
!       NSHIFT   :  For positioning number of rows in band storage
!       WK       :  Temporary storage for EXP evaluations
    
    
!   BAND STORAGE
    
!      LINPACK requires band matrices to be input in a special
!      form where the elements of each diagonal are moved up or
!      down (in their column) so that each diagonal becomes a row.
!      (The column locations of diagonal elements are unchanged.)
    
!      Example:  if the original matrix is
    
!          11 12 13  0  0  0
!          21 22 23 24  0  0
!           0 32 33 34 35  0
!           0  0 43 44 45 46
!           0  0  0 54 55 56
!           0  0  0  0 65 66
    
!      then its LINPACK input form would be:
    
!           *  *  *  +  +  +  , * = not used
!           *  * 13 24 35 46  , + = used for pivoting
!           * 12 23 34 45 56
!          11 22 33 44 55 66
!          21 32 43 54 65  *
    
!      If A is a band matrix, the following program segment
!      will convert it to the form (ABD) required by LINPACK
!      band-matrix routines:
    
!               N  = (column dimension of A, ABD)
!               ML = (band width below the diagonal)
!               MU = (band width above the diagonal)
!               M = ML + MU + 1
!               DO J = 1, N
!                  I1 = MAX(1, J-MU)
!                  I2 = MIN(N, J+ML)
!                  DO I = I1, I2
!                     K = I - J + M
!                     ABD(K,J) = A(I,J)
!                  END DO
!               END DO
    
!      This uses rows  ML+1  through  2*ML+MU+1  of ABD.
!      The total number of rows needed in ABD is  2*ML+MU+1 .
!      In the example above, N = 6, ML = 1, MU = 2, and the
!      row dimension of ABD must be >= 5.
    
    
!   Called by- DISORT, ALBTRN
!   Calls- ZEROIT
! +-------------------------------------------------------------------+
    
!     .. Scalar Arguments ..
    
    
    REAL, INTENT(IN)                         :: BDR( MI, 0:MI )
    REAL, INTENT(OUT)                        :: CBAND( MI9M2, NNLYRI )
    REAL, INTENT(IN)                         :: CMU( MXCMU )
    REAL, INTENT(IN)                         :: CWT( MXCMU )
    REAL, INTENT(IN)                         :: DELM0
    REAL, INTENT(IN)                         :: DTAUCP( * )
    REAL, INTENT(IN)                         :: GC( MXCMU, MXCMU, * )
    REAL, INTENT(IN)                         :: KK( MXCMU, * )
    LOGICAL, INTENT(IN)                      :: LAMBER
    LOGICAL, INTENT(IN)                      :: LYRCUT
    INTEGER, INTENT(IN OUT)                  :: MI
    INTEGER, INTENT(IN OUT)                  :: MI9M2
    INTEGER, INTENT(IN OUT)                  :: MXCMU
    INTEGER, INTENT(OUT)                     :: NCOL
    INTEGER, INTENT(IN)                      :: NCUT
    INTEGER, INTENT(IN OUT)                  :: NNLYRI
    INTEGER, INTENT(IN)                      :: NN
    INTEGER, INTENT(IN OUT)                  :: NSTR
    REAL, INTENT(IN)                         :: TAUCPR( 0:* )
    REAL, INTENT(OUT)                        :: WK( MXCMU )
    
    
    
!     ..
!     .. Array Arguments ..
    
    
!     ..
!     .. Local Scalars ..
    
    INTEGER :: IQ, IROW, JCOL, JQ, K, LC, LDA, NCD, NNCOL, NSHIFT
    REAL :: EXPA, SUM
!     ..
!     .. External Subroutines ..
    
    EXTERNAL  ZEROIT
!     ..
!     .. Intrinsic Functions ..
    
    INTRINSIC EXP
!     ..
    
    
    CALL ZEROIT( CBAND, MI9M2*NNLYRI )
    
    NCD    = 3*NN - 1
    LDA    = 3*NCD + 1
    NSHIFT = LDA - 2*NSTR + 1
    NCOL   = 0
!                         ** Use continuity conditions of Eq. STWJ(17)
!                         ** to form coefficient matrix in STWJ(20);
!                         ** employ scaling transformation STWJ(22)
    DO  LC = 1, NCUT
      
      DO  IQ = 1, NN
        WK( IQ ) = EXP( KK( IQ,LC )*DTAUCP( LC ) )
      END DO
      
      JCOL  = 0
      
      DO  IQ = 1, NN
        
        NCOL  = NCOL + 1
        IROW  = NSHIFT - JCOL
        
        DO  JQ = 1, NSTR
          CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )
          CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )*WK( IQ )
          IROW  = IROW + 1
        END DO
        
        JCOL  = JCOL + 1
        
      END DO
      
      
      DO  IQ = NN + 1, NSTR
        
        NCOL  = NCOL + 1
        IROW  = NSHIFT - JCOL
        
        DO  JQ = 1, NSTR
          CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )* WK( NSTR + 1 - IQ )
          CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )
          IROW  = IROW + 1
        END DO
        
        JCOL  = JCOL + 1
        
      END DO
      
    END DO
!                  ** Use top boundary condition of STWJ(20a) for
!                  ** first layer
    JCOL  = 0
    
    DO  IQ = 1, NN
      
      EXPA  = EXP( KK( IQ,1 )*TAUCPR( 1 ) )
      IROW  = NSHIFT - JCOL + NN
      
      DO  JQ = NN, 1, -1
        CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )*EXPA
        IROW  = IROW + 1
      END DO
      
      JCOL  = JCOL + 1
      
    END DO
    
    
    DO  IQ = NN + 1, NSTR
      
      IROW  = NSHIFT - JCOL + NN
      
      DO  JQ = NN, 1, -1
        CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )
        IROW  = IROW + 1
      END DO
      
      JCOL  = JCOL + 1
      
    END DO
!                           ** Use bottom boundary condition of
!                           ** STWJ(20c) for last layer
    
    NNCOL = NCOL - NSTR
    JCOL  = 0
    
    DO  IQ = 1, NN
      
      NNCOL  = NNCOL + 1
      IROW   = NSHIFT - JCOL + NSTR
      
      DO  JQ = NN + 1, NSTR
        
        IF( LYRCUT .OR. ( LAMBER .AND. DELM0 == 0 ) ) THEN
          
!                          ** No azimuthal-dependent intensity if Lam-
!                          ** bert surface; no intensity component if
!                          ** truncated bottom layer
          
          CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )
          
        ELSE
          
          SUM  = 0.0
          
          DO  K = 1, NN
            SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )*  &
                GC( NN + 1 - K, IQ, NCUT )
          END DO
          
          CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT ) - ( 1.+ DELM0 )*SUM
        END IF
        
        IROW  = IROW + 1
        
      END DO
      
      JCOL  = JCOL + 1
      
    END DO
    
    
    DO  IQ = NN + 1, NSTR
      
      NNCOL  = NNCOL + 1
      IROW   = NSHIFT - JCOL + NSTR
      EXPA   = WK( NSTR + 1 - IQ )
      
      DO  JQ = NN + 1, NSTR
        
        IF( LYRCUT .OR. ( LAMBER .AND. DELM0 == 0 ) ) THEN
          
          CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )*EXPA
          
        ELSE
          
          SUM  = 0.0
          
          DO  K = 1, NN
            SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )*  &
                GC( NN + 1 - K, IQ, NCUT )
          END DO
          
          CBAND( IROW, NNCOL ) = ( GC( JQ,IQ,NCUT ) - ( 1.+ DELM0 )*SUM )*EXPA
        END IF
        
        IROW  = IROW + 1
        
      END DO
      
      JCOL  = JCOL + 1
      
    END DO
    
    
    RETURN
  END SUBROUTINE SETMTX
  
!************************************************************************
  
  REAL FUNCTION  SINSCA( DITHER, LAYRU, NLYR, PHASE, OMEGA, TAU,  &
      UMU, UMU0, UTAU, FBEAM, PI )
  
!        Calculates single-scattered intensity from EQS. STWL (65b,d,e)
  
!                I N P U T   V A R I A B L E S
  
!        DITHER   10 times machine precision
  
!        LAYRU    index of UTAU in multi-layered system
  
!        NLYR     number of sublayers
  
!        PHASE    phase functions of sublayers
  
!        OMEGA    single scattering albedos of sublayers
  
!        TAU      optical thicknesses of sublayers
  
!        UMU      cosine of emergent angle
  
!        UMU0     cosine of incident zenith angle
  
!        UTAU     user defined optical depth for output intensity
  
!        FBEAM   incident beam radiation at top
  
!        PI       3.1415...
  
!   Called by- INTCOR
! +-------------------------------------------------------------------+
  
!     .. Scalar Arguments ..
  
  
  REAL, INTENT(IN)                         :: DITHER
  INTEGER, INTENT(IN)                      :: LAYRU
  INTEGER, INTENT(IN)                      :: NLYR
  REAL, INTENT(IN)                         :: PHASE( * )
  REAL, INTENT(IN)                         :: OMEGA( * )
  REAL, INTENT(IN)                         :: TAU( 0:* )
  REAL, INTENT(IN)                         :: UMU
  REAL, INTENT(IN)                         :: UMU0
  REAL, INTENT(IN)                         :: UTAU
  REAL, INTENT(IN)                         :: FBEAM
  REAL, INTENT(IN)                         :: PI
  
  
!     ..
!     .. Array Arguments ..
  
  
!     ..
!     .. Local Scalars ..
  
  INTEGER :: LYR
  REAL :: EXP0, EXP1
!     ..
!     .. Intrinsic Functions ..
  
  INTRINSIC ABS, EXP
!     ..
  
  
  SINSCA = 0.
  EXP0 = EXP( -UTAU/UMU0 )
  
  IF( ABS( UMU+UMU0 ) <= DITHER ) THEN
    
!                                 ** Calculate downward intensity when
!                                 ** UMU=UMU0, Eq. STWL (65e)
    
    DO  LYR = 1, LAYRU - 1
      SINSCA = SINSCA + OMEGA( LYR ) * PHASE( LYR ) *  &
          ( TAU( LYR ) - TAU( LYR-1 ) )
    END DO
    
    SINSCA = FBEAM / ( 4.*PI * UMU0 ) * EXP0 * ( SINSCA +  &
        OMEGA( LAYRU )*PHASE( LAYRU )*( UTAU-TAU(LAYRU-1) ) )
    
    RETURN
    
  END IF
  
  
  IF( UMU > 0. ) THEN
!                                 ** Upward intensity, Eq. STWL (65b)
    DO  LYR = LAYRU, NLYR
      
      EXP1 = EXP( -( ( TAU( LYR )-UTAU )/UMU + TAU( LYR )/UMU0 ) )
      SINSCA = SINSCA + OMEGA( LYR )*PHASE( LYR )*( EXP0 - EXP1 )
      EXP0 = EXP1
      
    END DO
    
  ELSE
!                                 ** Downward intensity, Eq. STWL (65d)
    DO  LYR = LAYRU, 1, -1
      
      EXP1 = EXP( -( ( TAU(LYR-1)-UTAU )/UMU + TAU(LYR-1)/UMU0 ) )
      SINSCA = SINSCA + OMEGA( LYR )*PHASE( LYR )*( EXP0 - EXP1 )
      EXP0 = EXP1
      
    END DO
    
  END IF
  
  SINSCA = FBEAM / ( 4.*PI * ( 1. + UMU/UMU0 ) ) * SINSCA
  
  
  RETURN
END FUNCTION  SINSCA
!************************************************************************

SUBROUTINE SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MI, MAZIM,  &
    MXCMU, NN, NSTR, YLMC, CC, EVECC, EVAL, KK, GC, WKD )


REAL, INTENT(OUT)                        :: AMB( MI, MI )
REAL, INTENT(OUT)                        :: APB( MI, MI )
REAL, INTENT(OUT)                        :: ARRAY( MI, * )
REAL, INTENT(IN)                         :: CMU( MXCMU )
REAL, INTENT(IN)                         :: CWT( MXCMU )
REAL, INTENT(IN)                         :: GL( 0:MXCMU )
INTEGER, INTENT(IN OUT)                  :: MI
INTEGER, INTENT(IN)                      :: MAZIM
INTEGER, INTENT(IN OUT)                  :: MXCMU
INTEGER, INTENT(IN OUT)                  :: NN
INTEGER, INTENT(IN)                      :: NSTR
REAL, INTENT(IN)                         :: YLMC( 0:MXCMU, MXCMU )
REAL, INTENT(OUT)                        :: CC( MXCMU, MXCMU )
REAL, INTENT(IN OUT)                     :: EVECC( MXCMU, MXCMU )
REAL, INTENT(OUT)                        :: EVAL( MI )
REAL, INTENT(OUT)                        :: KK( MXCMU )
REAL, INTENT(OUT)                        :: GC( MXCMU, MXCMU )
REAL, INTENT(IN OUT)                     :: WKD( MXCMU )
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

!         Solves eigenvalue/vector problem necessary to construct
!         homogeneous part of discrete ordinate solution; STWJ(8b),
!         STWL(23f)
!         ** NOTE ** Eigenvalue problem is degenerate when single
!                    scattering albedo = 1;  present way of doing it
!                    seems numerically more stable than alternative
!                    methods that we tried


!   I N P U T     V A R I A B L E S:

!       GL     :  Delta-M scaled Legendre coefficients of phase function
!                 (including factors 2l+1 and single-scatter albedo)

!       CMU    :  Computational polar angle cosines

!       CWT    :  Weights for quadrature over polar angle cosine

!       MAZIM  :  Order of azimuthal component

!       NN     :  Half the total number of streams

!       YLMC   :  Normalized associated Legendre polynomial
!                 at the quadrature angles CMU

!       (remainder are DISORT input variables)


!   O U T P U T    V A R I A B L E S:

!       CC     :  C-sub-ij in Eq. SS(5); needed in SS(15&18)

!       EVAL   :  NN eigenvalues of Eq. SS(12), STWL(23f) on return
!                 from ASYMTX but then square roots taken

!       EVECC  :  NN eigenvectors  (G+) - (G-)  on return
!                 from ASYMTX ( column j corresponds to EVAL(j) )
!                 but then  (G+) + (G-)  is calculated from SS(10),
!                 G+  and  G-  are separated, and  G+  is stacked on
!                 top of  G-  to form NSTR eigenvectors of SS(7)

!       GC     :  Permanent storage for all NSTR eigenvectors, but
!                 in an order corresponding to KK

!       KK     :  Permanent storage for all NSTR eigenvalues of SS(7),
!                 but re-ordered with negative values first ( square
!                 roots of EVAL taken and negatives added )


!   I N T E R N A L   V A R I A B L E S:

!       AMB,APB :  Matrices (alpha-beta), (alpha+beta) in reduced
!                    eigenvalue problem
!       ARRAY   :  Complete coefficient matrix of reduced eigenvalue
!                    problem: (alfa+beta)*(alfa-beta)
!       GPPLGM  :  (G+) + (G-) (cf. Eqs. SS(10-11))
!       GPMIGM  :  (G+) - (G-) (cf. Eqs. SS(10-11))
!       WKD     :  Scratch array required by ASYMTX

!   Called by- DISORT, ALBTRN
!   Calls- ASYMTX, ERRMSG
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..


!     ..
!     .. Array Arguments ..




!     ..
!     .. Local Scalars ..

INTEGER :: IER, IQ, JQ, KQ, L
REAL :: ALPHA, BETA, GPMIGM, GPPLGM, SUM
!     ..
!     .. External Subroutines ..

EXTERNAL  ASYMTX, ERRMSG
!     ..
!     .. Intrinsic Functions ..

INTRINSIC ABS, SQRT
!     ..

!                             ** Calculate quantities in Eqs. SS(5-6),
!                             ** STWL(8b,15,23f)
DO  IQ = 1, NN
  
  DO  JQ = 1, NSTR
    
    SUM  = 0.0
    DO  L = MAZIM, NSTR - 1
      SUM  = SUM + GL( L )*YLMC( L, IQ )*YLMC( L, JQ )
    END DO
    
    CC( IQ, JQ ) = 0.5*SUM*CWT( JQ )
    
  END DO
  
  DO  JQ = 1, NN
!                             ** Fill remainder of array using symmetry
!                             ** relations  C(-mui,muj) = C(mui,-muj)
!                             ** and        C(-mui,-muj) = C(mui,muj)
    
    CC( IQ + NN, JQ ) = CC( IQ, JQ + NN )
    CC( IQ + NN, JQ + NN ) = CC( IQ, JQ )
    
!                                       ** Get factors of coeff. matrix
!                                       ** of reduced eigenvalue problem
    
    ALPHA  = CC( IQ, JQ ) / CMU( IQ )
    BETA   = CC( IQ, JQ + NN ) / CMU( IQ )
    AMB( IQ, JQ ) = ALPHA - BETA
    APB( IQ, JQ ) = ALPHA + BETA
    
  END DO
  
  AMB( IQ, IQ ) = AMB( IQ, IQ ) - 1.0 / CMU( IQ )
  APB( IQ, IQ ) = APB( IQ, IQ ) - 1.0 / CMU( IQ )
  
END DO
!                      ** Finish calculation of coefficient matrix of
!                      ** reduced eigenvalue problem:  get matrix
!                      ** product (alfa+beta)*(alfa-beta); SS(12),
!                      ** STWL(23f)
DO  IQ = 1, NN
  
  DO  JQ = 1, NN
    
    SUM  = 0.
    DO  KQ = 1, NN
      SUM  = SUM + APB( IQ, KQ )*AMB( KQ, JQ )
    END DO
    
    ARRAY( IQ, JQ ) = SUM
    
  END DO
  
END DO
!                      ** Find (real) eigenvalues and eigenvectors

CALL ASYMTX( ARRAY, EVECC, EVAL, NN, MI, MXCMU, IER, WKD )

IF( IER > 0 ) THEN
  
  WRITE( KSTDWARN, '(//,A,I4,A)' ) ' ASYMTX--eigenvalue no. ',  &
      IER, '  didnt converge.  Lower-numbered eigenvalues wrong.'
  
! commented this out ----- sergio
  CALL ERRMSG( 'ASYMTX--convergence problems',.True.)
  
END IF


DO  IQ = 1, NN
  EVAL( IQ )    = SQRT( ABS( EVAL( IQ ) ) )
  KK( IQ + NN ) = EVAL( IQ )
!                                      ** Add negative eigenvalue
  KK( NN + 1 - IQ ) = -EVAL( IQ )
END DO

!                          ** Find eigenvectors (G+) + (G-) from SS(10)
!                          ** and store temporarily in APB array
DO  JQ = 1, NN
  
  DO  IQ = 1, NN
    
    SUM  = 0.
    DO  KQ = 1, NN
      SUM  = SUM + AMB( IQ, KQ )*EVECC( KQ, JQ )
    END DO
    
    APB( IQ, JQ ) = SUM / EVAL( JQ )
    
  END DO
  
END DO


DO  JQ = 1, NN
  
  DO  IQ = 1, NN
    
    GPPLGM = APB( IQ, JQ )
    GPMIGM = EVECC( IQ, JQ )
!                                ** Recover eigenvectors G+,G- from
!                                ** their sum and difference; stack them
!                                ** to get eigenvectors of full system
!                                ** SS(7) (JQ = eigenvector number)
    
    EVECC( IQ,      JQ ) = 0.5*( GPPLGM + GPMIGM )
    EVECC( IQ + NN, JQ ) = 0.5*( GPPLGM - GPMIGM )
    
!                                ** Eigenvectors corresponding to
!                                ** negative eigenvalues (corresp. to
!                                ** reversing sign of 'k' in SS(10) )
    GPPLGM = - GPPLGM
    EVECC(IQ,   JQ+NN) = 0.5 * ( GPPLGM + GPMIGM )
    EVECC(IQ+NN,JQ+NN) = 0.5 * ( GPPLGM - GPMIGM )
    GC( IQ+NN,   JQ+NN )   = EVECC( IQ,    JQ )
    GC( NN+1-IQ, JQ+NN )   = EVECC( IQ+NN, JQ )
    GC( IQ+NN,   NN+1-JQ ) = EVECC( IQ,    JQ+NN )
    GC( NN+1-IQ, NN+1-JQ ) = EVECC( IQ+NN, JQ+NN )
    
  END DO
  
END DO


RETURN
END SUBROUTINE SOLEIG

!************************************************************************

SUBROUTINE SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,  &
    FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM,  &
    MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, NNLYRI,  &
    PI, TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )

!        Construct right-hand side vector B for general boundary
!        conditions STWJ(17) and solve system of equations obtained
!        from the boundary conditions and the continuity-of-
!        intensity-at-layer-interface equations.
!        Thermal emission contributes only in azimuthal independence.


!    I N P U T      V A R I A B L E S:

!       BDR      :  Surface bidirectional reflectivity

!       BEM      :  Surface bidirectional emissivity

!       BPLANK   :  Bottom boundary thermal emission

!       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
!                   scaled by Eq. SC(12); in banded form required
!                   by LINPACK solution routines

!       CMU,CWT  :  Abscissae, weights for Gauss quadrature
!                   over angle cosine

!       EXPBEA   :  Transmission of incident beam, EXP(-TAUCPR/UMU0)

!       LYRCUT   :  Logical flag for truncation of computational layers

!       MAZIM    :  Order of azimuthal component

!       NCOL     :  Number of columns in CBAND

!       NN       :  Order of double-Gauss quadrature (NSTR/2)

!       NCUT     :  Total number of computational layers considered

!       TPLANK   :  Top boundary thermal emission

!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)

!       ZZ       :  Beam source vectors in Eq. SS(19), STWL(24b)

!       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16),
!                   Y0 in STWL(26b)

!       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16),
!                   Y1 in STWL(26a)

!       (remainder are DISORT input variables)


!    O U T P U T     V A R I A B L E S:

!       B        :  Right-hand side vector of Eq. SC(5) going into
!                   SGBSL; returns as solution vector of Eq. SC(12),
!                   constants of integration without exponential term

!      LL        :  Permanent storage for B, but re-ordered


!   I N T E R N A L    V A R I A B L E S:

!       IPVT     :  Integer vector of pivot indices
!       IT       :  Pointer for position in  B
!       NCD      :  Number of diagonals below or above main diagonal
!       RCOND    :  Indicator of singularity for CBAND
!       Z        :  Scratch array required by SGBCO

!   Called by- DISORT
!   Calls- ZEROIT, SGBCO, ERRMSG, SGBSL
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..


REAL, INTENT(OUT)                        :: B( NNLYRI )
REAL, INTENT(IN)                         :: BDR( MI, 0:MI )
REAL, INTENT(IN)                         :: BEM( MI )
REAL, INTENT(IN)                         :: BPLANK
REAL, INTENT(IN OUT)                     :: CBAND( MI9M2, NNLYRI )
REAL, INTENT(IN)                         :: CMU( MXCMU )
REAL, INTENT(IN)                         :: CWT( MXCMU )
REAL, INTENT(IN)                         :: EXPBEA( 0:* )
REAL, INTENT(IN)                         :: FBEAM
REAL, INTENT(IN)                         :: FISOT
INTEGER, INTENT(IN OUT)                  :: IPVT( * )
LOGICAL, INTENT(IN)                      :: LAMBER
REAL, INTENT(OUT)                        :: LL( MXCMU, * )
LOGICAL, INTENT(IN)                      :: LYRCUT
INTEGER, INTENT(IN)                      :: MAZIM
INTEGER, INTENT(IN OUT)                  :: MI
INTEGER, INTENT(IN OUT)                  :: MI9M2
INTEGER, INTENT(IN OUT)                  :: MXCMU
INTEGER, INTENT(OUT)                     :: NCOL
INTEGER, INTENT(IN)                      :: NCUT
INTEGER, INTENT(IN OUT)                  :: NN
INTEGER, INTENT(IN)                      :: NSTR
INTEGER, INTENT(IN OUT)                  :: NNLYRI
REAL, INTENT(IN)                         :: PI
REAL, INTENT(IN)                         :: TPLANK
REAL, INTENT(IN OUT)                     :: TAUCPR( 0:* )
REAL, INTENT(IN)                         :: UMU0
REAL, INTENT(IN OUT)                     :: Z( NNLYRI )
REAL, INTENT(IN)                         :: ZZ( MXCMU, * )
REAL, INTENT(IN)                         :: ZPLK0( MXCMU, * )
REAL, INTENT(IN OUT)                     :: ZPLK1( MXCMU, * )



!     ..
!     .. Array Arguments ..



!     ..
!     .. Local Scalars ..

INTEGER :: IPNT, IQ, IT, JQ, LC, NCD
REAL :: RCOND, SUM
!     ..
!     .. External Subroutines ..

EXTERNAL  ERRMSG, SGBCO, SGBSL, ZEROIT
!     ..
!     .. Intrinsic Functions ..

INTRINSIC EXP
!     ..


CALL ZEROIT( B, NNLYRI )
!                              ** Construct B,  STWJ(20a,c) for
!                              ** parallel beam + bottom reflection +
!                              ** thermal emission at top and/or bottom

IF( MAZIM > 0 .AND. FBEAM > 0.0 ) THEN
  
!                                         ** Azimuth-dependent case
!                                         ** (never called if FBEAM = 0)
  IF( LYRCUT .OR. LAMBER ) THEN
    
!               ** No azimuthal-dependent intensity for Lambert surface;
!               ** no intensity component for truncated bottom layer
    
    DO  IQ = 1, NN
!                                                  ** Top boundary
      B( IQ ) = -ZZ( NN + 1 - IQ, 1 )
!                                                  ** Bottom boundary
      
      B( NCOL - NN + IQ ) = -ZZ( IQ + NN, NCUT )*EXPBEA( NCUT )
      
    END DO
    
    
  ELSE
    
    DO  IQ = 1, NN
      
      B( IQ ) = -ZZ( NN + 1 - IQ, 1 )
      
      SUM  = 0.
      DO  JQ = 1, NN
        SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ )*  &
            ZZ( NN + 1 - JQ, NCUT )*EXPBEA( NCUT )
      END DO
      
      B( NCOL - NN + IQ ) = SUM
      IF( FBEAM > 0.0 ) B( NCOL - NN + IQ ) = SUM +  &
          ( BDR( IQ,0 )*UMU0*FBEAM/PI - ZZ( IQ+NN,NCUT ) )* EXPBEA( NCUT )
      
    END DO
    
  END IF
!                             ** Continuity condition for layer
!                             ** interfaces of Eq. STWJ(20b)
  IT  = NN
  
  DO  LC = 1, NCUT - 1
    
    DO  IQ = 1, NSTR
      IT  = IT + 1
      B( IT ) = ( ZZ( IQ,LC+1 ) - ZZ( IQ,LC ) )*EXPBEA( LC )
    END DO
    
  END DO
  
  
ELSE
!                                   ** Azimuth-independent case
  
  IF( FBEAM == 0.0 ) THEN
    
    DO  IQ = 1, NN
!                                      ** Top boundary
      
      B( IQ ) = -ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK
      
    END DO
    
    
    IF( LYRCUT ) THEN
!                               ** No intensity component for truncated
!                               ** bottom layer
      DO  IQ = 1, NN
!                                      ** Bottom boundary
        
        B( NCOL - NN + IQ ) = - ZPLK0( IQ + NN, NCUT ) -  &
            ZPLK1( IQ + NN, NCUT ) * TAUCPR( NCUT )
      END DO
      
      
    ELSE
      
      DO  IQ = 1, NN
        
        SUM  = 0.
        DO  JQ = 1, NN
          SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ )*  &
              ( ZPLK0( NN+1-JQ, NCUT ) + ZPLK1( NN+1-JQ, NCUT ) *TAUCPR( NCUT ) )
        END DO
        
        B( NCOL - NN + IQ ) = 2.*SUM + BEM( IQ )*BPLANK -  &
            ZPLK0( IQ + NN, NCUT ) - ZPLK1( IQ + NN, NCUT ) *  &
            TAUCPR( NCUT )
      END DO
      
    END IF
!                             ** Continuity condition for layer
!                             ** interfaces, STWJ(20b)
    IT  = NN
    DO  LC = 1, NCUT - 1
      
      DO  IQ = 1, NSTR
        IT  = IT + 1
        B( IT ) =   ZPLK0( IQ, LC + 1 ) - ZPLK0( IQ, LC ) +  &
            ( ZPLK1( IQ, LC + 1 ) - ZPLK1( IQ, LC ) )* TAUCPR( LC )
      END DO
      
    END DO
    
    
  ELSE
    
    DO  IQ = 1, NN
      B( IQ ) = -ZZ( NN + 1 - IQ, 1 ) -  &
          ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK
    END DO
    
    IF( LYRCUT ) THEN
      
      DO  IQ = 1, NN
        B( NCOL-NN+IQ ) = - ZZ(IQ+NN, NCUT) * EXPBEA(NCUT)  &
            - ZPLK0(IQ+NN, NCUT) - ZPLK1(IQ+NN, NCUT) * TAUCPR(NCUT)
      END DO
      
      
    ELSE
      
      DO  IQ = 1, NN
        
        SUM  = 0.
        DO  JQ = 1, NN
          SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)  &
              * ( ZZ(NN+1-JQ, NCUT) * EXPBEA(NCUT) + ZPLK0(NN+1-JQ, NCUT)  &
              + ZPLK1(NN+1-JQ, NCUT) * TAUCPR(NCUT))
        END DO
        
        B(NCOL-NN+IQ) = 2.*SUM + ( BDR(IQ,0) * UMU0*FBEAM/PI  &
            - ZZ(IQ+NN, NCUT) ) * EXPBEA(NCUT) + BEM(IQ) * BPLANK  &
            - ZPLK0(IQ+NN, NCUT) - ZPLK1(IQ+NN, NCUT) * TAUCPR(NCUT)
      END DO
      
    END IF
    
    
    IT  = NN
    
    DO  LC = 1, NCUT - 1
      
      DO  IQ = 1, NSTR
        
        IT  = IT + 1
        B(IT) = ( ZZ(IQ,LC+1) - ZZ(IQ,LC) ) * EXPBEA(LC)  &
            + ZPLK0(IQ,LC+1) - ZPLK0(IQ,LC) +  &
            ( ZPLK1(IQ,LC+1) - ZPLK1(IQ,LC) ) * TAUCPR(LC)
      END DO
      
    END DO
    
  END IF
  
END IF
!                     ** Find L-U (lower/upper triangular) decomposition
!                     ** of band matrix CBAND and test if it is nearly
!                     ** singular (note: CBAND is destroyed)
!                     ** (CBAND is in LINPACK packed format)
RCOND  = 0.0
NCD    = 3*NN - 1

CALL SGBCO( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z )

IF( 1.0 + RCOND == 1.0 )  &
    CALL ERRMSG('SOLVE0--SGBCO says matrix near singular',.FALSE.)

!                   ** Solve linear system with coeff matrix CBAND
!                   ** and R.H. side(s) B after CBAND has been L-U
!                   ** decomposed.  Solution is returned in B.

CALL SGBSL( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, B, 0 )

!                   ** Zero CBAND (it may contain 'foreign'
!                   ** elements upon returning from LINPACK);
!                   ** necessary to prevent errors

CALL ZEROIT( CBAND, MI9M2*NNLYRI )

DO  LC = 1, NCUT
  
  IPNT  = LC*NSTR - NN
  
  DO  IQ = 1, NN
    LL( NN + 1 - IQ, LC ) = B( IPNT + 1 - IQ )
    LL( IQ + NN,     LC ) = B( IQ + IPNT )
  END DO
  
END DO


RETURN
END SUBROUTINE SOLVE0

!************************************************************************

SUBROUTINE SURFAC( ALBEDO  DELM0, CMU, FBEAM, LAMBER, MI, MAZIM,  &
      MXUMU, NN, NUMU, ONLYFL, PI, UMU, UMU0,  &
      USRANG, WVNMLO, WVNMHI, BDR, EMU, BEM, RMU )
  
!       Computes user's surface bidirectional properties, STWL(41)
  
!   I N P U T     V A R I A B L E S:
  
!       CMU    :  Computational polar angle cosines (Gaussian)
  
!       DELM0  :  Kronecker delta, delta-sub-m0
  
!       MAZIM  :  Order of azimuthal component
  
!       NN     :  Order of Double-Gauss quadrature (NSTR/2)
  
!       (Remainder are 'DISORT' input variables)
  
!    O U T P U T     V A R I A B L E S:
  
!       BDR :  Fourier expansion coefficient of surface bidirectional
!                 reflectivity (computational angles)
  
!       RMU :  Surface bidirectional reflectivity (user angles)
  
!       BEM :  Surface directional emissivity (computational angles)
  
!       EMU :  Surface directional emissivity (user angles)
  
!    I N T E R N A L     V A R I A B L E S:
  
!       DREF   :  Directional reflectivity
  
!       NMUG   :  Number of angle cosine quadrature points on (-1,1)
!                 for integrating bidirectional reflectivity to get
!                 directional emissivity (it is necessary to use a
!                 quadrature set distinct from the computational angles,
!                 because the computational angles may not be dense
!                 enough -- i.e. 'NSTR' may be too small-- to give an
!                 accurate approximation for the integration).
  
!       GMU    :  The 'NMUG' angle cosine quadrature points on (0,1)
  
!       GWT    :  The 'NMUG' angle cosine quadrature weights on (0,1)
  
!   Called by- DISORT
!   Calls- QGAUSN, BDREF, ZEROIT
!+---------------------------------------------------------------------+
  
!     .. Parameters ..
  
  
  REAL, INTENT(IN OUT)                     :: ALBEDO  DE
    REAL, INTENT(IN OUT)                     :: CMU( * )
    REAL, INTENT(IN)                         :: FBEAM
    LOGICAL, INTENT(IN)                      :: LAMBER
    INTEGER, INTENT(IN OUT)                  :: MI
    INTEGER, INTENT(IN)                      :: MAZIM
    INTEGER, INTENT(IN OUT)                  :: MXUMU
    INTEGER, INTENT(IN)                      :: NN
    INTEGER, INTENT(IN)                      :: NUMU
    LOGICAL, INTENT(IN OUT)                  :: ONLYFL
    REAL, INTENT(IN OUT)                     :: PI
    REAL, INTENT(IN)                         :: UMU( * )
    REAL, INTENT(IN OUT)                     :: UMU0
    LOGICAL, INTENT(IN)                      :: USRANG
    REAL, INTENT(IN OUT)                     :: WVNMLO
    REAL, INTENT(IN OUT)                     :: WVNMHI
    REAL, INTENT(OUT)                        :: BDR( MI, 0:MI )
    REAL, INTENT(OUT)                        :: EMU( MXUMU )
    REAL, INTENT(OUT)                        :: BEM( MI )
    REAL, INTENT(OUT)                        :: RMU( MXUMU, 0:MI )
    
    INTEGER, PARAMETER :: NMUG = 50
!     ..
!     .. Scalar Arguments ..
    
    
    
    REAL :: ALBEDO  DELM0
!     ..
!     .. Array Arguments ..
      
!     ..
!     .. Local Scalars ..
      
      LOGICAL :: PASS1
      INTEGER :: IQ, IU, JG, JQ, K
      REAL :: DREF, SUM
!     ..
!     .. Local Arrays ..
      
      REAL :: GMU( NMUG ), GWT( NMUG )
!     ..
!     .. External Functions ..
      
      REAL :: BDREF
      EXTERNAL  BDREF
!     ..
!     .. External Subroutines ..
      
      EXTERNAL  QGAUSN, ZEROIT
!     ..
!     .. Intrinsic Functions ..
      
      INTRINSIC COS
!     ..
      SAVE      PASS1, GMU, GWT
      DATA      PASS1 / .True. /
      
      
      IF( PASS1 ) THEN
        
        PASS1  = .FALSE.
        
        CALL QGAUSN( NMUG/2, GMU, GWT )
        
        DO  K = 1, NMUG / 2
          GMU( K + NMUG/2 ) = -GMU( K )
          GWT( K + NMUG/2 ) = GWT( K )
        END DO
        
      END IF
      
      
      CALL ZEROIT( BDR, MI*( MI+1 ) )
      CALL ZEROIT( BEM, MI )
      
!                             ** Compute Fourier expansion coefficient
!                             ** of surface bidirectional reflectance
!                             ** at computational angles Eq. STWL (41)
      
      IF( LAMBER .AND. MAZIM == 0 ) THEN
        
        DO  IQ = 1, NN
          
          BEM( IQ ) = 1.0 - ALBEDO
            
            DO  JQ = 0, NN
              BDR( IQ, JQ ) = ALBEDO
              END DO
              
            END DO
            
          ELSE IF( .NOT.LAMBER ) THEN
            
            DO  IQ = 1, NN
              
              DO  JQ = 1, NN
                
                SUM  = 0.0
                DO  K = 1, NMUG
                  SUM  = SUM + GWT( K ) *  &
                      BDREF( WVNMLO, WVNMHI, CMU(IQ), CMU(JQ),  &
                      PI*GMU(K) ) * COS( MAZIM*PI*GMU( K ) )
                END DO
                
                BDR( IQ, JQ ) = 0.5 * ( 2. - DELM0 ) * SUM
                
              END DO
              
              
              IF( FBEAM > 0.0 ) THEN
                
                SUM  = 0.0
                DO  K = 1, NMUG
                  SUM  = SUM + GWT( K ) *  &
                      BDREF( WVNMLO, WVNMHI, CMU(IQ), UMU0,  &
                      PI*GMU(K) ) * COS( MAZIM*PI*GMU( K ) )
                END DO
                
                BDR( IQ, 0 ) = 0.5 * ( 2. - DELM0 ) * SUM
                
              END IF
              
            END DO
            
            
            IF( MAZIM == 0 ) THEN
              
!                             ** Integrate bidirectional reflectivity
!                             ** at reflection polar angle cosines -CMU-
!                             ** and incident angle cosines -GMU- to get
!                             ** directional emissivity at computational
!                             ** angle cosines -CMU-.
              DO  IQ = 1, NN
                
                DREF  = 0.0
                
                DO  JG = 1, NMUG
                  
                  SUM  = 0.0
                  DO  K = 1, NMUG / 2
                    SUM  = SUM + GWT( K ) * GMU( K ) *  &
                        BDREF( WVNMLO, WVNMHI, CMU(IQ), GMU(K), PI*GMU(JG) )
                  END DO
                  
                  DREF  = DREF + GWT( JG )*SUM
                  
                END DO
                
                BEM( IQ ) = 1.0 - DREF
                
              END DO
              
            END IF
            
          END IF
!                             ** Compute Fourier expansion coefficient
!                             ** of surface bidirectional reflectance
!                             ** at user angles Eq. STWL (41)
          
          IF( .NOT.ONLYFL .AND. USRANG ) THEN
            
            CALL ZEROIT( EMU, MXUMU )
            CALL ZEROIT( RMU, MXUMU*( MI+1 ) )
            
            DO  IU = 1, NUMU
              
              IF( UMU( IU ) > 0.0 ) THEN
                
                IF( LAMBER .AND. MAZIM == 0 ) THEN
                  
                  DO  IQ = 0, NN
                    RMU( IU, IQ ) = ALBEDO
                    END DO
                    
                    EMU( IU ) = 1.0 - ALBEDO
                      
                    ELSE IF( .NOT.LAMBER ) THEN
                      
                      DO  IQ = 1, NN
                        
                        SUM  = 0.0
                        DO  K = 1, NMUG
                          SUM  = SUM + GWT( K ) *  &
                              BDREF( WVNMLO, WVNMHI, UMU(IU), CMU(IQ),  &
                              PI*GMU(K) ) * COS( MAZIM*PI*GMU( K ) )
                        END DO
                        
                        RMU( IU, IQ ) = 0.5 * ( 2. - DELM0 ) * SUM
                        
                      END DO
                      
                      IF( FBEAM > 0.0 ) THEN
                        
                        SUM  = 0.0
                        DO  K = 1, NMUG
                          SUM  = SUM + GWT( K ) *  &
                              BDREF( WVNMLO, WVNMHI, UMU(IU), UMU0,  &
                              PI*GMU(K) ) * COS( MAZIM*PI*GMU( K ) )
                        END DO
                        
                        RMU( IU, 0 ) = 0.5 * ( 2. - DELM0 ) * SUM
                        
                      END IF
                      
                      
                      IF( MAZIM == 0 ) THEN
                        
!                               ** Integrate bidirectional reflectivity
!                               ** at reflection angle cosines -UMU- and
!                               ** incident angle cosines -GMU- to get
!                               ** directional emissivity at
!                               ** user angle cosines -UMU-.
                        DREF  = 0.0
                        
                        DO  JG = 1, NMUG
                          
                          SUM  = 0.0
                          DO  K = 1, NMUG / 2
                            SUM  = SUM + GWT( K )*GMU( K )*  &
                                BDREF( WVNMLO, WVNMHI, UMU(IU),  &
                                GMU(K), PI*GMU(JG) )
                          END DO
                          
                          DREF  = DREF + GWT( JG ) * SUM
                          
                        END DO
                        
                        EMU( IU ) = 1.0 - DREF
                        
                      END IF
                      
                    END IF
                    
                  END IF
                  
                END DO
                
              END IF
              
              
              RETURN
            END SUBROUTINE SURFAC
            
!************************************************************************
            
            SUBROUTINE TERPEV( CWT, EVECC, GL, GU, MAZIM, MXCMU, MXUMU, NN,  &
                NSTR, NUMU, WK, YLMC, YLMU )
            
!         Interpolate eigenvectors to user angles; Eq SD(8)
            
!   Called by- DISORT, ALBTRN
! --------------------------------------------------------------------+
            
!     .. Scalar Arguments ..
            
            
            REAL, INTENT(IN)                         :: CWT( MXCMU )
            REAL, INTENT(IN)                         :: EVECC( MXCMU, MXCMU )
            REAL, INTENT(IN)                         :: GL( 0:MXCMU )
            REAL, INTENT(OUT)                        :: GU( MXUMU, MXCMU )
            INTEGER, INTENT(IN)                      :: MAZIM
            INTEGER, INTENT(IN OUT)                  :: MXCMU
            INTEGER, INTENT(IN OUT)                  :: MXUMU
            INTEGER, INTENT(IN OUT)                  :: NN
            INTEGER, INTENT(IN OUT)                  :: NSTR
            INTEGER, INTENT(IN)                      :: NUMU
            REAL, INTENT(OUT)                        :: WK( MXCMU )
            REAL, INTENT(IN)                         :: YLMC( 0:MXCMU, MXCMU )
            REAL, INTENT(IN)                         :: YLMU( 0:MXCMU, MXUMU )
            
!     ..
!     .. Array Arguments ..
            
            
!     ..
!     .. Local Scalars ..
            
            INTEGER :: IQ, IU, JQ, L
            REAL :: SUM
!     ..
            
            
            DO  IQ = 1, NSTR
              
              DO  L = MAZIM, NSTR - 1
!                                   ** Inner sum in SD(8) times all
!                                   ** factors in outer sum but PLM(mu)
                SUM  = 0.0
                DO  JQ = 1, NSTR
                  SUM  = SUM + CWT( JQ )*YLMC( L, JQ )*EVECC( JQ, IQ )
                END DO
                
                WK( L + 1 ) = 0.5*GL( L )*SUM
                
              END DO
!                                    ** Finish outer sum in SD(8)
!                                    ** and store eigenvectors
              DO  IU = 1, NUMU
                
                SUM  = 0.
                DO  L = MAZIM, NSTR - 1
                  SUM  = SUM + WK( L + 1 )*YLMU( L, IU )
                END DO
                
                IF( IQ <= NN ) GU( IU, IQ + NN )       = SUM
                IF( IQ > NN ) GU( IU, NSTR + 1 - IQ ) = SUM
                
              END DO
              
            END DO
            
            
            RETURN
          END SUBROUTINE TERPEV
          
!************************************************************************
          
          SUBROUTINE TERPSO( CWT, DELM0, FBEAM, GL, MAZIM, MXCMU, PLANK,  &
              NUMU, NSTR, OPRIM, PI, YLM0, YLMC, YLMU, PSI0,  &
              PSI1, XR0, XR1, Z0, Z1, ZJ, ZBEAM, Z0U, Z1U )
          
!         Interpolates source functions to user angles, Eq. STWL(30)
          
          
!    I N P U T      V A R I A B L E S:
          
!       CWT    :  Weights for Gauss quadrature over angle cosine
          
!       DELM0  :  Kronecker delta, delta-sub-m0
          
!       GL     :  Delta-M scaled Legendre coefficients of phase function
!                 (including factors 2L+1 and single-scatter albedo)
          
!       MAZIM  :  Order of azimuthal component
          
!       OPRIM  :  Single scattering albedo
          
!       XR0    :  Expansion of thermal source function, Eq. STWL(24d)
          
!       XR1    :  Expansion of thermal source function Eq. STWL(24d)
          
!       YLM0   :  Normalized associated Legendre polynomial
!                 at the beam angle
          
!       YLMC   :  Normalized associated Legendre polynomial
!                 at the quadrature angles
          
!       YLMU   :  Normalized associated Legendre polynomial
!                 at the user angles
          
!       Z0     :  Solution vectors Z-sub-zero of Eq. SS(16), STWL(26a)
          
!       Z1     :  Solution vectors Z-sub-one  of Eq. SS(16), STWL(26b)
          
!       ZJ     :  Solution vector Z-sub-zero after solving Eq. SS(19),
!                 STWL(24b)
          
!       (remainder are DISORT input variables)
          
          
!    O U T P U T     V A R I A B L E S:
          
!       ZBEAM  :  Incident-beam source function at user angles
          
!       Z0U,Z1U:  Components of a linear-in-optical-depth-dependent
!                 source (approximating the Planck emission source)
          
          
!   I N T E R N A L    V A R I A B L E S:
          
!       PSI0  :  Sum just after square bracket in  Eq. SD(9)
!       PSI1  :  Sum in Eq. STWL(31d)
          
!   Called by- DISORT
! +-------------------------------------------------------------------+
          
!     .. Scalar Arguments ..
          
          
          REAL, INTENT(IN)                         :: CWT( MXCMU )
          REAL, INTENT(IN)                         :: DELM0
          REAL, INTENT(IN)                         :: FBEAM
          REAL, INTENT(IN)                         :: GL( 0:MXCMU )
          INTEGER, INTENT(IN)                      :: MAZIM
          INTEGER, INTENT(IN OUT)                  :: MXCMU
          LOGICAL, INTENT(IN)                      :: PLANK
          INTEGER, INTENT(IN)                      :: NUMU
          INTEGER, INTENT(IN)                      :: NSTR
          REAL, INTENT(IN)                         :: OPRIM
          REAL, INTENT(IN)                         :: PI
          REAL, INTENT(IN OUT)                     :: YLM0( 0:MXCMU )
          REAL, INTENT(IN)                         :: YLMC( 0:MXCMU, MXCMU )
          REAL, INTENT(IN)                         :: YLMU( 0:MXCMU, * )
          REAL, INTENT(OUT)                        :: PSI0( MXCMU )
          REAL, INTENT(OUT)                        :: PSI1( MXCMU )
          REAL, INTENT(IN)                         :: XR0
          REAL, INTENT(IN)                         :: XR1
          REAL, INTENT(IN)                         :: Z0( MXCMU )
          REAL, INTENT(IN)                         :: Z1( MXCMU )
          REAL, INTENT(IN)                         :: ZJ( MXCMU )
          REAL, INTENT(OUT)                        :: ZBEAM( * )
          REAL, INTENT(OUT)                        :: Z0U( * )
          REAL, INTENT(OUT)                        :: Z1U( * )
          
          
          
!     ..
!     .. Array Arguments ..
          
          
!     ..
!     .. Local Scalars ..
          
          INTEGER :: IQ, IU, JQ
          REAL :: FACT, PSUM, PSUM0, PSUM1, SUM, SUM0, SUM1
!     ..
          
          
          IF( FBEAM > 0.0 ) THEN
!                                  ** Beam source terms; Eq. SD(9)
            
            DO  IQ = MAZIM, NSTR - 1
              
              PSUM   = 0.
              DO  JQ = 1, NSTR
                PSUM  = PSUM + CWT( JQ )*YLMC( IQ, JQ )*ZJ( JQ )
              END DO
              
              PSI0( IQ + 1 ) = 0.5*GL( IQ )*PSUM
              
            END DO
            
            FACT   = ( 2. - DELM0 )*FBEAM / ( 4.0*PI )
            
            DO  IU = 1, NUMU
              
              SUM    = 0.
              DO  IQ = MAZIM, NSTR - 1
                SUM  = SUM + YLMU( IQ, IU )*  &
                    ( PSI0( IQ+1 ) + FACT*GL( IQ )*YLM0( IQ ) )
              END DO
              
              ZBEAM( IU ) = SUM
              
            END DO
            
          END IF
          
          
          IF( PLANK .AND. MAZIM == 0 ) THEN
            
!                          ** Thermal source terms, STWJ(27c), STWL(31c)
            
            DO  IQ = MAZIM, NSTR - 1
              
              PSUM0  = 0.0
              PSUM1  = 0.0
              DO  JQ = 1, NSTR
                PSUM0  = PSUM0 + CWT( JQ )*YLMC( IQ, JQ )*Z0( JQ )
                PSUM1  = PSUM1 + CWT( JQ )*YLMC( IQ, JQ )*Z1( JQ )
              END DO
              
              PSI0( IQ + 1 ) = 0.5*GL( IQ ) * PSUM0
              PSI1( IQ + 1 ) = 0.5*GL( IQ ) * PSUM1
              
            END DO
            
            DO  IU = 1, NUMU
              
              SUM0   = 0.0
              SUM1   = 0.0
              DO  IQ = MAZIM, NSTR - 1
                SUM0  = SUM0 + YLMU( IQ, IU ) * PSI0( IQ + 1 )
                SUM1  = SUM1 + YLMU( IQ, IU ) * PSI1( IQ + 1 )
              END DO
              
              Z0U( IU ) = SUM0 + ( 1. - OPRIM ) * XR0
              Z1U( IU ) = SUM1 + ( 1. - OPRIM ) * XR1
              
            END DO
            
          END IF
          
          
          RETURN
        END SUBROUTINE TERPSO
        
!************************************************************************
        
        SUBROUTINE UPBEAM( ARRAY, CC, CMU, DELM0, FBEAM, GL, IPVT, MAZIM,  &
            MXCMU, NN, NSTR, PI, UMU0, WK, YLM0, YLMC, ZJ, ZZ )
        
!         Finds the incident-beam particular solution of SS(18),
!         STWL(24a)
        
!   I N P U T    V A R I A B L E S:
        
!       CC     :  C-sub-ij in Eq. SS(5)
        
!       CMU    :  Abscissae for Gauss quadrature over angle cosine
        
!       DELM0  :  Kronecker delta, delta-sub-m0
        
!       GL     :  Delta-M scaled Legendre coefficients of phase function
!                 (including factors 2L+1 and single-scatter albedo)
        
!       MAZIM  :  Order of azimuthal component
        
!       YLM0   :  Normalized associated Legendre polynomial
!                 at the beam angle
        
!       YLMC   :  Normalized associated Legendre polynomial
!                 at the quadrature angles
        
!       (remainder are DISORT input variables)
        
        
!   O U T P U T    V A R I A B L E S:
        
!       ZJ     :  Right-hand side vector X-sub-zero in SS(19),STWL(24b);
!                 also the solution vector Z-sub-zero after solving
!                 that system
        
!       ZZ     :  Permanent storage for ZJ, but re-ordered
        
        
!   I N T E R N A L    V A R I A B L E S:
        
!       ARRAY  :  Coefficient matrix in left-hand side of Eq. SS(19),
!                   STWL(24b)
!       IPVT   :  Integer vector of pivot indices required by LINPACK
!       WK     :  Scratch array required by LINPACK
        
!   Called by- DISORT
!   Calls- SGECO, ERRMSG, SGESL
! +-------------------------------------------------------------------+
        
!     .. Scalar Arguments ..
        
        
        REAL, INTENT(OUT)                        :: ARRAY( MXCMU, MXCMU )
        REAL, INTENT(IN)                         :: CC( MXCMU, MXCMU )
        REAL, INTENT(IN)                         :: CMU( MXCMU )
        REAL, INTENT(IN)                         :: DELM0
        REAL, INTENT(IN)                         :: FBEAM
        REAL, INTENT(IN)                         :: GL( 0:MXCMU )
        INTEGER, INTENT(IN OUT)                  :: IPVT( * )
        INTEGER, INTENT(IN)                      :: MAZIM
        INTEGER, INTENT(IN OUT)                  :: MXCMU
        INTEGER, INTENT(IN OUT)                  :: NN
        INTEGER, INTENT(IN)                      :: NSTR
        REAL, INTENT(IN)                         :: PI
        REAL, INTENT(IN)                         :: UMU0
        REAL, INTENT(IN OUT)                     :: WK( MXCMU )
        REAL, INTENT(IN)                         :: YLM0( 0:MXCMU )
        REAL, INTENT(IN)                         :: YLMC( 0:MXCMU, * )
        REAL, INTENT(OUT)                        :: ZJ( MXCMU )
        REAL, INTENT(OUT)                        :: ZZ( MXCMU )
        
        
!     ..
!     .. Array Arguments ..
        
        
        
!     ..
!     .. Local Scalars ..
        
        INTEGER :: IQ, JOB, JQ, K
        REAL :: RCOND, SUM
!     ..
!     .. External Subroutines ..
        
        EXTERNAL  ERRMSG, SGECO, SGESL
!     ..
        
        
        DO  IQ = 1, NSTR
          
          DO  JQ = 1, NSTR
            ARRAY( IQ, JQ ) = -CC( IQ, JQ )
          END DO
          
          ARRAY( IQ, IQ ) = 1.+ CMU( IQ ) / UMU0 + ARRAY( IQ, IQ )
          
          SUM  = 0.
          DO  K = MAZIM, NSTR - 1
            SUM  = SUM + GL( K )*YLMC( K, IQ )*YLM0( K )
          END DO
          
          ZJ( IQ ) = ( 2.- DELM0 )*FBEAM*SUM / ( 4.*PI )
        END DO
        
!                  ** Find L-U (lower/upper triangular) decomposition
!                  ** of ARRAY and see if it is nearly singular
!                  ** (NOTE:  ARRAY is altered)
        RCOND  = 0.0
        
        CALL SGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )
        
        IF( 1.0 + RCOND == 1.0 )  &
            CALL ERRMSG('UPBEAM--SGECO says matrix near singular',.FALSE.)
        
!                ** Solve linear system with coeff matrix ARRAY
!                ** (assumed already L-U decomposed) and R.H. side(s)
!                ** ZJ;  return solution(s) in ZJ
        JOB  = 0
        
        CALL SGESL( ARRAY, MXCMU, NSTR, IPVT, ZJ, JOB )
        
        
        DO  IQ = 1, NN
          ZZ( IQ + NN )     = ZJ( IQ )
          ZZ( NN + 1 - IQ ) = ZJ( IQ + NN )
        END DO
        
        
        RETURN
      END SUBROUTINE UPBEAM
      
!************************************************************************
      
      SUBROUTINE UPISOT( ARRAY, CC, CMU, IPVT, MXCMU, NN, NSTR, OPRIM,  &
          WK, XR0, XR1, Z0, Z1, ZPLK0, ZPLK1 )
      
!       Finds the particular solution of thermal radiation of STWL(25)
      
      
      
!    I N P U T     V A R I A B L E S:
      
!       CC     :  C-sub-ij in Eq. SS(5), STWL(8b)
      
!       CMU    :  Abscissae for Gauss quadrature over angle cosine
      
!       OPRIM  :  Delta-M scaled single scattering albedo
      
!       XR0    :  Expansion coefficient b-sub-zero of thermal source
!                   function, Eq. STWL(24c)
      
!       XR1    :  Expansion coefficient b-sub-one of thermal source
!                   function Eq. STWL(24c)
      
!       (remainder are DISORT input variables)
      
      
!    O U T P U T    V A R I A B L E S:
      
!       Z0     :  Solution vectors Z-sub-zero of Eq. SS(16), STWL(26a)
      
!       Z1     :  Solution vectors Z-sub-one  of Eq. SS(16), STWL(26b)
      
!       ZPLK0, :  Permanent storage for Z0,Z1, but re-ordered
!        ZPLK1
      
      
!   I N T E R N A L    V A R I A B L E S:
      
!       ARRAY  :  Coefficient matrix in left-hand side of EQ. SS(16)
!       IPVT   :  Integer vector of pivot indices required by LINPACK
!       WK     :  Scratch array required by LINPACK
      
!   Called by- DISORT
!   Calls- SGECO, ERRMSG, SGESL
! +-------------------------------------------------------------------+
      
!     .. Scalar Arguments ..
      
      
      REAL, INTENT(OUT)                        :: ARRAY( MXCMU, MXCMU )
      REAL, INTENT(IN)                         :: CC( MXCMU, MXCMU )
      REAL, INTENT(IN)                         :: CMU( MXCMU )
      INTEGER, INTENT(IN OUT)                  :: IPVT( * )
      INTEGER, INTENT(IN OUT)                  :: MXCMU
      INTEGER, INTENT(IN OUT)                  :: NN
      INTEGER, INTENT(IN)                      :: NSTR
      REAL, INTENT(IN)                         :: OPRIM
      REAL, INTENT(IN OUT)                     :: WK( MXCMU )
      REAL, INTENT(IN)                         :: XR0
      REAL, INTENT(IN)                         :: XR1
      REAL, INTENT(OUT)                        :: Z0( MXCMU )
      REAL, INTENT(OUT)                        :: Z1( MXCMU )
      REAL, INTENT(OUT)                        :: ZPLK0( MXCMU )
      REAL, INTENT(OUT)                        :: ZPLK1( MXCMU )
      
      
!     ..
!     .. Array Arguments ..
      
      
      
!     ..
!     .. Local Scalars ..
      
      INTEGER :: IQ, JQ
      REAL :: RCOND
!     ..
!     .. External Subroutines ..
      
      EXTERNAL  ERRMSG, SGECO, SGESL
!     ..
      
      
      DO  IQ = 1, NSTR
        
        DO  JQ = 1, NSTR
          ARRAY( IQ, JQ ) = -CC( IQ, JQ )
        END DO
        
        ARRAY( IQ, IQ ) = 1.0 + ARRAY( IQ, IQ )
        
        Z1( IQ ) = ( 1. - OPRIM ) * XR1
        
      END DO
!                       ** Solve linear equations: same as in UPBEAM,
!                       ** except ZJ replaced by Z1 and Z0
      RCOND  = 0.0
      
      CALL SGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )
      
      IF( 1.0 + RCOND == 1.0 )  &
          CALL ERRMSG('UPISOT--SGECO says matrix near singular',.False.)
      
      CALL SGESL( ARRAY, MXCMU, NSTR, IPVT, Z1, 0 )
      
      DO  IQ = 1, NSTR
        Z0( IQ ) = ( 1. - OPRIM ) * XR0 + CMU( IQ ) * Z1( IQ )
      END DO
      
      CALL SGESL( ARRAY, MXCMU, NSTR, IPVT, Z0, 0 )
      
      DO  IQ = 1, NN
        ZPLK0( IQ + NN ) = Z0( IQ )
        ZPLK1( IQ + NN ) = Z1( IQ )
        ZPLK0( NN + 1 - IQ ) = Z0( IQ + NN )
        ZPLK1( NN + 1 - IQ ) = Z1( IQ + NN )
      END DO
      
      
      RETURN
    END SUBROUTINE UPISOT
    
!************************************************************************
    
    SUBROUTINE USRINT( BPLANK, CMU, CWT, DELM0, DTAUCP, EMU, EXPBEA,  &
        FBEAM, FISOT, GC, GU, KK, LAMBER, LAYRU, LL,  &
        LYRCUT, MAZIM, MXCMU, MXULV, MXUMU, NCUT, NLYR,  &
        NN, NSTR, PLANK, NUMU, NTAU, PI, RMU, TAUCPR,  &
        TPLANK, UMU, UMU0, UTAUPR, WK, ZBEAM, Z0U, Z1U, ZZ, ZPLK0, ZPLK1, UUM )
    
!       Computes intensity components at user output angles
!       for azimuthal expansion terms in Eq. SD(2), STWL(6)
    
    
!   I N P U T    V A R I A B L E S:
    
!       BPLANK :  Integrated Planck function for emission from
!                 bottom boundary
    
!       CMU    :  Abscissae for Gauss quadrature over angle cosine
    
!       CWT    :  Weights for Gauss quadrature over angle cosine
    
!       DELM0  :  Kronecker delta, delta-sub-M0
    
!       EMU    :  Surface directional emissivity (user angles)
    
!       EXPBEA :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
    
!       GC     :  Eigenvectors at polar quadrature angles, SC(1)
    
!       GU     :  Eigenvectors interpolated to user polar angles
!                    (i.e., G in Eq. SC(1) )
    
!       KK     :  Eigenvalues of coeff. matrix in Eq. SS(7), STWL(23b)
    
!       LAYRU  :  Layer number of user level UTAU
    
!       LL     :  Constants of integration in Eq. SC(1), obtained
!                 by solving scaled version of Eq. SC(5);
!                 exponential term of Eq. SC(12) not included
    
!       LYRCUT :  Logical flag for truncation of computational layer
    
!       MAZIM  :  Order of azimuthal component
    
!       NCUT   :  Total number of computational layers considered
    
!       NN     :  Order of double-Gauss quadrature (NSTR/2)
    
!       RMU    :  Surface bidirectional reflectivity (user angles)
    
!       TAUCPR :  Cumulative optical depth (delta-M-Scaled)
    
!       TPLANK :  Integrated Planck function for emission from
!                 top boundary
    
!       UTAUPR :  Optical depths of user output levels in delta-M
!                 coordinates;  equal to UTAU if no delta-M
    
!       Z0U    :  Z-sub-zero in Eq. SS(16) interpolated to user
!                 angles from an equation derived from SS(16),
!                 Y-sub-zero on STWL(26b)
    
!       Z1U    :  Z-sub-one in Eq. SS(16) interpolated to user
!                 angles from an equation derived from SS(16),
!                 Y-sub-one in STWL(26a)
    
!       ZZ     :  Beam source vectors in Eq. SS(19), STWL(24b)
    
!       ZPLK0  :  Thermal source vectors Z0, by solving Eq. SS(16),
!                 Y-sub-zero in STWL(26)
    
!       ZPLK1  :  Thermal source vectors Z1, by solving Eq. SS(16),
!                 Y-sub-one in STWL(26)
    
!       ZBEAM  :  Incident-beam source vectors
    
!       (Remainder are DISORT input variables)
    
    
!    O U T P U T    V A R I A B L E S:
    
!       UUM    :  Azimuthal components of the intensity in EQ. STWJ(5),
!                 STWL(6)
    
    
!    I N T E R N A L    V A R I A B L E S:
    
!       BNDDIR :  Direct intensity down at the bottom boundary
!       BNDDFU :  Diffuse intensity down at the bottom boundary
!       BNDINT :  Intensity attenuated at both boundaries, STWJ(25-6)
!       DTAU   :  Optical depth of a computational layer
!       LYREND :  End layer of integration
!       LYRSTR :  Start layer of integration
!       PALINT :  Intensity component from parallel beam
!       PLKINT :  Intensity component from planck source
!       WK     :  Scratch vector for saving EXP evaluations
    
!       All the exponential factors ( EXP1, EXPN,... etc.)
!       come from the substitution of constants of integration in
!       Eq. SC(12) into Eqs. S1(8-9).  They all have negative
!       arguments so there should never be overflow problems.
    
!   Called by- DISORT
! +-------------------------------------------------------------------+
    
!     .. Scalar Arguments ..
    
    
    REAL, INTENT(IN)                         :: BPLANK
    REAL, INTENT(IN OUT)                     :: CMU( MXCMU )
    REAL, INTENT(IN OUT)                     :: CWT( MXCMU )
    REAL, INTENT(IN)                         :: DELM0
    REAL, INTENT(IN)                         :: DTAUCP( * )
    REAL, INTENT(IN)                         :: EMU( MXUMU )
    REAL, INTENT(IN)                         :: EXPBEA( 0:* )
    REAL, INTENT(IN)                         :: FBEAM
    REAL, INTENT(IN OUT)                     :: FISOT
    REAL, INTENT(IN)                         :: GC( MXCMU, MXCMU, * )
    REAL, INTENT(OUT)                        :: GU( MXUMU, MXCMU, * )
    REAL, INTENT(IN)                         :: KK( MXCMU, * )
    LOGICAL, INTENT(IN OUT)                  :: LAMBER
    INTEGER, INTENT(IN)                      :: LAYRU( * )
    REAL, INTENT(IN)                         :: LL( MXCMU, * )
    LOGICAL, INTENT(IN)                      :: LYRCUT
    INTEGER, INTENT(IN)                      :: MAZIM
    INTEGER, INTENT(IN OUT)                  :: MXCMU
    INTEGER, INTENT(IN OUT)                  :: MXULV
    INTEGER, INTENT(IN OUT)                  :: MXUMU
    INTEGER, INTENT(IN)                      :: NCUT
    INTEGER, INTENT(IN)                      :: NLYR
    INTEGER, INTENT(IN)                      :: NN
    INTEGER, INTENT(IN)                      :: NSTR
    LOGICAL, INTENT(IN)                      :: PLANK
    INTEGER, INTENT(IN)                      :: NUMU
    INTEGER, INTENT(IN)                      :: NTAU
    REAL, INTENT(IN)                         :: PI
    REAL, INTENT(IN)                         :: RMU( MXUMU, 0:* )
    REAL, INTENT(IN)                         :: TAUCPR( 0:* )
    REAL, INTENT(IN OUT)                     :: TPLANK
    REAL, INTENT(IN)                         :: UMU( * )
    REAL, INTENT(IN)                         :: UMU0
    REAL, INTENT(IN)                         :: UTAUPR( MXULV )
    REAL, INTENT(OUT)                        :: WK( MXCMU )
    REAL, INTENT(IN)                         :: ZBEAM( MXUMU, * )
    REAL, INTENT(IN)                         :: Z0U( MXUMU, * )
    REAL, INTENT(IN)                         :: Z1U( MXUMU, * )
    REAL, INTENT(IN OUT)                     :: ZZ( MXCMU, * )
    REAL, INTENT(IN)                         :: ZPLK0( MXCMU, * )
    REAL, INTENT(IN OUT)                     :: ZPLK1( MXCMU, * )
    REAL, INTENT(OUT)                        :: UUM( MXUMU, MXULV )
    
    
    
!     ..
!     .. Array Arguments ..
    
    
    
!     ..
!     .. Local Scalars ..
    
    LOGICAL :: NEGUMU
    INTEGER :: IQ, IU, JQ, LC, LU, LYREND, LYRSTR, LYU
    REAL :: BNDDFU, BNDDIR, BNDINT, DENOM, DFUINT, DTAU, DTAU1,  &
        DTAU2, EXP0, EXP1, EXP2, EXPN, F0N, F1N, FACT, PALINT, PLKINT, SGN
!     ..
!     .. Intrinsic Functions ..
    
    INTRINSIC ABS, EXP
!     ..
    
!                          ** Incorporate constants of integration into
!                          ** interpolated eigenvectors
    DO  LC = 1, NCUT
      
      DO  IQ = 1, NSTR
        
        DO  IU = 1, NUMU
          GU( IU, IQ, LC ) = GU( IU, IQ, LC ) * LL( IQ, LC )
        END DO
        
      END DO
      
    END DO
!                           ** Loop over levels at which intensities
!                           ** are desired ('user output levels')
    DO  LU = 1, NTAU
      
      IF( FBEAM > 0.0 ) EXP0  = EXP( -UTAUPR( LU ) / UMU0 )
      LYU  = LAYRU( LU )
!                              ** Loop over polar angles at which
!                              ** intensities are desired
      DO  IU = 1, NUMU
        
        IF( LYRCUT .AND. LYU > NCUT ) CYCLE
        
        NEGUMU = UMU( IU ) < 0.0
        
        IF( NEGUMU ) THEN
          
          LYRSTR = 1
        LYREND = LYU - 1
        SGN    = -1.0
        
      ELSE
        
        LYRSTR = LYU + 1
      LYREND = NCUT
      SGN    = 1.0
      
    END IF
!                          ** For downward intensity, integrate from top
!                          ** to LYU-1 in Eq. S1(8); for upward,
!                          ** integrate from bottom to LYU+1 in S1(9)
    PALINT = 0.0
    PLKINT = 0.0
    
  DO  LC = LYRSTR, LYREND
  
  DTAU = DTAUCP( LC )
  EXP1 = EXP( ( UTAUPR(LU) - TAUCPR(LC-1) ) / UMU( IU ) )
  EXP2 = EXP( ( UTAUPR(LU) - TAUCPR(LC)   ) / UMU( IU ) )
  
  IF( PLANK .AND. MAZIM == 0 ) THEN
    
!                          ** Eqs. STWL(36b,c, 37b,c)
    
    F0N = SGN * ( EXP1 - EXP2 )
    
    F1N = SGN * ( ( TAUCPR( LC-1 ) + UMU( IU ) ) * EXP1 -  &
        ( TAUCPR( LC )   + UMU( IU ) ) * EXP2 )
    
    PLKINT = PLKINT + Z0U( IU,LC )*F0N + Z1U( IU,LC )*F1N
    
  END IF
  
  
  IF( FBEAM > 0.0 ) THEN
    
    DENOM  = 1. + UMU( IU ) / UMU0
    
    IF( ABS( DENOM ) < 0.0001 ) THEN
!                                                   ** L'Hospital limit
      EXPN   = ( DTAU / UMU0 )*EXP0
      
    ELSE
      
      EXPN   = ( EXP1*EXPBEA( LC-1 ) - EXP2*EXPBEA( LC ) ) * SGN / DENOM
      
    END IF
    
    PALINT = PALINT + ZBEAM( IU, LC )*EXPN
    
  END IF
  
!                                                   ** KK is negative
  DO  IQ = 1, NN
    
    WK( IQ ) = EXP( KK( IQ,LC )*DTAU )
    DENOM  = 1.0 + UMU( IU )*KK( IQ, LC )
    
    IF( ABS( DENOM ) < 0.0001 ) THEN
!                                                   ** L'Hospital limit
      EXPN   = DTAU / UMU( IU )*EXP2
      
    ELSE
      
      EXPN   = SGN*( EXP1*WK( IQ ) - EXP2 ) / DENOM
      
    END IF
    
    PALINT = PALINT + GU( IU, IQ, LC )*EXPN
    
  END DO
  
!                                                   ** KK is positive
  DO  IQ = NN + 1, NSTR
    
    DENOM  = 1.0 + UMU( IU )*KK( IQ, LC )
    
    IF( ABS( DENOM ) < 0.0001 ) THEN
!                                                   ** L'Hospital limit
      EXPN  = -DTAU / UMU( IU )*EXP1
      
    ELSE
      
      EXPN  = SGN*( EXP1 - EXP2*WK( NSTR+1-IQ ) ) / DENOM
      
    END IF
    
    PALINT = PALINT + GU( IU, IQ, LC )*EXPN
    
  END DO
  
  
END DO
!                           ** Calculate contribution from user
!                           ** output level to next computational level

DTAU1  = UTAUPR( LU ) - TAUCPR( LYU - 1 )
DTAU2  = UTAUPR( LU ) - TAUCPR( LYU )

IF( ABS( DTAU1 ) < 1.E-6 .AND. NEGUMU ) GO TO  90
IF( ABS( DTAU2 ) < 1.E-6 .AND. (.NOT.NEGUMU ) ) GO TO  90

IF( NEGUMU )      EXP1  = EXP( DTAU1/UMU( IU ) )
IF( .NOT.NEGUMU ) EXP2  = EXP( DTAU2/UMU( IU ) )

IF( FBEAM > 0.0 ) THEN
  
  DENOM  = 1. + UMU( IU ) / UMU0
  
  IF( ABS( DENOM ) < 0.0001 ) THEN
    
    EXPN   = ( DTAU1 / UMU0 )*EXP0
    
  ELSE IF( NEGUMU ) THEN
    
    EXPN  = ( EXP0 - EXPBEA( LYU-1 )*EXP1 ) / DENOM
    
  ELSE
    
    EXPN  = ( EXP0 - EXPBEA( LYU )*EXP2 ) / DENOM
    
  END IF
  
  PALINT = PALINT + ZBEAM( IU, LYU )*EXPN
  
END IF

!                                                   ** KK is negative
DTAU  = DTAUCP( LYU )

DO  IQ = 1, NN
  
  DENOM  = 1. + UMU( IU )*KK( IQ, LYU )
  
  IF( ABS( DENOM ) < 0.0001 ) THEN
    
    EXPN = -DTAU2 / UMU( IU )*EXP2
    
  ELSE IF( NEGUMU ) THEN
    
    EXPN = ( EXP( -KK( IQ,LYU ) * DTAU2 ) -  &
        EXP(  KK( IQ,LYU ) * DTAU  ) * EXP1 ) / DENOM
    
  ELSE
    
    EXPN = ( EXP( -KK( IQ,LYU ) * DTAU2 ) - EXP2 ) / DENOM
    
  END IF
  
  PALINT = PALINT + GU( IU, IQ, LYU )*EXPN
  
END DO

!                                                   ** KK is positive
DO  IQ = NN + 1, NSTR
  
  DENOM  = 1. + UMU( IU )*KK( IQ, LYU )
  
  IF( ABS( DENOM ) < 0.0001 ) THEN
    
    EXPN   = -DTAU1 / UMU( IU )*EXP1
    
  ELSE IF( NEGUMU ) THEN
    
    EXPN = ( EXP( -KK( IQ,LYU ) * DTAU1 ) - EXP1 ) / DENOM
    
  ELSE
    
    EXPN = ( EXP( -KK( IQ,LYU ) * DTAU1 ) -  &
        EXP( -KK( IQ,LYU ) * DTAU  ) * EXP2 ) / DENOM
    
  END IF
  
  PALINT = PALINT + GU( IU, IQ, LYU )*EXPN
  
END DO


IF( PLANK .AND. MAZIM == 0 ) THEN
  
!                            ** Eqs. STWL (35-37) with tau-sub-n-1
!                            ** replaced by tau for upward, and
!                            ** tau-sub-n replaced by tau for downward
!                            ** directions
  
  IF( NEGUMU ) THEN
    
    EXPN  = EXP1
    FACT  = TAUCPR( LYU - 1 ) + UMU( IU )
    
  ELSE
    
    EXPN  = EXP2
    FACT  = TAUCPR( LYU ) + UMU( IU )
    
  END IF
  
  F0N  = 1. - EXPN
  F1N  = UTAUPR( LU ) + UMU( IU ) - FACT * EXPN
  
  PLKINT = PLKINT + Z0U( IU, LYU )*F0N + Z1U( IU, LYU )*F1N
  
END IF

!                            ** Calculate intensity components
!                            ** attenuated at both boundaries.
!                            ** NOTE: no azimuthal intensity
!                            ** component for isotropic surface
90       CONTINUE
BNDINT = 0.0

IF( NEGUMU .AND. MAZIM == 0 ) THEN
  
  BNDINT = (FISOT + TPLANK) * EXP( UTAUPR(LU ) / UMU(IU) )
  
  
ELSE IF( .NOT.NEGUMU ) THEN
  
  IF( LYRCUT .OR. ( LAMBER.AND.MAZIM > 0 ) ) GO TO  140
  
  DO  JQ = NN + 1, NSTR
    WK( JQ ) = EXP( -KK( JQ,NLYR )*DTAUCP( NLYR ) )
  END DO
  
  BNDDFU = 0.0
  
  DO  IQ = NN, 1, -1
    
    DFUINT = 0.0
    DO  JQ = 1, NN
      DFUINT = DFUINT + GC( IQ, JQ, NLYR )*LL( JQ, NLYR )
    END DO
    
    DO  JQ = NN + 1, NSTR
      DFUINT = DFUINT + GC( IQ, JQ, NLYR )* LL( JQ, NLYR )*WK( JQ )
    END DO
    
    IF( FBEAM > 0.0 ) DFUINT = DFUINT + ZZ( IQ, NLYR )*EXPBEA( NLYR )
    
    DFUINT = DFUINT + DELM0 * ( ZPLK0( IQ, NLYR ) +  &
        ZPLK1( IQ,NLYR ) *TAUCPR( NLYR ) )
    BNDDFU = BNDDFU + ( 1.+DELM0 ) * RMU(IU,NN+1-IQ)  &
        * CMU(NN+1-IQ) * CWT(NN+1-IQ)* DFUINT
  END DO
  
  BNDDIR = 0.0
  IF( FBEAM > 0.0 ) BNDDIR = UMU0*FBEAM / PI*RMU( IU, 0 )* EXPBEA( NLYR )
  
  BNDINT = ( BNDDFU + BNDDIR + DELM0 * EMU(IU) * BPLANK )  &
      * EXP( (UTAUPR(LU)-TAUCPR(NLYR)) / UMU(IU) )
  
END IF

140       CONTINUE

UUM( IU, LU ) = PALINT + PLKINT + BNDINT

END DO

END DO


RETURN
END SUBROUTINE USRINT

!************************************************************************

REAL FUNCTION  XIFUNC( UMU1, UMU2, UMU3, TAU )

!          Calculates Xi function of EQ. STWL (72)

!                    I N P U T   V A R I A B L E S

!        TAU         optical thickness of the layer

!        UMU1,2,3    cosine of zenith angle_1, _2, _3

!   Called by- SECSCA
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..


REAL, INTENT(IN)                         :: UMU1
REAL, INTENT(IN)                         :: UMU2
REAL, INTENT(IN)                         :: UMU3
REAL, INTENT(IN)                         :: TAU

!     ..
!     .. Local Scalars ..

REAL :: EXP1, X1, X2
!     ..
!     .. Intrinsic Functions ..

INTRINSIC EXP
!     ..


X1     = 1. / UMU1 - 1. / UMU2
X2     = 1. / UMU1 - 1. / UMU3

EXP1 = EXP( -TAU/UMU1 )

IF( UMU2 == UMU3 .AND. UMU1 == UMU2 ) THEN
  
  XIFUNC = TAU*TAU * EXP1 / ( 2.*UMU1*UMU2 )
  
ELSE IF( UMU2 == UMU3 .AND. UMU1 /= UMU2 ) THEN
  
  XIFUNC = ( ( TAU - 1./X1 ) * EXP( -TAU/UMU2 ) + EXP1 / X1 )  &
      / ( X1*UMU1*UMU2 )
  
ELSE IF( UMU2 /= UMU3 .AND. UMU1 == UMU2 ) THEN
  
  XIFUNC = ( ( EXP( -TAU/UMU3 ) - EXP1 ) / X2 - TAU * EXP1 )  &
      / ( X2*UMU1*UMU2 )
  
ELSE IF( UMU2 /= UMU3 .AND. UMU1 == UMU3 ) THEN
  
  XIFUNC = ( ( EXP( -TAU/UMU2 ) - EXP1 ) / X1 - TAU * EXP1 )  &
      / ( X1*UMU1*UMU2 )
  
ELSE
  
  XIFUNC = ( ( EXP( -TAU/UMU3 ) - EXP1 ) / X2 -  &
      (   EXP( -TAU/UMU2 ) - EXP1 ) / X1 ) / ( X2*UMU1*UMU2 )
  
END IF


RETURN
END FUNCTION  XIFUNC

! ******************************************************************
! ********** DISORT service routines ************************
! ******************************************************************

SUBROUTINE CHEKIN( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,  &
    WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,  &
    NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,  &
PHI0, FISOT, LAMBER, ALBEDO  BTEMP, TTEMP,  &
      TEMIS, PLANK, ONLYFL, DELTAM, CORINT, ACCUR,  &
      TAUC, MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM,  &
      MXCLY, MXULV, MXUMU, MXCMU, MXPHI, MXSQT )
  
!           Checks the input dimensions and variables
  
!   Calls- WRTBAD, WRTDIM, DREF, ERRMSG
!   Called by- DISORT
! --------------------------------------------------------------------
  
!     .. Scalar Arguments ..
  
  
  INTEGER, INTENT(IN)                      :: NLYR
  REAL, INTENT(IN)                         :: DTAUC( MAXCLY )
  REAL, INTENT(IN)                         :: SSALB( MAXCLY )
  INTEGER, INTENT(IN)                      :: NMOM
  REAL, INTENT(IN)                         :: PMOM( 0:MAXMOM, MAXCLY )
  REAL, INTENT(IN)                         :: TEMPER( 0:MAXCLY )
  REAL, INTENT(IN)                         :: WVNMLO
  REAL, INTENT(IN)                         :: WVNMHI
  LOGICAL, INTENT(IN)                      :: USRTAU
  INTEGER, INTENT(IN)                      :: NTAU
  REAL, INTENT(IN OUT)                     :: UTAU( MAXULV )
  INTEGER, INTENT(IN)                      :: NSTR
  LOGICAL, INTENT(IN)                      :: USRANG
  INTEGER, INTENT(IN)                      :: NUMU
  REAL, INTENT(IN OUT)                     :: UMU( MAXUMU )
  INTEGER, INTENT(IN)                      :: NPHI
  REAL, INTENT(IN)                         :: PHI( MAXPHI )
  INTEGER, INTENT(IN)                      :: IBCND
  REAL, INTENT(IN)                         :: FBEAM
  REAL, INTENT(IN)                         :: UMU0
  REAL, INTENT(IN)                         :: PHI0
  REAL, INTENT(IN)                         :: FISOT
  LOGICAL, INTENT(IN)                      :: LAMBER
  REAL, INTENT(IN OUT)                     :: ALBEDO  BT
    REAL, INTENT(IN)                         :: TTEMP
    REAL, INTENT(IN)                         :: TEMIS
    LOGICAL, INTENT(IN)                      :: PLANK
    LOGICAL, INTENT(IN OUT)                  :: ONLYFL
    LOGICAL, INTENT(IN OUT)                  :: DELTAM
    LOGICAL, INTENT(IN OUT)                  :: CORINT
    REAL, INTENT(IN)                         :: ACCUR
    REAL, INTENT(IN)                         :: TAUC( 0:MXCLY )
    INTEGER, INTENT(IN)                      :: MAXCLY
    INTEGER, INTENT(IN)                      :: MAXULV
    INTEGER, INTENT(IN)                      :: MAXUMU
    INTEGER, INTENT(IN)                      :: MAXPHI
    INTEGER, INTENT(IN)                      :: MAXMOM
    INTEGER, INTENT(IN)                      :: MXCLY
    INTEGER, INTENT(IN)                      :: MXULV
    INTEGER, INTENT(IN)                      :: MXUMU
    INTEGER, INTENT(IN)                      :: MXCMU
    INTEGER, INTENT(IN)                      :: MXPHI
    INTEGER, INTENT(IN)                      :: MXSQT
    
    
    REAL :: ALBEDO  BTEMP,  &
          
!     ..
!     .. Array Arguments ..
      
      
!     ..
!     .. Local Scalars ..
      
      LOGICAL :: INPERR
      INTEGER :: IRMU, IU, J, K, LC, LU, NUMSQT
      REAL :: FLXALB, RMU, YESSCT
!     ..
!     .. External Functions ..
      
      LOGICAL :: WRTBAD, WRTDIM
      REAL :: DREF
      EXTERNAL  WRTBAD, WRTDIM, DREF
!     ..
!     .. External Subroutines ..
      
      EXTERNAL  ERRMSG
!     ..
!     .. Intrinsic Functions ..
      
      INTRINSIC ABS, MAX, MOD
!     ..
      
      
      INPERR = .FALSE.
      
      IF( NSTR < 2 .OR. MOD( NSTR,2 ) /= 0 ) INPERR = WRTBAD( 'NSTR' )
      
      IF( NSTR == 2 ) CALL ERRMSG( 'CHEKIN--2 streams not recommended; '//  &
          'use specialized 2-stream code TWOSTR instead', .True.)
      
      IF( NLYR < 1 ) INPERR = WRTBAD( 'NLYR' )
      
      IF( NLYR > MAXCLY ) INPERR = WRTBAD( 'MAXCLY' )
      
      YESSCT = 0.0
      
      DO  LC = 1, NLYR
        
        IF( DTAUC( LC ) < 0.0 ) INPERR = WRTBAD( 'DTAUC' )
        
        IF( SSALB( LC ) < 0.0 .OR. SSALB( LC ) > 1.0 )  &
            INPERR = WRTBAD( 'SSALB' )
        
        YESSCT = YESSCT + SSALB( LC )
        
        IF( PLANK .AND. IBCND /= 1 ) THEN
          
          IF( LC == 1 .AND. TEMPER( 0 ) < 0.0 ) INPERR = WRTBAD( 'TEMPER' )
          
          IF( TEMPER( LC ) < 0.0 ) INPERR = WRTBAD( 'TEMPER' )
          
        END IF
        
      END DO
      
      IF( NMOM < 0 .OR. ( YESSCT > 0.0 .AND. NMOM < NSTR ) )  &
          INPERR = WRTBAD( 'NMOM' )
      
      IF( MAXMOM < NMOM ) INPERR = WRTBAD( 'MAXMOM' )
      
      DO  LC = 1, NLYR
        
        DO  K = 0, NMOM
          
          IF( PMOM( K,LC ) < -1.0 .OR. PMOM( K,LC ) > 1.0 )  &
              INPERR = WRTBAD( 'PMOM' )
          
        END DO
        
      END DO
      
      IF( IBCND == 1 ) THEN
        
        IF( MAXULV < 2 ) INPERR = WRTBAD( 'MAXULV' )
        
      ELSE IF( USRTAU ) THEN
        
        IF( NTAU < 1 ) INPERR = WRTBAD( 'NTAU' )
        
        IF( MAXULV < NTAU ) INPERR = WRTBAD( 'MAXULV' )
        
        DO  LU = 1, NTAU
          
          IF( ABS( UTAU( LU )-TAUC( NLYR ) ) <= 1.E-4 )  &
              UTAU( LU ) = TAUC( NLYR )
          
          IF( UTAU( LU ) < 0.0 .OR. UTAU( LU ) > TAUC( NLYR ) ) THEN
            INPERR = WRTBAD( 'UTAU' )
          END IF
          
        END DO
        
      ELSE
        
        IF( MAXULV < NLYR + 1 ) INPERR = WRTBAD( 'MAXULV' )
        
      END IF
      
      
      IF( USRANG ) THEN
        
        IF( NUMU < 0 ) INPERR = WRTBAD( 'NUMU' )
        
        IF( .NOT.ONLYFL .AND. NUMU == 0 ) INPERR = WRTBAD( 'NUMU' )
        
        IF( NUMU > MAXUMU ) INPERR = WRTBAD( 'MAXUMU' )
        
        IF( IBCND == 1 .AND. 2*NUMU > MAXUMU ) INPERR = WRTBAD( 'MAXUMU' )
        
        DO  IU = 1, NUMU
          
          IF( UMU( IU ) < -1.0 .OR. UMU( IU ) > 1.0 .OR.  &
              UMU( IU ) == 0.0 ) INPERR = WRTBAD( 'UMU' )
          
          IF( IBCND == 1 .AND. UMU( IU ) < 0.0 ) INPERR = WRTBAD( 'UMU' )
          
          IF( IU > 1 ) THEN
            
            IF( UMU( IU ) < UMU( IU-1 ) ) INPERR = WRTBAD( 'UMU' )
            
          END IF
          
        END DO
        
      ELSE
        
        IF( MAXUMU < NSTR ) INPERR = WRTBAD( 'MAXUMU' )
        
      END IF
      
      
      IF( .NOT.ONLYFL .AND. IBCND /= 1 ) THEN
        
        IF( NPHI <= 0 ) INPERR = WRTBAD( 'NPHI' )
        
        IF( NPHI > MAXPHI ) INPERR = WRTBAD( 'MAXPHI' )
        
        DO  J = 1, NPHI
          
          IF( PHI( J ) < 0.0 .OR. PHI( J ) > 360.0 ) INPERR = WRTBAD( 'PHI' )
          
        END DO
        
      END IF
      
      
      IF( IBCND < 0 .OR. IBCND > 1 ) INPERR = WRTBAD( 'IBCND' )
      
      IF( IBCND == 0 ) THEN
        
        IF( FBEAM < 0.0 ) INPERR = WRTBAD( 'FBEAM' )
        
        IF( FBEAM > 0.0 .AND. ( UMU0 <= 0.0 .OR. UMU0 > 1.0 ) )  &
            INPERR = WRTBAD( 'UMU0' )
        
        IF( FBEAM > 0.0 .AND. ( PHI0 < 0.0 .OR. PHI0 > 360.0 ) )  &
            INPERR = WRTBAD( 'PHI0' )
        
        IF( FISOT < 0.0 ) INPERR = WRTBAD( 'FISOT' )
        
        IF( LAMBER ) THEN
          
          IF( ALBEDO < 0.0 .OR. ALBEDO > 1.0 )  &
                INPERR = WRTBAD( 'ALBEDO' )
            
          ELSE
!                    ** Make sure flux albedo at dense mesh of incident
!                    ** angles does not assume unphysical values
            
            DO  IRMU = 0, 100
              
              RMU  = IRMU*0.01
              FLXALB = DREF( WVNMLO, WVNMHI, RMU )
              
              IF( FLXALB < 0.0 .OR. FLXALB > 1.0 )  &
                  INPERR = WRTBAD( 'FUNCTION BDREF' )
              
            END DO
            
          END IF
          
          
        ELSE IF( IBCND == 1 ) THEN
          
          IF( ALBEDO < 0.0 .OR. ALBEDO > 1.0 )  &
                INPERR = WRTBAD( 'ALBEDO' )
            
          END IF
          
          
          IF( PLANK .AND. IBCND /= 1 ) THEN
            
            IF( WVNMLO < 0.0 .OR. WVNMHI <= WVNMLO )  &
                INPERR = WRTBAD( 'WVNMLO,HI' )
            
            IF( TEMIS < 0.0 .OR. TEMIS > 1.0 ) INPERR = WRTBAD( 'TEMIS' )
            
            IF( BTEMP < 0.0 ) INPERR = WRTBAD( 'BTEMP' )
            
            IF( TTEMP < 0.0 ) INPERR = WRTBAD( 'TTEMP' )
            
          END IF
          
          
          IF( ACCUR < 0.0 .OR. ACCUR > 1.E-2 ) INPERR = WRTBAD( 'ACCUR' )
          
          IF( MXCLY < NLYR ) INPERR = WRTDIM( 'MXCLY', NLYR )
          
          IF( IBCND /= 1 ) THEN
            
            IF( USRTAU .AND. MXULV < NTAU ) INPERR = WRTDIM( 'MXULV', NTAU )
            
            IF( .NOT.USRTAU .AND. MXULV < NLYR + 1 )  &
                INPERR = WRTDIM( 'MXULV', NLYR + 1 )
            
          ELSE
            
            IF( MXULV < 2 ) INPERR = WRTDIM( 'MXULV', 2 )
            
          END IF
          
          IF( MXCMU < NSTR ) INPERR = WRTDIM( 'MXCMU', NSTR )
          
          IF( USRANG .AND. MXUMU < NUMU ) INPERR = WRTDIM( 'MXUMU', NUMU )
          
          IF( USRANG .AND. IBCND == 1 .AND. MXUMU < 2*NUMU )  &
              INPERR = WRTDIM( 'MXUMU', 2*NUMU )
          
          IF( .NOT.USRANG .AND. MXUMU < NSTR )  &
              INPERR = WRTDIM( 'MXUMU', NSTR )
          
          IF( .NOT.ONLYFL .AND. IBCND /= 1 .AND. MXPHI < NPHI )  &
              INPERR = WRTDIM( 'MXPHI', NPHI )
          
          NUMSQT = 2*MAX( 100, NSTR )
          IF( MXSQT < NUMSQT ) INPERR = WRTDIM( 'MXSQT', NUMSQT )
          
          IF( INPERR )  &
              CALL ERRMSG( 'DISORT--input and/or dimension errors', .True. )
          
          IF( PLANK ) THEN
            
            DO  LC = 1, NLYR
              
              IF( ABS( TEMPER( LC )-TEMPER( LC-1 ) ) > 10.0 )  &
                  CALL ERRMSG('CHEKIN--vertical temperature step may'  &
                  //' be too large for good accuracy', .False. )
            END DO
            
          END IF
          
          IF( .NOT.CORINT .AND. .NOT.ONLYFL .AND. FBEAM > 0.0 .AND.  &
              YESSCT > 0.0 .AND. DELTAM )  &
              CALL ERRMSG( 'CHEKIN--intensity correction is off; '//  &
              'intensities may be less accurate', .False. )
          
          
          RETURN
        END SUBROUTINE CHEKIN
        
!************************************************************************
        
        REAL FUNCTION  DREF( WVNMLO, WVNMHI, MU )
        
!        Flux albedo for given angle of incidence, given
!        a bidirectional reflectivity.
        
!  INPUT :   MU      Cosine of incidence angle
        
!            WVNMLO  Lower wavenumber (inv-cm) of spectral interval
        
!            WVNMHI  Upper wavenumber (inv-cm) of spectral interval
        
        
!  INTERNAL VARIABLES :
        
!       NMUG   :  Number of angle cosine quadrature points on (-1,1)
!                 for integrating bidirectional reflectivity to get
!                 directional emissivity (it is necessary to use a
!                 quadrature set distinct from the computational angles,
!                 because the computational angles may not be dense
!                 enough -- i.e. 'NSTR' may be too small -- to give an
!                 accurate approximation for the integration).
        
!       GMU    :  The 'NMUG' angle cosine quadrature points on (0,1)
        
!       GWT    :  The 'NMUG' angle cosine quadrature weights on (0,1)
        
!   Called by- CHEKIN
!   Calls- QGAUSN, ERRMSG, BDREF
! +--------------------------------------------------------------------+
        
!     .. Parameters ..
        
        
        REAL, INTENT(IN OUT)                     :: WVNMLO
        REAL, INTENT(IN OUT)                     :: WVNMHI
        REAL, INTENT(IN)                         :: MU
        
        INTEGER, PARAMETER :: NMUG = 50
!     ..
!     .. Scalar Arguments ..
        
        
!     ..
!     .. Local Scalars ..
        
        LOGICAL :: PASS1
        INTEGER :: JG, K
        REAL :: PI, SUM
!     ..
!     .. Local Arrays ..
        
        REAL :: GMU( NMUG ), GWT( NMUG )
!     ..
!     .. External Functions ..
        
        REAL :: BDREF
        EXTERNAL  BDREF
!     ..
!     .. External Subroutines ..
        
        EXTERNAL  ERRMSG, QGAUSN
!     ..
!     .. Intrinsic Functions ..
        
        INTRINSIC ABS, ASIN
!     ..
        SAVE      PASS1, GMU, GWT, PI
        DATA      PASS1 / .True. /
        
        
        IF( PASS1 ) THEN
          
          PASS1 = .FALSE.
          PI   = 2.*ASIN( 1.0 )
          
          CALL QGAUSN( NMUG/2, GMU, GWT )
          
          DO  K = 1, NMUG / 2
            GMU( K + NMUG/2 ) = -GMU( K )
            GWT( K + NMUG/2 ) = GWT( K )
          END DO
          
        END IF
        
        IF( ABS( MU ) > 1.0 )  &
            CALL ERRMSG( 'DREF--input argument error(s)',.True. )
        
        DREF = 0.0
        
!                       ** Loop over azimuth angle difference
        DO  JG = 1, NMUG
          
          SUM  = 0.0
!                       ** Loop over angle of reflection
          DO  K = 1, NMUG / 2
            SUM  = SUM + GWT( K )*GMU( K )*  &
                BDREF( WVNMLO, WVNMHI, GMU( K ), MU, PI*GMU( JG ) )
          END DO
          
          DREF = DREF + GWT( JG )*SUM
          
        END DO
        
        IF( DREF < 0.0 .OR. DREF > 1.0 )  &
            CALL ERRMSG( 'DREF--albedo value not in (0,1)',.False. )
        
        
        RETURN
      END FUNCTION  DREF
      
!************************************************************************
      
      SUBROUTINE LEPOLY( NMU, M, MAXMU, TWONM1, MU, SQT, YLM )
      
!       Computes the normalized associated Legendre polynomial,
!       defined in terms of the associated Legendre polynomial
!       Plm = P-sub-l-super-m as
      
!             Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)
      
!       for fixed order m and all degrees from l = m to TWONM1.
!       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available
!       from a prior call to the routine.
      
!       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
!                  High-Order Associated Legendre Polynomials,
!                  J. Quant. Spectrosc. Radiat. Transfer 10,
!                  557-562, 1970.  (hereafter D/A)
      
!       METHOD: Varying degree recurrence relationship.
      
!       NOTES:
!       (1) The D/A formulas are transformed by setting M=n-1; L=k-1.
!       (2) Assumes that routine is called first with  M = 0, then with
!           M = 1, etc. up to  M = TWONM1.
      
      
!  I N P U T     V A R I A B L E S:
      
!       NMU    :  Number of arguments of YLM
      
!       M      :  Order of YLM
      
!       MAXMU  :  First dimension of YLM
      
!       TWONM1 :  Max degree of YLM
      
!       MU(i)  :  Arguments of YLM (i = 1 to NMU)
      
!       SQT(k) :  Square root of k
      
!       If M.GT.0, YLM(M-1,i) for i = 1 to NMU is assumed to exist
!       from a prior call.
      
      
!  O U T P U T     V A R I A B L E:
      
!       YLM(l,i) :  l = M to TWONM1, normalized associated Legendre
!                   polynomials evaluated at argument MU(i)
      
!   Called by- DISORT, ALBTRN
! +-------------------------------------------------------------------+
      
!     .. Scalar Arguments ..
      
      
      INTEGER, INTENT(IN)                      :: NMU
      INTEGER, INTENT(IN OUT)                  :: M
      INTEGER, INTENT(IN OUT)                  :: MAXMU
      INTEGER, INTENT(IN)                      :: TWONM1
      REAL, INTENT(IN)                         :: MU( * )
      REAL, INTENT(IN)                         :: SQT( * )
      REAL, INTENT(OUT)                        :: YLM( 0:MAXMU, * )
      
!     ..
!     .. Array Arguments ..
      
      
!     ..
!     .. Local Scalars ..
      
      INTEGER :: I, L
      REAL :: TMP1, TMP2
!     ..
      
      
      IF( M == 0 ) THEN
!                             ** Upward recurrence for ordinary
!                             ** Legendre polynomials
        DO  I = 1, NMU
          YLM( 0, I ) = 1.0
          YLM( 1, I ) = MU( I )
        END DO
        
        
        DO  L = 2, TWONM1
          
          DO  I = 1, NMU
            YLM( L, I ) = ( ( 2*L - 1 )*MU( I )*YLM( L-1, I ) -  &
                ( L - 1 )*YLM( L-2, I ) ) / L
          END DO
          
        END DO
        
        
      ELSE
        
        DO  I = 1, NMU
!                               ** Y-sub-m-super-m; derived from
!                               ** D/A Eqs. (11,12), STWL(58c)
          
          YLM( M, I ) = - SQT( 2*M - 1 ) / SQT( 2*M )*  &
              SQRT( 1.- MU(I)**2 )*YLM( M-1, I )
          
!                              ** Y-sub-(m+1)-super-m; derived from
!                              ** D/A Eqs.(13,14) using Eqs.(11,12),
!                              ** STWL(58f)
          
          YLM( M+1, I ) = SQT( 2*M + 1 )*MU( I )*YLM( M, I )
          
        END DO
        
!                                   ** Upward recurrence; D/A EQ.(10),
!                                   ** STWL(58a)
        DO  L = M + 2, TWONM1
          
          TMP1  = SQT( L - M )*SQT( L + M )
          TMP2  = SQT( L - M - 1 )*SQT( L + M - 1 )
          
          DO  I = 1, NMU
            YLM( L, I ) = ( ( 2*L - 1 )*MU( I )*YLM( L-1, I ) -  &
                TMP2*YLM( L-2, I ) ) / TMP1
          END DO
          
        END DO
        
      END IF
      
      
      RETURN
    END SUBROUTINE LEPOLY
    
!************************************************************************
!sergio this has been usurped; however, call this for the SLFTST
!      REAL FUNCTION PLKAVG( WNUMLO, WNUMHI, T )
    
    REAL FUNCTION PLKAVG_ORIG( WNUMLO, WNUMHI, T )
    
!        Computes Planck function integrated between two wavenumbers
    
!  INPUT :  WNUMLO : Lower wavenumber (inv cm) of spectral interval
    
!           WNUMHI : Upper wavenumber
    
!           T      : Temperature (K)
    
!  OUTPUT : PLKAVG : Integrated Planck function ( Watts/sq m )
!                      = Integral (WNUMLO to WNUMHI) of
!                        2h c**2  nu**3 / ( EXP(hc nu/kT) - 1)
!                        (where h=Plancks constant, c=speed of
!                         light, nu=wavenumber, T=temperature,
!                         and k = Boltzmann constant)
    
!  Reference : Specifications of the Physical World: New Value
!                 of the Fundamental Constants, Dimensions/N.B.S.,
!                 Jan. 1974
    
!  Method :  For WNUMLO close to WNUMHI, a Simpson-rule quadrature
!            is done to avoid ill-conditioning; otherwise
    
!            (1)  For WNUMLO or WNUMHI small,
!                 integral(0 to WNUMLO/HI) is calculated by expanding
!                 the integrand in a power series and integrating
!                 term by term;
    
!            (2)  Otherwise, integral(WNUMLO/HI to INFINITY) is
!                 calculated by expanding the denominator of the
!                 integrand in powers of the exponential and
!                 integrating term by term.
    
!  Accuracy :  At least 6 significant digits, assuming the
!              physical constants are infinitely accurate
    
!  ERRORS WHICH ARE NOT TRAPPED:
    
!      * power or exponential series may underflow, giving no
!        significant digits.  This may or may not be of concern,
!        depending on the application.
    
!      * Simpson-rule special case is skipped when denominator of
!        integrand will cause overflow.  In that case the normal
!        procedure is used, which may be inaccurate if the
!        wavenumber limits (WNUMLO, WNUMHI) are close together.
    
!  LOCAL VARIABLES
    
!        A1,2,... :  Power series coefficients
!        C2       :  h * c / k, in units cm*K (h = Plancks constant,
!                      c = speed of light, k = Boltzmann constant)
!        D(I)     :  Exponential series expansion of integral of
!                       Planck function from WNUMLO (i=1) or WNUMHI
!                       (i=2) to infinity
!        EPSIL    :  Smallest number such that 1+EPSIL .GT. 1 on
!                       computer
!        EX       :  EXP( - V(I) )
!        EXM      :  EX**M
!        MMAX     :  No. of terms to take in exponential series
!        MV       :  Multiples of V(I)
!        P(I)     :  Power series expansion of integral of
!                       Planck function from zero to WNUMLO (I=1) or
!                       WNUMHI (I=2)
!        PI       :  3.14159...
!        SIGMA    :  Stefan-Boltzmann constant (W/m**2/K**4)
!        SIGDPI   :  SIGMA / PI
!        SMALLV   :  Number of times the power series is used (0,1,2)
!        V(I)     :  C2 * (WNUMLO(I=1) or WNUMHI(I=2)) / temperature
!        VCUT     :  Power-series cutoff point
!        VCP      :  Exponential series cutoff points
!        VMAX     :  Largest allowable argument of EXP function
    
!   Called by- DISORT
!   Calls- D1MACH, ERRMSG
! ----------------------------------------------------------------------
    
!     .. Parameters ..
    
!      REAL plkavg_orig
    
    REAL, INTENT(IN)                         :: WNUMLO
    REAL, INTENT(IN)                         :: WNUMHI
    REAL, INTENT(IN)                         :: T
    REAL :: plkavg
    
    REAL, PARAMETER :: A1 = 1. / 3.
    REAL, PARAMETER :: A2 = -1. / 8.
    REAL, PARAMETER :: A3 = 1. / 60.
    REAL, PARAMETER :: A4 = -1. / 5040.
    REAL, PARAMETER :: A5 = 1. / 272160.
    REAL, PARAMETER :: A6 = -1. / 13305600.
!     ..
!     .. Scalar Arguments ..
    
    
!     ..
!     .. Local Scalars ..
    
    INTEGER :: I, K, M, MMAX, N, SMALLV
    REAL :: C2, CONC, DEL, EPSIL, EX, EXM, HH, MV, OLDVAL, PI,  &
        SIGDPI, SIGMA, VAL, VAL0, VCUT, VMAX, VSQ, X
!     ..
!     .. Local Arrays ..
    
    REAL :: D( 2 ), P( 2 ), V( 2 ), VCP( 7 )
!     ..
!     .. External Functions ..
    
    REAL :: D1MACH
    EXTERNAL  D1MACH
!     ..
!     .. External Subroutines ..
    
    EXTERNAL  ERRMSG
!     ..
!     .. Intrinsic Functions ..
    
    INTRINSIC ABS, ASIN, EXP, LOG, MOD
!     ..
!     .. Statement Functions ..
    
    REAL :: PLKF
!     ..
    SAVE      PI, CONC, VMAX, EPSIL, SIGDPI
    
    DATA      C2 / 1.438786 / , SIGMA / 5.67032E-8 / , VCUT / 1.5 / ,  &
        VCP / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
    DATA      PI / 0.0 /
    
!     .. Statement Function definitions ..
    
    PLKF( X ) = X**3 / ( EXP( X ) - 1 )
!     ..
    
    
    IF( PI == 0.0 ) THEN
      
      PI     = 2.*ASIN( 1.0 )
      VMAX   = LOG( D1MACH( 2 ) )
      EPSIL  = D1MACH( 4 )
      SIGDPI = SIGMA / PI
      CONC   = 15. / PI**4
      
    END IF
    
    
    IF( T < 0.0 .OR. WNUMHI <= WNUMLO .OR. WNUMLO < 0. )  &
        CALL ERRMSG('PLKAVG--temperature or wavenums. wrong',.TRUE.)
    
    
    IF( T < 1.E-4 ) THEN
      
      PLKAVG = 0.0
      plkavg_orig=plkavg
      RETURN
      
    END IF
    
    
    V( 1 ) = C2*WNUMLO / T
    V( 2 ) = C2*WNUMHI / T
    
    IF( V( 1 ) > EPSIL .AND. V( 2 ) < VMAX .AND.  &
          ( WNUMHI - WNUMLO ) / WNUMHI < 1.E-2 ) THEN
      
!                          ** Wavenumbers are very close.  Get integral
!                          ** by iterating Simpson rule to convergence.
      
      HH     = V( 2 ) - V( 1 )
      OLDVAL = 0.0
      VAL0   = PLKF( V( 1 ) ) + PLKF( V( 2 ) )
      
      DO  N = 1, 10
        
        DEL  = HH / ( 2*N )
        VAL  = VAL0
        
        DO  K = 1, 2*N - 1
          VAL  = VAL + 2*( 1 + MOD( K,2 ) )* PLKF( V( 1 ) + K*DEL )
        END DO
        
        VAL  = DEL / 3.*VAL
        IF( ABS( ( VAL - OLDVAL ) / VAL ) <= 1.E-6 ) GO TO  30
        OLDVAL = VAL
        
      END DO
      
      CALL ERRMSG( 'PLKAVG--Simpson rule didnt converge',.FALSE.)
      
      30    CONTINUE
      
      PLKAVG = SIGDPI * T**4 * CONC * VAL
      plkavg_orig=plkavg
      RETURN
      
    END IF
    
!                          *** General case ***
    SMALLV = 0
    
    DO  I = 1, 2
      
      IF( V( I ) < VCUT ) THEN
!                                   ** Use power series
        SMALLV = SMALLV + 1
        VSQ    = V( I )**2
        P( I ) = CONC*VSQ*V( I )*( A1 +  &
            V( I )*( A2 + V( I )*( A3 + VSQ*( A4 + VSQ*( A5 + VSQ*A6 ) ) ) ) )
        
      ELSE
!                      ** Use exponential series
        MMAX  = 0
!                                ** Find upper limit of series
        40       CONTINUE
        MMAX  = MMAX + 1
        
        IF( V(I) < VCP( MMAX ) ) GO TO  40
        
        EX     = EXP( - V(I) )
        EXM    = 1.0
        D( I ) = 0.0
        
        DO  M = 1, MMAX
          MV     = M*V( I )
          EXM    = EX*EXM
          D( I ) = D( I ) + EXM*( 6.+ MV*( 6.+ MV*( 3.+ MV ) ) ) / M**4
        END DO
        
        D( I ) = CONC*D( I )
        
      END IF
      
    END DO
    
!                              ** Handle ill-conditioning
    IF( SMALLV == 2 ) THEN
!                                    ** WNUMLO and WNUMHI both small
      PLKAVG = P( 2 ) - P( 1 )
      
    ELSE IF( SMALLV == 1 ) THEN
!                                    ** WNUMLO small, WNUMHI large
      PLKAVG = 1.- P( 1 ) - D( 2 )
      
    ELSE
!                                    ** WNUMLO and WNUMHI both large
      PLKAVG = D( 1 ) - D( 2 )
      
    END IF
    
    PLKAVG = SIGDPI * T**4 * PLKAVG
    
    IF( PLKAVG == 0.0 )  &
        CALL ERRMSG('PLKAVG--returns zero; possible underflow', .FALSE.)
    
    plkavg_orig=plkavg
    
    RETURN
  END FUNCTION PLKAVG_ORIG
  
!************************************************************************
  
  
  SUBROUTINE PRAVIN( UMU, NUMU, MXUMU, UTAU, NTAU, U0U )
  
!        Print azimuthally averaged intensities at user angles
  
!   Called by- DISORT
  
!     LENFMT   Max number of polar angle cosines UMU that can be
!              printed on one line, as set in FORMAT statement
! --------------------------------------------------------------------
  
  
  REAL, INTENT(OUT)                        :: UMU( NUMU )
  INTEGER, INTENT(IN)                      :: NUMU
  INTEGER, INTENT(IN OUT)                  :: MXUMU
  REAL, INTENT(IN OUT)                     :: UTAU( NTAU )
  INTEGER, INTENT(IN)                      :: NTAU
  REAL, INTENT(OUT)                        :: U0U( MXUMU, * )
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/scatterparam.f90'
  
!     .. Scalar Arguments ..
  
  
!     ..
!     .. Array Arguments ..
  
  
!     ..
!     .. Local Scalars ..
  
  INTEGER :: IU, IUMAX, IUMIN, LENFMT, LU, NP, NPASS
!     ..
!     .. Intrinsic Functions ..
  
  INTRINSIC MIN
!     ..
  
  
  IF( NUMU < 1 )  RETURN
  
  WRITE( KSTDWARN, '(//,A)' )  &
      ' *******  AZIMUTHALLY AVERAGED INTENSITIES ' //  &
      '(at user polar angles)  ********'
  
  LENFMT = 8
  NPASS  = 1 + (NUMU-1) / LENFMT
  
  WRITE( KSTDWARN,'(/,A,/,A)') '   Optical   Polar Angle Cosines',  &
      '     Depth'
  
  DO  NP = 1, NPASS
    
    IUMIN  = 1 + LENFMT * ( NP - 1 )
    IUMAX  = MIN( LENFMT*NP, NUMU )
    WRITE( KSTDWARN,'(/,10X,8F14.5)') ( UMU(IU), IU = IUMIN, IUMAX )
    
    DO  LU = 1, NTAU
      WRITE( KSTDWARN, '(0P,F10.4,1P,8E14.4)' ) UTAU( LU ),  &
          ( U0U( IU,LU ), IU = IUMIN, IUMAX )
    END DO
    
  END DO
  
  
  RETURN
END SUBROUTINE PRAVIN

!************************************************************************
!sergio got rid of maxmom
!      SUBROUTINE PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER,
!     &                   WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,
!     &                   NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT,
!     &                   LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM,
!     &                   PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT,
!     &                   OPRIM, TAUC, TAUCPR, MAXMOM, PRTMOM )

SUBROUTINE PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER,  &
    WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,  &
    NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT,  &
LAMBER, ALBEDO  BTEMP, TTEMP, TEMIS, DELTAM,  &
      PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, OPRIM, TAUC, TAUCPR, PRTMOM )
  
!        Print values of input variables
  
!   Called by- DISORT
! --------------------------------------------------------------------
  
!     .. Scalar Arguments ..
  
  
  INTEGER, INTENT(IN)                      :: NLYR
  REAL, INTENT(IN OUT)                     :: DTAUC( * )
  REAL, INTENT(IN OUT)                     :: DTAUCP( * )
  REAL, INTENT(IN)                         :: SSALB( * )
  INTEGER, INTENT(IN)                      :: NMOM
  REAL, INTENT(OUT)                        :: PMOM( 0:MAXMOM, * )
  REAL, INTENT(IN OUT)                     :: TEMPER( 0:* )
  REAL, INTENT(IN OUT)                     :: WVNMLO
  REAL, INTENT(IN OUT)                     :: WVNMHI
  INTEGER, INTENT(OUT)                     :: NTAU
  REAL, INTENT(OUT)                        :: UTAU( * )
  INTEGER, INTENT(IN)                      :: NSTR
  INTEGER, INTENT(OUT)                     :: NUMU
  REAL, INTENT(OUT)                        :: UMU( * )
  INTEGER, INTENT(OUT)                     :: NPHI
  REAL, INTENT(OUT)                        :: PHI( * )
  INTEGER, INTENT(IN OUT)                  :: IBCND
  REAL, INTENT(IN)                         :: FBEAM
  REAL, INTENT(IN)                         :: UMU0
  REAL, INTENT(IN)                         :: PHI0
  REAL, INTENT(IN)                         :: FISOT
  LOGICAL, INTENT(IN)                      :: LAMBER
  NO TYPE, INTENT(IN OUT)                  :: ALBEDO  BT
    REAL, INTENT(IN)                         :: TTEMP
    REAL, INTENT(IN)                         :: TEMIS
    LOGICAL, INTENT(IN)                      :: DELTAM
    LOGICAL, INTENT(IN)                      :: PLANK
    LOGICAL, INTENT(IN OUT)                  :: ONLYFL
    LOGICAL, INTENT(IN)                      :: CORINT
    REAL, INTENT(IN OUT)                     :: ACCUR
    REAL, INTENT(IN OUT)                     :: FLYR( * )
    LOGICAL, INTENT(IN)                      :: LYRCUT
    REAL, INTENT(IN OUT)                     :: OPRIM( * )
    REAL, INTENT(IN OUT)                     :: TAUC( 0:* )
    REAL, INTENT(IN OUT)                     :: TAUCPR( 0:* )
    LOGICAL, INTENT(IN)                      :: PRTMOM
    IMPLICIT NONE
    
    INCLUDE '../INCLUDE/scatterparam.f90'
    
    
!sergio got rid of maxmom
!      INTEGER   IBCND, MAXMOM, NLYR, NMOM, NPHI, NSTR, NTAU, NUMU
    
    REAL :: ALBEDO  BTEMP,  &
          
!     ..
!     .. Array Arguments ..
      
      
!     ..
!     .. Local Scalars ..
      
      INTEGER :: IU, J, K, LC, LU
      REAL :: YESSCT
!     ..
      
      
      WRITE( KSTDWARN, '(/,A,I4,A,I4)' ) ' No. streams =', NSTR,  &
          '     No. computational layers =', NLYR
      
      IF( IBCND /= 1 ) WRITE( KSTDWARN, '(I4,A,10F10.4,/,(26X,10F10.4))' )  &
          NTAU, ' User optical depths :', ( UTAU(LU), LU = 1, NTAU )
      
      IF( .NOT.ONLYFL ) WRITE( KSTDWARN, '(I4,A,10F9.5,/,(31X,10F9.5))' )  &
          NUMU, ' User polar angle cosines :', ( UMU(IU), IU = 1, NUMU )
      
      IF( .NOT.ONLYFL .AND. IBCND /= 1 )  &
          WRITE( KSTDWARN, '(I4,A,10F9.2,/,(28X,10F9.2))' )  &
          NPHI,' User azimuthal angles :',( PHI(J), J = 1, NPHI )
      
      IF( .NOT.PLANK .OR. IBCND == 1 )  &
          WRITE( KSTDWARN, '(A)' ) ' No thermal emission'
      
      
      WRITE( KSTDWARN, '(A,I2)' ) ' Boundary condition flag: IBCND =', IBCND
      
      IF( IBCND == 0 ) THEN
        
        WRITE( KSTDWARN, '(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P,E11.3)' )  &
            '    Incident beam with intensity =', FBEAM,  &
            ' and polar angle cosine = ', UMU0, '  and azimuth angle =', PHI0,  &
            '    plus isotropic incident intensity =', FISOT
        
        IF( LAMBER ) WRITE( KSTDWARN, '(A,0P,F8.4)' )  &
        '    Bottom albedo (Lambertian) =', ALBEDO
          
          IF( .NOT.LAMBER ) WRITE( KSTDWARN, '(A)' )  &
              '    Bidirectional reflectivity at bottom'
          
          IF( PLANK ) WRITE( KSTDWARN, '(A,2F14.4,/,A,F10.2,A,F10.2,A,F8.4)' )  &
              '    Thermal emission in wavenumber interval :', WVNMLO, WVNMHI,  &
              '    Bottom temperature =', BTEMP,  &
              '    Top temperature =', TTEMP, '    Top emissivity =', TEMIS
          
        ELSE IF( IBCND == 1 ) THEN
          
          WRITE( KSTDWARN, '(A)' )  &
              '    Isotropic illumination from top and bottom'
          WRITE( KSTDWARN, '(A,0P,F8.4)' )  &
          '    Bottom albedo (Lambertian) =', ALBEDO
            
          END IF
          
          
          IF( DELTAM ) WRITE( KSTDWARN, '(A)' ) ' Uses delta-M method'
          IF( .NOT.DELTAM ) WRITE( KSTDWARN, '(A)' ) ' Does not use delta-M method'
          
          IF( CORINT ) WRITE( KSTDWARN, '(A)' ) ' Uses TMS/IMS method'
          IF( .NOT.CORINT ) WRITE( KSTDWARN,'(A)' ) ' Does not use TMS/IMS method'
          
          
          IF( IBCND == 1 ) THEN
            
            WRITE( KSTDWARN, '(A)' ) ' Calculate albedo and transmissivity of'//  &
                ' medium vs. incident beam angle'
            
          ELSE IF( ONLYFL ) THEN
            
            WRITE( KSTDWARN, '(A)' ) ' Calculate fluxes only'
            
          ELSE
            
            WRITE( KSTDWARN, '(A)' ) ' Calculate fluxes and intensities'
            
          END IF
          
          WRITE( KSTDWARN, '(A,1P,E11.2)' )  &
              ' Relative convergence criterion for azimuth series =', ACCUR
          
          IF( LYRCUT ) WRITE( KSTDWARN, '(A)' )  &
              ' Sets radiation = 0 below absorption optical depth 10'
          
          
!                                    ** Print layer variables
!                                    ** (to read, skip every other line)
          
          IF( PLANK ) WRITE( KSTDWARN, '(/,37X,A,3(/,2A))' )  &
              '<------------- Delta-M --------------->',  &
              '                   Total    Single                           ',  &
              'Total    Single',  &
              '       Optical   Optical   Scatter   Separated   ',  &
              'Optical   Optical   Scatter    Asymm',  &
              '         Depth     Depth    Albedo    Fraction     ',  &
              'Depth     Depth    Albedo   Factor   Temperature'
          
          IF( .NOT.PLANK ) WRITE( KSTDWARN, '(/,37X,A,3(/,2A))' )  &
              '<------------- Delta-M --------------->',  &
              '                   Total    Single                           ',  &
              'Total    Single',  &
              '       Optical   Optical   Scatter   Separated   ',  &
              'Optical   Optical   Scatter    Asymm',  &
              '         Depth     Depth    Albedo    Fraction     ',  &
              'Depth     Depth    Albedo   Factor'
          
          
          YESSCT = 0.0
          
          DO  LC = 1, NLYR
            
            YESSCT = YESSCT + SSALB( LC )
!                                       ** f90 nonadvancing I/O would
!                                       ** simplify this a lot (also the
!                                       ** two WRITEs above)
            IF( PLANK )  &
                WRITE( KSTDWARN,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4,F14.3)')  &
                LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),  &
                DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ),  &
                PMOM( 1,LC ), TEMPER( LC-1 )
            
            IF( .NOT.PLANK )  &
                WRITE( KSTDWARN,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4)' )  &
                LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),  &
                DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM( 1,LC )
          END DO
          
          IF( PLANK ) WRITE( KSTDWARN, '(85X,F14.3)' ) TEMPER( NLYR )
          
          
          IF( PRTMOM .AND. YESSCT > 0.0 ) THEN
            
            WRITE( KSTDWARN, '(/,A,I5)' ) ' Number of Phase Function Moments = ',  &
                NMOM + 1
            WRITE( KSTDWARN, '(A)' ) ' Layer   Phase Function Moments'
            
            DO  LC = 1, NLYR
              
              IF( SSALB( LC ) > 0.0 )  &
                  WRITE( KSTDWARN, '(I6,10F11.6,/,(6X,10F11.6))' )  &
                  LC, ( PMOM( K, LC ), K = 0, NMOM )
            END DO
            
          END IF
          
          
          RETURN
        END SUBROUTINE PRTINP
        
!************************************************************************
!sergio comment out maxulv,maxumu
!      SUBROUTINE PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
!     &                   MAXUMU )
        
        SUBROUTINE PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI)
        
        
        REAL, INTENT(OUT)                        :: UU( MAXUMU, MAXULV, * )
        REAL, INTENT(OUT)                        :: UTAU( * )
        INTEGER, INTENT(IN)                      :: NTAU
        REAL, INTENT(OUT)                        :: UMU( * )
        INTEGER, INTENT(IN)                      :: NUMU
        REAL, INTENT(OUT)                        :: PHI( * )
        INTEGER, INTENT(IN)                      :: NPHI
        IMPLICIT NONE
        INCLUDE '../INCLUDE/scatterparam.f90'
        
!         Prints the intensity at user polar and azimuthal angles
        
!     All arguments are DISORT input or output variables
        
!   Called by- DISORT
        
!     LENFMT   Max number of azimuth angles PHI that can be printed
!                on one line, as set in FORMAT statement
! +-------------------------------------------------------------------+
        
!     .. Scalar Arguments ..
        
!sergio comment out maxulv,maxumu
!      INTEGER   MAXULV, MAXUMU, NPHI, NTAU, NUMU
        
!     ..
!     .. Array Arguments ..
        
        
!     ..
!     .. Local Scalars ..
        
        INTEGER :: IU, J, JMAX, JMIN, LENFMT, LU, NP, NPASS
!     ..
!     .. Intrinsic Functions ..
        
        INTRINSIC MIN
!     ..
        
        
        IF( NPHI < 1 )  RETURN
        
        WRITE( KSTDWARN, '(//,A)' )  &
            ' *********  I N T E N S I T I E S  *********'
        
        LENFMT = 10
        NPASS  = 1 + (NPHI-1) / LENFMT
        
        WRITE( KSTDWARN, '(/,A,/,A,/,A)' )  &
            '             Polar   Azimuth angles (degrees)',  &
            '   Optical   Angle', '    Depth   Cosine'
        
        DO  LU = 1, NTAU
          
          DO  NP = 1, NPASS
            
            JMIN   = 1 + LENFMT * ( NP - 1 )
            JMAX   = MIN( LENFMT*NP, NPHI )
            
            WRITE( KSTDWARN, '(/,18X,10F11.2)' ) ( PHI(J), J = JMIN, JMAX )
            
            IF( NP == 1 ) WRITE( KSTDWARN, '(F10.4,F8.4,1P,10E11.3)' )  &
                UTAU(LU), UMU(1), (UU(1, LU, J), J = JMIN, JMAX)
            IF( NP > 1 ) WRITE( KSTDWARN, '(10X,F8.4,1P,10E11.3)' )  &
                UMU(1), (UU(1, LU, J), J = JMIN, JMAX)
            
            DO  IU = 2, NUMU
              WRITE( KSTDWARN, '(10X,F8.4,1P,10E11.3)' )  &
                  UMU( IU ), ( UU( IU, LU, J ), J = JMIN, JMAX )
            END DO
            
          END DO
          
        END DO
        
        
        RETURN
      END SUBROUTINE PRTINT
      
!************************************************************************
      
      SUBROUTINE QGAUSN( M, GMU, GWT )
      
!       Compute weights and abscissae for ordinary Gaussian quadrature
!       on the interval (0,1);  that is, such that
      
!           sum(i=1 to M) ( GWT(i) f(GMU(i)) )
      
!       is a good approximation to
      
!           integral(0 to 1) ( f(x) dx )
      
!   INPUT :    M       order of quadrature rule
      
!   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M)
!             GWT(I)   array of weights (I = 1 TO M)
      
!   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
!                   Integration, Academic Press, New York, pp. 87, 1975
      
!   METHOD:  Compute the abscissae as roots of the Legendre
!            polynomial P-sub-M using a cubically convergent
!            refinement of Newton's method.  Compute the
!            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
!            that Newton's method can very easily diverge; only a
!            very good initial guess can guarantee convergence.
!            The initial guess used here has never led to divergence
!            even for M up to 1000.
      
!   ACCURACY:  relative error no better than TOL or computer
!              precision (machine epsilon), whichever is larger
      
!   INTERNAL VARIABLES:
      
!    ITER      : number of Newton Method iterations
!    MAXIT     : maximum allowed iterations of Newton Method
!    PM2,PM1,P : 3 successive Legendre polynomials
!    PPR       : derivative of Legendre polynomial
!    P2PRI     : 2nd derivative of Legendre polynomial
!    TOL       : convergence criterion for Legendre poly root iteration
!    X,XI      : successive iterates in cubically-convergent version
!                of Newtons Method (seeking roots of Legendre poly.)
      
!   Called by- DREF, SETDIS, SURFAC
!   Calls- D1MACH, ERRMSG
! +-------------------------------------------------------------------+
      
!     .. Scalar Arguments ..
      
      
      INTEGER, INTENT(IN)                      :: M
      REAL, INTENT(OUT)                        :: GMU( M )
      REAL, INTENT(OUT)                        :: GWT( M )
      
!     ..
!     .. Array Arguments ..
      
      
!     ..
!     .. Local Scalars ..
      
      INTEGER :: ITER, K, LIM, MAXIT, NN, NP1
      REAL :: CONA, PI, T
      REAL :: EN, NNP1, ONE, P, P2PRI, PM1, PM2, PPR, PROD,  &
          TMP, TOL, TWO, X, XI
!     ..
!     .. External Functions ..
      
      REAL :: D1MACH
      EXTERNAL  D1MACH
!     ..
!     .. External Subroutines ..
      
      EXTERNAL  ERRMSG
!     ..
!     .. Intrinsic Functions ..
      
      INTRINSIC ABS, ASIN, COS, FLOAT, MOD, TAN
!     ..
      SAVE      PI, TOL
      
      DATA      PI / 0.0 /, MAXIT / 1000 /, ONE / 1.0 /, TWO / 2.0 /
      
      IF( PI == 0.0 ) THEN
        
        PI   = 2.*ASIN( 1.0 )
        TOL  = 10.*D1MACH( 4 )
        
      END IF
      
      
      IF( M < 1 ) CALL ERRMSG( 'QGAUSN--Bad value of M',.True.)
      
      IF( M == 1 ) THEN
        
        GMU( 1 ) = 0.5
        GWT( 1 ) = 1.0
        RETURN
        
      END IF
      
      EN   = M
      NP1  = M + 1
      NNP1 = M*NP1
      CONA = FLOAT( M - 1 ) / ( 8*M**3 )
      
      LIM  = M / 2
      
      DO  K = 1, LIM
!                                        ** Initial guess for k-th root
!                                        ** of Legendre polynomial, from
!                                        ** Davis/Rabinowitz (2.7.3.3a)
        T  = ( 4*K - 1 )*PI / ( 4*M + 2 )
        X  = COS( T + CONA / TAN( T ) )
        ITER = 0
!                                        ** Upward recurrence for
!                                        ** Legendre polynomials
        10    CONTINUE
        ITER   = ITER + 1
        PM2    = ONE
        PM1    = X
        
        DO  NN = 2, M
          P    = ( ( 2*NN - 1 )*X*PM1 - ( NN - 1 )*PM2 ) / NN
          PM2  = PM1
          PM1  = P
        END DO
!                                              ** Newton Method
        TMP    = ONE / ( ONE - X**2 )
        PPR    = EN*( PM2 - X*P )*TMP
        P2PRI  = ( TWO*X*PPR - NNP1*P )*TMP
        XI     = X - ( P / PPR )*( ONE + ( P / PPR )*P2PRI / ( TWO*PPR ) )
        
!                                              ** Check for convergence
        IF( ABS( XI - X ) > TOL ) THEN
          IF( ITER > MAXIT )  &
              CALL ERRMSG( 'QGAUSN--max iteration count',.True.)
          
          X  = XI
          GO TO  10
          
        END IF
!                             ** Iteration finished--calculate weights,
!                             ** abscissae for (-1,1)
        GMU( K ) = -X
        GWT( K ) = TWO / ( TMP*( EN*PM2 )**2 )
        GMU( NP1 - K ) = -GMU( K )
        GWT( NP1 - K ) = GWT( K )
      END DO
!                                    ** Set middle abscissa and weight
!                                    ** for rules of odd order
      IF( MOD( M,2 ) /= 0 ) THEN
        
        GMU( LIM + 1 ) = 0.0
        PROD   = ONE
        
        DO  K = 3, M, 2
          PROD   = PROD * K / ( K - 1 )
        END DO
        
        GWT( LIM + 1 ) = TWO / PROD**2
      END IF
      
!                                        ** Convert from (-1,1) to (0,1)
      DO  K = 1, M
        GMU( K ) = 0.5*GMU( K ) + 0.5
        GWT( K ) = 0.5*GWT( K )
      END DO
      
      
      RETURN
    END SUBROUTINE QGAUSN
    
!************************************************************************
    
    REAL FUNCTION RATIO( A, B )
    
!        Calculate ratio  A/B  with over- and under-flow protection
!        (thanks to Prof. Jeff Dozier for some suggestions here).
!        Since this routine takes two logs, it is no speed demon,
!        but it is invaluable for comparing results from two runs
!        of a program under development.
    
!        NOTE:  In Fortran90, built-in functions TINY and HUGE
!               can replace the D1MACH calls.
    
!   Called by- DISORT
!   Calls- D1MACH
! +-------------------------------------------------------------------+
    
!     .. Scalar Arguments ..
    
    
    REAL, INTENT(IN OUT)                     :: A
    REAL, INTENT(IN OUT)                     :: B
    
!     ..
!     .. Local Scalars ..
    
    LOGICAL :: PASS1
    REAL :: ABSA, ABSB, HUGE, POWA, POWB, POWMAX, POWMIN, TINY
!     ..
!     .. External Functions ..
    
    REAL :: D1MACH
    EXTERNAL  D1MACH
!     ..
!     .. Intrinsic Functions ..
    
    INTRINSIC ABS, LOG10, SIGN
!     ..
!     .. Save statement ..
    
    SAVE      PASS1, TINY, HUGE, POWMAX, POWMIN
!     ..
!     .. Data statements ..
    
    DATA      PASS1 / .TRUE. /
!     ..
    
    
    IF( PASS1 ) THEN
      
      TINY   = D1MACH( 1 )
      HUGE   = D1MACH( 2 )
      POWMAX = LOG10( HUGE )
      POWMIN = LOG10( TINY )
      PASS1  = .FALSE.
      
    END IF
    
    
    IF( A == 0.0 ) THEN
      
      IF( B == 0.0 ) THEN
        
        RATIO  = 1.0
        
      ELSE
        
        RATIO  = 0.0
        
      END IF
      
      
    ELSE IF( B == 0.0 ) THEN
      
      RATIO  = SIGN( HUGE, A )
      
    ELSE
      
      ABSA   = ABS( A )
      ABSB   = ABS( B )
      POWA   = LOG10( ABSA )
      POWB   = LOG10( ABSB )
      
      IF( ABSA < TINY .AND. ABSB < TINY ) THEN
        
        RATIO  = 1.0
        
      ELSE IF( POWA - POWB >= POWMAX ) THEN
        
        RATIO  = HUGE
        
      ELSE IF( POWA - POWB <= POWMIN ) THEN
        
        RATIO  = TINY
        
      ELSE
        
        RATIO  = ABSA / ABSB
        
      END IF
!                      ** DONT use old trick of determining sign
!                      ** from A*B because A*B may (over/under)flow
      
      IF( ( A > 0.0 .AND. B < 0.0 ) .OR.  &
          ( A < 0.0 .AND. B > 0.0 ) ) RATIO = -RATIO
      
    END IF
    
    
    RETURN
  END FUNCTION RATIO
  
!************************************************************************
  
  SUBROUTINE SLFTST( CORINT, ACCUR, ALBEDO  BTEMP, DELTAM, DTAUC,  &
        FBEAM, FISOT, IBCND, LAMBER, NLYR, PLANK, NPHI,  &
        NUMU, NSTR, NTAU, ONLYFL, PHI, PHI0, NMOM,  &
        PMOM, PRNT, PRNTU0, SSALB, TEMIS, TEMPER,  &
        TTEMP, UMU, USRANG, USRTAU, UTAU, UMU0, WVNMHI,  &
        WVNMLO, COMPAR, FLUP, RFLDIR, RFLDN, UU )
    
!       If  COMPAR = FALSE, save user input values that would otherwise
!       be destroyed and replace them with input values for self-test.
!       If  COMPAR = TRUE, compare self-test case results with correct
!       answers and restore user input values if test is passed.
    
!       (See file 'DISORT.doc' for variable definitions.)
    
    
!     I N T E R N A L    V A R I A B L E S:
    
!         ACC     Relative accuracy required for passing self-test
    
!         ERRORn  Relative errors in DISORT output variables
    
!         OK      Logical variable for determining failure of self-test
    
!         All variables ending in 'S' are temporary 'S'torage for input
    
!   Called by- DISORT
!   Calls- TSTBAD, ERRMSG
! +-------------------------------------------------------------------+
    
!     .. Scalar Arguments ..
    
    
    LOGICAL, INTENT(IN OUT)                  :: CORINT
    REAL, INTENT(IN OUT)                     :: ACCUR
    REAL, INTENT(IN OUT)                     :: ALBEDO  BT
      LOGICAL, INTENT(IN OUT)                  :: DELTAM
      REAL, INTENT(IN OUT)                     :: DTAUC
      REAL, INTENT(IN OUT)                     :: FBEAM
      REAL, INTENT(IN OUT)                     :: FISOT
      INTEGER, INTENT(IN OUT)                  :: IBCND
      LOGICAL, INTENT(IN OUT)                  :: LAMBER
      INTEGER, INTENT(IN OUT)                  :: NLYR
      LOGICAL, INTENT(IN OUT)                  :: PLANK
      INTEGER, INTENT(IN OUT)                  :: NPHI
      INTEGER, INTENT(IN OUT)                  :: NUMU
      INTEGER, INTENT(IN OUT)                  :: NSTR
      INTEGER, INTENT(IN OUT)                  :: NTAU
      LOGICAL, INTENT(IN OUT)                  :: ONLYFL
      REAL, INTENT(IN OUT)                     :: PHI
      REAL, INTENT(IN OUT)                     :: PHI0
      INTEGER, INTENT(IN OUT)                  :: NMOM
      REAL, INTENT(IN OUT)                     :: PMOM( 0:* )
      LOGICAL, INTENT(IN OUT)                  :: PRNT( * )
      LOGICAL, INTENT(IN OUT)                  :: PRNTU0( * )
      REAL, INTENT(IN OUT)                     :: SSALB
      REAL, INTENT(IN OUT)                     :: TEMIS
      REAL, INTENT(IN OUT)                     :: TEMPER( 0:* )
      REAL, INTENT(IN OUT)                     :: TTEMP
      REAL, INTENT(IN OUT)                     :: UMU
      LOGICAL, INTENT(IN OUT)                  :: USRANG
      LOGICAL, INTENT(IN OUT)                  :: USRTAU
      REAL, INTENT(IN OUT)                     :: UTAU
      REAL, INTENT(IN OUT)                     :: UMU0
      REAL, INTENT(IN OUT)                     :: WVNMHI
      REAL, INTENT(IN OUT)                     :: WVNMLO
      LOGICAL, INTENT(IN OUT)                  :: COMPAR
      REAL, INTENT(IN)                         :: FLUP
      REAL, INTENT(IN)                         :: RFLDIR
      REAL, INTENT(IN)                         :: RFLDN
      REAL, INTENT(IN)                         :: UU
      
      
      REAL :: ALBEDO  BTEMP,  &
            
!     ..
!     .. Array Arguments ..
        
        
        
!     ..
!     .. Local Scalars ..
        
        LOGICAL :: CORINS, DELTAS, LAMBES, OK, ONLYFS, PLANKS, USRANS, USRTAS
        INTEGER :: I, IBCNDS, N, NLYRS, NMOMS, NPHIS, NSTRS, NTAUS, NUMUS
        REAL :: ACC, ACCURS, ALBEDS, BTEMPS, DTAUCS, ERROR1, ERROR2,  &
            ERROR3, ERROR4, FBEAMS, FISOTS, PHI0S, PHIS, SSALBS,  &
            TEMISS, TTEMPS, UMU0S, UMUS, UTAUS, WVNMHS, WVNMLS
!     ..
!     .. Local Arrays ..
        
        LOGICAL :: PRNTS( 5 ), PRNU0S( 2 )
        REAL :: PMOMS( 0:4 ), TEMPES( 0:1 )
!     ..
!     .. External Functions ..
        
        LOGICAL :: TSTBAD
        EXTERNAL  TSTBAD
!     ..
!     .. External Subroutines ..
        
        EXTERNAL  ERRMSG
!     ..
!     .. Intrinsic Functions ..
        
        INTRINSIC ABS
!     ..
        SAVE
        DATA      ACC / 1.E-4 /
        
        IF( .NOT.COMPAR ) THEN
!                                     ** Save user input values
          NLYRS  = NLYR
          DTAUCS = DTAUC
          SSALBS = SSALB
          
          DO  N = 0, 4
            PMOMS( N ) = PMOM( N )
          END DO
          
          NSTRS  = NSTR
          NMOMS  = NMOM
          USRANS = USRANG
          NUMUS  = NUMU
          UMUS   = UMU
          USRTAS = USRTAU
          NTAUS  = NTAU
          UTAUS  = UTAU
          NPHIS  = NPHI
          PHIS   = PHI
          IBCNDS = IBCND
          FBEAMS = FBEAM
          UMU0S  = UMU0
          PHI0S  = PHI0
          FISOTS = FISOT
          LAMBES = LAMBER
          ALBEDS = ALBEDO
            DELTAS = DELTAM
            ONLYFS = ONLYFL
            CORINS = CORINT
            ACCURS = ACCUR
            PLANKS = PLANK
            WVNMLS = WVNMLO
            WVNMHS = WVNMHI
            BTEMPS = BTEMP
            TTEMPS = TTEMP
            TEMISS = TEMIS
            TEMPES( 0 ) = TEMPER( 0 )
            TEMPES( 1 ) = TEMPER( 1 )
            
            DO  I = 1, 5
              PRNTS( I ) = PRNT( I )
            END DO
            
            DO  I = 1, 2
              PRNU0S( I ) = PRNTU0( I )
            END DO
            
!                                     ** Set input values for self-test
            NSTR   = 4
            NLYR   = 1
            DTAUC  = 1.0
            SSALB  = 0.9
            NMOM   = 4
!                          ** Haze L moments
            PMOM( 0 ) = 1.0
            PMOM( 1 ) = 0.8042
            PMOM( 2 ) = 0.646094
            PMOM( 3 ) = 0.481851
            PMOM( 4 ) = 0.359056
            USRANG = .TRUE.
            NUMU   = 1
            UMU    = 0.5
            USRTAU = .TRUE.
            NTAU   = 1
            UTAU   = 0.5
            NPHI   = 1
            PHI    = 90.0
            IBCND  = 0
            FBEAM  = 3.14159265
            UMU0   = 0.866
            PHI0   = 0.0
            FISOT  = 1.0
            LAMBER = .TRUE.
            ALBEDO = 0.7
              DELTAM = .TRUE.
              ONLYFL = .FALSE.
              CORINT = .TRUE.
              ACCUR  = 1.E-4
              PLANK  = .TRUE.
              WVNMLO = 0.0
              WVNMHI = 50000.
              BTEMP  = 300.0
              TTEMP  = 100.0
              TEMIS  = 0.8
              TEMPER( 0 ) = 210.0
              TEMPER( 1 ) = 200.0
              
              DO  I = 1, 5
                PRNT( I ) = .FALSE.
              END DO
              
              DO  I = 1, 2
                PRNTU0( I ) = .FALSE.
              END DO
              
              
            ELSE
!                                    ** Compare test case results with
!                                    ** correct answers and abort if bad
              OK     = .TRUE.
              
!  for testing purposes
!      print *,   COMPAR, CORINT, DELTAM, LAMBER, ONLYFL, PLANK, USRANG,
!     &           USRTAU
!      print *,   IBCND, NLYR, NMOM, NPHI, NSTR, NTAU, NUMU
!      print *,   ACCUR, ALBEDO, BTEMP, DTAUC, FBEAM, FISOT, FLUP, PHI
!      print *,    PHI0, RFLDIR, RFLDN, SSALB, TEMIS, TTEMP, UMU, UMU0
!      print *,    UTAU, UU, WVNMHI, WVNMLO
!         stop
              
              ERROR1 = ( UU - 47.865571 ) / 47.865571
              ERROR2 = ( RFLDIR - 1.527286 ) / 1.527286
              ERROR3 = ( RFLDN - 28.372225 ) / 28.372225
              ERROR4 = ( FLUP - 152.585284 ) / 152.585284
              
              IF( ABS( ERROR1 ) > ACC ) OK = TSTBAD( 'UU', ERROR1 )
              
              IF( ABS( ERROR2 ) > ACC ) OK = TSTBAD( 'RFLDIR', ERROR2 )
              
              IF( ABS( ERROR3 ) > ACC ) OK = TSTBAD( 'RFLDN', ERROR3 )
              
              IF( ABS( ERROR4 ) > ACC ) OK = TSTBAD( 'FLUP', ERROR4 )
              
              IF( .NOT.OK ) CALL ERRMSG( 'DISORT--self-test failed', .True. )
              
!                                      ** Restore user input values
              NLYR   = NLYRS
              DTAUC  = DTAUCS
              SSALB  = SSALBS
              
              DO  N = 0, 4
                PMOM( N ) = PMOMS( N )
              END DO
              
              NSTR   = NSTRS
              NMOM   = NMOMS
              USRANG = USRANS
              NUMU   = NUMUS
              UMU    = UMUS
              USRTAU = USRTAS
              NTAU   = NTAUS
              UTAU   = UTAUS
              NPHI   = NPHIS
              PHI    = PHIS
              IBCND  = IBCNDS
              FBEAM  = FBEAMS
              UMU0   = UMU0S
              PHI0   = PHI0S
              FISOT  = FISOTS
              LAMBER = LAMBES
              ALBEDO = ALBEDS
                DELTAM = DELTAS
                ONLYFL = ONLYFS
                CORINT = CORINS
                ACCUR  = ACCURS
                PLANK  = PLANKS
                WVNMLO = WVNMLS
                WVNMHI = WVNMHS
                BTEMP  = BTEMPS
                TTEMP  = TTEMPS
                TEMIS  = TEMISS
                TEMPER( 0 ) = TEMPES( 0 )
                TEMPER( 1 ) = TEMPES( 1 )
                
                DO  I = 1, 5
                  PRNT( I ) = PRNTS( I )
                END DO
                
                DO  I = 1, 2
                  PRNTU0( I ) = PRNU0S( I )
                END DO
                
              END IF
              
              RETURN
            END SUBROUTINE SLFTST
            
!************************************************************************
            
            SUBROUTINE ZEROAL( ND1, EXPBEA, FLYR, OPRIM, PHASA, PHAST, PHASM,  &
                TAUCPR, XR0, XR1, ND2, CMU, CWT, PSI0, PSI1, WK, Z0, Z1, ZJ,  &
                ND3, YLM0, ND4, ARRAY, CC, EVECC,  &
                ND5, GL, ND6, YLMC,  &
                ND7, YLMU, ND8, KK, LL, ZZ, ZPLK0, ZPLK1,  &
                ND9, GC, ND10, LAYRU, UTAUPR,  &
                ND11, GU, ND12, Z0U, Z1U, ZBEAM,  &
                ND13, EVAL, ND14, AMB, APB,  &
                ND15, IPVT, Z, ND16, RFLDIR, RFLDN, FLUP, UAVG, DFDT,  &
                ND17, ALBMED, TRNMED, ND18, U0U,  &
                ND19, UU )
            
!         ZERO ARRAYS; NDn is dimension of all arrays following
!         it in the argument list
            
!   Called by- DISORT
! --------------------------------------------------------------------
            
!     .. Scalar Arguments ..
            
            
            INTEGER, INTENT(IN)                      :: ND1
            REAL, INTENT(OUT)                        :: EXPBEA( * )
            REAL, INTENT(OUT)                        :: FLYR( * )
            REAL, INTENT(OUT)                        :: OPRIM( * )
            REAL, INTENT(OUT)                        :: PHASA( * )
            REAL, INTENT(OUT)                        :: PHAST( * )
            REAL, INTENT(OUT)                        :: PHASM( * )
            REAL, INTENT(OUT)                        :: TAUCPR( * )
            REAL, INTENT(OUT)                        :: XR0( * )
            REAL, INTENT(OUT)                        :: XR1( * )
            INTEGER, INTENT(IN)                      :: ND2
            REAL, INTENT(OUT)                        :: CMU( * )
            REAL, INTENT(OUT)                        :: CWT( * )
            REAL, INTENT(OUT)                        :: PSI0( * )
            REAL, INTENT(OUT)                        :: PSI1( * )
            REAL, INTENT(OUT)                        :: WK( * )
            REAL, INTENT(OUT)                        :: Z0( * )
            REAL, INTENT(OUT)                        :: Z1( * )
            REAL, INTENT(OUT)                        :: ZJ( * )
            INTEGER, INTENT(IN)                      :: ND3
            REAL, INTENT(OUT)                        :: YLM0( * )
            INTEGER, INTENT(IN)                      :: ND4
            REAL, INTENT(OUT)                        :: ARRAY( * )
            REAL, INTENT(OUT)                        :: CC( * )
            REAL, INTENT(OUT)                        :: EVECC( * )
            INTEGER, INTENT(IN)                      :: ND5
            REAL, INTENT(OUT)                        :: GL( * )
            INTEGER, INTENT(IN)                      :: ND6
            REAL, INTENT(OUT)                        :: YLMC( * )
            INTEGER, INTENT(IN)                      :: ND7
            REAL, INTENT(OUT)                        :: YLMU( * )
            INTEGER, INTENT(IN)                      :: ND8
            REAL, INTENT(OUT)                        :: KK( * )
            REAL, INTENT(OUT)                        :: LL( * )
            REAL, INTENT(OUT)                        :: ZZ( * )
            REAL, INTENT(OUT)                        :: ZPLK0( * )
            REAL, INTENT(OUT)                        :: ZPLK1( * )
            INTEGER, INTENT(IN)                      :: ND9
            REAL, INTENT(OUT)                        :: GC( * )
            INTEGER, INTENT(IN)                      :: ND10
            INTEGER, INTENT(OUT)                     :: LAYRU( * )
            REAL, INTENT(OUT)                        :: UTAUPR( * )
            INTEGER, INTENT(IN)                      :: ND11
            REAL, INTENT(OUT)                        :: GU( * )
            INTEGER, INTENT(IN)                      :: ND12
            REAL, INTENT(OUT)                        :: Z0U( * )
            REAL, INTENT(OUT)                        :: Z1U( * )
            REAL, INTENT(OUT)                        :: ZBEAM( * )
            INTEGER, INTENT(IN)                      :: ND13
            REAL, INTENT(OUT)                        :: EVAL( * )
            INTEGER, INTENT(IN)                      :: ND14
            REAL, INTENT(OUT)                        :: AMB( * )
            REAL, INTENT(OUT)                        :: APB( * )
            INTEGER, INTENT(IN)                      :: ND15
            INTEGER, INTENT(OUT)                     :: IPVT( * )
            REAL, INTENT(OUT)                        :: Z( * )
            INTEGER, INTENT(IN)                      :: ND16
            REAL, INTENT(OUT)                        :: RFLDIR( * )
            REAL, INTENT(OUT)                        :: RFLDN( * )
            REAL, INTENT(OUT)                        :: FLUP( * )
            REAL, INTENT(OUT)                        :: UAVG( * )
            REAL, INTENT(OUT)                        :: DFDT( * )
            INTEGER, INTENT(IN)                      :: ND17
            REAL, INTENT(OUT)                        :: ALBMED( * )
            REAL, INTENT(OUT)                        :: TRNMED( * )
            INTEGER, INTENT(IN)                      :: ND18
            REAL, INTENT(OUT)                        :: U0U( * )
            INTEGER, INTENT(IN)                      :: ND19
            REAL, INTENT(OUT)                        :: UU( * )
            
!     ..
!     .. Array Arguments ..
            
            
            
!     ..
!     .. Local Scalars ..
            
            INTEGER :: N
!     ..
            
            
            DO  N = 1, ND1
              EXPBEA( N ) = 0.0
              FLYR( N )   = 0.0
              OPRIM( N )  = 0.0
              PHASA( N )  = 0.0
              PHAST( N )  = 0.0
              PHASM( N )  = 0.0
              TAUCPR( N ) = 0.0
              XR0( N )    = 0.0
              XR1( N )    = 0.0
            END DO
            
            DO  N = 1, ND2
              CMU( N )  = 0.0
              CWT( N )  = 0.0
              PSI0( N ) = 0.0
              PSI1( N ) = 0.0
              WK( N )   = 0.0
              Z0( N )   = 0.0
              Z1( N )   = 0.0
              ZJ( N )   = 0.0
            END DO
            
            DO  N = 1, ND3
              YLM0( N ) = 0.0
            END DO
            
            DO  N = 1, ND4
              ARRAY( N ) = 0.0
              CC( N )    = 0.0
              EVECC( N ) = 0.0
            END DO
            
            DO  N = 1, ND5
              GL( N ) = 0.0
            END DO
            
            DO  N = 1, ND6
              YLMC( N ) = 0.0
            END DO
            
            DO  N = 1, ND7
              YLMU( N ) = 0.0
            END DO
            
            DO  N = 1, ND8
              KK( N )    = 0.0
              LL( N )    = 0.0
              ZZ( N )    = 0.0
              ZPLK0( N ) = 0.0
              ZPLK1( N ) = 0.0
            END DO
            
            DO  N = 1, ND9
              GC( N ) = 0.0
            END DO
            
            DO  N = 1, ND10
              LAYRU( N )  = 0
              UTAUPR( N ) = 0.0
            END DO
            
            DO  N = 1, ND11
              GU( N ) = 0.0
            END DO
            
            DO  N = 1, ND12
              Z0U( N )   = 0.0
              Z1U( N )   = 0.0
              ZBEAM( N ) = 0.0
            END DO
            
            DO  N = 1, ND13
              EVAL( N ) = 0.0
            END DO
            
            DO  N = 1, ND14
              AMB( N ) = 0.0
              APB( N ) = 0.0
            END DO
            
            DO  N = 1, ND15
              IPVT( N ) = 0
              Z( N )    = 0.0
            END DO
            
            DO  N = 1, ND16
              RFLDIR( N ) = 0.
              RFLDN( N )  = 0.
              FLUP( N )   = 0.
              UAVG( N )   = 0.
              DFDT( N )   = 0.
            END DO
            
            DO  N = 1, ND17
              ALBMED( N ) = 0.
              TRNMED( N ) = 0.
            END DO
            
            DO  N = 1, ND18
              U0U( N ) = 0.
            END DO
            
            DO  N = 1, ND19
              UU( N ) = 0.
            END DO
            
            
            RETURN
          END SUBROUTINE ZEROAL
          
!************************************************************************
          
          SUBROUTINE ZEROIT( A, LENGTH )
          
!         Zeros a real array A having LENGTH elements
          
!   Called by- DISORT, ALBTRN, SOLVE1, SURFAC, SETMTX, SOLVE0, FLUXES
! --------------------------------------------------------------------
          
!     .. Scalar Arguments ..
          
          
          REAL, INTENT(OUT)                        :: A( LENGTH )
          INTEGER, INTENT(IN)                      :: LENGTH
          
!     ..
!     .. Array Arguments ..
          
          
!     ..
!     .. Local Scalars ..
          
          INTEGER :: L
!     ..
          
          
          DO  L = 1, LENGTH
            A( L ) = 0.0
          END DO
          
          
          RETURN
        END SUBROUTINE ZEROIT
        
! ******************************************************************
! ********** end of DISORT service routines ************************
! ******************************************************************
        
! ******************************************************************
! ********** IBCND=1 special case routines *************************
! ******************************************************************
        
        SUBROUTINE ALBTRN( ALBEDO  AMB, APB, ARRAY, B, BDR, CBAND, CC,  &
              CMU, CWT, DTAUCP, EVAL, EVECC, GL, GC, GU,  &
              IPVT, KK, LL, NLYR, NN, NSTR, NUMU, PRNT,  &
              TAUCPR, UMU, U0U, WK, YLMC, YLMU, Z, WKD, MI,  &
              MI9M2, MAXUMU, MXCMU, MXUMU, NNLYRI, SQT, ALBMED, TRNMED )
          
!    DISORT special case to get only albedo and transmissivity
!    of entire medium as a function of incident beam angle
!    (many simplifications because boundary condition is just
!    isotropic illumination, there are no thermal sources, and
!    particular solutions do not need to be computed).  See
!    Ref. S2 and references therein for details.
          
!    The basic idea is as follows.  The reciprocity principle leads to
!    the following relationships for a plane-parallel, vertically
!    inhomogeneous medium lacking thermal (or other internal) sources:
          
!       albedo(theta) = u_0(theta) for unit-intensity isotropic
!                       illumination at *top* boundary
          
!       trans(theta) =  u_0(theta) for unit-intensity isotropic
!                       illumination at *bottom* boundary
          
!    where
          
!       albedo(theta) = albedo for beam incidence at angle theta
!       trans(theta) = transmissivity for beam incidence at angle theta
!       u_0(theta) = upward azim-avg intensity at top boundary
!                    at angle theta
          
          
!   O U T P U T    V A R I A B L E S:
          
!       ALBMED(IU)   Albedo of the medium as a function of incident
!                    beam angle cosine UMU(IU)
          
!       TRNMED(IU)   Transmissivity of the medium as a function of
!                    incident beam angle cosine UMU(IU)
          
          
!    I N T E R N A L   V A R I A B L E S:
          
!       NCD         number of diagonals below/above main diagonal
          
!       RCOND       estimate of the reciprocal condition of matrix
!                   CBAND; for system  CBAND*X = B, relative
!                   perturbations in CBAND and B of size epsilon may
!                   cause relative perturbations in X of size
!                   epsilon/RCOND.  If RCOND is so small that
!                          1.0 + RCOND .EQ. 1.0
!                   is true, then CBAND may be singular to working
!                   precision.
          
!       CBAND       Left-hand side matrix of linear system Eq. SC(5),
!                   scaled by Eq. SC(12); in banded form required
!                   by LINPACK solution routines
          
!       NCOL        number of columns in CBAND matrix
          
!       IPVT        INTEGER vector of pivot indices
          
!       (most others documented in DISORT)
          
!   Called by- DISORT
!   Calls- LEPOLY, ZEROIT, SGBCO, SOLEIG, TERPEV, SETMTX, SOLVE1,
!          ALTRIN, SPALTR, PRALTR
! +-------------------------------------------------------------------+
          
!     .. Scalar Arguments ..
          
          
          REAL, INTENT(IN OUT)                     :: ALBEDO  AM
            REAL, INTENT(IN OUT)                     :: APB( MI, MI )
            REAL, INTENT(IN OUT)                     :: ARRAY( MXCMU, MXCMU )
            REAL, INTENT(IN OUT)                     :: B( NNLYRI )
            REAL, INTENT(IN OUT)                     :: BDR( MI, 0:MI )
            REAL, INTENT(IN OUT)                     :: CBAND( MI9M2, NNLYRI )
            REAL, INTENT(IN OUT)                     :: CC( MXCMU, MXCMU )
            REAL, INTENT(IN OUT)                     :: CMU( MXCMU )
            REAL, INTENT(IN OUT)                     :: CWT( MXCMU )
            REAL, INTENT(IN OUT)                     :: DTAUCP( * )
            REAL, INTENT(IN OUT)                     :: EVAL( MI )
            REAL, INTENT(IN OUT)                     :: EVECC( MXCMU, MXCMU )
            REAL, INTENT(IN OUT)                     :: GL( 0:MXCMU, * )
            REAL, INTENT(IN OUT)                     :: GC( MXCMU, MXCMU, * )
            REAL, INTENT(IN OUT)                     :: GU( MXUMU, MXCMU, * )
            INTEGER, INTENT(IN OUT)                  :: IPVT( * )
            REAL, INTENT(IN OUT)                     :: KK( MXCMU, * )
            REAL, INTENT(IN OUT)                     :: LL( MXCMU, * )
            INTEGER, INTENT(IN)                      :: NLYR
            INTEGER, INTENT(IN)                      :: NN
            INTEGER, INTENT(IN)                      :: NSTR
            INTEGER, INTENT(IN OUT)                  :: NUMU
            LOGICAL, INTENT(IN)                      :: PRNT( * )
            REAL, INTENT(IN OUT)                     :: TAUCPR( 0:* )
            REAL, INTENT(OUT)                        :: UMU( MAXUMU )
            REAL, INTENT(IN)                         :: U0U( MXUMU, * )
            REAL, INTENT(IN OUT)                     :: WK( MXCMU )
            REAL, INTENT(OUT)                        :: YLMC( 0:MXCMU, MXCMU )
            REAL, INTENT(IN OUT)                     :: YLMU( 0:MXCMU, * )
            REAL, INTENT(IN OUT)                     :: Z( NNLYRI )
            REAL, INTENT(IN OUT)                     :: WKD( MXCMU )
            INTEGER, INTENT(IN OUT)                  :: MI
            INTEGER, INTENT(IN OUT)                  :: MI9M2
            INTEGER, INTENT(IN OUT)                  :: MAXUMU
            INTEGER, INTENT(IN OUT)                  :: MXCMU
            INTEGER, INTENT(IN OUT)                  :: MXUMU
            INTEGER, INTENT(IN OUT)                  :: NNLYRI
            REAL, INTENT(IN OUT)                     :: SQT( * )
            REAL, INTENT(OUT)                        :: ALBMED( MAXUMU )
            REAL, INTENT(OUT)                        :: TRNMED( MAXUMU )
            
            REAL :: ALBEDO
!     ..
!     .. Array Arguments ..
              
              
              
              REAL :: AMB( MI, MI ),   &
                     &
                     &
                     &
                  
              
              
!     ..
!     .. Local Scalars ..
              
              LOGICAL :: LAMBER, LYRCUT
              INTEGER :: IQ, IU, L, LC, MAZIM, NCD, NCOL, NCUT
              REAL :: DELM0, FISOT, RCOND, SGN, SPHALB, SPHTRN
!     ..
!     .. External Subroutines ..
              
              EXTERNAL  ALTRIN, ERRMSG, LEPOLY, PRALTR, SETMTX, SGBCO, SOLEIG,  &
                  SOLVE1, SPALTR, TERPEV, ZEROIT
!     ..
!     .. Intrinsic Functions ..
              
              INTRINSIC EXP
!     ..
              
              MAZIM  = 0
              DELM0  = 1.0
!                    ** Set DISORT variables that are ignored in this
!                    ** special case but are needed below in argument
!                    ** lists of subroutines shared with general case
              NCUT   = NLYR
              LYRCUT = .FALSE.
              FISOT  = 1.0
              LAMBER = .TRUE.
!                          ** Get Legendre polynomials for computational
!                          ** and user polar angle cosines
              
              CALL LEPOLY( NUMU, MAZIM, MXCMU, NSTR-1, UMU, SQT, YLMU )
              
              CALL LEPOLY( NN, MAZIM, MXCMU, NSTR-1, CMU, SQT, YLMC )
              
!                       ** Evaluate Legendre polynomials with negative
!                       ** arguments from those with positive arguments;
!                       ** Dave/Armstrong Eq. (15), STWL(59)
              SGN  = -1.0
              
              DO  L = MAZIM, NSTR - 1
                
                SGN  = -SGN
                
                DO  IQ = NN + 1, NSTR
                  YLMC( L, IQ ) = SGN*YLMC( L, IQ - NN )
                END DO
                
              END DO
!                                  ** Zero out bottom reflectivity
!                                  ** (ALBEDO is used only in analytic
!                                  ** formulae involving ALBEDO = 0
!                                  ** solutions; Eqs 16-17 of Ref S2)
              
              CALL ZEROIT( BDR, MI*( MI+1 ) )
              
              
! ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============
              
              DO  LC = 1, NLYR
                
!                                       ** Solve eigenfunction problem
!                                       ** in Eq. STWJ(8b), STWL(23f)
                
                CALL SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL( 0,LC ), MI, MAZIM,  &
                    MXCMU, NN, NSTR, YLMC, CC, EVECC, EVAL,  &
                    KK( 1,LC ), GC( 1,1,LC ), WKD )
                
!                          ** Interpolate eigenvectors to user angles
                
                CALL TERPEV( CWT, EVECC, GL( 0,LC ), GU( 1,1,LC ), MAZIM,  &
                    MXCMU, MXUMU, NN, NSTR, NUMU, WK, YLMC, YLMU )
                
              END DO
              
! ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============
              
              
!                      ** Set coefficient matrix (CBAND) of equations
!                      ** combining boundary and layer interface
!                      ** conditions (in band-storage mode required by
!                      ** LINPACK routines)
              
              CALL SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,  &
                  LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT,  &
                  NNLYRI, NN, NSTR, TAUCPR, WK )
              
!                      ** LU-decompose the coeff. matrix (LINPACK)
              
              NCD  = 3*NN - 1
              CALL SGBCO( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z )
              IF( 1.0+RCOND == 1.0 )  &
                  CALL ERRMSG('ALBTRN--SGBCO says matrix near singular',.FALSE.)
              
!                             ** First, illuminate from top; if only
!                             ** one layer, this will give us everything
              
!                             ** Solve for constants of integration in
!                             ** homogeneous solution
              
              CALL SOLVE1( B, CBAND, FISOT, 1, IPVT, LL, MI9M2, MXCMU,  &
                  NCOL, NLYR, NN, NNLYRI, NSTR )
              
!                             ** Compute azimuthally-averaged intensity
!                             ** at user angles; gives albedo if multi-
!                             ** layer (Eq. 9 of Ref S2); gives both
!                             ** albedo and transmissivity if single
!                             ** layer (Eqs. 3-4 of Ref S2)
              
              CALL ALTRIN( GU, KK, LL, MXCMU, MXUMU, MAXUMU, NLYR, NN, NSTR,  &
                  NUMU, TAUCPR, UMU, U0U, WK )
              
!                               ** Get beam-incidence albedos from
!                               ** reciprocity principle
              DO  IU = 1, NUMU / 2
                ALBMED( IU ) = U0U( IU + NUMU/2, 1 )
              END DO
              
              
              IF( NLYR == 1 ) THEN
                
                DO  IU = 1, NUMU / 2
!                               ** Get beam-incidence transmissivities
!                               ** from reciprocity principle (1 layer);
!                               ** flip them end over end to correspond
!                               ** to positive UMU instead of negative
                  
                  TRNMED( IU ) = U0U( NUMU/2 + 1 - IU, 2 )  &
                      + EXP( -TAUCPR( NLYR ) / UMU( IU + NUMU/2 ) )
                  
                END DO
                
              ELSE
!                             ** Second, illuminate from bottom
!                             ** (if multiple layers)
                
                CALL SOLVE1( B, CBAND, FISOT, 2, IPVT, LL, MI9M2, MXCMU,  &
                    NCOL, NLYR, NN, NNLYRI, NSTR )
                
                CALL ALTRIN( GU, KK, LL, MXCMU, MXUMU, MAXUMU, NLYR, NN, NSTR,  &
                    NUMU, TAUCPR, UMU, U0U, WK )
                
!                               ** Get beam-incidence transmissivities
!                               ** from reciprocity principle
                DO  IU = 1, NUMU / 2
                  TRNMED( IU ) = U0U( IU + NUMU/2, 1 )  &
                      + EXP( -TAUCPR( NLYR ) / UMU( IU + NUMU/2 ) )
                END DO
                
              END IF
              
              
              IF( ALBEDO > 0.0 ) THEN
                
!                             ** Get spherical albedo and transmissivity
                IF( NLYR == 1 ) THEN
                  
                  CALL SPALTR( CMU, CWT, GC, KK, LL, MXCMU, NLYR,  &
                      NN, NSTR, TAUCPR, SPHALB, SPHTRN )
                ELSE
                  
                  CALL SPALTR( CMU, CWT, GC, KK, LL, MXCMU, NLYR,  &
                      NN, NSTR, TAUCPR, SPHTRN, SPHALB )
                END IF
                
!                                ** Ref. S2, Eqs. 16-17 (these eqs. have
!                                ** a simple physical interpretation
!                                ** like that of adding-doubling eqs.)
                DO  IU = 1, NUMU
                  
                  ALBMED(IU) = ALBMED(IU) + ( ALBEDO / (1.-ALBEDO*SPHALB) )  &
                        * SPHTRN * TRNMED(IU)
                    
                    TRNMED(IU) = TRNMED(IU) + ( ALBEDO / (1.-ALBEDO*SPHALB) )  &
                          * SPHALB * TRNMED(IU)
                    END DO
                    
                  END IF
!                          ** Return UMU to all positive values, to
!                          ** agree with ordering in ALBMED, TRNMED
                  NUMU  = NUMU / 2
                  DO  IU = 1, NUMU
                    UMU( IU ) = UMU( IU + NUMU )
                  END DO
                  
                  IF( PRNT(4) ) CALL PRALTR( UMU, NUMU, ALBMED, TRNMED )
                  
                  
                  RETURN
                END SUBROUTINE ALBTRN
                
!************************************************************************
                
                SUBROUTINE ALTRIN( GU, KK, LL, MXCMU, MXUMU, MAXUMU, NLYR, NN,  &
                    NSTR, NUMU, TAUCPR, UMU, U0U, WK )
                
!       Computes azimuthally-averaged intensity at top and bottom
!       of medium (related to albedo and transmission of medium by
!       reciprocity principles; see Ref S2).  User polar angles are
!       used as incident beam angles. (This is a very specialized
!       version of USRINT)
                
!       ** NOTE **  User input values of UMU (assumed positive) are
!                   temporarily in upper locations of  UMU  and
!                   corresponding negatives are in lower locations
!                   (this makes GU come out right).  I.e. the contents
!                   of the temporary UMU array are:
                
!                     -UMU(NUMU),..., -UMU(1), UMU(1),..., UMU(NUMU)
                
                
!   I N P U T    V A R I A B L E S:
                
!       GU     :  Eigenvectors interpolated to user polar angles
!                   (i.e., g in Eq. SC(1), STWL(31ab))
                
!       KK     :  Eigenvalues of coeff. matrix in Eq. SS(7), STWL(23b)
                
!       LL     :  Constants of integration in Eq. SC(1), obtained
!                   by solving scaled version of Eq. SC(5);
!                   exponential term of Eq. SC(12) not included
                
!       NN     :  Order of double-Gauss quadrature (NSTR/2)
                
!       TAUCPR :  Cumulative optical depth (delta-M-scaled)
                
!       (remainder are DISORT input variables)
                
                
!   O U T P U T    V A R I A B L E:
                
!       U0U  :    Diffuse azimuthally-averaged intensity at top and
!                 bottom of medium (directly transmitted component,
!                 corresponding to BNDINT in USRINT, is omitted).
                
                
!   I N T E R N A L    V A R I A B L E S:
                
!       DTAU   :  Optical depth of a computational layer
!       PALINT :  Non-boundary-forced intensity component
!       UTAUPR :  Optical depths of user output levels (delta-M scaled)
!       WK     :  Scratch vector for saving 'EXP' evaluations
!       All the exponential factors (i.e., EXP1, EXPN,... etc.)
!       come from the substitution of constants of integration in
!       Eq. SC(12) into Eqs. S1(8-9).  All have negative arguments.
                
!   Called by- ALBTRN
! +-------------------------------------------------------------------+
                
!     .. Scalar Arguments ..
                
                
                REAL, INTENT(IN)                         :: GU( MXUMU, MXCMU, * )
                REAL, INTENT(IN)                         :: KK( MXCMU, * )
                REAL, INTENT(IN)                         :: LL( MXCMU, * )
                INTEGER, INTENT(IN OUT)                  :: MXCMU
                INTEGER, INTENT(IN OUT)                  :: MXUMU
                INTEGER, INTENT(IN OUT)                  :: MAXUMU
                INTEGER, INTENT(IN)                      :: NLYR
                INTEGER, INTENT(IN)                      :: NN
                INTEGER, INTENT(IN)                      :: NSTR
                INTEGER, INTENT(IN)                      :: NUMU
                REAL, INTENT(IN)                         :: TAUCPR( 0:* )
                REAL, INTENT(IN)                         :: UMU( MAXUMU )
                REAL, INTENT(OUT)                        :: U0U( MXUMU, * )
                REAL, INTENT(OUT)                        :: WK( MXCMU )
                
!     ..
!     .. Array Arguments ..
                
                
!     ..
!     .. Local Scalars ..
                
                INTEGER :: IQ, IU, IUMAX, IUMIN, LC, LU
                REAL :: DENOM, DTAU, EXP1, EXP2, EXPN, MU, PALINT, SGN
!     ..
!     .. Local Arrays ..
                
                REAL :: UTAUPR( 2 )
!     ..
!     .. Intrinsic Functions ..
                
                INTRINSIC ABS, EXP
!     ..
                
                
                UTAUPR( 1 ) = 0.0
                UTAUPR( 2 ) = TAUCPR( NLYR )
                
                DO  LU = 1, 2
                  
                  IF( LU == 1 ) THEN
                    
                    IUMIN  = NUMU / 2 + 1
                    IUMAX  = NUMU
                    SGN    = 1.0
                    
                  ELSE
                    
                    IUMIN  = 1
                    IUMAX  = NUMU / 2
                    SGN    = - 1.0
                    
                  END IF
!                                   ** Loop over polar angles at which
!                                   ** albedos/transmissivities desired
!                                   ** ( upward angles at top boundary,
!                                   ** downward angles at bottom )
                  DO  IU = IUMIN, IUMAX
                    
                    MU   = UMU( IU )
!                                     ** Integrate from top to bottom
!                                     ** computational layer
                    PALINT = 0.0
                    
                    DO  LC = 1, NLYR
                      
                      DTAU   = TAUCPR( LC ) - TAUCPR( LC - 1 )
                      EXP1   = EXP( ( UTAUPR( LU ) - TAUCPR( LC - 1 ) ) / MU )
                      EXP2   = EXP( ( UTAUPR( LU ) - TAUCPR( LC ) ) / MU )
                      
!                                      ** KK is negative
                      DO  IQ = 1, NN
                        
                        WK( IQ ) = EXP( KK( IQ,LC )*DTAU )
                        DENOM  = 1.0 + MU*KK( IQ, LC )
                        
                        IF( ABS( DENOM ) < 0.0001 ) THEN
!                                                   ** L'Hospital limit
                          EXPN   = DTAU / MU*EXP2
                          
                        ELSE
                          
                          EXPN   = ( EXP1*WK( IQ ) - EXP2 )*SGN / DENOM
                          
                        END IF
                        
                        PALINT = PALINT + GU( IU, IQ, LC )*LL( IQ, LC )*EXPN
                        
                      END DO
                      
!                                        ** KK is positive
                      DO  IQ = NN + 1, NSTR
                        
                        DENOM  = 1.0 + MU*KK( IQ, LC )
                        
                        IF( ABS( DENOM ) < 0.0001 ) THEN
                          
                          EXPN   = - DTAU / MU * EXP1
                          
                        ELSE
                          
                          EXPN = ( EXP1 - EXP2 * WK(NSTR+1-IQ) ) *SGN / DENOM
                          
                        END IF
                        
                        PALINT = PALINT + GU( IU, IQ, LC )*LL( IQ, LC )*EXPN
                        
                      END DO
                      
                    END DO
                    
                    U0U( IU, LU ) = PALINT
                    
                  END DO
                  
                END DO
                
                
                RETURN
              END SUBROUTINE ALTRIN
              
!************************************************************************
              
              SUBROUTINE PRALTR( UMU, NUMU, ALBMED, TRNMED )
              
              
              REAL, INTENT(IN OUT)                     :: UMU( NUMU )
              INTEGER, INTENT(IN)                      :: NUMU
              REAL, INTENT(IN OUT)                     :: ALBMED( NUMU )
              REAL, INTENT(IN OUT)                     :: TRNMED( NUMU )
              IMPLICIT NONE
              INCLUDE '../INCLUDE/scatterparam.f90'
              
!        Print planar albedo and transmissivity of medium
!        as a function of incident beam angle
              
!   Called by- ALBTRN
! --------------------------------------------------------------------
              
!     .. Parameters ..
              
              
              REAL, PARAMETER :: DPR = 180.0 / 3.14159265
!     ..
!     .. Scalar Arguments ..
              
              
!     ..
!     .. Array Arguments ..
              
              
!     ..
!     .. Local Scalars ..
              
              INTEGER :: IU
!     ..
!     .. Intrinsic Functions ..
              
              INTRINSIC ACOS
!     ..
              
              
              WRITE( KSTDWARN, '(///,A,//,A)' )  &
                  ' *******  Flux Albedo and/or Transmissivity of ' //  &
                  'entire medium  ********',  &
                  ' Beam Zen Ang   cos(Beam Zen Ang)      Albedo   Transmissivity'
              
              DO  IU = 1, NUMU
                WRITE( KSTDWARN, '(0P,F13.4,F20.6,F12.5,1P,E17.4)' )  &
                    DPR*ACOS( UMU( IU ) ), UMU( IU ), ALBMED( IU ), TRNMED( IU )
              END DO
              
              
              RETURN
            END SUBROUTINE PRALTR
            
!************************************************************************
            
            SUBROUTINE SOLVE1( B, CBAND, FISOT, IHOM, IPVT, LL, MI9M2, MXCMU,  &
                NCOL, NCUT, NN, NNLYRI, NSTR )
            
!        Construct right-hand side vector B for isotropic incidence
!        (only) on either top or bottom boundary and solve system
!        of equations obtained from the boundary conditions and the
!        continuity-of-intensity-at-layer-interface equations
            
            
!     I N P U T      V A R I A B L E S:
            
!       CBAND    :  Left-hand side matrix of banded linear system
!                   Eq. SC(5), scaled by Eq. SC(12); assumed already
!                   in LU-decomposed form, ready for LINPACK solver
            
!       IHOM     :  Direction of illumination flag (1, top; 2, bottom)
            
!       NCOL     :  Number of columns in CBAND
            
!       NN       :  Order of double-Gauss quadrature (NSTR/2)
            
!       (remainder are DISORT input variables)
            
            
!    O U T P U T     V A R I A B L E S:
            
!       B        :  Right-hand side vector of Eq. SC(5) going into
!                   SGBSL; returns as solution vector of Eq.
!                   SC(12), constants of integration without
!                   exponential term
            
!       LL      :   permanent storage for B, but re-ordered
            
            
!    I N T E R N A L    V A R I A B L E S:
            
!       IPVT     :  INTEGER vector of pivot indices
!       NCD      :  Number of diagonals below or above main diagonal
            
!   Called by- ALBTRN
!   Calls- ZEROIT, SGBSL
! +-------------------------------------------------------------------+
            
!     .. Scalar Arguments ..
            
            
            REAL, INTENT(OUT)                        :: B( NNLYRI )
            REAL, INTENT(IN OUT)                     :: CBAND( MI9M2, NNLYRI )
            REAL, INTENT(IN)                         :: FISOT
            INTEGER, INTENT(IN)                      :: IHOM
            INTEGER, INTENT(IN OUT)                  :: IPVT( NNLYRI )
            REAL, INTENT(OUT)                        :: LL( MXCMU, * )
            INTEGER, INTENT(IN OUT)                  :: MI9M2
            INTEGER, INTENT(IN OUT)                  :: MXCMU
            INTEGER, INTENT(OUT)                     :: NCOL
            INTEGER, INTENT(IN)                      :: NCUT
            INTEGER, INTENT(IN OUT)                  :: NN
            INTEGER, INTENT(IN OUT)                  :: NNLYRI
            INTEGER, INTENT(IN)                      :: NSTR
            
            
!     ..
!     .. Array Arguments ..
            
            
            
!     ..
!     .. Local Scalars ..
            
            INTEGER :: I, IPNT, IQ, LC, NCD
!     ..
!     .. External Subroutines ..
            
            EXTERNAL  SGBSL, ZEROIT
!     ..
            
            
            CALL ZEROIT( B, NNLYRI )
            
            IF( IHOM == 1 ) THEN
!                             ** Because there are no beam or emission
!                             ** sources, remainder of B array is zero
              DO  I = 1, NN
                B( I )             = FISOT
                B( NCOL - NN + I ) = 0.0
              END DO
              
            ELSE IF( IHOM == 2 ) THEN
              
              DO  I = 1, NN
                B( I )             = 0.0
                B( NCOL - NN + I ) = FISOT
              END DO
              
            END IF
            
            
            NCD  = 3*NN - 1
            CALL SGBSL( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, B, 0 )
            
            DO  LC = 1, NCUT
              
              IPNT  = LC*NSTR - NN
              
              DO  IQ = 1, NN
                LL( NN + 1 - IQ, LC ) = B( IPNT + 1 - IQ )
                LL( IQ + NN,     LC ) = B( IQ + IPNT )
              END DO
              
            END DO
            
            
            RETURN
          END SUBROUTINE SOLVE1
          
!************************************************************************
          
          SUBROUTINE SPALTR( CMU, CWT, GC, KK, LL, MXCMU, NLYR, NN, NSTR,  &
              TAUCPR, SFLUP, SFLDN )
          
!       Calculates spherical albedo and transmissivity for the entire
!       medium from the m=0 intensity components
!       (this is a very specialized version of FLUXES)
          
          
!    I N P U T    V A R I A B L E S:
          
!       CMU,CWT    Abscissae, weights for Gauss quadrature
!                  over angle cosine
          
!       KK      :  Eigenvalues of coeff. matrix in eq. SS(7)
          
!       GC      :  Eigenvectors at polar quadrature angles, SC(1)
          
!       LL      :  Constants of integration in eq. SC(1), obtained
!                  by solving scaled version of Eq. SC(5);
!                  exponential term of Eq. SC(12) not included
          
!       NN      :  Order of double-Gauss quadrature (NSTR/2)
          
!       (remainder are DISORT input variables)
          
          
!    O U T P U T   V A R I A B L E S:
          
!       SFLUP   :  Up-flux at top (equivalent to spherical albedo due to
!                  reciprocity).  For illumination from below it gives
!                  spherical transmissivity
          
!       SFLDN   :  Down-flux at bottom (for single layer, equivalent to
!                  spherical transmissivity due to reciprocity)
          
          
!    I N T E R N A L   V A R I A B L E S:
          
!       ZINT    :  Intensity of m=0 case, in Eq. SC(1)
          
!   Called by- ALBTRN
! +--------------------------------------------------------------------
          
!     .. Scalar Arguments ..
          
          
          REAL, INTENT(IN)                         :: CMU( MXCMU )
          REAL, INTENT(IN)                         :: CWT( MXCMU )
          REAL, INTENT(IN)                         :: GC( MXCMU, MXCMU, * )
          REAL, INTENT(IN OUT)                     :: KK( MXCMU, * )
          REAL, INTENT(IN)                         :: LL( MXCMU, * )
          INTEGER, INTENT(IN OUT)                  :: MXCMU
          INTEGER, INTENT(IN)                      :: NLYR
          INTEGER, INTENT(IN)                      :: NN
          INTEGER, INTENT(IN)                      :: NSTR
          REAL, INTENT(IN OUT)                     :: TAUCPR( 0:* )
          REAL, INTENT(OUT)                        :: SFLUP
          REAL, INTENT(OUT)                        :: SFLDN
          
          
!     ..
!     .. Array Arguments ..
          
          
!     ..
!     .. Local Scalars ..
          
          INTEGER :: IQ, JQ
          REAL :: ZINT
!     ..
!     .. Intrinsic Functions ..
          
          INTRINSIC EXP
!     ..
          
          
          SFLUP  = 0.0
          
          DO  IQ = NN + 1, NSTR
            
            ZINT   = 0.0
            DO  JQ = 1, NN
              ZINT  = ZINT + GC( IQ, JQ, 1 )*LL( JQ, 1 )*  &
                  EXP( KK( JQ,1 )*TAUCPR( 1 ) )
            END DO
            
            DO  JQ = NN + 1, NSTR
              ZINT  = ZINT + GC( IQ, JQ, 1 )*LL( JQ, 1 )
            END DO
            
            SFLUP  = SFLUP + CWT( IQ - NN )*CMU( IQ - NN )*ZINT
            
          END DO
          
          
          SFLDN  = 0.0
          
          DO  IQ = 1, NN
            
            ZINT   = 0.0
            DO  JQ = 1, NN
              ZINT  = ZINT + GC( IQ, JQ, NLYR )*LL( JQ, NLYR )
            END DO
            
            DO  JQ = NN + 1, NSTR
              ZINT  = ZINT + GC( IQ, JQ, NLYR )*LL( JQ, NLYR )*  &
                  EXP( - KK( JQ,NLYR ) * ( TAUCPR( NLYR ) - TAUCPR( NLYR-1 ) ) )
            END DO
            
            SFLDN  = SFLDN + CWT( NN + 1 - IQ )*CMU( NN + 1 - IQ )*ZINT
            
          END DO
          
          SFLUP  = 2.0*SFLUP
          SFLDN  = 2.0*SFLDN
          
          
          RETURN
        END SUBROUTINE SPALTR
        
! ******************************************************************
! ********** End of IBCND=1 special case routines ******************
! ******************************************************************
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! RCS version control information:
! $Header: ErrPack.f,v 2.1 2000/03/27 21:40:49 laszlo Exp $
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        SUBROUTINE  ErrMsg( MESSAG, FATAL )
        
        
        CHARACTER (LEN=*), INTENT(IN OUT)        :: MESSAG
        LOGICAL, INTENT(IN)                      :: FATAL
        IMPLICIT NONE
        INCLUDE '../INCLUDE/scatterparam.f90'
        
!        Print out a warning or error message;  abort if error
        
        LOGICAL :: MsgLim
        
        INTEGER :: MaxMsg, NumMsg
        SAVE          MaxMsg, NumMsg, MsgLim
        DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /
        
        
        IF ( FATAL )  THEN
          WRITE ( KSTDERR, '(/,2A,/)' )  ' ******* ERROR >>>>>>  ', MESSAG
          STOP
        END IF
        
        NumMsg = NumMsg + 1
        IF( MsgLim )  RETURN
        
        IF ( NumMsg <= MaxMsg )  THEN
          WRITE ( KSTDWARN, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
        ELSE
          WRITE ( KSTDWARN,99 )
          MsgLim = .True.
        END IF
        
        RETURN
        
        99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',  &
            'They will no longer be printed  <<<<<<<', // )
      END SUBROUTINE  ErrMsg
      
!************************************************************************
      
      LOGICAL FUNCTION  WrtBad ( VarNam )
      
      
      CHARACTER (LEN=*), INTENT(IN OUT)        :: VarNam
      IMPLICIT NONE
      INCLUDE '../INCLUDE/scatterparam.f90'
      
!          Write names of erroneous variables and return 'TRUE'
      
!      INPUT :   VarNam = Name of erroneous variable to be written
!                         ( CHARACTER, any length )
      
      
      INTEGER :: MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 50 /
      
      
      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( KSTDWARN, '(3A)' )  ' ****  Input variable  ', VarNam,  &
          '  in error  ****'
      IF ( NumMsg == MaxMsg )  &
          CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )
      
      RETURN
    END FUNCTION  WrtBad
    
!************************************************************************
    
    LOGICAL FUNCTION  WrtDim ( DimNam, MinVal )
    
    
    CHARACTER (LEN=*), INTENT(IN OUT)        :: DimNam
    INTEGER, INTENT(IN OUT)                  :: MinVal
    IMPLICIT NONE
    INCLUDE '../INCLUDE/scatterparam.f90'
    
!          Write name of too-small symbolic dimension and
!          the value it should be increased to;  return 'TRUE'
    
!      INPUT :  DimNam = Name of symbolic dimension which is too small
!                        ( CHARACTER, any length )
!               Minval = Value to which that dimension should be
!                        increased (at least)
    
    
    
    
    
    WRITE ( KSTDWARN, '(/,3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,  &
        '  should be increased to at least ', MinVal
    
    WrtDim = .TRUE.
    
    RETURN
  END FUNCTION  WrtDim
  
!************************************************************************
  
  LOGICAL FUNCTION  TstBad( VarNam, RelErr )
  
  
  CHARACTER (LEN=*), INTENT(IN OUT)        :: VarNam
  REAL, INTENT(IN OUT)                     :: RelErr
  IMPLICIT NONE
  INCLUDE '../INCLUDE/scatterparam.f90'
  
!       Write name (VarNam) of variable failing self-test and its
!       percent error from the correct value;  return  'FALSE'.
  
  
  
  
  
  TstBad = .FALSE.
  WRITE( KSTDWARN, '(/,3A,1P,E11.2,A)' )  &
      ' Output variable ', VarNam,' differed by ', 100.*RelErr,  &
      ' per cent from correct value.  Self-test failed.'
  
  RETURN
END FUNCTION  TstBad

!************************************************************************
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! RCS version control information:
! $Header: D1MACH.f,v 1.2 97/03/18 17:05:25 wiscombe Exp $
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DOUBLE PRECISION FUNCTION D1MACH(I)

!  Double-precision machine constants (see R1MACH for documentation).

!  By default, returns values appropriate for a computer with IEEE
!  arithmetic.  This is an abbreviated version of a routine widely
!  used for 20+ years by numerical analysts.  Most of the values in
!  the original version pertain to computers which went to computer
!  heaven years ago and are of little if any interest.

!  If the values herein do not work for any reason, just look in
!  your Fortran manual for the correct values (usually in the part
!  discussing representations of numbers) and insert them. The exact
!  values are not that important; they can be a factor of 2-3 off
!  without causing any harm.

!  Only I = 1,2,4 is actually used by DISORT.

!  This routine is superseded in Fortran-90 by the intrinsic numeric
!  inquiry functions HUGE(1.D0), TINY(1.D0), and EPSILON(1.D0).

!  The original version can be found on NetLib (search by name):
!      http://www.netlib.org/
! ====================================================================


INTEGER, INTENT(IN)                      :: I

EXTERNAL  ERRMSG

IF( I == 1 )  THEN
  D1MACH = 2.3D-308
!        D1MACH = TINY(1.D0)
ELSE IF( I == 2 )  THEN
  D1MACH = 1.7D+308
!        D1MACH = HUGE(1.D0)
ELSE IF( I == 4 )  THEN
  D1MACH = 2.3D-16
!        D1MACH = EPSILON(1.D0)
ELSE
  CALL ERRMSG( 'D1MACH--argument incorrect', .TRUE.)
END IF

RETURN
END FUNCTION D1MACH

!************************************************************************

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! RCS version control information:
! $Header: LINPAK.f,v 2.1 2000/03/27 21:40:49 laszlo Exp $
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Call tree:

!    SGBCO
!       SASUM_DIS
!       SDOT_DIS
!       SAXPY_DIS
!       SGBFA
!           ISAMAX_DIS
!           SAXPY_DIS
!           SSCAL_DIS
!       SSCAL_DIS
!   SGBSL
!       SDOT_DIS
!       SAXPY_DIS
!   SGECO
!       SASUM_DIS
!       SDOT_DIS
!       SAXPY_DIS
!       SGEFA
!           ISAMAX_DIS
!           SAXPY_DIS
!           SSCAL_DIS
!       SSCAL_DIS
!   SGESL
!       SDOT_DIS
!       SAXPY_DIS
!   SSWAP_DIS
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


SUBROUTINE SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )

!         Factors a real band matrix by Gaussian elimination
!         and estimates the condition of the matrix.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!     If  RCOND  is not needed, SGBFA is slightly faster.
!     To solve  A*X = B , follow SBGCO by SGBSL.

!     input:

!        ABD     REAL(LDA, N)
!                contains the matrix in band storage.  The columns
!                of the matrix are stored in the columns of  ABD  and
!                the diagonals of the matrix are stored in rows
!                ML+1 through 2*ML+MU+1 of  ABD .
!                See the comments below for details.

!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be .GE. 2*ML + MU + 1 .

!        N       INTEGER
!                the order of the original matrix.

!        ML      INTEGER
!                number of diagonals below the main diagonal.
!                0 .LE. ML .LT. N .

!        MU      INTEGER
!                number of diagonals above the main diagonal.
!                0 .LE. MU .LT. N .
!                more efficient if  ML .LE. MU .

!     on return

!        ABD     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.

!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.

!        RCOND   REAL
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  epsilon  may cause
!                relative perturbations in  X  of size  epsilon/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND .EQ. 1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.

!        Z       REAL(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                norm(a*z) = rcond*norm(a)*norm(z) .

!     Band storage

!           If  A  is a band matrix, the following program segment
!           will set up the input.

!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-MU)
!                      I2 = MIN(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE

!           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
!           In addition, the first  ML  rows in  ABD  are used for
!           elements generated during the triangularization.
!           The total number of rows needed in  ABD  is  2*ML+MU+1 .
!           The  ML+MU by ML+MU  upper left triangle and the
!           ML by ML  lower right triangle are not referenced.

!     Example:  if the original matrix is

!           11 12 13  0  0  0
!           21 22 23 24  0  0
!            0 32 33 34 35  0
!            0  0 43 44 45 46
!            0  0  0 54 55 56
!            0  0  0  0 65 66

!      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABD should contain

!            *  *  *  +  +  +  , * = not used
!            *  * 13 24 35 46  , + = used for pivoting
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!           21 32 43 54 65  *

! --------------------------------------------------------------------


!     .. Scalar Arguments ..


REAL, INTENT(IN)                         :: ABD( LDA, * )
INTEGER, INTENT(IN OUT)                  :: LDA
INTEGER, INTENT(IN)                      :: N
INTEGER, INTENT(IN)                      :: ML
INTEGER, INTENT(IN)                      :: MU
INTEGER, INTENT(IN)                      :: IPVT( * )
REAL, INTENT(OUT)                        :: RCOND
REAL, INTENT(OUT)                        :: Z( * )


!     ..
!     .. Array Arguments ..



!     ..
!     .. Local Scalars ..

INTEGER :: INFO, IS, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM
REAL :: ANORM, EK, S, SM, T, WK, WKM, YNORM
!     ..
!     .. External Functions ..

REAL :: SASUM_DIS, SDOT_DIS
EXTERNAL  SASUM_DIS, SDOT_DIS
!     ..
!     .. External Subroutines ..

EXTERNAL  SAXPY_DIS, SGBFA, SSCAL_DIS
!     ..
!     .. Intrinsic Functions ..

INTRINSIC ABS, MAX, MIN, SIGN
!     ..


!                       ** compute 1-norm of A
ANORM  = 0.0E0
L  = ML + 1
IS = L + MU

DO  J = 1, N
  
  ANORM  = MAX( ANORM, SASUM_DIS( L,ABD( IS,J ),1 ) )
  
  IF( IS > ML + 1 ) IS = IS - 1
  
  IF( J <= MU ) L  = L + 1
  
  IF( J >= N - ML ) L  = L - 1
  
END DO
!                                               ** factor

CALL SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

!     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
!     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E.
!     trans(A) is the transpose of A.  The components of E  are
!     chosen to cause maximum local growth in the elements of W  where
!     trans(U)*W = E.  The vectors are frequently rescaled to avoid
!     overflow.

!                     ** solve trans(U)*W = E
EK = 1.0E0

DO  J = 1, N
  Z( J ) = 0.0E0
END DO


M  = ML + MU + 1
JU = 0

DO  K = 1, N
  
  IF( Z( K ) /= 0.0E0 ) EK = SIGN( EK, -Z( K ) )
  
  IF( ABS( EK - Z( K ) ) > ABS( ABD( M,K ) ) ) THEN
    
    S  = ABS( ABD( M,K ) ) / ABS( EK - Z( K ) )
    
    CALL SSCAL_DIS( N, S, Z, 1 )
    
    EK = S*EK
    
  END IF
  
  WK   = EK - Z( K )
  WKM  = -EK - Z( K )
  S    = ABS( WK )
  SM   = ABS( WKM )
  
  IF( ABD( M,K ) /= 0.0E0 ) THEN
    
    WK   = WK / ABD( M, K )
    WKM  = WKM / ABD( M, K )
    
  ELSE
    
    WK   = 1.0E0
    WKM  = 1.0E0
    
  END IF
  
  KP1  = K + 1
  JU   = MIN( MAX( JU,MU + IPVT( K ) ), N )
  MM   = M
  
  IF( KP1 <= JU ) THEN
    
    DO  J = KP1, JU
      MM     = MM - 1
      SM     = SM + ABS( Z( J ) + WKM*ABD( MM,J ) )
      Z( J ) = Z( J ) + WK*ABD( MM, J )
      S      = S + ABS( Z( J ) )
    END DO
    
    IF( S < SM ) THEN
      
      T  = WKM - WK
      WK = WKM
      MM = M
      
      DO  J = KP1, JU
        MM = MM - 1
        Z( J ) = Z( J ) + T*ABD( MM, J )
      END DO
      
    END IF
    
  END IF
  
  Z( K ) = WK
  
END DO


S  = 1.0E0 / SASUM_DIS( N, Z, 1 )

CALL SSCAL_DIS( N, S, Z, 1 )

!                         ** solve trans(L)*Y = W
DO  KB = 1, N
  K  = N + 1 - KB
  LM = MIN( ML, N - K )
  
  IF( K < N ) Z( K ) = Z( K ) + SDOT_DIS( LM, ABD( M+1, K ), 1, Z( K+1 ), 1 )
  
  IF( ABS( Z( K ) ) > 1.0E0 ) THEN
    
    S  = 1.0E0 / ABS( Z( K ) )
    
    CALL SSCAL_DIS( N, S, Z, 1 )
    
  END IF
  
  L      = IPVT( K )
  T      = Z( L )
  Z( L ) = Z( K )
  Z( K ) = T
  
END DO


S  = 1.0E0 / SASUM_DIS( N, Z, 1 )

CALL SSCAL_DIS( N, S, Z, 1 )

YNORM  = 1.0E0
!                         ** solve L*V = Y
DO  K = 1, N
  
  L      = IPVT( K )
  T      = Z( L )
  Z( L ) = Z( K )
  Z( K ) = T
  LM     = MIN( ML, N - K )
  
  IF( K < N ) CALL SAXPY_DIS( LM, T, ABD( M+1, K ), 1, Z( K+1 ), 1 )
  
  IF( ABS( Z(K) ) > 1.0E0 ) THEN
    
    S  = 1.0E0 / ABS( Z(K) )
    
    CALL SSCAL_DIS( N, S, Z, 1 )
    
    YNORM  = S*YNORM
    
  END IF
  
END DO


S  = 1.0E0 / SASUM_DIS( N, Z, 1 )

CALL SSCAL_DIS( N, S, Z, 1 )

YNORM  = S*YNORM

!                           ** solve  U*Z = W
DO  KB = 1, N
  
  K  = N + 1 - KB
  
  IF( ABS( Z( K ) ) > ABS( ABD( M,K ) ) ) THEN
    
    S  = ABS( ABD( M,K ) ) / ABS( Z( K ) )
    
    CALL SSCAL_DIS( N, S, Z, 1 )
    
    YNORM  = S*YNORM
    
  END IF
  
  IF( ABD( M,K ) /= 0.0E0 ) Z( K ) = Z( K ) / ABD( M, K )
  IF( ABD( M,K ) == 0.0E0 ) Z( K ) = 1.0E0
  
  LM = MIN( K, M ) - 1
  LA = M - LM
  LZ = K - LM
  T  = -Z( K )
  
  CALL SAXPY_DIS( LM, T, ABD( LA,K ), 1, Z( LZ ), 1 )
  
END DO
!                              ** make znorm = 1.0

S  = 1.0E0 / SASUM_DIS( N, Z, 1 )

CALL SSCAL_DIS( N, S, Z, 1 )

YNORM  = S*YNORM
IF( ANORM /= 0.0E0 ) RCOND  = YNORM / ANORM
IF( ANORM == 0.0E0 ) RCOND  = 0.0E0

END SUBROUTINE SGBCO

SUBROUTINE SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

!         Factors a real band matrix by elimination.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!     SGBFA is usually called by SBGCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.

!     Input:  same as SGBCO

!     On return:

!        ABD,IPVT    same as SGBCO

!        INFO    INTEGER
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that SGBSL will divide by zero if
!                     called.  Use  RCOND  in SBGCO for a reliable
!                     indication of singularity.

!     (see SGBCO for description of band storage mode)

! ----------------------------------------------------------------


!     .. Scalar Arguments ..


REAL, INTENT(OUT)                        :: ABD( LDA, * )
INTEGER, INTENT(IN OUT)                  :: LDA
INTEGER, INTENT(IN OUT)                  :: N
INTEGER, INTENT(IN)                      :: ML
INTEGER, INTENT(IN)                      :: MU
INTEGER, INTENT(OUT)                     :: IPVT( * )
INTEGER, INTENT(OUT)                     :: INFO

!     ..
!     .. Array Arguments ..



!     ..
!     .. Local Scalars ..

INTEGER :: I, I0, J, J0, J1, JU, JZ, K, KP1, L, LM, M, MM, NM1
REAL :: T
!     ..
!     .. External Functions ..

INTEGER :: ISAMAX_DIS
EXTERNAL  ISAMAX_DIS
!     ..
!     .. External Subroutines ..

EXTERNAL  SAXPY_DIS, SSCAL_DIS
!     ..
!     .. Intrinsic Functions ..

INTRINSIC MAX, MIN
!     ..


M    = ML + MU + 1
INFO = 0
!                        ** zero initial fill-in columns
J0 = MU + 2
J1 = MIN( N, M ) - 1

DO  JZ = J0, J1
  
  I0 = M + 1 - JZ
  
  DO  I = I0, ML
    ABD( I, JZ ) = 0.0E0
  END DO
  
END DO

JZ = J1
JU = 0
!                       ** Gaussian elimination with partial pivoting
NM1  = N - 1

DO  K = 1, NM1
  
  KP1 = K + 1
!                                  ** zero next fill-in column
  JZ = JZ + 1
  
  IF( JZ <= N ) THEN
    
    DO  I = 1, ML
      ABD( I, JZ ) = 0.0E0
    END DO
    
  END IF
!                                  ** find L = pivot index
  LM  = MIN( ML, N - K )
  L   = ISAMAX_DIS( LM + 1, ABD( M, K ), 1 ) + M - 1
  IPVT( K ) = L + K - M
  
  IF( ABD( L,K ) == 0.0E0 ) THEN
!                                      ** zero pivot implies this column
!                                      ** already triangularized
    INFO = K
    
  ELSE
!                                ** interchange if necessary
    IF( L /= M ) THEN
      
      T           = ABD( L, K )
      ABD( L, K ) = ABD( M, K )
      ABD( M, K ) = T
    END IF
!                                      ** compute multipliers
    T  = - 1.0E0 / ABD( M, K )
    
    CALL SSCAL_DIS( LM, T, ABD( M + 1,K ), 1 )
    
!                               ** row elimination with column indexing
    
    JU = MIN( MAX( JU,MU + IPVT( K ) ), N )
    MM = M
    
    DO  J = KP1, JU
      
      L  = L - 1
      MM = MM - 1
      T  = ABD( L, J )
      
      IF( L /= MM ) THEN
        
        ABD( L, J ) = ABD( MM, J )
        ABD( MM, J ) = T
        
      END IF
      
      CALL SAXPY_DIS( LM, T, ABD( M+1, K ), 1, ABD( MM+1, J ), 1)
      
    END DO
    
  END IF
  
END DO


IPVT( N ) = N
IF( ABD( M,N ) == 0.0E0 ) INFO = N

END SUBROUTINE SGBFA

SUBROUTINE SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )

!         Solves the real band system
!            A * X = B  or  transpose(A) * X = B
!         using the factors computed by SBGCO or SGBFA.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!     Input:

!        ABD     REAL(LDA, N)
!                the output from SBGCO or SGBFA.

!        LDA     INTEGER
!                the leading dimension of the array  ABD .

!        N       INTEGER
!                the order of the original matrix.

!        ML      INTEGER
!                number of diagonals below the main diagonal.

!        MU      INTEGER
!                number of diagonals above the main diagonal.

!        IPVT    INTEGER(N)
!                the pivot vector from SBGCO or SGBFA.

!        B       REAL(N)
!                the right hand side vector.

!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  transpose(A)*X = B

!     On return

!        B       the solution vector  X

!     Error condition

!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically, this indicates singularity,
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if SBGCO has set RCOND .GT. 0.0
!        or SGBFA has set INFO .EQ. 0 .

!     To compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!        10 continue

! --------------------------------------------------------

!     .. Scalar Arguments ..


REAL, INTENT(IN)                         :: ABD( LDA, * )
INTEGER, INTENT(IN OUT)                  :: LDA
INTEGER, INTENT(IN)                      :: N
INTEGER, INTENT(IN)                      :: ML
INTEGER, INTENT(IN)                      :: MU
INTEGER, INTENT(IN)                      :: IPVT( * )
REAL, INTENT(IN OUT)                     :: B( * )
INTEGER, INTENT(IN)                      :: JOB

!     ..
!     .. Array Arguments ..



!     ..
!     .. Local Scalars ..

INTEGER :: K, KB, L, LA, LB, LM, M, NM1
REAL :: T
!     ..
!     .. External Functions ..

REAL :: SDOT_DIS
EXTERNAL  SDOT_DIS
!     ..
!     .. External Subroutines ..

EXTERNAL  SAXPY_DIS
!     ..
!     .. Intrinsic Functions ..

INTRINSIC MIN
!     ..


M   = MU + ML + 1
NM1 = N - 1

IF( JOB == 0 ) THEN
!                           ** solve  A * X = B
  
!                               ** first solve L*Y = B
  IF( ML /= 0 ) THEN
    
    DO  K = 1, NM1
      
      LM = MIN( ML, N - K )
      L  = IPVT( K )
      T  = B( L )
      
      IF( L /= K ) THEN
        
        B( L ) = B( K )
        B( K ) = T
        
      END IF
      
      CALL SAXPY_DIS( LM, T, ABD( M + 1,K ), 1, B( K + 1 ), 1 )
      
    END DO
    
  END IF
  
!                           ** now solve  U*X = Y
  DO  KB = 1, N
    
    K      = N + 1 - KB
    B( K ) = B( K ) / ABD( M, K )
    LM     = MIN( K, M ) - 1
    LA     = M - LM
    LB     = K - LM
    T      = -B( K )
    
    CALL SAXPY_DIS( LM, T, ABD( LA,K ), 1, B( LB ), 1 )
    
  END DO
  
  
ELSE
!                          ** solve  trans(A) * X = B
  
!                                  ** first solve  trans(U)*Y = B
  DO  K = 1, N
    
    LM     = MIN( K, M ) - 1
    LA     = M - LM
    LB     = K - LM
    T      = SDOT_DIS( LM, ABD( LA,K ), 1, B( LB ), 1 )
    B( K ) = ( B( K ) - T ) / ABD( M, K )
    
  END DO
  
!                                  ** now solve trans(L)*X = Y
  IF( ML /= 0 ) THEN
    
    DO  KB = 1, NM1
      
      K      = N - KB
      LM     = MIN( ML, N - K )
      B( K ) = B( K ) + SDOT_DIS( LM, ABD( M+1, K ), 1, B( K+1 ), 1 )
      L      = IPVT( K )
      
      IF( L /= K ) THEN
        
        T    = B( L )
        B( L ) = B( K )
        B( K ) = T
        
      END IF
      
    END DO
    
  END IF
  
END IF

END SUBROUTINE SGBSL

SUBROUTINE SGECO( A, LDA, N, IPVT, RCOND, Z )

!         Factors a real matrix by Gaussian elimination
!         and estimates the condition of the matrix.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!         If  RCOND  is not needed, SGEFA is slightly faster.
!         To solve  A*X = B , follow SGECO by SGESL.

!     On entry

!        A       REAL(LDA, N)
!                the matrix to be factored.

!        LDA     INTEGER
!                the leading dimension of the array  A .

!        N       INTEGER
!                the order of the matrix  A .

!     On return

!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U , where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.

!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.

!        RCOND   REAL
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  epsilon  may cause
!                relative perturbations in  X  of size  epsilon/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND .EQ. 1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.

!        Z       REAL(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                norm(A*Z) = RCOND*norm(A)*norm(Z) .

! ------------------------------------------------------------------

!     .. Scalar Arguments ..


REAL, INTENT(IN)                         :: A( LDA, * )
INTEGER, INTENT(IN OUT)                  :: LDA
INTEGER, INTENT(IN)                      :: N
INTEGER, INTENT(IN)                      :: IPVT( * )
REAL, INTENT(OUT)                        :: RCOND
REAL, INTENT(OUT)                        :: Z( * )


!     ..
!     .. Array Arguments ..



!     ..
!     .. Local Scalars ..

INTEGER :: INFO, J, K, KB, KP1, L
REAL :: ANORM, EK, S, SM, T, WK, WKM, YNORM
!     ..
!     .. External Functions ..

REAL :: SASUM_DIS, SDOT_DIS
EXTERNAL  SASUM_DIS, SDOT_DIS
!     ..
!     .. External Subroutines ..

EXTERNAL  SAXPY_DIS, SGEFA, SSCAL_DIS
!     ..
!     .. Intrinsic Functions ..

INTRINSIC ABS, MAX, SIGN
!     ..


!                        ** compute 1-norm of A
ANORM  = 0.0E0
DO  J = 1, N
  ANORM  = MAX( ANORM, SASUM_DIS( N,A( 1,J ),1 ) )
END DO
!                                      ** factor

CALL SGEFA( A, LDA, N, IPVT, INFO )

!     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
!     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E .
!     trans(A) is the transpose of A.  The components of E  are
!     chosen to cause maximum local growth in the elements of W  where
!     trans(U)*W = E.  The vectors are frequently rescaled to avoid
!     overflow.

!                        ** solve trans(U)*W = E
EK = 1.0E0

DO  J = 1, N
  Z( J ) = 0.0E0
END DO


DO  K = 1, N
  
  IF( Z( K ) /= 0.0E0 ) EK = SIGN( EK, -Z( K ) )
  
  IF( ABS( EK - Z( K ) ) > ABS( A( K,K ) ) ) THEN
    
    S  = ABS( A( K,K ) ) / ABS( EK - Z( K ) )
    
    CALL SSCAL_DIS( N, S, Z, 1 )
    
    EK = S*EK
    
  END IF
  
  WK   = EK - Z( K )
  WKM  = -EK - Z( K )
  S    = ABS( WK )
  SM   = ABS( WKM )
  
  IF( A( K,K ) /= 0.0E0 ) THEN
    
    WK   = WK / A( K, K )
    WKM  = WKM / A( K, K )
    
  ELSE
    
    WK   = 1.0E0
    WKM  = 1.0E0
    
  END IF
  
  KP1  = K + 1
  
  IF( KP1 <= N ) THEN
    
    DO  J = KP1, N
      SM     = SM + ABS( Z( J ) + WKM*A( K,J ) )
      Z( J ) = Z( J ) + WK*A( K, J )
      S      = S + ABS( Z( J ) )
    END DO
    
    IF( S < SM ) THEN
      
      T  = WKM - WK
      WK = WKM
      
      DO  J = KP1, N
        Z( J ) = Z( J ) + T*A( K, J )
      END DO
      
    END IF
    
  END IF
  
  Z( K ) = WK
  
END DO


S  = 1.0E0 / SASUM_DIS( N, Z, 1 )

CALL SSCAL_DIS( N, S, Z, 1 )
!                                ** solve trans(L)*Y = W
DO  KB = 1, N
  K  = N + 1 - KB
  
  IF( K < N ) Z( K ) = Z( K ) + SDOT_DIS( N - K, A( K+1, K ), 1, Z( K+1 ), 1)
  
  IF( ABS( Z( K ) ) > 1.0E0 ) THEN
    
    S  = 1.0E0 / ABS( Z( K ) )
    
    CALL SSCAL_DIS( N, S, Z, 1 )
    
  END IF
  
  L      = IPVT( K )
  T      = Z( L )
  Z( L ) = Z( K )
  Z( K ) = T
END DO


S  = 1.0E0 / SASUM_DIS( N, Z, 1 )

CALL SSCAL_DIS( N, S, Z, 1 )
!                                 ** solve L*V = Y
YNORM  = 1.0E0

DO  K = 1, N
  L      = IPVT( K )
  T      = Z( L )
  Z( L ) = Z( K )
  Z( K ) = T
  
  IF( K < N ) CALL SAXPY_DIS( N - K, T, A( K + 1,K ), 1, Z( K + 1 ), 1 )
  
  IF( ABS( Z( K ) ) > 1.0E0 ) THEN
    
    S  = 1.0E0 / ABS( Z( K ) )
    
    CALL SSCAL_DIS( N, S, Z, 1 )
    
    YNORM  = S*YNORM
  END IF
  
END DO


S  = 1.0E0 / SASUM_DIS( N, Z, 1 )

CALL SSCAL_DIS( N, S, Z, 1 )
!                                  ** solve  U*Z = V
YNORM  = S*YNORM

DO  KB = 1, N
  
  K  = N + 1 - KB
  
  IF( ABS( Z( K ) ) > ABS( A( K,K ) ) ) THEN
    
    S  = ABS( A( K,K ) ) / ABS( Z( K ) )
    
    CALL SSCAL_DIS( N, S, Z, 1 )
    
    YNORM  = S*YNORM
    
  END IF
  
  IF( A( K,K ) /= 0.0E0 ) Z( K ) = Z( K ) / A( K, K )
  
  IF( A( K,K ) == 0.0E0 ) Z( K ) = 1.0E0
  
  T  = -Z( K )
  
  CALL SAXPY_DIS( K - 1, T, A( 1,K ), 1, Z( 1 ), 1 )
  
END DO
!                                   ** make znorm = 1.0
S  = 1.0E0 / SASUM_DIS( N, Z, 1 )

CALL SSCAL_DIS( N, S, Z, 1 )

YNORM  = S*YNORM

IF( ANORM /= 0.0E0 ) RCOND = YNORM / ANORM
IF( ANORM == 0.0E0 ) RCOND = 0.0E0

END SUBROUTINE SGECO

SUBROUTINE SGEFA( A, LDA, N, IPVT, INFO )

!         Factors a real matrix by Gaussian elimination.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!     SGEFA is usually called by SGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (time for SGECO) = (1 + 9/N) * (time for SGEFA) .

!     Input:  same as SGECO

!     On return:

!        A,IPVT  same as SGECO

!        INFO    INTEGER
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that SGESL or SGEDI will divide by zero
!                     if called.  Use  RCOND  in SGECO for a reliable
!                     indication of singularity.

! ---------------------------------------------------------------------

!     .. Scalar Arguments ..


REAL, INTENT(IN OUT)                     :: A( LDA, * )
INTEGER, INTENT(IN OUT)                  :: LDA
INTEGER, INTENT(IN OUT)                  :: N
INTEGER, INTENT(OUT)                     :: IPVT( * )
INTEGER, INTENT(OUT)                     :: INFO

!     ..
!     .. Array Arguments ..



!     ..
!     .. Local Scalars ..

INTEGER :: J, K, KP1, L, NM1
REAL :: T
!     ..
!     .. External Functions ..

INTEGER :: ISAMAX_DIS
EXTERNAL  ISAMAX_DIS
!     ..
!     .. External Subroutines ..

EXTERNAL  SAXPY_DIS, SSCAL_DIS
!     ..


!                      ** Gaussian elimination with partial pivoting
INFO = 0
NM1  = N - 1

DO  K = 1, NM1
  
  KP1  = K + 1
!                                            ** find L = pivot index
  
  L  = ISAMAX_DIS( N - K + 1, A( K,K ), 1 ) + K - 1
  IPVT( K ) = L
  
  IF( A( L,K ) == 0.0E0 ) THEN
!                                     ** zero pivot implies this column
!                                     ** already triangularized
    INFO = K
    
  ELSE
!                                     ** interchange if necessary
    IF( L /= K ) THEN
      
      T         = A( L, K )
      A( L, K ) = A( K, K )
      A( K, K ) = T
      
    END IF
!                                     ** compute multipliers
    T  = -1.0E0 / A( K, K )
    
    CALL SSCAL_DIS( N - K, T, A( K + 1,K ), 1 )
    
!                              ** row elimination with column indexing
    DO  J = KP1, N
      
      T  = A( L, J )
      
      IF( L /= K ) THEN
        
        A( L, J ) = A( K, J )
        A( K, J ) = T
        
      END IF
      
      CALL SAXPY_DIS( N-K, T, A( K+1, K ), 1, A( K+1, J ), 1 )
      
    END DO
    
  END IF
  
END DO


IPVT( N ) = N
IF( A( N,N ) == 0.0E0 ) INFO = N

END SUBROUTINE SGEFA

SUBROUTINE SGESL( A, LDA, N, IPVT, B, JOB )

!         Solves the real system
!            A * X = B  or  transpose(A) * X = B
!         using the factors computed by SGECO or SGEFA.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!     On entry

!        A       REAL(LDA, N)
!                the output from SGECO or SGEFA.

!        LDA     INTEGER
!                the leading dimension of the array  A

!        N       INTEGER
!                the order of the matrix  A

!        IPVT    INTEGER(N)
!                the pivot vector from SGECO or SGEFA.

!        B       REAL(N)
!                the right hand side vector.

!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  transpose(A)*X = B

!     On return

!        B       the solution vector  X

!     Error condition

!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically, this indicates singularity,
!        but it is often caused by improper arguments or improper
!        setting of LDA.  It will not occur if the subroutines are
!        called correctly and if SGECO has set RCOND .GT. 0.0
!        or SGEFA has set INFO .EQ. 0 .

!     To compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call sgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue

! ---------------------------------------------------------------------

!     .. Scalar Arguments ..


REAL, INTENT(IN)                         :: A( LDA, * )
INTEGER, INTENT(IN OUT)                  :: LDA
INTEGER, INTENT(IN)                      :: N
INTEGER, INTENT(IN)                      :: IPVT( * )
REAL, INTENT(IN OUT)                     :: B( * )
INTEGER, INTENT(IN)                      :: JOB

!     ..
!     .. Array Arguments ..



!     ..
!     .. Local Scalars ..

INTEGER :: K, KB, L, NM1
REAL :: T
!     ..
!     .. External Functions ..

REAL :: SDOT_DIS
EXTERNAL  SDOT_DIS
!     ..
!     .. External Subroutines ..

EXTERNAL  SAXPY_DIS
!     ..


NM1  = N - 1

IF( JOB == 0 ) THEN
!                                 ** solve  A * X = B
  
!                                     ** first solve  L*Y = B
  DO  K = 1, NM1
    
    L  = IPVT( K )
    T  = B( L )
    
    IF( L /= K ) THEN
      
      B( L ) = B( K )
      B( K ) = T
      
    END IF
    
    CALL SAXPY_DIS( N - K, T, A( K+1, K ), 1, B( K+1 ), 1 )
    
  END DO
!                                    ** now solve  U*X = Y
  DO  KB = 1, N
    
    K      = N + 1 - KB
    B( K ) = B( K ) / A( K, K )
    T      = - B( K )
    
    CALL SAXPY_DIS( K-1, T, A( 1, K ), 1, B(1), 1 )
    
  END DO
  
  
ELSE
!                         ** solve  trans(A) * X = B
  
!                                    ** first solve  trans(U)*Y = B
  DO  K = 1, N
    
    T      = SDOT_DIS( K - 1, A( 1,K ), 1, B( 1 ), 1 )
    B( K ) = ( B( K ) - T ) / A( K, K )
    
  END DO
  
!                                    ** now solve  trans(l)*x = y
  DO  KB = 1, NM1
    
    K      = N - KB
    B( K ) = B( K ) + SDOT_DIS( N - K, A( K+1, K ), 1, B( K+1 ), 1)
    L      = IPVT( K )
    
    IF( L /= K ) THEN
      
      T      = B( L )
      B( L ) = B( K )
      B( K ) = T
      
    END IF
    
  END DO
  
END IF

END SUBROUTINE SGESL

REAL FUNCTION SASUM_DIS( N, SX, INCX )

!  INPUT--    N  Number of elements in vector to be summed
!            SX  Sing-prec array, length 1+(N-1)*INCX, containing vector
!          INCX  Spacing of vector elements in SX

!  OUTPUT-- SASUM_DIS   Sum from 0 to N-1 of  ABS(SX(1+I*INCX))
! ----------------------------------------------------------

!     .. Scalar Arguments ..


INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN)                         :: SX( * )
INTEGER, INTENT(IN)                      :: INCX

!     ..
!     .. Array Arguments ..


!     ..
!     .. Local Scalars ..

INTEGER :: I, M
!     ..
!     .. Intrinsic Functions ..

INTRINSIC ABS, MOD
!     ..

SASUM_DIS  = 0.0

IF( N <= 0 ) RETURN

IF( INCX /= 1 ) THEN
!                                          ** non-unit increments
  DO  I = 1, 1 + ( N - 1 )*INCX, INCX
    SASUM_DIS  = SASUM_DIS + ABS( SX( I ) )
  END DO
  
ELSE
!                                          ** unit increments
  M  = MOD( N, 6 )
  
  IF( M /= 0 ) THEN
!                             ** clean-up loop so remaining vector
!                             ** length is a multiple of 6.
    DO  I = 1, M
      SASUM_DIS  = SASUM_DIS + ABS( SX( I ) )
    END DO
    
  END IF
!                              ** unroll loop for speed
  DO  I = M + 1, N, 6
    SASUM_DIS  = SASUM_DIS + ABS( SX( I ) ) + ABS( SX( I + 1 ) ) +  &
        ABS( SX( I + 2 ) ) + ABS( SX( I + 3 ) ) +  &
        ABS( SX( I + 4 ) ) + ABS( SX( I + 5 ) )
  END DO
  
END IF

END FUNCTION SASUM_DIS

SUBROUTINE SAXPY_DIS( N, SA, SX, INCX, SY, INCY )

!          Y = A*X + Y  (X, Y = vectors, A = scalar)

!  INPUT--
!        N  Number of elements in input vectors X and Y
!       SA  Single precision scalar multiplier A
!       SX  Sing-prec array containing vector X
!     INCX  Spacing of elements of vector X in SX
!       SY  Sing-prec array containing vector Y
!     INCY  Spacing of elements of vector Y in SY

! OUTPUT--
!       SY   For I = 0 to N-1, overwrite  SY(LY+I*INCY) with
!                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
!            where LX = 1          if INCX .GE. 0,
!                     = (-INCX)*N  if INCX .LT. 0
!            and LY is defined analogously using INCY.
! ------------------------------------------------------------

!     .. Scalar Arguments ..


INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN)                         :: SA
REAL, INTENT(IN)                         :: SX( * )
INTEGER, INTENT(IN)                      :: INCX
REAL, INTENT(OUT)                        :: SY( * )
INTEGER, INTENT(IN)                      :: INCY


!     ..
!     .. Array Arguments ..


!     ..
!     .. Local Scalars ..

INTEGER :: I, IX, IY, M
!     ..
!     .. Intrinsic Functions ..

INTRINSIC MOD
!     ..


IF( N <= 0 .OR. SA == 0.0 ) RETURN

IF( INCX == INCY .AND. INCX > 1 ) THEN
  
  DO  I = 1, 1 + ( N - 1 )*INCX, INCX
    SY( I ) = SY( I ) + SA*SX( I )
  END DO
  
ELSE IF( INCX == INCY .AND. INCX == 1 ) THEN
  
!                                        ** equal, unit increments
  M  = MOD( N, 4 )
  
  IF( M /= 0 ) THEN
!                            ** clean-up loop so remaining vector length
!                            ** is a multiple of 4.
    DO  I = 1, M
      SY( I ) = SY( I ) + SA*SX( I )
    END DO
    
  END IF
!                              ** unroll loop for speed
  DO  I = M + 1, N, 4
    SY( I ) = SY( I ) + SA*SX( I )
    SY( I + 1 ) = SY( I + 1 ) + SA*SX( I + 1 )
    SY( I + 2 ) = SY( I + 2 ) + SA*SX( I + 2 )
    SY( I + 3 ) = SY( I + 3 ) + SA*SX( I + 3 )
  END DO
  
  
ELSE
!               ** nonequal or nonpositive increments.
  IX = 1
  IY = 1
  IF( INCX < 0 ) IX = 1 + ( N - 1 )*( -INCX )
  IF( INCY < 0 ) IY = 1 + ( N - 1 )*( -INCY )
  
  DO  I = 1, N
    SY( IY ) = SY( IY ) + SA*SX( IX )
    IX = IX + INCX
    IY = IY + INCY
  END DO
  
END IF

END SUBROUTINE SAXPY_DIS

REAL FUNCTION SDOT_DIS( N, SX, INCX, SY, INCY )

!        Single-prec dot product of vectors  X  and  Y

!  INPUT--
!        N  Number of elements in input vectors X and Y
!       SX  Sing-prec array containing vector X
!     INCX  Spacing of elements of vector X in SX
!       SY  Sing-prec array containing vector Y
!     INCY  Spacing of elements of vector Y in SY

! OUTPUT--
!     SDOT_DIS   Sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
!            where  LX = 1          if INCX .GE. 0,
!                      = (-INCX)*N  if INCX .LT. 0,
!            and LY is defined analogously using INCY.
! ------------------------------------------------------------------

!     .. Scalar Arguments ..


INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN)                         :: SX( * )
INTEGER, INTENT(IN)                      :: INCX
REAL, INTENT(IN)                         :: SY( * )
INTEGER, INTENT(IN)                      :: INCY

!     ..
!     .. Array Arguments ..


!     ..
!     .. Local Scalars ..

INTEGER :: I, IX, IY, M
!     ..
!     .. Intrinsic Functions ..

INTRINSIC MOD
!     ..


SDOT_DIS = 0.0

IF( N <= 0 ) RETURN

IF( INCX == INCY .AND. INCX > 1 ) THEN
  
  DO  I = 1, 1 + ( N - 1 )*INCX, INCX
    SDOT_DIS = SDOT_DIS + SX( I )*SY( I )
  END DO
  
  
ELSE IF( INCX == INCY .AND. INCX == 1 ) THEN
  
!                                        ** equal, unit increments
  M  = MOD( N, 5 )
  
  IF( M /= 0 ) THEN
!                            ** clean-up loop so remaining vector length
!                            ** is a multiple of 4.
    DO  I = 1, M
      SDOT_DIS = SDOT_DIS + SX( I )*SY( I )
    END DO
    
  END IF
!                              ** unroll loop for speed
  DO  I = M + 1, N, 5
    SDOT_DIS = SDOT_DIS + SX( I )*SY( I ) + SX( I + 1 )*SY( I + 1 ) +  &
        SX( I + 2 )*SY( I + 2 ) + SX( I + 3 )*SY( I + 3 ) +  &
        SX( I + 4 )*SY( I + 4 )
  END DO
  
ELSE
!               ** nonequal or nonpositive increments.
  IX = 1
  IY = 1
  
  IF( INCX < 0 ) IX = 1 + ( N - 1 )*( -INCX )
  IF( INCY < 0 ) IY = 1 + ( N - 1 )*( -INCY )
  
  DO  I = 1, N
    SDOT_DIS = SDOT_DIS + SX( IX )*SY( IY )
    IX   = IX + INCX
    IY   = IY + INCY
  END DO
  
END IF

END FUNCTION SDOT_DIS

SUBROUTINE SSCAL_DIS( N, SA, SX, INCX )

!         Multiply vector SX by scalar SA

!  INPUT--  N  Number of elements in vector
!          SA  Single precision scale factor
!          SX  Sing-prec array, length 1+(N-1)*INCX, containing vector
!        INCX  Spacing of vector elements in SX

! OUTPUT-- SX  Replace  SX(1+I*INCX)  with  SA * SX(1+I*INCX)
!                for I = 0 to N-1
! ---------------------------------------------------------------------

!     .. Scalar Arguments ..


INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN)                         :: SA
REAL, INTENT(OUT)                        :: SX( * )
INTEGER, INTENT(IN)                      :: INCX


!     ..
!     .. Array Arguments ..


!     ..
!     .. Local Scalars ..

INTEGER :: I, M
!     ..
!     .. Intrinsic Functions ..

INTRINSIC MOD
!     ..


IF( N <= 0 ) RETURN

IF( INCX /= 1 ) THEN
  
  DO  I = 1, 1 + ( N - 1 )*INCX, INCX
    SX( I ) = SA*SX( I )
  END DO
  
  
ELSE
  
  M  = MOD( N, 5 )
  
  IF( M /= 0 ) THEN
!                           ** clean-up loop so remaining vector length
!                           ** is a multiple of 5.
    DO  I = 1, M
      SX( I ) = SA*SX( I )
    END DO
    
  END IF
!                             ** unroll loop for speed
  DO  I = M + 1, N, 5
    SX( I ) = SA*SX( I )
    SX( I + 1 ) = SA*SX( I + 1 )
    SX( I + 2 ) = SA*SX( I + 2 )
    SX( I + 3 ) = SA*SX( I + 3 )
    SX( I + 4 ) = SA*SX( I + 4 )
  END DO
  
END IF

END SUBROUTINE SSCAL_DIS

SUBROUTINE SSWAP_DIS( N, SX, INCX, SY, INCY )

!          Interchange s.p vectors  X  and  Y, as follows:

!     For I = 0 to N-1, interchange  SX(LX+I*INCX) and SY(LY+I*INCY),
!     where LX = 1          if INCX .GE. 0,
!              = (-INCX)*N  if INCX .LT. 0
!     and LY is defined analogously using INCY.


!  INPUT--
!        N  Number of elements in input vectors X and Y
!       SX  Sing-prec array containing vector X
!     INCX  Spacing of elements of vector X in SX
!       SY  Sing-prec array containing vector Y
!     INCY  Spacing of elements of vector Y in SY

! OUTPUT--
!       SX  Input vector SY (unchanged if N .LE. 0)
!       SY  Input vector SX (unchanged IF N .LE. 0)
! --------------------------------------------------------------

!     .. Scalar Arguments ..


INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN OUT)                     :: SX( * )
INTEGER, INTENT(IN)                      :: INCX
REAL, INTENT(IN OUT)                     :: SY( * )
INTEGER, INTENT(IN)                      :: INCY

!     ..
!     .. Array Arguments ..


!     ..
!     .. Local Scalars ..

INTEGER :: I, IX, IY, M
REAL :: STEMP1, STEMP2, STEMP3
!     ..
!     .. Intrinsic Functions ..

INTRINSIC MOD
!     ..


IF( N <= 0 ) RETURN

IF( INCX == INCY .AND. INCX > 1 ) THEN
  
  DO  I = 1, 1 + ( N-1 )*INCX, INCX
    STEMP1 = SX( I )
    SX( I ) = SY( I )
    SY( I ) = STEMP1
  END DO
  
  
ELSE IF( INCX == INCY .AND. INCX == 1 ) THEN
  
!                                        ** equal, unit increments
  M  = MOD( N, 3 )
  
  IF( M /= 0 ) THEN
!                            ** clean-up loop so remaining vector length
!                            ** is a multiple of 3.
    DO  I = 1, M
      STEMP1 = SX( I )
      SX( I ) = SY( I )
      SY( I ) = STEMP1
    END DO
    
  END IF
!                              ** unroll loop for speed
  DO  I = M + 1, N, 3
    STEMP1 = SX( I )
    STEMP2 = SX( I + 1 )
    STEMP3 = SX( I + 2 )
    SX( I ) = SY( I )
    SX( I + 1 ) = SY( I + 1 )
    SX( I + 2 ) = SY( I + 2 )
    SY( I ) = STEMP1
    SY( I + 1 ) = STEMP2
    SY( I + 2 ) = STEMP3
  END DO
  
  
ELSE
!               ** nonequal or nonpositive increments.
  IX = 1
  IY = 1
  
  IF( INCX < 0 ) IX = 1 + ( N - 1 )*( -INCX )
  IF( INCY < 0 ) IY = 1 + ( N - 1 )*( -INCY )
  
  DO  I = 1, N
    STEMP1 = SX( IX )
    SX( IX ) = SY( IY )
    SY( IY ) = STEMP1
    IX   = IX + INCX
    IY   = IY + INCY
  END DO
  
END IF

END SUBROUTINE SSWAP_DIS

INTEGER FUNCTION ISAMAX_DIS( N, SX, INCX )

! INPUT--  N     Number of elements in vector of interest
!          SX    Sing-prec array, length 1+(N-1)*INCX, containing vector
!          INCX  Spacing of vector elements in SX

! OUTPUT-- ISAMAX_DIS   First I, I = 1 to N, to maximize
!                         ABS(SX(1+(I-1)*INCX))
! ---------------------------------------------------------------------

!     .. Scalar Arguments ..


INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN)                         :: SX( * )
INTEGER, INTENT(IN)                      :: INCX

!     ..
!     .. Array Arguments ..


!     ..
!     .. Local Scalars ..

INTEGER :: I, II
REAL :: SMAX, XMAG
!     ..
!     .. Intrinsic Functions ..

INTRINSIC ABS
!     ..


IF( N <= 0 ) THEN
  
  ISAMAX_DIS = 0
  
ELSE IF( N == 1 ) THEN
  
  ISAMAX_DIS = 1
  
ELSE
  
  SMAX = 0.0
  II   = 1
  
  DO  I = 1, 1 + ( N-1 )*INCX, INCX
    
    XMAG = ABS( SX( I ) )
    
    IF( SMAX < XMAG ) THEN
      
      SMAX   = XMAG
      ISAMAX_DIS = II
      
    END IF
    
    II = II + 1
    
  END DO
  
END IF

END FUNCTION ISAMAX_DIS

!************************************************************************
!ccccccc sergio shorter faster version
!     COMPUTES PLANCK FUNCTION

REAL FUNCTION PLKAVG ( WNUMLO, WNUMHI, T )


REAL, INTENT(IN OUT)                     :: WNUMLO
REAL, INTENT(IN)                         :: WNUMHI
REAL, INTENT(IN OUT)                     :: T
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'


REAL :: PLKAVG_ORIG

REAL :: wvn,ttorad

IF (T < 1.0E-4) THEN
  PLKAVG = 0.0
ELSE IF (wnumhi-wnumlo < 0) THEN
  WRITE(*,*) 'what!!!! wnumhi < wnmlo!!!!'
  STOP
ELSE IF (wnumhi-wnumlo > 1.0) THEN
!!!self test purposes; do an integral
  PLKAVG = PLKAVG_ORIG(WNUMLO, WNUMHI, T)
ELSE   !!!kCARTA purposes : monochromatic
  wvn  =  (wnumlo+wnumhi)/2.0
  plkavg = ttorad(wvn,T)
END IF

RETURN
END FUNCTION PLKAVG

!************************************************************************
! slftest : when working
! in disort nana  0.000000000000000  50000.0000000000
!  T T T T F T T T
! 0 1 4 1 4 1 1
! 1.000E-004  0.700  300.  1.00  3.1415926500  1.00  152.585283866739  90.0
! 0.000  1.52728600717066  28.3722254400388  0.900  0.800  100.  0.500  0.866
! 0.500  47.8655714383180  5.00  0.000
