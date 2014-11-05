** MemSys4 original source modules
** mem4.for

      SUBROUTINE MEM4
     * (ST, MEMRUN,SX,TESTX,CHISQX,SCALEX,PLOWX,PHIGHX,PDEVX,
     *             GLOWX,GHIGHX,GDEVX,OMEGAX,ALPHAX,ISTATX,NTRNSX)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*                                                                      *
*           -----------------------------------------------            *
*          |                                               |           *
*          | "MemSys4"  Quantified Maximum Entropy System  |           *
*          |                                               |           *
*           -----------------------------------------------            *
*                                                                      *
*                 Conceptualised and programmed by                     *
*            Stephen Gull, Mark Charter and John Skilling              *
*                                                                      *
*         (C)   Copyright 1990,1991                                    *
*               Maximum Entropy Data Consultants Ltd                   *
*               33 North End, Meldreth, Royston, England               *
*                                                                      *
*               Version 0.92            15 March 1991                  *
************************************************************************
*
* Purpose:
*        Perform one iterate of hidden-space maximum entropy,
*        and find the new distribution, with statistical
*        information such as
*               Evidence = log Prob(Data)     (base e)
*
*   Enter with:
*               <1> = hidden distribution h   (unless MEMRUN=1)
*               <3> = default model           (if DEF <= 0)
*              <21> = data D
*              <22> = accuracies = 1/sigma    (if ACC <= 0)
*
*    Exit with:
*               <1> = new hidden distribution h
*              <11> = new visible distribution f = ICF * h
*              <24> = normalised residuals = (D-Transform(old f))/sigma
*     and optionally
*              <22> = Poisson accuracies = sqrt(D+1)/Transform(old f)
*
*    Workspace:
*               <2> , <5>,
*              <23>, <25> , <26> , <27> , <28>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    MEMRUN  I     I      -         1 is  start,  > 1 is continue
*    SX      R       O    -         Entropy of old distribution h
*    TESTX   R       O    -         1-cos(angle( w, D-Transform(old f)))
*    CHISQX  R       O    -         Chisquared of old distribution f
*    SCALEX  R       O    -         estimate of overall scaling of
*                                   sigma(data)  (=1 if not being found)
*    PLOWX   R       O    -         Numerical lower limit on Evidence
*    PHIGHX  R       O    -         Numerical upper limit on Evidence
*    PDEVX   R       O    -         Random vector std dev of Evidence
*    GLOWX   R       O    -         Numerical lower limit on Good
*    GHIGHX  R       O    -         Numerical upper limit on Good
*    GDEVX   R       O    -         Random vector std dev of Good
*    OMEGAX  R       O    -         value of stopping criterion,
*                                   rising to 1.00
*    ALPHAX  R       O    -         alpha used to calculate new h
*    ISTATX  I       O    -         Composite status code
*    NTRNSX  I       O    -         Accumulated number of transforms
* Local externals:
*    LMAX    I     I      -         Maximum number of iterations
*    METHD0  I     I      -         Switch for stopping criterion
*    METHD1  I     I      -         Switch for entropy type
*    METHD2  I     I      -         Switch for likelihood type
*    AIM     R     I      -         Ultimate target of constraint >= 0
*    UTOL    R     I      -         Tolerance (dimensionless O(.1)) used
*                                   for maximisations & to define CODEs
*    NTRANS  I     I O    -         Accumulated number of transforms
*    ALPHA   R     I O    -         Regularisation parameter
*    SCALE   R       O    -         Estimate of scaling of sigma(data)
*    HHIGH   R     I O    -         Upper limit for bubble mismatch
*    ISTAT   I     I O    -         Composite status code
*    TOL     R     I      -         Arithmetical tolerance
*    RATE    R     I      -         Dimensionless distance limit, O(1)
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
*
* Areas:
*     <1>          I                Hidden distribution h
*                                     (not input if MEMRUN=1)
*     <2>           (O)             Workspace
*     <3>          I                Default model m (if DEF <= 0)
*     <5>            O              Hidden distribution h
*    <11>            O              New visible distribution f = ICF * h
*
*    <21>          I                Data D
*    <22>          I(O)             Accuracies = 1/sigma (if ACC <= 0)
*                                     (Poisson = sqrt(D+1)/Transform(f))
*    <24>            O              Normalised residuals =
*                                       (D-Transform(old f))/sigma
*    <25>           (O)             Workspace
*    <26>           (O)             Workspace
*    <27>           (O)             Workspace
*    <28>           (O)             Workspace
*
* External calls:
*       Hidden-Visible-Data
*    MECOPY               Copy area
*    MEMCGH    h  f  d    Conjugate gradient in hidden (dual) space
*    MEMCGP    h  f  d    Conjugate gradient in data space
*    MEMCGY    h  f  d    Conjugate gradient in hidden (dual) space
*    MEMCHI          d    Likelihood and gradient (Gaussian errors)
*    MEMDET               Eigenstructure of symmetric tridiagonal matrix
*    MEMDOT               Dot product
*    MEMENT    h          Entropy, gradS and SQRT(metric)
*    MEMICF    h  f       Apply ICF
*    MEMLA                Alpha control
*    MEMLA0    h          Alpha control if MEMRUN = 1
*    MEMLAR               Interpolate omega(alpha) table
*    MEMLAW               Insert new entry into omega(alpha) table
*    MEMLB                Beta control
*    MEMOP     h  f  d    Apply scaled forward transform
*    MEMPOI          d    Likelihood and gradient (Poisson errors)
*    MEMSMA               Area scalar multiply and add
*    MEMSMD               Area scalar multiply, add and dot product
*    MEMTR     h  f  d    Apply scaled backward transform
*    MEMUPD    h          Update area within limits
*    MEPRB0    h  f  d    alpha*good for "Classic" if MEMRUN = 1
*    MEPROB    h  f  d    LogProb and Good degrees of freedom
*    MESCAL    h  f  d    Statistics for input h
*    MESMUL               Area scalar multiply
*    MESURE    h  f  d    Conjugate gradient determination of t
*    METEST          d    Compares gradS, gradL
*    MEVMUL               Area vector multiply
*    MEVRND    h     d    Random vector
*    MEZERO               Initialise area to zero
*
*   Called directly from this routine:
*    MCHECK      Check control parameters
*    MCLEAR      Clear read/write flags
*    MECOPY      Copy area
*    MEHSET      Set up initial h if MEMRUN = 1
*    MEMCGY      Conjugate gradient in hidden space
*    MEMICF      Apply ICF
*    MEMLA       Alpha control
*    MEMLA0      Alpha control if MEMRUN=1
*    MEMSMA      Area scalar multiply and add
*    MESCAL      Statistics for input h
*    MEZERO      Initialise area to zero
*    UINFO       User's diagnostics handler
* 
* 
* Parameter interpretations:
*
* Method(2) = 1 is "Gaussian likelihood"
* Method(2) = 2 is "Poisson likelihood"
*
* Method(1) = 1 is "Standard entropy"
* Method(1) = 2 is "Positive/negative entropy"
* Method(1) = 3 is "Fermi-Dirac entropy"
* Method(1) = 4 is "Quadratic regularisation"
*
* Method(0) = 1 is "Classic with known noise"
* Method(0) = 2 is "Classic with unknown noise"
* Method(0) = 3 is "Alpha given"
* Method(0) = 4 is "Historic Chisquared=Ndata"
*
*    ISTATX =  64*code6 + 32*code5 + 16*code4 + 8*code3 + 4*code2
*                                             + 2*code1 +   code0
*
*     code6 = 0 is "conjugate gradient returns GOOD to within UTOL"
*     code6 = 1 is "conjgrad failure to maximise accurately"
*
*     code5 = 0 is "conjugate gradient returns dy,H,|dr| within UTOL"
*     code5 = 1 is "conjgrad failure to maximise accurately"
*
*     code4 = 0 is "test <= 1"
*     code4 = 1 is "test >  1 , unacceptable"
*
*     code3 = 0 is "Alpha control is satisfied by current alpha"
*     code3 = 1 is "Alpha control changes alpha"
*
*     code2 = 0 is "Beta control does not impose distance penalty"
*     code2 = 1 is "Beta control imposes distance penalty"
*
*     code1 = 0 is "stopping criterion omega = 1 reached within UTOL"
*     code1 = 1 is "stopping criterion on omega not reached"
*
*     code0 = 0 is "bubble overlaps true bubble(alpha) to within UTOL"
*     code0 = 1 is "bubble not overlapping true bubble(alpha)"
*
*  but if Omega > 1 on setup, ISTATX = 0 regardless of other code values
*
*
* History:
*    MKC/JS          10 Dec 1990     First release
*    JS               8 Mar 1991     RATE in beta control
*    MKC             15 Mar 1991     In line with MemSys5 v1.11
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER MEMRUN,NTRNSX,ISTATX
      REAL ST(0:*)
      REAL ALPHAX,SX,TESTX,CHISQX,SCALEX,OMEGAX,
     * PLOWX,PHIGHX,PDEVX,GLOWX,GHIGHX,GDEVX
      INCLUDE 'memsys.inc'

      INTEGER METHD0,METHD1,METHD2
      INTEGER NTRANS,ISTAT
      REAL UTOL,TOL,RATE,AIM,ALPHA,SCALE,HHIGH
      EQUIVALENCE (METHD0,ICOM(KMETH0)),(METHD1,ICOM(KMETH1))
      EQUIVALENCE (METHD2,ICOM(KMETH2))
      EQUIVALENCE (NTRANS,ICOM(KNTRNS)),(ISTAT, ICOM(KISTAT))
      EQUIVALENCE (UTOL,  RCOM(KUTOL )),(TOL,   RCOM(KTOL  ))
      EQUIVALENCE (RATE,  RCOM(KRATE )),(AIM,   RCOM(KAIM  ))
      EQUIVALENCE (ALPHA, RCOM(KALPHA)),(SCALE, RCOM(KSCALE))
      EQUIVALENCE (HHIGH, RCOM(KHHIGH))

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL LCODE0,LCODE1,LCODE2,LCODE3,LCODE4,LCODE5,LCODE6
      REAL H,HLOW,DIST
      REAL S,TEST,CHISQ,OMEGA,SUMMET,GRADL,VAR,
     * GLOW,GHIGH,GDEV,PLOW,PHIGH,PDEV
      REAL ONE,TWO
      PARAMETER (ONE=1.0D0,TWO=2.0D0)

      WRITE (INFO,'(''  MemSys4 Version 0.92 15 Mar 1991   Method'',
     *              3(1X,I1))') METHD2,METHD1,METHD0
      CALL UINFO(INFO,1)
      CALL MCHECK
      CALL MCLEAR
      IF ((MEMRUN.LT.1).OR.(4.LT.MEMRUN)) STOP ' Illegal MEMRUN value'
      IF (MEMRUN.EQ.1) THEN
        ISTAT=0
        NTRANS=0
        CALL MEHSET(ST)
        IF(METHD2.EQ.2) CALL MEZERO(ST,22)
      ELSE
        CALL MECOPY(ST,1,5)
      ENDIF
      CALL MESCAL(ST,MEMRUN,METHD0,METHD2,AIM,UTOL,ALPHA,SCALE,
     *            S,TEST,CHISQ,PLOW,PHIGH,PDEV,GLOW,GHIGH,GDEV,HHIGH,
     *            OMEGA,VAR,SUMMET,GRADL,LCODE0,LCODE1,LCODE4,LCODE6)
      LCODE2=(MOD(ISTAT/4,2).EQ.0)
      ISTAT=0
      IF (.NOT.LCODE0) ISTAT=ISTAT+ 1
      IF (.NOT.LCODE1) ISTAT=ISTAT+ 2
      IF (.NOT.LCODE4) ISTAT=ISTAT+16
      IF (.NOT.LCODE6) ISTAT=ISTAT+64
      IF ((MEMRUN.EQ.1).AND.(OMEGA.GE.ONE)) THEN
        OMEGA=0.5
      ENDIF
      IF (MEMRUN.EQ.4) THEN
        WRITE (INFO,'(''  Omega   === '',0PF10.6)') OMEGA
        CALL UINFO(INFO,10)
      ELSE
        IF (MEMRUN.EQ.1) THEN
          CALL MEMLA0(RATE,SUMMET,GRADL,OMEGA,VAR,ALPHA,LCODE3)
        ELSEIF (MEMRUN.EQ.2) THEN
          CALL MEMLA(ALPHA,RATE,SUMMET,OMEGA,TOL,VAR,GRADL,.FALSE.,
     *                  LCODE2,LCODE3,LCODE4)
        ELSEIF (MEMRUN.EQ.3) THEN
          CALL MEMLA(ALPHA,RATE,SUMMET,OMEGA,TOL,VAR,GRADL,.TRUE.,
     *                  LCODE2,LCODE3,LCODE4)
        ENDIF
        CALL MEMSMA(ST,2,-ALPHA,4,4)
        CALL MEMCGY(ST,LMAX,ALPHA,UTOL,TOL,RATE,SUMMET,
     *                             HLOW,HHIGH,DIST,LCODE2,LCODE5)
        H=(HLOW+HHIGH)/TWO
        WRITE (INFO,'('' H      ='',1PG14.7,'' = ['',1PG14.7,'' ,'',
     *              1PG14.7,'' ]'')') H,HLOW,HHIGH
        CALL UINFO(INFO,20)
        WRITE (INFO,'(''  Omega   === '',0PF10.6,
     *             ''      dist ==='',0PF7.4,
     *              ''    Alpha ==='',1PE12.4)') OMEGA,DIST,ALPHA
        CALL UINFO(INFO,10)
        IF (.NOT.LCODE2) ISTAT=ISTAT+ 4
        IF (.NOT.LCODE3) ISTAT=ISTAT+ 8
        IF (.NOT.LCODE5) ISTAT=ISTAT+32
      ENDIF
      CALL MECOPY(ST,5,1)
      CALL MEMICF(ST,1,11)
      WRITE (INFO,'(''  Ntrans  === '',I6,25X,
     *              ''    Code  ===  '',6I1.1)') NTRANS,
     *                 MOD(ISTAT/32,2),MOD(ISTAT/16,2),MOD(ISTAT/ 8,2),
     *                 MOD(ISTAT/ 4,2),MOD(ISTAT/ 2,2),MOD(ISTAT   ,2)
      CALL UINFO(INFO,10)
      NTRNSX=NTRANS
      ISTATX=ISTAT
      ALPHAX=ALPHA
      SX=S
      TESTX=TEST
      CHISQX=CHISQ
      SCALEX=SCALE
      PLOWX=PLOW
      PHIGHX=PHIGH
      PDEVX=PDEV
      GLOWX=GLOW
      GHIGHX=GHIGH
      GDEVX=GDEV
      OMEGAX=OMEGA
      END

      SUBROUTINE MESCAL(ST,MEMRUN,METHD0,METHD2,AIM,UTOL,ALPHA,SCALE,
     *           S,TEST,CHISQ,PLOW,PHIGH,PDEV,GLOW,GHIGH,GDEV,HHIGH,
     *           OMEGA,VAR,SUMMET,GRADL,LCODEH,LCODEO,LCODET,LCODEG)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*          Calculate scalar statistics for current h,
*          and set up metric, normalised residuals, -gradS and -gradL.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    MEMRUN  I     I      -         1 is  start,  > 1 is continue
*    METHD0  I     I      -         Switch for stopping criterion
*    METHD2  I     I      -         Switch for type of errors
*    AIM     R     I      -         Ultimate target of constraint >= 0
*    UTOL    R     I      -         Tolerance (dimensionless O(.1))
*    ALPHA   R     I      -         Regularisation parameter
*    SCALE   R       O    -         Estimate of scaling of sigma(data)
*    S       R       O    -         Entropy of old distribution h
*    TEST    R       O    -         1-cos(angle( w, D-Transform(old f)))
*    CHISQ   R       O    -         Chisquared of old distribution
*    PLOW    R       O    -         Numerical lower limit on Evidence
*    PHIGH   R       O    -         Numerical upper limit on Evidence
*    PDEV    R       O    -         Random vector std dev of Evidence
*    GLOW    R       O    -         Numerical lower limit on Good
*    GHIGH   R       O    -         Numerical upper limit on Good
*    GDEV    R       O    -         Random vector std dev of Good
*    HHIGH   R     I      -         Upper limit for bubble mismatch
*    OMEGA   R       O    -         value of stopping criterion,
*                                   rising to 1.00
*    VAR     R       O    -         Estimated variance of OMEGA
*    SUMMET  R       O    -         SUM([metric])
*    GRADL   R       O    -         gradL
*    LCODEH  L       O    -         Flag for bubble match
*    LCODEO  L       O    -         Flag for OMEGA = 1
*    LCODET  L       O    -         Flag for TEST <= 1
*    LCODEG  L       O    -         Flag for GOOD accuracy
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <1>            O              SQRT([metric])
*     <2>            O              -gradS
*     <3>          I                Default model m
*     <4>            O              -gradL
*     <5>          I                Hidden-space distribution h
*    <11>           (O)                      workspace
*    <21>          I                Data D
*    <22>          I(O)             [acc] (output if Poisson)
*    <24>            O              Normalised residuals
*    <25>           (O)                      workspace
*    <26>           (O)                      workspace
*    <28>           (O)                      workspace
*
* External calls:
*    MEMCHI      Likelihood and gradient (Gaussian errors)
*    MEMDOT      Dot product
*    MEMENT      Entropy, gradS and SQRT(metric)
*    MEMOP       Mock data from h
*    MEMPOI      Likelihood and gradient (Poisson errors)
*    MEMTR       Initialising backward transform
*    MEPRB0      alpha*good for "Classic" if alpha = infinity
*    MEPROB      LogProb and Good degrees of freedom
*    METEST      Compares gradS, gradL
*    MEVMUL      Area vector multiply
*    MEZERO      Initialise area to zero
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS         10 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL LCODEH,LCODEO,LCODET,LCODEG
      INTEGER MEMRUN,METHD0,METHD2
      REAL ST(0:*)
      REAL AIM,UTOL,ALPHA,SCALE,HHIGH,S,SUMMET,GRADL,
     * CHISQ,TEST,GLOW,GHIGH,GDEV,PLOW,PHIGH,PDEV,OMEGA,VAR

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      INTEGER ICODEG
      REAL ALHOOD,UNITS,DIV,XXX,YYY,PROBL,GOOD,
     * DETL,DLOW,DHIGH,DDEV,ALFGD,ALFS,ALF2S,DATA
      REAL ZERO,EPS,ONE,TWO,TWELVE,PILOG
* PILOG = log(sqrt(2*pi))
      PARAMETER (ZERO=0.0D0,EPS=1.0D-20,ONE=1.0D0,TWO=2.0D0,
     * TWELVE=12.0D0,PILOG=0.9189385332046727417803296D0)

      DIV(XXX,YYY) = XXX/MAX(YYY,EPS)

      CALL UINFO('   MeScal',1)
      CALL MEMOP(ST,5,25,.FALSE.)
      IF (METHD2.EQ.1) THEN
        CALL MEMCHI(ST,ALHOOD,DATA,UNITS)
      ELSEIF (METHD2.EQ.2) THEN
        CALL MEMPOI(ST,ALHOOD,DATA,UNITS)
      ENDIF
	      UNITS=UNITS-DATA*PILOG
      CALL MEMTR(ST,24,4,.FALSE.)
      CALL MEMENT(ST,S,SUMMET)
      IF(MEMRUN.EQ.1) THEN
        S=ZERO
        CALL MEZERO(ST,2)
      ENDIF
      CALL METEST(ST,TEST)
      LCODET=(TEST.LE.ONE)
      CALL MEVMUL(ST,4,1,2)
      CALL MEMDOT(ST,2,2,ALF2S)
      GRADL=SQRT(ALF2S)
      IF(MEMRUN.EQ.1) THEN
        ALF2S=-ALF2S/TWO
        ALFS=ZERO
      ELSE
        ALFS=ALPHA*S
      ENDIF
      IF (METHD0.EQ.2) THEN
        SCALE=SQRT(TWO*(ALHOOD-ALFS)/DATA)
      ELSE
        SCALE=ONE
      ENDIF
      CHISQ=DIV(TWO*ALHOOD,SCALE*SCALE)
      WRITE(INFO,'(''  Entropy === '',1PE12.4,
     *              ''    Test ==='',0PF7.4,
     *             ''    Chisq ==='',1PE12.4)') S,TEST,CHISQ
      CALL UINFO(INFO,10)
      IF (METHD0.EQ.2) THEN
        WRITE (INFO,'(49X,''Scale ==='',1PE12.4)') SCALE
        CALL UINFO(INFO,10)
      ENDIF
* Stopping criterion OMEGA
      LCODEH=.TRUE.
      LCODEG=.TRUE.
      GLOW=ZERO
      GHIGH=ZERO
      GDEV=ZERO
      DLOW=ZERO
      DHIGH=ZERO
      DDEV=ZERO
      PLOW=ZERO
      PHIGH=ZERO
      PDEV=ZERO
      IF ((METHD0.EQ.1).OR.(METHD0.EQ.2).OR.(METHD0.EQ.3)) THEN
        IF(MEMRUN.EQ.1) THEN
          CALL MEPRB0(ST, ALFGD)
        ELSE
          CALL MEPROB(ST,ALPHA,UTOL,
     *                GLOW,GHIGH,GDEV,DLOW,DHIGH,DDEV,LCODEG)
          LCODEH=(HHIGH/(SCALE*SCALE).LE.UTOL*GLOW)
        ENDIF
        PLOW=UNITS-DATA*LOG(SCALE)-DHIGH/TWO+
     *                         ALFS/(SCALE*SCALE)-CHISQ/TWO
        PHIGH=UNITS-DATA*LOG(SCALE)-DLOW/TWO+
     *                         ALFS/(SCALE*SCALE)-CHISQ/TWO
        PROBL=(PLOW+PHIGH)/TWO
        PDEV=DDEV/TWO
        GOOD=(GLOW+GHIGH)/TWO
        DETL=(DLOW+DHIGH)/TWO
        WRITE (INFO,'('' Good   ='',1PG14.7,'' = ['',1PG14.7,'' ,'',
     *               1PG14.7,'' ] +-'',1PG14.7)') GOOD,GLOW,GHIGH,GDEV
        CALL UINFO(INFO,20)
        WRITE (INFO,'('' LogDet ='',1PG14.7,'' = ['',1PG14.7,'' ,'',
     *               1PG14.7,'' ] +-'',1PG14.7)') DETL,DLOW,DHIGH,DDEV
        CALL UINFO(INFO,20)
        IF(LCODEG) THEN
          ICODEG=0
        ELSE
          ICODEG=1
        ENDIF
        WRITE (INFO,'(''  LogProb ==='',1PE13.4,''    Code ===   '',I1,
     *             ''       Good  ==='',1PE12.4)') PROBL,ICODEG,GOOD   
        CALL UINFO(INFO,10)
        IF (METHD0.EQ.3) THEN
          IF(MEMRUN.EQ.1) THEN
            OMEGA=ZERO
          ELSE  
            OMEGA=DIV(AIM,ALPHA)
          ENDIF
          VAR=UTOL*UTOL/TWELVE
        ELSE
	    IF(MEMRUN.EQ.1) THEN
            OMEGA=DIV(ALFGD*SCALE*SCALE*AIM,-TWO*ALF2S)
            VAR=DIV(TEST,ONE-TEST)+UTOL*UTOL/TWELVE
          ELSEIF(MEMRUN.EQ.2) THEN
            OMEGA=DIV(GOOD*SCALE*SCALE*AIM,-TWO*ALFS)
            VAR=DIV(TEST,ONE-TEST)+UTOL*UTOL/TWELVE
          ENDIF
        ENDIF
      ELSEIF (METHD0.EQ.4) THEN
        OMEGA=DIV(DATA*AIM,CHISQ)
        VAR=DIV(TEST,ONE-TEST)+UTOL*UTOL/TWELVE
      ENDIF
      LCODEO=(ABS(OMEGA-ONE).LE.UTOL)
      IF(MEMRUN.EQ.1) THEN
        S=ZERO
        CALL MEZERO(ST,2)
      ELSE
        CALL MEMENT(ST,S,SUMMET)
      ENDIF
      END

      SUBROUTINE MEPROB(ST,ALPHA,UTOL,
     *                  GLOW,GHIGH,GDEV,DLOW,DHIGH,DDEV,LCODEG)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*   Purpose:
*        Calculate      Number of good degrees of freedom
*              and      log( prob ( Data | form of transform ) )
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    ALPHA   R     I      -         alpha
*    UTOL    R     I      -         User's termination criterion
*    GLOW    R       O    -         Numerical lower limit on Good
*    GHIGH   R       O    -         Numerical upper limit on Good
*    GDEV    R       O    -         Random vector std dev of Good
*    DLOW    R       O    -         Numerical lower limit on LogDet
*    DHIGH   R       O    -         Numerical upper limit on LogDet
*    DDEV    R       O    -         Random vector std dev of LogDet
*    LCODEG  L       O    -         Return code
* Local externals:
*    LMAX    I     I      -         Maximum number of iterations
*    NRAND   I     I      -         Number of random vectors to be used
*    ISEED   I     I      -         Seed for reproducible random vectors
*    TOL     R     I      -         Arithmetical tolerance
* Areas:
*     <1>          I                SQRT([metric])
*     <2>           (O)               (transform workspace)
*    <11>           (O)               (transform workspace)
*    <22>         (I)               [acc]
*    <25>           (O)               (transform workspace)
*    <26>           (O)             random vector
*    <28>           (O)             workspace
*
* External calls:
*    MEMCGP      Conjugate gradient in data space
*    MEMDET      Eigenstructure of symmetric tridiagonal matrix
*    MEVRND      Random sign generator
*    UINFO       User's diagnostics handler
*    VRAND0      Initialise random sign generator
*
* History:
*    MKC/JS         19 Dec 1990     First release
*                    4 Feb 1991     Overflow risk reduced
*
*  Notes:
*  (1) Each call to MEPROB can be as expensive as a MaxEnt iterate
*  (2) Inaccurate estimates are weighted down severely
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL LCODEG
      REAL ST(0:*)
      REAL ALPHA,UTOL,GLOW,GHIGH,GDEV,DLOW,DHIGH,DDEV
      INCLUDE 'memsys.inc'

      INTEGER NRAND,ISEED
      REAL TOL
      EQUIVALENCE (NRAND, ICOM(KNRAND)),(ISEED, ICOM(KISEED))
      EQUIVALENCE (TOL,   RCOM(KTOL  ))

      LOGICAL LG
      INTEGER I,J,M,N
      REAL WEIGHT,W,X,Y,TRACE,TOL2
      REAL GOODLO,GOODHI,GOODEV
      REAL DETLLO,DETLHI,DETDEV
      REAL GAM(0:LMAX),DEL(0:LMAX),
     *                 VAL(0:LMAX),VEC(0:LMAX),OFF(0:LMAX)
      REAL ZERO,SMALL,ONE,TWO
      PARAMETER (ZERO=0.0D0,SMALL=1.0D-6,ONE=1.0D0,TWO=2.0D0)

*  GOOD = [GLOW,GHIGH]+or-GDEV = < [GOODLO,GOODHI]+or-GOODEV >AvrgeNRAND
*  DETL = [DLOW,DHIGH]+or-DDEV = < [DETLLO,DETLHI]+or-DETDEV >AvrgeNRAND
      CALL UINFO('   MeProb',1)
      CALL VRAND0(ISEED)
      TOL2=SQRT(TOL)
      WEIGHT=ZERO
      GLOW=ZERO
      GHIGH=ZERO
      GDEV=ZERO
      DLOW=ZERO
      DHIGH=ZERO
      DDEV=ZERO
      LCODEG=.FALSE.
      DO 6 I=1,NRAND
        CALL MEVRND(ST,26,26,0)
        GOODLO=ZERO
        GOODEV=ZERO
        DETLLO=ZERO
        DETDEV=ZERO
        GOODHI=ZERO
        DETLHI=ZERO
        CALL MEMCGP(ST,LMAX,ALPHA,UTOL,TOL, GAM,DEL,M,N,LG)
        IF(M.GT.0) THEN
* Lower limit and variance for Good and LogDet
          TRACE=ZERO
          DO 2 J=0,M-1
            VAL(J)=GAM(J)*DEL(J)
            OFF(J+1)=-GAM(J+1)*DEL(J)
            TRACE=TRACE+VAL(J)
    2     CONTINUE
          VAL(M)=ZERO
          CALL MEMDET(TOL,VAL,VEC,OFF,M+1,1,-1)
          DO 3 J=0,M-1
            X=VAL(J)/ALPHA
            Y=VEC(J)**2/ALPHA
            GOODLO=GOODLO+Y/(ONE+X)
            GOODEV=GOODEV+(Y/(ONE+X))*(X/(ONE+X))
            IF(X.GT.TOL2) THEN
              DETLLO=DETLLO+Y*LOG(ONE+X)/X
              DETDEV=DETDEV+Y*(LOG(ONE+X))**2/X
            ELSE
              DETLLO=DETLLO+Y
              DETDEV=DETDEV+Y*X
            ENDIF
    3     CONTINUE
          X=GAM(0)*DEL(0)
          X=(GAM(0)*X)**2
          GOODLO=X*GOODLO
          GOODEV=TWO*X*GOODEV
          DETLLO=X*DETLLO
          DETDEV=TWO*X*DETDEV
        ENDIF
* Upper limit for Good and LogDet
* (if VAL(M) unavailable because M>N, worst case is VAL(M) big)
        VAL(M)=TRACE/TOL
        VAL(0)=GAM(0)*DEL(0)
        DO 4 J=1,N
          VAL(J)=GAM(J)*DEL(J)
          OFF(J)=-GAM(J)*DEL(J-1)
    4   CONTINUE
        CALL MEMDET(TOL,VAL,VEC,OFF,M+1,1,1)
        DO 5 J=0,M
          X=VAL(J)/ALPHA
          GOODHI=GOODHI+VEC(J)**2*X/(ONE+X)
          IF(X.GT.TOL2) THEN
            DETLHI=DETLHI+VEC(J)**2*LOG(ONE+X)
          ELSE
            DETLHI=DETLHI+VEC(J)**2*X
          ENDIF
    5   CONTINUE
        GOODHI=GAM(0)**2*GOODHI
        DETLHI=GAM(0)**2*DETLHI
* Increment totals
        IF (LG) THEN
          W=ONE
        ELSE
          W=SMALL
        ENDIF
        GLOW=GLOW+GOODLO*W
        GHIGH=GHIGH+GOODHI*W
        GDEV=GDEV+GOODEV*W
        DLOW=DLOW+DETLLO*W
        DHIGH=DHIGH+DETLHI*W
        DDEV=DDEV+DETDEV*W
        WEIGHT=WEIGHT+W
        LCODEG=(LCODEG.OR.LG)
    6 CONTINUE

      GLOW=GLOW/WEIGHT
      GHIGH=GHIGH/WEIGHT
      GDEV=SQRT(GDEV)/WEIGHT

      DLOW=DLOW/WEIGHT
      DHIGH=DHIGH/WEIGHT
      DDEV=SQRT(DDEV)/WEIGHT
      END

      SUBROUTINE MEPRB0(ST,ALFGD)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*  Purpose:
*           Calculate      alpha * good
*           for "Classic" MaxEnt if alpha = infinity
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    ALFGD   R       O    -         alpha*good
* Local externals:
*    NRAND   I     I      -         Number of random vectors to be used
*    ISEED   I     I      -         Seed for reproducible random vectors
* Areas:
*     <1>          I                SQRT([metric])
*     <2>           (O)               (transform workspace)
*    <11>           (O)               (transform workspace)
*    <22>         (I)               [acc]
*    <25>           (O)               (transform workspace)
*    <26>           (O)             random vector
*
* External calls:
*    MEMDOT      Area scalar product
*    MEMTR       Apply SQRT([metric]) * TrICF * Tropus * [acc]
*    MEVRND      Random sign generator
*    UINFO       User's diagnostics handler
*    VRAND0      Initialise random sign generator
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
*  Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      REAL ALFGD
      INCLUDE 'memsys.inc'

      INTEGER NRAND,ISEED
      EQUIVALENCE (NRAND, ICOM(KNRAND)),(ISEED, ICOM(KISEED))

      INTEGER I
      REAL AG,ZERO
      PARAMETER (ZERO=0.0D0)

      CALL UINFO('   MePrb0',1)
      CALL VRAND0(ISEED)
      ALFGD=ZERO
      DO 3 I=1,NRAND
        CALL MEVRND(ST,26,26,0)
        CALL MEMTR(ST,26,2,.TRUE.)
        CALL MEMDOT(ST,2,2,AG)
        ALFGD=ALFGD+AG
    3 CONTINUE
      ALFGD=ALFGD/NRAND
      END

** extra.for

      SUBROUTINE MEMSET(METHDX,NRANDX,ISEEDX,AIMX,RATEX,DEFX,ACCX,UTOLX)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                Set MemSys control parameters
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    METHDX  I     I     0:2        Switch array for stopping criterion,
*                                   default model, noise type
*    NRANDX  I     I      -         Number of random vectors to be used
*    ISEEDX  I     I      -         Seed for reproducible random vectors
*    AIMX    R     I      -         Ultimate target of constraint >= 0
*    RATEX   R     I      -         Dimensionless distance limit, O(1)
*    DEFX    R     I      -         Default level (used only if > 0)
*    ACCX    R     I      -         Accuracies (used only if > 0)
*    UTOLX   R     I      -         Tolerance (dimensionless O(.1)) used
*                                   for maximisations & to define CODEs
* Local externals:
*    METHD0  I       O    -         Switch for stopping criterion
*    METHD1  I       O    -         Switch for entropy type
*    METHD2  I       O    -         Switch for likelihood type
*    NRAND   I       O    -         Number of random vectors to be used
*    ISEED   I       O    -         Seed for reproducible random vectors
*    AIM     R       O    -         Ultimate target of constraint >= 0
*    RATE    R       O    -         Dimensionless distance limit, O(1)
*    DEF     R       O    -         Default level (used only if > 0)
*    ACC     R       O    -         Accuracies (used only if > 0)
*    UTOL    R       O    -         Tolerance (dimensionless O(.1)) used
*                                   for maximisations & to define CODEs
*
* External calls:
*    MCHECK      Check control parameters
*    UINFO       User's diagnostics handler
*    VRAND0      Re-initialise random number generator to ISEED
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER METHDX(0:2),NRANDX,ISEEDX
      REAL RATEX,AIMX,DEFX,ACCX,UTOLX
      INCLUDE 'memsys.inc'

      INTEGER METHD0,METHD1,METHD2,NRAND,ISEED
      REAL RATE,AIM,DEF,ACC,UTOL
      EQUIVALENCE (METHD0,ICOM(KMETH0)),(METHD1,ICOM(KMETH1))
      EQUIVALENCE (METHD2,ICOM(KMETH2))
      EQUIVALENCE (NRAND, ICOM(KNRAND)),(ISEED, ICOM(KISEED))
      EQUIVALENCE (RATE,  RCOM(KRATE )),(AIM,   RCOM(KAIM  ))
      EQUIVALENCE (DEF,   RCOM(KDEF  )),(ACC,   RCOM(KACC  ))
      EQUIVALENCE (UTOL,  RCOM(KUTOL ))

      CALL UINFO('   MemSet',1)
      METHD0=METHDX(0)
      METHD1=METHDX(1)
      METHD2=METHDX(2)
      NRAND=NRANDX
      ISEED=ISEEDX
      RATE=RATEX
      AIM=AIMX
      DEF=DEFX
      ACC=ACCX
      UTOL=UTOLX
      CALL MCHECK
      CALL VRAND0(ISEED)
      END

      SUBROUTINE MEMGET(METHDX,NRANDX,ISEEDX,AIMX,RATEX,DEFX,ACCX,UTOLX)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*           Return current values of MemSys control parameters
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    METHDX  I       O   0:2        Switch array for stopping criterion,
*                                   default model, noise type
*    NRANDX  I       O    -         Number random vectors to be used
*    ISEEDX  I       O    -         Seed for reproducible random vectors
*    AIMX    R       O    -         Ultimate target of constraint >= 0
*    RATEX   R       O    -         Dimensionless distance limit, O(1)
*    DEFX    R       O    -         Default level (used only if > 0)
*    ACCX    R       O    -         Accuracies (used only if > 0)
*    UTOLX   R       O    -         Tolerance (dimensionless O(.1)) used
*                                   for maximisations & to define CODEs
* Local externals:
*    METHD0  I     I      -         Switch for stopping criterion
*    METHD1  I     I      -         Switch for entropy type
*    METHD2  I     I      -         Switch for likelihood type
*    NRAND   I     I      -         Number of random vectors to be used
*    ISEED   I     I      -         Seed for reproducible random vectors
*    AIM     R     I      -         Ultimate target of constraint >= 0
*    RATE    R     I      -         Dimensionless distance limit, O(1)
*    DEF     R     I      -         Default level (used only if > 0)
*    ACC     R     I      -         Accuracies (used only if > 0)
*    UTOL    R     I      -         Tolerance (dimensionless O(.1)) used
*                                   for maximisations & to define CODEs
*
* External calls:
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER METHDX(0:2),NRANDX,ISEEDX
      REAL RATEX,AIMX,DEFX,ACCX,UTOLX
      INCLUDE 'memsys.inc'

      INTEGER METHD0,METHD1,METHD2,NRAND,ISEED
      REAL RATE,AIM,DEF,ACC,UTOL
      EQUIVALENCE (METHD0,ICOM(KMETH0)),(METHD1,ICOM(KMETH1))
      EQUIVALENCE (METHD2,ICOM(KMETH2))
      EQUIVALENCE (NRAND, ICOM(KNRAND)),(ISEED, ICOM(KISEED))
      EQUIVALENCE (RATE,  RCOM(KRATE )),(AIM,   RCOM(KAIM  ))
      EQUIVALENCE (DEF,   RCOM(KDEF  )),(ACC,   RCOM(KACC  ))
      EQUIVALENCE (UTOL,  RCOM(KUTOL ))

      CALL UINFO('   MemGet',1)
      METHDX(0)=METHD0
      METHDX(1)=METHD1
      METHDX(2)=METHD2
      NRANDX=NRAND
      ISEEDX=ISEED
      AIMX=AIM
      RATEX=RATE
      DEFX=DEF
      ACCX=ACC
      UTOLX=UTOL
      END

      SUBROUTINE MCHECK
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                Check MemSys control parameters
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Local externals:
*    METHD0  I     I      -         Switch for stopping criterion
*    METHD1  I     I      -         Switch for entropy type
*    METHD2  I     I      -         Switch for likelihood type
*    NRAND   I     I      -         Number of random vectors to be used
*    AIM     R     I      -         Ultimate target of constraint >= 0
*    RATE    R     I      -         Dimensionless distance limit, O(1)
*    UTOL    R     I      -         Tolerance (dimensionless O(.1)) used
*                                   for maximisations & to define CODEs
*
* External calls:
*     -
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INCLUDE 'memsys.inc'

      INTEGER METHD0,METHD1,METHD2,NRAND
      REAL RATE,AIM,UTOL
      EQUIVALENCE (METHD0,ICOM(KMETH0)),(METHD1,ICOM(KMETH1))
      EQUIVALENCE (METHD2,ICOM(KMETH2))
      EQUIVALENCE (NRAND, ICOM(KNRAND))
      EQUIVALENCE (RATE,  RCOM(KRATE )),(AIM,   RCOM(KAIM  ))
      EQUIVALENCE (UTOL,  RCOM(KUTOL ))

      REAL ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)

      IF((METHD0.LT.1).OR.(4.LT.METHD0)) STOP ' Illegal METHOD(0) value'
      IF((METHD1.LT.1).OR.(4.LT.METHD1)) STOP ' Illegal METHOD(1) value'
      IF((METHD2.LT.1).OR.(2.LT.METHD2)) STOP ' Illegal METHOD(2) value'
      IF((METHD2.EQ.2) .AND. (METHD0.EQ.2))
     *               STOP ' Method incompatible with Poisson statistics'
      IF(((METHD0.EQ.1).OR.(METHD0.EQ.2)) .AND. (NRAND.LE.0))
     *                                     STOP ' Illegal NRAND value'
      IF (AIM.LT.ZERO)                     STOP ' Illegal AIM value'
      IF (RATE.LT.ZERO)                    STOP ' Illegal RATE value'
      IF ((UTOL.LT.ZERO).OR.(ONE.LT.UTOL)) STOP ' Illegal UTOL value'
      END

      SUBROUTINE MOVIE4(ST,NCORR, KSTAT,NTRNSX)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*   Purpose:
*         Obtain visible and hidden sample from posterior probability
*         bubble, optionally correlated with previous sample(s).
*
*   IF( NCORR > 0 ) THEN
*          set random normal vector as source v
*   ELSE
*          use <4> as supplied externally for source v
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    NCORR   I     I      -         Sample overlap (none if NCORR=1)
*    KSTAT   I       O    -         Return code
*    NTRNSX  I       O    -         Accumulated number of transforms
* Local externals:
*    UTOL    R     I      -         User's termination criterion
*    TOL     R     I      -         Arithmetic tolerance
*    ALPHA   R     I      -         Regularisation parameter
*    SCALE   R     I      -         Estimate of scaling of sigma(data)
*    NTRANS  I     I O    -         Accumulated number of transforms
* Areas:
*     <1>            O              Hidden bubble sample = h+(B**-1/2)v
*     <2>           (O)                      workspace
*     <3>          I                Default model m (if DEF <= 0)
*     <4>         (I)               (Source vector v, if NCORR=0)
*                    O              Hidden offset        = (B**-1/2)v
*     <5>          I                Hidden distribution h
*    <11>            O              Visible bubble sample = ICF*<1>
*    <22>         (I)               [acc]    (if ACC <= 0)
*    <25>           (O)                      workspace
*    <26>           (O)                      workspace
*    <27>           (O)                      workspace
*    <28>           (O)                      workspace
*
* External calls:
*    MCLEAR      Clear read/write flags (on entry)
*    MEMENT      Set up h, SQRT([metric]) and associated scalars
*    MEMICF      Apply ICF
*    MEMSMA      Area scalar multiply add
*    MESMUL      Area scalar multiply
*    MESURE      Conjugate gradient determination of t
*    MEVMUL      Area multiply
*    MEVRND      Area Gaussian
*    UINFO       User's diagnostics handler
*
* The return codes KSTAT are
*    0 is "offset vector was estimated accurately"
*    1 is "conjugate gradient failed in some way"
*    2 is "offset vector has wrong magnitude"
*    3 is "conjugate gradient failed and offset has wrong magnitude"
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
*  Notes:
*  (1) MEMICF may also be used to obtain visible-space distribution
*      corresponding to <4>
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER NCORR,KSTAT,NTRNSX
      REAL ST(0:*)
      INCLUDE 'memsys.inc'

      INTEGER NTRANS
      REAL UTOL,TOL,ALPHA,SCALE
      EQUIVALENCE (NTRANS,ICOM(KNTRNS))
      EQUIVALENCE (UTOL,  RCOM(KUTOL )),(TOL,   RCOM(KTOL  ))
      EQUIVALENCE (ALPHA, RCOM(KALPHA)),(SCALE, RCOM(KSCALE))

      LOGICAL LCODEC,LCODET
      REAL S,SUMMET,QLO,QHI,ONE
      PARAMETER (ONE=1.0D0)

      CALL UINFO('   Movie4',1)
      CALL MCLEAR
* <4> := random input vector, possibly correlated between samples
      IF( NCORR.GT.0 ) CALL MEVRND(ST,4,2, NCORR)
* <1> := SQRT([metric])
      CALL MEMENT(ST,S,SUMMET)
* <4> := B**(-1/2) <4>     ,    <2> = workspace
      CALL MESURE(ST,.TRUE.,UTOL,TOL,ALPHA, QLO,QHI,LCODEC,LCODET)
      KSTAT=0
      IF(.NOT.LCODEC) KSTAT=KSTAT+1
      IF(.NOT.LCODET) KSTAT=KSTAT+2
*   Rescale to include final SQRT([metric]) factor from <1>
      CALL MESMUL(ST,2,SCALE,2)
      CALL MEVMUL(ST,1,2,4)
* <1> := <5> + <4>  =  h + (B**-1/2) random = hidden bubble sample
      CALL MEMSMA(ST,4,ONE,5,1)
* <11>:= ICF * <1>
      CALL MEMICF(ST,1,11)
      NTRNSX=NTRANS
      END

      SUBROUTINE MASK4(ST,AMEANX,ERRLOX,ERRHIX,JSTAT,NTRNSX)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*   Purpose:
*               Estimate  SUM  h * (hidden mask)    with error bars
*                          i    i               i
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    AMEANX  R       O    -         mean of required sum
*    ERRLOX  R       O    -         lower limit on standard deviation
*    ERRHIX  R       O    -         upper limit on standard deviation 
*    JSTAT   I       O    -         Return code
*    NTRNSX  I       O    -         Accumulated number of transforms
* Local externals:
*    UTOL    R     I      -         User's termination criterion
*    TOL     R     I      -         Arithmetic tolerance
*    ALPHA   R     I      -         Regularisation parameter
*    SCALE   R     I      -         Estimate of scaling of sigma(data)
*    NTRANS  I     I O    -         Accumulated number of transforms
* Areas:
*     <1>           (O)                      workspace
*     <2>           (O)                      workspace
*     <3>          I                Default model m (if DEF <= 0)
*     <4>          I                Input hidden mask
*     <5>          I                Hidden distribution h
*    <11>           (O)                      workspace
*    <22>         (I)               [acc]    (if ACC <= 0)
*    <25>           (O)                      workspace
*    <26>           (O)                      workspace
*    <28>           (O)                      workspace
*
* External calls:
*    MCLEAR      Clear read/write flags (on entry)
*    MEMDOT      Dot product
*    MEMENT      Set up h, SQRT([metric]) and associated scalars
*    MESURE      Conjugate gradient determination of p.B**(-1).p
*    MEVMUL      Area multiply
*    UINFO       User's diagnostics handler
*
* The return codes are
*    0 is "standard deviation was estimated accurately"
*    1 is "standard deviation may not have been estimated accurately"
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
*  Notes:
*  (1) A call to MASK4 may be more expensive than a MaxEnt iterate.
*  (2) Calculation of standard deviation can be very ill-conditioned,
*      involving inverse of B without a factor A, and is prone to
*      serious rounding error.  The value of the return code should be
*      checked (zero indicates correct operation of the procedure).
*      The only remedy for serious error may be full double precision
*      operation.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER JSTAT,NTRNSX
      REAL ST(0:*)
      REAL AMEANX,ERRLOX,ERRHIX
      INCLUDE 'memsys.inc'

      INTEGER NTRANS
      REAL UTOL,TOL,ALPHA,SCALE
      EQUIVALENCE (NTRANS,ICOM(KNTRNS))
      EQUIVALENCE (UTOL,  RCOM(KUTOL )),(TOL,   RCOM(KTOL  ))
      EQUIVALENCE (ALPHA, RCOM(KALPHA)),(SCALE, RCOM(KSCALE))

      LOGICAL LCODEC,LCODET
      REAL AMEAN,ERRLO,ERRHI,VARLO,VARHI
      REAL S,SUMMET

      CALL UINFO('   Mask4',1)
      CALL MCLEAR
* <1> := SQRT([metric])
      CALL MEMENT(ST,S,SUMMET)
      CALL MEMDOT(ST,5,4,AMEAN)
      CALL MEVMUL(ST,1,4,4)
      CALL MESURE(ST,.FALSE.,UTOL,TOL,ALPHA, VARLO,VARHI,LCODEC,LCODET)
      JSTAT=0
      IF(.NOT.LCODEC) JSTAT=JSTAT+1
      ERRLO=SQRT(VARLO)*SCALE
      ERRHI=SQRT(VARHI)*SCALE
      AMEANX=AMEAN
      ERRLOX=ERRLO
      ERRHIX=ERRHI
      NTRNSX=NTRANS
      END

      SUBROUTINE MESURE(ST,BUBBLE,UTOL,TOL,ALPHA, QLO,QHI,LCODEC,LCODET)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*              ------------------------------------------ 
*             |          T                        -1     |
*         Set |   Q :=  b . ( Alpha*Identity + A )  . b  |
*             |         -                             -  |
*              ------------------------------------------ 
*    where      A = Memtr.Memop ,
*               Memop = ( [acc] *  Opus  *  ICF  * SQRT([metric]) *)
*               Memtr = ( SQRT([metric]) * TrICF * Tropus * [acc] *)
*
* IF ( BUBBLE ) then additionally
*              ------------------------------------------ 
*             |                               -1/2       |
*         set |   t :=  ( Alpha*Identity + A )    . b    |
*             |   -                                 -    |
*              ------------------------------------------ 
*         When b is a random normal vector, the statistics of t are
*                 T                          -1
*            < t t > = ( Alpha*Identity + A )
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    BUBBLE  L     I      -         Switch for "bubble" option
*    UTOL    R     I      -         User's termination criterion
*    TOL     R     I      -         Arithmetic tolerance
*    ALPHA   R     I      -         Regularisation parameter
*    QLO     R       O    -         Numerical lower limit on Q
*    QHI     R       O    -         Numerical upper limit on Q
*    LCODEC  L       O    -         Status code for conj. gradients
*    LCODET  L       O    -         Status code for bubble sample
* Local externals:
*    LMAX    I     I      -         Maximum number of iterations
* Areas:
*     <1>          I                SQRT([metric])
*     <2>            O              t   (only if BUBBLE)
*     <4>          I                b = input vector
*    <11>           (O)                      workspace
*    <22>         (I)               [acc]    (if ACC <= 0)
*    <25>           (O)                      workspace
*    <26>           (O)                      workspace
*    <27>           (O)                      workspace (only if BUBBLE)
*    <28>           (O)                      workspace
*
* External calls:
*    MEMCGH      Quantify
*    MEMDET      Diagonalisation
*
* The status codes are
*     LCODEC   .TRUE. : Conjugate gradient returned OK
*              .FALSE.: Conjugate gradient not satisfactory
* and LCODET   .TRUE. : Bubble not called or bubble sample OK
*              .FALSE.: Bubble offset has detectably incorrect length
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL BUBBLE,LCODEC,LCODET
      REAL ST(0:*)
      REAL UTOL,TOL,ALPHA,QLO,QHI
      INCLUDE 'memsys.inc'

      INTEGER I,J,K,M,LMAX2,N
      PARAMETER (LMAX2=LMAX*(LMAX+2))
      REAL GAM(0:LMAX),DEL(0:LMAX)
      REAL OFF(0:LMAX),VAL(0:LMAX),VEC(0:LMAX2)
      REAL TEMP,ZERO
      PARAMETER (ZERO=0.0D0)

* Calculate quantification scalars
      CALL MEMCGH(ST,LMAX,ALPHA,UTOL,TOL,.FALSE.,
     *            QLO,QHI,LCODEC,N,GAM,DEL,OFF)
      LCODET=.TRUE.
      IF (BUBBLE) THEN
* Calculate coefficients OFF() of gradient vectors g for random sample
        VAL(0)=GAM(0)*DEL(0)
        DO 1 I=1,N
          VAL(I)=GAM(I)*DEL(I)
          OFF(I)=-GAM(I)*DEL(I-1)
    1   CONTINUE
        M=N+1
        CALL MEMDET(TOL,VAL,VEC,OFF,M,M,1)
        DO 2 I=0,N
          OFF(I)=ZERO
    2   CONTINUE
        K=0
        DO 4 J=0,N
          TEMP=VEC(K)/SQRT(ALPHA+VAL(J))
          DO 3 I=0,N
            OFF(I)=OFF(I)+VEC(K)*TEMP*GAM(0)/GAM(I)
            K=K+1
    3     CONTINUE
    4   CONTINUE
* Build up output vector
        CALL MEMCGH(ST,LMAX,ALPHA,UTOL,TOL,.TRUE.,
     *              QLO,QHI,LCODET,N,GAM,DEL,OFF)
      ENDIF
      END

      SUBROUTINE MEMTRQ(ST,ISEED, ERRX)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*         Check for   OPUS*ICF / TRICF*TROPUS   inconsistencies.
*
*         Dimensionless return parameter ERR should be O(rounding error)
*         if OPUS*ICF and TRICF*TROPUS are each other's transposes.
*
* Method:
*         IF ( ISEED > 0 ) THEN
*           <1> := random   ,   <26> := random
*         ELSE
*           use <1> and <26> as set by user.
*
*                ABS(  <1> . TrICF*Tropus*<26>  -  <26> . Opus*ICF*<1> )
*         ERR := -------------------------------------------------------
*                SQRT( |<1>| |TrICF*Tropus*<26>| |<26>| |Opus*ICF*<1>| )
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    ISEED   I     I      -         Seed for random vector (or 0)
*    ERRX    R       O    -         Fractional error
* Areas:
*     <1>          I(O)             input random vector
*     <2>           (O)               (transform workspace)
*    <11>           (O)               (transform workspace)
*    <25>           (O)               (transform workspace)
*    <26>          I(O)             input random vector
*
* External calls:
*    ICF         Intrinsic correlation function (ICF)
*    MCLEAR      Clear read/write flags (on entry)
*    MECOPY      Copy area
*    MEFLAG      Set read/write flag
*    MEMDOT      Dot product
*    MEVRND      Random sign generator
*    OPUS        Visible-to-Data differential transform
*    TRICF       Transpose ICF operator
*    TROPUS      Data-to-Visible differential transform
*    UINFO       User's diagnostics handler
*    VRAND0      Initialise random sign generator
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER ISEED
      REAL ST(0:*)
      REAL ERRX

      REAL U,V,UU,VV,UTV,VRU,URRU,VTTV,ZERO,ERR
      PARAMETER (ZERO=0.0D0)

      CALL UINFO('   MemTrQ',1)
      CALL MCLEAR
      IF (ISEED.GT.0) THEN
*   Optional setup source areas <1> and <26>
        CALL VRAND0(ISEED)
        CALL MEVRND(ST,1,1,0)
        CALL MEVRND(ST,26,26,0)
      ENDIF
* <2> := TrICF * Tropus <26>
*  VV := <26>.<26>
      CALL MEMDOT(ST,26,26,VV)
      CALL MECOPY(ST,26,25)
      CALL TROPUS(ST,25,11)
      CALL MEFLAG(25, 1)
      CALL MEFLAG(11, 2)
      CALL MCLEAR
      CALL TRICF(ST,11,2)
      CALL MEFLAG(11, 1)
      CALL MEFLAG(2,  2)
      CALL MCLEAR
* UTV  := <1>.<2>
* VTTV := <2>.<2>
      CALL MEMDOT(ST,1,2,UTV)
      CALL MEMDOT(ST,2,2,VTTV)
* <25> := Opus * ICF <1>
*   UU := <1>.<1>
      CALL MEMDOT(ST,1,1,UU)
      CALL MECOPY(ST,1,2)
      CALL ICF(ST,2,11)
      CALL MEFLAG(2,  1)
      CALL MEFLAG(11, 2)
      CALL MCLEAR
      CALL OPUS(ST,11,25)
      CALL MEFLAG(11, 1)
      CALL MEFLAG(25, 2)
      CALL MCLEAR
* VRU  := <26>.<25>
* URRU := <25>.<25>
      CALL MEMDOT(ST,25,26,VRU)
      CALL MEMDOT(ST,25,25,URRU)
* Fractional error
      U=SQRT(UU)*SQRT(VV)
      V=SQRT(URRU)*SQRT(VTTV)
      ERR=ABS(UTV-VRU)/SQRT(U*V)
      ERRX=ERR
      END

** transfrm.for

      SUBROUTINE MEMICF(ST,M,N)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                      <N> :=  ICF * <M>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    M       I     I      -         Input hidden area number
*    N       I     I      -         Output visible area number
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <M>          I                Input hidden area
*     <N>            O              Output visible area
*     <2>           (O)               (transform workspace)
*    <11>           (O)               (transform workspace)
*
* External calls:
*    ICF         Intrinsic correlation function
*    MCLEAR      Clear read/write flags
*    MECOPY      Copy area
*    MEFLAG      Set read/write flag
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER M,N
      REAL ST(0:*)

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      WRITE (INFO,'(''   MemICF('',I2,'','',I2,'')'')') M,N
      CALL UINFO(INFO,1)
      CALL MECOPY(ST,M,2)
      CALL ICF(ST,2,11)
      CALL MEFLAG(2,  1)
      CALL MEFLAG(11, 2)
      CALL MCLEAR
      CALL MECOPY(ST,11,N)
      END

      SUBROUTINE MTRICF(ST,M,N)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                                 T
*                      <N> :=  ICF  * <M>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    M       I     I      -         Input visible area number
*    N       I     I      -         Output hidden area number
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <M>          I                Input visible area
*     <N>            O              Output hidden area
*     <2>           (O)               (transform workspace)
*    <11>           (O)               (transform workspace)
*
* External calls:
*    MCLEAR      Clear read/write flags
*    MECOPY      Copy area
*    MEFLAG      Set read/write flag
*    TRICF       Transpose ICF operator
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER M,N
      REAL ST(0:*)

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      WRITE (INFO,'(''   MTrICF('',I2,'','',I2,'')'')') M,N
      CALL UINFO(INFO,1)
      CALL MECOPY(ST,M,11)
      CALL TRICF(ST,11,2)
      CALL MEFLAG(11, 1)
      CALL MEFLAG(2,  2)
      CALL MCLEAR
      CALL MECOPY(ST,2,N)
      END

      SUBROUTINE MEMOP(ST,J,K,FLAG)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*         Differential forward transform with full scaling,
*          or generation of mock data
*
*     IF( FLAG ) THEN
*                  <K> :=  [acc] * Opus  * ICF * SQRT([metric]) * <J>
*     ELSE
*                  <K> :=          MemEx * ICF * <J>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    J       I     I      -         Input hidden area number
*    K       I     I      -         Output data area number
*    FLAG    L     I      -         Opus / MockData switch
* Local externals:
*    METHD2  I     I      -         Switch for likelihood type
*    NTRANS  I     I O    -         Transform counter
*    ACC     R     I      -         Accuracies (used only if > 0)
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <J>          I                Input hidden area
*     <K>            O              Output data area
*     <1>         (I)                 (metric)
*     <2>           (O)               (transform workspace)
*    <11>           (O)               (transform workspace)
*    <22>          I                [acc]
*    <25>           (O)               (transform workspace)
*
* External calls:
*    ICF         Intrinsic correlation function
*    MCLEAR      Clear read/write flags
*    MECOPY      Copy area
*    MEFLAG      Set read/write flag
*    MEMEX       Generate mock data
*    MESMUL      Area scalar multiply
*    MEVMUL      Area multiply
*    OPUS        Visible-to-Data differential transform
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*   (1)  ICF and Opus must NOT change frequently between iterates
*   (2)  The actual ICF transform is performed from <2> to <11>
*   (3)  The actual Opus transform is performed from <11> to <25>
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL FLAG
      INTEGER J,K
      REAL ST(0:*)
      INCLUDE 'memsys.inc'

      INTEGER METHD2,NTRANS
      EQUIVALENCE (METHD2,ICOM(KMETH2)),(NTRANS,ICOM(KNTRNS))
      REAL ACC
      EQUIVALENCE (ACC,   RCOM(KACC  ))

      REAL ZERO
      PARAMETER (ZERO=0.0D0)
      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      WRITE (INFO,'(''   MemOp ('',I2,'','',I2,'')'')') J,K
      CALL UINFO(INFO,1)
      CALL MECOPY(ST,J,2)
      IF(FLAG) CALL MEVMUL(ST,2,1,2)
      CALL ICF(ST,2,11)
      CALL MEFLAG(2,  1)
      CALL MEFLAG(11, 2)
      CALL MCLEAR
      IF(FLAG) THEN
        CALL OPUS(ST,11,25)
      ELSE
        CALL MEMEX(ST,11,25)
      ENDIF
      CALL MEFLAG(11, 1)
      CALL MEFLAG(25, 2)
      CALL MCLEAR
      IF(FLAG) THEN
        IF (METHD2.EQ.1 .AND. ACC.GT.ZERO) THEN
          CALL MESMUL(ST,25,ACC,25)
        ELSE
          CALL MEVMUL(ST,25,22,25)
        ENDIF
      ENDIF
      CALL MECOPY(ST,25,K)
      NTRANS=NTRANS+1
      END

      SUBROUTINE MEMTR(ST,K,J,FLAG)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*         Differential backward transform, with or without full scaling
*
*     IF( FLAG ) THEN
*                  <J> := SQRT([metric]) * TrICF * Tropus * [acc] * <K>
*     ELSE
*                  <J> :=                  TrICF * Tropus * [acc] * <K>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    K       I     I      -         Input data area number
*    J       I     I      -         Output hidden area number
*    FLAG    L     I      -         Tropus / SetDistribution switch
* Local externals:
*    METHD2  I     I      -         Switch for likelihood type
*    NTRANS  I     I O    -         Transform counter
*    ACC     R     I      -         Accuracies (used only if > 0)
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <K>          I                Input data area
*     <J>            O              Output hidden area
*     <1>         (I)                 (metric)
*     <2>           (O)               (transform workspace)
*    <11>           (O)               (transform workspace)
*    <22>          I                [acc]
*    <25>           (O)               (transform workspace)
*
* External calls:
*    MCLEAR      Clear read/write flags
*    MECOPY      Copy area
*    MEFLAG      Set read/write flag
*    MESMUL      Area scalar multiply
*    MEVMUL      Area multiply
*    TRICF       Transpose ICF operator
*    TROPUS      Data-to-Visible differential transform
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*   (1)  TrICF and Tropus must NOT change frequently between iterates
*        (see Opus)
*   (2)  The actual Tropus transform is performed from <25> to <11>
*   (3)  The actual TrICF transform is performed from <11> to <2>
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL FLAG
      INTEGER J,K
      REAL ST(0:*)
      INCLUDE 'memsys.inc'

      INTEGER METHD2,NTRANS
      EQUIVALENCE (METHD2,ICOM(KMETH2)),(NTRANS,ICOM(KNTRNS))
      REAL ACC
      EQUIVALENCE (ACC,   RCOM(KACC  ))

      REAL ZERO
      PARAMETER (ZERO=0.0D0)
      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      WRITE (INFO,'(''   MemTr ('',I2,'','',I2,'')'')') K,J
      CALL UINFO(INFO,1)
      CALL MECOPY(ST,K,25)
      IF (METHD2.EQ.1 .AND. ACC.GT.ZERO) THEN
        CALL MESMUL(ST,25,ACC,25)
      ELSE
        CALL MEVMUL(ST,25,22,25)
      ENDIF
      CALL TROPUS(ST,25,11)
      CALL MEFLAG(25, 1)
      CALL MEFLAG(11, 2)
      CALL MCLEAR
      CALL TRICF(ST,11,2)
      CALL MEFLAG(11, 1)
      CALL MEFLAG(2,  2)
      CALL MCLEAR
      IF(FLAG) CALL MEVMUL(ST,2,1,2)
      CALL MECOPY(ST,2,J)
      NTRANS=NTRANS+1
      END

** conjgrad.for

      SUBROUTINE MEMCGH( ST,LMAX,ALPHA,UTOL,TOL,FLAG,
     *                   QLO,QHI,LCODE, N,GAM,DEL,COEFF)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Hidden-space conjugate gradient
* -------------------------------
*
*   Estimate Q = [ QLO , QHI ]  and set up GAMma and DELta scalars
*     by maximising
*                         QB =  2 y.b  -  y.B.y    (scalars only)
*                        QAB = 2 z.A.b - z.A.B.z   (scalars only)
*
*   Slave these to master conjugate gradient maximisation of
*                         QA =  2 x.b  -  x.A.x    (scalars only)
*    where      B = alpha*I + A    ,  A = Memtr.Memop ,
*               Memop = ( [acc] *  Opus  *  ICF  * SQRT([metric]) *)
*               Memtr = ( SQRT([metric]) * TrICF * Tropus * [acc] *)
*
*   IF( FLAG ) then additionally construct the vector
*                       N
*           t(.)  :=   Sum COEFF(i) * g(.,i)
*                      i=0
*       repeating conjugate gradient to regenerate the gradients g,
*       (as simulated in dataspace).
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    LMAX    I     I      -         Maximum number of iterations
*    ALPHA   R     I      -         Regularisation parameter
*    UTOL    R     I      -         User's termination criterion
*    TOL     R     I      -         Arithmetic tolerance
*    FLAG    L     I      -         Output vector switch
*    QLO     R       O    -         Numerical lower limit on Q
*    QHI     R       O    -         Numerical upper limit on Q
*    LCODE   L       O    -         Status code
*    N       I       O    -         GAMma(i),DELta(i), for i=0,N
*    GAM     R       O   0:LMAX     Gamma from each iteration
*    DEL     R       O   0:LMAX     Delta from each iteration
*    COEFF   R     I     0:N        Coefficients of g (only if FLAG)
* Areas:
*     <1>          I                SQRT([metric])
*     <2>            O              t = output (if FLAG), and workspace
*     <4>          I                b = input vector
*    <22>         (I)               [acc]
*    <25>           (O)             uu (workspace for transforms)
*    <26>           (O)             gg = simulated gradient  (workspace)
*    <27>           (O)             tt = simulated output t  (if FLAG)
*    <28>           (O)             hh = simulated conjugate (workspace)
*
* External calls:
*    MEMDOT      Dot product
*    MEMOP       Apply ( [acc] *  Opus  *  ICF  * SQRT([metric]) *)
*    MEMSMA      Update area
*    MEMSMD      Update area with scalar product
*    MEMTR       Apply ( SQRT([metric]) * TrICF * Tropus * [acc] *)
*    MESMUL      Multiply area by scalar
*    MEZERO      Initialise area to zero
*    UINFO       User's diagnostics handler
*
* IF( .NOT.FLAG ) the status code relates to the conjugate gradient: 
*     LCODE = ( ICODE <= 2 )
* where the internal codes are
*     ICODE = 0: User's termination criterion satisfied      (OK)
*     ICODE = 1: Convergence to within arithmetic TOLerance  (OK)
*     ICODE = 2: Zero results from TINY input vector b     (Warning)
*     ICODE = 3: Iteration limit LMAX exceeded             (Warning)
*     ICODE = 4: Probable error in transform routines       (Error)
*
* IF( FLAG ) the status code relates to the output vector length
*     LCODE = ( t.t is within a fraction UTOL of QLO )
* this being a consistency test when called from subroutine MESURE.
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*  (1) Conjugate gradient termination is either when the user's
*      criterion is satisfied, by
*         (1.+UTOL) * alpha * QB + QAB   passing   (length of b)**2
*      or when L passes
*         LMAX = maximum number of matrix applications
*      or when gamma underflows the internal arithmetic TOLerance.
*  (2) A repeat run (with FLAG=.TRUE.) to calculate the vector output
*      must reproduce the first run exactly, otherwise serious rounding
*      errors can accumulate.  This is why the same code is re-used.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL FLAG,LCODE
      INTEGER N,LMAX
      REAL ST(0:*)
      REAL ALPHA,UTOL,TOL,QLO,QHI
      REAL GAM(0:LMAX),DEL(0:LMAX),COEFF(0:LMAX)

      INTEGER L,ICODE
      REAL GAMMA0,GAMMA,DELTA, QB,EPSB,PHIB,DELB,T0
      REAL QAB,EPSAB,PHIAB,DELAB, ZERO,ONE,TINY,TEMP,SCALE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY=1.0D-30)

      CALL UINFO('   MemCGH',1)
* Initialise
      L=0
      QB=ZERO
      QAB=ZERO
      N=0
      GAM(0)=ZERO
      DEL(0)=ZERO
* SCALE arbitrarily (make input vector unit) to reduce overflow risk
      CALL MEMDOT(ST,4,4,TEMP)
      IF (TEMP.LT.TINY) THEN
        GAMMA0=ZERO
        T0=ZERO
        SCALE=ONE
        ICODE=2
        GOTO 9
      ENDIF
      SCALE=SQRT(TEMP)
      GAMMA0=(TEMP/SCALE)/SCALE
      GAMMA=GAMMA0
      GAM(L)=SQRT(GAMMA)*SCALE
      TEMP=ONE/SCALE
      CALL MESMUL(ST,4,TEMP,2)
      CALL MEZERO(ST,26)
      CALL MEZERO(ST,28)
      IF( FLAG ) THEN
        CALL MEZERO(ST,27)
        T0=COEFF(0)
      ENDIF
      EPSB=ONE
      PHIB=ALPHA/GAMMA
      EPSAB=GAMMA
      PHIAB=ZERO
      DELAB=ZERO
* Loop
    1 CONTINUE
        TEMP=ONE/GAMMA
        CALL MEMOP(ST,2,25,.TRUE.)
        CALL MEMSMD(ST,25,TEMP,28,DELTA)
        N=L
        DEL(L)=SQRT(DELTA)/SCALE
        L=L+1
        DELB=PHIB+(DELTA/EPSB)/EPSB
        QB=QB+ONE/DELB
        IF (ALPHA*QB*UTOL+ALPHA*QB+QAB.GE.GAMMA0) THEN
          ICODE=0
          GOTO 9
        ENDIF
        IF (DELTA.LE.TINY) THEN
          ICODE=4
          GOTO 9
        ENDIF
        IF (L.GT.LMAX) THEN
          ICODE=3
          GOTO 9
        ENDIF
        TEMP=ONE/EPSAB-DELAB*EPSAB/GAMMA
        PHIAB=PHIAB+ALPHA/(EPSAB*DELTA*EPSAB)+TEMP*GAMMA*TEMP
        TEMP=-ONE/DELTA
        CALL MEMSMA(ST,28,TEMP,26,26)
        CALL MEMTR(ST,26,2,.TRUE.)
        TEMP=ONE/SCALE
        CALL MEMSMD(ST,4,TEMP,2,GAMMA)
        GAM(L)=SQRT(GAMMA)*SCALE
        DELAB=PHIAB+(GAMMA/EPSAB)/EPSAB
        EPSAB=GAMMA/(EPSAB*DELAB)
        QAB=QAB+ONE/DELAB
        IF (ALPHA*QB*UTOL+ALPHA*QB+QAB.GE.GAMMA0) THEN
          ICODE=0
          GOTO 9
        ENDIF
        IF (GAMMA.LE.GAMMA0*TOL) THEN
          ICODE=1
          GOTO 9
        ENDIF
        TEMP=ONE/EPSB
        EPSB=DELTA/(EPSB*DELB)
        PHIB=PHIB+ALPHA/(EPSB*GAMMA*EPSB)+
     *            (ONE/EPSB-TEMP)*DELTA*(ONE/EPSB-TEMP)
        IF( FLAG ) THEN
          CALL MEMSMA(ST,26,COEFF(L),27,27)
          T0=T0+COEFF(L)
        ENDIF
      GOTO 1
* Rescale and exit
    9 CONTINUE
      QLO=SCALE*QB*SCALE
      QHI=SCALE*(GAMMA0-QAB)*SCALE/ALPHA
      LCODE=(ICODE.LE.2)
      IF( FLAG ) THEN
        CALL MEMTR(ST,27,2,.TRUE.)
        TEMP=T0/SCALE
        CALL MEMSMD(ST,4,TEMP,2,TEMP)
        CALL MESMUL(ST,2,SCALE,2)
        TEMP=SCALE*TEMP*SCALE
        LCODE = ABS(TEMP-QLO).LE.UTOL*QLO
      ENDIF
      END

      SUBROUTINE MEMCGY(ST,LMAX,ALPHA,UTOL,TOL,RATE,SUMMET,
     *                   HLOW,HHIGH,DIST,LCODEB,LCODEC)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Hidden-space conjugate gradient
* -------------------------------
*                  ------------------------------------
*   Find          |                      -1            |
*                 |    y = ( beta*I + A )  b           |
*   by maximising |                                    |
*                 |   QB = 2 y.b - y.(beta*I+A).y      |
*                  ------------------------------------
*           where  beta >= alpha  ensures that
*                   DIST = dimensionless length of y
*                        = sqrt( y.y / SUMMET )      <=   RATE
*   Return  DIST  and  H = [ HLOW , HHIGH ] = QB/2 = cross-entropy
*
*           where      A = Memtr.Memop
*                  Memtr = ( SQRT([metric]) * TrICF * Tropus * [acc] *)
*                  Memop = ( [acc] *  Opus  *  ICF  * SQRT([metric]) *)
*
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    LMAX    I     I      -         Maximum number of iterations
*    ALPHA   R     I      -         Regularisation parameter
*    UTOL    R     I      -         User's termination criterion
*    TOL     R     I      -         Arithmetic tolerance
*    RATE    R     I      -         Dimensionless distance limit
*    SUMMET  R     I      -         SUM([metric])
*    HLOW    R       O    -         Lower limit on cross-entropy
*    HHIGH   R       O    -         Upper limit on cross-entropy
*    DIST    R       O    -         Dimensionless distance moved
*    LCODEB  L       O    -         Status code for beta control
*    LCODEC  L       O    -         Status code for conj. gradients
* Areas:
*     <1>          I                SQRT([metric])
*     <2>           (O)             u (transform workspace)
*     <4>          I(O)             b input vector
*     <5>          I O              Hidden-space distribution
*    <22>         (I)               [acc]
*    <23>           (O)             yy = simulated result
*    <25>           (O)             uu (transform workspace)
*    <26>           (O)             gg = simulated gradient  (workspace)
*    <27>           (O)             vv = simulated conjugate for yy
*    <28>           (O)             hh = simulated conjugate (workspace)
*
* External calls:
*    MEMDOT      Dot product
*    MEMLB       Beta control
*    MEMOP       Apply ( [acc] *  Opus  *  ICF  * SQRT([metric]) *)
*    MEMSMA      Update area
*    MEMSMD      Update area with scalar product
*    MEMTR       Apply ( SQRT([metric]) * TrICF * Tropus * [acc] *)
*    MEMUPD      Update vector, keeping result within limits
*    MESMUL      Multiply area by scalar
*    MEVMUL      Area vector multiply
*    MEZERO      Initialise area to zero
*    UINFO       User's diagnostics handler
*
* The status codes are
*     LCODEB = ( beta == alpha , so that no distance penalty is invoked)
* and
*     LCODEC = ( ICODE <= 2 )
* where the internal codes are
*      ICODE = 0: User's termination criterion satisfied      (OK)
*      ICODE = 1: Convergence to within arithmetic TOLerance  (OK)
*      ICODE = 2: Zero results from TINY input vector b     (Warning)
*      ICODE = 3: Iteration limit LMAX exceeded             (Warning)
*      ICODE = 4: Probable error in transform routines       (Error)
*
* History:
*    MKC/JS          9 Feb 1991     First release
*    JS              8 Mar 1991     RATE included in beta control
*
* Notes:
*  (1) Maximise QB and  QAB = 2 z.A.b - z.A.(beta*I+A).z
*      by slaving to master conjugate gradient maximisation of
*                        QA = 2 x.b - x.A.x    (scalars only)
*  (2) Termination is either when QB has TOLerably converged, by
*         (1.+UTOL) * alpha * QB + QAB   passing   (length of b)**2
*      or when L passes
*         LMAX = maximum number of matrix applications
*      or when gamma underflows the internal arithmetic TOLerance.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL LCODEB,LCODEC
      INTEGER LMAX
      REAL ST(0:*)
      REAL  ALPHA,UTOL,TOL,RATE,SUMMET,HLOW,HHIGH,DIST

      INTEGER L,ICODE
      REAL GAMMA0,GAMMA,DELTA,TEMP,SCALE,BETA,YY
      REAL QB,EPSB,PHIB,DELB,V0,Y0,QAB,EPSAB,PHIAB,DELAB
      REAL ZERO,TINY,ONE,TWO
      PARAMETER (ZERO=0.0D0,TINY=1.0D-30,ONE=1.0D0,TWO=2.0D0)

      CALL UINFO('   MemCGY',1)
* Initialise
      L=0
      EPSB=ONE
      QB=ZERO
      CALL MEZERO(ST,23)
      CALL MEZERO(ST,27)
      QAB=ZERO

* SCALE arbitrarily (make input vector unit) to reduce overflow risk 
      CALL MEVMUL(ST,4,1,4)
      CALL MEMDOT(ST,4,4,TEMP)
      IF (TEMP.LT.TINY) THEN
        GAMMA0=ZERO
        SCALE=ONE
        ICODE=2
        GOTO 9
      ENDIF
      SCALE=SQRT(TEMP)
      GAMMA0=(TEMP/SCALE)/SCALE
      GAMMA=GAMMA0
      
* Beta control
      CALL MEMLB(ALPHA,RATE,TEMP/(RATE*SUMMET*RATE),BETA,LCODEB)

      CALL MEMOP(ST,4,28,.TRUE.)
      TEMP=ONE/SCALE
      CALL MESMUL(ST,28,TEMP,28)
      CALL MEMDOT(ST,28,28,DELTA)
      CALL MEZERO(ST,26)
      PHIB=BETA
      EPSB=ONE
      DELB=PHIB+DELTA
      QB=ONE/DELB
      V0=ONE
      Y0=V0/DELB
      EPSAB=ONE
      PHIAB=ZERO
      DELAB=ZERO
      QAB=ZERO
* Loop
    1 L=L+1
        IF (BETA*QB*UTOL+BETA*QB+QAB.GE.GAMMA0) THEN
          ICODE=0
          GOTO 9
        ENDIF
        IF (DELTA.LE.TINY) THEN
          ICODE=4
          GOTO 9
        ENDIF
        IF (L.GT.LMAX) THEN
          ICODE=3
          GOTO 9
        ENDIF
        TEMP=ONE/EPSAB-DELAB*EPSAB/GAMMA
        PHIAB=PHIAB+BETA/(EPSAB*DELTA*EPSAB)+TEMP*GAMMA*TEMP
        TEMP=-ONE/DELTA
        CALL MEMSMA(ST,28,TEMP,26,26)
        TEMP=ONE/SCALE
        CALL MEMTR(ST,26,2,.TRUE.)
        CALL MEMSMD(ST,4,TEMP,2,GAMMA)
        DELAB=PHIAB+(GAMMA/EPSAB)/EPSAB
        EPSAB=GAMMA/(EPSAB*DELAB)
        QAB=QAB+ONE/DELAB
        IF (BETA*QB*UTOL+BETA*QB+QAB.GE.GAMMA0) THEN
          ICODE=0
          GOTO 9
        ENDIF
        IF (GAMMA.LE.GAMMA0*TOL) THEN
          ICODE=1
          GOTO 9
        ENDIF
        TEMP=ONE/EPSB
        EPSB=DELTA/(EPSB*DELB)
        PHIB=PHIB+BETA/(EPSB*GAMMA*EPSB)+
     *            (ONE/EPSB-TEMP)*DELTA*(ONE/EPSB-TEMP)
        TEMP=ONE/GAMMA
        CALL MEMOP(ST,2,25,.TRUE.)
        CALL MEMSMD(ST,25,TEMP,28,DELTA)
        DELB=PHIB+(DELTA/EPSB)/EPSB
        QB=QB+ONE/DELB
        TEMP=ONE/(EPSB*GAMMA)
        V0=V0+TEMP
        CALL MEMSMA(ST,26,TEMP,27,27)
        TEMP=ONE/DELB
        Y0=Y0+V0*TEMP
        CALL MEMSMA(ST,27,TEMP,23,23)
      GOTO 1
* Rescale and exit
    9 CONTINUE
      LCODEC=(ICODE.LE.2)
      CALL MESMUL(ST,23,SCALE,23)
      CALL MEMTR(ST,23,2,.TRUE.)
      CALL MEMSMA(ST,4,Y0,2,4)
      CALL MEMDOT(ST,4,4,YY)
      HLOW=SCALE*QB*SCALE/TWO
      HHIGH=SCALE*(GAMMA0-QAB)*SCALE/(TWO*BETA)
      DIST=SQRT(YY/SUMMET)
      CALL MEVMUL(ST,4,1,4)
      CALL MEMUPD(ST,5,4,5)
      END

      SUBROUTINE MEMCGP(ST,LMAX,ALPHA,UTOL,TOL, GAM,DEL,M,N,LCODEC)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Data-space conjugate gradient
* -----------------------------
*                  -------------------------------------
*   Find          |                       -1            |
*                 |    z = ( alpha*I + A )  b           |
*   by maximising |                                     |
*                 |  QAB = 2 z.A.b - z.A.(alpha*I+A).z  |
*                  -------------------------------------
*
*                      ---------------------------
*   and set up        |  GAMma and DELta scalars  |    defining A
*                      ---------------------------
*
*           where      A = Memop.Memtr
*                  Memtr = ( SQRT([metric]) * TrICF * Tropus * [acc] *)
*                  Memop = ( [acc] *  Opus  *  ICF  * SQRT([metric]) *)
*
*    In terms of normalised (and orthogonal) gradient vectors
*      --------       --------  -------- 
*     |        |     |        ||\       |
*     |   ^    |     |        ||  1/GAM |
*     |   g    |  =  |   g    ||       \|
*     |        |     |        | -------- 
*     |        |     |        |
*      --------       -------- 
*    the matrix A becomes represented by the factorised tridiagonal form
*      --------       --------  --------  --------  --------  --------
*     | ^t   ^ |     |\       || 1      ||\       || 1 -1   ||\       |
*     | g .A.g |  =  |  GAM   ||-1  1   || DEL**2 ||    1 -1||  GAM   |
*     |        |     |       \||   -1  1||       \||       1||       \|
*      --------       --------  --------  --------  --------  --------
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    LMAX    I     I      -         Maximum number of iterations
*    ALPHA   R     I      -         Regularisation parameter
*    UTOL    R     I      -         User's termination criterion
*    TOL     R     I      -         Arithmetic tolerance
*    GAM     R       O   0:LMAX     Gamma from each iteration
*    DEL     R       O   0:LMAX     Delta from each iteration
*    M       I       O    -         GAM(0..M) set on exit
*    N       I       O    -         DEL(0..N) set on exit
*    LCODEC  L       O    -         Status code for conj. gradients
* Areas:
*     <1>          I                SQRT([metric])
*     <2>           (O)             u (transform workspace)
*    <22>         (I)               [acc]
*    <25>           (O)             u (transform workspace)
*    <26>          I(O)             b input vector (& g gradient for x)
*    <28>           (O)             h conjugate vector for x
*
* External calls:
*    MEMDOT      Dot product
*    MEMOP       Apply ( [acc] *  Opus  *  ICF  * SQRT([metric]) *)
*    MEMSMA      Update area
*    MEMSMD      Update area with scalar product
*    MEMTR       Apply ( SQRT([metric]) * TrICF * Tropus * [acc] *)
*    MESMUL      Multiply area by scalar
*    UINFO       User's diagnostics handler
*
* The status code is
*     LCODEC = ( ICODE <= 2 )
* where the internal codes are
*      ICODE = 0: User's termination criterion satisfied      (OK)
*      ICODE = 1: Convergence to within arithmetic TOLerance  (OK)
*      ICODE = 2: Zero results from TINY input vector b     (Warning)
*      ICODE = 3: Iteration limit LMAX exceeded             (Warning)
*      ICODE = 4: Probable error in transform routines       (Error)
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*  (1) Maximise QAB and  QB = 2 y.b - y.(alpha*I+A).y  and set GAM,DEL
*      by slaving them to master conjugate gradient maximisation of
*                        QA = 2 x.b - x.A.x    (scalars only)
*  (2) Termination is either when QAB has TOLerably converged, by
*         (1.+UTOL) * QAB + alpha * QB   passing   (length of b)**2
*      or when L passes
*         LMAX = maximum number of matrix applications
*      or when gamma underflows the internal arithmetic TOLerance.
*  (3) On exit, either M = N or M = N+1 .
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL LCODEC
      INTEGER LMAX,M,N
      REAL ST(0:*),GAM(0:LMAX),DEL(0:LMAX)
      REAL ALPHA,UTOL,TOL

      INTEGER L,ICODE
      REAL GAMMA0,GAMMA,DELTA,TEMP,SCALE
      REAL QB,EPSB,PHIB,DELB,QAB,EPSAB,PHIAB,DELAB
      REAL ZERO,TINY,ONE,TWO
      PARAMETER (ZERO=0.0D0,TINY=1.0D-30,ONE=1.0D0,TWO=2.0D0)

      CALL UINFO('   MemCGP',1)
* Initialise
      L=0
      M=0
      N=0
      EPSB=ONE
      QB=ZERO
      QAB=ZERO

* SCALE arbitrarily (make input vector unit) to reduce overflow risk 
* Internal variables rescale proportionally to powers of SCALE:
* U +1,               GAMMA -2, DELTA +2,     G -1,   H +1
*              QB -2,   EPSB 0,  PHIB +2,  DELB +2
* W +1, Z -1, QAB -2,  EPSAB-2, PHIAB +2, DELAB +2
      CALL MEMDOT(ST,26,26,TEMP)
      IF (TEMP.LT.TINY) THEN
        ICODE=2
        SCALE=ONE
        GAMMA0=ZERO
        GAM(0)=ZERO
        DEL(0)=ZERO
        GOTO 9
      ENDIF
      SCALE=SQRT(TEMP)
      GAMMA0=(TEMP/SCALE)/SCALE
      GAMMA=GAMMA0
      GAM(0)=SQRT(GAMMA)*SCALE
      TEMP=ONE/SCALE
      CALL MESMUL(ST,26,TEMP,26)
      EPSAB=GAMMA
      PHIAB=ZERO
      DELAB=ZERO
      TEMP=ONE/GAMMA
      CALL MESMUL(ST,26,TEMP,28)
      CALL MEMTR(ST,28,2,.TRUE.)
      CALL MEMDOT(ST,2,2,DELTA)
      CALL MEMOP(ST,2,25,.TRUE.)
      DEL(0)=SQRT(DELTA)/SCALE
      PHIB=ALPHA/GAMMA
* Loop
    1 L=L+1
        DELB=PHIB+(DELTA/EPSB)/EPSB
        QB=QB+ONE/DELB
        IF (UTOL*QAB+QAB+ALPHA*QB.GE.GAMMA0) THEN
          ICODE=0
          GOTO 9
        ENDIF
        IF (DELTA.LE.TINY) THEN
          ICODE=4
          GOTO 9
        ENDIF
        IF (L.GT.LMAX) THEN
          ICODE=3
          GOTO 9
        ENDIF
        IF (L.GT.1) CALL MEMOP(ST,2,25,.TRUE.)
        TEMP=ONE/EPSAB-DELAB*EPSAB/GAMMA
        PHIAB=PHIAB+ALPHA/(EPSAB*DELTA*EPSAB)+TEMP*GAMMA*TEMP
        TEMP=-ONE/DELTA
        CALL MEMSMD(ST,25,TEMP,26,GAMMA)
        DELAB=PHIAB+(GAMMA/EPSAB)/EPSAB
        M=L
        GAM(M)=SQRT(GAMMA)*SCALE
        EPSAB=GAMMA/(EPSAB*DELAB)
        QAB=QAB+ONE/DELAB
        IF (UTOL*QAB+QAB+ALPHA*QB.GE.GAMMA0) THEN
          ICODE=0
          GOTO 9
        ENDIF
        IF (GAMMA.LE.GAMMA0*TOL) THEN
          ICODE=1
          GOTO 9
        ENDIF
        TEMP=ONE/GAMMA
        CALL MEMSMA(ST,26,TEMP,28,28)
        TEMP=ONE/EPSB
        EPSB=DELTA/(EPSB*DELB)
        PHIB=PHIB+ALPHA/(EPSB*GAMMA*EPSB)+
     *            (ONE/EPSB-TEMP)*DELTA*(ONE/EPSB-TEMP)
        CALL MEMTR(ST,28,2,.TRUE.)
        CALL MEMDOT(ST,2,2,DELTA)
        N=L
        DEL(N)=SQRT(DELTA)/SCALE
      GOTO 1
* Rescale and exit
    9 CONTINUE
      LCODEC=(ICODE.LE.2)
      END

      SUBROUTINE MEMDET(TOL,VAL,VEC,OFF,N,M,IFLAG)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                                                                
*         Eigenstructure of symmetric tridiagonal matrix
*                t            t
*         A = L L   or   A = L  L    where L is lower bidiagonal.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    TOL     R     I       -         Arithmetic tolerance
*    VAL     R     I      0:N-1      Diagonal elements    0..N-1 of L
*                    O               Eigenvalues of A
*    VEC     R       O  0:M-1,0:N-1  VEC(.,K) is 1st part of evec K of A
*    OFF     R     I      0:N-1      Subdiagonal elements 1..N-1 of L
*                    O               (workspace)
*    N       I     I       -         Dimension
*    M       I     I       -         Number of evector components needed
*    IFLAG   I     I       -         +ve is A = L Lt , -ve is A = Lt L
*
* External calls:
*     -
*
* History:
*    MKC/JS          1 Dec 1990     First release
*    MKC/JS         17 Feb 1991     N <= 1 return fixed
*
* Notes:
*    (1)  VAL must be strictly positive, except that VAL(N-1) may be 0.
*         OFF must be strictly negative, except that OFF(N-1) may be 0.
*         The internal signs in this routine are specifically chosen for
*         this case, to ensure that the rotations preserve this form.
*    (2)  The final run, after the input angle has been found by binary
*         chop, must repeat the final simulation exactly, which is why
*         the same code is re-used.
*    (3)  DO NOT try to edit this - almost all changes will cause at
*         least occasional damage!
*    (4)  The eigenvalues return in decreasing order, unless an early
*         one is TOLerably small - see loop at instruction 3.
*    (5)  Running time proportional to  N * N *( M + log[2](1/TOL) )
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER IFLAG,M,N
      REAL TOL,VAL(0:N-1),VEC(0:M-1,0:N-1),OFF(0:N-1)

      LOGICAL FINISH,EXIT
      INTEGER I,J,K,L
* Register......
      REAL C,S,R,X,D0,D1,E0,E1
* .....variables
      REAL CLO,CEQ,CHI,DC,DCOLD,SLO,SEQ,SHI,DS,DSOLD
      REAL ZERO,ONE,TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)

      OFF(0)=ZERO
* Set vectors VEC to be rotated into eigenvector coords.
      DO 2 J=0,M-1
        DO 1 I=0,N-1
          VEC(J,I)=ZERO
    1   CONTINUE
        VEC(J,J)=ONE
    2 CONTINUE
* Exceptional return
      IF(N.LE.1) GOTO 9
* Iterate each successive OFF(K) ( K = N-1,N-2,...) towards 0.
* Ignore any early OFF(L) ( L = 0,1,2,...) which are TOLerably small.
      K=N-1
      L=0
    3 IF( OFF(L+1)+TOL*(VAL(L)+VAL(L+1)) .GT. ZERO ) THEN
        L=L+1
        IF(L.LT.K-1) GOTO 3
      ENDIF
* Initialise binary chop, using both Cos and Sin for full accuracy.
    4 FINISH=.FALSE.
        CLO=ONE
        SLO=ZERO
        CEQ=ONE/TWO
        SEQ=-ONE/TWO
        CHI=ZERO
        SHI=-ONE
        DC=ONE
        DS=ONE
* Chop the initial rotation angle (Q) between 0 and -pi/2,
* until the last off-diagonal element vanishes.
    5   EXIT=FINISH
          R=SQRT(CEQ**2+SEQ**2)
          C=CEQ/R
          S=SEQ/R
*     Enter with lower bidiagonal matrix L and rotation Q (C=cos,S=sin).
*     Apply Q to L,
*              -----------   -----------
*             | C -S      | |d0         |
*             | S  C      | |e1 d1      |   (e.g. dimension n=4)
*      Q L =  |       1   | |   e2 d2   |   (d in VAL, e in OFF)
*             |          1| |      e3 d3|
*              -----------   ----------- 
*     with such other right and left rotations as are needed to preserve
*     the lower bidiagonal form of A.  Successively (left to right):-
*  -------   -------   -------   -------   -------   -------   ------- 
* |d      | |d X    | |d      | |d      | |d      | |d      | |d      |
* |e d    | |e d    | |e d    | |e d X  | |e d    | |e d    | |e d    |
* |  e d  | |  e d  | |X e d  | |  e d  | |  e d  | |  e d X| |  e d  |
* |    e d| |    e d| |    e d| |    e d| |  X e d| |    e d| |    e d|
*  -------   -------   -------   -------   -------   -------   ------- 
          I=L
          D1=VAL(L)
          E1=OFF(L+1)
    6     I=I+1
            D0=D1
            D1=VAL(I)
            R =C*E1-S*D0
            D0=S*E1+C*D0
            E1=R
            X =S*D1
            D1=C*D1
            IF( EXIT ) THEN
              IF (IFLAG.GE.0) THEN
                DO 7 J=0,M-1
                  R         =C*VEC(J,I)-S*VEC(J,I-1)
                  VEC(J,I-1)=S*VEC(J,I)+C*VEC(J,I-1)
                  VEC(J,I)  =R
    7           CONTINUE
              ENDIF
            ENDIF
*     X is now overflow element on upper-right
*                ----
*     Apply  i-1|C -S|    on right to eliminate X
*             i |S  C|
*                ----
            R=ABS(D0)+ABS(X)
            R=R*SQRT((D0/R)**2+(X/R)**2)
            C=D0/R  
            S=X/R
            D0=S*X+C*D0
            R =C*D1-S*E1
            E1=S*D1+C*E1
            D1=R 
            IF( EXIT ) THEN
              VAL(I-1)=D0
              IF (IFLAG.LT.0) THEN
                DO 8 J=0,M-1
                  R         =C*VEC(J,I)-S*VEC(J,I-1)
                  VEC(J,I-1)=S*VEC(J,I)+C*VEC(J,I-1)
                  VEC(J,I)  =R
    8           CONTINUE
              ENDIF
            ENDIF
* The zeros of the last off-diagonal terms OFF(I) at successive
* dimensions I interleave.  Hence, if any OFF(I) would change sign,
* then reduce the initial angle to approach the first zero of OFF(K).
            IF(E1.LT.ZERO) THEN
              IF (I.LT.K) THEN
                E0=E1
                E1=OFF(I+1)
                X =S*E1
                E1=C*E1
*     X is now overflow element on lower-left
*                ----
*     Apply   i | C S|    on left to eliminate X
*            i+1|-S C|
*                ----
                R=ABS(E0)+ABS(X)
                R=-R*SQRT((E0/R)**2+(X/R)**2)
                C=E0/R
                S=X/R
                IF( EXIT ) OFF(I)=R
                GOTO 6
              ENDIF
              CLO=CEQ
              SLO=SEQ
* Exit if last off-diagonal element OFF(K) is TOLerably small.
              FINISH = FINISH .OR.( E1+TOL*(D0+D1) .GT. ZERO )
            ELSE
              CHI=CEQ
              SHI=SEQ
            ENDIF
            DCOLD=DC
            DSOLD=DS
            DC=CLO-CHI
            DS=SLO-SHI
            FINISH = FINISH .OR. ( DCOLD.LE.DC .AND. DSOLD.LE.DS )
* Repeat the final successful (low) iterate to get output
          IF( FINISH ) THEN
            CEQ=CLO
            SEQ=SLO
          ELSE
            CEQ=CHI+DC/TWO
            SEQ=SLO-DS/TWO
          ENDIF
        IF(.NOT.EXIT) GOTO 5
        VAL(I)=D1
* Exit having re-used the rotation angle which kills the last off-
* diagonal term OFF(K), effectively reducing the dimension by 1.
* E1 contains new OFF(K), which should be small negative diagnostic.
        K=K-1
      IF(K.GT.L) GOTO 4
    9 CONTINUE
* Square up eigenvalues on exit
      DO 10 I=0,N-1
        VAL(I)=VAL(I)**2
   10 CONTINUE
      END

** control.for

      SUBROUTINE MEMLA0(RATE,SUMMET,GRADL,OMEGA,VAR, ALPHA,LCODE)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*          Alpha control if MEMRUN=1
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    RATE    R     I      -         Dimensionless distance limit
*    SUMMET  R     I      -         SUM([metric])
*    GRADL   R     I      -         gradL
*    OMEGA   R     I      -         value of stopping criterion
*    VAR     R     I      -         Intrinsic variance of (OMEGA)
*    ALPHA   R       O    -         Output alpha
*    LCODE   L       O    -         Alpha control diagnostic
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
*
* External calls:
*    MEMLAW      Insert new entry into omega(alpha) table
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL LCODE
      REAL RATE,SUMMET,GRADL,OMEGA,VAR,ALPHA

      REAL EPS,ONE,BIG
      PARAMETER (EPS=1.0D-20,ONE=1.0D0,BIG=1.0D20)
      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      CALL UINFO('   MemLA0',1)
      CALL MEMLAW(BIG,OMEGA,VAR,.TRUE.)
      IF (OMEGA.LT.ONE) THEN
        ALPHA=(GRADL/SQRT(SUMMET))/RATE
        LCODE=.FALSE.
      ELSE
        ALPHA=BIG
        LCODE=.TRUE.
      ENDIF
      WRITE (INFO,'(''     alpha  '',1PE16.7)') ALPHA
      CALL UINFO(INFO,20)
      END

      SUBROUTINE MEMLA(ALPHA,RATE,SUMMET,OMEGA,TOL,VAR,GRADL,INIT,
     *                 LCODEB,LCODEA,LCODET)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*          Alpha control, iterating towards Omega rising to 1.00
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ALPHA   R     I O    -         Regularisation parameter
*    RATE    R     I      -         Dimensionless distance limit
*    SUMMET  R     I      -         SUM([metric])
*    OMEGA   R     I      -         value of stopping criterion
*    TOL     R     I      -         Arithmetic tolerance
*    VAR     R     I      -         Intrinsic variance of (OMEGA)
*    GRADL   R     I      -         gradL
*    INIT    L     I      -         Flag for initialisation of table
*    LCODEB  L     I      -         Beta control diagnostic
*    LCODEA  L       O    -         Alpha control diagnostic
*    LCODET  L     I      -         TEST diagnostic
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
*
* External calls:
*    MEMLAR      Interpolate omega(alpha) table
*    MEMLAW      Insert new entry into omega(alpha) table
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL INIT,LCODEB,LCODEA,LCODET
      REAL ALPHA,RATE,SUMMET,OMEGA,TOL,VAR,GRADL

      REAL R,YNEW,SIGMA,ALF1,ALF2,Y1,Y2
      REAL ZERO,ONE,TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      CALL UINFO('   Control',1)
      CALL MEMLAW(ALPHA,OMEGA,VAR,INIT)
      WRITE (INFO,'(''     rate   '',1PE16.7)') RATE
      CALL UINFO(INFO,20)
      WRITE (INFO,'(''     summet '',1PE16.7)') SUMMET
      CALL UINFO(INFO,20)
      WRITE (INFO,'(''     gradL  '',1PE16.7)') GRADL
      CALL UINFO(INFO,20)
      WRITE (INFO,'(''     alpha  '',1PE16.7)') ALPHA
      CALL UINFO(INFO,20)
      IF (LCODEB.AND.LCODET) THEN
        CALL MEMLAR(ALPHA,YNEW,SIGMA)
        IF(ABS(YNEW).LT.SIGMA) THEN
          LCODEA=.TRUE.
        ELSE
          LCODEA=.FALSE.
          R=ONE+ALPHA*SQRT(SUMMET)*RATE/GRADL
          ALF1=ALPHA*R
          ALF2=ALPHA/R
          ALF1=MAX(ALF1,ALPHA)
          ALF2=MIN(ALF2,ALPHA)
          CALL MEMLAR(ALF1,Y1,SIGMA)
          IF (Y1.GT.ZERO) THEN
            ALPHA=ALF1
          ELSE
            CALL MEMLAR(ALF2,Y2,SIGMA)
            IF (Y2.LE.ZERO) THEN
              ALPHA=ALF2
            ELSE
              R=ONE
    1           R=R/TWO
                ALPHA=(ALF1+ALF2)/TWO
                CALL MEMLAR(ALPHA,YNEW,SIGMA)
                IF (YNEW.LE.ZERO) THEN
                  ALF1=ALPHA
                  Y1=YNEW
                ELSE
                  ALF2=ALPHA
                  Y2=YNEW
                ENDIF
              IF (R.GT.TOL) GOTO 1
              IF (Y1-Y2.NE.ZERO) THEN
                R=Y1/(Y1-Y2)
                ALPHA=ALF1+(ALF2-ALF1)*R
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ELSE
        LCODEA=.TRUE.
      ENDIF
      WRITE (INFO,'(''     alpha  '',1PE16.7)') ALPHA
      CALL UINFO(INFO,20)
      END

      SUBROUTINE MEMLAW(ALPHA,OMEGA,VAR,INIT)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*         Insert new entry into table of OMEGA(alpha) with variances,
*         discarding least significant earlier entry if table is full.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ALPHA   R     I      -         New alpha
*    OMEGA   R     I      -         New omega(alpha)
*    VAR     R     I      -         Intrinsic variance of log(OMEGA)
*    INIT    L     I      -         Initialisation flag
* Local externals:
*    NSIZE   I     I      -         Dimensioned length of table
*    NTABLE  I     I O    -         Current length of table
*    XTABLE  R     I O   0:NSIZE-1  Values of log(alpha)
*    YTABLE  R     I O   0:NSIZE-1  Values of log(omega)
*    VTABLE  R     I O   0:NSIZE-1  Values of intrinsic variance
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
*
* External calls:
*    UINFO       User's diagnostics handler
*
* History:
*    John Skilling   7 Dec 1988     MEMSYS3 version
*    MKC/JS          1 Dec 1990     First release
*    JS              8 Mar 1991     Overflow risk reduced
*
* Notes:
*   (1) Deletion algorithm is ad hoc.
*   (2) "Development code" can be used to test control of the tables.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL INIT
      REAL ALPHA,OMEGA,VAR
      INCLUDE 'memsys.inc'

      INTEGER NTABLE
      REAL XTABLE(0:NSIZE-1),YTABLE(0:NSIZE-1),
     *     VTABLE(0:NSIZE-1)
      EQUIVALENCE (NTABLE,ICOM(KNTABL))
      EQUIVALENCE (XTABLE,RCOM(KXTABL)),(YTABLE,RCOM(KYTABL))
      EQUIVALENCE (VTABLE,RCOM(KVTABL))

      INTEGER I,J
      REAL D,S,X,XNEW,YNEW,VNEW,EPS
      PARAMETER (EPS=1.0D-20)
      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      XNEW=LOG(ALPHA)
      YNEW=LOG(MAX(OMEGA,EPS))
      VNEW=MAX(VAR,EPS)
      IF (INIT) NTABLE=0
      DO 1 I=0,NTABLE-1
        IF(XNEW.EQ.XTABLE(I)) THEN
          S=VNEW/(VNEW+VTABLE(I))
          X=VTABLE(I)/(VNEW+VTABLE(I))
          YTABLE(I)=YNEW*X+YTABLE(I)*S
          VTABLE(I)=VNEW*X
          RETURN
        ENDIF
    1 CONTINUE
      IF(NTABLE.GE.NSIZE) THEN
* Delete worst alpha, preserving order
        S=MAX(VAR,EPS)
        J=-1
        DO 2 I=0,NTABLE-1
          X=LOG(ALPHA)-XTABLE(I)
          D=VTABLE(I)+X*X*X*X
          IF( D.GT.S ) THEN
            S=D
            J=I
          ENDIF
    2   CONTINUE
        IF(J.EQ.-1) RETURN
        NTABLE=NTABLE-1
        DO 3 I=J,NTABLE-1
          XTABLE(I)=XTABLE(I+1)
          YTABLE(I)=YTABLE(I+1)
          VTABLE(I)=VTABLE(I+1)
    3   CONTINUE
      ENDIF
* Insert alpha
      J=0
      DO 4 I=0,NTABLE-1
        IF(XNEW.GT.XTABLE(I)) J=I+1
    4 CONTINUE
      DO 5 I=NTABLE-1,J,-1
        XTABLE(I+1)=XTABLE(I)
        YTABLE(I+1)=YTABLE(I)
        VTABLE(I+1)=VTABLE(I)
    5 CONTINUE
      XTABLE(J)=XNEW
      YTABLE(J)=YNEW
      VTABLE(J)=VNEW
      NTABLE=NTABLE+1
      CALL UINFO('   Control on  log alpha'//
     *           '       log omega       var omega',20)

      DO 6 I=0,NTABLE-1
        WRITE(INFO,'('' '',I4,1PE20.7,2E16.7)')
     *                      I,XTABLE(I),YTABLE(I),VTABLE(I)
        CALL UINFO(INFO,20)
    6 CONTINUE
      END

      SUBROUTINE MEMLAR(ALPHA,YVAL,SIGMA)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*         Interpolate table to estimate log(omega(alpha)) and its error
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ALPHA   R     I      -         alpha
*    YVAL    R       O    -         Estimate of log(omega(alpha))
*    SIGMA   R       O    -         Standard deviation of estimate
* Local externals:
*    NSIZE   I     I      -         Dimensioned length of table
*    NTABLE  I     I      -         Current length of table
*    XTABLE  R     I     0:NSIZE-1  Values of log(alpha)
*    YTABLE  R     I     0:NSIZE-1  Values of log(omega)
*    VTABLE  R     I     0:NSIZE-1  Values of intrinsic variance
*    TOL     R     I      -         Arithmetic tolerance
*
* External calls:
*     -
*
* History:
*    John Skilling   7 Dec 1988     MEMSYS3 version
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ALPHA,YVAL,SIGMA
      INCLUDE 'memsys.inc'

      INTEGER NTABLE
      REAL XTABLE(0:NSIZE-1),YTABLE(0:NSIZE-1),
     *     VTABLE(0:NSIZE-1),TOL
      EQUIVALENCE (NTABLE,ICOM(KNTABL))
      EQUIVALENCE (XTABLE,RCOM(KXTABL)),(YTABLE,RCOM(KYTABL))
      EQUIVALENCE (VTABLE,RCOM(KVTABL)),(TOL,   RCOM(KTOL  ))

      INTEGER I
      REAL X,Y,WT,W,WX,WXX,WY,WXY,XBAR,YBAR,DIV,XXX,YYY
      REAL ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)

      DIV(XXX,YYY) = XXX/MAX(YYY,TOL)

      W=ZERO
      WX=ZERO
      WXX=ZERO
      WY=ZERO
      WXY=ZERO
      IF(NTABLE.EQ.1) THEN
        YVAL=YTABLE(0)
        SIGMA=ZERO
      ELSE
        DO 1 I=0,NTABLE-1
          X=LOG(ALPHA)-XTABLE(I)
          Y=YTABLE(I)
          WT=DIV(ONE,VTABLE(I)+X*X*X*X)
          W=W+WT
          WX=WX+WT*X
          WY=WY+WT*Y
    1   CONTINUE
        XBAR=WX/W
        YBAR=WY/W
        DO 2 I=0,NTABLE-1
          X=LOG(ALPHA)-XTABLE(I)
          Y=YTABLE(I)
          WT=DIV(ONE,VTABLE(I)+X*X*X*X)
          X=X-XBAR
          Y=Y-YBAR
          WXX=WXX+WT*X*X
          WXY=WXY+WT*X*Y
    2   CONTINUE
        YVAL=YBAR-XBAR*DIV(WXY,WXX)
        SIGMA=ONE/SQRT(W)
      ENDIF
      END

      SUBROUTINE MEMLB(ALPHA,RATE,D0,BETA,LCODE)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*          Beta control.
*
*          Given g.g , choose beta >= alpha/rate such that
*              -2
*           g.B  .g  <  SUMMET*RATE**2  where  B = beta*I + A
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ALPHA   R     I      -         Regularisation parameter
*    RATE    R     I      -         Dimensionless distance limit
*    D0      R     I      -         g.g / (SUMMET*RATE**2)
*    BETA    R       O    -         Revised alpha
*    LCODE   L       O    -         Return code
*                                   .FALSE. = distance penalty invoked
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
*
* External calls:
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*    JS              8 Mar 1991     RATE included in beta control
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL LCODE
      REAL ALPHA,RATE,D0,BETA

      REAL BMIN
      REAL ONE
      PARAMETER (ONE=1.0D0)
      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      CALL UINFO('   MemLB',1)
      BMIN=SQRT(D0)
*      IF (.FALSE.) THEN
      IF (BMIN.GT.ALPHA/MIN(ONE,RATE)) THEN
        BETA=BMIN
        LCODE=.FALSE.
      ELSE
        BETA=ALPHA/MIN(ONE,RATE)
        LCODE=.TRUE.
      ENDIF
      WRITE (INFO,'(''      beta  '',1PE16.7)') BETA
      CALL UINFO(INFO,20)
      END

** utility.for

      SUBROUTINE MCLEAR
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*   Clear read/write flags
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Local externals:
*    KF      I     I O    40        Read/Write flags, set to -1
*    LENGTH    O     O    -         Old IOFF = length last area accessed
*    IOFF    I O     O    -         Offset, set to zero
*    MORE    I       O    -         Any more elements available?
*
* External calls:
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INCLUDE 'memsys.inc'

      INTEGER KF(40),IOFF,LENGTH,MORE
      EQUIVALENCE (IOFF,  ICOM(KIOFF )),(KF,    ICOM(KKF   ))
      EQUIVALENCE (LENGTH,ICOM(KLNGTH)),(MORE,  ICOM(KMORE ))

      INTEGER J
      CHARACTER CHART(-1:3)*1,A(4)*10,B(40)*1
      EQUIVALENCE(A,B)
      SAVE CHART
      DATA CHART/ '.' , ':' , 'R' , 'W' , 'U' /

      DO 1 J=1,40
        B(J)=CHART(KF(J))
        KF(J)=-1
    1 CONTINUE
      CALL UINFO('                     '//
     *           A(1)//'  '//A(2)//'  '//A(3)//'  '//A(4) , 2 )
      LENGTH=IOFF
      IOFF=0
      MORE=0
      END

      SUBROUTINE MFETCH(ST,J,IFLAG, KCORE,LCORE)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*             Prepare to read or write a block of area J.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I      *         Vector arithmetic workspace
*    J       I     I      -         Area number
*    IFLAG   I     I      -         Flag to supply contents of block
*                                   0=':'    1='R'    2='W'    3='U'
*    KCORE   I       O    -         Core address for use by VLIB
*    LCORE   I       O    -         Block length for use by VLIB
* Local externals:
*    KC(J)   I       O    -         Core address
*    KF(J)   I     I O    -         Read/Write flag
*    IOFF    I     I      -         Relative address within area
*    LENGTH  I       O    -         Block length
*    MORE    I     I(O)   -         Any more elements available?
*
* External calls:
*    MEFLAG      Set read/write/update flag
*    VFETCH      Prepare block and return address in ST
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*    (1) A chain of MFETCH calls is designed to be terminated by MWHILE
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER J,IFLAG,KCORE,LCORE
      REAL ST(0:*)
      INCLUDE 'memsys.inc'

      INTEGER KC(40),KF(40),IOFF,LENGTH,MORE
      EQUIVALENCE (KC,    ICOM(KKC   )),(KF,    ICOM(KKF   ))
      EQUIVALENCE (IOFF,  ICOM(KIOFF )),(LENGTH,ICOM(KLNGTH))
      EQUIVALENCE (MORE,  ICOM(KMORE ))

      LOGICAL LREAD(0:3)
      SAVE LREAD
      DATA LREAD/.FALSE.,.TRUE.,.FALSE.,.TRUE./

      IF (J.LT.1 .OR. J.GT.40) STOP ' Illegal area number'
      IF (KF(J).LT.0)
     *   CALL VFETCH(ST,J,IOFF,LREAD(IFLAG), KC(J),LENGTH,MORE)
      LCORE=LENGTH
      KCORE=KC(J)
      CALL MEFLAG(J,IFLAG)
      END

      LOGICAL FUNCTION MWHILE(ST)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*   Deal with current block and flag while more blocks remain.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
* Local externals:
*    KC      I     I      40        Core address
*    KF      I     I O    40        Read/Write flags
*    IOFF    I     I      -         Relative address within area
*    LENGTH  I     I      -         Block length
*    MORE    I    (I)O    -         Any more elements available?
*
* External calls:
*    MCLEAR      Clear read/write flags
*    VSTORE      Signal end of block use
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*    (1) MWHILE is designed to terminate a chain of MFETCH calls
*    (2) After each MWHILE, IOFF = total number of vector elements used
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      INCLUDE 'memsys.inc'

      INTEGER KC(40),KF(40),IOFF,LENGTH,MORE
      EQUIVALENCE (KC,    ICOM(KKC   )),(KF,    ICOM(KKF   ))
      EQUIVALENCE (IOFF,  ICOM(KIOFF )),(LENGTH,ICOM(KLNGTH))
      EQUIVALENCE (MORE,  ICOM(KMORE ))

      INTEGER J
      LOGICAL LWRITE(0:3)
      SAVE LWRITE
      DATA LWRITE/.FALSE.,.FALSE.,.TRUE.,.TRUE./

      DO 1 J=1,40
        IF (KF(J).GE.0)
     *      CALL VSTORE(ST,J,IOFF,LWRITE(KF(J)),KC(J),LENGTH, MORE)
    1 CONTINUE
      IOFF=IOFF+LENGTH
      MWHILE=(MORE.NE.0)
      IF (.NOT.MWHILE) CALL MCLEAR
      DO 2 J=1,40
        KF(J)=-1
    2 CONTINUE
      END

      SUBROUTINE MEFLAG(J,IFLAG)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*           Reset the Read/Write flag for area J
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    J       I     I      -         Area number
*    IFLAG   I     I      -         Read or write indicator
*                                   0='.'    1='R'    2='W'    3='U'
* Local externals:
*    KF(J)   I     I O    -         Read/Write flag
*
* External calls:
*    -
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*   (1) KF flag is reset as follows:
*                    .   R   W   U
*          IFLAG =   0   1   2   3
*                 ..................
*     Input KF=-1:   0   1   2   3  : (set each supplied bit)
*     Input KF=0 :   0   1   2   3  : (set each supplied bit)
*     Input KF=1 :   1   1   3   3  : (set read bit)
*     Input KF=2 :   2   2   2   2  : (if already  writing, leave alone)
*     Input KF=3 :   3   3   3   3  : (if already updating, leave alone)
*                 ..................
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER J,IFLAG
      INCLUDE 'memsys.inc'

      INTEGER KF(40)
      EQUIVALENCE (KF,    ICOM(KKF   ))

      IF (KF(J).LE.0) THEN
        KF(J)=IFLAG
      ELSEIF (KF(J).EQ.1) THEN
        IF (IFLAG.GE.2) KF(J)=3
      ENDIF
      END

      SUBROUTINE MESAVE
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*    Provide arrays of MemSys internal variables for the user to save.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Local externals:
*    NINTS   I     I      -         Number of INTEGERs to be preserved
*    NREALS  I     I      -         Number of REALs to be preserved
*
* External calls:
*    UINFO       User's diagnostics handler
*    USAVE       User's save routine
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INCLUDE 'memsys.inc'

      CALL UINFO('   MeSave',1)
      CALL USAVE(ICOM,NINTS,RCOM,NREALS)
      END

      SUBROUTINE MEREST
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*    Restore MemSys internal variables from the user.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Local externals:
*    NINTS   I     I      -         Number of INTEGERs to be preserved
*    NREALS  I     I      -         Number of REALs to be preserved
*
* External calls:
*    MCLEAR      Clear read/write flags
*    UINFO       User's diagnostics handler
*    UREST       User's restore routine
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INCLUDE 'memsys.inc'

      CALL UINFO('   MeRest',1)
      CALL UREST(ICOM,NINTS,RCOM,NREALS)
      CALL MCLEAR
      END

      SUBROUTINE MEINIT
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*   Initialise vector space dimensions and block size to impossible
*     values (to check initialisation), allocate areas to spaces,
*     and initialise the read/write flags.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Local externals:
*    TOL     R       O    -         Arithmetic tolerance  (1.+TOL>1.)
*    KF      I       O    40        Read/Write flags (-1)
*    IOFF    I       O    -         Relative address within area (0.)
*    MORE    I       O    -         Any more elements available? (0.)
*
* External calls:
*    VRAND0      Random generator initialisation
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*   (1) In the unusual event that automatic tolerance determination
*       fails, set TOL explicitly, e.g.
*          TOL = 1.0E-7  for single precision,
*          TOL = 1.0D-16 for double precision.
*       The exact value is not critical.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INCLUDE 'memsys.inc'

      SAVE /MEMCOM/
      INTEGER KF(40),IOFF,MORE
      REAL TOL
      EQUIVALENCE (KF,    ICOM(KKF   ))
      EQUIVALENCE (IOFF,  ICOM(KIOFF )),(MORE,  ICOM(KMORE ))
      EQUIVALENCE (TOL,   RCOM(KTOL  ))

      REAL ZERO,ONE,THREE,FOUR,X13,X43
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,THREE=3.0D0,FOUR=4.0D0,
     *           X13=ONE/THREE,X43=FOUR/THREE)

      INTEGER I
* Arithmetic tolerance = smallest number for which    (1+TOL)  >  1
      TOL=ABS((X43-ONE)-X13)
      IF (TOL.LE.ZERO)
     *  STOP ' TOLerance determination failed - see MEINIT source code'
* Initialise random number generator
      CALL VRAND0(0)

* Read/write flags and vector space allocation
      DO 1 I=1,40
        KF(I)=-1
    1 CONTINUE
      IOFF=0
      MORE=0
      END

** library.for

      SUBROUTINE MEZERO(ST,M)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*
*            <M> := zero
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    M       I     I      -         Output area number
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <M>            O              Output area
*
* External calls:
*    MFILL       Fill one block
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      INTEGER M

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL MWHILE
      REAL ZERO
      PARAMETER (ZERO=0.0D0)
      EXTERNAL MWHILE

      WRITE (INFO,'(''   MeZero('',I2,'')'')') M
      CALL UINFO(INFO,1)
    1   CALL MFILL(ST,M,ZERO)
      IF (MWHILE(ST)) GOTO 1
      END

      SUBROUTINE MEHSET(ST)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*
*            <5> := h
*
*      when alpha = infinity
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
* Local externals:
*    METHD1  I     I      -         Switch for entropy type
*    DEF     R     I      -         Default level (used only if > 0)
* Areas:
*     <3>          I                Default model m (if DEF <= 0)
*     <5>            O              Hidden distribution h
*
* External calls:
*    MDIV        Vector divide
*    MFILL       Fill one block
*    MMOV        Copy one block
*    MSADD       Scalar add
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      INCLUDE 'memsys.inc'

      INTEGER METHD1
      REAL DEF
      EQUIVALENCE (METHD1,ICOM(KMETH1))
      EQUIVALENCE (DEF,   RCOM(KDEF))

      LOGICAL MWHILE
      REAL ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      EXTERNAL MWHILE

      CALL UINFO('   MeHSet',1)
*    1   IF (METHD1.EQ.1) THEN
*          IF (DEF.GT.ZERO) THEN
*            CALL MFILL(ST,5,DEF)
*          ELSE
 1    CALL MMOV(ST,3,5)
*          ENDIF
*        ELSEIF (METHD1.EQ.2) THEN
*          CALL MFILL(ST,5,ZERO)
*        ELSEIF (METHD1.EQ.3) THEN
*          IF (DEF.GT.ZERO) THEN
*            CALL MFILL(ST,5,DEF/(ONE+DEF))
*          ELSE
*            CALL MSADD(ST,3,ONE,5)
*            CALL MDIV(ST,3,5,5)
*          ENDIF 
*        ELSEIF (METHD1.EQ.4) THEN
*          CALL MFILL(ST,5,ZERO)
*        ENDIF
      IF (MWHILE(ST)) GOTO 1
      END

      SUBROUTINE MESMUL(ST,M,A,N)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*             <N> := A * <M>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    M       I     I      -         Input area number
*    A       R     I      -         Scalar
*    N       I     I      -         Output area number
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <M>          I                Input area
*     <N>            O              Output area
*
* External calls:
*    MSMUL       Scalar multiply
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER M,N
      REAL ST(0:*)
      REAL A

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL MWHILE
      EXTERNAL MWHILE

      WRITE (INFO,'(''   MeSmul('',I2,'','',I2,'')'')') M,N
      CALL UINFO(INFO,1)
    1   CALL MSMUL(ST,M,A,N)
      IF (MWHILE(ST)) GOTO 1
      END

      SUBROUTINE MEVMUL(ST,L,M,N)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*             <N> := <L> * <M>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    L       I     I      -         Input area number
*    M       I     I      -         Input area number
*    N       I     I      -         Output area number
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <L>          I                Input area
*     <M>          I                Input area
*     <N>            O              Output area
*
* External calls:
*    MMUL        Multiply
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER L,M,N
      REAL ST(0:*)

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL MWHILE
      EXTERNAL MWHILE

      WRITE (INFO,'(''   MeVmul('',I2,'','',I2,'','',I2,'')'')') L,M,N
      CALL UINFO(INFO,1)
    1   CALL MMUL(ST,L,M,N)
      IF (MWHILE(ST)) GOTO 1
      END

      SUBROUTINE MEMSMA(ST,L,X,M,N)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                      <N> :=  <L> * X + <M>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    L       I     I      -         Input area number
*    X       R     I      -         Scalar
*    M       I     I      -         Input area number
*    N       I     I      -         Output area number
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <L>          I                Input area
*     <M>          I                Input area
*     <N>            O              Output area
*   
*
* External calls:
*    MSMULA      Scalar multiply and add
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      INTEGER L,M,N
      REAL X

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL MWHILE
      EXTERNAL MWHILE

      WRITE (INFO,'(''   MemSma('',I2,'','',I2,'','',I2,'')'')') L,M,N
      CALL UINFO(INFO,1)
    1   CALL MSMULA(ST,L,X,M,N)
      IF (MWHILE(ST)) GOTO 1
      END

      SUBROUTINE MEMUPD(ST,L,M,N)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                      <N> :=  <L> + <M>
*
*              keeping elements of <N> within allowed limits.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    L       I     I      -         Input area number
*    M       I     I      -         Input area number
*    N       I     I      -         Output area number
* Local externals:
*    METHD1  I     I      -         Switch for entropy type
*    TOL     R     I      -         Arithmetical tolerance
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <L>          I                Input area
*     <M>          I                Input area
*     <N>            O              Output area
*
* External calls:
*    MADD        Vector add
*    MUPDT1      Vector add (ensuring 0.0 < result)
*    MUPDT3      Vector add (ensuring 0.0 < result < 1.0)
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      INTEGER L,M,N
      INCLUDE 'memsys.inc'

      INTEGER METHD1
      EQUIVALENCE (METHD1,ICOM(KMETH1))
      REAL TOL
      EQUIVALENCE (TOL,   RCOM(KTOL)  )

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL MWHILE
      EXTERNAL MWHILE

      WRITE (INFO,'(''   MemUpd('',I2,'','',I2,'','',I2,'')'')') L,M,N
      CALL UINFO(INFO,1)
    1   IF (METHD1.EQ.1) THEN
          CALL MUPDT1(ST,L,M,N) 
        ELSEIF (METHD1.EQ.3) THEN
          CALL MUPDT1(ST,L,M,N)
        ELSEIF (METHD1.EQ.4) THEN
          CALL MUPDT1(ST,L,M,N) 
*          CALL MUPDT3(ST,L,M,N,TOL) 
        ELSE
          CALL MADD(ST,L,M,N)
        ENDIF
      IF (MWHILE(ST)) GOTO 1
      END

      SUBROUTINE MEMSMD(ST,L,X,M,PROD)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                      <M> :=  <L> * X + <M>
*                     PROD :=  <M>.<M>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    L       I     I      -         Input area number
*    X       R     I      -         Scalar
*    M       I     I      -         Input/Output area number
*    PROD    R     O      -         Dot product
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <L>          I                Input area
*     <M>          I O              Input/Output area
*   
*
* External calls:
*    MDOT        Scalar product
*    MSMULA      Scalar multiply and add
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      INTEGER L,M
      REAL X,PROD

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL MWHILE
      REAL A,PPROD,ZERO
      PARAMETER (ZERO=0.0D0)
      EXTERNAL MWHILE

      WRITE (INFO,'(''   MemSmd('',I2,'','',I2,'')'')') L,M
      CALL UINFO(INFO,1)
      A=X
      PROD=ZERO
    1   CALL MSMULA(ST,L,A,M,M)
        CALL MDOT(ST,M,M,PPROD)
        PROD=PROD+PPROD
      IF (MWHILE(ST)) GOTO 1
      END

      SUBROUTINE MEMDOT(ST,M,N,PROD)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                      PROD :=  <M>.<N>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    M       I     I      -         Input area number
*    N       I     I      -         Input area number
*    PROD    R     I      -         Scalar product
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <M>          I                Input area
*     <N>          I                Input area
*
* External calls:
*    MDOT        Dot product
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER M,N
      REAL ST(0:*)
      REAL PROD

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL MWHILE
      REAL PPROD,ZERO
      PARAMETER (ZERO=0.0D0)
      EXTERNAL MWHILE

      WRITE (INFO,'(''   MemDot('',I2,'','',I2,'')'')') M,N
      CALL UINFO(INFO,1)
      PROD=ZERO
    1   CALL MDOT(ST,M,N,PPROD)
        PROD=PROD+PPROD
      IF (MWHILE(ST)) GOTO 1
      END

      SUBROUTINE MECOPY(ST,M,N)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*          Copy area M to N   ,   <N> := <M>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    M       I     I      -         Input area number
*    N       I     I      -         Output area number
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <M>          I                Input area
*     <N>            O              Output area
*
* External calls:
*    MMOV        Copy one block
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
* (1)  Replaces MemSys3 routines MEMJ and MEMK
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER M,N
      REAL ST(0:*)

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL MWHILE
      EXTERNAL MWHILE

      IF (M.NE.N) THEN
        WRITE (INFO,'(''   MeCopy('',I2,'','',I2,'')'')') M,N
        CALL UINFO(INFO,1)
    1     CALL MMOV(ST,M,N)
        IF (MWHILE(ST)) GOTO 1
      ENDIF
      END

      SUBROUTINE MESWAP(ST,M,N)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*          Interchange areas <M> and <N>
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    M       I     I      -         First area number
*    N       I     I      -         Second area number
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <M>          I O              First area
*     <N>          I O              Second area
*
* External calls:
*    MSWAP       Copy one block
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER M,N
      REAL ST(0:*)

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL MWHILE
      EXTERNAL MWHILE

      IF (M.NE.N) THEN
        WRITE (INFO,'(''   MeSwap('',I2,'','',I2,'')'')') M,N
        CALL UINFO(INFO,1)
    1     CALL MSWAP(ST,M,N)
        IF (MWHILE(ST)) GOTO 1
      ENDIF
      END

      SUBROUTINE MEMENT(ST,S,SUMMET)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*          Generate metric, entropy and gradS from h
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    S       R       O    -         Entropy
*    SUMMET  R       O    -         SUM([metric])
* Local externals:
*    METHD1  I     I      -         Switch for entropy type
*    DEF     R     I      -         Default level (used only if > 0)
* Areas:
*     <1>            O              SQRT([metric])
*     <2>            O              -gradS
*     <3>          I                Default model m
*     <5>          I                Hidden distribution h
*
* External calls:
*    MENT1       Vectors and scalars for standard entropy
*    MENT2       Vectors and scalars for positive/negative entropy
*    MENT3       Vectors and scalars for Fermi-Dirac entropy
*    MENT4       Vectors and scalars for quadratic 'entropy'
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*   (1)  This routine contains the entire definitions of entropy and
*        related quantities.
*        The metric is  [metric] = -1/(d2S/dh2).
*   (2)  The entropy S must be uncorrelated and strictly convex, so that
*        [metric] is diagonal and positive definite.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      REAL S,SUMMET
      INCLUDE 'memsys.inc'

      INTEGER METHD1
      EQUIVALENCE (METHD1,ICOM(KMETH1))
      REAL DEF
      EQUIVALENCE (DEF,   RCOM(KDEF  ))

      LOGICAL MWHILE
      REAL PS,PSUM
      REAL ZERO
      PARAMETER (ZERO=0.0D0)
      EXTERNAL MWHILE

      CALL UINFO('   MemEnt',1)
      S=ZERO
      SUMMET=ZERO
    1   IF (METHD1.EQ.1) THEN
          CALL MENT1(ST,DEF, PS,PSUM)
        ELSEIF (METHD1.EQ.2) THEN
          CALL MENT2(ST,DEF, PS,PSUM)
        ELSEIF (METHD1.EQ.3) THEN
          CALL MENT3(ST,DEF, PS,PSUM)
        ELSEIF (METHD1.EQ.4) THEN
          CALL MENT4(ST,DEF, PS,PSUM)
        ENDIF
        S=S+PS
        SUMMET=SUMMET+PSUM
      IF (MWHILE(ST)) GOTO 1
      END

      SUBROUTINE MENT1(ST,DEF, S,SUM)
* One block of standard entropy
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*),DEF,S,SUM
      REAL ZERO,ONE,A,B
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      IF (DEF.GT.ZERO) THEN
        CALL MFILL(ST,2,DEF)
        CALL MSUM(ST,2,A)
        CALL MSMUL(ST,5,ONE/DEF,2)
      ELSE
        CALL MSUM(ST,3,A)
        CALL MDIV(ST,5,3,2)
      ENDIF
      CALL MSUM(ST,5,SUM)
      CALL MSQRT(ST,5,1)
      CALL MLOG(ST,2,2)
      CALL MDOT(ST,5,2,B)
      S=SUM-A-B
      END

      SUBROUTINE MENT2(ST,DEF, S,SUM)
* One block of positive/negative entropy
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*),DEF,S,SUM
      INTEGER KC1,KC2,KC3,KC4,N
      REAL ZERO
      PARAMETER (ZERO=0.0D0)
      CALL MFETCH(ST,5,1, KC1,N)
      IF (DEF.LE.ZERO) CALL MFETCH(ST,3,1, KC2,N)
      CALL MFETCH(ST,2,2, KC3,N)
      CALL MFETCH(ST,1,2, KC4,N)
      CALL VENT2(ST,KC1,KC2,KC3,KC4,DEF,S,SUM,N)
      END


      SUBROUTINE MENT3(ST,DEF, S,SUM)
* L1L2w entropy S(h)= m log( 1 + h / m ) - h 
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*),DEF,S,SUM
      REAL ZERO,ONE,A,B
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)  
      CALL MSUM(ST,5,A)
      CALL MADD(ST, 3, 5, 1)
      CALL MDIV(ST, 1, 3, 2)
      CALL MLOG( ST, 2, 2)
      CALL MDOT( ST, 3, 2, B)
      S = B - A
      CALL MMUL( ST, 1,  1 , 1) 
      CALL MDIV( ST, 1, 3 , 1)
      CALL MSUM(ST, 1, SUM)
      CALL MSQRT(ST, 1 , 1)
      CALL MADD(ST, 3, 5, 2)
      CALL MDIV(ST, 5, 2, 2)
      END

      SUBROUTINE MENT4(ST,DEF, S,SUM)
* One block of quadratic 'entropy'
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*),DEF,S,SUM
      REAL ZERO,ONE,TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      IF (DEF.GT.ZERO) THEN
        CALL MFILL(ST,2,DEF)
        CALL MSUM(ST,2,SUM)
        CALL MFILL(ST,1,SQRT(DEF))
        CALL MSMUL(ST,5,ONE/DEF,2)
      ELSE
        CALL MSUM(ST,3,SUM)
        CALL MSQRT(ST,3,1)
        CALL MDIV(ST,5,3,2)
      ENDIF
      CALL MDOT(ST,2,5,S)
      S=-S/TWO
      END

      SUBROUTINE MENT5(ST,DEF, S,SUM)
* One block of positive/negative entropy - old bugged routine
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*),DEF,S,SUM
      REAL ZERO,ONE,TWO,A,B,DEF2
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      IF (DEF.GT.ZERO) THEN
        DEF2=TWO*DEF
        CALL MFILL(ST,2,DEF2)
        CALL MSUM(ST,2,A)
        CALL MSMUL(ST,5,ONE/DEF2,2)
      ELSE
        CALL MSMUL(ST,3,TWO,2)
        CALL MSUM(ST,2,A)
        CALL MDIV(ST,5,2,2)
      ENDIF
      CALL MMUL(ST,2,2,1)
      CALL MSADD(ST,1,ONE,1)
      CALL MSQRT(ST,1,1)
      CALL MADD(ST,1,2,2)
      IF (DEF.GT.ZERO) THEN
        CALL MSMUL(ST,1,DEF2,1)
      ELSE
        CALL MMUL(ST,1,3,1)
        CALL MSMUL(ST,1,TWO,1)
      ENDIF
      CALL MSUM(ST,1,SUM)
      CALL MSQRT(ST,1,1)
      CALL MLOG(ST,2,2)
      CALL MDOT(ST,2,5,B)
      S=SUM-A-B
      END

      SUBROUTINE MEMCHI(ST,ALHOOD,DATA,UNITS)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*         Generate residuals and likelihood from mock Gaussian data, for
*              pr(Data|Mock) = exp( UNITS - ALHOOD ) ,
*
*             <24>  := ( <21> - <25> ) * <22>
*                    = ( Data - Mockdata ) /sigma = Scaled residuals
*            ALHOOD := <24>.<24> / 2          ( >= 0)
*
*         Also generate dimensionalities
*             DATA  := Number of nonzero accuracies
*             UNITS := Sum log ( nonzero accuracies )
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    ALHOOD  R       O    -         L = - log likelihood
*    DATA    R       O    -         Number of data
*    UNITS   R       O    -         Sum log accuracy
* Local externals:
*    ACC     R     I      -         Accuracies (used only if > 0)
*    LENGTH  I     I      -         Dimension of data space
* Areas:
*    <21>          I                Data
*    <22>         (I)               [acc]
*    <24>            O              Normalised residuals
*    <25>          I                Mock data
*
* External calls:
*    MDIV        Vector divide
*    MDOT        Dot product
*    MLOG        Vector logarithm
*    MMUL        Vector multiply
*    MRECIP      Vector reciprocal
*    MSADD       Scalar add
*    MSMUL       Scalar multiply
*    MSUB        Vector subtract
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    John Skilling   7 Dec 1988     First MEMSYS3 version
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*   (1) This routine contains the entire definition of likelihood and
*       related quantities, with zero likelihood minimum at F=D.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      REAL ALHOOD,DATA,UNITS
      INCLUDE 'memsys.inc'

      INTEGER LENGTH
      REAL ACC
      EQUIVALENCE (LENGTH,ICOM(KLNGTH))
      EQUIVALENCE (ACC,   RCOM(KACC  ))

      LOGICAL MWHILE
      REAL ZERO,EPS,TWO,PI2,PLHOOD,D,PSUML
      PARAMETER (ZERO=0.0D0,EPS=1.0D-20,TWO=2.0D0)
*  PI2 = 1/(2*pi)
      PARAMETER (PI2  =0.1591549430918953357688838D0)
      EXTERNAL MWHILE

      CALL UINFO('   MemChi',1)
      ALHOOD=ZERO
      DATA=ZERO
      UNITS=ZERO
      IF (ACC.GT.ZERO) THEN
    3     CALL MSUB(ST,21,25,24)
          CALL MSMUL(ST,24,ACC,24)
          CALL MDOT(ST,24,24,PLHOOD)
          ALHOOD=ALHOOD+PLHOOD
        IF (MWHILE(ST)) GOTO 3
        DATA=FLOAT(LENGTH)
        UNITS=DATA*LOG(ACC)
      ELSE
    4     CALL MSADD(ST,22,EPS,26)
          CALL MLOG(ST,26,24)
          CALL MDIV(ST,24,26,24)
          CALL MDOT(ST,24,22, PSUML)
          CALL MRECIP(ST,26,24)
          CALL MDOT(ST,24,22, D)
          DATA=DATA+D
          UNITS=UNITS+PSUML
        IF (MWHILE(ST)) GOTO 4
    5     CALL MSUB(ST,21,25,24)
          CALL MMUL(ST,24,22,24)
          CALL MDOT(ST,24,24,PLHOOD)
          ALHOOD=ALHOOD+PLHOOD
        IF (MWHILE(ST)) GOTO 5
      ENDIF
      ALHOOD=ALHOOD/TWO
      END

      SUBROUTINE MEMPOI(ST,ALHOOD,DATA,UNITS)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*         Generate residuals and likelihood from mock Poisson data, for
*              pr(Data|Mock) = exp( UNITS - ALHOOD ) ,
*
*             <22>  := sqrt(<21>+1) / <25>
*                    = sqrt(gradgrad(-log(lhood))
*                    = sqrt(data) / mockdata = accuracies
*             <24>  := ( <21> - <25> ) / sqrt(<21>+1)
*                    = - grad (-log Lhood) / sqrt(gradgrad(-log Lhood))
*                    = Scaled residuals
*            ALHOOD := SUM ( <25> - <21> + <21> log(<21>/<25>) )
*
*         Also generate dimensionalities
*             DATA  := Number of nonzero accuracies
*             UNITS := Sum log ( nonzero accuracies )
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    ALHOOD  R       O    -         L = - log likelihood
*    DATA    R       O    -         Number of data
*    UNITS   R       O    -         Sum log accuracy
* Local externals:
*    LENGTH  I     I      -         Dimension of data space
* Areas:
*    <21>          I                Data
*    <22>          I O              [acc]**-2
*    <24>            O              Normalised residuals
*    <25>          I                Mock data
*
* External calls:
*    MADD        Vector add
*    MDOT        Dot product
*    MLOG        Vector logarithm
*    MMUL        Vector multiply
*    MSADD       Scalar add
*    MSQRT       Vector square root
*    MSUB        Vector subtract
*    MSUM        Vector sum
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    John Skilling   7 Dec 1988     First MEMSYS3 version
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*   (1) This routine contains the entire definition of likelihood and
*       related quantities, with zero likelihood minimum at F=D.
*   (2) Proper operation of this routine requires positive mock data.
*   (3) Adding 1 to D when setting Poisson accuracies can be justified
*       by averaging gradgradL over a reasonable range of values of F.
*   (4) The Poisson contribution ( D log D - D - log D! ) to 
*       UNITS is here approximated by  -(1/2) log( 1 + 2*pi*D ) .
*       Importantly, error = 0 for D=0,  error -> 0 as D -> infinity.
*       Max error = 0.007 (less than 1%) at D=1.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      REAL ALHOOD,DATA,UNITS
      INCLUDE 'memsys.inc'

      INTEGER LENGTH
      EQUIVALENCE (LENGTH,ICOM(KLNGTH))

      LOGICAL MWHILE
      REAL ZERO,EPS,ONE,TWO,PI2,PLHOOD,PSUML
      PARAMETER (ZERO=0.0D0,EPS=1.0D-20,ONE=1.0D0,TWO=2.0D0)
*  PI2 = 1/(2*pi)
      PARAMETER (PI2  =0.1591549430918953357688838D0)
      EXTERNAL MWHILE

      CALL UINFO('   MemPoi',1)
      ALHOOD=ZERO
      UNITS=ZERO
    3   CALL MDIV(ST,21,25,22)
        CALL MSADD(ST,22,EPS,22)
        CALL MLOG(ST,22,22)
        CALL MMUL(ST,21,22,22)
        CALL MSUB(ST,22,21,22)
        CALL MADD(ST,25,22,22)
        CALL MSUM(ST,22,PLHOOD)
        CALL MSADD(ST,21,PI2,22)
        CALL MLOG(ST,22,24)
        CALL MSUM(ST,24,PSUML)
        CALL MSADD(ST,21,ONE,22)
        CALL MSQRT(ST,22,22)
        CALL MSUB(ST,21,25,24)
        CALL MDIV(ST,24,22,24)
        CALL MDIV(ST,22,25,22)
        ALHOOD=ALHOOD+PLHOOD
        UNITS=UNITS+PSUML
      IF (MWHILE(ST)) GOTO 3
      DATA=FLOAT(LENGTH)
      UNITS=-UNITS/TWO        
      END

      SUBROUTINE METEST(ST,TEST)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*               TEST = 1 - cos(angle(<2>,<4>))
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    TEST    R       O    -         1 - cos(angle( gradS, gradL ))
* Areas:
*     <2>          I                -gradS
*     <4>          I                -gradL
*
* External calls:
*    MDOT        Dot product
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      REAL TEST

      LOGICAL MWHILE
      REAL ZERO,EPS,ONE,SS,SC,CC,PSS,PSC,PCC,DIV,XXX,YYY
      PARAMETER (ZERO=0.0D0,EPS=1.0D-20,ONE=1.0D0)
      EXTERNAL MWHILE

      DIV(XXX,YYY) = XXX/MAX(YYY,EPS)

      CALL UINFO('   MeTest',1)
      SS=ZERO
      SC=ZERO
      CC=ZERO
    1   CALL MDOT(ST,2,2,PSS)
        CALL MDOT(ST,2,4,PSC)
        CALL MDOT(ST,4,4,PCC)
        SS=SS+PSS
        SC=SC+PSC
        CC=CC+PCC
      IF (MWHILE(ST)) GOTO 1
      TEST=ONE-DIV(SC,SQRT(SS)*SQRT(CC))
      END

      SUBROUTINE MEVRND(ST,M,N,NCORR)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                       <M> := random vector
*
* IF    ( NCORR = 0 ) THEN
*                   components are uncorrelated random signs
* ELSEIF( NCORR > 0 ) THEN
*                   components are from N(0,1),
*                   with possible correlation with previous samples.
*
* During several successive calls to this routine with the same value of
* NCORR (strictly positive), let the m'th and n'th vector results be
*     (m)                             (n)
*    x    (i=0,1,2,....)     and     x    (i=0,1,2,....)
*     i                               i
* All values are individually from N(0,1), but have covariance structure
*                        -
*       (m)  (n)        |  max( 1 - abs(m-n)/NCORR , 0 )    if  i=j
*    < x    x   >   =   |
*       i    j          |    0   if i not equal to j
*                        -
* Examples: NCORR = 3  gives (for each i=j) covariance
*
*                       | 3  2  1  0  0 |  
*       (m)  (n)        | 2  3  2  1  0 |
*    < x    x   >   =   | 1  2  3  2  1 | / 3
*                       | 0  1  2  3  2 |
*                       | 0  0  1  2  3 |
*
*           between the first five samples.
*           NCORR = 1 would have given the uncorrelated identity.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    M       I     I      -         Output area number
*    N       I     I      -         Work area (not M if NCORR>1)
*    NCORR   I     I      -         =0 for random signs, or
*                                   >0 for normal, possibly correlated
* Globals:
*    INFO    C*78   (O)   -         Diagnostics workspace
* Areas:
*     <M>            O              Output area
*     <N>            O              Work area (if NCORR>1)
*
* External calls:
*    MRAND       Normal deviates
*    MSIGN       Vector signs
*    MSMUL       Area scalar multiply
*    MSMULA      Area scalar multiply and add
*    MWHILE      Manage blocking     
*    UINFO       User's diagnostics handler
*    VRAND1      Save random generator
*    VRAND2      Restore random generator
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*  (1) For correlated Gaussian, CPU time proportional to NCORR
*  (2) In this implementation, the vector samples depend on addressing
*      details such as the block-size.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      REAL ST(0:*)
      INTEGER M,N,NCORR

      CHARACTER INFO*78
      COMMON /MEINFO/ INFO

      LOGICAL MWHILE
      INTEGER I
      REAL ZERO,ONE,TEMP
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      EXTERNAL MWHILE

      WRITE (INFO,'(''   MeVrnd('',I2,'')'')') M
      CALL UINFO(INFO,1)
* Random signs
      IF( NCORR.LE.0 ) THEN
    1     CALL MRAND(ST,M)
          CALL MSIGN(ST,M,M)
        IF (MWHILE(ST)) GOTO 1
* Normal deviates
      ELSEIF( NCORR.EQ.1 ) THEN
* Independent samples
    2     CALL MRAND(ST,M)
        IF (MWHILE(ST)) GOTO 2
      ELSE
* Correlated samples
        IF( N.EQ.M ) STOP ' MEVRND needs separate output and work areas'
* Normalise the output correctly.
        TEMP=SQRT(ONE/FLOAT(NCORR))
    3     CALL MRAND(ST,N)
          CALL MSMUL(ST,N,TEMP,M)
        IF (MWHILE(ST)) GOTO 3
* After the first pass, store the random generator
        CALL VRAND1
        DO 5 I=2,NCORR
    4       CALL MRAND(ST,N)
            CALL MSMULA(ST,N,TEMP,M,M)
          IF (MWHILE(ST)) GOTO 4
    5   CONTINUE
* After all NCORR passes are done, recover the random generator so
* that the last NCORR-1 samples will be repeated at the start of the
* next call.  This is the source of the correlation.
        CALL VRAND2
      ENDIF
      END

** msubs.for

      SUBROUTINE MFILL(ST,J1, A)
* One block of  <J1> := A
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1
      REAL ST(0:*)
      REAL A
      INTEGER KC1,N
      CALL MFETCH(ST,J1,2, KC1,N)
      CALL VFILL(ST,KC1, A,N)
      END

      SUBROUTINE MMOV(ST,J1, J2)
* One block of  <J2> := <J1>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2
      REAL ST(0:*)
      INTEGER KC1,KC2,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,2, KC2,N)
      CALL VMOV(ST,KC1, KC2,N)
      END

      SUBROUTINE MSWAP(ST,J1,J2)
* Interchange one block of <J1> and <J2>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2
      REAL ST(0:*)
      INTEGER KC1,KC2,N
      CALL MFETCH(ST,J1,3, KC1,N)
      CALL MFETCH(ST,J2,3, KC2,N)
      CALL VSWAP(ST,KC1,KC2, N)
      END

      SUBROUTINE MSIGN(ST,J1, J2)
* One block of  <J2> := sign<J1>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2
      REAL ST(0:*)
      INTEGER KC1,KC2,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,2, KC2,N)
      CALL VSIGN(ST,KC1, KC2,N)
      END

      SUBROUTINE MADD(ST,J1,J2, J3)
* One block of  <J3> := <J1> + <J2>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2,J3
      REAL ST(0:*)
      INTEGER KC1,KC2,KC3,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,1, KC2,N)
      CALL MFETCH(ST,J3,2, KC3,N)
      CALL VADD(ST,KC1,KC2, KC3,N)
      END

      SUBROUTINE MSUB(ST,J1,J2, J3)
* One block of  <J3> := <J1> - <J2>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2,J3
      REAL ST(0:*)
      INTEGER KC1,KC2,KC3,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,1, KC2,N)
      CALL MFETCH(ST,J3,2, KC3,N)
      CALL VSUB(ST,KC1,KC2, KC3,N)
      END

      SUBROUTINE MMUL(ST,J1,J2, J3)
* One block of  <J3> := <J1> * <J2>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2,J3
      REAL ST(0:*)
      INTEGER KC1,KC2,KC3,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,1, KC2,N)
      CALL MFETCH(ST,J3,2, KC3,N)
      CALL VMUL(ST,KC1,KC2, KC3,N)
      END

      SUBROUTINE MDIV(ST,J1,J2, J3)
* One block of  <J3> := <J1> / <J2>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2,J3
      REAL ST(0:*)
      INTEGER KC1,KC2,KC3,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,1, KC2,N)
      CALL MFETCH(ST,J3,2, KC3,N)
      CALL VDIV(ST,KC1,KC2, KC3,N)
      END

      SUBROUTINE MRECIP(ST,J1, J2)
* One block of  <J2> := 1. / <J1>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2
      REAL ST(0:*)
      INTEGER KC1,KC2,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,2, KC2,N)
      CALL VRECIP(ST,KC1, KC2,N)
      END

      SUBROUTINE MLOG(ST,J1,J2)
* One block of  <J2> := log <J1>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2
      REAL ST(0:*)
      INTEGER KC1,KC2,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,2, KC2,N)
      CALL VLOG(ST,KC1,KC2, N)
      END

      SUBROUTINE MEXP(ST,J1,J2)
* One block of  <J2> := exp <J1>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2
      REAL ST(0:*)
      INTEGER KC1,KC2,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,2, KC2,N)
      CALL VEXP(ST,KC1,KC2, N)
      END

      SUBROUTINE MSQRT(ST,J1, J2)
* One block of  <J2> := sqrt <J1>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2
      REAL ST(0:*)
      INTEGER KC1,KC2,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,2, KC2,N)
      CALL VSQRT(ST,KC1, KC2,N)
      END

      SUBROUTINE MSADD(ST,J1,A, J2)
* One block of  <J2> := <J1> + A
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2
      REAL ST(0:*)
      REAL A
      INTEGER KC1,KC2,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,2, KC2,N)
      CALL VSADD(ST,KC1,A, KC2,N)
      END

      SUBROUTINE MSMUL(ST,J1,A, J2)
* One block of  <J2> := A * <J1>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2
      REAL ST(0:*)
      REAL A
      INTEGER KC1,KC2,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,2, KC2,N)
      CALL VSMUL(ST,KC1,A, KC2,N)
      END

      SUBROUTINE MSMULA(ST,J1,A,J2, J3)
* One block of  <J3> := A*<J1> + <J2>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2,J3
      REAL ST(0:*)
      REAL A
      INTEGER KC1,KC2,KC3,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,1, KC2,N)
      CALL MFETCH(ST,J3,2, KC3,N)
      CALL VSMULA(ST,KC1,A,KC2, KC3,N)
      END

      SUBROUTINE MUPDT1(ST,J1,J2, J3)
* One block of  <J3> := MAX(<J1> + <J2>, <J3>/5, FLOOR)
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2,J3
      REAL ST(0:*)
      INTEGER KC1,KC2,KC3,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,1, KC2,N)
      CALL MFETCH(ST,J3,3, KC3,N)
      CALL VUPDT1(ST,KC1,KC2, KC3,N)
      END

      SUBROUTINE MUPDT3(ST,J1,J2, J3, TOL)
* One block of  <J3> := TOL <= <J1> + <J2> <= 1.0 - TOL
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2,J3
      REAL ST(0:*),TOL
      INTEGER KC1,KC2,KC3,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,1, KC2,N)
      CALL MFETCH(ST,J3,3, KC3,N)
      CALL VUPDT3(ST,KC1,KC2, KC3, TOL, N)
      END

      SUBROUTINE MSUM(ST,J1, SUM)
* One block of  SUM := SUM( <J1> )
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1
      REAL ST(0:*)
      REAL SUM
      INTEGER KC1,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL VSUM(ST,KC1, SUM,N)
      END

      SUBROUTINE MDOT(ST,J1,J2, SUM)
* One block of  SUM := <J1>.<J2>
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1,J2
      REAL ST(0:*)
      REAL SUM
      INTEGER KC1,KC2,N
      CALL MFETCH(ST,J1,1, KC1,N)
      CALL MFETCH(ST,J2,1, KC2,N)
      CALL VDOT(ST,KC1,KC2, SUM,N)
      END

      SUBROUTINE MRAND(ST,J1)
* One block of  <J1> := Random normal
      IMPLICIT CHARACTER (A-Z)
      INTEGER J1
      REAL ST(0:*)
      INTEGER KC1,N
      CALL MFETCH(ST,J1,2, KC1,N)
      CALL VRAND(ST,KC1,N)
      END

