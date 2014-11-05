*                          Vector routines
*                          ---------------
*
* Vector Library
*  VFILL (ST,JU,A,N)       U[] = A
*  VMOV  (ST,JU,JV,N)      V[] = U[]
*  VSWAP (ST,JU,JV,N)      Interchange U[] and V[]
*  VSIGN (ST,JU,JV,N)      V[] = sign(U[])
*  VADD  (ST,JU,JV,JW,N)   W[] = U[] + V[]
*  VSUB  (ST,JU,JV,JW,N)   W[] = U[] - V[]
*  VMUL  (ST,JU,JV,JW,N)   W[] = U[] * V[]
*  VDIV  (ST,JU,JV,JW,N)   W[] = U[] / V[]
*  VRECIP(ST,JU,JV,N)      V[] = 1 / U[]
*  VLOG  (ST,JU,JV,N)      V[] = log(U[])
*  VEXP  (ST,JU,JV,N)      V[] = exp(U[])
*  VSQRT (ST,JU,JV,N)      V[] = sqrt(U[])
*  VSADD (ST,JU,A,JV,N)    V[] = U[] + A
*  VSMUL (ST,JU,A,JV,N)    V[] = U[] * A
*  VSMULA(ST,JU,A,JV,JW,N) W[] = U[] * A + V[]
*  VSUM  (ST,JU,A,N)         A = Sum(U[])
*  VDOT  (ST,JU,JV,A,N)      A = U[].V[]
*  VENT2(ST,JU,JV,JW,JX,DEF,S,SUM,N) Pos/neg entropy
*
*  VRAND0(ISEED)           Initialise random generator
*  VRAND1                  Save       random generator
*  VRAND2                  Restore    random generator
*  VRAND (ST,JU,N)         U[] = random normal samples
*
* Storage addressing
*  VFETCH(ST,J,IOFF,ACTION, KCORE,LENGTH,MORE)  Pull a block in
*  VSTORE(ST,J,IOFF,ACTION, KCORE,LENGTH,MORE)  Push a block out
*
      SUBROUTINE VFILL(ST,JU,A,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,N
      REAL A
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JU+I)=A
    1 CONTINUE
      END

      SUBROUTINE VMOV(ST,JU,JV,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,N
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JV+I)=ST(JU+I)
    1 CONTINUE
      END

      SUBROUTINE VSWAP(ST,JU,JV,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,N
      REAL ST(0:*)
      INTEGER I
      REAL SWAP
      DO 1 I=0,N-1
        SWAP=ST(JU+I)
        ST(JU+I)=ST(JV+I)
        ST(JV+I)=SWAP
    1 CONTINUE
      END

      SUBROUTINE VSIGN(ST,JU,JV,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,N
      REAL ST(0:*)
      INTEGER I
      REAL ONE
      PARAMETER (ONE=1.0D0)
      DO 1 I=0,N-1
        ST(JV+I)=SIGN(ONE,ST(JU+I))
    1 CONTINUE
      END

      SUBROUTINE VADD(ST,JU,JV,JW,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,JW,N
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JW+I)=ST(JU+I)+ST(JV+I)
    1 CONTINUE
      END

      SUBROUTINE VSUB(ST,JU,JV,JW,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,JW,N
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JW+I)=ST(JU+I)-ST(JV+I)
    1 CONTINUE
      END

      SUBROUTINE VMUL(ST,JU,JV,JW,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,JW,N
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JW+I)=ST(JU+I)*ST(JV+I)
    1 CONTINUE
      END

      SUBROUTINE VDIV(ST,JU,JV,JW,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,JW,N
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JW+I)=ST(JU+I)/ST(JV+I)
    1 CONTINUE
      END


      SUBROUTINE VRECIP(ST,JU,JV,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,N
      REAL ST(0:*)
      INTEGER I
      REAL ONE
      PARAMETER (ONE=1.0D0)
      DO 1 I=0,N-1
        ST(JV+I)=ONE/ST(JU+I)
    1 CONTINUE
      END

      SUBROUTINE VLOG(ST,JU,JV,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,N
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JV+I)=LOG(ST(JU+I))
    1 CONTINUE
      END

      SUBROUTINE VSQRT(ST,JU,JV,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,N
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JV+I)=SQRT(ST(JU+I))
    1 CONTINUE
      END

      SUBROUTINE VEXP(ST,JU,JV,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,N
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JV+I)=EXP(ST(JU+I))
    1 CONTINUE
      END

      SUBROUTINE VSADD(ST,JU,A,JV,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,N
      REAL A
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JV+I)=ST(JU+I)+A
    1 CONTINUE
      END

      SUBROUTINE VSMUL(ST,JU,A,JV,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,N
      REAL A
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JV+I)=A*ST(JU+I)
    1 CONTINUE
      END

      SUBROUTINE VSMULA(ST,JU,A,JV,JW,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,JW,N
      REAL A
      REAL ST(0:*)
      INTEGER I
      DO 1 I=0,N-1
        ST(JW+I)=A*ST(JU+I)+ST(JV+I)
    1 CONTINUE
      END

      SUBROUTINE VUPDT1(ST,JU,JV,JW,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,JW,N
      REAL ST(0:*)
      REAL FLOOR,FIVE
      PARAMETER (FLOOR=1.0D-12,FIVE=5.0D0)
      INTEGER I
      DO 1 I=0,N-1
        ST(JW+I)=MAX(ST(JU+I)+ST(JV+I),ST(JW+I)/FIVE,FLOOR)
    1 CONTINUE
      END

      SUBROUTINE VUPDT3(ST,JU,JV,JW,TOL,N)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,JW,N
      REAL ST(0:*),TOL
      INTEGER I
      REAL CEILNG
      REAL ZERO,ONE,FOUR,FIVE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,FOUR=4.0D0,FIVE=5.0D0)
      CEILNG=ONE-TOL
      DO 1 I=0,N-1
        IF (ST(JV+I).LT.ZERO) THEN
          ST(JW+I)=MAX(ST(JU+I)+ST(JV+I),ST(JW+I)/FIVE,TOL)
        ELSE
          ST(JW+I)=MIN(ST(JU+I)+ST(JV+I),(FOUR+ST(JW+I))/FIVE,CEILNG)
        ENDIF
    1 CONTINUE
      END

      SUBROUTINE VSUM(ST,JU, A,N)
* Vector sum
*                        A := Sum U[.]
* Sum the product as binary tree, e.g. A=((a+b)+(c+d))+((e+f)+(g+h)),
*  effectively guaranteeing full accuracy.
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,N
      REAL A
      REAL ST(0:*)
      INTEGER I
      REAL ZERO
      PARAMETER (ZERO=0.0D0)
* This routine is a slower but more accurate version of the
*  following more easily vectorisable code ................
*      A=ZERO
*      DO 1 I=0,N-1
*        A=A+ST(JU+I)
*    1 CONTINUE
* .........................................................
* Allow for addresses up to 32 bits. Z stores the branch sums
      INTEGER J,NBIT
      PARAMETER (NBIT=32)
      REAL Z(NBIT)
      DO 1 J=1,NBIT
        Z(J)=ZERO
    1 CONTINUE
      DO 3 I=0,N-1
        A=ST(JU+I)
        J=1
    2   IF(Z(J).NE.ZERO) THEN
          A=A+Z(J)
          Z(J)=ZERO
          J=J+1
          GOTO 2
        ENDIF
        Z(J)=A
    3 CONTINUE
* Accumulate all the branches, in increasing order,
      A=Z(1)
      DO 4 J=2,NBIT
        A=A+Z(J)
    4 CONTINUE
      END

      SUBROUTINE VDOT(ST,JU,JV, A,N)
* Vector dot
*                        A := U[.].V[.]
* Sum the product as binary tree, e.g. A=((a+b)+(c+d))+((e+f)+(g+h)),
*  effectively guaranteeing full accuracy.
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,N
      REAL A
      REAL ST(0:*)
      INTEGER I
      REAL ZERO
      PARAMETER (ZERO=0.0D0)
* This routine is a slower but more accurate version of the
*  following more easily vectorisable code ................
*      A=ZERO
*      DO 1 I=0,N-1
*        A=A+ST(JU+I)*ST(JV+I)
*    1 CONTINUE
* .........................................................
* Allow for addresses up to 32 bits. Z stores the branch sums
      INTEGER J,NBIT
      PARAMETER (NBIT=32)
      REAL Z(NBIT)
      DO 1 J=1,NBIT
        Z(J)=ZERO
    1 CONTINUE
      DO 3 I=0,N-1
        A=ST(JU+I)*ST(JV+I)
        J=1
    2   IF(Z(J).NE.ZERO) THEN
          A=A+Z(J)
          Z(J)=ZERO
          J=J+1
          GOTO 2
        ENDIF
        Z(J)=A
    3 CONTINUE
* Accumulate all the branches, in increasing order,
      A=Z(1)
      DO 4 J=2,NBIT
        A=A+Z(J)
    4 CONTINUE
      END

      SUBROUTINE VENT2(ST,JU,JV,JW,JX,DEF,S,SUM,N)
* U = h, V = m, W = -gradS, X = sqrt(metric)
      IMPLICIT CHARACTER (A-Z)
      INTEGER JU,JV,JW,JX,N
      REAL ST(0:*),DEF,S,SUM

      INTEGER I,ITERM,NTERM
      PARAMETER (NTERM=6)
      REAL H,M,A,X,Y,Z,RX2,LOGZ
      REAL ZERO,ONE,TWO,TEN
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,TEN=10.0D0)

      S=ZERO
      SUM=ZERO
      DO 1 I=0,N-1
        H=ST(JU+I)
        IF (DEF.GT.ZERO) THEN
          M=DEF
        ELSE
          M=ST(JV+I)
        ENDIF
        A=TWO*M
        X=H/A
        Y=SQRT(ONE+X**2)
        IF (X.GT.-TEN) THEN
          Z=Y+X
        ELSE
          RX2=ONE/X**2
          Z=ONE
          DO 10 ITERM=NTERM,2,-1
            Z=ONE+(1.5-ITERM)*RX2*Z/ITERM
   10     CONTINUE
          Z=-Z/(TWO*X)
        ENDIF
        LOGZ=LOG(Z)
        SUM=SUM+A*Y
        ST(JW+I)=LOGZ
        S=S+A*(Y-ONE-X*LOGZ)
        ST(JX+I)=SQRT(A*Y)
    1 CONTINUE
      END

      SUBROUTINE VRAND0(ISEED)
* with     ENTRY VRAND1
* and      ENTRY VRAND2
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*           VRAND0: Initialise random number generator
*           VRAND1: Save       random number generator
*           VRAND2: Restore    random number generator
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ISEED   I     I      *         Seed
* Local externals:
*    P1-P6   R       O    -         Random generator state
* Saved:
*    Q1-Q6   R     I O    -         Stored generator state
*
* External calls:
*     -
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*   (1) Call to VRAND0 is essential before any call to VRAND.
*   (2) This package approximates a Gaussian by the sum of 12
*       uniform deviates: it returns a sample from the 12th order
*       B-spline with cross-entropy -0.00023 relative to N(0,1).
*   (3) This package effectively uses integer arithmetic, but all
*       integers are less than 2**22.  Hence the routines will
*       (with most hardware) work identically in fully floating-point
*       arithmetic.
*   (4) Each individual generator limited to fairly small integers 
*       is rather crude, so the package uses six different generators.
*       The period is 15395410800000.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER ISEED

      INTEGER M1,M2,M3,M4,M5,M6
      PARAMETER (M1=6075,M2=7875,M3=53125,M4=11979,M5=6655,M6=31104)
      REAL P1,P2,P3,P4,P5,P6,Q1,Q2,Q3,Q4,Q5,Q6
      COMMON /VRANDC/ P1,P2,P3,P4,P5,P6
      SAVE Q1,Q2,Q3,Q4,Q5,Q6
      P1=MOD(ISEED+1,M1)
      P2=MOD(ISEED+2,M2)
      P3=MOD(ISEED+3,M3)
      P4=MOD(ISEED+4,M4)
      P5=MOD(ISEED+5,M5)
      P6=MOD(ISEED+6,M6)
      RETURN

      ENTRY VRAND1
      Q1=P1
      Q2=P2
      Q3=P3
      Q4=P4
      Q5=P5
      Q6=P6
      RETURN

      ENTRY VRAND2
      P1=Q1
      P2=Q2
      P3=Q3
      P4=Q4
      P5=Q5
      P6=Q6
      END

      SUBROUTINE VRAND(ST,JU,N)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*                 Normal random number generator
*        ST(JU),ST(JU+1),...,ST(JU+N-1) := samples from N(0,1)
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    JU      I     I      -         Start address
*    N       I     I      -         Length of block
* Local externals:
*    P1-P6   R       O    -         Random generator state
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*   (1) Requires previous call to VRAND0 to initialise the generator.
*   (2) Code is floating-point, and easily vectorisable.
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      INTEGER N,JU
      REAL ST(0:*)

      REAL P1,P2,P3,P4,P5,P6
      COMMON /VRANDC/ P1,P2,P3,P4,P5,P6

      INTEGER I,KU
      REAL A1,A2,A31,A32,A41,A42,A51,A52,A61,A62
      REAL C1,C2,C31,C32,C41,C42,C51,C52,C61,C62
      REAL X1,X2,X3,X4,X5,X6,R1,R2,R3,R4,R5,R6
      REAL ZERO,FRPT,XXX   
      PARAMETER (ZERO=0.0D0)
* In these parameter definitions for the six generators,
* X = individual period, as for M in VRAND0
* A = multiplier, factorised if necessary to avoid overflow
* C = offset, with 1/2 included to ensure that INT rounds correctly down
* R = 1/X, evaluated correctly as 1/sqrt(X*X-1) to get exact variance.
      PARAMETER (X1= 6075.0D0, R1= 1.646090557280774D-4)
      PARAMETER (A1=  106.0D0/ 6075.0D0, C1= 1283.5D0/ 6075.0D0)
      PARAMETER (X2= 7875.0D0, R2= 1.269841280079345D-4)
      PARAMETER (A2=  211.0D0/ 7875.0D0, C2= 1663.5D0/ 7875.0D0)
      PARAMETER (X3=53125.0D0, R3= 1.882352941509953D-5)
      PARAMETER (A31=   9.0D0/53125.0D0, C31=    0.5D0/53125.0D0)
      PARAMETER (A32=  19.0D0/53125.0D0, C32=11213.5D0/53125.0D0)
      PARAMETER (X4=11979.0D0, R4= 8.347942261327381D-5)
      PARAMETER (A41=  10.0D0/11979.0D0, C41=    0.5D0/11979.0D0)
      PARAMETER (A42=  43.0D0/11979.0D0, C42= 2531.5D0/11979.0D0)
      PARAMETER (X5= 6655.0D0, R5= 1.502629618767060D-4)
      PARAMETER (A51=  26.0D0/ 6655.0D0, C51=    0.5D0/ 6655.0D0)
      PARAMETER (A52=  36.0D0/ 6655.0D0, C52= 1399.5D0/ 6655.0D0)
      PARAMETER (X6=31104.0D0, R6= 3.215020577793267D-5)
      PARAMETER (A61=  25.0D0/31104.0D0, C61=    0.5D0/31104.0D0)
      PARAMETER (A62=  25.0D0/31104.0D0, C62= 6571.5D0/31104.0D0)
      FRPT(XXX)=XXX-FLOAT(INT(XXX))

      KU=JU
      DO 1 I=0,N-1
        ST(KU)=ZERO
        P1=INT(FRPT(A1*P1+C1)*X1)
        ST(KU)=ST(KU)+P1*R1
        P1=INT(FRPT(A1*P1+C1)*X1)
        ST(KU)=ST(KU)-P1*R1
        P2=INT(FRPT(A2*P2+C2)*X2)
        ST(KU)=ST(KU)+P2*R2
        P2=INT(FRPT(A2*P2+C2)*X2)
        ST(KU)=ST(KU)-P2*R2
        P3=INT(FRPT(A31*P3+C31)*X3)
        P3=INT(FRPT(A32*P3+C32)*X3)
        ST(KU)=ST(KU)+P3*R3
        P3=INT(FRPT(A31*P3+C31)*X3)
        P3=INT(FRPT(A32*P3+C32)*X3)
        ST(KU)=ST(KU)-P3*R3
        P4=INT(FRPT(A41*P4+C41)*X4)
        P4=INT(FRPT(A42*P4+C42)*X4)
        ST(KU)=ST(KU)+P4*R4
        P4=INT(FRPT(A41*P4+C41)*X4)
        P4=INT(FRPT(A42*P4+C42)*X4)
        ST(KU)=ST(KU)-P4*R4
        P5=INT(FRPT(A51*P5+C51)*X5)
        P5=INT(FRPT(A52*P5+C52)*X5)
        ST(KU)=ST(KU)+P5*R5
        P5=INT(FRPT(A51*P5+C51)*X5)
        P5=INT(FRPT(A52*P5+C52)*X5)
        ST(KU)=ST(KU)-P5*R5
        P6=INT(FRPT(A61*P6+C61)*X6)
        P6=INT(FRPT(A62*P6+C62)*X6)
        ST(KU)=ST(KU)+P6*R6
        P6=INT(FRPT(A61*P6+C61)*X6)
        P6=INT(FRPT(A62*P6+C62)*X6)
        ST(KU)=ST(KU)-P6*R6
        KU=KU+1
    1 CONTINUE
      END

************************************************************************
*
*    Storage addressing routines
*
************************************************************************

      SUBROUTINE VFETCH(ST,J,IOFF,ACTION, KCORE,LENGTH,MORE)
* with     ENTRY VSTORE()
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*         Prepare block of <J> and return address in ST
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I      *         Vector arithmetic workspace
*    J       I     I      -         Area number
*    IOFF    I     I      -         Relative address within area
*    ACTION  L     I      -         Flag to supply contents of block
*    KCORE   I       O    -         Core address            (VFETCH)
*                  I      -         Core address            (VSTORE)
*    LENGTH  I       O    -         Block length            (VFETCH)
*                  I      -         Block length            (VSTORE)
*    MORE    I     I(O)   -         Any more elements left? (VFETCH)
*                 (I)O    -         Any more elements left? (VSTORE)
* Saved:
*    NBUF    I     I O    -         Number of buffers
*    LAREA   I     I O    -         Dimension of current vector space
*    LBLOCK  I     I O    -         Local block size
* Globals:
*    KA      I     I      40        Allocation flags
*    KB      I     I      40        Base addresses in core
*    KL      I     I      40        Lengths
*    LWORK   I     I      -         Length of workspace in core
*    LENREC  I     I      -         Record length for disc storage
*
* External calls:
*    VSTACK      Stack manager
*    UFETCH      Read from external device
*    USTORE      Write to  external device
*
* History:
*    MKC/JS          1 Dec 1990     First release
*
* Notes:
*  (1)  This implementation uses a dummy run through the vector chain
*       in order to establish the optimal (largest) block size.  The
*       package always enters a new chain with MORE=0, so VFETCH can
*       use this as a flag to start counting the number of buffers it
*       needs (NBUF).  VFETCH immediately resets MORE=1 to avoid re-
*       starting the count, and LENGTH=0 to ensure that vector routines
*       do nothing cleanly.  Finally, VSTORE reads NBUF (once only) and
*       works out the allowed block size (LBLOCK) on the first (dummy)
*       time it is called.  This enables VFETCH to set LENGTH properly
*       to a positive value in subsequent loops in the chain.  These
*       then behave normally, picking up blocks of size up to LENGTH
*       until the areas are exhausted, when VSTORE returns MORE=0.
*  (2)  Local use of MORE is as follows.
*       MORE = 0  is start new chain or end old chain
*       MORE = 1  is counting buffers in VFETCH
*       MORE = 2  is setting LBLOCK in VSTORE
*       MORE = 3  is rolling along chain
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER (A-Z)
      LOGICAL ACTION
      INTEGER J,IOFF,LENGTH,KCORE,MORE
      REAL ST(0:*)

      INTEGER KA,KB,KL,KWORK,LWORK,LENREC
      COMMON /VPOINT/ KL(40),KB(40),KA(40),KWORK,LWORK,LENREC

      INTEGER NBUF,LAREA,LBLOCK
      SAVE    NBUF,LAREA,LBLOCK
	
      IF (MORE.EQ.0) THEN
* Initialise and count disc buffers if using dynamic block sizes .....
        MORE=1
        NBUF=0
        LAREA=KL(J)
        LENGTH=0
      ENDIF
      IF (MORE.EQ.1) THEN
        IF (KL(J).NE.LAREA) STOP ' Inconsistent lengths'
        IF (KA(J).EQ.1)  NBUF=NBUF+1
      ELSE
        MORE=3
* ....................................................................
* Else set up as much as possible of LBLOCK
        LENGTH=MIN(KL(J)-IOFF,LBLOCK)
        IF (KA(J).EQ.0) THEN
*   held in core
          KCORE=KB(J)+IOFF
        ELSEIF (KA(J).EQ.1) THEN
*   held on disc
          CALL VSTACK(1,KCORE)
          IF (ACTION) CALL UFETCH(ST,KCORE,KB(J)+IOFF,LENGTH)
        ENDIF
      ENDIF
      RETURN

      ENTRY VSTORE(ST,J,IOFF,ACTION,KCORE,LENGTH,MORE)
* If using dynamic block sizes instead of fixed LBLOCK ...............
      IF (MORE.EQ.1) THEN
* Calculate block size, and set up stack
        LBLOCK=KL(J)
        IF (NBUF.GT.0) THEN
          LBLOCK=MIN(LBLOCK,LWORK/NBUF)
          IF (LBLOCK.LT.KL(J)) LBLOCK=LENREC*(LBLOCK/LENREC)
          IF (LBLOCK.LE.0) STOP ' Not enough space for buffers'
          CALL VSTACK(0,LBLOCK)
        ENDIF
        MORE=2
      ENDIF
      IF(MORE.NE.2) THEN
* ....................................................................
*   if held in core, no operation
        IF (KA(J).EQ.1) THEN
*   if held on disc
          IF (ACTION) CALL USTORE(ST,KCORE,KB(J)+IOFF,LENGTH)
          CALL VSTACK(-1,KCORE)
        ENDIF
      ENDIF
* Any more elements?
      IF (IOFF+LENGTH.GE.KL(J)) MORE=0
      END

      SUBROUTINE VSTACK(I,K)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*         Manage stack of pointers
*
*         I > 0  :  Pop pointer
*         I = 0  :  Set pointers
*         I < 0  :  Push pointer
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    I       I     I      -         Required stack operation
*    K       I     I      -         Pointer (I>0)
*                  I      -         Buffer size (I=0)
*                    O    -         Pointer (I<0)
* Saved:
*    KSTACK  I     I O    4         Pointer list
*    ISTACK  I     I O    -         Position in list
* Globals:
*    KWORK   I     I      -         Start of workspace in core
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
      INTEGER I,K

      INTEGER KA,KB,KL,KWORK,LWORK,LENREC
      COMMON /VPOINT/ KL(40),KB(40),KA(40),KWORK,LWORK,LENREC

      INTEGER MAXBUF
      PARAMETER (MAXBUF=5)
      INTEGER KSTACK(MAXBUF),ISTACK,J
      SAVE KSTACK,ISTACK

      IF (I.GT.0) THEN
* Pop a buffer
        IF (ISTACK.GE.MAXBUF) STOP ' Stack overflow'
        ISTACK=ISTACK+1
        K=KSTACK(ISTACK)
      ELSEIF (I.EQ.0) THEN
* Initialise stack
        DO 1 J=1,MAXBUF
          KSTACK(J)=KWORK+(J-1)*K
    1   CONTINUE
        ISTACK=0
      ELSE
* Push a buffer
        IF (ISTACK.LE.0) STOP ' Stack underflow'
        KSTACK(ISTACK)=K
        ISTACK=ISTACK-1
      ENDIF
      END
