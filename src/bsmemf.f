* Copyright (C) 1990,1991 Maximum Entropy Data Consultants Ltd

      SUBROUTINE UAREA(L,M,N,NSTORE,LEVEL1,KBC)
      INTEGER KBC(*)
* Set up storage.
      COMMON /OUTPUT/ IOUT,LEVEL
      COMMON /VPOINT/ KL(40),KB(40),KA(40),KWORK,LWORK,LENREC

      IOUT=6
	LEVEL=LEVEL1
** (L,M,N) = length of (hidden,visible,data) space
* Record length for disc I/O
      LENREC=0
*  Set up transform areas <2>, <11>, <25> in core (<2> and <11> 
*  overlaid), and the rest on disc
      KCORE=0
      KDISC=0
      KA(1)=0
      KB(1)=KCORE
      KL(1)=L
      KCORE=KCORE+L

      KA(2)=0
      KB(2)=KCORE
      KL(2)=L
      KCORE=KCORE+MAX(L,M)

      KA(3)=0
      KB(3)=KCORE
      KL(3)=L
      KCORE=KCORE+L

      KA(4)=0
      KB(4)=KCORE
      KL(4)=L
      KCORE=KCORE+L

      KA(5)=0
      KB(5)=KCORE
      KL(5)=L
      KCORE=KCORE+L

      KA(11)=0
      KB(11)=KB(2)
      KL(11)=M

      KA(21)=0
      KB(21)=KCORE
      KL(21)=N
      KCORE=KCORE+N

      KA(22)=0
      KB(22)=KCORE
      KL(22)=N
      KCORE=KCORE+N

      KA(23)=0
      KB(23)=KCORE
      KL(23)=N
      KCORE=KCORE+N

      KA(24)=0
      KB(24)=KCORE
      KL(24)=N
      KCORE=KCORE+N

      KA(25)=0
      KB(25)=KCORE
      KL(25)=N
      KCORE=KCORE+N

      KA(26)=0
      KB(26)=KCORE
      KL(26)=N
      KCORE=KCORE+N

      KA(27)=0
      KB(27)=KCORE
      KL(27)=N
      KCORE=KCORE+N

      KA(28)=0
      KB(28)=KCORE
      KL(28)=N
	DO I=1,40
        KBC(I) = KB(I)
      ENDDO

      IF (KCORE.GT.NSTORE) STOP ' Not enough Storage available'
      KWORK=KCORE
      LWORK=NSTORE-KCORE
      END

      SUBROUTINE SETMSK(ST,J,K)
* Set <11> := square wave from cell j to k
      REAL ST(0:*)

      COMMON /VPOINT/ KL(40),KB(40),KA(40),KWORK,LWORK,LENREC

      KC=KB(11)
      M=KL(11)
      DO 1 I=0,M-1
        IF( I.LT.J-1 .OR. K-1.LT.I ) THEN
          ST(KC+I)=0.
        ELSE
          ST(KC+I)=1.
        ENDIF
    1 CONTINUE
      END

      SUBROUTINE UINFO(STRING,MLEVEL)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*   Manage numerical and progress diagnostics output
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    STRING  C *   I      -         Diagnostic string
*    MLEVEL  I     I      -         Composite diagnostic type
* Globals:
*    LEVEL   I     I      -         Composite switch for diagnostics
*    IOUT    I     I      -         Output channel
*
* External calls:
*     -
*
* MLEVEL is interpreted as:-
*   10 = Major numerical diagnostic   1 = Name of major subroutine
*   20 = Technical numerical          2 = Flowchart (Read/Write/Update)
*
* History:
*    MKC/JS           1 Dec 1990     First release
*
* Notes:
*  (1)  In this example, LEVEL is a 2-decimal-digit integer like MLEVEL
*       Its "tens"  column defines the  tens(MLEVEL) above which the
*                                 numerical diagnostics will be ignored.
*       Its "units" column defines the units(MLEVEL) above which the
*                                  progress diagnostics will be ignored.
*   Example: LEVEL=12 writes according to the following values of MLEVEL
*   10 = Numerical    YES             1 = Routine     YES
*   20 = Technical    no              2 = Flowchart   YES
*-----------------------------------------------------------------------
      CHARACTER*(*) STRING
      COMMON /OUTPUT/ IOUT,LEVEL
      IF( MLEVEL .GE. 10 ) THEN
        IF( MLEVEL .GT. LEVEL ) RETURN
      ELSE
        IF( MOD(MLEVEL,10) .GT. MOD(LEVEL,10) ) RETURN
      ENDIF
*      WRITE(IOUT,*) STRING
      END

      SUBROUTINE UDIAG(ST,STRING)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Purpose:
*   Diagnostic routine, to inspect any part of any area in core.
*
*     NAME   TYPE  I/O  DIMENSION   DESCRIPTION
* Arguments:
*    ST      R     I O    *         Vector arithmetic workspace
*    STRING  C *   I      -         Message to show calling location
* Globals:
*    KA      I     I      40        Allocation flags
*    KB      I     I      40        Base addresses in core
*
* External calls:
*    MESWAP      Interchange areas
*
* History:
*    MKC/JS           1 Nov 1990     First release
*        JS          20 Feb 1991     Sum incorporated
*
* Notes:
*    (1)  It can be used to investigate areas held on disc, by
*         swapping them with <2>, <11> and <25> (assumed to be in core)
*         for hidden, visible and data space areas respectively.
*         The previous contents of these core areas are restored before
*         returning.
*    (2)  This means that UDIAG for external areas should ONLY be used
*         when all the areas have been completely read or written,
*         and NOT part way through a vector chain.
*    (3)  This implementation is numerical, and interrogates the user
*         via the default Fortran stream (*).
*-----------------------------------------------------------------------
      REAL ST(0:*)
      CHARACTER*(*) STRING

      COMMON /VPOINT/ KL(40),KB(40),KA(40),KWORK,LWORK,LENREC

      WRITE(*,*) ' UDIAG:- ',STRING
    1 CONTINUE
        WRITE(*,*) ' Which area would you like?  (0 exits)  '
        READ(*,*,END=1,ERR=1) J
        IF (J.LT.0 .OR. J.GT.40) GOTO 1
        IF (J.EQ.0) RETURN
        IF (KA(J).EQ.1) THEN
          WRITE (*,'(''      (swapping external area in)  '')')
          IF (J.LE.10) THEN
            K=2
          ELSEIF (J.LE.20) THEN
            K=11
          ELSEIF (J.LE.30) THEN
            K=25
          ELSE
            WRITE (*,*) ' Unknown area.'
            GOTO 1
          ENDIF
          CALL MESWAP(ST,J,K)
	ELSE
          WRITE (*,'(''               (core area)         '',
     *               '' KB ='',I7,'' = base address'')') KB(J)
          K=J
        ENDIF
        WRITE(*,*) ' First element (from 1) ?  '
        READ(*,*) I1
        WRITE(*,*) '  Last element ?  '
        READ(*,*) I2
        I1=KB(K)+I1-1
        I2=KB(K)+I2-1
        SUM=0.
        DO 2 I=I1,I2
          SUM=SUM+ST(I)
    2   CONTINUE
        WRITE(*,'(''  Sum ='',1PE13.4)') SUM
        WRITE(*,'(1PE14.4,3E14.4)') (ST(I),I=I1,I2)
        IF(KA(J).EQ.1) THEN
          WRITE (*,'(''      (swapping external area out) '')')
          CALL MESWAP(ST,J,K)
        ENDIF
      GOTO 1
      END

      SUBROUTINE UFETCH(ST,KCORE,KDISC,LENGTH)
      REAL ST(0:*)
      STOP ' Disc Storage no longer implemented'
      END

      SUBROUTINE USTORE(ST,KCORE,KDISC,LENGTH)
      REAL ST(0:*)
      STOP ' Disc Storage no longer implemented'
      END


***********************************************************************
*  Application-specific OPUS and ICF operators, with their transposes.
***********************************************************************

      SUBROUTINE OPUS(ST,JM, JN)
* Forward differential operator:    <JN> := Opus<JM>
      REAL ST(0:*)
      COMMON /VPOINT/ KL(40),KB(40),KA(40),KWORK,LWORK,LENREC
      CALL VOPUS(ST(KB(JM)),ST(KB(JN)))
      END

      SUBROUTINE TROPUS(ST,JN, JM)
* Transpose differential operator:  <JM> := Tropus<JN>
      REAL ST(0:*)
      COMMON /VPOINT/ KL(40),KB(40),KA(40),KWORK,LWORK,LENREC
      CALL VTROP(ST(KB(JM)),ST(KB(JN)))
      END

      SUBROUTINE ICF(ST,JL, JM)
* Forward intrinsic correlation function operator:    <JM> := ICF<JL>
      REAL ST(0:*)
      CALL MECOPY(ST,JL,JM)
      END

      SUBROUTINE TRICF(ST,JM, JL)
* Transpose intrinsic correlation function operator:  <JL> := TrICF<JM>
      REAL ST(0:*)
      CALL MECOPY(ST,JM,JL)
      END

***********************************************************************
*  Application-specific GRADL, GRDGRL operators
***********************************************************************

      SUBROUTINE MEMEX(ST,JM, JN)
* Model the ACTUAL behaviour of the experiment: <JN> := MemEx<JM>
* LINEAR experiments are treated by a call to the
*  differential response routine VOPUS.
* NONLINEAR experiments should transform a distribution on file <JM>
*  to mock data on file <JN>.
      REAL ST(0:*)

      COMMON /VPOINT/ KL(40),KB(40),KA(40),KWORK,LWORK,LENREC

      CALL VMEMEX(ST(KB(JM)),ST(KB(JN)))
      END

