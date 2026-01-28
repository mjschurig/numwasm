      SUBROUTINE SODEX(N,FCN,IFCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JAC ,IJAC,MLJAC,MUJAC,DFX,IDFX,
     &                  MAS ,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,IDID)
C ----------------------------------------------------------
C     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
C     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y).
C     THIS IS AN EXTRAPOLATION-ALGORITHM, BASED ON THE
C     LINEARLY IMPLICIT MID-POINT RULE, DUE TO BADER-DEUFLHARD
C     (WITH STEP SIZE CONTROL AND ORDER SELECTION).
C     C.F. SECTION IV.9 
C
C     AUTHORS: E. HAIRER AND G. WANNER
C              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
C              CH-1211 GENEVE 24, SWITZERLAND 
C              E-MAIL:  Ernst.Hairer@math.unige.ch
C                       Gerhard.Wanner@math.unige.ch
C     
C     THIS CODE IS PART OF THE BOOK:
C         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
C         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
C         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
C         SPRINGER-VERLAG (1990)               
C      
C     VERSION OF NOVEMBER 17, 1992
C         SMALL CORRECTIONS ON JUNE 11, 1999
C
C     INPUT PARAMETERS  
C     ----------------  
C     N           DIMENSION OF THE SYSTEM 
C
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
C                 VALUE OF F(X,Y):
C                    SUBROUTINE FCN(N,X,Y,F)
C                    REAL*8 X,Y(N),F(N)
C                    F(1)=...   ETC.
C
C     IFCN        GIVES INFORMATION ON FCN:
C                    IFCN=0: F(X,Y) INDEPENDENT OF X (AUTONOMOUS)
C                    IFCN=1: F(X,Y) MAY DEPEND ON X (NON-AUTONOMOUS)
C
C     X           INITIAL X-VALUE
C
C     Y(N)        INITIAL VALUES FOR Y
C
C     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
C
C     H           INITIAL STEP SIZE GUESS;
C                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT, 
C                 H=1.D0/(NORM OF F'), USUALLY 1.D-2 OR 1.D-3, IS GOOD.
C                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
C                 ADAPTS ITS STEP SIZE. STUDY THE CHOSEN VALUES FOR A FEW
C                 STEPS IN SUBROUTINE "SOLOUT", WHEN YOU ARE NOT SURE.
C                 (IF H=0.D0, THE CODE PUTS H=1.D-6).
C
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
C                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
C
C     ITOL        SWITCH FOR RTOL AND ATOL:
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
C                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
C                     RTOL(I)*ABS(Y(I))+ATOL(I).
C
C     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
C                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
C                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
C                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
C                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
C                    SUBROUTINE JAC(N,X,Y,DFY,LDFY)
C                    REAL*8 X,Y(N),DFY(LDFY,N)
C                    DFY(1,1)= ...
C                 LDFY, THE COLOMN-LENGTH OF THE ARRAY, IS
C                 FURNISHED BY THE CALLING PROGRAM.
C                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO
C                    BE FULL AND THE PARTIAL DERIVATIVES ARE
C                    STORED IN DFY AS
C                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
C                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
C                    THE PARTIAL DERIVATIVES ARE STORED
C                    DIAGONAL-WISE AS
C                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
C
C     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
C                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
C                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
C                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
C
C     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
C                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
C                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN 
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
C                       THE MAIN DIAGONAL).
C
C     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
C                 NEED NOT BE DEFINED IF MLJAC=N.
C
C     DFX         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
C                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO X
C                 (THIS ROUTINE IS ONLY CALLED IF IDFX=1 AND IFCN=1;
C                 SUPPLY A DUMMY SUBROUTINE IN THE CASE IDFX=0 OR IFCN=0).
C                 OTHERWISE, THIS SUBROUTINE MUST HAVE THE FORM
C                    SUBROUTINE DFX(N,X,Y,FX)
C                    REAL*8 X,Y(N),FX(N)
C                    FX(1)= ...
C                
C     IDFX        SWITCH FOR THE COMPUTATION OF THE DF/DX:
C                    IDFX=0: DF/DX IS COMPUTED INTERNALLY BY FINITE
C                       DIFFERENCES, SUBROUTINE "DFX" IS NEVER CALLED.
C                    IDFX=1: DF/DX IS SUPPLIED BY SUBROUTINE DFX.
C
C     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      -----
C     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
C
C     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
C                 MATRIX M.
C                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
C                 MATRIX AND NEEDS NOT TO BE DEFINED;
C                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
C                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
C                    SUBROUTINE MAS(N,AM,LMAS)
C                    REAL*8 AM(LMAS,N)
C                    AM(1,1)= ....
C                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
C                    AS FULL MATRIX LIKE
C                         AM(I,J) = M(I,J)
C                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
C                    DIAGONAL-WISE AS
C                         AM(I-J+MUMAS+1,J) = M(I,J).
C
C     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
C                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
C                       MATRIX, MAS IS NEVER CALLED.
C                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
C
C     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
C                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
C                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
C                       THE MAIN DIAGONAL).
C                 MLMAS IS SUPPOSED TO BE .LE. MLJAC.
C
C     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON-
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
C                 NEED NOT BE DEFINED IF MLMAS=N.
C                 MUMAS IS SUPPOSED TO BE .LE. MUJAC.
C
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION. 
C                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
C                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. 
C                 IT MUST HAVE THE FORM
C                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,IRTRN)
C                    REAL*8 X,Y(N)
C                    ....  
C                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
C                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
C                    THE FIRST GRID-POINT).
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, SODEX RETURNS TO THE CALLING PROGRAM.
C           
C     IOUT        GIVES INFORMATION ON THE SUBROUTINE SOLOUT:
C                    IOUT=0: SUBROUTINE IS NEVER CALLED
C                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
C
C     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
C                 SERVES AS WORKING SPACE FOR ALL VECTORS AND MATRICES.
C                 "LWORK" MUST BE AT LEAST
C                             N*(LJAC+LMAS+LE1+KM+9)+3*KM+13
C                 WHERE
C                    LJAC=N              IF MLJAC=N (FULL JACOBIAN)
C                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.)
C                 AND                  
C                    LMAS=0              IF IMAS=0
C                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL)
C                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.)
C                 AND
C                    LE1=N               IF MLJAC=N (FULL JACOBIAN)
C                    LE1=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.).
C                 AND                                 
C                    KM=6                IF IWORK(3)=0
C                    KM=IWORK(3)         IF IWORK(3).GT.0
C
C                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE
C                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM
C                 STORAGE REQUIREMENT IS 
C                             LWORK = 2*N*N+(KM+9)*N+3*KM+13.
C
C     LWORK       DECLARED LENGHT OF ARRAY "WORK".
C
C     IWORK       INTEGER WORKING SPACE OF LENGHT "LIWORK".
C                 "LIWORK" MUST BE AT LEAST  2*N+KM+4.
C
C     LIWORK      DECLARED LENGHT OF ARRAY "IWORK".
C
C ----------------------------------------------------------------------
C 
C     SOPHISTICATED SETTING OF PARAMETERS
C     -----------------------------------
C              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK 
C              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),..,WORK(13)
C              AS WELL AS IWORK(1),..,IWORK(4) DIFFERENT FROM ZERO.
C              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
C
C    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN
C              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY
C              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN.
C              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N)
C              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1). IT IS
C              ALSO NOT GOOD FOR SPARSE JACOBIANS.
C
C    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
C
C    IWORK(3)  THE MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION 
C              TABLE. THE DEFAULT VALUE (FOR IWORK(3)=0) IS 6.
C              IF IWORK(3).NE.0 THEN IWORK(3) SHOULD BE .GE.3.
C
C    IWORK(4)  SWITCH FOR THE STEP SIZE SEQUENCE
C              IF IWORK(4).EQ.1 THEN  2,6,10,14,22,34,50,...
C              THE DEFAULT VALUE (FOR IWORK(4)=0) IS IWORK(4)=1.
C
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT  1.D-16.
C
C    WORK(2)   MAXIMAL STEP SIZE, DEFAULT  XEND-X.
C
C    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
C              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS
C              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER 
C              (0.001D0, SAY).     DEFAULT RTOL(1).
C
C    WORK(4), WORK(5)   PARAMETERS FOR STEP SIZE SELECTION
C              THE NEW STEP SIZE FOR THE J-TH DIAGONAL ENTRY IS
C              CHOSEN SUBJECT TO THE RESTRICTION
C                 FACMIN/WORK(5) <= HNEW(J)/HOLD <= 1/FACMIN
C              WHERE FACMIN=WORK(4)**(1/(J-1)) 
C              DEFAULT VALUES: WORK(4)=0.1D0, WORK(5)=4.D0
C
C    WORK(6), WORK(7)   PARAMETERS FOR THE ORDER SELECTION
C              STEP SIZE IS DECREASED IF    W(K-1) <= W(K)*WORK(6)
C              STEP SIZE IS INCREASED IF    W(K) <= W(K-1)*WORK(7)
C              DEFAULT VALUES: WORK(6)=0.9D0, WORK(7)=0.9D0
C
C    WORK(8), WORK(9)   SAFETY FACTORS FOR STEP CONTROL ALGORITHM
C             HNEW=H*WORK(9)*(WORK(8)*TOL/ERR)**(1/(J-1))
C              DEFAULT VALUES: WORK(8)=0.8D0, WORK(9)=0.93D0
C
C    WORK(10), WORK(11), WORK(12), WORK(13)   ESTIMATED WORKS FOR
C             A CALL TO FCN, JAC, DEC, SOL, RESPECTIVELY.
C             DEFAULT VALUES ARE: WORK(10)=1.D0, WORK(11)=5.D0,
C             WORK(12)=1.D0, WORK(13)=1.D0.
C
C-----------------------------------------------------------------------
C
C     OUTPUT PARAMETERS 
C     ----------------- 
C     X           X-VALUE WHERE THE SOLUTION IS COMPUTED
C                 (AFTER SUCCESSFUL RETURN X=XEND)
C
C     Y(N)        SOLUTION AT X
C  
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
C
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
C                   IDID=1  COMPUTATION SUCCESSFUL,
C                   IDID=-1 COMPUTATION UNSUCCESSFUL.
C
C --------------------------------------------------------- 
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C          DECLARATIONS 
C *** *** *** *** *** *** *** *** *** *** *** *** ***
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(N),ATOL(1),RTOL(1),WORK(LWORK),IWORK(LIWORK)
      LOGICAL AUTNMS,IMPLCT,ARRET,JBAND
      EXTERNAL FCN,JAC,DFX,MAS,SOLOUT
      COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
C --- COMMON STAT CAN BE USED FOR STATISTICS
C ---    NFCN      NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
C                  EVALUATION OF THE JACOBIAN ARE NOT COUNTED)  
C ---    NJAC      NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
C                  OR NUMERICALLY), COMPUTATION OF DFX IS INCLUDED
C ---    NSTEP     NUMBER OF COMPUTED STEPS
C ---    NACCPT    NUMBER OF ACCEPTED STEPS
C ---    NREJCT    NUMBER OF REJECTED STEPS (AFTER AT LEAST ONE STEP
C                  HAS BEEN ACCEPTED)
C ---    NDEC      NUMBER OF LU-DECOMPOSITIONS (N-DIMENSIONAL MATRIX)
C ---    NSOL      NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS
C *** *** *** *** *** *** ***
C        SETTING THE PARAMETERS 
C *** *** *** *** *** *** ***
       NFCN=0
       NJAC=0
       NSTEP=0
       NACCPT=0
       NREJCT=0
       NDEC=0
       NSOL=0
      ARRET=.FALSE.
C -------- SWITCH FOR TRANSFORMATION OF JACOBIAN TO HESSIAN FORM ---
      NHESS=IWORK(1)
      IF (N.LE.2) NHESS=0
C -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
      IF(IWORK(2).EQ.0)THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(2)
         IF(NMAX.LE.0)THEN
            WRITE(6,*)' WRONG INPUT IWORK(2)=',IWORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C -------- KM     MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION 
      IF(IWORK(3).EQ.0)THEN
         KM=6
      ELSE
         KM=IWORK(3)
         IF(KM.LE.2)THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(3)=',IWORK(3)
            ARRET=.TRUE.
         END IF
      END IF
C -------- NSEQU     MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION 
      NSEQU=IWORK(4)
      IF(IWORK(4).EQ.0) NSEQU=1
      IF(NSEQU.LE.0.OR.NSEQU.GE.2)THEN
         WRITE(6,*)' CURIOUS INPUT IWORK(4)=',IWORK(4)
         ARRET=.TRUE.
      END IF
C -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0  
      IF(WORK(1).EQ.0.D0)THEN
         UROUND=1.D-16
      ELSE
         UROUND=WORK(1)
         IF(UROUND.LE.1.D-19.OR.UROUND.GE.1.D0)THEN
            WRITE(6,*)' COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1)
            ARRET=.TRUE.
         END IF
      END IF
C -------- MAXIMAL STEP SIZE
      IF(WORK(2).EQ.0.D0)THEN
         HMAX=XEND-X
      ELSE
         HMAX=WORK(2)
      END IF
C ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
      IF(WORK(3).EQ.0.D0)THEN
         THET=RTOL(1)
      ELSE
         THET=WORK(3)
      END IF
C -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
      IF(WORK(4).EQ.0.D0)THEN
         FAC1=0.1D0
      ELSE
         FAC1=WORK(4)
      END IF
      IF(WORK(5).EQ.0.D0)THEN
         FAC2=4.0D0
      ELSE
         FAC2=WORK(5)
      END IF
C -------  FAC3, FAC4   PARAMETERS FOR THE ORDER SELECTION
      IF(WORK(6).EQ.0.D0)THEN
         FAC3=0.9D0
      ELSE
         FAC3=WORK(6)
      END IF
      IF(WORK(7).EQ.0.D0)THEN
         FAC4=0.9D0
      ELSE
         FAC4=WORK(7)
      END IF
C ------- SAFE1, SAFE2 SAFETY FACTORS FOR STEP SIZE PREDICTION
      IF(WORK(8).EQ.0.D0)THEN
         SAFE1=0.8D0
      ELSE
         SAFE1=WORK(8)
      END IF
      IF(WORK(9).EQ.0.D0)THEN
         SAFE2=0.93D0
      ELSE
         SAFE2=WORK(9)
      END IF
C ------- WKFCN,WKJAC,WKDEC,WKSOL  ESTIMATED WORK FOR  FCN,JAC,DEC,SOL
      IF(WORK(10).EQ.0.D0)THEN
         WKFCN=1.D0
      ELSE
         WKFCN=WORK(10)
      END IF
      IF(WORK(11).EQ.0.D0)THEN
         WKJAC=5.D0
      ELSE
         WKJAC=WORK(11)
      END IF
      IF(WORK(12).EQ.0.D0)THEN
         WKDEC=1.D0
      ELSE
         WKDEC=WORK(12)
      END IF
      IF(WORK(13).EQ.0.D0)THEN
         WKSOL=1.D0
      ELSE
         WKSOL=WORK(13)
      END IF
      WKROW=WKFCN+WKSOL
C --------- CHECK IF TOLERANCES ARE O.K.
      IF (ITOL.EQ.0) THEN
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES ARE TOO SMALL'
              ARRET=.TRUE.
          END IF
      ELSE
          DO 15 I=1,N
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL'
              ARRET=.TRUE.
          END IF
  15      CONTINUE
      END IF
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C         COMPUTATION OF ARRAY ENTRIES
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C ---- AUTONOMOUS, IMPLICIT, BANDED OR NOT ?
      AUTNMS=IFCN.EQ.0
      IMPLCT=IMAS.NE.0
      JBAND=MLJAC.NE.N
      ARRET=.FALSE.
C -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS ---
C -- JACOBIAN 
      IF(JBAND)THEN
         LDJAC=MLJAC+MUJAC+1
      ELSE
         LDJAC=N
      END IF
C -- MATRIX E FOR LINEAR ALGEBRA
      IF(JBAND)THEN
         LDE=2*MLJAC+MUJAC+1
      ELSE
         LDE=N
      END IF
C -- MASS MATRIX
      IF (IMPLCT) THEN
          IF (MLMAS.NE.N) THEN
              LDMAS=MLMAS+MUMAS+1
          ELSE
              LDMAS=N
          END IF
C ------ BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC"
          IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN
              WRITE (6,*) 'BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF
     & "JAC"'
            ARRET=.TRUE.
          END IF
      ELSE
          LDMAS=0
      END IF
      LDMAS2=MAX(1,LDMAS)
C ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN
      IF ((IMPLCT.OR.JBAND).AND.NHESS.NE.0) THEN
         WRITE(6,*)' HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS WITH 
     &FULL JACOBIAN'
         ARRET=.TRUE.
      END IF
C ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
      IEYH=14
      IEDY=IEYH+N
      IEFX=IEDY+N
      IEYHH=IEFX+N
      IEDYH=IEYHH+N
      IEDEL=IEDYH+N
      IEFXH=IEDEL+N
      IEWH =IEFXH+N
      IESCAL=IEWH+N
      IEHH =IESCAL+N
      IEW  =IEHH+KM
      IEA  =IEW+KM
      IEJAC=IEA+KM
      IEE  =IEJAC+N*LDJAC
      IEMAS=IEE+N*LDE
      IET=IEMAS+N*LDMAS
C ------ TOTAL STORAGE REQUIREMENT -----------
      ISTORE=IET+N*KM-1
      IF(ISTORE.GT.LWORK)THEN
         WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------- ENTRY POINTS FOR INTEGER WORKSPACE -----
      IEIP=5
      IENJ=IEIP+N
      IEIPH=IENJ+KM
C --------- TOTAL REQUIREMENT ---------------
      ISTORE=IEIPH+N-1
      IF(ISTORE.GT.LIWORK)THEN
         WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF
C -------- CALL TO CORE INTEGRATOR ------------
      CALL SODCOR(N,FCN,X,Y,XEND,HMAX,H,KM,RTOL,ATOL,ITOL,JAC,IJAC,
     &   MLJAC,MUJAC,DFX,IDFX,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,
     &   NMAX,UROUND,NSEQU,NHESS,AUTNMS,IMPLCT,JBAND,LDJAC,LDE,LDMAS2,
     &   WORK(IEYH),WORK(IEDY),WORK(IEFX),WORK(IEYHH),WORK(IEDYH),
     &   WORK(IEDEL),WORK(IEFXH),WORK(IEWH),WORK(IESCAL),WORK(IEHH),
     &   WORK(IEW),WORK(IEA),WORK(IEJAC),WORK(IEE),WORK(IEMAS),
     &   WORK(IET),IWORK(IEIP),IWORK(IENJ),IWORK(IEIPH),FAC1,FAC2,FAC3,
     &   FAC4,THET,SAFE1,SAFE2,WKJAC,WKDEC,WKROW)
C ----------- RETURN -----------
      RETURN
      END
C
C
C
C  ----- ... AND HERE IS THE CORE INTEGRATOR  ----------
C
      SUBROUTINE SODCOR(N,FCN,X,Y,XEND,HMAX,H,KM,RTOL,ATOL,ITOL,
     &  JAC,IJAC,MLJAC,MUJAC,DFX,IDFX,MAS,MLB,MUB,SOLOUT,IOUT,IDID,
     &  NMAX,UROUND,NSEQU,NHESS,AUTNMS,IMPLCT,BANDED,LFJAC,LE,
     &  LDMAS,YH,DY,FX,YHH,DYH,DEL,FXH,WH,SCAL,HH,W,A,FJAC,E,FMAS,T,IP,
     &  NJ,IPHES,FAC1,FAC2,FAC3,FAC4,THET,SAFE1,SAFE2,WKJAC,WKDEC,WKROW)
C ----------------------------------------------------------
C     CORE INTEGRATOR FOR SODEX
C     PARAMETERS SAME AS IN SODEX WITH WORKSPACE ADDED 
C ---------------------------------------------------------- 
C         DECLARATIONS 
C ---------------------------------------------------------- 
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION Y(N),YH(N),DY(N),FX(N),YHH(N),DYH(N),DEL(N),FXH(N)
       DIMENSION WH(N),SCAL(N),HH(KM),W(KM),A(KM),FJAC(LFJAC,N),E(LE,N)
       DIMENSION FMAS(LDMAS,N),T(KM,N),IP(N),NJ(KM),RTOL(1),ATOL(1)
       DIMENSION IPHES(N)
       LOGICAL REJECT,LAST,ATOV,CALJAC,AUTNMS,IMPLCT,BANDED
       EXTERNAL FCN
       COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
C ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
      IF (IMPLCT) THEN
        CALL MAS(N,FMAS,LDMAS)
        MBDIAG=MUB+1
      END IF
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** ***
       IF (BANDED) THEN
          MLE=MLJAC
          MUE=MUJAC
          MBJAC=MLJAC+MUJAC+1
          MBB=MLB+MUB+1
          MDIAG=MLE+MUE+1
          MDIFF=MLE+MUE-MUB
       END IF
C --- DEFINE THE STEP SIZE SEQUENCE
       IF (NSEQU.EQ.1) THEN
           NJ(1)=2
           NJ(2)=6
           IF (KM.GE.3) NJ(3)=10
           IF (KM.GE.4) NJ(4)=14
           IF (KM.GE.5) NJ(5)=22
           IF (KM.GE.6) NJ(6)=34
           IF (KM.GE.7) NJ(7)=50
           IF (KM.GE.8) THEN
              DO 1 I=8,KM
   1          NJ(I)=2*NJ(I-2)+NJ(1) 
           END IF
       END IF
       A(1)=WKJAC+(NJ(1)+1)*WKROW+WKDEC
       DO 4 I=2,KM
   4   A(I)=A(I-1)+NJ(I)*WKROW+WKDEC
       POSNEG=SIGN(1.D0,XEND-X) 
       K=MAX0(2,MIN0(KM-1,INT(-DLOG10(RTOL(1)+ATOL(1))*.6D0+1.5D0))) 
       HMAXN=MIN(ABS(HMAX),ABS(XEND-X)) 
       H=MAX(ABS(H),1.D-6) 
       H=POSNEG*MIN(H,HMAXN)
       THETA=2*ABS(THET)
       IF (IOUT.NE.0) THEN
          IRTRN=1
          NRSOL=1
          XSOL=X 
          XOLD=X
          DO 7 I=1,N
  7       YH(I)=Y(I)
          NSOLU=N
          CALL SOLOUT(NRSOL,XOLD,XSOL,YH,NSOLU,IRTRN)
          IF (IRTRN.LT.0) GOTO 120
       END IF
       ERR=0.D0
       W(1)=1.D30  
       DO 8 I=1,N
       IF (ITOL.EQ.0) THEN
         SCAL(I)=ATOL(1)+RTOL(1)*DABS(Y(I))
       ELSE
         SCAL(I)=ATOL(I)+RTOL(I)*DABS(Y(I))
       END IF
   8   CONTINUE
       CALJAC=.FALSE.
       REJECT=.FALSE.
       LAST=.FALSE.
  10   CONTINUE
       IF (REJECT) THETA=2*ABS(THET)
       ATOV=.FALSE.
C *** *** *** *** *** *** ***
C --- IS XEND REACHED IN THE NEXT STEP?
C *** *** *** *** *** *** ***
       IF (0.1D0*ABS(XEND-X).LE.ABS(X)*UROUND) GO TO 110
       HOPT=H
       H=POSNEG*MIN(ABS(H),ABS(XEND-X),HMAXN)
       IF ((X+1.01D0*H-XEND)*POSNEG.GT.0.D0) THEN
          H=XEND-X 
          LAST=.TRUE.
       END IF
       CALL FCN(N,X,Y,DY)
       NFCN=NFCN+1
       IF (THETA.GT.THET.AND..NOT.CALJAC) THEN
C *** *** *** *** *** *** ***
C  COMPUTATION OF THE JACOBIAN
C *** *** *** *** *** *** ***
          NJAC=NJAC+1
          IF (IJAC.EQ.0) THEN
C --- COMPUTE JACOBIAN MATRIX NUMERICALLY
              IF (BANDED) THEN
C --- JACOBIAN IS BANDED
                  MUJACP=MUJAC+1
                  MD=MIN(MBJAC,N)
                  DO 16 KK=1,MD
                  J=KK
 12               YHH(J)=Y(J)
                  DEL(J)=DSQRT(UROUND*MAX(1.D-5,ABS(Y(J))))
                  Y(J)=Y(J)+DEL(J)
                  J=J+MD
                  IF (J.LE.N) GOTO 12 
                  CALL FCN(N,X,Y,YH)
                  J=KK
                  LBEG=MAX(1,J-MUJAC)
 14               LEND=MIN(N,J+MLJAC)
                  Y(J)=YHH(J)
                  MUJACJ=MUJACP-J
                  DO 15 L=LBEG,LEND
 15               FJAC(L+MUJACJ,J)=(YH(L)-DY(L))/DEL(J) 
                  J=J+MD
                  LBEG=LEND+1
                  IF (J.LE.N) GOTO 14
 16               CONTINUE
              ELSE
C --- JACOBIAN IS FULL
                  DO 18 I=1,N
                  YSAFE=Y(I)
                  DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE)))
                  Y(I)=YSAFE+DELT
                  CALL FCN(N,X,Y,YH)
                  DO 17 J=1,N
  17              FJAC(J,I)=(YH(J)-DY(J))/DELT
  18              Y(I)=YSAFE
                  MLJAC=N
              END IF
          ELSE
C --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
              CALL JAC(N,X,Y,FJAC,LFJAC)
          END IF
          CALJAC=.TRUE.
          IF (NHESS.NE.0) CALL ELMHES (LFJAC,N,1,N,FJAC,IPHES) 
       END IF
       IF (.NOT.AUTNMS.AND..NOT.REJECT) THEN
           IF (IDFX.EQ.0) THEN
C --- COMPUTE NUMERICALLY THE DERIVATIVE WITH RESPECT TO X
               DELT=DSQRT(UROUND*MAX(1.D-5,ABS(X)))
               XDELT=X+DELT
               CALL FCN(N,XDELT,Y,YH)
               DO 19 J=1,N
  19           FX(J)=(YH(J)-DY(J))/DELT
           ELSE
C --- COMPUTE ANALYTICALLY THE DERIVATIVE WITH RESPECT TO X
               CALL DFX(N,X,Y,FX)
           END IF
       END IF
C *** *** *** *** *** *** ***
C --- THE FIRST AND LAST STEP 
C *** *** *** *** *** *** ***
       IF (NSTEP.EQ.0.OR.LAST) THEN 
          NSTEP=NSTEP+1 
          DO 20 J=1,K
          KC=J
       CALL SIMEX(J,N,FCN,X,Y,DY,FX,FJAC,LFJAC,FMAS,LDMAS,E,LE,IP,H,KM,
     &              HMAXN,T,SCAL,NJ,HH,W,A,YHH,DYH,DEL,FXH,WH,ERR,SAFE1,
     &              FAC,FAC1,FAC2,SAFE2,THETA,NHESS,MLJAC,MUJAC,MBJAC,
     &              MLB,MUB,MBB,MBDIAG,MDIFF,MLE,MUE,MDIAG,ERROLD,IPHES,
     &              AUTNMS,IMPLCT,BANDED,REJECT,ATOV)
          IF (ATOV) GOTO 10
          IF (J.GT.1.AND.ERR.LE.1.D0) GO TO 60
  20      CONTINUE
          GO TO 55
       END IF
C --- BASIC INTEGRATION STEP  
  30   CONTINUE
       NSTEP=NSTEP+1
       IF (NSTEP.GE.NMAX) GO TO 120 
       KC=K-1
       DO 40 J=1,KC
       CALL SIMEX(J,N,FCN,X,Y,DY,FX,FJAC,LFJAC,FMAS,LDMAS,E,LE,IP,H,KM,
     &              HMAXN,T,SCAL,NJ,HH,W,A,YHH,DYH,DEL,FXH,WH,ERR,SAFE1,
     &              FAC,FAC1,FAC2,SAFE2,THETA,NHESS,MLJAC,MUJAC,MBJAC,
     &              MLB,MUB,MBB,MBDIAG,MDIFF,MLE,MUE,MDIAG,ERROLD,IPHES,
     &              AUTNMS,IMPLCT,BANDED,REJECT,ATOV)
       IF (ATOV) GOTO 10
  40   CONTINUE
C *** *** *** *** *** *** ***
C --- CONVERGENCE MONITOR
C *** *** *** *** *** *** ***
       IF (K.EQ.2.OR.REJECT) GO TO 50
       IF (ERR.LE.1.D0) GO TO 60
       IF (ERR.GT.(DBLE(NJ(K+1)*NJ(K))/4.D0)**2) GO TO 100  
 50    CALL SIMEX(K,N,FCN,X,Y,DY,FX,FJAC,LFJAC,FMAS,LDMAS,E,LE,IP,H,KM,
     &              HMAXN,T,SCAL,NJ,HH,W,A,YHH,DYH,DEL,FXH,WH,ERR,SAFE1,
     &              FAC,FAC1,FAC2,SAFE2,THETA,NHESS,MLJAC,MUJAC,MBJAC,
     &              MLB,MUB,MBB,MBDIAG,MDIFF,MLE,MUE,MDIAG,ERROLD,IPHES,
     &              AUTNMS,IMPLCT,BANDED,REJECT,ATOV)
       IF (ATOV) GOTO 10
       KC=K 
       IF (ERR.LE.1.D0) GO TO 60
C --- HOPE FOR CONVERGENCE IN LINE K+1
  55   IF (ERR.GT.(DBLE(NJ(K+1))/2.D0)**2) GO TO 100  
       KC=K+1
       CALL SIMEX(KC,N,FCN,X,Y,DY,FX,FJAC,LFJAC,FMAS,LDMAS,E,LE,IP,H,KM,
     &              HMAXN,T,SCAL,NJ,HH,W,A,YHH,DYH,DEL,FXH,WH,ERR,SAFE1,
     &              FAC,FAC1,FAC2,SAFE2,THETA,NHESS,MLJAC,MUJAC,MBJAC,
     &              MLB,MUB,MBB,MBDIAG,MDIFF,MLE,MUE,MDIAG,ERROLD,IPHES,
     &              AUTNMS,IMPLCT,BANDED,REJECT,ATOV)
       IF (ATOV) GOTO 10
       IF (ERR.GT.1.D0) GO TO 100
C *** *** *** *** *** *** ***
C --- STEP IS ACCEPTED  
C *** *** *** *** *** *** ***
  60   XOLD=X
       X=X+H
       DO 70 I=1,N
       T1I=T(1,I)
       IF (ITOL.EQ.0) THEN
         SCAL(I)=ATOL(1)+RTOL(1)*DABS(T1I)
       ELSE
         SCAL(I)=ATOL(I)+RTOL(I)*DABS(T1I)
       END IF
  70   Y(I)=T1I
       NACCPT=NACCPT+1  
       CALJAC=.FALSE.
       IF (IOUT.NE.0) THEN
          IRTRN=1
          NRSOL=NACCPT+1
          XSOL=X
          DO 75 I=1,N
  75      YH(I)=Y(I)
          NSOLU=N
          CALL SOLOUT(NRSOL,XOLD,XSOL,YH,NSOLU,IRTRN)
          IF (IRTRN.LT.0) GOTO 120
       END IF
C --- COMPUTE OPTIMAL ORDER
       IF (KC.EQ.2) THEN
          KOPT=3  
          IF (REJECT) KOPT=2  
          GO TO 80
       END IF
       IF (KC.LE.K) THEN
          KOPT=KC 
          IF (W(KC-1).LT.W(KC)*FAC3) KOPT=KC-1  
          IF (W(KC).LT.W(KC-1)*FAC4) KOPT=MIN0(KC+1,KM-1)
       ELSE 
          KOPT=KC-1
          IF (KC.GT.3.AND.W(KC-2).LT.W(KC-1)*FAC3) KOPT=KC-2
          IF (W(KC).LT.W(KOPT)*FAC4) KOPT=MIN0(KC,KM-1)
       END IF
C --- AFTER A REJECTED STEP
  80   IF (REJECT) THEN 
          K=MIN0(KOPT,KC)
          H=POSNEG*MIN(ABS(H),ABS(HH(K)))
          REJECT=.FALSE.
          GO TO 10
       END IF
C --- COMPUTE STEP SIZE FOR NEXT STEP
       IF (KOPT.LE.KC) THEN
          H=HH(KOPT)
       ELSE 
          IF (KC.LT.K.AND.W(KC).LT.W(KC-1)*FAC4) THEN 
             H=HH(KC)*A(KOPT+1)/A(KC)
          ELSE
             H=HH(KC)*A(KOPT)/A(KC) 
          END IF  
       END IF
       K=KOPT
       H=POSNEG*ABS(H)
       GO TO 10
C *** *** *** *** *** *** ***
C --- STEP IS REJECTED  
C *** *** *** *** *** *** ***
 100   K=MIN(K,KC)
       IF (K.GT.2.AND.W(K-1).LT.W(K)*FAC3) K=K-1
       NREJCT=NREJCT+1  
       H=POSNEG*HH(K)
       LAST=.FALSE.
       REJECT=.TRUE.
       IF (CALJAC) GOTO 30
       GO TO 10
C --- SOLUTION EXIT
 110   CONTINUE
       H=HOPT
       IDID=1
       RETURN
C --- FAIL EXIT
 120   WRITE (6,979) X,H
 979   FORMAT(' EXIT OF SODEX AT X=',D14.7,'   H=',D14.7)
       IDID=-1
       RETURN
       END  
C
C
C *** *** *** *** *** *** ***
C     S U B R O U T I N E    S I M E X
C *** *** *** *** *** *** ***
C
      SUBROUTINE SIMEX(JJ,N,FCN,X,Y,DY,FX,FJAC,LFJAC,FMAS,LDMAS,E,LE,IP,
     &          H,KM,HMAXN,T,SCAL,NJ,HH,W,A,YH,DYH,DEL,FXH,WH,ERR,SAFE1,
     &          FAC,FAC1,FAC2,SAFE2,THETA,NHESS,MLJAC,MUJAC,MBJAC,
     &          MLB,MUB,MBB,MBDIAG,MDIFF,MLE,MUE,MDIAG,ERROLD,IPHES,
     &          AUTNMS,IMPLCT,BANDED,REJECT,ATOV)
C --- THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE
C --- EXTRAPOLATION TABLE AND PROVIDES AN ESTIMATE  
C --- OF THE OPTIMAL STEP SIZE 
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION Y(N),YH(N),DY(N),FX(N),DYH(N),DEL(N),FXH(N)
       DIMENSION WH(N),SCAL(N),HH(KM),W(KM),A(KM),FJAC(LFJAC,N),E(LE,N)
       DIMENSION FMAS(LDMAS,N),T(KM,N),IP(N),NJ(KM),IPHES(N)
       LOGICAL ATOV,REJECT,AUTNMS,IMPLCT,BANDED
       COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
C *** *** *** *** *** *** ***
C  COMPUTE THE MATRIX E AND ITS DECOMPOSITION
C *** *** *** *** *** *** ***
      HJ=H/NJ(JJ)
      HJI=1.D0/HJ
      IF (IMPLCT) THEN
          IF (BANDED) THEN
C --- THE MATRIX E (B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX)
              DO 127 J=1,N
              I1J=MAX0(1,MUJAC+2-J)
              I2J=MIN(MBJAC,MUJAC+1-J+N)
              DO 125 I=I1J,I2J
 125          E(I+MLE,J)=-FJAC(I,J)
              I1B=MAX0(1,MUB+2-J)
              I2B=MIN0(MBB,MUB+1-J+N)
              DO 126 I=I1B,I2B
              IB=I+MDIFF
 126          E(IB,J)=E(IB,J)+HJI*FMAS(I,J)
 127          CONTINUE
              CALL DECB(N,LE,E,MLE,MUE,IP,IER)
              IF (IER.NE.0) GOTO 79
          ELSE
              IF (MLB.NE.N) THEN
C --- THE MATRIX E (B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX)
                  DO 225 J=1,N
                  DO 225 I=1,N
 225              E(I,J)=-FJAC(I,J)
                  DO 226 J=1,N
                  I1=MAX0(1,J-MUB)
                  I2=MIN0(N,J+MLB)
                  DO 226 I=I1,I2
 226              E(I,J)=E(I,J)+HJI*FMAS(I-J+MBDIAG,J)
                  CALL DEC(N,LE,E,IP,IER)
                  IF (IER.NE.0) GOTO 79
              ELSE
C --- THE MATRIX E (B IS A FULL MATRIX, JACOBIAN A FULL MATRIX)
                  IF (MLJAC.EQ.N) THEN
                      DO 324 J=1,N
                      DO 324 I=1,N
 324                  E(I,J)=FMAS(I,J)*HJI-FJAC(I,J)
                      CALL DEC(N,LE,E,IP,IER)
                      IF (IER.NE.0) GOTO 79
                  ELSE
C --- THE MATRIX E (B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX)
                      MADD=MUJAC+1
                      DO 405 J=1,N
                      DO 405 I=1,N
  405                 E(I,J)=FMAS(I,J)*HJI
                      DO 406 J=1,N
                      I1=MAX0(1,J-MUJAC)
                      I2=MIN0(N,J+MLJAC)
                      DO 406 I=I1,I2
  406                 E(I,J)=E(I,J)-FJAC(I-J+MADD,J)
                      CALL DEC(N,LE,E,IP,IER)
                      IF (IER.NE.0) GOTO 79
                  END IF
              END IF
          END IF
      ELSE
          IF (BANDED) THEN
C --- THE MATRIX E (B=IDENTITY, JACOBIAN A BANDED MATRIX)
              DO 427 J=1,N
              I1J=MAX0(1,MUJAC+2-J)
              I2J=MIN(MBJAC,MUJAC+1-J+N)
              DO 425 I=I1J,I2J
 425          E(I+MLE,J)=-FJAC(I,J)
 427          E(MDIAG,J)=E(MDIAG,J)+HJI
              CALL DECB(N,LE,E,MLE,MUE,IP,IER)
              IF (IER.NE.0) GOTO 79
          ELSE
C --- THE MATRIX E (B=IDENTITY, JACOBIAN A FULL MATRIX)
              IF (NHESS.EQ.0) THEN
                  DO 526 J=1,N
                  DO 525 I=1,N
 525              E(I,J)=-FJAC(I,J)
 526              E(J,J)=E(J,J)+HJI
                  CALL DEC(N,LE,E,IP,IER)
                  IF (IER.NE.0) GOTO 79 
              ELSE
                  DO 624 J=1,N-1
                  J1=J+1
 624              E(J1,J)=-FJAC(J1,J)
                  DO 626 J=1,N
                  DO 625 I=1,J
 625              E(I,J)=-FJAC(I,J)
 626              E(J,J)=E(J,J)+HJI
                  CALL DECH(N,LE,E,1,IP,IER)
                  IF (IER.NE.0) GOTO 79 
              END IF
          END IF
      END IF
      NDEC=NDEC+1
C *** *** *** *** *** *** ***
C --- STARTING PROCEDURE
C *** *** *** *** *** *** ***
       DO 20 I=1,N
       YH(I)=Y(I)
  20   DEL(I)=DY(I)
       IF (.NOT.AUTNMS) THEN
           DO 634 I=1,N
 634       DEL(I)=DEL(I)+HJ*FX(I)
       END IF
       IF (BANDED) THEN
           CALL SOLB(N,LE,E,MLE,MUE,DEL,IP)
       ELSE
           IF (NHESS.EQ.0) THEN
               CALL SOL(N,LE,E,DEL,IP)
           ELSE
               DO 140 MMM=N-2,1,-1
               MP=N-MMM
               MP1=MP-1
               I=IPHES(MP)
               IF (I.EQ.MP) GOTO 110
               ZSAFE=DEL(MP)
               DEL(MP)=DEL(I)
               DEL(I)=ZSAFE
 110           CONTINUE
               DO 100 I=MP+1,N 
 100           DEL(I)=DEL(I)-FJAC(I,MP1)*DEL(MP)
 140           CONTINUE
                   CALL SOLH(N,LE,E,1,DEL,IP)
               DO 240 MMM=1,N-2
               MP=N-MMM
               MP1=MP-1
               DO 200 I=MP+1,N 
 200           DEL(I)=DEL(I)+FJAC(I,MP1)*DEL(MP)
               I=IPHES(MP)
               IF (I.EQ.MP) GOTO 240
               ZSAFE=DEL(MP)
               DEL(MP)=DEL(I)
               DEL(I)=ZSAFE
 240           CONTINUE
           END IF
       END IF
       NSOL=NSOL+1
       M=NJ(JJ)
C *** *** *** *** *** *** ***
C --- SEMI-IMPLICIT MID-POINT RULE
C *** *** *** *** *** *** ***
       DO 30 MM=1,M
       DO 23 I=1,N
  23   YH(I)=YH(I)+DEL(I)
       CALL FCN(N,X+HJ*MM,YH,DYH)
       NFCN=NFCN+1
       IF (IMPLCT) THEN
           IF (MLB.EQ.N) THEN
               DO 361 I=1,N
               SUM=0.D0
               DO 360 J=1,N
 360           SUM=SUM+FMAS(I,J)*DEL(J)
 361           WH(I)=SUM
           ELSE
               DO 363 I=1,N
               SUM=0.D0
               J1B=MAX0(1,I-MLB)
               J2B=MIN0(N,I+MUB)
               DO 362 J=J1B,J2B
 362           SUM=SUM+FMAS(I-J+MBDIAG,J)*DEL(J)
 363           WH(I)=SUM
           END IF
           IF (MM.EQ.M) GOTO 30
           DO 24 I=1,N
  24       DYH(I)=2.D0*(DYH(I)-HJI*WH(I))
       ELSE
           IF (MM.EQ.M) GOTO 30
           DO 28 I=1,N
  28       DYH(I)=2.D0*(DYH(I)-HJI*DEL(I))
       END IF
       IF (BANDED) THEN
           CALL SOLB(N,LE,E,MLE,MUE,DYH,IP)
       ELSE
                 IF (NHESS.EQ.0) THEN
                     CALL SOL(N,LE,E,DYH,IP)
                 ELSE
                     DO 340 MMM=N-2,1,-1
                     MP=N-MMM
                     MP1=MP-1
                     I=IPHES(MP)
                     IF (I.EQ.MP) GOTO 310
                     ZSAFE=DYH(MP)
                     DYH(MP)=DYH(I)
                     DYH(I)=ZSAFE
 310                 CONTINUE
                     DO 300 I=MP+1,N 
 300                 DYH(I)=DYH(I)-FJAC(I,MP1)*DYH(MP)
 340                 CONTINUE
                         CALL SOLH(N,LE,E,1,DYH,IP)
                     DO 440 MMM=1,N-2
                     MP=N-MMM
                     MP1=MP-1
                     DO 400 I=MP+1,N 
 400                 DYH(I)=DYH(I)+FJAC(I,MP1)*DYH(MP)
                     I=IPHES(MP)
                     IF (I.EQ.MP) GOTO 440
                     ZSAFE=DYH(MP)
                     DYH(MP)=DYH(I)
                     DYH(I)=ZSAFE
 440                 CONTINUE
                 END IF
       END IF
       NSOL=NSOL+1
       IF (MM.EQ.1.AND.JJ.LE.2) THEN
C --- STABILITY CHECK
          DEL1=0.D0
          DO 21 I=1,N
  21      DEL1=DEL1+(DEL(I)/SCAL(I))**2
          DEL1=DSQRT(DEL1)
          DEL2=0.D0
          DO 26 I=1,N
  26      DEL2=DEL2+(DYH(I)/SCAL(I))**2
          DEL2=.5D0*DSQRT(DEL2)
          THETA=DEL2/MAX(1.D0,DEL1)
          IF (THETA.GT.1.D0) GOTO 79
       END IF
       DO 25 I=1,N
  25   DEL(I)=DEL(I)+DYH(I)
  30   CONTINUE
C --- FINAL STEP (DUE TO BADER)
       IF (IMPLCT) THEN
           DO 34 I=1,N
  34       DYH(I)=DYH(I)-HJI*WH(I)
       ELSE
           DO 36 I=1,N
  36       DYH(I)=DYH(I)-HJI*DEL(I)
       END IF
       IF (BANDED) THEN
           CALL SOLB(N,LE,E,MLE,MUE,DYH,IP)
       ELSE
               IF (NHESS.EQ.0) THEN
                   CALL SOL(N,LE,E,DYH,IP)
                ELSE
                   DO 540 MMM=N-2,1,-1
                   MP=N-MMM
                   MP1=MP-1
                   I=IPHES(MP)
                   IF (I.EQ.MP) GOTO 510
                   ZSAFE=DYH(MP)
                   DYH(MP)=DYH(I)
                   DYH(I)=ZSAFE
 510               CONTINUE
                   DO 500 I=MP+1,N 
 500               DYH(I)=DYH(I)-FJAC(I,MP1)*DYH(MP)
 540               CONTINUE
                       CALL SOLH(N,LE,E,1,DYH,IP)
                   DO 640 MMM=1,N-2
                   MP=N-MMM
                   MP1=MP-1
                   DO 600 I=MP+1,N 
 600               DYH(I)=DYH(I)+FJAC(I,MP1)*DYH(MP)
                   I=IPHES(MP)
                   IF (I.EQ.MP) GOTO 640
                   ZSAFE=DYH(MP)
                   DYH(MP)=DYH(I)
                   DYH(I)=ZSAFE
 640               CONTINUE
               END IF
       END IF
       NSOL=NSOL+1
       DO 35 I=1,N
  35   T(JJ,I)=YH(I)+DYH(I)
C *** *** *** *** *** *** ***
C --- POLYNOMIAL EXTRAPOLATION
C *** *** *** *** *** *** ***
       IF (JJ.EQ.1) RETURN
       DO 60 L=JJ,2,-1
       FAC=(DBLE(NJ(JJ))/DBLE(NJ(L-1)))**2-1.D0
       DO 60 I=1,N
       T(L-1,I)=T(L,I)+(T(L,I)-T(L-1,I))/FAC
  60   CONTINUE
       ERR=0.D0
       DO 65 I=1,N
  65   ERR=ERR+MIN(ABS((T(1,I)-T(2,I)))/SCAL(I),1.D10)**2  
       IF (ERR.GE.1.D20) GOTO 79
       ERR=DSQRT(ERR/DBLE(N))
       IF (JJ.GT.2.AND.ERR.GE.ERROLD) GOTO 79
       ERROLD=DMAX1(4*ERR,1.D0)
C --- COMPUTE OPTIMAL STEP SIZES
       EXPO=1.D0/(2*JJ-2)
       FACMIN=FAC1**EXPO
       FAC=MIN(FAC2/FACMIN,MAX(FACMIN,(ERR/SAFE1)**EXPO/SAFE2))
       FAC=1.D0/FAC
       HH(JJ)=MIN(ABS(H)*FAC,HMAXN)
       W(JJ)=A(JJ)/HH(JJ)  
       RETURN
  79   ATOV=.TRUE.
       H=H*0.5D0
       REJECT=.TRUE.
       RETURN
       END  
