C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RADAR5 
C * * * * * * * * * * * * * * *
        include '../radar5.f'
        include '../decsol.f'
        include '../dc_decdel.f'

        IMPLICIT REAL*8 (A-H,O-Z)

        REAL*4 start,finish
        INTEGER, PARAMETER :: DP=kind(1D0)
C --->  PARAMETERS FOR RADAR5 (FULL JACOBIAN) <---
        INTEGER, PARAMETER :: ND=1
        INTEGER, PARAMETER :: NRDENS=1
        INTEGER, PARAMETER :: NGRID=30
        INTEGER, PARAMETER :: NLAGS=1
        INTEGER, PARAMETER :: NJACL=1
        INTEGER, PARAMETER :: MXST=10000
        INTEGER, PARAMETER :: LWORK=30
        INTEGER, PARAMETER :: LIWORK=30
        REAL(kind=DP), dimension(ND) :: Y
        REAL(kind=DP), dimension(NGRID+1) :: GRID
        REAL(kind=DP), dimension(LWORK) :: WORK
        INTEGER, dimension(LIWORK) :: IWORK
        INTEGER, dimension(NRDENS+2*ND) :: IPAST
        DIMENSION IPAR(1),RPAR(5)
        DIMENSION ATOL(1),RTOL(1)
        REAL(kind=DP), PARAMETER :: A=5.D0
        REAL(kind=DP), PARAMETER :: B=-4.0D0
        REAL(kind=DP), PARAMETER :: TAU=1.D0
        REAL(kind=DP), PARAMETER :: EPS=1.D-6
        REAL(kind=DP), PARAMETER :: PI=3.1415926535897932385D0
        EXTERNAL  FCN,PHI,ARGLAG,JFCN,JACLAG,QFUN,SOLOUT

C ------ FILE TO OPEN ----------
        OPEN(9,FILE='sol.out')
        OPEN(10,FILE='cont.out')
        REWIND 9
        REWIND 10
  
C --- PARAMETERS IN THE DIFFERENTIAL EQUATION
        RPAR(1)=A
        RPAR(2)=B
        RPAR(3)=TAU
        RPAR(4)=EPS
        RPAR(5)=PI
C -----------------------------------------------------------------------
C
C --- DIMENSION OF THE SYSTEM
        N=ND
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
      IMAS=1
      MLMAS=N
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES 
        X=0.0D0
        Y(1)=0.0D0
C       Consistent with initial function
C --- ENDPOINT OF INTEGRATION
        XEND=100.5D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.D-6 
        ATOL=RTOL
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- DEFAULT VALUES FOR PARAMETERS
        IWORK=0
        WORK=0.D0  
        IWORK(2)=1000000
C --- PARAMETERS FOR ERROR ESTIMATION
        WORK(10)=0.D0
        IWORK(11)=1
C --- WORKSPACE FOR PAST
        IWORK(12)=MXST
C --- THE COMPONENT USES RETARDED ARGUMENT
        IWORK(15)=NRDENS
        IPAST(1)=1
C --- SET THE PRESCRIBED GRID-POINTS
        DO I=1,NGRID
          GRID(I)=I*TAU  
        END DO
        LGRID = NGRID+1
C --- WORKSPACE FOR GRID
        IWORK(13)=NGRID
C _____________________________________________________________________
C --- CALL OF THE SUBROUTINE RADAR5   
        CALL cpu_time(start)
        CALL RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JFCN,IJAC,MLJAC,MUJAC,
     &                  JACLAG,NLAGS,NJACL,
     &                  QFUN,IMAS,MLMAS,MUMAS,SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID,
     &                  GRID,LGRID,IPAST,NRDENS)   
        CALL cpu_time(finish)
C --- PRINT FINAL SOLUTION SOLUTION
        WRITE (6,90) X,Y(1)
C --- PRINT STATISTICS
        WRITE(6,*)' *** TOL=',RTOL,' ***','   Time = ',finish-start
        WRITE (6,91) (IWORK(J),J=14,20),IWORK(13)
 90     FORMAT(1X,'X =',F8.2,'    Y =',4E18.10)
 91     FORMAT(' fcn=',I7,' jac=',I6,' step=',I6,
     &        ' accpt=',I6,' rejct=',I6,' dec=',I6,
     &        ' sol=',I7,' fullits =',I7)
        WRITE(6,*) 'SOLUTION IS TABULATED IN FILES: sol.out & cont.out'
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,HSOL,Y,CONT,LRC,N,
     &                     RPAR,IPAR,IRTRN)

C ----- PRINTS THE DISCRETE OUTPUT AND THE CONTINUOUS OUTPUT
C       AT EQUIDISTANT OUTPUT-POINTS
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), PARAMETER :: XSTEP=0.01D0
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(LRC) :: CONT
        DIMENSION IPAR(1),RPAR(5)
C       XOUT IS USED FOR THE DENSE OUTPUT
        COMMON /INTERN/XOUT

        WRITE (9,99) X,Y(1)
C
        IF (NR.EQ.1) THEN
           WRITE (10,99) X,Y(1)
           XOUT=XSTEP
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (10,99) XOUT,CONTR5(1,N,XOUT,CONT,X,HSOL)
              XOUT=XOUT+XSTEP
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F12.8,'    Y =',E18.10)
        RETURN
        END
C
        FUNCTION ARGLAG(IL,X,N,Y,RPAR,IPAR,PHI,PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(NRDS+2*N) :: IPAST
        DIMENSION IPAR(1),RPAR(5)

        TAU=RPAR(3)
        ARGLAG=X-TAU     
        RETURN
        END
C
        SUBROUTINE FCN(N,X,Y,F,ARGLAG,PHI,RPAR,IPAR,
     &                 PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,K,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(N) :: F
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(NRDS+2*N) :: IPAST
        DIMENSION IPAR(1),RPAR(5)
        EXTERNAL PHI,ARGLAG

        A=RPAR(1)
        B=RPAR(2)
        TAU=RPAR(3)
        EPS=RPAR(4)
        PI=RPAR(5)

C       ONE DELAY
        CALL LAGR5(1,X,N,Y,ARGLAG,PAST,ALPHA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
        Y1LT=YLAGR5(1,ALPHA1,IPOS1,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)

        F(1)=EPS*(1.D0/(1.D0 + X*X)) - Y(1) + (B/A)*Y1LT + 
     &            ((A-B)/A)*PI/2.D0
        RETURN
        END
C
        SUBROUTINE JFCN(N,X,Y,DFY,LDFY,ARGLAG,PHI,RPAR,IPAR,
     &                  PAST,IPAST,NRDS)
C ----- STANDARD JACOBIAN OF THE EQUATION
        IMPLICIT REAL*8 (A-H,K,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(LDFY,N) :: DFY
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(NRDS+2*N) :: IPAST
        DIMENSION IPAR(1),RPAR(5)
        EXTERNAL PHI,ARGLAG

        DFY(1,1)=-1
        RETURN
        END
C
        SUBROUTINE JACLAG(N,X,Y,DFYL,ARGLAG,PHI,IVE,IVC,IVL,
     &                    RPAR,IPAR,PAST,IPAST,NRDS)
C ----- JACOBIAN OF DELAY TERMS IN THE EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(1) :: DFYL
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(NRDS+2*N) :: IPAST
        DIMENSION IPAR(1),RPAR(5)
        INTEGER, dimension(1) :: IVE,IVC,IVL
        EXTERNAL PHI,ARGLAG

        A=RPAR(1)
        B=RPAR(2)

        IVL(1)=1
        IVE(1)=1
        IVC(1)=1
        DFYL(1)=B/A    

        RETURN
        END

C
        FUNCTION PHI(I,X,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION IPAR(1),RPAR(5)
        PHI=ATAN(X)   
        RETURN
        END

        SUBROUTINE QFUN(N,Q,LQ,RPAR,IPAR)
C --- MATRIX "M" FOR THE TEST PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        DIMENSION IPAR(1),RPAR(5)
        REAL(kind=DP), dimension(LQ,N) :: Q
        Q(1,1)=RPAR(4)
        RETURN
        END
