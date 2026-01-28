C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RADAR5 at Example 2 of "Applying codes for ..."
C * * * * * * * * * * * * * * *
        include 'radar5.f'
        include '/Users/hairer/fortran/programme/dgln_stiff/decsol.f'
        include 'dc_sumexpdel.f'
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*4 start,finish
C ---  PARAMETERS FOR RADAR5       
        INTEGER, PARAMETER :: NDIM=200,NRDENS=2,NGRID=10,NLAGS=2
        INTEGER, PARAMETER :: NJACL=NDIM,MXST=4000,LWORK=30,LIWORK=30
        DIMENSION Y(NDIM),RPAR(4),IPAR(4),WORK(LWORK),IWORK(LIWORK)
        DIMENSION RTOL(NDIM),ATOL(NDIM),IPASTS(NRDENS),GRID(NGRID+1)
        DIMENSION IPAST(NRDENS+2*NDIM)
        EXTERNAL  FCN,PHI,ARGLAG,JFCN,JFCNEX,JACLAG,MAS,SOLOUT
C --- PARAMETERS OF THE PROBLEM
        alpha=0.5D0
        beta = 1.0D0
        RPAR(1)=alpha
        RPAR(2)=beta
C --- DESIRED ACCURACY
        TOL=1.0D-6
        EPSIL=TOL
C --- END OF INTEGRATION
        XEND=10.D0
        CALL HMNPAR (XEND,alpha,beta,epsil,T,HEXP,MEXP,NEXP)
        write (6,*) '  '
        WRITE(*,*) 'HEXP, MEXP, NEXP ', HEXP, MEXP, NEXP
        RPAR(4)=HEXP
        IPAR(1)=MEXP
        IPAR(2)=NEXP
C --- FULL LINEAR ALGEBRA (ILIN = 1) AND ADAPTED TO SUMEXP (ILIN = 2)
        DO ILIN = 1,2
        WRITE (6,*) '  '
        WRITE (6,*) '   ILIN = ',ILIN
           MI=1
           ND=2
           N=ND+MI*(NEXP-MEXP)
           IF (ILIN.EQ.2) N=N+2
C --- INITIAL VALUES 
           X=0.0D0
           Y=0.0D0
           Y(1)=PHI(1,X,RPAR,IPAR)
C --- DIFFERENTIAL EQUATION IS IN IMPLICIT FORM ?
           IMAS=1
	       MLMAS=0
           MUMAS=0
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
           IOUT=0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
           ITOL=0
           RTOL=TOL
           ATOL=TOL
C --- INITIAL STEP SIZE
           H=TOL
C --- DEFAULT VALUES FOR PARAMETERS
           IWORK = 0
           WORK = 0.0D0  
C --- WORKSPACE FOR PAST 
           IWORK(12)=MXST
C --- THE FIRST TWO COMPONENTS USE RETARDED ARGUMENT
           IWORK(15)=NRDENS
           IPAST(1)=1
           IPAST(2)=2
C --- SET THE PRESCRIBED GRID-POINTS
           TAU=BETA    
           PI=4.D0*ATAN(1.D0)
           TAU2=PI/4.D0
           GRID=0.D0
           GRID(1)=TAU2
           GRID(2)=TAU
           GRID(3)=2*TAU2
           GRID(4)=TAU2+TAU
           GRID(5)=2*TAU
           GRID(6)=3*TAU2
           GRID(7)=2*TAU2+TAU
           GRID(8)=2*TAU+TAU2
           GRID(9)=3*TAU
           GRID(10)=4*TAU2
           IWORK(13)=NGRID
           LGRID=NGRID+1
C --- CALL OF TIME INTEGRATOR
           IF (ILIN.EQ.1) THEN
             IJAC=1
             MLJAC=N
             MUJAC=N
             CALL cpu_time(start)
             CALL RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JFCN,IJAC,MLJAC,MUJAC,
     &                  JACLAG,NLAGS,NJACL,
     &                  MAS,IMAS,MLMAS,MUMAS,SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID,
     &                  GRID,LGRID,IPAST,NRDENS)   
             CALL cpu_time(finish)
           END IF
           IF (ILIN.EQ.2) THEN
             IJAC=1
             MLJAC=MAX(3-ND,0)
             MUJAC=ND-1
             CALL cpu_time(start)
             CALL RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JFCNEX,IJAC,MLJAC,MUJAC,
     &                  JACLAG,NLAGS,NJACL,
     &                  MAS,IMAS,MLMAS,MUMAS,SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID,
     &                  GRID,LGRID,IPAST,NRDENS)   
             CALL cpu_time(finish)
           END IF
C --- PRINT FINAL SOLUTION
        YEX = 0.570525788119D0
        WRITE (*,99) X,Y(1),ABS(Y(1)-YEX)/YEX
 99     FORMAT(1X,'X =',F5.2,'   Y = ',E18.12,'   error = ',E18.12)
C --- PRINT STATISTICS
        WRITE(6,*)' *** TOL=',TOL,'  TIME=',finish-start,' ***'
        WRITE (6,91) (IWORK(J),J=14,20)
 91     FORMAT(' fcn=',I5,' jac=',I5,' step=',I5,
     &        ' accpt=',I5,' rejct=',I5,' dec=',I5,
     &        ' sol=',I5)
        END DO
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,HSOL,Y,CONT,LRC,N,
     &                     RPAR,IPAR,IRTRN)
C ----- PRINTS THE SOLUTION AT EQUIDISTANT OUTPUT-POINTS
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N), CONT(LRC), RPAR(4), IPAR(4)
        COMMON /INTERN/XOUT
        XSTEP=0.5D0
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1)
           XOUT=XSTEP
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (6,99) XOUT,CONTR5(1,N,XOUT,CONT,X,HSOL)
              XOUT=XOUT+XSTEP
              GOTO 10
           END IF
        END IF
   99   FORMAT(1X,'X =',F12.6,'    Y =',F18.12)
        RETURN
        END
C
        FUNCTION ARGLAG(IL,X,N,Y,RPAR,IPAR,PHI,PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),RPAR(4),IPAR(4),PAST(*),IPAST(*)
        EXTERNAL PHI
        PI=4.D0*ATAN(1.D0)
        IF (IL.EQ.1) ARGLAG=X-PI/4.D0  
        IF (IL.EQ.2) ARGLAG=X-RPAR(2) 
        RETURN
        END
C
        SUBROUTINE FCN(N,X,Y,F,ARGLAG,PHI,RPAR,IPAR,
     &                  PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,K,O-Z)
        DIMENSION Y(N),F(N),RPAR(4),IPAR(4),PAST(N),IPAST(N)
        EXTERNAL PHI,ARGLAG
        F = 0.0D0
        MEXP=IPAR(1)
        NEXP=IPAR(2)
        HEXP=RPAR(4)
        alpha=RPAR(1)
        beta=RPAR(2)
        LEXP=NEXP-MEXP
        CALL LAGR5(1,X,N,Y,ARGLAG,PAST,THETA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
        Y1L1=YLAGR5(1,THETA1,IPOS1,PHI,RPAR,IPAR,PAST,IPAST,NRDS)
        CALL LAGR5(2,X,N,Y,ARGLAG,PAST,THETA2,IPOS2,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
        Y2L2=YLAGR5(2,THETA2,IPOS2,PHI,RPAR,IPAR,PAST,IPAST,NRDS)
        SUMZ=0.D0
        DO I=3,LEXP+2
           IND=MEXP+(I-3)
           GAMMAI =  EXP(IND*HEXP)
           CI =      EXP((alpha+1.D0)*IND*HEXP-beta*GAMMAI)
           F(I)   = - GAMMAI*Y(I) + Y(1)
           SUMZ =  SUMZ + CI*Y(I)
        END DO
        FAC=alpha*(beta**alpha)*HEXP/GAMMA(alpha+1)
        F(2) = Y(2) - FAC*SUMZ
        F(1) = -5.D0*Y2L2 - (Y1L1-2.D0)/(Y(1)+1.D0)
        RETURN
        END
C
        SUBROUTINE JFCN(N,X,Y,DFY,LDFY,ARGLAG,PHI,RPAR,IPAR,
     &                  PAST,IPAST,NRDS)
C ----- STANDARD JACOBIAN OF THE EQUATION
        IMPLICIT REAL*8 (A-H,K,O-Z)
        DIMENSION Y(N),DFY(LDFY,N),PAST(N),IPAST(N),RPAR(4),IPAR(4)
        EXTERNAL PHI,ARGLAG
        MEXP=IPAR(1)
        NEXP=IPAR(2)
        HEXP=RPAR(4)
        alpha=RPAR(1)
        beta=RPAR(2)
        LEXP=NEXP-MEXP
        DFY = 0.0D0    
        CALL LAGR5(1,X,N,Y,ARGLAG,PAST,THETA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
        Y1L1=YLAGR5(1,THETA1,IPOS1,PHI,RPAR,IPAR,PAST,IPAST,NRDS)
        DFY(1,1)=(Y1L1-2.D0)/(Y(1)+1.D0)**2
        DFY(2,2)=1.0D0
        FAC=alpha*(beta**alpha)*HEXP/GAMMA(alpha+1)
        DO I=3,LEXP+2
          IND=MEXP+(I-3)
          GAMMAI =  EXP(IND*HEXP)
          CI =      EXP((alpha+1.D0)*IND*HEXP-beta*GAMMAI)
          DFY(2,I) = - CI*FAC
          DFY(I,I) = - GAMMAI
          DFY(I,1) = 1.0D0
        END DO
        RETURN
        END
C
        SUBROUTINE JFCNEX(N,X,Y,DFY,LDFY,ARGLAG,PHI,RPAR,IPAR,
     &                  PAST,IPAST,NRDS)
C --- JACOBIAN for the option "sumexp" in the linear algebra
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N),PAST(N),IPAST(N),RPAR(4),IPAR(4)
        EXTERNAL PHI,ARGLAG
        MEXP=IPAR(1)
        NEXP=IPAR(2)
        HEXP=RPAR(4)
        alpha=RPAR(1)
        beta=RPAR(2)
        LEXP=NEXP-MEXP
        DFY = 0.0D0    
        CALL LAGR5(1,X,N,Y,ARGLAG,PAST,THETA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
        Y1L1=YLAGR5(1,THETA1,IPOS1,PHI,RPAR,IPAR,PAST,IPAST,NRDS)
        DFY(1,1)=(Y1L1-2.D0)/(Y(1)+1.D0)**2
        DFY(2,2)=1.0D0
        DFY(1,3)=0.0D0
        DFY(2,3)=-1.0D0
        DFY(1,4)=1.0D0
        DFY(2,4)=0.0D0
        FAC=alpha*(beta**alpha)*HEXP/GAMMA(alpha+1)
        II=5
        DO I=1,LEXP
          IND=MEXP+(I-1)
          GAMMAI =  EXP(IND*HEXP)
          CI =      EXP((alpha+1.D0)*IND*HEXP-beta*GAMMAI)
          DFY(1,II) = CI*FAC
          DFY(2,II) = - GAMMAI
          DFY(3,II) = 0.0D0
          II=II+1
        END DO
        RETURN
        END
C    
        SUBROUTINE JACLAG(N,X,Y,DFYL,ARGLAG,PHI,IVE,IVC,IVL,
     &                    RPAR,IPAR,PAST,IPAST,NRDS)
C ----- JACOBIAN OF DELAY TERMS IN THE EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFYL(2),PAST(1),IPAST(1),RPAR(4),IPAR(4)
        DIMENSION IVE(2),IVC(2),IVL(2)
        EXTERNAL PHI,ARGLAG
        IVE(1)=1
        IVC(1)=1
        IVL(1)=1
        DFYL(1) = - 1.D0/(Y(1)+1.D0)
        IVE(2)=1
        IVC(2)=2
        IVL(2)=2
        DFYL(2) = -5.D0   
        RETURN
        END
C
        FUNCTION PHI(I,X,RPAR,IPAR)
C --- INITIAL FUNCTIONS FOR THE TEST PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RPAR(4),IPAR(4)
        IF (I.EQ.1) PHI = X
        IF (I.EQ.2) PHI = 0.0D0
        RETURN
        END
C
        SUBROUTINE MAS(N,Q,LQ,RPAR,IPAR)
C --- MATRIX "M" FOR THE TEST PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RPAR(4),IPAR(4),Q(LQ,N)
        DO I=1,N
           Q(1,I)=1.0D0
        END DO
        Q(1,2) = 0.0D0
        RETURN
        END
C        
        SUBROUTINE HMNPAR (XEND,alpha,beta,epsil,T,HEXP,MEXP,NEXP)
        IMPLICIT REAL*8 (A-H,O-Z)
C --- COMPUTES THE VALUES OF HEXP,MEXP,NEXP FOR THE PARETO DISTRIBUTION
        PI=4.D0*ATAN(1.D0)
        T=beta*EPSIL**(-1.D0/alpha)
        T=DMIN1(T,XEND)    
        a = (1.D0-(alpha+1.D0)/((alpha+2.D0)*LOG(1.D0/EPSIL)))*PI/2.D0
        HEXP = 2.D0*PI*a/LOG(1.D0+2.D0*COS(a)**(-(alpha+1))/EPSIL)
        xl=  GAMMA(alpha+2.D0)*epsil
        xs= -LOG(GAMMA(alpha+1.D0)*epsil)
        NEXP = INT(LOG(xS/beta)/HEXP)+1
        MEXP = INT(LOG(xL/T)/HEXP)-1
        RETURN
        END