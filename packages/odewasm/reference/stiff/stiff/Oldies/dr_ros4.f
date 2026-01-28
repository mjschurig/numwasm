C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR ROSENBROCK CODES ROS4 AND ROSA AT VAN DER POL
C * * * * * * * * * * * * * * * * * * * * * * * * *
clink dr_ros4 ros4 decsol
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR ROS4 (FULL JACOBIAN)
        PARAMETER (ND=2,LWORK=2*ND*ND+8*ND+5,LIWORK=ND+2)
C --- DECLARATIONS
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
        COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
        EXTERNAL FVPOL,JVPOL,SOLOUT
C --- DIMENSION OF THE SYSTEM
        N=2
C --- PROBLEM IS AUTONOMOUS
        IFCN=0
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=2.0D0
        Y(2)=-0.66D0
C --- ENDPOINT OF INTEGRATION
        XEND=2.0D0
C --- REQUIRED TOLERANCE
        RTOL=1.0D-5
        ATOL=1.0D-6*RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6 
C --- SET DEFAULT VALUES (NPAR=5 FOR ROS4, =4 FOR ROSA)
        NPAR=5
C        NPAR=4
        DO 10 I=1,NPAR
  10    WORK(I)=0.D0
        IWORK(1)=0
        IWORK(2)=0
C --- CALL OF THE SUBROUTINE ROS4
        CALL ROS4(N,FVPOL,IFCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JVPOL,IJAC,MLJAC,MUJAC,FVPOL,IDFX,
     &                  FVPOL,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,IDID)
C --- PRINT FINAL SOLUTION
        WRITE (6,99) X,Y(1),Y(2)
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10)
C --- PRINT STATISTICS
        WRITE (6,90) RTOL
 90     FORMAT('       rtol=',D8.2)
        WRITE (6,91) NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,
     &        ' accpt=',I4,' rejct=',I3,' dec=',I4,
     &        ' sol=',I5)
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,IRTRN)
C --- PRINTS SOLUTION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N)
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),Y(2),NR-1
           XOUT=0.1D0
        ELSE
           IF (X.GE.XOUT) THEN
              WRITE (6,99) X,Y(1),Y(2),NR-1
              XOUT=DMAX1(XOUT+0.1D0,X)
           END IF
        END IF
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10,'    NSTEP =',I4)
        RETURN
        END
C
C
        SUBROUTINE FVPOL(N,X,Y,F)
C --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        EPS=1.0D-6
        F(1)=Y(2)
        F(2)=((1-Y(1)**2)*Y(2)-Y(1))/EPS
        RETURN
        END 
C
C
        SUBROUTINE JVPOL(N,X,Y,DFY,LDFY)
C --- JACOBIAN OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        EPS=1.0D-6
        DFY(1,1)=0.0D0
        DFY(1,2)=1.0D0
        DFY(2,1)=(-2.0D0*Y(1)*Y(2)-1.0D0)/EPS
        DFY(2,2)=(1.0D0-Y(1)**2)/EPS
        RETURN
        END

