C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR EXTRAPOLATION CODES SODEX AT VAN DER POL
C * * * * * * * * * * * * * * * * * * * * * * * * *
clink dr_seulex sodex decsol
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR SODEX (FULL JACOBIAN)
        PARAMETER (ND=2,KM=6,LWORK=2*ND*ND+(KM+9)*ND+3*KM+13)
        PARAMETER (LIWORK=2*ND+KM+4)
C --- DECLARATIONS
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
        COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
        EXTERNAL FVPOL,JVPOL,SOLOUT 
C --- DIMENSION OF THE SYSTEM
        N=2
C --- PROBLEM IS AUTONOMOUS
        IFCN=1
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=0
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
        ATOL=RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6 
C --- SET DEFAULT VALUES
        DO 10 I=1,13
  10    WORK(I)=0.D0
        DO 12 I=1,4
  12    IWORK(I)=0
C --- CALL OF THE SUBROUTINE SODEX
        CALL SODEX(N,FVPOL,IFCN,X,Y,XEND,H,
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
              XOUT=MAX(XOUT+0.1D0,X)
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

