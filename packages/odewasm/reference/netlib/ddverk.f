      SUBROUTINE APFORM(N,NW,MU,OMEGA,LOLD,S,T,Y,H,TOLD,KOLD,TOL,W,C,
     +     FCN,SIGMA,YINIT,DYINIT)

C**********************************************************************
C                                                                     *
C     Purpose - This routine integrates one step by RK-methods        *
c               for delay DEs.                                        *
C                                                                     *
C**********************************************************************
C                                                                     *
C     N      Number of equations                                      *
C     NW     First dimension of workspace W & queue KOLD              *
C     MU     Number of solution delay terms                           *
C     OMEGA  Number of derivative delay terms                         *
C     LOLD   The length of queue                                      *
C     S      The number of stages for RK-formula                      *
C     T      Current independent value                                *
C     Y      Current solution values, i.e., at T                      *
C     H      Trial stepsize                                           *
C     TOLD   Queue of T values                                        *
C     KOLD   Queue oF K values and solution values                    *
C     TOL    Tolerance                                                *
C     W      Workspace                                                *
C     C      Communication vector                                     *
C            ...USE C(13) : Minimum stepsize                          *
C                   C(17) : T trial                                   *
C                   C(18) : H trial                                   *
C                   C(24) : # of function evaluations                 *
C                   C(25) : Initial point                             *
C                   C(26) : The first index of queue                  *
C                   C(27) : The last index of queue                   *
C                   C(28) : Flag for computation of K1 in the next    *
C                           step(if = 0.0, need to compute)           *
C                   C(29) : Flag indicating which of the              *
c                           extrapolation or the special interpolant  *
c                           must be used.                             *
C                           If C(29) = 0.0 then use extrapolation     *
C     FCN    Name of derivative routine                               *
C     SIGMA  Name of subroutine for computing delay argument          *
C     YINIT  Name of subroutine for computing initial values at       *
C            initial intervals.                                       *
C     DYINIT Name of subroutine for computing derivative of initial   *
C            values at initial intervals.                             *
C                                                                     *
C**********************************************************************

C**********************************************************************
C Note:                                                               *
C       When using extrapolation, we use global intrp routine.        *
C       And use told(c(27)-1). But c(27) may be one.                  *
C       In this case, we have to use told(lold - 1)                   *
C**********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,MU,OMEGA,LOLD,S
      DOUBLE PRECISION T,H,TOL

C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N),TOLD(LOLD),KOLD(NW,7,LOLD),W(NW,*),C(*)

C     .. Subroutine Arguments ..
      EXTERNAL FCN,SIGMA,YINIT,DYINIT

C     .. Arrays in Common ..
      DOUBLE PRECISION AA(11,11),BB(11),CC(11)


C     .. Local Scalars ..
      INTEGER I,J,K,L,INDEX,IW,MAXNIT,ND,NIT,TMP
      DOUBLE PRECISION ERRKI,SIGN,TEMP,TEMP2

C     .. External Subroutines ..
      EXTERNAL EVALAGS,EVALAGD,INTRPG,DERIVG

C     .. Intrisic Functions ..
      INTRINSIC DMAX1,DFLOAT,IDINT

C     .. Common blocks ..
      COMMON/METH/AA,BB,CC

C     .. Data statements ..
      DATA MAXNIT/5/

C     .. Set up the nubmer of delay terms
      ND = MU + OMEGA

C     ... LOOP for Advanced Delay ...
 1000 CONTINUE

C     Set flag to be 1
      C(28) = 1.D0
C
C***** Beginning of the 1st phase *****

C     Set index to 0
C ... INDEX: The first vanishing delay occures in the INDEX-th
C     stage if not 0.
      INDEX = 0

C     The 2nd stage
C ... Compute Y values
      DO 110 I = 1, N
         W(I,16) = Y(I) + H*AA(2,1)*W(I,1)
 110  CONTINUE

C ... Compute argument for delay values
      TEMP = T + CC(2)*H
      CALL SIGMA(N,ND,TEMP,W(1,16),W(1,15))

C ... Compute solution delay values
      DO 120 I = 1, MU

         IW = 20 + I

         IF (W(I,15) .GT. TEMP) THEN
            PRINT *, '**** WARNING:APFORM *****'
            PRINT *, 'Advanced Delay at Stage  2'
            PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
            PRINT *

C ... Cut the stepseize to half and integrate again
            IF (W(I,15) - TEMP .GE. 1.D0) THEN
               PRINT *, 'ADVANCED TOO MUCH !!!'
               PRINT *
               C(18) = C(18) / 2.D0
               IF (C(18) .LT. C(13)) THEN
                  PRINT *, 'IND = -3,  At T = ', T
                  PRINT *
                  STOP
               END IF
               C(17) = T + C(18)
               H = C(18)
               GOTO 1000
            END IF
         END IF
C
         SIGN = W(I,15) - T
C
C ... In case there is no vanishing delay
         IF (SIGN .LE. 0.D0) THEN
            CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +           TOLD,KOLD,W,YINIT)

C ... In case there is vanishing delay
         ELSE
            INDEX = 2

C ... Compute delay value using extrapolation
            IF (C(29) .EQ. 0.D0) THEN
	       TMP = IDINT(C(27))
               IF (TMP .EQ. 1) TMP = LOLD
               CALL INTRPG(N,NW,TMP,TOLD(TMP-1),W(I,15),W(1,IW),
     +              T-TOLD(TMP-1),KOLD)

C ... Compute delay value using the special interpolant
            ELSE 
               DO 130 J = 1, N
                  W(J,IW) = Y(J) +  SIGN*W(J,1)
 130           CONTINUE
            END IF
         END IF

 120  CONTINUE

C ... Compute derivative delay values      

      DO 125 I = MU + 1, ND

         IW = 20 + I

         IF (W(I,15) .GT. TEMP) THEN
            PRINT *, '**** WARNING:APFORM *****'
            PRINT *, 'Advanced Delay at Stage  2'
            PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
            PRINT *

C ... Cut the stepseize to half and integrate again
            IF (W(I,15) - TEMP .GE. 1.D0) THEN
               PRINT *, 'ADVANCED TOO MUCH !!!'
               PRINT *
               C(18) = C(18) / 2.D0
               IF (C(18) .LT. C(13)) THEN
                  PRINT *, 'IND = -3,  At T = ', T
                  PRINT *
                  STOP
               END IF
               C(17) = T + C(18)
               H = C(18)
               GOTO 1000
            END IF
         END IF
C
         SIGN = W(I,15) - T
C
C ... In case there is no vanishing delay
         IF (SIGN .LE. 0.D0) THEN
            CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +           TOLD,KOLD,W,DYINIT)

C ... In case there is vanishing delay
         ELSE
            INDEX = 2

C ... Compute delay value using extrapolation
            IF (C(29) .EQ. 0.D0) THEN
	       TMP = IDINT(C(27))
               IF (TMP .EQ. 1) TMP = LOLD
               CALL DERIVG(N,NW,TMP,TOLD(TMP-1),W(I,15),W(1,IW),
     +              T-TOLD(TMP-1),KOLD)

C ... Compute delay value using the special interpolant
            ELSE 
               DO 135 J = 1, N
                  W(J,IW) = W(J,1)
 135           CONTINUE
            END IF
         END IF

 125  CONTINUE

C ... Compute K2

      CALL FCN(N,ND,NW,TEMP,W(1,16),W(1,21),W(1,2))

C ... Compute K3 to KS
      DO 140 K = 3, S

         TEMP = T + CC(K)*H

C ... Compute Y values
         DO 150 I = 1, N
            W(I,16) = AA(K,1)*W(I,1)
            DO 160 J = 2, K - 1
               W(I,16) = W(I,16) + AA(K,J)*W(I,J)
 160        CONTINUE
            W(I,16) = Y(I) + H*W(I,16)
 150     CONTINUE

C ... For consistency with asymptotically valid estimate
         if (k. eq. 8) then
            do 165 i = 1, n
               w(i,12) = w(i,16)
 165        continue
         end if

C ... Compute arguments for delay values
         CALL SIGMA(N,ND,TEMP,W(1,16),W(1,15))

C ... Compute solution delay values
         DO 170 I = 1, MU

            IW = 20 + I

            IF (W(I,15) .GT. TEMP) THEN
               PRINT *, '**** WARNING:APFORM *****'
               PRINT *, 'AdvanceD Delay at Stage', K
               PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
               PRINT *

C ... Cut the stepseize to half and integrate again
               IF (W(I,15)-TEMP .GE. 1.D0) THEN
                  PRINT *, 'ADVANCED TOO MUCH !!!'
                  PRINT *
                  C(18) = C(18) / 2.D0
                  IF (C(18) .LT. C(13)) THEN
                     PRINT *, 'IND = -3,  At T = ', T
                     STOP
                  END IF
                  C(17) = T + C(18)
                  C(24) = C(24) + DFLOAT(K - 2)
                  H = C(18)
                  GOTO 1000
               END IF
            END IF

            SIGN = W(I,15) - T
            TEMP2 = SIGN**2 / (CC(K)*H)

C ... In case there is no vanishing delay
            IF (SIGN .LE. 0.D0) THEN
               CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +              TOLD,KOLD,W,YINIT)

C ... In case there is vanishing delay
            ELSE
               IF (INDEX .EQ. 0) INDEX = K

C ... compute delay values using extrapolation
               IF (C(29) .EQ. 0.D0) THEN
	          TMP = IDINT(C(27))
                  IF (TMP .EQ. 1) TMP = LOLD
                  CALL INTRPG(N,NW,TMP,TOLD(TMP-1),W(I,15),W(1,IW),
     +                 T-TOLD(TMP-1),KOLD)

C ... Compute delay value using the special interpolant
               ELSE 
                  DO 180 L = 1, N
                     W(L,IW) = AA(K,1) * W(L,1)
                     DO 190 J = 2, K-1
                        W(L,IW) = W(L,IW) + AA(K,J) * W(L,J)
 190                 CONTINUE
                     W(L,IW) = Y(L) + (SIGN - TEMP2) * W(L,1) 
     +                    + TEMP2 * W(L,IW) / CC(K)

 180              CONTINUE
               END IF  
            END IF
 170     CONTINUE

C ... Compute derivative delay terms
         DO 175 I = MU + 1, ND

            IW = 20 + I

            IF (W(I,15) .GT. TEMP) THEN
               PRINT *, '**** WARNING:APFORM *****'
               PRINT *, 'AdvanceD Delay at Stage', K
               PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
               PRINT *

C ... Cut the stepseize to half and integrate again
               IF (W(I,15)-TEMP .GE. 1.D0) THEN
                  PRINT *, 'ADVANCED TOO MUCH !!!'
                  PRINT *
                  C(18) = C(18) / 2.D0
                  IF (C(18) .LT. C(13)) THEN
                     PRINT *, 'IND = -3,  At T = ', T
                     STOP
                  END IF
                  C(17) = T + C(18)
                  C(24) = C(24) + DFLOAT(K - 2)
                  H = C(18)
                  GOTO 1000
               END IF
            END IF

            SIGN = W(I,15) - T

C ... In case there is no vanishing delay
            IF (SIGN .LE. 0.D0) THEN
               CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +              TOLD,KOLD,W,DYINIT)

C ... In case there is vanishing delay
            ELSE
               IF (INDEX .EQ. 0) INDEX = K

C ... compute delay values using extrapolation
               IF (C(29) .EQ. 0.D0) THEN
	          TMP = IDINT(C(27))
                  IF (TMP .EQ. 1) TMP = LOLD
                  CALL DERIVG(N,NW,TMP,TOLD(TMP-1),W(I,15),W(1,IW),
     +                 T-TOLD(TMP-1),KOLD)

C ... Compute delay value using the special interpolant
               ELSE 

                  DO 185 L = 1, N
                     W(L,IW) = AA(K,1) * W(L,1)
                     DO 195 J = 2, K-1
                        W(L,IW) = W(L,IW) + AA(K,J) * W(L,J)
 195                 CONTINUE
                     W(L,IW) = W(L,1) - (W(L,1) - W(L,IW) / CC(K)) * 
     +                    (2.D0 * SIGN / (H * CC(K)))
 185              CONTINUE
               END IF  
            END IF
 175        CONTINUE

C ... Compute K-value

         CALL FCN(N,ND,NW,TEMP,W(1,16),W(1,21),W(1,K))

 140  CONTINUE

C***** End of the 1st phase *****

C***** Beginning of the 2nd phase *****

C ... Set # of iteration to be 0
      NIT = 0

C ... If index is not 0, then there is a vanishg delay 
      IF (INDEX .NE. 0) THEN

C ... Loop for iterative scheme ...
 2000    CONTINUE

         NIT = NIT + 1
         ERRKI = 0.D0

C ... Recompute K-values for stage INDEX+1 to s
         DO 210 K = INDEX, S

            DO 220 I = 1, N
               W(I,13) = W(I,K)
 220        CONTINUE

            TEMP = T + CC(K)*H

C ... Recompute y values
            DO 230 I = 1, N
               W(I,16) = AA(K,1)*W(I,1)
               DO 240 J = 2, K-1
                  W(I,16) = W(I,16) + AA(K,J)*W(I,J)
 240           CONTINUE
               W(I,16) = Y(I) + H*W(I,16)
 230        CONTINUE

C ... For consistency with Asymptotically valid estimate
            if (k .eq. 8) then
               do 235 i = 1, n
                  w(i,12) = w(i,16)
 235           continue
            end if

C ... Recompute arguments for delay values
            CALL SIGMA(N,ND,TEMP,W(1,16),W(1,15))

C ... Recompute solution delay values using intrp routine
            DO 250 I = 1, MU

               IW = 20 + I

               IF (W(I,15) .GT. TEMP) THEN
                  PRINT *, '**** WARNING:APFORM *****'
                  PRINT *, '    Advanced Delay at Stage ', K
                  PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
                  PRINT *

C ... Cut the stepseize to half and integrate again
                  IF (W(I,15)-TEMP .GE. 1.D0) THEN
                     PRINT *, 'ADVANCED TOO MUCH !!!'
                     PRINT *
                     C(18) = C(18) / 2.D0
                     IF (C(18) .LT. C(13)) THEN
                        PRINT *, 'IND = -3,  At T = ', T
                        STOP
                     END IF
                     C(17) = T + C(18)
                     C(24) = C(24) + DFLOAT(S-1 + K-INDEX)
                     H = C(18)
                     GOTO 1000
                  END IF
               END IF

               SIGN = W(I,15) - T

C ... In case there is no vanisning delay
               IF (SIGN .LE. 0.D0) THEN
                  CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +                 TOLD,KOLD,W,YINIT)
               ELSE 
                  CALL INTRP(N,NW,T,Y,W(I,15),W(1,IW),H,W)
               END IF
 250        CONTINUE

C ... Recompute derivative delay values using intrp routine
            DO 255 I = MU + 1, ND

               IW = 20 + I

               IF (W(I,15) .GT. TEMP) THEN
                  PRINT *, '**** WARNING:APFORM *****'
                  PRINT *, '    Advanced Delay at Stage ', K
                  PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
                  PRINT *

C ... Cut the stepseize to half and integrate again
                  IF (W(I,15)-TEMP .GE. 1.D0) THEN
                     PRINT *, 'ADVANCED TOO MUCH !!!'
                     PRINT *
                     C(18) = C(18) / 2.D0
                     IF (C(18) .LT. C(13)) THEN
                        PRINT *, 'IND = -3,  At T = ', T
                        STOP
                     END IF
                     C(17) = T + C(18)
                     C(24) = C(24) + DFLOAT(S-1 + K-INDEX)
                     H = C(18)
                     GOTO 1000
                  END IF
               END IF

               SIGN = W(I,15) - T

C ... In case there is no vanisning delay
               IF (SIGN .LE. 0.D0) THEN
                  CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +                 TOLD,KOLD,W,DYINIT)
               ELSE 

                  CALL DERIV(N,NW,T,Y,W(I,15),W(1,IW),H,W)
               END IF
 255        CONTINUE

C ... Recompute K-value

            CALL FCN(N,ND,NW,TEMP,W(1,16),W(1,21),W(1,K))

            DO 260 I = 1, N
               ERRKI = DMAX1(ERRKI, DABS(W(I,13) - W(I,K)))
 260        CONTINUE

 210     CONTINUE
         C(24) = C(24) + DFLOAT(S - INDEX + 1)

         IF ((ERRKI .GT. TOL) .AND. (NIT .LT. MAXNIT)) GOTO 2000 

      END IF
      
c      PRINT *, 'INDEX =', INDEX
c      PRINT *, 'NIT =', NIT
c      PRINT *, 'ERRKI =', ERRKI
c      PRINT *

C***** End of the the 2nd phase *****

C***** Beginning of the 3rd phase *****

C ... Improve derivative values on the mesh point for the interpolation
C ... to archieve C1 continuity
      IF ((INDEX .GT. 0) .AND. (INDEX .LE. 8)) THEN
         TEMP = T + CC(8)*H

C ... Compute Y values
         DO 310 I = 1, N
            W(I,12) = AA(8,1)*W(I,1)
            DO 320 J = 2, 7
               W(I,12) = W(I,12) + AA(8,J)*W(I,J)
 320        CONTINUE
            W(I,12) = Y(I) + H*W(I,12)
 310     CONTINUE

C ... Compute arguments for delay values
         CALL SIGMA(N,ND,TEMP,W(1,12),W(1,15))

C ... Compute solution delay values
         DO 330 I = 1, MU

            IW = 20 + I
            SIGN = W(I,15) - T
            IF (SIGN .LE. 0.D0) THEN
               CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +                      TOLD,KOLD,W,YINIT)
            ELSE 
               CALL INTRP(N,NW,T,Y,W(I,15),W(1,IW),H,W)
            END IF

 330     CONTINUE

C ... Compute derivative delay values
         DO 335 I = MU + 1, ND

            IW = 20 + I
            SIGN = W(I,15) - T
            IF (SIGN .LE. 0.D0) THEN
               CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +                      TOLD,KOLD,W,DYINIT)
            ELSE 

               CALL DERIV(N,NW,T,Y,W(I,15),W(1,IW),H,W)
            END IF
 335     CONTINUE

C ... Recompute K-value
         CALL FCN(N,ND,NW,TEMP,W(1,12),W(1,21),W(1,8))

         C(24) = C(24) + 1.D0
         C(28) = 0.D0
      END IF

C***** End of the 3nd phase *****

C ... Finally compute solution at T + H in W(*,12)
      DO 410 I = 1, N
         W(I,12) = BB(1)*W(I,1)
         DO 420 J = 2, S
            W(I,12) = W(I,12) + BB(J)*W(I,J)
 420     CONTINUE
         W(I,12) = Y(I) + H*W(I,12)
 410  CONTINUE

c      print *, 'index =', index

      RETURN
      END


C *** End of APFORM ***






      SUBROUTINE APFORM_ASM(N,NW,MU,OMEGA,LOLD,S,T,Y,H,TOLD,KOLD,TOL,W,C
     +     ,FCN,SIGMA,YINIT,DYINIT)

C**********************************************************************
C                                                                     *
C     Purpose - This routine integrates one step by RK-methods        *
c               for delay DEs.                                        *
C                                                                     *
C**********************************************************************
C                                                                     *
C     N      Number of equations                                      *
C     NW     First dimension of workspace W & queue KOLD              *
C     MU     Number of solution delay terms                           *
C     OMEGA  Number of derivative delay terms                         *
C     LOLD   The length of queue                                      *
C     S      The number of stages for RK-formula                      *
C     T      Current independent value                                *
C     Y      Current solution values, i.e., at T                      *
C     H      Trial stepsize                                           *
C     TOLD   Queue of T values                                        *
C     KOLD   Queue oF K values and solution values                    *
C     TOL    Tolerance                                                *
C     W      Workspace                                                *
C     C      Communication vector                                     *
C            ...USE C(13) : Minimum stepsize                          *
C                   C(17) : T trial                                   *
C                   C(18) : H trial                                   *
C                   C(24) : # of function evaluations                 *
C                   C(25) : Initial point                             *
C                   C(26) : The first index of queue                  *
C                   C(27) : The last index of queue                   *
C                   C(28) : Flag for computation of K1 in the next    *
C                           step(if = 0.0, need to compute)           *
C                   C(29) : Flag indicating which of the              *
c                           extrapolation or the special interpolant  *
c                           must be used.                             *
C                           If C(29) = 0.0 then use extrapolation     *
C     FCN    Name of derivative routine                               *
C     SIGMA  Name of subroutine for computing delay argument          *
C     YINIT  Name of subroutine for computing initial values at       *
C            initial intervals.                                       *
C     DYINIT Name of subroutine for computing derivative of initial   *
C            values at initial intervals.                             *
C                                                                     *
C**********************************************************************

C**********************************************************************
C Note:                                                               *
C       When using extrapolation, we use global intrp routine.        *
C       And use told(c(27)-1). But c(27) may be one.                  *
C       In this case, we have to use told(lold - 1)                   *
C**********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,MU,OMEGA,LOLD,S
      DOUBLE PRECISION T,H,TOL

C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N),TOLD(LOLD),KOLD(NW,7,LOLD),W(NW,*),C(*)

C     .. Subroutine Arguments ..
      EXTERNAL FCN,SIGMA,YINIT,DYINIT

C     .. Arrays in Common ..
      DOUBLE PRECISION AA(11,11),BB(11),CC(11)


C     .. Local Scalars ..
      INTEGER I,J,K,L,INDEX,IW,MAXNIT,ND,NIT,TMP
      DOUBLE PRECISION ERRKI,SIGN,TEMP,TEMP2

C     .. External Subroutines ..
      EXTERNAL EVALAGS,EVALAGD,INTRPG,DERIVG

C     .. Intrisic Functions ..
      INTRINSIC DMAX1,DFLOAT,IDINT

C     .. Common blocks ..
      COMMON/METH/AA,BB,CC

C     .. Data statements ..
      DATA MAXNIT/11/
c      DATA MAXNIT/5/

C     .. Set up the nubmer of delay terms
      ND = MU + OMEGA

C     ... LOOP for Advanced Delay ...
 1000 CONTINUE

C     Set flag to be 1
      C(28) = 1.D0
C
C***** Beginning of the 1st phase *****

C     Set index to 0
C ... INDEX: The first vanishing delay occures in the INDEX-th
C     stage if not 0.
      INDEX = 0

C     The 2nd stage
C ... Compute Y values
      DO 110 I = 1, N
         W(I,16) = Y(I) + H*AA(2,1)*W(I,1)
 110  CONTINUE

C ... Compute argument for delay values
      TEMP = T + CC(2)*H
      CALL SIGMA(N,ND,TEMP,W(1,16),W(1,15))

C ... Compute solution delay values
      DO 120 I = 1, MU

         IW = 20 + I

         IF (W(I,15) .GT. TEMP) THEN
            PRINT *, '**** WARNING:APFORM *****'
            PRINT *, 'Advanced Delay at Stage  2'
            PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
            PRINT *

C ... Cut the stepseize to half and integrate again
            IF (W(I,15) - TEMP .GE. 1.D0) THEN
               PRINT *, 'ADVANCED TOO MUCH !!!'
               PRINT *
               C(18) = C(18) / 2.D0
               IF (C(18) .LT. C(13)) THEN
                  PRINT *, 'IND = -3,  At T = ', T
                  PRINT *
                  STOP
               END IF
               C(17) = T + C(18)
               H = C(18)
               GOTO 1000
            END IF
         END IF
C
         SIGN = W(I,15) - T
C
C ... In case there is no vanishing delay
         IF (SIGN .LE. 0.D0) THEN
            CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +           TOLD,KOLD,W,YINIT)

C ... In case there is vanishing delay
         ELSE
            INDEX = 2

C ... Compute delay value using extrapolation
            IF (C(29) .EQ. 0.D0) THEN
	       TMP = IDINT(C(27))
               IF (TMP .EQ. 1) TMP = LOLD
               CALL INTRPG(N,NW,TMP,TOLD(TMP-1),W(I,15),W(1,IW),
     +              T-TOLD(TMP-1),KOLD)

C ... Compute delay value using the special interpolant
            ELSE 
               DO 130 J = 1, N
                  W(J,IW) = Y(J) +  SIGN*W(J,1)
 130           CONTINUE
            END IF
         END IF

 120  CONTINUE

C ... Compute derivative delay values      

      DO 125 I = MU + 1, ND

         IW = 20 + I

         IF (W(I,15) .GT. TEMP) THEN
            PRINT *, '**** WARNING:APFORM *****'
            PRINT *, 'Advanced Delay at Stage  2'
            PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
            PRINT *

C ... Cut the stepseize to half and integrate again
            IF (W(I,15) - TEMP .GE. 1.D0) THEN
               PRINT *, 'ADVANCED TOO MUCH !!!'
               PRINT *
               C(18) = C(18) / 2.D0
               IF (C(18) .LT. C(13)) THEN
                  PRINT *, 'IND = -3,  At T = ', T
                  PRINT *
                  STOP
               END IF
               C(17) = T + C(18)
               H = C(18)
               GOTO 1000
            END IF
         END IF
C
         SIGN = W(I,15) - T
C
C ... In case there is no vanishing delay
         IF (SIGN .LE. 0.D0) THEN
            CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +           TOLD,KOLD,W,DYINIT)

C ... In case there is vanishing delay
         ELSE
            INDEX = 2

C ... Compute delay value using extrapolation
            IF (C(29) .EQ. 0.D0) THEN
	       TMP = IDINT(C(27))
               IF (TMP .EQ. 1) TMP = LOLD
               CALL DERIVG(N,NW,TMP,TOLD(TMP-1),W(I,15),W(1,IW),
     +              T-TOLD(TMP-1),KOLD)

C ... Compute delay value using the special interpolant
            ELSE 
               DO 135 J = 1, N
                  W(J,IW) = W(J,1)
 135           CONTINUE
            END IF
         END IF

 125  CONTINUE

C ... Compute K2

      CALL FCN(N,ND,NW,TEMP,W(1,16),W(1,21),W(1,2))

C ... Compute K3 to KS
      DO 140 K = 3, S

         TEMP = T + CC(K)*H

C ... Compute Y values
         DO 150 I = 1, N
            W(I,16) = AA(K,1)*W(I,1)
            DO 160 J = 2, K - 1
               W(I,16) = W(I,16) + AA(K,J)*W(I,J)
 160        CONTINUE
            W(I,16) = Y(I) + H*W(I,16)
 150     CONTINUE

c .. for assmyptocally correct defect control
         if (k .eq. 8) then
            do 165 i = 1, n
               w(i,12) = w(i,16)
 165        continue
         end if

C ... Compute arguments for delay values
         CALL SIGMA(N,ND,TEMP,W(1,16),W(1,15))

C ... Compute solution delay values
         DO 170 I = 1, MU

            IW = 20 + I

            IF (W(I,15) .GT. TEMP) THEN
               PRINT *, '**** WARNING:APFORM *****'
               PRINT *, 'AdvanceD Delay at Stage', K
               PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
               PRINT *

C ... Cut the stepseize to half and integrate again
               IF (W(I,15)-TEMP .GE. 1.D0) THEN
                  PRINT *, 'ADVANCED TOO MUCH !!!'
                  PRINT *
                  C(18) = C(18) / 2.D0
                  IF (C(18) .LT. C(13)) THEN
                     PRINT *, 'IND = -3,  At T = ', T
                     STOP
                  END IF
                  C(17) = T + C(18)
                  C(24) = C(24) + DFLOAT(K - 2)
                  H = C(18)
                  GOTO 1000
               END IF
            END IF

            SIGN = W(I,15) - T
            TEMP2 = SIGN**2 / (CC(K)*H)

C ... In case there is no vanishing delay
            IF (SIGN .LE. 0.D0) THEN
               CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +              TOLD,KOLD,W,YINIT)

C ... In case there is vanishing delay
            ELSE
               IF (INDEX .EQ. 0) INDEX = K

C ... compute delay values using extrapolation
               IF (C(29) .EQ. 0.D0) THEN
	          TMP = IDINT(C(27))
                  IF (TMP .EQ. 1) TMP = LOLD
                  CALL INTRPG(N,NW,TMP,TOLD(TMP-1),W(I,15),W(1,IW),
     +                 T-TOLD(TMP-1),KOLD)

C ... Compute delay value using the special interpolant
               ELSE 
                  DO 180 L = 1, N
                     W(L,IW) = AA(K,1) * W(L,1)
                     DO 190 J = 2, K-1
                        W(L,IW) = W(L,IW) + AA(K,J) * W(L,J)
 190                 CONTINUE
                     W(L,IW) = Y(L) + (SIGN - TEMP2) * W(L,1) 
     +                    + TEMP2 * W(L,IW) / CC(K)
 180              CONTINUE
               END IF  
            END IF
 170     CONTINUE

C ... Compute derivative delay terms
         DO 175 I = MU + 1, ND

            IW = 20 + I

            IF (W(I,15) .GT. TEMP) THEN
               PRINT *, '**** WARNING:APFORM *****'
               PRINT *, 'AdvanceD Delay at Stage', K
               PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
               PRINT *

C ... Cut the stepseize to half and integrate again
               IF (W(I,15)-TEMP .GE. 1.D0) THEN
                  PRINT *, 'ADVANCED TOO MUCH !!!'
                  PRINT *
                  C(18) = C(18) / 2.D0
                  IF (C(18) .LT. C(13)) THEN
                     PRINT *, 'IND = -3,  At T = ', T
                     STOP
                  END IF
                  C(17) = T + C(18)
                  C(24) = C(24) + DFLOAT(K - 2)
                  H = C(18)
                  GOTO 1000
               END IF
            END IF

            SIGN = W(I,15) - T

C ... In case there is no vanishing delay
            IF (SIGN .LE. 0.D0) THEN
               CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +              TOLD,KOLD,W,DYINIT)

C ... In case there is vanishing delay
            ELSE
               IF (INDEX .EQ. 0) INDEX = K

C ... compute delay values using extrapolation
               IF (C(29) .EQ. 0.D0) THEN
	          TMP = IDINT(C(27))
                  IF (TMP .EQ. 1) TMP = LOLD
                  CALL DERIVG(N,NW,TMP,TOLD(TMP-1),W(I,15),W(1,IW),
     +                 T-TOLD(TMP-1),KOLD)

C ... Compute delay value using the special interpolant
               ELSE 

                  DO 185 L = 1, N
                     W(L,IW) = AA(K,1) * W(L,1)
                     DO 195 J = 2, K-1
                        W(L,IW) = W(L,IW) + AA(K,J) * W(L,J)
 195                 CONTINUE
                     W(L,IW) = W(L,1) - (W(L,1) - W(L,IW) / CC(K)) * 
     +                    (2.D0 * SIGN / (H * CC(K)))
 185              CONTINUE
               END IF  
            END IF
 175        CONTINUE

C ... Compute K-value
         CALL FCN(N,ND,NW,TEMP,W(1,16),W(1,21),W(1,K))

 140  CONTINUE

C***** End of the 1st phase *****

C***** Beginning of the 2nd phase *****

C ... Set # of iteration to be 0
      NIT = 0

C ... If index is not 0, then there is a vanishg delay 
      IF (INDEX .NE. 0) THEN

C ... Loop for iterative scheme ...
 2000    CONTINUE

         NIT = NIT + 1
         ERRKI = 0.D0

C ... Recompute K-values for stage INDEX+1 to s
         DO 210 K = INDEX, S

            DO 220 I = 1, N
               W(I,13) = W(I,K)
 220        CONTINUE

            TEMP = T + CC(K)*H

C ... Recompute y values
            DO 230 I = 1, N
               W(I,16) = AA(K,1)*W(I,1)
               DO 240 J = 2, K-1
                  W(I,16) = W(I,16) + AA(K,J)*W(I,J)
 240           CONTINUE
               W(I,16) = Y(I) + H*W(I,16)
 230        CONTINUE

c .. for assmyptocally correct defect control

            if (k .eq. 8) then
               do 235 i = 1, n
                  w(i,12) = w(i,16)
 235           continue
            end if

C ... Recompute arguments for delay values
            CALL SIGMA(N,ND,TEMP,W(1,16),W(1,15))

C ... Recompute solution delay values using intrp routine
            DO 250 I = 1, MU

               IW = 20 + I

               IF (W(I,15) .GT. TEMP) THEN
                  PRINT *, '**** WARNING:APFORM *****'
                  PRINT *, '    Advanced Delay at Stage ', K
                  PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
                  PRINT *

C ... Cut the stepseize to half and integrate again
                  IF (W(I,15)-TEMP .GE. 1.D0) THEN
                     PRINT *, 'ADVANCED TOO MUCH !!!'
                     PRINT *
                     C(18) = C(18) / 2.D0
                     IF (C(18) .LT. C(13)) THEN
                        PRINT *, 'IND = -3,  At T = ', T
                        STOP
                     END IF
                     C(17) = T + C(18)
                     C(24) = C(24) + DFLOAT(S-1 + K-INDEX)
                     H = C(18)
                     GOTO 1000
                  END IF
               END IF

               SIGN = W(I,15) - T

C ... In case there is no vanisning delay
               IF (SIGN .LE. 0.D0) THEN
                  CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +                 TOLD,KOLD,W,YINIT)
               ELSE 
                  CALL INTRP(N,NW,T,Y,W(I,15),W(1,IW),H,W)
               END IF
 250        CONTINUE

C ... Recompute derivative delay values using intrp routine
            DO 255 I = MU + 1, ND

               IW = 20 + I

               IF (W(I,15) .GT. TEMP) THEN
                  PRINT *, '**** WARNING:APFORM *****'
                  PRINT *, '    Advanced Delay at Stage ', K
                  PRINT *, 'Ti =', TEMP, ' Delay Argument =', W(I,15)
                  PRINT *

C ... Cut the stepseize to half and integrate again
                  IF (W(I,15)-TEMP .GE. 1.D0) THEN
                     PRINT *, 'ADVANCED TOO MUCH !!!'
                     PRINT *
                     C(18) = C(18) / 2.D0
                     IF (C(18) .LT. C(13)) THEN
                        PRINT *, 'IND = -3,  At T = ', T
                        STOP
                     END IF
                     C(17) = T + C(18)
                     C(24) = C(24) + DFLOAT(S-1 + K-INDEX)
                     H = C(18)
                     GOTO 1000
                  END IF
               END IF

               SIGN = W(I,15) - T

C ... In case there is no vanisning delay
               IF (SIGN .LE. 0.D0) THEN
                  CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +                 TOLD,KOLD,W,DYINIT)
               ELSE 

                  CALL DERIV(N,NW,T,Y,W(I,15),W(1,IW),H,W)
               END IF
 255        CONTINUE

C ... Recompute K-value

            CALL FCN(N,ND,NW,TEMP,W(1,16),W(1,21),W(1,K))

            DO 260 I = 1, N
               ERRKI = DMAX1(ERRKI, DABS(W(I,13) - W(I,K)))
 260        CONTINUE

 210     CONTINUE
         C(24) = C(24) + DFLOAT(S - INDEX + 1)

        IF ((ERRKI .GT. TOL) .AND. (NIT .LT. MAXNIT)) GOTO 2000 

      END IF

c      PRINT *, 'INDEX =', INDEX
c      PRINT *, 'NIT =', NIT
c      PRINT *, 'ERRKI =', ERRKI
c      PRINT *

C***** End of the the 2nd phase *****

C***** Beginning of the 3rd phase *****

C ... Improve derivative values on the mesh point for the interpolation
C ... to archieve C1 continuity
      IF ((INDEX .GT. 0) .AND. (INDEX .LE. 8)) THEN
         TEMP = T + CC(8)*H

C ... Compute Y values
         DO 310 I = 1, N
            W(I,12) = AA(8,1)*W(I,1)
            DO 320 J = 2, 7
               W(I,12) = W(I,12) + AA(8,J)*W(I,J)
 320        CONTINUE
            W(I,12) = Y(I) + H*W(I,12)
 310     CONTINUE


C ... Compute arguments for delay values
         CALL SIGMA(N,ND,TEMP,W(1,12),W(1,15))

C ... Compute solution delay values
         DO 330 I = 1, MU

            IW = 20 + I
            SIGN = W(I,15) - T

            IF (SIGN .LE. 0.D0) THEN
               CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +                      TOLD,KOLD,W,YINIT)
            ELSE 
               CALL INTRP(N,NW,T,Y,W(I,15),W(1,IW),H,W)

            END IF
 330     CONTINUE

C ... Compute derivative delay values
         DO 335 I = MU + 1, ND

            IW = 20 + I
            SIGN = W(I,15) - T

            IF (SIGN .LE. 0.D0) THEN
               CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +                      TOLD,KOLD,W,DYINIT)
            ELSE 

               CALL DERIV(N,NW,T,Y,W(I,15),W(1,IW),H,W)
            END IF
 335     CONTINUE

C ... Recompute K-value
         CALL FCN(N,ND,NW,TEMP,W(1,12),W(1,21),W(1,8))

         C(24) = C(24) + 1.D0
         C(28) = 0.D0
      END IF

C***** End of the 3nd phase *****



C***** Beginning of 4th phase (Improve intrpolant for defect) *****
      do 500 k = 9, s
         temp = t + cc(k)*h

C ... Compute Y values
         call intrp(n,nw,t,y,temp,w(1,16),h,w)

C ... Compute arguments for delay values
         call sigma(n,nd,temp,w(1,16),w(1,15))

C ... Compute solution delay values
         do 510 i = 1, mu

            iw = 20 + i
            sign = w(i,15) - t

            if (sign .le. 0.d0) then
               call evalags(n,nw,iw,lold,c(25),c(26),c(27),w(i,15),
     +                      told,kold,w,yinit)
            else
               call intrp(n,nw,t,y,w(i,15),w(1,iw),h,w)
            end if
 510     continue

C ... Compute derivative delay values
         do 520 i = mu + 1, nd

            iw = 20 + i
            sign = w(i,15) - t

            if (sign .le. 0.d0) then
               call evalagd(n,nw,iw,lold,c(25),c(26),c(27),w(i,15),
     +                      told,kold,w,dyinit)
            else 
               call deriv(n,nw,t,y,w(i,15),w(1,iw),h,w)
            end if
 520     continue

C ... Recompute K-value
         call fcn(n,nd,nw,temp,w(1,16),w(1,21),w(1,k))

         c(24) = c(24) + 1.d0
 500  continue

C ... Finally compute solution at T + H in W(*,12)
      DO 410 I = 1, N
         W(I,12) = BB(1)*W(I,1)
         DO 420 J = 2, S
            W(I,12) = W(I,12) + BB(J)*W(I,J)
 420     CONTINUE
         W(I,12) = Y(I) + H*W(I,12)
 410  CONTINUE

c      print *, 'index =', index

c      print '(3d25.14)', (w(1,i),i=1,12)
c      print *
c      print '(3d25.14)', (w(2,i),i=1,12)
c      print *
c      print '(3d25.14)', (w(3,i),i=1,12)      
c      print *

      RETURN
      END


C *** End of APFORM_ASM ***








      SUBROUTINE DDVERK(N,NW,MU,OMEGA,IND,LOLD,T,TEND,Y,TOL,TOLD,KOLD
     +     ,W,C,FCN,SIGMA,YINIT,DYINIT)

C************************************************************************
C                                                                       *
C            Description of DDVERK (Experimental Code)                  *
C                                                                       *
C     DDVERK is a driver which invokes delay and discontinuity routines.*
C     This driver first initializes parameters and then calls the DVERK *
C     integration routine, RDMETH. After each step of the integration,  *
C     RDMETH returns to this driver and this driver checks if a         *
C     discontinuity exists in the step. If the existence of a           *
C     discontinuity is suspected, this driver calls DSLOC and           *
C     DSLOC detects the discontinuity point. Finally, if the step is    *
C     accepted, DDVERK calls STORE or STORED to store information       *
C     necessary to construct an interpolant in future.                  *
C                                                                       *
C************************************************************************
C                                                                       *
C                          Calling sequence                             *
C                                                                       *
C     CALL DDVERK(N,NW,MU,OMEGA,IND,LOLD,T,TEND,Y,TOL,TOLD,KOLD,W,C,    *
C    +            FCN,SIGMA,YINIT,DYINIT)                               *
C                                                                       *
C     INTEGER N,NW,IND,LOLD,MU,OMEGA                                    *
C     DOUBLE PRECISION T,TEND,Y(N),TOL,TOLD(LOLD),KOLD(NW,7,LOLD)       *
C     DOUBLE PRECISION W(NW,20+MU+OMEGA),C(*)                          *
C     EXTERNAL FCN,SIGMA,YINIT,DYINIT                                   *
C                                                                       *
C************************************************************************
C                                                                       *
C                       Arguments for DDERK                             *
C                                                                       *
C    N      Number of equations                                         *
C    NW     First dimension of workspace W and queue KOLD               *
C           NW must be at least MAX(N,MU+OMEGA)                         *
C    MU     The number of solution delay terms                          *
C    OMEGA  The number of derivative delay terms                        *
C    IND    Indicator --- same as DVERK except for IND = 7:             *
C           IND is set to be 7 after a discontinuity is found.          *
C    LOLD   Length of queue                                             *
C    T      Independent value                                           *
C    TEND   Value of T to which integration is to be carried out        *
C    Y      Solution vector at T                                        *
C    TOL    Tolerance                                                   *
C    TOLD   Queue of independent values T                               *
C    KOLD   Queue of derivative and solution approximations             *
C           The first dimension of KOLD must be NW                      *
C           The second dimension of KOLD must be at least 7             *
C           The third dimension of KOLD must be LOLD                    *
C    W      Workspace array                                             *
C           The first dimension of W must be NW                         *
C           The second dimension of W must be at least 20 + ND          *
C    C      Communication vector - The dimension must be greater than   *
C           or equal to 33. C(1) to C(24) are the same as DVERK.        *
C    FCN    Name of subroutine for computing derivative                 *
C           (must be supplied by a user)                                *
C    SIGMA  Name of subroutine for computing delay argument             *
C           (must be supplied by a user)                                *
C    YINIT  Name of subroutine for computing initial values on the      *
C           initial interval. (must be supplied by a user)              *
C    DYINIT Name of subroutine for computing derivative of initial      *
C           values on the  initial interval. (must be supplied by       *
C           a user)                                                     *
C                                                                       *
C************************************************************************
C                                                                       *
C                          Communication vector                         *
C                                                                       *
C     C(25) Initial time                                                *
C     C(26) The first index of queue                                    *
C     C(27) The last index of queue                                     *
C     C(28) Flag indicating if K1 needs to be computed in the next      *
C           step: if C(28) = 0.0 then we need to compute K1 in the next *
C           step                                                        *
C     C(29) Flag indicating which of the extrapolation or the special   *
C           interpolant must be used: if C(29) = 0.0 then we use        *
C           extrapolation                                               *
C     C(30) TL values used in discontinuity locating routine DSLOC.     *
C     C(31) TH - TL used in DSLOC                                       *
C     C(32) Flag indicating which of the defect estimate must be used:  *
C              1: 1pt, 2: 2pt, 3: asymptotically valid                  *
C              The default value is 1.                                  *
c     C(33) TS value (location of disconinuous point) used in ddverk    *
C                                                                       *
C                                                                       *
C************************************************************************
C                                                                       *
C                               Work array                              *
C                                                                       *
C     W(*, 1:11)           K1 to K11 of CRK formula                     *
C     W(*, 12)             Solution at the right end point of the       *
C                          each step                                    *
C     W(*, 13)             Work space                                   *
C     W(*, 14)             Estimate of Defect                           *
C     W(*, 15)             Delay argument                               *
C     W(*, 16)             Work space                                   *
C     W(*, 17)             Work space                                   *
C     W(*, 18)             Work space                                   *
C     W(*, 21:20+ND)       Delay values                                 *
C                                                                       *
C************************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,MU,OMEGA,IND,LOLD
      DOUBLE PRECISION T,TEND,TOL

C     .. Array Arguments ..
      DOUBLE PRECISION Y(N),TOLD(LOLD),KOLD(NW,7,LOLD)
      DOUBLE PRECISION W(NW,*),C(*)

C     .. Subroutine Arguments ..
      EXTERNAL FCN,SIGMA,YINIT,DYINIT

C     .. Local Scalars ..
      INTEGER I,IND2,P
      DOUBLE PRECISION TS,K, h2, h3

C    .. Local Logicals ..
      LOGICAL GOBACK

C     .. External Subroutines ..
      EXTERNAL RDMETH,DSLOC,INTRP,STORED,STORE

c .. scalr in common
      integer disc_count

C     .. Common blocks ..
      COMMON/GOBACK/GOBACK
      common/disc_c/disc_count

C     .. Data Statements .. 
C ... Set the order of methods to be 6 
      DATA P/6/
C ... GOBACK: Logical variable 
C             --- if true, then return to main routine in every step
      DATA GOBACK/.FALSE./

C ... Initialize communication array

      IF (IND .EQ. 1) THEN
         IND = 2
         DO 110 I = 1, 9
            C(I) = 0.D0
 110     CONTINUE
         C(32) = 1
      END IF


      IF (IND .EQ. 2) THEN
C ... Set up initial value
         CALL YINIT(N,T,Y)

C ... In the first step, the extrapolation cannot be used. Hence set
C     C(29) to be 1.
         C(29) = 1.D0

         GOBACK = C(9) .EQ. 1.D0
C ... Has to be returned to this routine in every step
         C(9) = 1.D0

C ... Substitue initial time to C(25) and TOLD(1)
         C(25) = T
         TOLD(1) = T

C ... Set the first index of the history queue to be 1
         C(26) = 1.D0

C ... Set the last index of the history queue to be 1
         C(27) = 1.D0

C ... Set flag, which determins if K1 (in the current step) is 
C     equal to k8 (in the previous step), to 0.d0.
C     If FLAG(=C(28)) is 0, then K1 != K8.
         C(28) = 0.D0
      
      END IF

C ... Loop for integration till T = Tend
 1000 CONTINUE

C ... Change the information when IND = 7
      IF (IND .EQ. 7) THEN
         T = C(17)
         DO 210 I = 1, N
            Y(I) = W(I,12)
 210     CONTINUE

C ... Add 1 to the number of successful steps
         C(22) = C(22) + 1.D0

C ... Set the number of successive failures to be 0
         C(23) = 0.D0

C ... If T >= Tend, then end the integration
         IF (T .GE. TEND) THEN
            IND = 3
            C(20) = TEND
            C(21) = 1.D0         
            RETURN
         END IF
      END IF

 2000 continue
      CALL RDMETH(N,NW,MU,OMEGA,IND,LOLD,T,TEND,Y,TOLD,KOLD,TOL,W,C,
     +            FCN,SIGMA,YINIT,DYINIT) 

C ... Exit successfully if IND = 3
      IF (IND .EQ. 3) RETURN

C ... Error exit if IND < 1 or IND > 8
      IF ((IND .LT. 1) .OR. (IND .GT. 9)) GOTO 500

C ... Check for signal of possible discontinuity -- the signal by
C     an attempt to 1.6 times the step when the immediately preceding 
C     step was a failure

      IF ((C(23) .GT. 0.D0) .AND. (IND .EQ. 5)
     +     .AND. (TOL .GE. 1.6d0**P*C(19)/0.9D0)) THEN

C ... Use defect to locate and cross discontinuity if one exist
         CALL DSLOC(N,NW,MU,OMEGA,IND2,LOLD,T,TS,Y,W(1,13),TOLD,KOLD,
     +              TOL,W,C,FCN,SIGMA,YINIT,DYINIT)

C ... Now cross the discontinuity if one has been found (IND2 > 0)
C     and prepare for a restart of the integrator.

         print *, 'ind2 =', ind2

         IF (IND2 .GT. 0) THEN
            PRINT *
            PRINT *, 'DETECT DISCONTINUITY AT T = ', TS

            IF (TS .LT. C(17)) THEN
               PRINT *, '***** ERROR AT DDVERK *****'
               PRINT *, 'TS is less than Tn-1 + H'
               PRINT *
               STOP
            END IF

C ... In case the discontinous point is greater than or equal to Tend,
C     we have to inerpolate tend and store informaiton and then return
C     to main routine.

            IF ((TS .GE. TEND) .OR. (DABS(TEND - TS) .LE. C(10))) THEN
               CALL INTRP(N,NW,T,Y,TEND,W(1,13),C(18),W)
               CALL STORED(N,NW,LOLD,T,TEND,Y,W(1,13),TOLD,KOLD,W,C)
               C(17) = TEND

               DO 310 I = 1, N
                  W(I,12) = W(I,13)
 310           CONTINUE

               write(13, 10000) ts
               disc_count = disc_count + 1

               RETURN
            END IF

c ... Try not to use extration, so extend the stepsize and integrate
c ... Note: c(30) = XL and c(33) = XH - XL
            h3 = ts - t
            h2 = dmin1(h3/1.1d0, c(30) - t)

c            print *, 'h3/h1=', h3/c(18), ', h2-h1/h1 =', 
c     +           (h2 - c(18))/c(18), ', hl =', c(30)
            if ( (h3/c(18) .gt. 1.2d0) .and. (dabs(h2 - c(18))/c(18)
     +           .gt. 0.1d0)) then

c ... Save info in case the next step is rejected

               CALL STORE(N,NW,LOLD,C(17),Y,C(26),C(27),TOLD,KOLD,W)
               c(27) = c(27) - 1.d0
               ind = 8
               c(18) = h2
               c(17) = t + c(18)
               c(14) = c(18)
               c(33) = ts
               go to 2000

            end if
 
C ... If there is a discontinuity, the extrapolation cannot be used.
C     Hence set C(29) to be 1.
            C(29) = 1.D0

            IND = 7
            K = C(23) - 1.D0
            C(4) = C(18) * 2.D0**K
            C(28) = 0.D0         
            c(23) = 0.d0
   
            CALL STORED(N,NW,LOLD,T,TS,Y,W(1,13),TOLD,KOLD,W,C)
            C(17) = TS

            DO 320 I = 1,N
               W(I,12) = W(I,13)
 320        CONTINUE

            write(13, 10000) ts
            disc_count = disc_count + 1

         ELSE 

            CALL STORE(N,NW,LOLD,C(17),Y,C(26),C(27),TOLD,KOLD,W) 

         END IF

      ELSE IF (IND .EQ. 5) THEN

         CALL STORE(N,NW,LOLD,C(17),Y,C(26),C(27),TOLD,KOLD,W)

      ELSE IF (IND .EQ. 8) THEN

         c(29) = 1.d0
         ind = 7
         k = c(23) - 1.d0
         c(4) = c(18)
c         c(4) = c(18) * 2.d0 ** k
         c(28) = 0.d0
         c(23) = 0.d0
         call intrp(n,nw,t,y,c(33),w(1,13),c(18),w)
         call stored(n,nw,lold,t,c(33),y,w(1,13),told,kold,w,c)
         c(17) = c(33)
         do 350 i = 1, n
            w(i,12) = w(i,13)
 350     continue

         write(13, 10000) c(33)
         disc_count = disc_count + 1
         
      else if (ind .eq. 9) then
         
         c(27) = c(27) + 1.d0
         c(23) = 0.d0
         ind = 5
         do 370 i = 1, n
            w(i,12) = kold(i,3,idint(c(27)))
            w(i,8) = kold(i,4,idint(c(27)))
 370     continue
         c(17) = told(idint(c(27)))
         c(18) = c(17) - told(idint(c(27)) - 1.d0)
         c(14) = c(18)
         c(29) = 0.d0

      end if

 500  CONTINUE

      IF (GOBACK) RETURN

      GO TO 1000

10000 format(f20.12)

      END


C *** End of DDVERK ***





      SUBROUTINE DEFECT(N,NW,MU,OMEGA,LOLD,T,TS,Y,DFCT,TOLD,KOLD,W,C,
     +                  FCN,SIGMA,YINIT,DYINIT)

C********************************************************************
C                                                                   *
C     Purpose - This subroutine computes weighted max norm of       *
c               defect at the point TS.                             *
C                                                                   *
C********************************************************************
C                                                                   *
C    N      Number of equation                                      *
C    NW     First dimension of workspace w & queue kold             *
C    MU     Number of solution delay terms                          *
C    OMEGA  Number of derivative delay terms                        *
C    LOLD   Length of queue                                         *
C    T      Independent value (left end point)                      * 
C    TS     Independent value at which defect is to be computed     *
C    Y      Dependent values (at T)                                 *
C    DFCT   Max norm of defect at ts (output)                       *
C    TOLD   Queue of independent values                             *
C    KOLD   Queue of k-values and solutions                         *
C    W      Workspace                                               *
C    C      Communication array                                     *
C    FCN    Name of subroutine for computing derivative             *
C    SIGMA  Name of subroutine for computing delay argument         *
C    YINIT  Name of subroutine for computing initial values at      *
C           the initial intervals.                                  *
C    DYINIT Name of subroutine for computing derivative of initial  *
C           values at the initial intervals.                        *
C                                                                   *
C********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,MU,OMEGA,LOLD
      DOUBLE PRECISION T,TS,DFCT

C     .. Array Arguments ..
      DOUBLE PRECISION Y(N),TOLD(LOLD),KOLD(NW,7,LOLD),W(NW,*),C(*)

C     .. Subroutine Arguments ..
      EXTERNAL FCN,SIGMA,YINIT,DYINIT

C     .. Local Scalars ..
      INTEGER I,K,IW,ND
      DOUBLE PRECISION TEMP,SIGN,H

C     .. External Subroutines ..
      EXTERNAL INTRP,EVALAGS,EVALAGD,DERIV

C     .. Intrinsic Functions ..
      INTRINSIC DMAX1,DABS

c      print *, 't =', t, ' ts =', ts

C ... Error check
      IF (TS .LT. T) THEN
         PRINT *, '***** ERROR AT DEFECT *****'
         PRINT *, 'T=', T, ' TS=', TS
         PRINT *, 'TS is less than T !!!'
         PRINT *
         STOP
      END IF

      ND = MU + OMEGA
      H = C(18)

C ... Compute P(TS)
      CALL INTRP(N,NW,T,Y,TS,W(1,2),H,W)

C ... Compute argument of delay value
      CALL SIGMA(N,ND,TS,W(1,2),W(1,15))

C     Compute P(Delay Arugment)
      DO 110 I = 1, MU
         IW = 20 + I

         IF (W(I,15) .GT. TS) THEN
            PRINT *, '**** Warning:DEFECT *****'
            PRINT *, '    ADVANCED DELAY in DEFECT'
            PRINT *, 'T=',TS,'  Delay Argument=',W(I,15)
            PRINT *
         END IF

         SIGN = W(I,15) - T

C ... In case there is no vanishing delay
         IF (SIGN .LE. 0.D0) THEN
            CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +           TOLD,KOLD,W,YINIT)

C ... In case there is vanishing delay
         ELSE
            CALL INTRP(N,NW,T,Y,W(I,15),W(1,IW),H,W)
         END IF
 110  CONTINUE
       
C ...  Compute P'(Delay Argument)
      DO 115 I = MU + 1, ND

         IW = 20 + I

         IF (W(I,15) .GT. TS) THEN
            PRINT *, '**** Warning:DEFECT *****'
            PRINT *, '    ADVANCED DELAY in DEFECT'
            PRINT *, 'T=',TS,'  Delay Argument=',W(I,15)
            PRINT *
         END IF

         SIGN = W(I,15) - T

C ... In case there is no vanishing delay
         IF (SIGN .LE. 0.D0) THEN
            CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +           TOLD,KOLD,W,DYINIT)

C ... In case there is vanishing delay
         ELSE
            CALL DERIV(N,NW,T,Y,W(I,15),W(1,IW),H,W)
         END IF
 115  CONTINUE

C ... Compute F(TS, P(TS), P(DEL ARG), P'(DELAY ARG))
      CALL FCN(N,ND,NW,TS,W(1,2),W(1,21),W(1,18))
      C(24) = C(24) + 1.D0

C ... Compute P'(TS) and then defect
      CALL DERIV(N,NW,T,Y,TS,W(1,2),H,W)

c      print *, 'derv =', (w(i,2), i = 1, 3) 

      DO 120 I = 1, N
         W(I,14) = W(I,18) - W(I,2)
 120  CONTINUE

c      print *, 'def =', (w(i,14), i=1,3)

C ... Calculate the weighted max norm of W(*,14) as specified by
C     the error control indicator C(1)
      TEMP = 0.D0
      IF (C(1) .NE. 1.D0) GO TO 210

C ... Absolute error control
      DO 200 K = 1, N
         TEMP = DMAX1(TEMP,DABS(W(K,14)))
 200  CONTINUE
      GO TO 500
 210  IF (C(1) .NE. 2.D0) GO TO 240

C ... Relative error control
      DO 230 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)/W(K,18)))
 230  CONTINUE
      GO TO 500
 240  IF (C(1) .NE. 3.D0) GO TO 260

C ... Weights are 1/MAX(C(2),ABS(W(K,18)))
      DO 250 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)) / DMAX1(C(2), 
     +        DABS(W(K,18))))
 250  CONTINUE
      GO TO 500
 260  IF (C(1) .NE. 4.D0) GO TO 280

C ... Weights are 1/MAX(C(K+30),ABS(W(K,18)))
      DO 270 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)) / DMAX1(C(K+30), 
     +        DABS(W(K,18))))
 270  CONTINUE
      GO TO 500
 280  IF (C(1) .NE. 5.D0) GO TO 300

C ... Weights are 1/C(K+30)
      DO 290 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)/C(K+30)))
 290  CONTINUE
      GO TO 500
 300  CONTINUE

C ... Default case - Weights are 1/MAX(1,ABS(Y'(K)))
      DO 310 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)) / DMAX1(1.D0, 
     +        DABS(W(K,18))) )
 310  CONTINUE
 500  CONTINUE

      DFCT = TEMP

      RETURN
      END


C *** End of DEFECT ***





      SUBROUTINE DEFECT2(N,NW,MU,OMEGA,LOLD,T,TS,TS2,Y,DFCT,TOLD,KOLD,W,
     +                  C,FCN,SIGMA,YINIT,DYINIT)

C********************************************************************
C                                                                   *
C     Purpose - This subroutine computes weighted max norm of       *
c               defects at the point TS and TS2.                    *
C                                                                   *
C********************************************************************
C                                                                   *
C    N      Number of equation                                      *
C    NW     First dimension of workspace w & queue kold             *
C    MU     Number of solution delay terms                          *
C    OMEGA  Number of derivative delay terms                        *
C    LOLD   Length of queue                                         *
C    T      Independent value (left end point)                      * 
C    TS     Independent value at which defect is to be computed     *
C    TS2    Independent value at which defect is to be computed     *
C    Y      Dependent values (at T)                                 *
C    DFCT   Max norm of defect at ts (output)                       *
C    TOLD   Queue of independent values                             *
C    KOLD   Queue of k-values and solutions                         *
C    W      Workspace                                               *
C    C      Communication array                                     *
C    FCN    Name of subroutine for computing derivative             *
C    SIGMA  Name of subroutine for computing delay argument         *
C    YINIT  Name of subroutine for computing initial values at      *
C           the initial intervals.                                  *
C    DYINIT Name of subroutine for computing derivative of initial  *
C           values at the initial intervals.                        *
C                                                                   *
C********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,MU,OMEGA,LOLD
      DOUBLE PRECISION T,TS,TS2,DFCT

C     .. Array Arguments ..
      DOUBLE PRECISION Y(N),TOLD(LOLD),KOLD(NW,7,LOLD),W(NW,*),C(*)

C     .. Subroutine Arguments ..
      EXTERNAL FCN,SIGMA,YINIT,DYINIT

C     .. Local Scalars ..
      INTEGER I,K,IW,ND
      DOUBLE PRECISION TEMP,SIGN,H

C     .. External Subroutines ..
      EXTERNAL INTRP,EVALAGS,EVALAGD,DERIV

C     .. Intrinsic Functions ..
      INTRINSIC DMAX1,DABS

C ... Error check
      IF ((TS .LT. T) .or. (ts2 .lt. t)) THEN
         PRINT *, '***** ERROR AT DEFECT *****'
         PRINT *, 'T=', T, ' TS=', TS, ' TS2 =', TS2
         PRINT *, 'TS (or TS2) is less than T !!!'
         PRINT *
         STOP
      END IF

      ND = MU + OMEGA
      H = C(18)

C ... Compute P(TS)
      CALL INTRP(N,NW,T,Y,TS,W(1,2),H,W)

C ... Compute argument of delay value
      CALL SIGMA(N,ND,TS,W(1,2),W(1,15))

C     Compute P(Delay Arugment)
      DO 110 I = 1, MU
         IW = 20 + I

         IF (W(I,15) .GT. TS) THEN
            PRINT *, '**** Warning:DEFECT *****'
            PRINT *, '    ADVANCED DELAY in DEFECT'
            PRINT *, 'T=',TS,'  Delay Argument=',W(I,15)
            PRINT *
         END IF

         SIGN = W(I,15) - T

C ... In case there is no vanishing delay
         IF (SIGN .LE. 0.D0) THEN
            CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +           TOLD,KOLD,W,YINIT)

C ... In case there is vanishing delay
         ELSE
            CALL INTRP(N,NW,T,Y,W(I,15),W(1,IW),H,W)
         END IF
 110  CONTINUE
       
C ...  Compute P'(Delay Argument)
      DO 115 I = MU + 1, ND

         IW = 20 + I

         IF (W(I,15) .GT. TS) THEN
            PRINT *, '**** Warning:DEFECT *****'
            PRINT *, '    ADVANCED DELAY in DEFECT'
            PRINT *, 'T=',TS,'  Delay Argument=',W(I,15)
            PRINT *
         END IF

         SIGN = W(I,15) - T

C ... In case there is no vanishing delay
         IF (SIGN .LE. 0.D0) THEN
            CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +           TOLD,KOLD,W,DYINIT)

C ... In case there is vanishing delay
         ELSE
            CALL DERIV(N,NW,T,Y,W(I,15),W(1,IW),H,W)
         END IF
 115  CONTINUE

C ... Compute F(TS, P(TS), P(DEL ARG), P'(DELAY ARG))
      CALL FCN(N,ND,NW,TS,W(1,2),W(1,21),W(1,18))
      C(24) = C(24) + 1.D0

C ... Compute P'(TS) and then defect
      CALL DERIV(N,NW,T,Y,TS,W(1,2),H,W)

      DO 120 I = 1, N
         W(I,14) = W(I,18) - W(I,2)
 120  CONTINUE

C ... Calculate the weighted max norm of W(*,14) as specified by
C     the error control indicator C(1)
      TEMP = 0.D0
      IF (C(1) .NE. 1.D0) GO TO 210

C ... Absolute error control
      DO 200 K = 1, N
         TEMP = DMAX1(TEMP,DABS(W(K,14)))
 200  CONTINUE
      GO TO 500
 210  IF (C(1) .NE. 2.D0) GO TO 240

C ... Relative error control
      DO 230 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)/W(K,18)))
 230  CONTINUE
      GO TO 500
 240  IF (C(1) .NE. 3.D0) GO TO 260

C ... Weights are 1/MAX(C(2),ABS(W(K,18)))
      DO 250 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)) / DMAX1(C(2), 
     +        DABS(W(K,18))))
 250  CONTINUE
      GO TO 500
 260  IF (C(1) .NE. 4.D0) GO TO 280

C ... Weights are 1/MAX(C(K+30),ABS(W(K,18)))
      DO 270 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)) / DMAX1(C(K+30), 
     +        DABS(W(K,18))))
 270  CONTINUE
      GO TO 500
 280  IF (C(1) .NE. 5.D0) GO TO 300

C ... Weights are 1/C(K+30)
      DO 290 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)/C(K+30)))
 290  CONTINUE
      GO TO 500
 300  CONTINUE

C ... Default case - Weights are 1/MAX(1,ABS(Y'(K)))
      DO 310 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)) / DMAX1(1.D0, 
     +        DABS(W(K,18))) )
 310  CONTINUE
 500  CONTINUE

      DFCT = TEMP


C ... Compute P(TS2)
      CALL INTRP(N,NW,T,Y,TS2,W(1,2),H,W)

C ... Compute argument of delay value
      CALL SIGMA(N,ND,TS2,W(1,2),W(1,15))

C     Compute P(Delay Arugment)
      DO 610 I = 1, MU
         IW = 20 + I

         IF (W(I,15) .GT. TS2) THEN
            PRINT *, '**** Warning:DEFECT *****'
            PRINT *, '    ADVANCED DELAY in DEFECT'
            PRINT *, 'T=',TS2,'  Delay Argument=',W(I,15)
            PRINT *
         END IF

         SIGN = W(I,15) - T

C ... In case there is no vanishing delay
         IF (SIGN .LE. 0.D0) THEN
            CALL EVALAGS(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +           TOLD,KOLD,W,YINIT)

C ... In case there is vanishing delay
         ELSE
            CALL INTRP(N,NW,T,Y,W(I,15),W(1,IW),H,W)
         END IF
 610  CONTINUE
       
C ...  Compute P'(Delay Argument)
      DO 615 I = MU + 1, ND

         IW = 20 + I

         IF (W(I,15) .GT. TS2) THEN
            PRINT *, '**** Warning:DEFECT *****'
            PRINT *, '    ADVANCED DELAY in DEFECT'
            PRINT *, 'T=',TS2,'  Delay Argument=',W(I,15)
            PRINT *
         END IF

         SIGN = W(I,15) - T

C ... In case there is no vanishing delay
         IF (SIGN .LE. 0.D0) THEN
            CALL EVALAGD(N,NW,IW,LOLD,C(25),C(26),C(27),W(I,15),
     +           TOLD,KOLD,W,DYINIT)

C ... In case there is vanishing delay
         ELSE
            CALL DERIV(N,NW,T,Y,W(I,15),W(1,IW),H,W)
         END IF
 615  CONTINUE

C ... Compute F(TS2, P(TS2), P(DEL ARG), P'(DELAY ARG))
      CALL FCN(N,ND,NW,TS2,W(1,2),W(1,21),W(1,18))
      C(24) = C(24) + 1.D0

C ... Compute P'(TS2) and then defect
      CALL DERIV(N,NW,T,Y,TS2,W(1,2),H,W)

      DO 620 I = 1, N
         W(I,14) = W(I,18) - W(I,2)
 620  CONTINUE

C ... Calculate the weighted max norm of W(*,14) as specified by
C     the error control indicator C(1)
      TEMP = 0.D0
      IF (C(1) .NE. 1.D0) GO TO 710

C ... Absolute error control
      DO 700 K = 1, N
         TEMP = DMAX1(TEMP,DABS(W(K,14)))
 700  CONTINUE
      GO TO 900
 710  IF (C(1) .NE. 2.D0) GO TO 740

C ... Relative error control
      DO 730 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)/W(K,18)))
 730  CONTINUE
      GO TO 900
 740  IF (C(1) .NE. 3.D0) GO TO 760

C ... Weights are 1/MAX(C(2),ABS(W(K,18)))
      DO 750 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)) / DMAX1(C(2), 
     +        DABS(W(K,18))))
 750  CONTINUE
      GO TO 900
 760  IF (C(1) .NE. 4.D0) GO TO 780

C ... Weights are 1/MAX(C(K+30),ABS(W(K,18)))
      DO 770 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)) / DMAX1(C(K+30), 
     +        DABS(W(K,18))))
 770  CONTINUE
      GO TO 900
 780  IF (C(1) .NE. 5.D0) GO TO 800

C ... Weights are 1/C(K+30)
      DO 790 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)/C(K+30)))
 790  CONTINUE
      GO TO 900
 800  CONTINUE

C ... Default case - Weights are 1/MAX(1,ABS(Y'(K)))
      DO 810 K = 1, N
         TEMP = DMAX1(TEMP, DABS(W(K,14)) / DMAX1(1.D0, 
     +        DABS(W(K,18))) )
 810  CONTINUE
 900  CONTINUE

c      print *, 'ts =', ts, ' ts2 =', ts2
c      print *, 'dfct1 =', dfct, ', defct2 =', temp

      dfct = dmax1(dfct,temp)

      RETURN
      END


C *** End of DEFECT2 ***




      SUBROUTINE DEFECTG(N,NW,MU,OMEGA,LOLD,M,T1,TS,H,TOLD,KOLD,W,C,
     +                   FCN,SIGMA,YINIT,DYINIT)


C**********************************************************************
C                                                                     *
C     Purpose - This subroutine computes defect at XS and put them in *
C               W(*,14). The DEFECT computes the defect only on the   *
C               current step using INTRP and DERIV while DEFECTG      *
C               computes the defect in any step using INTRPG and      *
C               DERIVG.                                               *
C                                                                     *
C**********************************************************************
C                                                                     *
C    N      Number of equation                                        *
C    NW     First dimension of workspace W & queue KOLD               *
C    MU     Number of solution delay terms                            *
C    OMEGA  Number of derivative delay terms                          *
C    LOLD   Length of queue                                           *
C    M      Index for the queue, i.e, TOLD(M-1) <= TS <= TOLD(M)      *
C    T1     Independent value (left end point)                        *  
C    TS     Independent value at which defect is to be computed       *
C    H      Stepsize, i.e, TOLD(M) - TOLD(M-1)                        *
C    W      Work space                                                *
C    C      Communication array                                       *
C    TOLD   Queue of independent values                               *
C    KOLD   Queue of k-values and solutions                           *
C    FCN    Name of subroutine for computing derivative               *
C    SIGMA  Name of subroutine for computing delay argument           *
C    YINIT  Name of subroutine for computing initial values at        *
C           initial intervals                                         *
C    DYINIT Name of subroutine for computing derivative of initial    *
C           values at initial intervals                               *
C                                                                     *
C**********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,LOLD,M,MU,OMEGA
      DOUBLE PRECISION T1,TS,H

C     .. Array Arguments ..
      DOUBLE PRECISION TOLD(LOLD),KOLD(NW,7,LOLD),W(NW,*),C(*)

C     .. Subroutine Arguments ..
      EXTERNAL FCN,SIGMA,YINIT,DYINIT

C     .. Local Scalars ..
      INTEGER I,ND

C     .. External Subroutines ..
      EXTERNAL INTRPG,EVALAGS,EVALAGD,DERIVG

C ... Error check      
      IF (TS .LT. T1) THEN
         PRINT *, '***** ERROR AT DEFECT *****'
         PRINT *, 'X=', T1, ' TS=', TS
         PRINT *, 'TS is less than X !!!'
         PRINT *
         STOP
      END IF

      ND = MU + OMEGA

C ... Compute P(TS)
      CALL INTRPG(N,NW,M,T1,TS,W(1,2),H,KOLD)

C ... Compute argument of delay value
      CALL SIGMA(N,ND,TS,W(1,2),W(1,15))

C ... Compute P(Delay Arguemnt)
      DO 20 I = 1, MU

         IF (W(I,15) .GT. TS) THEN
            PRINT *, '**** WARNING AT DEFECTG *****'
            PRINT *, '    ADVANCED DELAY in DEFECTG'
            PRINT *, 'X=',TS,' Delay Argument =', W(I,15)
            PRINT *
         END IF

         CALL EVALAGS(N,NW,20+I,LOLD,C(25),C(26),C(27),W(I,15),
     +        TOLD,KOLD,W,YINIT)

 20   CONTINUE

C ... Compute P'(Delay Arugment)
      DO 25 I = MU + 1, ND

         IF (W(I,15) .GT. TS) THEN
            PRINT *, '**** WARNING AT DEFECTG *****'
            PRINT *, '    ADVANCED DELAY in DEFECTG'
            PRINT *, 'X=',TS,' Delay Argument =', W(I,15)
            PRINT *
         END IF

         CALL EVALAGD(N,NW,20+I,LOLD,C(25),C(26),C(27),W(I,15),
     +        TOLD,KOLD,W,DYINIT)

 25   CONTINUE

C ... Compute F(TS, P(TS), P(DEL ARG), P'(DELAY ARG))
      CALL FCN(N,ND,NW,TS,W(1,2),W(1,21),W(1,18))
      C(24) = C(24) + 1.D0

C ... Compute P'(TS) and then defect
      CALL DERIVG(N,NW,M,T1,TS,W(1,2),H,KOLD)

      DO 30 I = 1, N
         W(I,14) = W(I,18) - W(I,2)
 30   CONTINUE

      RETURN
      END


C *** End of DEFECTG ***




      SUBROUTINE DERIV(N,NW,T,Y,TOUT,YPOUT,HTRIAL,W)     

C**********************************************************************
C                                                                     *
C   Purpose - This subroutine computes derivatives at TOUT            *
C                                                                     *
C**********************************************************************
C                                                                     *
C     N      Number of equations                                      *
C     NW     First dimension of workspace W & queue KOLD              *
C     T      Current independent value                                *
C     Y      Current solution values, i.e., at T                      *
C     TOUT   Independent value at which derivative is computed        *
C     YPOUT  Derivative value at TOUT (OUTPUT)                        *
C     HTRIAL Stepsize                                                 *
C     W      Workspace                                                *
C                                                                     *
C**********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW
      DOUBLE PRECISION T,TOUT,HTRIAL

C     .. Array Arguments ..
      DOUBLE PRECISION Y(N),YPOUT(N),W(NW,12)

C     .. Local Scalars ..
      INTEGER I
      DOUBLE PRECISION TAU,TAU2,TAU3,TAU4,TAU5
      double precision DD0,DD1,DD2,DD3,DD9,DD10,DD11


      TAU = (TOUT - T)/HTRIAL
      TAU2 = TAU * TAU
      TAU3 = TAU2 * TAU
      TAU4 = TAU3 * TAU
      TAU5 = TAU4 * TAU


      DD0 = 4.285714285716394d0 * TAU - 57.85714285714626d0 * TAU2 +
     +     248.5714285713638d0 * TAU3 - 404.9999999997968d0 * TAU4  +
     +     209.9999999998538d0 * TAU5
      DD1 = 1.d0 - 14.53571428571327d0 * TAU + 71.98214285713236d0 *
     +     TAU2 - 154.5714285713923d0 * TAU3  + 146.8749999999520d0 * 
     +     TAU4 - 50.74999999997890d0 * TAU5 
      DD2 = - 4.285714285716394d0 * TAU + 57.85714285714626d0 * TAU2 -
     +     248.5714285713638d0 * TAU3  + 404.9999999997968d0 * TAU4 -
     +     209.9999999998538d0 * TAU5
      DD3 = 0.5309523809522572D0 * TAU - 7.234523809521525D0 
     +  * TAU2  + 31.62857142856110D0 * TAU3 - 53.20833333331760D0 
     +  * TAU4  + 29.28333333332581D0 * TAU5
      DD9 =  6.857142857142692D0 * TAU - 88.30476190475958D0 
     +  * TAU2  + 348.6476190476093D0 * TAU3 - 498.6666666666516D0
     +  * TAU4  + 231.4666666666593D0 * TAU5
      DD10 =  28.58333333332976D0 * TAU - 225.8083333333001D0 
     +  * TAU2  + 617.3999999998936D0 * TAU3 - 700.2916666665310D0 
     +  * TAU4  + 280.1166666666077D0 * TAU5
      DD11 = - 17.14999999999986D0 * TAU + 191.5083333333260D0
     +  * TAU2  - 594.5333333332940D0 * TAU3 + 700.2916666666015D0 
     +  * TAU4  - 280.1166666666346D0 * TAU5

      DO 110 I = 1, N
         YPOUT(I) = (DD0 * y(i) + DD2 * w(i,12)) / HTRIAL +
     +        DD1 * w(i,1) + DD3 * w(i,8) + DD9 * 
     +        w(i,9) + DD10 * w(i,10) + DD11 * w(i,11)
 110  CONTINUE

      RETURN
      END


C *** End of DERIV ***






      SUBROUTINE DERIVG(N,NW,M,T,TOUT,YPOUT,HTRIAL,KOLD)     

C**********************************************************************
C                                                                     *
C   Purpose - This subroutine computes derivatives at TOUT            *
C             The DERIV computes the derivative only on the current   *
C             step while DERIVG computes the derivative in any step   *
C                                                                     *
C**********************************************************************
C                                                                     *
C     N      Number of equations                                      *
C     NW     First dimension of workspace W & queue KOLD              *
C     M      Index for the queue, i.e. TOLD(M-1) <= TOUT <= TOLD(M)   *
C     T      Current independent value                                *
C     TOUT   Independent value at which derivative is computed        *
C     YPOUT  Derivative value at TOUT (OUTPUT)                        *
C     HTRIAL Stepsize                                                 *
C     KOLD   Queue of K values and solution values                    *
C                                                                     *
C**********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,M
      DOUBLE PRECISION T,TOUT,HTRIAL

C     .. Array Arguments ..
      DOUBLE PRECISION YPOUT(N),KOLD(NW,7,*)

C     .. Local Scalars ..
      INTEGER I
      DOUBLE PRECISION TAU,TAU2,TAU3,TAU4,TAU5
      DOUBLE PRECISION DD0,DD1,DD2,DD3,DD9,DD10,DD11

      TAU = (TOUT - T)/HTRIAL
      TAU2 = TAU * TAU
      TAU3 = TAU2 * TAU
      TAU4 = TAU3 * TAU
      TAU5 = TAU4 * TAU

      DD0 = 4.285714285716394d0 * TAU - 57.85714285714626d0 * TAU2 +
     +     248.5714285713638d0 * TAU3 - 404.9999999997968d0 * TAU4  +
     +     209.9999999998538d0 * TAU5
      DD1 = 1.d0 - 14.53571428571327d0 * TAU + 71.98214285713236d0 *
     +     TAU2 - 154.5714285713923d0 * TAU3  + 146.8749999999520d0 * 
     +     TAU4 - 50.74999999997890d0 * TAU5 
      DD2 = - 4.285714285716394d0 * TAU + 57.85714285714626d0 * TAU2 -
     +     248.5714285713638d0 * TAU3  + 404.9999999997968d0 * TAU4 -
     +     209.9999999998538d0 * TAU5
      DD3 = 0.5309523809522572D0 * TAU - 7.234523809521525D0 
     +  * TAU2  + 31.62857142856110D0 * TAU3 - 53.20833333331760D0 
     +  * TAU4  + 29.28333333332581D0 * TAU5
      DD9 =  6.857142857142692D0 * TAU - 88.30476190475958D0 
     +  * TAU2  + 348.6476190476093D0 * TAU3 - 498.6666666666516D0
     +  * TAU4  + 231.4666666666593D0 * TAU5
      DD10 =  28.58333333332976D0 * TAU - 225.8083333333001D0 
     +  * TAU2  + 617.3999999998936D0 * TAU3 - 700.2916666665310D0 
     +  * TAU4  + 280.1166666666077D0 * TAU5
      DD11 = - 17.14999999999986D0 * TAU + 191.5083333333260D0
     +  * TAU2  - 594.5333333332940D0 * TAU3 + 700.2916666666015D0 
     +  * TAU4  - 280.1166666666346D0 * TAU5

      DO 110 I = 1, N
         YPOUT(I) = (DD0 * KOLD(I,1,M) + DD2 * KOLD(I,3,M)) / HTRIAL +
     +        DD1 * KOLD(I,2,M) + DD3 * KOLD(I,4,M) + DD9 * 
     +        KOLD(I,5,M) + DD10 * KOLD(I,6,M) + DD11 * KOLD(I,7,M)
 110  CONTINUE

      RETURN
      END


C *** End of DERIVG ***





      SUBROUTINE DSLOC(N,NW,MU,OMEGA,IND2,LOLD,T,TS,Y,YS,TOLD,KOLD,
     +                 TOL,W,C,FCN,SIGMA,YINIT,DYINIT)

C********************************************************************
C                                                                   *
C     Purpose - This routine searches discontinuity point           *
c               using the techiniques proposed by                   *
c               enright et al                                       *
C                                                                   *
C********************************************************************
C                                                                   *
C    N      Number of equation                                      *
C    NW     First dimension of workspace W & queue KOLD             *
C    MU     Number of solution delay terms                          *
C    OMEGA  Number of derivative delay terms                        *
C    IND2   Flag for the discontinuity                              *
C           ... if found IND2 > 0                                   *
C    LOLD   Length of queue                                         *
C    T      Independent value (Tn-1)                                * 
C    TS     Discontinous point to be found (output)                 *
C    Y      Solution at T                                           *
C    YS     Solution at TS (Output)                                 *
C    TOLD   Queue of independent values                             *
C    KOLD   Queue of K-values and solutions                         *
C    TOL    Tolerance                                               *
C    W      Workspace                                               *
C    C      Communication array                                     *
C    FCN    Name of subroutine for computing derivative             *
C    SIGMA  Name of subroutine for computing delay argument         *
C    YINIT  Name of subroutine for computing initial values at      *
C           the initial interval                                    *
C    DYINIT Name of subroutine for computing initial values at      *
C           the initial interval                                    *
C                                                                   *
C********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,IND2,LOLD,MU,OMEGA
      DOUBLE PRECISION T,TS,TOL

C     .. Array Arguments ..
      DOUBLE PRECISION Y(N),YS(N),TOLD(LOLD),KOLD(NW,7,LOLD)
      DOUBLE PRECISION W(NW,*),C(*)

C     .. Subroutine Arguments ..
      EXTERNAL FCN,SIGMA,YINIT,DYINIT

C     .. Local Scalars ..
      DOUBLE PRECISION K2,K21,K22,DFCT,EPS,DELTA,TL,TH

C     .. External Subroutines ..
      EXTERNAL DEFECT,INTRP

C     .. Intrisic Functions ..
      INTRINSIC DMAX1,DMIN1

C     .. Data statements ..      
C        SET THE VALUES FOR DECIDING IF WE ARE TO THE RIGHT OF Td.
      DATA K21,K22/400.D0, 50.D0/

C ... Set Tl to be Tn-1 + H and Th to be Tn-1 + 2H
      TL = C(17)
      TH = TL + C(18)

C ... To check if there is discontinuity, first we have to check 
C     The scale of defect at TH. Thus, seT K2 = 400.
      K2 = K21
      DFCT = C(19)

C ... Set DELTA = 400 * DEFECT at Tstar, WHERE Tn-1 < Tstar < Tn-1 + H
C     DELTA SHOULD BE 10*TOL < DELTA < 1
      DELTA = DMAX1(10.D0*TOL, K2*DFCT)
      DELTA = DMIN1(1.D0, DELTA)

      TS = TH
      CALL DEFECT(N,NW,MU,OMEGA,LOLD,T,TS,Y,DFCT,TOLD,KOLD,W,C,
     +            FCN,SIGMA,YINIT,DYINIT)

C ... Set EPS to be TOL / DEFECT at Tn-1 + 2H for stopping criterion
c      EPS = DMIN1(C(18), TOL / DFCT)
      EPS = DMIN1(C(18)/20.d0, TOL / DFCT)

C ... If DEFECT at Tn-1 + 2H is less than 400 * DEFECT at Tstar,
C     where Tn-1 < Tstar < Tn-1 + H, there is no discontinuity.
      IND2 = -1
      IF (DFCT .LT. DELTA) RETURN

C ... Check the DEFECT at Tn-1 + 3/2 H and then decide which K2 value
C     should be used.
      TS = TL + 0.5D0 * C(18)
      CALL DEFECT(N,NW,MU,OMEGA,LOLD,T,TS,Y,DFCT,TOLD,KOLD,W,C,
     +            FCN,SIGMA,YINIT,DYINIT)

C ... If DEFECT at Tn-1 + 3/2 H > 400 * DEFECT at Tstar,
C     then Tn-1 + H < Td < Tn-1 + 3/2 H
      IF (DFCT .GT. DELTA) THEN
         TH = TS
         EPS = DMIN1(C(18)/20.d0, TOL / DFCT)

C     If DEFECT at Tn-1 + 3/2 H <= 400 * DEFECT at Tstar,
C     then Tn-1 + 3/2 H < Td < Tn-1 + 2H
      ELSE
         TL = TS
         K2 = K22
         DELTA = DMAX1(10.D0*TOL, K2*DFCT)
         DELTA = DMIN1(1.D0, DELTA)
      END IF

      IND2 = 2

C ... Loop --- Bisection search
 1000 CONTINUE

         IF ((TH - TL) .LT. EPS) THEN
            TS = TH
            CALL INTRP(N,NW,T,Y,TS,YS,C(18),W)
            C(30) = TL
            C(31) = TH - TL
            RETURN
         END IF

         TS = (TH + TL) / 2.D0
         CALL DEFECT(N,NW,MU,OMEGA,LOLD,T,TS,Y,DFCT,TOLD,KOLD,W,C,
     +               FCN,SIGMA,YINIT,DYINIT)

         IND2 = IND2 + 1

         IF (DFCT .LT. DELTA) THEN
            TL = TS
         ELSE
            TH = TS
            EPS = DMIN1(C(18)/20.d0, TOL / DFCT)
         END IF

         GOTO 1000

C     End loop --- Bisection search

      END


C *** End of DSLOC ***





      SUBROUTINE EVALAGS(N,NW,IW,LOLD,X0,IFIRST,ILAST,DELARG,TOLD,KOLD,
     +                  W,YINIT)

C********************************************************************
C                                                                   *
C     Purpose - This is the routine that computes delay values      *
c               by using interpolation routine.                     *
C                                                                   *
C********************************************************************
C                                                                   *
C     N       Number of equations                                   *
C     NW      First dimension of workspace W & queue KOLD           *
C     IW      Index for strage space                                *
C     LOLD    Length of the queue                                   *
C     X0      Initial independet value                              *
C     IFIRST  The first index of the queue                          *
C     ILAST   The last index of the queue                           *
C     DELARG  Delay argement                                        *
C     TOLD    Queue of the past independent variables               *
C     KOLD    Queue of the past K-values and soutions at mesh       *
C     W       Workspace matrix (W(*,21:20+I) is used to store delay *
C             values and W(*,21+I:20+ND+I) is used to store         *
C             derivative delay values)                              *
C     YINIT   Name of subroutine for computing initial values       *
C             at initial interval.                                  *
C                                                                   *
C********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,IW,LOLD
      DOUBLE PRECISION X0,IFIRST,ILAST,DELARG

C     .. Array Arguments ..
      DOUBLE PRECISION TOLD(LOLD),KOLD(NW,7,LOLD),W(NW,*)

C     .. Subroutine Arguments ..
      EXTERNAL YINIT


C     .. Local Scalars ..
      INTEGER J,M
      DOUBLE PRECISION H

C     .. External Subroutines ..
      EXTERNAL TFIND,INTRPG

C ... Check whthere the delay is before initial point
      IF (DELARG .LE. X0) THEN
         CALL YINIT(N,DELARG,W(1,IW))

C ... Otherwise find the index in the queue for the delay function
      ELSE 

         CALL TFIND(LOLD,TOLD,DELARG,IFIRST,ILAST,M)

C ... In case the data is on the mesh
         IF (TOLD(M) .EQ. DELARG) THEN
            DO 120 J = 1, N
               W(J,IW) = KOLD(J,3,M)
 120        CONTINUE

C ... Otherwise evaluate delay value using intrp routine
C ... In case the index is equal to the last index
         ELSE IF (M .EQ. LOLD) THEN
c         IF (M .EQ. LOLD) THEN
            H = TOLD(1) - TOLD(M)
            CALL INTRPG(N,NW,1,TOLD(M),DELARG,W(1,IW),H,KOLD)

C ... Otherwise
         ELSE
            H = TOLD(M+1) - TOLD(M)
            CALL INTRPG(N,NW,M+1,TOLD(M),DELARG,W(1,IW),H,KOLD)
         END IF 
      END IF

      RETURN
      END


C *** End of EVALAGS ***


      SUBROUTINE EVALAGD(N,NW,IW,LOLD,X0,IFIRST,ILAST,DELARG,TOLD,KOLD,
     +                  W,DYINIT)

C********************************************************************
C                                                                   *
C     Purpose - This is the routine that computes delay values      *
c               by using interpolation routine.                     *
C                                                                   *
C********************************************************************
C                                                                   *
C     N       Number of equations                                   *
C     NW      First dimension of workspace W & queue KOLD           *
C     IW      Index of storage space                                *
C     LOLD    Length of the queue                                   *
C     X0      Initial independet value                              *
C     IFIRST  The first index of the queue                          *
C     ILAST   The last index of the queue                           *
C     DELARG  Delay argement                                        *
C     TOLD    Queue of the past independent variables               *
C     KOLD    Queue of the past K-values and soutions at mesh       *
C     W       Workspace matrix (W(*,21:20+I) is used to store delay *
C             values and W(*,21+I:20+ND+I) is used to store         *
C             derivative delay values)                              *
C     DYINIT  Name of subroutine for computing derivative initial   *
C             values at initial interval.                           *
C                                                                   *
C********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,IW,LOLD
      DOUBLE PRECISION X0,IFIRST,ILAST,DELARG

C     .. Array Arguments ..
      DOUBLE PRECISION TOLD(LOLD),KOLD(NW,7,LOLD),W(NW,*)

C     .. Subroutine Arguments ..
      EXTERNAL DYINIT

C     .. Local Scalars ..
      INTEGER J,M
      DOUBLE PRECISION H

c     .. arrays in common      
c      double precision c(40)
c      common /comm/c

C     .. External Subroutines ..
      EXTERNAL TFIND,DERIVG

C ... Check whthere the delay is before initial point
      IF (DELARG .LE. X0) THEN
         CALL DYINIT(N,DELARG,W(1,IW))

c      if ((c(40) .eq. 1.d0) .and. (delarg .le. x0)) then
c         call dyinit(n,delarg,w(1,iw))
c      else if ((c(40) .eq. 2.d0) .and. (delarg .lt. x0)) then 
c         call dyinit(n,delarg,w(1,iw))

C ... Otherwise find the index in the queue for the delay function
      ELSE 

         CALL TFIND(LOLD,TOLD,DELARG,IFIRST,ILAST,M)


C ... In case the data is on the mesh
         IF (TOLD(M) .EQ. DELARG) THEN
            DO 120 J = 1, N
               W(J,IW) = KOLD(J,4,M)
 120        CONTINUE

C ... Otherwise evaluate delay value using intrp routine
c ... In case the index is equal to the last index
         ELSE IF (M .EQ. LOLD) THEN
c         IF (M .EQ. LOLD) THEN
            H = TOLD(1) - TOLD(M)
            CALL DERIVG(N,NW,1,TOLD(M),DELARG,W(1,IW),H,KOLD)

C ... Otherwise
         ELSE
            H = TOLD(M+1) - TOLD(M)
            CALL DERIVG(N,NW,M+1,TOLD(M),DELARG,W(1,IW),H,KOLD)
         END IF 
      END IF

      RETURN
      END


C *** End of EVALAGD ***




      SUBROUTINE INTRP(N,NW,T,Y,TOUT,YPOUT,HTRIAL,W)

C**********************************************************************
C                                                                     *
C   Purpose - This subroutine computes solutions at TOUT              *
C                                                                     *
C**********************************************************************
C                                                                     *
C     N      Number of equations                                      *
C     NW     First dimension of workspace W & queue KOLD              *
C     T      Current independent value                                *
C     Y      Current solution values, i.e., at T                      *
C     TOUT   Independent value at which solution is computed          *
C     YPOUT  Solution  at TOUT (Output)                               *
C     HTRIAL Stepsize                                                 *
C     W      Workspace                                                *
C                                                                     *
C**********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW
      DOUBLE PRECISION T,TOUT,HTRIAL

C     .. Array Arguments ..
      DOUBLE PRECISION Y(N),YPOUT(N),W(NW,12)

C     .. Local Scalars ..
      INTEGER I
      DOUBLE PRECISION TAU,TAU2,TAU3,TAU4,TAU5,TAU6
      double precision d0,d1,d2,d3,d9,d10,d11
c     double precision bb(11)

      TAU = (TOUT - T)/HTRIAL
      TAU2 = TAU * TAU
      TAU3 = TAU2 * TAU
      TAU4 = TAU3 * TAU
      TAU5 = TAU4 * TAU
      TAU6 = TAU5 * TAU


      D0 = 1.D0 + 2.142857142858197D0 * TAU2 - 19.28571428571542D0
     + * TAU3 + 62.14285714284096D0 * TAU4 - 80.99999999995936D0
     + * TAU5 + 34.99999999997563D0 * TAU6
      D1 = TAU - 7.267857142856635D0 * TAU2 + 23.99404761904412D0
     + * TAU3 - 38.64285714284807D0 * TAU4 + 29.37499999999040D0
     + * TAU5 - 8.458333333329817D0 * TAU6
      D2 =  - 2.142857142858197D0 * TAU2 + 19.28571428571542D0 * TAU3
     + - 62.14285714284096D0 * TAU4 + 80.99999999995936D0 * TAU5 - 
     + 34.99999999997563D0 * TAU6
      D3 = 0.2654761904761286D0 * TAU2 - 2.411507936507175D0 * TAU3
     + + 7.907142857140276D0 * TAU4 - 10.64166666666352D0 * TAU5
     + + 4.880555555554301D0 * TAU6
      D9 = 3.428571428571346D0 * TAU2  - 29.43492063491986D0 
     + * TAU3  + 87.16190476190233D0 * TAU4 - 99.73333333333032D0 
     + * TAU5  + 38.57777777777655D0 * TAU6 
      D10 = 14.29166666666488D0 * TAU2  - 75.26944444443335D0 
     + * TAU3  + 154.3499999999734D0 * TAU4 - 140.0583333333062D0 
     + * TAU5  + 46.68611111110129D0 * TAU6
      D11 = - 8.574999999999928D0 * TAU2  + 63.83611111110867D0 
     + * TAU3  - 148.6333333333235D0 * TAU4 + 140.0583333333203D0 
     + * TAU5  - 46.68611111110577D0 * TAU6

      DO 110 I = 1, N
         YPOUT(I) = D0 * y(i) + D2 * w(i,12) + HTRIAL * (D1
     +        * w(i,1) + D3 * w(i,8) + D9 * w(i,9) + D10
     +        * w(i,10) + D11 * w(i,11))
 110  CONTINUE


      RETURN
      END


C *** End of INTRP ***






      SUBROUTINE INTRPG(N,NW,M,T,TOUT,YPOUT,HTRIAL,KOLD)


C**********************************************************************
C                                                                     *
C   Purpose - This subroutine computes solutions at TOUT              *
C             The INTRP computes the solution only on the current     *
C             step while INTRPG computes the solution in any step     *
C                                                                     *
C**********************************************************************
C                                                                     *
C     N      Number of equations                                      *
C     NW     First dimension of workspace W & queue KOLD              *
C     M      Index for the queue, i.e. TOLD(M-1) <= TOUT <= TOLD(M)   *
C     T      Current independent value                                *
C     TOUT   Independent value at which solution is computed          *
C     YPOUT  Solution at TOUT (OUTPUT)                                *
C     HTRIAL Stepsize                                                 *
C     KOLD   Queue of K values and solution values                    *
C                                                                     *
C**********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,M
      DOUBLE PRECISION T,TOUT,HTRIAL

C     .. Array Arguments ..
      DOUBLE PRECISION YPOUT(N),KOLD(NW,7,*)

C     .. Local Scalars ..
      INTEGER I
      DOUBLE PRECISION TAU,TAU2,TAU3,TAU4,TAU5,TAU6
      DOUBLE PRECISION D0,D1,D2,D3,D9,D10,D11
C
      TAU = (TOUT - T)/HTRIAL
      TAU2 = TAU * TAU
      TAU3 = TAU2 * TAU
      TAU4 = TAU3 * TAU
      TAU5 = TAU4 * TAU
      TAU6 = TAU5 * TAU

      D0 = 1.D0 + 2.142857142858197D0 * TAU2 - 19.28571428571542D0
     + * TAU3 + 62.14285714284096D0 * TAU4 - 80.99999999995936D0
     + * TAU5 + 34.99999999997563D0 * TAU6
      D1 = TAU - 7.267857142856635D0 * TAU2 + 23.99404761904412D0
     + * TAU3 - 38.64285714284807D0 * TAU4 + 29.37499999999040D0
     + * TAU5 - 8.458333333329817D0 * TAU6
      D2 =  - 2.142857142858197D0 * TAU2 + 19.28571428571542D0 * TAU3
     + - 62.14285714284096D0 * TAU4 + 80.99999999995936D0 * TAU5 - 
     + 34.99999999997563D0 * TAU6
      D3 = 0.2654761904761286D0 * TAU2 - 2.411507936507175D0 * TAU3
     + + 7.907142857140276D0 * TAU4 - 10.64166666666352D0 * TAU5
     + + 4.880555555554301D0 * TAU6
      D9 = 3.428571428571346D0 * TAU2  - 29.43492063491986D0 
     + * TAU3  + 87.16190476190233D0 * TAU4 - 99.73333333333032D0 
     + * TAU5  + 38.57777777777655D0 * TAU6 
      D10 = 14.29166666666488D0 * TAU2  - 75.26944444443335D0 
     + * TAU3  + 154.3499999999734D0 * TAU4 - 140.0583333333062D0 
     + * TAU5  + 46.68611111110129D0 * TAU6
      D11 = - 8.574999999999928D0 * TAU2  + 63.83611111110867D0 
     + * TAU3  - 148.6333333333235D0 * TAU4 + 140.0583333333203D0 
     + * TAU5  - 46.68611111110577D0 * TAU6

      DO 110 I = 1, N
         YPOUT(I) = D0 * KOLD(I,1,M) + D2 * KOLD(I,3,M) + HTRIAL * (D1
     +        * KOLD(I,2,M) + D3 * KOLD(I,4,M) + D9 * KOLD(I,5,M) + D10
     +        * KOLD(I,6,M) + D11 * KOLD(I,7,M))
 110  CONTINUE

      RETURN
      END




      SUBROUTINE RDMETH(N,NW,MU,OMEGA,IND,LOLD,T,TEND,Y,TOLD,KOLD,TOL,
     +                  W,C,FCN,SIGMA,YINIT,DYINIT)

C**********************************************************************
C                                                                     *
C     Purpose - This routine sets up necessary data such as           *
C               stepsize, minimum stepsize, etc. Then call APFORM to  *
C               integrates one step by RK-methods for delay DEs.      *
C               This subroutine is based on DVERK.                    *
C                                                                     *
C**********************************************************************
C                                                                     *
C     N      Number of equations                                      *
C     NW     First dimension of workspace W & queue KOLD              *
C     MU     Number of solution delay terms                           *
C     OMEGA  Number of derivative delay terms                         *
C     IND    Indicator --- same as DVERK except IND = 7:              *
C            IND is set to be 7 after a discontinuity is found.       *
C     LOLD   The length of queue                                      *
C     T      Current independent value                                *
C     Y      Current solution, i.e., at T                             *
C     TOLD   Queue of T values                                        *
C     KOLD   Queue oF K values and solution values                    *
C     TOL    Tolerance                                                *
C     W      Workspace                                                *
C     C      Communication vector                                     *
C     FCN    Name of derivative routine                               *
C     SIGMA  Name of subroutine for computing delay argument          *
C     YINIT  Name of subroutine for computing initial values at       *
C            initial intervals.                                       *
C     DYINIT Name of subroutine for computing derivative of initial   *
C            values at initial intervals.                             *
C                                                                     *
C**********************************************************************


C     .. Scalar Arguments ..
      INTEGER N,NW,IND,LOLD,MU,OMEGA
      DOUBLE PRECISION T,TEND,TOL

C     .. Array Arguments ..
      DOUBLE PRECISION Y(N),TOLD(LOLD),KOLD(NW,7,LOLD),W(NW,*),C(*)

C     .. Subroutine Arguments ..
      EXTERNAL FCN,SIGMA,YINIT,DYINIT

C     .. Arrays in Common ..
      DOUBLE PRECISION AA(11,11),BB(11),CC(11)

C     .. Local Scalars ..
      INTEGER I,J,K,P,SSTAR,ND
      DOUBLE PRECISION H,HSCALE,TAUS,TAUS2,TEMP

C     .. External Subroutines ..
      EXTERNAL APFORM,DEFECT,EVALAGS,EVALAGD

C     .. Intrisic Functions ..
      INTRINSIC DABS,DFLOAT,DMIN1,DMAX1,DSIGN

C     .. Common blocks ..
      COMMON/METH/AA,BB,CC

C     ******************************************************************
C     * BEGIN INITIALIZATION, PARAMETER CHECKING, INTERRUPT RE-ENTRIES *
C     ******************************************************************

C ... Abort if IND out of range 1 to 7
      IF (IND.LT.1 .OR. IND.GT.8) GO TO 3000

      SSTAR = 11
      P = 6
      TAUS = .64D0
      TAUS2 = 0.85d0
      if (C(32) .eq. 3) TAUS = .859d0

      ND = MU + OMEGA

C ... Cases - Initial entry, normal re-entry, interrupt re-entries
      GO TO (100, 100, 200, 1111, 2222, 2222, 200, 2500), IND

C ... Case 1 - Initial entry (IND .EQ. 1 OR 2)
C ... Abort if N  > NW or TOL <= 0
 100  IF (N.GT.NW .OR. TOL.LE.0.D0) GO TO 3000

C ... Set up coefficients of RK-formula
      CC(1) = 0.D0
      CC(2) = 0.1D0
      CC(3) = 0.16D0
      CC(4) = 0.2040816326530613D0
      CC(5) = 0.6111D0
      CC(6) = 0.7778D0
      CC(7) = 0.1D1
      CC(8) = 0.1D1
      CC(9) = 0.5D0
      CC(10) = 0.1428571428571428D0
      CC(11) = 0.2857142857142857D0

      AA(2,1) = 0.1D0
      AA(3,1) = 0.32D-1
      AA(3,2) = 0.128D0
      AA(4,1) = 0.388656087174562D-1
      AA(4,2) = 0.9349845727545482D-1
      AA(4,3) = 0.7171756666015026D-1
      AA(5,1) = 0.9182106096978206D0
      AA(5,2) = -0.2806349313292671D0
      AA(5,3) = -0.4995012311203205D1
      AA(5,4) = 0.4968536632834652D1
      AA(6,1) = -0.1172309621097296D1
      AA(6,2) = 0.2577443286829264D0
      AA(6,3) = 0.7019866335902897D1
      AA(6,4) = -0.5919067674431333D1
      AA(6,5) = 0.5915666309428068D0
      AA(7,1) = 0.2052791398807182D1
      AA(7,2) = -0.2063792647412788D0
      AA(7,3) = -0.1319762903121892D2
      AA(7,4) = 0.1236821894967594D2
      AA(7,5) = -0.7279931113442874D0
      AA(7,6) = 0.7109910588213573D0
      AA(8,1) = 0.9752588552523605D-1
      AA(8,2) = 0.D0
      AA(8,3) = -0.4526249411455864D0
      AA(8,4) = 0.8022821731202383D0
      AA(8,5) = 0.2301388037778522D0
      AA(8,6) = 0.24584593549251D0
      AA(8,7) = 0.7683214322974987D-1
      AA(9,1) = 0.1089434604424868D0
      AA(9,2) = 0.D0
      AA(9,3) = -0.5148608705527948D0
      AA(9,4) = 0.8462848943497683D0
      AA(9,5) = 0.6712861232697209D-1
      AA(9,6) = -0.5132033903415767D-2
      AA(9,7) = -0.3361406266302458D-1
      AA(9,8) = 0.312500000000111D-1
      AA(10,1) = 0.8592255366477703D-1
      AA(10,2) = 0.D0
      AA(10,3) = -0.1524656894033379D0
      AA(10,4) = 0.2543389082003915D0
      AA(10,5) = 0.3082416724957052D-1
      AA(10,6) = 0.1449348604485854D-1
      AA(10,7) = -0.3149720157394562D-2
      AA(10,8) = -0.1427976438387382D-2
      AA(10,9) = -0.8567858630333475D-1
      AA(11,1) = 0.4859699015987465D-1
      AA(11,2) = 0.D0
      AA(11,3) = -0.573655693811402D-1
      AA(11,4) = 0.9353712191970651D-1
      AA(11,5) = 0.5261500935769923D-2
      AA(11,6) = -0.3816683388062108D-2
      AA(11,7) = -0.5124089070117709D-2
      AA(11,8) = 0.5665928935000797D-2
      AA(11,9) = 0.4744196280052546D-2
      AA(11,10) = 0.1942148893232018D0

      DO 110 J = 1, 7
         BB(J) = AA(8,J)
 110  CONTINUE
      DO 120 J = 8, 11
         BB(J) = 0.D0
 120  CONTINUE

      IF (IND .EQ. 2) GO TO 130

C ... Initial entry without options (IND .EQ. 1)
C ... Set C(1) to C(9) equal to 0
      DO 140 K = 1, 9
         C(K) = 0.D0
 140  CONTINUE
      GO TO 150

 130  CONTINUE
C ... Initial entry with options (IND .EQ. 2)
C ... Make C(1) to C(9) non-negative
      DO 160 K = 1, 9
         C(K) = DABS(C(K))
 160   CONTINUE

C ... Make floor values non-negative if they are to be used
      IF (C(1).NE.4.D0 .AND. C(1).NE.5.D0) GO TO 180
      DO 170 K = 1, N
         C(K+30) = DABS(C(K+30))
 170  CONTINUE
 180  CONTINUE
 150  CONTINUE

C ... Initialize rreb, dwarf, prev tend, flag, counts
      C(10) = 16.D0**(-13)
      C(11) = 1.D-28

C ... Set previous TEND initially to initial value of T
      C(20) = T

      DO 190 K = 21, 24
         C(K) = 0.D0
 190  CONTINUE
      GO TO 300

C ... Case 2 - Normal re-entry (IND .EQ. 3) or restart after
C              discontinuity (IND .EQ. 7)
C ....Abort if TEND reached, and either T changed or TEND not
 200  IF (C(21) .NE. 0.D0 .AND.
     +     (T. NE. C(20) .OR. TEND .EQ. C(20))) GO TO 3000

C ... Re-initialize flag
      C(21) = 0.D0
      GO TO 300

C ... Case 3 - Re-entry following an interrupt (IND .EQ. 4 TO 6)
C              transfer control to the appropriate re-entry point
C              This has already been handled by the computed go to
C              end case
 300   CONTINUE

C ... End initialization, etc.

C     ******************************************************************
C     * LOOP THROUGH THE FOLLOWING 4 STAGES, ONCE FOR EACH TRIAL  STEP *
C     * UNTIL THE OCCURRENCE OF ONE OF THE FOLLOWING		       *
C     *    (A) THE NORMAL RETURN (WITH IND .EQ. 3) ON REACHING TEND IN *
C     *        STAGE 4						       *
C     *    (B) AN ERROR RETURN (WITH IND .LT. 0) IN STAGE 1 OR STAGE 4 *
C     *    (C) AN INTERRUPT RETURN (WITH IND  .EQ.  4,	5  OR  6),  IF *
C     *        REQUESTED, IN STAGE 1 OR STAGE 4 		       *
C     ******************************************************************
C
99999 CONTINUE

C
C	 ***************************************************************
C	 * STAGE 1 - PREPARE - DO CALCULATIONS OF  HMIN,  HMAX,  ETC., *
C	 * AND SOME PARAMETER  CHECKING,  AND  END  UP	WITH  SUITABLE *
C	 * VALUES OF HMAG, TTRIAL AND HTRIAL IN PREPARATION FOR TAKING *
C	 * AN INTEGRATION STEP. 				       *
C	 ***************************************************************
C
C***********ERROR RETURN (WITH IND=-1) IF NO OF FCN EVALS TOO GREAT
      IF (C(7).EQ.0.D0 .OR. C(24).LT.C(7)) GO TO 400
      IND = -1
      RETURN
 400  CONTINUE

C ... Calculate slope (adding 1 to no of fcn evals) if IND .NE. 6
      IF (IND .EQ. 6) GO TO 410

C ... If C(28) is equal to 1.D0, then K1(in this step) is equal to K8
C     (in the previous step). Thus, no need to compute K1 here.


      IF (C(28) .EQ. 1.D0) THEN
         DO 420 I = 1, N
            W(I,1) = W(I,8)
 420     CONTINUE
      ELSE

         CALL SIGMA(N,ND,T,Y,W(1,15))

         DO 430 I = 1, MU
            IF (W(I,15) .GT. T) THEN
               PRINT *, '***** Warning: RDMETH *****'
               PRINT *, ' ADVANCED DELAY AT STAGE  1'
               PRINT *, 'T= ', T, ' Delay Argument =', W(I,15)

            ELSE 
               CALL EVALAGS(N,NW,20+I,LOLD,C(25),C(26),C(27),W(I,15),
     +                     TOLD,KOLD,W,YINIT)
            END IF

 430     CONTINUE

         DO 435 I = MU + 1, ND
            IF (W(I,15) .GT. T) THEN
               PRINT *, '***** Warning: RDMETH *****'
               PRINT *, ' ADVANCED DELAY AT STAGE  1'
               PRINT *, 'T= ', T, ' Delay Argument =', W(I,15)

            ELSE 
               
               CALL EVALAGD(N,NW,20+I,LOLD,C(25),C(26),C(27),W(I,15),
     +                     TOLD,KOLD,W,DYINIT)
               
            END IF

 435     CONTINUE

         CALL FCN(N,ND,NW,T,Y,W(1,21),W(1,1))
         C(24) = C(24) + 1.D0
      END IF

 410  CONTINUE

C ... Calculate HMIN - Use default unless value prescribed
      C(13) = C(3)
      IF (C(3) .NE. 0.D0) GO TO 440
      C(13) = 1.D-8
c      C(13) = 1.D-6
 440  CONTINUE

C ... Calculate SCALE - Use default unless value prescribed
      C(15) = C(5)
      IF (C(5) .EQ. 0.D0) C(15) = 1.D0

C ... Calculate HMAX - Consider 4 cases
C ... Case 1 both HMAX and SCALE prescribed
      IF (C(6).NE.0.D0 .AND. C(5).NE.0.D0)
     +     C(16) = DMIN1(C(6), 2.D0/C(5))

C ... Case 2 - HMAX prescribed, but SCALE not
      IF (C(6).NE.0.D0 .AND. C(5).EQ.0.D0) C(16) = C(6)

C ... Case 3 - HMAX not prescribed, but SCALE is
      IF (C(6).EQ.0.D0 .AND. C(5).NE.0.D0) C(16) = 2.D0/C(5)

C ... Case 4 - Neither HMAX nor SCALE is provided
      IF (C(6).EQ.0.D0 .AND. C(5).EQ.0.D0) C(16) = 1.5D0

C ... Error return (with IND=-2) if HMIN > HMAX
      IF (C(13) .LE. C(16)) GO TO 450
      IND = -2
      RETURN
 450  CONTINUE

C ... Calculate preliminary HMAG - Consider 3 cases
      IF (IND .GT. 2 .AND. IND .LT. 7) GO TO 460

C ... Case 1 - Initial entry - use prescribed value of HSTART, if
C              any, else default
      C(14) = C(4)

C ... If the delay vanishes at the initial time, then set the stepsize
C     small.
      if (omega .ge. 1) then
         CALL SIGMA(N,ND,T,Y,W(1,15))
         DO 470 I = 1, ND
            IF (((W(I,15) - T) .GT. -C(10))) THEN
               IF ( C(4) .EQ.  0.D0)  C(14) = 1.D2*C(13) 
            END IF
 470     CONTINUE
      end if

      IF (C(14) .EQ. 0.D0) C(14) = C(16)*TOL**(1.D0/DFLOAT(P))

      if ((C(14) .gt. 0.d0) .and. (c(14) .le. 1.d-3)) then 
         hscale = 50.d0
      else 
         HSCALE = 2.D0
      end if

      GO TO 480

 460  IF (C(23) .GT. 1.D0) GO TO 490
      IF (C(23) .EQ. 1.D0) HSCALE = 2.D0

C ... Case 2 - After a successful step, or at most  one  failure,
C              Use MIN(2, .9*(TOL/EST)**(1/P))*HMAG, but avoid possible
C              overflow. Then avoid reduction by more than half.

      TEMP = HSCALE * C(14)
      IF (TOL .LT. (HSCALE/.9D0)**P*C(19))
     +     TEMP = .9D0*(TOL/C(19))**(1.D0/DFLOAT(P))*C(14)

      C(14) = DMAX1(TEMP, .5D0*C(14))
      GO TO 480
 490  CONTINUE

C ... Case 3 - After two or more successive failures
      C(14) = .5D0*C(14)
 480  CONTINUE

C ... Check against HMAX
      C(14) = DMIN1(C(14), C(16))

C ... Check against HMIN
      C(14) = DMAX1(C(14), C(13))

C***********INTERRUPT NO 1 (WITH IND=4) IF REQUESTED
      IF (C(8) .EQ. 0.D0) GO TO 1111
      IND = 4
      RETURN

C ... Resume here on re-entry with IND = 4   ........RE-ENTRY..
 1111 CONTINUE

C ... Calculate HMAG, TTRIAL - Depending on preliminary HMAG, TEND
      IF (C(14) .GE. DABS(TEND - T)) GO TO 500

C ... Do not step more than half way to TEND
      C(14) = DMIN1(C(14), .5D0*DABS(TEND - T))
      C(17) = T + DSIGN(C(14), TEND - T)
      GO TO 510
 500  CONTINUE

C ... Hit TEND exactly
      C(14) = DABS(TEND - T)
      C(17) = TEND
 510  CONTINUE

C ... Calculate HTRIAL
      C(18) = C(17) - T

C ... END STAGE 1

C*********************************************************************
C   Stage 2 - Calculate YTRIAL using the subroutine APFORM
 2500 continue

      H = C(18)

      if ((C(32) .eq. 1) .or. (C(32) .eq. 2)) then

         CALL APFORM(N,NW,MU,OMEGA,LOLD,SSTAR,T,Y,H,TOLD,KOLD,TOL,W,C, 
     +            FCN,SIGMA,YINIT,DYINIT) 

      else if (C(32) .eq. 3) then
         
         CALL APFORM_ASM(N,NW,MU,OMEGA,LOLD,SSTAR,T,Y,H,TOLD,KOLD,TOL,W
     +        ,C,FCN,SIGMA,YINIT,DYINIT) 

      else 
         
         print *, 'Parameter Error in C(32)'
         stop

      end if

         

C ... Add SSTAR to the no of fcn evals
          C(24) = C(24) + DFLOAT(SSTAR - 1)

C ... End stage 2
C*********************************************************************

C*********************************************************************
C   Stage 3 - Calculate the error estimate EST. First calculate 
C   the DEFECT at TAUS and use the norm of this vector to control
C   the stepsize. Note that W(*,2) and W(*,14) are used for temp
C   storage.

          IF ((C(32) .eq. 1) .or. (C(32) .eq. 3)) then
          
             CALL DEFECT(N,NW,MU,OMEGA,LOLD,T,T+TAUS*H,Y,C(19),TOLD,KOLD
     +            ,W,C,FCN,SIGMA,YINIT,DYINIT)

          else if (C(32) .eq. 2) then

             CALL DEFECT2(N,NW,MU,OMEGA,LOLD,T,T+TAUS*H,T+TAUS2*H,Y,
     +            C(19),TOLD,KOLD,W,C,FCN,SIGMA,YINIT,DYINIT)

          else

             print *, 'Parameter Error in C(32)'
             stop

          end if



C ... END Stage 3
C*********************************************************************

C*********************************************************************
C   Stage 4 - Make decisions.

C ... Set IND=5 if step acceptable, else set IND=6

         if (c(19) .le. tol) then
c         if ((dabs(w(1,14)) .le. tol) .and. 
c     +        (dmax1(dabs(w(2,14)), dabs(w(3,14))) .le. 1.d-4)) then

            if (ind .ne. 8) ind = 5
         else
            if (ind .ne. 8) then
               ind = 6
            else
               ind = 9
            end if
         end if



C ... If accepted then, we can use extrap in the next step,so C(29) = 0
      if ((ind .eq. 8) .or. (ind .eq. 9)) return   
      IF (IND .EQ. 5)  C(29) = 0.D0

C***********INTERRUPT NO 2 IF REQUESTED
      IF (C(9) .EQ. 0.D0) GO TO 2222
      RETURN

C ... Resume here on re-entry with IND = 5 OR 6   ...RE-ENTRY..
 2222 CONTINUE

      IF (IND .EQ. 6) GO TO 600

C ... Step accepted (IND .EQ. 5), so update T, Y from TTRIAL,
C     YTRIAL, add 1 to the no of successful steps, and set the no of
C     successive failures to zero
      T = C(17)
      DO 610 K = 1, N
         Y(K) = W(K,12)
 610  CONTINUE
      C(22) = C(22) + 1.D0
      C(23) = 0.D0

C**************RETURN(WITH IND=3, TEND SAVED, FLAG SET) IF T .EQ. TEND
      IF (T .NE. TEND) GO TO 620
      IND = 3
      C(20) = TEND
      C(21) = 1.D0
      RETURN
 620  CONTINUE
      GO TO 630
 600  CONTINUE

C ... Step not accepted (IND .EQ. 6), so add 1 to the no of
C     successive failures
      C(23) = C(23) + 1.D0

C**************ERROR RETURN (WITH IND=-3) IF HMAG .LE. HMIN
      IF (C(14) .GT. C(13)) GO TO 640
      IND = -3
      RETURN
 640  CONTINUE
 630  CONTINUE

C END STAGE 4
C***********************************************************************

      GO TO 99999
C ... End loop

C ... Begin abort action
 3000  CONTINUE

       WRITE(6,3010) IND, TOL, T, N, C(13), TEND, NW, C(16), C(20),
     +      C(22), C(23), C(24), (Y(K), K = 1, N)
 3010  FORMAT( /// 1H ,58HCOMPUTATION STOPPED IN DNORK WITH THE FOLLOWIN
     +G VALUES -
     +   / 1H , 5HIND =, I4, 5X, 6HTOL  =, 1PD13.6, 5X, 11HX         =,
     +          1PD22.15
     +   / 1H , 5HN   =, I4, 5X, 6HHMIN =, 1PD13.6, 5X, 11HTEND      =,
     +          1PD22.15
     +   / 1H , 5HNW  =, I4, 5X, 6HHMAX =, 1PD13.6, 5X, 11HPREV TEND =,
     +          1PD22.15
     +   / 1H , 14X, 27HNO OF SUCCESSFUL STEPS    =, 0PF8.0
     +   / 1H , 14X, 27HNO OF SUCCESSIVE FAILURES =, 0PF8.0
     +   / 1H , 14X, 27HNO OF FUNCTION EVALS      =,0PF8.0
     +   / 1H , 23HTHE COMPONENTS OF Y ARE
     +   // (1H , 1P5D24.15))

      STOP

C ... End abort action

      END


C *** End RDMETH ***





      SUBROUTINE STORE(N,NW,LOLD,T,Y,IFIRST,ILAST,TOLD,KOLD,W)

C*********************************************************************
C                                                                    *
C     Purpose - This routine is for storing T value, K-values, and   *
C               solutions to TOLD and KOLD queue.                    *
C               When K-values for Tn-1 to Tn and Yn have             *
C               computed, store                                      *
C               TOLD(ILAST+1) <--- Tn                                *
C               KOLD(*,1,ILAST+1) <--- Yn-1                          *
C               KOLD(*,2,ILAST+1) <--- K1                            *
C               KOLD(*,3,ILAST+1) <--- Yn                            *
C               KOLD(*,4,ILAST+1) <--- Ks+1(K8)                      *
C               KOLD(*,5,ILAST+1) <--- Ks+2(K9)                      *
C               KOLD(*,6,ILAST+1) <--- Ks+3(K10)                     *
C               KOLD(*,7,ILAST+1) <--- Ks+4(K11)                     *
C                                                                    *
C*********************************************************************
C                                                                    *
C    N      Number of equation                                       *
C    NW     First dimension of workspace W & queue KOLD              *
C    LOLD   Length of queue                                          *
C    T      Independent value (Tn)                                   *
C    Y      Solutions at Tn-1                                        *
C    IFIRST First index of queue                                     *
C    ILAST  Last index of queue                                      *
C    TOLD   Queue of independent values                              *
C    KOLD   Queue of k-values and solutions                          *
C    W      Workspace                                                *
C                                                                    *
C*********************************************************************


C     .. Scalar Arguments ..
      INTEGER N,NW,LOLD
      DOUBLE PRECISION T,IFIRST,ILAST

C     .. Array Arguments ..
      DOUBLE PRECISION Y(N),TOLD(LOLD),KOLD(NW,7,LOLD),W(NW,*)

C     .. Local Scalars ..
      INTEGER I,M,IFST,ILST

C     .. Intrisic Functions ..
      INTRINSIC IDINT,DFLOAT


      IFST = IDINT(IFIRST)
      ILST = IDINT(ILAST)

C ... CASE 1:  IFST <= ILST & ILST < LOLD
      IF ((ILST .GE. IFST) .AND. (ILST .LT. LOLD)) THEN
         M = ILST + 1
         ILAST = ILAST + 1.D0

C ... CASE 2: IFST < ILST = LOLD
      ELSE IF (ILST .EQ. LOLD) THEN
         print *, 'Warning:the length of the queue is short!' 
         M = 1
         IFIRST = 2.D0
         ILAST = 1.D0

C ... CASE 3: ILST < IFST < LOLD
      ELSE IF (IFST .LT. LOLD) THEN
         M = ILST + 1
         IFIRST = IFIRST + 1.D0
         ILAST = ILAST + 1.D0

C ... CASE 4: ILST < IFST = LOLD
      ELSE 
         M = LOLD
         IFIRST = 1.D0
         ILAST = DFLOAT(LOLD)
      END IF

C ... Store T value
      TOLD(M) = T

C ... Store K values and soltuion
      DO 100 I = 1, N
         KOLD(I,1,M) = Y(I)
         KOLD(I,2,M) = W(I,1)
         KOLD(I,3,M) = W(I,12)
         KOLD(I,4,M) = W(I,8) 
         KOLD(I,5,M) = W(I,9) 
         KOLD(I,6,M) = W(I,10) 
         KOLD(I,7,M) = W(I,11)   
 100  CONTINUE

c      print *, 'm =', m
c      print *, 'told =', told(m)
c      print *, 'kold1 =', (kold(1,i,m), i= 1, 7)
c      print *, 'kold2 =', (kold(2,i,m), i= 1, 7)
c      print *, 'kold3 =', (kold(3,i,m), i= 1, 7) 

      RETURN
      END



C *** End of STORE ***





      SUBROUTINE STORED(N,NW,LOLD,T,TS,Y,YS,TOLD,KOLD,W,C)

C*********************************************************************
C                                                                    *
C     Purpose - This routine is for storing T value, derivatives,    *
c               and solutions to TOLD and KOLD queue just after      *
c               finding discontinuity.                               *
C                                                                    *
C               TOLD(ILAST+1) <--- Ts                                *
C               KOLD(*,1,ILAST+1) <--- Yn-1                          *
C               KOLD(*,2,ILAST+1) <--- K1                            *
C               KOLD(*,3,ILAST+1) <--- Ys                            *
C               KOLD(*,4,ILAST+1) <--- Z'(Ts)                        *
C               KOLD(*,5,ILAST+1) <--- Z'(T+C9*Hs)                   *
C               KOLD(*,6,ILAST+1) <--- Z'(T+C10*Hs)                  *
C               KOLD(*,7,ILAST+1) <--- Z'(T+C11*Hs)                  *
C                                                                    *
C*********************************************************************
C                                                                    *
C    N      Number of equation                                       *
C    NW     First dimension of workspace W & queue KOLD              *
C    LOLD   Length of queue                                          *
C    T      Independent value (left mesh point)                      * 
C    TS     Independent value (right mesh point)                     *
C    Y      Solution at T                                            *
C    YS     Solution at TS                                           *
C    TOLD   QUEUE OF INDEPENDENT VALUES                              *
C    KOLD   QUEUE OF K-VALUES AND SOLUTIONS                          *
C    W      WORKSPACE                                                *
C    C      COMMUNICATION ARRAY                                      *
C                                                                    *
C*********************************************************************

C     .. Scalar Arguments ..
      INTEGER N,NW,LOLD
      DOUBLE PRECISION T,TS

C     .. Array Arguments ..
      DOUBLE PRECISION Y(N),YS(N),TOLD(LOLD),KOLD(NW,7,LOLD)
      DOUBLE PRECISION W(NW,*),C(*)

C     .. Local Scalars ..
      INTEGER I,IFST,ILST,M
      DOUBLE PRECISION  C9,C10,C11,H,HS

C     .. External Subroutines ..
      EXTERNAL DERIV

C     .. Intrisic Functions ..
      INTRINSIC IDINT,DFLOAT

C     .. Data statements ..
      DATA C9/0.5D0/
      DATA C10/0.1428571428571428D0/
      DATA C11/0.2857142857142857D0/

      H = C(18)
      HS = TS - T

      IFST = IDINT(C(26))
      ILST = IDINT(C(27))

C ... CASE 1:  IFST <= ILST & ILST < LOLD
      IF ((ILST .GE. IFST) .AND. (ILST .LT. LOLD)) THEN
         M = ILST + 1
         C(27) = C(27) + 1.D0

C ... CASE 2: IFST < ILST = LOLD
      ELSE IF (ILST .EQ. LOLD) THEN
         print *, 'Warning:the length of the queue is short!'
         M = 1
         C(26) = 2.D0
         C(27) = 1.D0

C ... CASE 3: ILST < IFST < LOLD
      ELSE IF (IFST .LT. LOLD) THEN
         M = ILST + 1
         C(26) = C(26) + 1.D0
         C(27) = C(27) + 1.D0

C ... CASE 4: ILST < IFST = LOLD
      ELSE 
         M = LOLD
         C(26) = 1.D0
         C(27) = DFLOAT(LOLD)
      END IF

C ... Store T value
      TOLD(M) = TS

C ... Store Yn-1, K1, Ys, F(Ts, Ys, DEL) to KOLD
      DO 20 I = 1, N
         KOLD(I,1,M) = Y(I)
         KOLD(I,2,M) = W(I,1)
         KOLD(I,3,M) = YS(I)
 20   CONTINUE

C ... Store Z'(T+C9*Hs), Z'(T+C10*Hs), Z'(T+C11*Hs) to KOLD
      CALL DERIV(N,NW,T,Y,TS,KOLD(1,4,M),H,W)
      CALL DERIV(N,NW,T,Y,T+C9*HS,KOLD(1,5,M),H,W)
      CALL DERIV(N,NW,T,Y,T+C10*HS,KOLD(1,6,M),H,W)
      CALL DERIV(N,NW,T,Y,T+C11*HS,KOLD(1,7,M),H,W) 

      RETURN
      END


C *** End of STORED ***





      SUBROUTINE TFIND(N,TOLD,T,IFIRST,ILAST,M)

C****************************************************************
C                                                               *
C     Purpose - This subroutine finds the index for the         *
C               independent valriable from the queue of the     *
C               past values by using binary search.             *
C                                                               *
C     If TOLD(M) <= T < TOLD(M+1), M is returned.               *
C                                                               *
C****************************************************************
C                                                               *
C     N       The length of queue                               *
C     TOLD    Queue of past independent variable                *
C     T       Independent variavle - wish to find the indet     *
C     IFIRST  The first index of the queue                      *
C     ILAST   The last index of the queue                       *
C     M       Index found (Output)                              *
C                                                               *
C****************************************************************

C     .. Scalar Arguments ..
      INTEGER N,M
      DOUBLE PRECISION T,IFIRST,ILAST

C     .. Array Arguments ..
      DOUBLE PRECISION TOLD(N)

C     .. Local Scalars ..
      INTEGER I,J,K,IFST,ILST,REALK,RK,RKP,RKM

C     .. Intrisic Functions ..
      INTRINSIC IDINT,DFLOAT


C ... Initialize variables
      IFST = IDINT(IFIRST)
      ILST = IDINT(ILAST)
      I = IFST

C ... Calculate the last index J for binary search
      IF (ILST .EQ. IFST) THEN
         PRINT *, '***** ERROR:TFIND ***** '
         PRINT *, 'No data is stored in the queue.:'
         Print *, 'T =', t
         PRINT *  
         STOP
      ELSE IF (IFST .LT. ILST) THEN
         J = ILST
      ELSE
         J = N + ILST
      END IF

C ... Loop - Binary search
 100  CONTINUE
C ... Calculate the middel of I and J
         K = (I + J) / 2

C ... Calculate real index for K, K+1, and K-1
         RK  = REALK(K, N)
         RKP = REALK(K+1, N)
         RKM = REALK(K-1, N)

C ... CASE 1:  TOLD(RK) = T
         IF (TOLD(RK) .EQ. T) THEN
            M = RK
            RETURN

C ... CASE 2:  TOLD(RK) < T
         ELSE IF (TOLD(RK) .LT. T) THEN
           
C ... CASE 2.1: IF RK IS ILST THEN MUST STOP
            IF (RK .EQ. ILST) THEN
               PRINT *, '***** ERROR:TFIND ***** '
               PRINT *, 'The delay is advanced than the history'
               print *, 'T =', t
               PRINT *
               STOP
           
C ... CASE 2.2: TOLD(RK) < T < TOLD(RK+1) THEN M = RK
            ELSE IF (T .LT. TOLD(RKP)) THEN
               M = RK
               RETURN

C ... CASE 2.3: TOLD(RK+1) < T , CHANGE I TO K + 1
            ELSE 
               I = K + 1
            END IF

C ... CASE 3: NOW T < TOLD(RK)
C ... CASE 3.1: IF K IS IFST THEN MUST STOP
         ELSE IF (K .EQ. IFST) THEN
            PRINT *, ' ***** ERROR:TFIND ***** '
            PRINT *, 'The length of the queue is short'
            PRINT *, 'T =', t, ', k =', k, ', Told(k) =', told(k)
            PRINT *
            STOP

C ... CASE 3.2: TOLD(RK-1) < T < TOLD(RK) THEM M = RK - 1
         ELSE IF (TOLD(RKM) .LT. T) THEN
            M = RKM
            RETURN

C ... CASE 3.3: T < TOLD(RK-1)
         ELSE
            J = K - 1
         END IF
      GO TO 100
C ... END LOOP - BINARY SEARCH

      END



C *** End of TFIND ***






      FUNCTION REALK(K,N)

C*****************************************************************
C                                                                *
C     Purpose - This function computes real index of the queue.  *
C                                                                *
C*****************************************************************
C                                                                *
C   K      The index                                             *
C   N      The length of queue                                   *
C   REALK  Real index (Output)                                   *
C                                                                *
C*****************************************************************

C     .. Scalar Arguments ..
      INTEGER K, N, REALK

C ... If the K is between 1 to N then K is the real index
      IF (K .LE. N) THEN
         REALK = K

C ... If the K is bigger than N then return K MOD N
      ELSE
         REALK = MOD(K, N)
      END IF

      RETURN
      END	



C *** End of REALK





