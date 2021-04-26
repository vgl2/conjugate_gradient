      FUNCTION FUNC(X)
      PARAMETER (NDIM=27)
      PARAMETER (NMONO=3)
      DIMENSION X(NDIM)
      X = X*0.52917721067
      CALL CALCPOT(NMONO,V,X)
      X = X/0.52917721067
      FUNC = DBLE(V/627.509474)
      RETURN
      END 
      SUBROUTINE DFUNC(X,DF) 
      PARAMETER (NDIM=27)
      PARAMETER (NMAX=50) 
      PARAMETER (NMONO=3)
      DIMENSION X(NDIM),DF(NMAX),COORD(NDIM),TEST(5)
      DX = 0.0001
      COORD = X
      DO 10 K = 1,5
        TEST(K) = 0.D0
10    CONTINUE
      DO 11 J=1,NDIM
C		Use below if statement for partial optimization
C        IF ((J.EQ.3).OR.(J.EQ.9).OR.(J.EQ.12).OR.(J.EQ.21)) THEN
          COORD(:) = COORD(:)*0.52917721067
          CALL CALCPOT(NMONO,V,COORD)
          TEST(1) = V/627.509474
          COORD(:) = COORD(:)/0.52917721067
          TEST(2) = TEST(1)
          TEST(3) = TEST(1)
          TEST(4) = TEST(1)
          TEST(5) = TEST(1)
        ELSE
          DO 12 K = 1,5
            IF (K.EQ.1) THEN
              COORD(J) = COORD(J) + 2*DX
              COORD(:) = COORD(:)*0.52917721067
              CALL CALCPOT(NMONO,V,COORD)
              TEST(1) = V/627.509474
              COORD(:) = COORD(:)/0.52917721067
              COORD(J) = COORD(J) - 2*DX
            ELSE IF (K.EQ.2) THEN
              COORD(J) = COORD(J) + DX
              COORD(:) = COORD(:)*0.52917721067
              CALL CALCPOT(NMONO,V,COORD)
              TEST(2) = V/627.509474
              COORD(:) = COORD(:)/0.52917721067
              COORD(J) = COORD(J) - DX
            ELSE IF (K.EQ.3) THEN
              COORD(J) = COORD(J) 
              COORD(:) = COORD(:)*0.52917721067
              CALL CALCPOT(NMONO,V,COORD)
              TEST(3) = V/627.509474
              COORD(:) = COORD(:)/0.52917721067
              COORD(J) = COORD(J) 
            ELSE IF (K.EQ.4) THEN
              COORD(J) = COORD(J) - DX
              COORD(:) = COORD(:)*0.52917721067
              CALL CALCPOT(NMONO,V,COORD)
              TEST(4) = V/627.509474
              COORD(:) = COORD(:)/0.52917721067
              COORD(J) = COORD(J) + DX
            ELSE IF (K.EQ.5) THEN
              COORD(J) = COORD(J) -2*DX
              COORD(:) = COORD(:)*0.52917721067
              CALL CALCPOT(NMONO,V,COORD)
              TEST(5) = V/627.509474
              COORD(:) = COORD(:)/0.52917721067
              COORD(J) = COORD(J) + 2*DX
            ENDIF
12        CONTINUE   
        ENDIF
        DF(J) =(((-1./12.)*TEST(1))+((2./3.)*TEST(2))+((-2./3.)*TEST(4))
     1  +((1./12.)*TEST(5)))/DX
        DO 13 K = 1,5
           TEST(K) = 0.D0
13      CONTINUE
11    CONTINUE
      RETURN 
      END 
      SUBROUTINE FRPRMN(P,N,FTOL,ITER,FRET) 
      PARAMETER (NMAX=100,ITMAX=100000,EPS=1.E-30) 
      DIMENSION P(N),G(NMAX),H(NMAX),XI(NMAX)
      FP=FUNC(P) 
      CALL DFUNC(P,XI) 
      DO 11 J=1,N 
        G(J)=-XI(J) 
        H(J)=G(J) 
        XI(J)=H(J) 
11    CONTINUE
      DO 14 ITS=1,ITMAX 
        ITER=ITS
        CALL LINMIN(P,XI,N,FRET) 
        IF(2.*ABS(FRET-FP).LE.FTOL*(ABS(FRET)+ABS(FP)+EPS))RETURN 
        FP=FUNC(P) 
        CALL DFUNC(P,XI) 
        GG=0. 
        DGG=0. 
        DO 12 J=1,N 
          GG=GG+G(J)**2 
         DGG=DGG+XI(J)**2 
12      CONTINUE 
        IF(GG.EQ.0.)RETURN 
        GAM=DGG/GG 
        DO 13 J=1,N 
          G(J)=-XI(J) 
          H(J)=G(J)+GAM*H(J) 
          XI(J)=H(J) 
13      CONTINUE 
14    CONTINUE 
      PAUSE 'FRPR MAXIMUM ITERATIONS EXCEEDED' 
      RETURN 
      END 
       SUBROUTINE LINMIN(P,XI,N,FRET) 
      PARAMETER (NMAX=50,TOL=1.E-4) 
      EXTERNAL F1DIM 
      DIMENSION P(N),XI(N) 
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX) 
      NCOM=N 
      DO 11 J=1,N 
        PCOM(J)=P(J) 
        XICOM(J)=XI(J) 
11    CONTINUE 
      AX=0. 
      XX=1. 
      BX=2. 
      CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM) 
      FRET=BRENT(AX,XX,BX,F1DIM,TOL,XMIN) 
      DO 12 J=1,N 
        XI(J)=XMIN*XI(J) 
        P(J)=P(J)+XI(J) 
12    CONTINUE 
      RETURN 
      END 
            SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC) 
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.E-20) 
      FA=FUNC(AX) 
      FB=FUNC(BX) 
      IF(FB.GT.FA)THEN 
        DUM=AX 
        AX=BX 
        BX=DUM 
        DUM=FB 
        FB=FA 
        FA=DUM 
      ENDIF 
      CX=BX+GOLD*(BX-AX) 
      FC=FUNC(CX) 
1     IF(FB.GE.FC)THEN 
        R=(BX-AX)*(FB-FC) 
        Q=(BX-CX)*(FB-FA) 
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q- 
     1        R)) 
        ULIM=BX+GLIMIT*(CX-BX) 
        IF((BX-U)*(U-CX).GT.0.)THEN 
          FU=FUNC(U) 
          IF(FU.LT.FC)THEN 
            AX=BX 
            FA=FB 
            BX=U 
            FB=FU 
            GO TO 1 
          ELSE IF(FU.GT.FB)THEN 
            CX=U 
            FC=FU 
            GO TO 1 
          ENDIF 
          U=CX+GOLD*(CX-BX) 
          FU=FUNC(U) 
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN 
          FU=FUNC(U) 
          IF(FU.LT.FC)THEN 
            BX=CX 
            CX=U 
            U=CX+GOLD*(CX-BX) 
            FB=FC 
            FC=FU 
            FU=FUNC(U) 
          ENDIF 
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN 
          U=ULIM 
          FU=FUNC(U) 
        ELSE 
          U=CX+GOLD*(CX-BX) 
          FU=FUNC(U) 
        ENDIF 
        AX=BX 
        BX=CX 
        CX=U 
        FA=FB 
        FB=FC 
        FC=FU 
        GO TO 1 
      ENDIF 
      RETURN 
      END 
        FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN) 
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10) 
      A=MIN(AX,CX) 
      B=MAX(AX,CX) 
      V=BX 
      W=V 
      X=V 
      E=0. 
      FX=F(X) 
      FV=FX 
      FW=FX 
      DO 11 ITER=1,ITMAX 
        XM=0.5*(A+B) 
        TOL1=TOL*ABS(X)+ZEPS 
        TOL2=2.*TOL1 
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3 
        IF(ABS(E).GT.TOL1) THEN 
          R=(X-W)*(FX-FV) 
          Q=(X-V)*(FX-FW) 
          P=(X-V)*Q-(X-W)*R 
          Q=2.*(Q-R) 
          IF(Q.GT.0.) P=-P 
          Q=ABS(Q) 
          ETEMP=E 
          E=D 
          IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR.  
     *        P.GE.Q*(B-X)) GOTO 1 
          D=P/Q 
          U=X+D 
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X) 
          GOTO 2 
        ENDIF 
1       IF(X.GE.XM) THEN 
          E=A-X 
        ELSE 
          E=B-X 
        ENDIF 
        D=CGOLD*E 
2       IF(ABS(D).GE.TOL1) THEN 
          U=X+D 
        ELSE 
          U=X+SIGN(TOL1,D) 
        ENDIF 
        FU=F(U) 
        IF(FU.LE.FX) THEN 
          IF(U.GE.X) THEN 
            A=X 
          ELSE 
            B=X 
          ENDIF 
          V=W 
          FV=FW 
          W=X 
          FW=FX 
          X=U 
          FX=FU 
        ELSE 
          IF(U.LT.X) THEN 
            A=U 
          ELSE 
            B=U 
          ENDIF 
          IF(FU.LE.FW .OR. W.EQ.X) THEN 
            V=W 
            FV=FW 
            W=U 
            FW=FU 
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN 
            V=U 
            FV=FU 
          ENDIF 
        ENDIF 
11    CONTINUE 
      PAUSE 'BRENT EXCEED MAXIMUM ITERATIONS.' 
3     XMIN=X 
      BRENT=FX 
      RETURN 
      END 
       FUNCTION F1DIM(X) 
      PARAMETER (NMAX=50) 
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX) 
      DIMENSION XT(NMAX) 
      DO 11 J=1,NCOM 
        XT(J)=PCOM(J)+X*XICOM(J) 
11    CONTINUE 
      F1DIM=FUNC(XT) 
      RETURN 
      END 
