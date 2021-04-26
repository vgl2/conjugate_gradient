      SUBROUTINE calc_grad(x,df) 
      PARAMETER (NDIM=27)
      PARAMETER (NMAX=50) 
      PARAMETER (NMONO=3)
      DIMENSION X(NDIM),DF(NMAX),COORD(NDIM),TEST(5)
C	  Calculates the gradient of the energy in respect
C	  to coordinates for a molecule and a given potential
C	  energy surface.
C	  Inputs:
C	  x = geometry of molecule
C     Outputs:
C     df = dE/dx_i of the coordinates of the molecule
      DX = 0.0001
      COORD = X
      DO 10 K = 1,5
        TEST(K) = 0.D0
10    CONTINUE
      DO 11 J=1,NDIM
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
        DF(J) =(((-1./12.)*TEST(1))+((2./3.)*TEST(2))+((-2./3.)*TEST(4))
     1  +((1./12.)*TEST(5)))/DX
        DO 13 K = 1,5
           TEST(K) = 0.D0
13      CONTINUE
11    CONTINUE
      RETURN 
