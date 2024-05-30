C          ----------------LAPLACE_3D_TEST.F-----------------

C   THIS FORTRAN CODE IS USED TO SOLVE THE LAPLACE EQUATION FOR THE
C   MAGNETIC SCALAR POTENTIAL IN A SPHERICAL COORDINATE SYSTEM. IT
C   CHAS BEEN DEVELOPED USING FISHPACK PACKAGE FROM HE CLASSIC 
C   LIBRARIES FOR GEOPHYSICS , WHICH WAS CREATED BY THE NATIONAL CENTER
C   FOR ATMOSPHERIC RESEARCH (NCAR). THE PACKAGE CONTAINS EFFICIENT
C   FORTRAN SUBPROGRAMS FOR SOLVING ELLIPTIC PARTIAL DIFFERENTIAL
C   EQUATIONS. SPECIFICALLY, WE USE THE SUBROUTINE HWSCSP, WHICH
C   SOLVES A FINITE DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
C   EQUATION IN SPHERICAL COORDINATES UNDER THE ASSUMPTION OF
C   AXISYMMETRY (NO DEPENDENCE ON LONGITUDE). WE HAVE MODIFIED THE
C   CODE TO SOLVE FOR A 3D SPHERICAL COORDINATE SYSTEM THAT IS
C   PERIODIC IN LONGITUDE.

C---------------------------------------------------------------------- 

C     THE DIMENSIONS REQUIRED FOR YOUR PROBLEM SHOULD BE SATISFY WITH 
C     CONDITIONS OF THE SUBRUTINE HWSCSP. YOU CAN REFER THE TECHNICAL 
C     NOTES.

C     FOR THIS CODE THE GRIDPOINTS ARE MENSION BELOW
C     THETA : 160   PHI : 288   R : 511   PHI_DOUBLED :: 575
C----------------------------------------------------------------------

	REAL :: WSAVE(2315)
        INTEGER :: I,J,K,H
        INTEGER, PARAMETER :: n_rows = 575
        INTEGER, PARAMETER :: n_cols = 160    
        REAL :: input_data(575*160) 
        REAL, DIMENSION(n_rows, n_cols) :: Brrs
        real, DIMENSION(n_cols,n_rows) ::  Br_rs
        COMPLEX :: A, B_r_THETA(575),B_K(160,575),B
        character(200)  f1
        real :: B_r(575,160), B_r_actual(160,575), PHI(575)   
        complex :: G(160,511,575) ,Gi(160,511,575)
        COMPLEX :: U(160,511,575) , FREQ_U(575)
        real :: BDTS(511), BDRS(160), n1, d1, n2, d2  
        DIMENSION      F(160,511), BDTF(511), W(13206), R(511), 
     &   THETA(160), BDRF(160)
     
C-----------------------------------------------------------------------

C    'f1' IS THE PATH OF THE INPUT FILE THAT CONTAINS THE RADIAL 
C     COMPONENT OF MAGNETIC FIELD AT THE LOWER BOUNDARY (Br_rs).
C     MAKE SURE THAT THE LOWER BOUNDARY FUNCTION IS SYMMETRIC ALONG 
C     LONGITUDE.IF NOT YOU CAN CREATE SYMMETRIC USING THE PYTHON CODE   
C     'Doubling_for_fourier.ipynb' IN THE FOLDER.

C-----------------------------------------------------------------------

C          ----------------SETTING UP THE FILES-----------------

        f1 = '/home/ijas/Fortran_IIA/testing_code/D_B/BrVAR74D.txt' 
        OPEN(UNIT=10, FILE=f1) 
        DO i = 1, n_rows*n_cols
            READ(10, *) input_data(i)
        END DO
        CLOSE(10)

        DO i = 1, n_rows	! the file was actually in phi,theta
            DO j = 1, n_cols	! so we need that file in theta,phi
                Brrs(i, j) = input_data((i - 1) * n_cols + j)
            END DO
        END DO
        
        DO i = 1, n_rows
            DO j = 1, n_cols
                Br_rs(j, i) = Brrs(i, j)
            END DO
        END DO
        PRINT*, '---Br_rs---'
        PRINT*, 'MAX_VAL', MAXVAL(Br_rs), 'MIN_VAL', MINVAL(Br_rs)

C          ---------------FORWARD FOURIER TRANSFORM----------------

	A = (1,0)
	B = (0,1)
        call CFFTI(575,WSAVE)
        DO I = 1, 160
         DO J=1,575
          B_r_THETA(J) = A * Br_rs(i, j)
         END DO
         !PRINT*, 1,J,B_r_THETA
         call CFFTF(575,B_r_THETA,WSAVE)
         DO K =1 , 575
          B_K(I,K) = B_r_THETA(K)/575	
         END DO
        end do
        PRINT*, '---REAL B_K---'
        PRINT*,'MAX_VAL',MAXVAL(REAL(B_K)), 'MIN_VAL',MINVAL(REAL(B_K))
        PRINT*, '---AIMAG B_K---'
        PRINT*,'MAX',MAXVAL(AIMAG(B_K)),'MIN',MINVAL(AIMAG(B_K))       

C######################################################################

C          ---------------HWSCSP FOR 3D_REAL----------------
 
        DO 1 K =1 , 575  ! THIS IS THE FOR EACH MODES

        PI = PIMACH(DUM)
        INTL = 0
        PS = -0.49087384
        PF = 0.49087384
        L = 574
        TS = 1.3089970
        TF = 1.8325958
        M = 159 
        MBDCND = 3
        RS = 1.
        RF = 6.
        N = 510 
        NBDCND = 4
        
c-----------------------------------------------------------------------      
c       DEFINING THE GRID POINTS IN PHI DIRECTION 
c-----------------------------------------------------------------------    
      
        LP1 = L+1
        DP = (PF-PS)/(287)
        DO 100 H=1,LP1
         PHI(H) = FLOAT(H-1)*DP + PS
  100   END DO
  
        ELMBDA=-(2*(1-COS(((2*PI)*(FLOAT(K)-1)/575))))/(DP**2)
        IDIMF = 160
        
c-----------------------------------------------------------------------      
c       DEFINING THE GRID POINTS AND VALUES ON THETA AND R
c-----------------------------------------------------------------------  

        MP1 = M+1
        DTHETA = (TF-TS)/FLOAT(M)   
        DO 101 I=1,MP1
         THETA(I) = FLOAT(I-1)*DTHETA + TS
  101   END DO
  
        NP1 = N+1
        DR = (RF-RS)/FLOAT(N)
        DO 102 J=1,NP1
         R(J) = FLOAT(J-1)*DR +RS
  102   END DO
  
c-----------------------------------------------------------------------  
c       NOW THE BOUNDARY CONDITIONS IN 2D IS DEFINED
c-----------------------------------------------------------------------     
      
      DO 3 J=1,NP1
         BDTS(J) = 0.
   3   END DO
     
      DO 4 J=1,NP1
         BDTF(J) = 0.
   4  END DO
      
      DO 5 I=1,MP1
         BDRS(I) = real(B_K(I,K))
   5   END DO
 
      
      DO 6 I=1,MP1
         F(I,N+1) = 0.
   6  END DO
   
c-----------------------------------------------------------------------  
c       DEFIENING THE RHS SIDE OF THE HELMHOLTZ EQUATION
c-----------------------------------------------------------------------  

      DO 106 I=1,MP1
         DO 105 J=1,N
            F(I,J) = 0.
  105    CONTINUE
  106 CONTINUE 
  
c----------------------------------------------------------------------- 
c       SOLVING FOR EACH PANEL FOR COEFFIENTS OF FOURIER TRANSFORM
c-----------------------------------------------------------------------  

      CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     1             BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)

      DO I = 1, 160
         DO J = 1,  511 
            G(I,J,K) = F(I,J)*A	
         end do		
      end do		! SOLVES AS REAL COMPONENT OF COMPLEX NUMBER

      
    1 end do ! K LOOP END
    
C######################################################################

C          ---------------HWSCSP FOR 3D_IMAGINARY----------------

        DO 60 K =1 , 575     ! THIS IS THE FOR EACH MODES
     
        PI = PIMACH(DUM)
        INTL = 0
        PS = -0.49087384
        PF = 0.49087384
        L = 574
        TS = 1.3089970
        TF = 1.8325958
        M = 159 
        MBDCND = 3
        RS = 1.
        RF = 6.
        N = 510 
        NBDCND = 4
         
c-----------------------------------------------------------------------      
c       DEFINING THE GRID POINTS IN PHI DIRECTION 
c-----------------------------------------------------------------------   
       
        LP1 = L+1
        DP = (PF-PS)/(287)
        DO 50 H=1,LP1
         PHI(H) = FLOAT(H-1)*DP + PS
   50   END DO
  
        ELMBDA=-(2*(1-COS(((2*PI)*(FLOAT(K)-1)/575))))/(DP**2)
        IDIMF = 160
        
c-----------------------------------------------------------------------      
c       DEFINING THE GRID POINTS AND VALUES ON THETA AND R
c-----------------------------------------------------------------------       

        MP1 = M+1
        DTHETA = (TF-TS)/FLOAT(M)   
        DO 51 I=1,MP1
         THETA(I) = FLOAT(I-1)*DTHETA + TS
   51   END DO
  
        NP1 = N+1
        DR = (RF-RS)/FLOAT(N)
        DO 52 J=1,NP1
         R(J) = FLOAT(J-1)*DR +RS
   52   END DO
  
c-----------------------------------------------------------------------  
c       NOW THE BOUNDARY CONDITIONS IN 2D IS DEFINED
c-----------------------------------------------------------------------        
      
      DO 53 J=1,NP1
         BDTS(J) = 0.
  53   END DO
     
      DO 54 J=1,NP1
         BDTF(J) = 0.
  54  END DO
      
      DO 55 I=1,MP1
         BDRS(I) = AIMAG(B_K(I,K))
  55   END DO

c-----------------------------------------------------------------------  
c       DEFIENING THE RHS SIDE OF THE HELMHOLTZ EQUATION
c-----------------------------------------------------------------------   
      
      DO 56 I=1,MP1
         F(I,N+1) = 0.
  56  END DO

c----------------------------------------------------------------------- 
c       SOLVING FOR EACH PANEL FOR COEFFIENTS OF FOURIER TRANSFORM
c-----------------------------------------------------------------------  

      DO 58 I=1,MP1
         DO 57 J=1,N
            F(I,J) = 0.
   57    CONTINUE
   58 CONTINUE 
   
c----------------------------------------------------------------------- 
c      SOLVING FOR EACH MODES FOR COEFFIENTS OF IMAG FOURIER TRANSFORM
c-----------------------------------------------------------------------  

      CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     1             BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)

      DO I = 1, 160
         DO J = 1,  511 
            Gi(I,J,K) = F(I,J)*B
         end do
      end do	! SOLVES AS IMAGINARY COMPONENT OF COMPLEX NUMBER
      
   60 end do ! K LOOP END

C######################################################################

c  ----SAVING THE FOURIER COEFFIENTS FOR POWER SPECTRUM ANALYSIS-------
    
      OPEN(UNIT=10, 
     & FILE='/home/ijas/Fortran_IIA/testing_code/res/coef_fourier.dat',
     & STATUS='replace')

      DO I = 1, MP1 
       DO K =1 ,575 !i only need the half
        DO J = 1, NP1
       WRITE(10,'(1F15.7,2E15.7)')theta(I),REAL(G(I,J,K))
     &       ,AIMAG(Gi(I,J,K))
        END DO
       END DO
      END DO 


C          ---------------BACKWARD FOURIER TRANSFORM----------------

      DO I = 1,160
       DO J = 1,511
        DO K = 1, 575
         FREQ_U(K) = G(I,J,K) + Gi(I,J,K)
        END DO         
        CALL CFFTB(575,FREQ_U,WSAVE)
        DO K = 1,575
         U(I,J,K) = FREQ_U(K) ! SAVING THE SCALAR MAGNETIC POTENTIAL
        END DO
       END DO
      END DO
    
      PRINT*, IERROR
      
C          ---------------SAVING THE RESULTS----------------   

      OPEN(UNIT=3, 
     & FILE='/home/ijas/Fortran_IIA/Sample_thwscsp3d/res/ac74.dat',
     & STATUS='replace')

      DO 110 I = 1, MP1 
       DO K =1 ,288 !i only need the half
        DO 109 J = 1, NP1
       WRITE(3,'(3F15.7,4E15.7)')theta(I),r(J),phi(K),U(I,J,K),B_K(I,K)
  109   END DO
       END DO
  110 END DO
   
        END
