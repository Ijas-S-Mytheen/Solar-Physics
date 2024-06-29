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

        INTEGER :: I,J,K,H,X,Z,Q,P
        INTEGER, PARAMETER :: phi_grid = 1151
        INTEGER, PARAMETER :: theta_grid = 320
        INTEGER, PARAMETER :: r_grid = 511
        
        REAL, ALLOCATABLE :: input_data(:) 
        REAL, ALLOCATABLE, DIMENSION(:,:) :: Brrs
        REAL, ALLOCATABLE, DIMENSION(:,:) :: Br_rs
        
        REAL :: WSAVE((4*phi_grid) + 15)
        COMPLEX :: A, B_r_THETA(phi_grid),B_K(theta_grid,phi_grid),B
        character(200)  f1
        REAL :: B_r(phi_grid,theta_grid)
        REAL :: B_r_actual(theta_grid,phi_grid)
        complex :: G(theta_grid,r_grid,phi_grid) 
     &    ,Gi(theta_grid,r_grid,phi_grid)
        COMPLEX ::  FREQ_U(phi_grid)
        COMPLEX :: U(theta_grid,r_grid,phi_grid)
        REAL :: BDTS(r_grid), BDRS(theta_grid), n1, d1, 
     &   n2, d2 , PHI(phi_grid)   
         
        DIMENSION     F(theta_grid,r_grid), BDTF(r_grid),  
     &   R(r_grid), THETA(theta_grid), BDRF(theta_grid), W(13259)


C        DIMENSION OF W CAN BE CALCULATED BY THIS EXPRESSION     
c        LOG_2 = LOG(REAL(r_grid - 1))/LOG(2.0)
c        X = INT(LOG_2) + 1
c        Y = 2**(X+1)
c        Z = (X-2)*Y+X+5*(theta_grid +r_grid - 2)
c        Q = MAX(2*(r_grid - 1),6*(theta_grid - 1))+23
c        P = Z + Q 
c        PRINT*, P    

        ALLOCATE(input_data(theta_grid*phi_grid))
        ALLOCATE(Brrs(phi_grid,theta_grid))
        ALLOCATE(Br_rs(theta_grid,phi_grid))          
     
C-----------------------------------------------------------------------

C    'f1' IS THE PATH OF THE INPUT FILE THAT CONTAINS THE RADIAL 
C     COMPONENT OF MAGNETIC FIELD AT THE LOWER BOUNDARY (Br_rs).
C     MAKE SURE THAT THE LOWER BOUNDARY FUNCTION IS SYMMETRIC ALONG 
C     LONGITUDE.IF NOT YOU CAN CREATE SYMMETRIC USING THE PYTHON CODE   
C     'Doubling_for_fourier.ipynb' IN THE FOLDER.

C-----------------------------------------------------------------------

C          ----------------SETTING UP THE FILES-----------------
c        f1 = '/home/ijas/Fortran_IIA/testing_code/D_B/BrVAR74D.txt' 
        f1='/home/ijas/Fortran_IIA/Sample_thwscsp3d/data/Br_wide_D.txt' 
        OPEN(UNIT=10, FILE=f1) 
        DO i = 1, phi_grid*theta_grid
            READ(10, *) input_data(i)
        END DO
        CLOSE(10)

        DO i = 1, phi_grid	! the file was actually in phi,theta
            DO j = 1, theta_grid	! so we need that file in theta,phi
                Brrs(i, j) = input_data((i - 1) * theta_grid + j)
            END DO
        END DO
        
        DO i = 1, phi_grid
            DO j = 1, theta_grid
                Br_rs(j, i) = Brrs(i, j)
            END DO
        END DO
        PRINT*, '---Br_rs---'
        PRINT*, 'MAX_VAL', MAXVAL(Br_rs), 'MIN_VAL', MINVAL(Br_rs)
        


C          ---------------FORWARD FOURIER TRANSFORM----------------

	A = (1,0)
	B = (0,1)
        call CFFTI(phi_grid,WSAVE)
        DO I = 1, theta_grid
         DO J=1,phi_grid
          B_r_THETA(J) = A * Br_rs(i, j)
         END DO
         !PRINT*, 1,J,B_r_THETA
         call CFFTF(phi_grid,B_r_THETA,WSAVE)
         DO K =1 , phi_grid
          B_K(I,K) = B_r_THETA(K)/phi_grid
         END DO
        end do
        PRINT*, '---REAL B_K---'
        PRINT*,'MAX_VAL',MAXVAL(REAL(B_K)), 'MIN_VAL',MINVAL(REAL(B_K))
        PRINT*, '---AIMAG B_K---'
        PRINT*,'MAX',MAXVAL(AIMAG(B_K)),'MIN',MINVAL(AIMAG(B_K))      
        
C######################################################################

C          ---------------HWSCSP FOR 3D_REAL----------------
 
        DO 1 K =1 , phi_grid ! THIS IS THE FOR EACH MODES
        !print*, K
        
        PI = PIMACH(DUM)
        INTL = 0
        PS = -0.98
        PF = 0.98
        L = phi_grid - 1 
        TS = 1.0471976
        TF = 1.8325958
        M = theta_grid - 1
        MBDCND = 3
        RS = 1.
        RF = 6.
        N = r_grid - 1
        NBDCND = 4
        
c-----------------------------------------------------------------------      
c       DEFINING THE GRID POINTS IN PHI DIRECTION 
c-----------------------------------------------------------------------    
      
        LP1 = L+1
        DP = (PF-PS)/((((LP1+1)/2) - 1))
        DO 100 H=1,LP1
         PHI(H) = FLOAT(H-1)*DP + PS
  100   END DO
  
        ELMBDA=-(2*(1-COS(((2*PI)*(FLOAT(K)-1)/phi_grid))))/(DP**2)
        IDIMF = M+1
        
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

      DO I = 1, theta_grid
         DO J = 1,  r_grid
            G(I,J,K) = F(I,J)*A	
         end do		
      end do		! SOLVES AS REAL COMPONENT OF COMPLEX NUMBER

      
    1 end do ! K LOOP END
    
C######################################################################

C          ---------------HWSCSP FOR 3D_IMAGINARY----------------

        DO 60 K =1 , phi_grid    ! THIS IS THE FOR EACH MODES
        ! print*, K
        
        PI = PIMACH(DUM)
        INTL = 0
        PS = -0.98
        PF = 0.98
        L = phi_grid - 1 
        TS = 1.0471976
        TF = 1.8325958
        M = theta_grid - 1
        MBDCND = 3
        RS = 1.
        RF = 6.
        N = r_grid - 1
        NBDCND = 4
         
c-----------------------------------------------------------------------      
c       DEFINING THE GRID POINTS IN PHI DIRECTION 
c-----------------------------------------------------------------------   
       
        LP1 = L+1
        DP = (PF-PS)/((((LP1+1)/2) - 1))
        DO 50 H=1,LP1
         PHI(H) = FLOAT(H-1)*DP + PS
   50   END DO
  
        ELMBDA=-(2*(1-COS(((2*PI)*(FLOAT(K)-1)/phi_grid))))/(DP**2)
        IDIMF = M + 1
        
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

      DO I = 1, theta_grid
         DO J = 1,  r_grid
            Gi(I,J,K) = F(I,J)*B
         end do
      end do	! SOLVES AS IMAGINARY COMPONENT OF COMPLEX NUMBER
      
   60 end do ! K LOOP END

C######################################################################

c  ----SAVING THE FOURIER COEFFIENTS FOR POWER SPECTRUM ANALYSIS-------
    
      OPEN(UNIT=10, 
     & FILE='coef_fourier.dat', 
     & STATUS='replace')

      DO I = 1, MP1 
      ! print*, I
       DO K =1 ,phi_grid !i only need the half
        DO J = 1, NP1
       WRITE(10,'(1F15.7,2E15.7)')theta(I),REAL(G(I,J,K))
     &       ,AIMAG(Gi(I,J,K))
        END DO
       END DO
      END DO 


C          ---------------BACKWARD FOURIER TRANSFORM----------------

      DO I = 1,theta_grid
       DO J = 1,r_grid
        DO K = 1, phi_grid
         FREQ_U(K) = G(I,J,K) + Gi(I,J,K)
        END DO         
        CALL CFFTB(phi_grid,FREQ_U,WSAVE)
        DO K = 1,phi_grid
         U(I,J,K) = FREQ_U(K) ! SAVING THE SCALAR MAGNETIC POTENTIAL
        END DO
       END DO
      END DO
    
      PRINT*, IERROR
      
C          ---------------SAVING THE RESULTS----------------   

      OPEN(UNIT=3, 
     & FILE='ac74.dat',   ! SAVE THE FILE
     & STATUS='replace')

      DO 110 I = 1, MP1
       DO K =1 ,(LP1+1)/2 !i only need the half
        DO 109 J = 1, NP1
       WRITE(3,'(3F15.7,4E15.7)')theta(I),r(J),phi(K),U(I,J,K),B_K(I,K)
  109   END DO
       END DO
  110 END DO

        DEALLOCATE(input_data)
        DEALLOCATE(Brrs)
        DEALLOCATE(Br_rs) 
   
        END
