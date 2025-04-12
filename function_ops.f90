!==================================================================================================================
! MODULE CONTAINS ROUTINES REPRESENTING OPERATIONS APPLIED TO FUNCTIONS
!
! Author: Jovan Zigic (inherited from Pritpal Matharu)                                   
! McMaster University                                 
! Date: 2024/12/06                                     
!
! CONTAINS:
! (*) initial_guess   - Creates/loads initial condition of system
! (*) vortJY          - Vorticity function used for creating the Jeong-Yoneda IC
! (*) phiIC           - Creates odd symmetric portion of the Jeong-Yoneda IC
! (*) bumpfn          - Creates radial bump function
! (*) vort2velFR      - Transform vorticity to velocity in Fourier space
! (*) vort_test       - Test the FFT and derivative operation
! (*) kinetic_energy  - Calculate the kinetic energy (performed in Fourier space)
! (*) enstrophy       - Calculate the enstrophy (performed in Fourier space)
! (*) palinstrophy    - Calculate the palinstrophy (performed in Fourier space)
! (*) palinstrophyreal- Calculate the palinstrophy (performed in physical space, ***more robust***)
! (*) cal_stream      - Determine the streamfunction from the vorticity in Fourier space
! (*) J_bilinear_form - Determine the convection term in 2D NS from vorticity and streamfunction in Fourier space
! (*) dealiasing      - Dealias field in Fourier space
! (*) inner_product   - Compute inner product between two functions
! (*) mean_vort       - Compute mean of the vorticity solution
! (*) trapzT          - Compute integral over time using trapezoidal rule
! (*) adj_stream      - Determine the adjoint streamfunction
! (*) adj_conv        - Determine the adjoint convection term
! (*) cal_deriv       - Determine x and y derivatives of a 2D field
! (*) SobolevGrad     - Compute the Sobolev gradient
! (*) ProjConstr      - Projection operator for constraint
!==================================================================================================================
MODULE function_ops
  IMPLICIT NONE ! Prevent using implicit typing rules for entire module

  CONTAINS
    !==========================================================================
    ! *** Generate Initial guess, assigns IC to vort0 ***
    ! Input:  mytype - String for type of IC
    !==========================================================================
    SUBROUTINE initial_guess(mytype)
      ! Load variables
      USE global_variables, ONLY: pr, PI, n_nse, dx, local_Ny, local_Nx, local_y_offset, vort0, BS_flag, Lx, Ly, MACH_EPSILON
      ! Load subroutines
      USE fftwfunction, ONLY: fftfwd, fftbwd                    ! FFT routines
      USE data_ops,     ONLY: read_IC, read_NS_Opt, read_BS_Opt ! Functions for reading IC from file
      ! Initialize variables
      CHARACTER(len = *),                INTENT(IN) :: mytype ! String for initial vorticity field
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx) :: w_hat  ! Fourier transform of complex vorticity and RK storage term
      REAL(pr)                                      :: X, Y   ! Scalars for storing x and y coordinates
      INTEGER                                       :: i, j   ! Temporary integers for loops
      REAL(pr)                                      :: noise(n_nse(1))  ! Scalars for random noise

      ! Determine initial condtion selected
      SELECT CASE (mytype)
        ! Initial condition from Jeong-Yoneda
        CASE ("jy")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            ! Set parameters in function call and determine initial vorticity field
            CALL vortJY( X, Y, 13, 2, vort0(i,j) )
          END DO
        END DO

        ! Filter initial condition (to get consistent values after dealiasing)
        CALL fftfwd(vort0,w_hat)  ! Forward Fourier transform
        CALL dealiasing(w_hat)    ! Filter in Fourier space
        CALL fftbwd(w_hat, vort0) ! Backward Fourier transform
        ! Normalize initial vorticity to have norm equal to 1 (dependent on global variable normconstr)
        CALL ProjConstr(vort0)

        ! Load optimal initial condition
        CASE ("optjy")
        ! Determine whether this has been bootstrappted before
        IF (BS_flag) THEN
          ! Call optimal initial condition from previous optimization run
          CALL read_NS_Opt( "jybs", vort0 )
        ELSE
          ! Call optimal initial condition from optimization
          CALL read_NS_Opt( "jy", vort0 )
        END IF

       ! Load optimal initial condition
        CASE ("optcontjy")
        ! Determine whether this has been bootstrappted before
        IF (BS_flag) THEN
          ! Call optimal initial condition from previous optimization run
          CALL read_NS_Opt( "jybs", vort0 )
        ELSE
          ! Call optimal initial condition from optimization
          CALL read_NS_Opt( "jy", vort0 )
        END IF

        ! Load optimal initial condition
        CASE ("jybs")
        ! Determine whether this is the first bootstrap or not
        IF (BS_flag) THEN
          ! Call optimal initial condition from more than one optimization
          CALL read_BS_Opt( "jybs", vort0 )
        ELSE
          ! Call optimal initial condition from single optimization
          CALL read_BS_Opt( "jy", vort0 )
        END IF

        ! Load optimal initial condition
        CASE ("optrand")
        ! Call random initial condition from single optimization
        CALL read_NS_Opt( "rand", vort0 )

        ! Load optimal initial condition
        CASE ("randbs")
        ! Determine whether this is the first bootstrap or not
        IF (BS_flag) THEN
          ! Call optimal initial condition from more than one optimization
          CALL read_BS_Opt( "randbs", vort0 )
        ELSE
          ! Call optimal initial condition from single optimization
          CALL read_BS_Opt( "rand", vort0 )
        END IF

        ! Sine initial condition
        CASE ("sine")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            vort0(i,j) = 2*sin(2.0_pr*PI*X)*sin(2.0_pr*PI*Y)
          END DO
        END DO

        ! Kuramoto-Sivashinsky fixed-scale sinusoidal initial condition
        CASE ("sineKS")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            vort0(i,j) = sin(X + Y) + sin(X) + sin(Y)
          END DO
        END DO

        ! Domain-scaled sinusoidal initial condition
        CASE ("sineL")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            vort0(i,j) = sin(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y) + sin(((2.0_pr*PI)/Lx)*X) + sin(((2.0_pr*PI)/Ly)*Y)
          END DO
        END DO

        ! Domain-scaled multimodal sinusoidal initial conditions of varying magnitudes
        CASE ("mn1")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            vort0(i,j) = ((10.0_pr)**(-1))*( sin(1.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(5.0*((2.0_pr*PI)/Lx)*X) + sin(10.0*((2.0_pr*PI)/Ly)*Y) &
                                          + sin(60.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(50.0*((2.0_pr*PI)/Lx)*X) + sin(30.0*((2.0_pr*PI)/Ly)*Y) &
                                          + sin(80.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(150.0*((2.0_pr*PI)/Lx)*X) + sin(90.0*((2.0_pr*PI)/Ly)*Y) )
          END DO
        END DO

        CASE ("mn5")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            vort0(i,j) = ((10.0_pr)**(-5))*( sin(1.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(5.0*((2.0_pr*PI)/Lx)*X) + sin(10.0*((2.0_pr*PI)/Ly)*Y) &
                                          + sin(60.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(50.0*((2.0_pr*PI)/Lx)*X) + sin(30.0*((2.0_pr*PI)/Ly)*Y) &
                                          + sin(80.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(150.0*((2.0_pr*PI)/Lx)*X) + sin(90.0*((2.0_pr*PI)/Ly)*Y) )
          END DO
        END DO

        CASE ("mn10")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            vort0(i,j) = ((10.0_pr)**(-10))*( sin(1.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(5.0*((2.0_pr*PI)/Lx)*X) + sin(10.0*((2.0_pr*PI)/Ly)*Y) &
                                           + sin(60.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(50.0*((2.0_pr*PI)/Lx)*X) + sin(30.0*((2.0_pr*PI)/Ly)*Y) &
                                           + sin(80.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(150.0*((2.0_pr*PI)/Lx)*X) + sin(90.0*((2.0_pr*PI)/Ly)*Y) )
          END DO
        END DO

        CASE ("mn14")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            vort0(i,j) = ((10.0_pr)**(-14))*( sin(1.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(5.0*((2.0_pr*PI)/Lx)*X) + sin(10.0*((2.0_pr*PI)/Ly)*Y) &
                                           + sin(60.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(50.0*((2.0_pr*PI)/Lx)*X) + sin(30.0*((2.0_pr*PI)/Ly)*Y) &
                                           + sin(80.0*(((2.0_pr*PI)/Lx)*X + ((2.0_pr*PI)/Ly)*Y)) + sin(150.0*((2.0_pr*PI)/Lx)*X) + sin(90.0*((2.0_pr*PI)/Ly)*Y) )
          END DO
        END DO

        ! Random noise initial condition
        CASE ("machepsnoise")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            call random_seed() ! Seed the random number generator
            call random_number(noise) ! Generate uniform random noise in the range [0,1)
            noise = 2.0 * noise - 1.0 ! Generate random noise in the range [-1, 1)
            vort0(i,j) = noise(i)*MACH_EPSILON*100 ! assign value at mach-eps magnitude to solution grid
          END DO
        END DO

        ! Kuramoto-Sivashinsky Gaussian initial condition
        CASE ("gaussianKS")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            vort0(i,j) = exp(-10*( (X - 0.5*(2.0_pr*PI)*10.0_pr)**(2) + (Y - 0.5*(2.0_pr*PI)*10.0_pr)**(2)))
          END DO
        END DO

        ! Taylor-Green initial condition
        CASE ("tg")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            ! Determine initial vorticity field
            vort0(i,j) = 2.0_pr*sin(2.0_pr*PI*10.0_pr*X)*sin(2.0_pr*PI*10.0_pr*Y)
          END DO
        END DO

        ! Filter initial condition (to get consistent values after dealiasing)
        CALL fftfwd(vort0,w_hat)  ! Forward Fourier transform
        CALL dealiasing(w_hat)    ! Filter in Fourier space
        CALL fftbwd(w_hat, vort0) ! Backward Fourier transform
        ! Normalize initial vorticity to have norm equal to 1 (dependent on global variable normconstr)
        CALL ProjConstr(vort0)

        ! Gaussian initial condition
        CASE ("gaussian")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            ! Determine initial vorticity field
            vort0(i,j) = EXP(-100.0_pr*(X-0.5_pr)**2 -100.0_pr*(Y-0.75_pr)**2 ) + &
                         EXP(-100.0_pr*(X-0.5_pr)**2 -100.0_pr*(Y-0.25_pr)**2 )
          END DO
        END DO

        ! Opposite gaussian initial condition
        CASE ("opposite_gaussian")
        DO j=1,local_Ny
          DO i=1,n_nse(1)
            ! x and y coordinates
            X = REAL(i-1,pr)*dx(1)
            Y = REAL(local_y_offset-1+j,pr)*dx(2)
            ! Determine initial vorticity field
            vort0(i,j) = EXP(-100.0_pr*(X-0.5_pr)**2 -100.0_pr*(Y-0.75_pr)**2 ) - &
                         EXP(-100.0_pr*(X-0.5_pr)**2 -100.0_pr*(Y-0.25_pr)**2 )
          END DO
        END DO

        ! Random initial condition from JY, from MATLAB
        CASE ("jyrand")
        ! Call random initial condition generated in MATLAB
        CALL read_IC( "jyrand", n_nse(1), local_Ny, vort0 )
        ! Filter initial condition (to get consistent values after dealiasing)
        CALL fftfwd(vort0,w_hat)  ! Forward Fourier transform
        CALL dealiasing(w_hat)    ! Filter in Fourier space
        CALL fftbwd(w_hat, vort0) ! Backward Fourier transform
        ! Normalize initial vorticity to have norm equal to 1 (dependent on global variable normconstr)
        CALL ProjConstr(vort0)

        ! Random initial condition from JY, from MATLAB
        CASE ("jyrand3")
        ! Call random initial condition generated in MATLAB
        CALL read_IC( mytype, n_nse(1), local_Ny, vort0 )
        ! Filter initial condition (to get consistent values after dealiasing)
        CALL fftfwd(vort0,w_hat)  ! Forward Fourier transform
        CALL dealiasing(w_hat)    ! Filter in Fourier space
        CALL fftbwd(w_hat, vort0) ! Backward Fourier transform
        ! Normalize initial vorticity to have norm equal to 1 (dependent on global variable normconstr)
        CALL ProjConstr(vort0)

        ! Random initial condition, from MATLAB
        CASE ("rand")
        ! Call random initial condition generated in MATLAB
        CALL read_IC( "rand", n_nse(1), local_Ny, vort0 )
        ! Filter initial condition (to get consistent values after dealiasing)
        CALL fftfwd(vort0,w_hat)  ! Forward Fourier transform
        CALL dealiasing(w_hat)    ! Filter in Fourier space
        CALL fftbwd(w_hat, vort0) ! Backward Fourier transform
        ! Normalize initial vorticity to have norm equal to 1 (dependent on global variable normconstr)
        CALL ProjConstr(vort0)
      END SELECT
    END SUBROUTINE initial_guess

    !==========================================================================
    ! *** Function for creating Jeong-Yoneda IC ***
    ! Input:  x - x coordinate
    !         y - y coordinate
    !         n - upper value for bubbles
    !         L - lower value for bubbles
    ! Output: w - Value of vorticity field at (x, y)
    !==========================================================================
    SUBROUTINE vortJY(x, y, n, L, w)
      ! Load variables
      USE global_variables, ONLY: pr
      ! Initialize variables
      REAL(pr), INTENT(IN)  :: x, y ! x and y coordinates, and parameters for creating the bump function
      INTEGER,  INTENT(IN)  :: n, L ! Upper and lower values for number of bubbles
      REAL(pr), INTENT(OUT) :: w    ! Value of the initial condition at coordinates
      REAL(pr)              :: tmp  ! Temporary storage value
      INTEGER               :: i    ! Temporary integer for loops

      ! Initialize initial condition with factor = 1
      CALL phiIC( x, y, 1.0_pr, w )
      ! Create Bourgain-Li bubbles for IC
      DO i = L, n
        ! Determine bubble
        CALL phiIC( x, y, 2.0_pr**i, tmp )
        ! Add to initial condition
        w = w + (REAL(n-i,pr))*tmp/20.0_pr
      END DO
    END SUBROUTINE vortJY

    !==========================================================================
    ! *** Function for creating odd symmetric portion of IC ***
    ! Input:   x  - x coordinate
    !          y  - y coordinate
    !         fac - Scaling factor
    ! Output: phi - Value of bump function
    !==========================================================================
    SUBROUTINE phiIC(x, y, fac, phi)
      ! Load variables
      USE global_variables, ONLY: pr
      ! Initialize variables
      REAL(pr), INTENT(IN)  :: x, y, fac  ! x and y coordinates, and scaling factor for creating the bump function
      REAL(pr), INTENT(OUT) :: phi        ! Value of the bump function at coordinates
      REAL(pr)              :: tmp        ! Temporary storage value
      INTEGER               :: i, j       ! Temporary integers for loops

      ! Start with zero value
      phi = 0.0_pr

      ! Create odd function of initial condition (loop through +/- 1)
      DO i = -1, 1, 2
        DO j = -1, 1, 2
          ! Determine appropritate bump function
          CALL bumpfn( fac*(x - 0.5_pr) + (REAL(i,pr)/4.0_pr), fac*(y - 0.5_pr) + (REAL(j,pr)/4.0_pr), tmp )
          ! Add to odd function
          phi = phi + REAL(i*j,pr)*tmp
        END DO
      END DO
    END SUBROUTINE phiIC

    !==========================================================================
    ! *** Function for creating radial bump function ***
    ! Input:   x  - x coordinate
    !          y  - y coordinate
    ! Output: phi - Value of bump function
    !==========================================================================
    SUBROUTINE bumpfn(x, y, phi)
      ! Load variables
      USE global_variables, ONLY: pr
      ! Initialize variables
      REAL(pr), INTENT(IN)  :: x, y         ! x and y coordinates for creating the bump function
      REAL(pr), INTENT(OUT) :: phi          ! Value of the bump function at coordinates
      REAL(pr)              :: r, rmax, fac ! Storage for the radius, max radius, and a scaling factor

      ! Scaling factor
      fac = 0.05_pr

      ! Max radius (1+1/8 since we have a 1/8 ball in support)
      rmax = SQRT(1.0_pr+(1.0_pr/8.0_pr))/(32.0_pr*fac)

      ! Determine the radius for bump function
      r = (x**2 + y**2)/fac

      ! Radial bump function
      IF (r >= rmax ) THEN
        phi = 0.0_pr ! Compact support in ball
      ELSE
        phi = EXP(-1.0_pr/(1.0_pr - r))
      END IF
    END SUBROUTINE bumpfn

    !==========================================================================
    ! *** Transform vorticity to velocity ***
    ! Input:  fw - vorticity in Fourier space
    ! Output: fU - velocity field in Fourier space
    !==========================================================================
    SUBROUTINE vort2velFR(fw, fU)
      ! Load variables
      USE global_variables, ONLY: pr, MACH_EPSILON, n_nse, K1, K2, ksq, local_Nx, local_x_offset
      ! Load subroutines
      USE fftwfunction, ONLY: fftfwd, fftbwd  ! FFT routines
      ! Initialize variables
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx),   INTENT(IN)  :: fw   ! Fourier transform of vorticity
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx,2), INTENT(OUT) :: fU   ! Fourier transform of velocity
      COMPLEX(pr), DIMENSION(:,:),   INTENT(IN)  :: fw   ! Fourier transform of vorticity
      COMPLEX(pr), DIMENSION(:,:,:), INTENT(OUT) :: fU   ! Fourier transform of velocity ! 2DNS
      !COMPLEX(pr), DIMENSION(:,:), INTENT(OUT) :: fU   ! Fourier transform of velocity ! 2DKS
      INTEGER :: i1, i2                                                    ! Integers for looping through values

      ! Compute the velocity field by using the Laplace operator and grad^perp (via the streamfunction)
      DO i2 = 1,local_Nx
        DO i1 = 1, n_nse(2)
          ! Compute each component of the velocity field
          IF (ksq(i1, i2) > MACH_EPSILON) THEN
            fU(i1,i2,1) = CMPLX(0.0_pr,K2(i1),pr)*fw(i1,i2)/ksq(i1, i2)
            fU(i1,i2,2) = -CMPLX(0.0_pr,K1(i2+local_x_offset),pr)*fw(i1,i2)/ksq(i1, i2)
          ELSE
            fU(i1,i2,1) = CMPLX(0.0_pr,0.0_pr,pr)
            fU(i1,i2,2) = CMPLX(0.0_pr,0.0_pr,pr)
          END IF
        END DO
      END DO
    END SUBROUTINE vort2velFR

    !==========================================================================
    ! *** Testing the transform and derivative operations ***
    ! Input: W - vorticity in physical space
    !==========================================================================
    SUBROUTINE vort_test(W)
      ! Load variables
      USE global_variables
      ! Load subroutines
      USE fftwfunction, ONLY: fftfwd, fftbwd  ! FFT routines
      USE data_ops, ONLY: save_NS_vorticity   ! Functions for saving files
      ! Initialize variables
!      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny), INTENT(IN) :: W             ! Vorticity field
      REAL(pr), DIMENSION(:,:), INTENT(IN) :: W             ! Vorticity field
      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny)             :: Wx, Wy, vortW ! Test values of vorticity field
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)          :: faux0, faux1  ! Storage for Fourier space
      INTEGER :: i1, i2

      ! Save original vorticity
      CALL save_NS_vorticity(W, Time_iter, "Original")
      ! Transform vorticity into Fourier space
      CALL fftfwd(W, faux0)
      ! Take Derivative with respect to x
      DO i2=1,local_Nx
        DO i1=1,n_nse(2)
          faux1(i1,i2) = CMPLX(0.0_pr,K1(i2+local_x_offset),pr)*faux0(i1,i2)
        END DO
      END DO
      ! Transform back to physical space
      CALL fftbwd(faux1, Wx)
      ! Save derivative with respect to x
      CALL save_NS_vorticity(Wx, Time_iter, "Testxderiv")

      ! Take derivative with respect to y
      DO i2=1,local_Nx
        DO i1=1,n_nse(2)
          faux1(i1,i2) = CMPLX(0.0_pr,K2(i1),pr)*faux0(i1,i2)
        END DO
      END DO
      ! Transform back to physical space
      CALL fftbwd(faux1, Wy)
      ! Save derivative with respect to y
      CALL save_NS_vorticity(Wy, Time_iter, "Testyderiv")
      ! Transform original back to physical space, for testing
      CALL fftbwd(faux0, vortW)
      ! Save vorticity (should be same as original)
      CALL save_NS_vorticity(vortW, Time_iter, "TestFFT")
    END SUBROUTINE vort_test

    !==========================================================================
    ! *** Calculate kinetic energy from function in Fourier space ***
    ! Input:     fu    - velocity in Fourier space
    ! Output: kin_ener - kinetic energy (local)
    !==========================================================================
    FUNCTION kinetic_energy(fu) RESULT (kin_ener)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, dV, local_Nx
      ! Initialize variables
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx,2), INTENT(IN) :: fu       ! Velocity in Fourier space
      COMPLEX(pr), DIMENSION(:,:,:), INTENT(IN) :: fu       ! Velocity in Fourier space
      REAL(pr)                                                    :: kin_ener ! Scalar value of kinetic energy
      INTEGER                                                     ::  i1, i2  ! Temporary integers for loops

      ! Start with zero value
      kin_ener = 0.0_pr
      ! Periodic domain, so guassian quadrature is equivalent to taking a double sum over the domain
      DO i2 = 1,local_Nx
        DO i1 = 1,n_nse(2)
          ! Due to conjugate symmetry, we only have half the Fourier modes, so we double the sum
          ! Hence, we remove the factor of 0.5_pr from the summation (except for along the x zeroth mode, R2C)
          IF (i2 == 1) THEN
            kin_ener = kin_ener + 0.5_pr*ABS(fu(i1,i2,1))**2*dV
            kin_ener = kin_ener + 0.5_pr*ABS(fu(i1,i2,2))**2*dV
          ELSE
            kin_ener = kin_ener + ABS(fu(i1,i2,1))**2*dV
            kin_ener = kin_ener + ABS(fu(i1,i2,2))**2*dV
          END IF
!          kin_ener = kin_ener + 0.5_pr*ABS(fu(i1,i2,1))**2*dV
!          kin_ener = kin_ener + 0.5_pr*ABS(fu(i1,i2,2))**2*dV
        END DO
      END DO
      ! Normalize kinetic energy
      ! FFTW computes an unnormalized DFT and we are not taking the inverse transform
      kin_ener = kin_ener/(PRODUCT(REAL(n_nse,pr)))
    END FUNCTION kinetic_energy

    !==========================================================================
    ! *** Calculate enstrophy from function in Fourier space ***
    ! Input:   fw - vorticity in Fourier space
    ! Output: ens - enstrophy (local)
    !==========================================================================
    FUNCTION enstrophy(fw) RESULT (ens)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, dV, local_Nx, ksq, Kcut
      ! Initialize variables
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx), INTENT(IN) :: fw     ! Vorticity in Fourier space
      COMPLEX(pr), DIMENSION(:,:), INTENT(IN) :: fw     ! Vorticity in Fourier space
      REAL(pr)                                                  :: ens    ! Scalar value of enstrophy
      REAL(pr)                                                  :: mode   ! Temporary value for Fourier mode
      INTEGER                                                   :: i1, i2 ! Temporary integers for loops

      ! Start with zero value
      ens = 0.0_pr
      ! Periodic domain, so gaussian quadrature is equivalent to taking a double sum over the domain
      DO i2 = 1,local_Nx
        DO i1 = 1,n_nse(2)
          ! Determine Fourier mode for dealiasing
          mode = SQRT( ksq(i1, i2) )
          ! Due to conjugate symmetry, we only have half the Fourier modes, so we double the sum
          ! Hence, we remove the factor of 0.5_pr from the summation (except for along the x zeroth mode, R2C)
          IF (i2 == 1) THEN
            ens = ens + 0.5_pr*ABS(fw(i1,i2)*EXP(-36.0_pr*(mode/Kcut)**36))**2*dV
          ELSE
            ens = ens + ABS(fw(i1,i2)*EXP(-36.0_pr*(mode/Kcut)**36))**2*dV
          END IF

        END DO
      END DO
      ! Normalize enstrophy
      ! FFTW computes an unnormalized DFT and we are not taking the inverse transform
      ens = ens/(PRODUCT(REAL(n_nse,pr)))
    END FUNCTION enstrophy

    !==========================================================================
    ! *** Calculate Palinstrophy from function in Fourier space ***
    ! Input:    fw   - vorticity in Fourier space
    ! Output: palins - palinstrophy (local)
    !==========================================================================
    FUNCTION palinstrophy(fw) RESULT (palins)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, dV, local_Nx, local_x_offset, K1, K2, ksq, Kcut
      ! Initialize variables
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx), INTENT(IN) :: fw     ! Vorticity in Fourier space
      COMPLEX(pr), DIMENSION(:,:), INTENT(IN) :: fw     ! Vorticity in Fourier space
      REAL(pr)                                                  :: palins ! Scalar value of palinstrophy
      REAL(pr)                                                  :: mode   ! Temporary value for Fourier mode
      INTEGER                                                   :: i1, i2 ! Temporary integers for loops

      ! Start with zero value
      palins = 0.0_pr
      ! Periodic domain, so gaussian quadrature is equivalent to taking a double sum over the domain
      DO i2=1,local_Nx
        DO i1=1,n_nse(2)
!        ! Determine Fourier mode for dealiasing
!        mode = SQRT( ksq(i1, i2) )
        ! Due to conjugate symmetry, we only have half the Fourier modes, so we double the sum
        ! Hence, we remove the factor of 0.5_pr from the summation (except for along the x zeroth mode, R2C)
        IF (i2 == 1) THEN
          palins = palins + 0.5_pr*ABS( CMPLX(0.0_pr,K1(i2+local_x_offset),pr)*fw(i1,i2) )**2*dV
          palins = palins + 0.5_pr*ABS( CMPLX(0.0_pr,K2(i1),pr)*fw(i1,i2) )**2*dV
        ELSE
          palins = palins + ABS( CMPLX(0.0_pr,K1(i2+local_x_offset),pr)*fw(i1,i2) )**2*dV
          palins = palins + ABS( CMPLX(0.0_pr,K2(i1),pr)*fw(i1,i2) )**2*dV
        END IF
        END DO
      END DO
      ! Normalize palinstrophy
      ! FFTW computes an unnormalized DFT and we are not taking the inverse transform
      palins = palins/(PRODUCT(REAL(n_nse,pr)))
    END FUNCTION palinstrophy

    !==========================================================================
    ! *** Calculate Palinstrophy from function in physical space ***
    ! !*Note: This version is more robust because this forces the summation to be
    ! of numbers that are relatively the same magnitude (where as in Fourier
    ! space/ using Parseval's identity is more susecptiable to round-off errors
    !, especially as number of discretization points increase)*!
    ! Input:    fw   - vorticity in Fourier space
    ! Output: palins - palinstrophy (local)
    !==========================================================================
    FUNCTION palinstrophyreal(fw) RESULT (palins)
      ! Load variables
      USE global_variables!, ONLY: pr, n_nse, dV, local_Ny, local_y_offset, K1, K2, ksq, Kcut
      ! Initialize variables
      COMPLEX(pr), DIMENSION(:,:), INTENT(IN) :: fw     ! Vorticity in Fourier space
      REAL(pr)                                                  :: palins ! Scalar value of enstrophy
      REAL(pr),    DIMENSION(1:n_nse(1),1:local_Ny)  :: wx, wy     ! Vorticity, partial derivatives, and velocity components from the forward solution

      ! Partial derivatives of vorticity solution (explicitly treated)
      CALL cal_deriv(fw, wx, wy)

      ! Calculate Palinstrophy
      palins = 0.5_pr*( inner_product(wx, wx, "L2") + inner_product(wy, wy, "L2") )

    END FUNCTION palinstrophyreal

    !==========================================================================
    ! *** Calculate stream function from vorticity ***
    ! Input:  faux - vorticity in Fouier space
    ! Output: fpsi - stream function in Fouier space
    !==========================================================================
    SUBROUTINE cal_stream(faux, fpsi)
      ! Load variables
      USE global_variables, ONLY: pr, MACH_EPSILON, n_nse, local_Nx, ksq
      ! Load subroutines
      USE fftwfunction, ONLY: fftfwd, fftbwd  ! FFT routines
      ! Initialize variables
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx), INTENT(IN)  :: faux   ! Vorticity field in Fourier space
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx), INTENT(OUT) :: fpsi   ! Streamfunction in Fourier space
      COMPLEX(pr), DIMENSION(:,:), INTENT(IN)  :: faux   ! Vorticity field in Fourier space
      COMPLEX(pr), DIMENSION(:,:), INTENT(OUT) :: fpsi   ! Streamfunction in Fourier space
      INTEGER                                                    :: i1, i2 ! Temporary integers for loops

      ! Loop through wavenumbers to determine streamfunction
      DO i2 = 1,local_Nx
        DO i1 = 1,n_nse(2)
          ! Determine stream function (negative cancels from Poisson operator)
          IF (ksq(i1, i2) > MACH_EPSILON) THEN
            fpsi(i1,i2) = faux(i1,i2)/ksq(i1, i2)
          ELSE
            fpsi(i1,i2) = 0.0_pr
          END IF
        END DO
      END DO
    END SUBROUTINE cal_stream

    !==========================================================================
    ! *** Calculate negative bilinear form J(f,g) = -det(f_x f_y; g_x g_y) ***
    ! Convection term in 2D NS
    ! Input:  fhat - vorticity in Fouier space
    !         ghat - stream function in Fouier space
    ! Output: Jhat - Convection term in Fouier space
    !==========================================================================
    SUBROUTINE J_bilinear_form(fhat, ghat, Jhat)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Nx, local_Ny, local_x_offset, K1, K2
      ! Load subroutines
      USE fftwfunction, ONLY: fftfwd, fftbwd  ! FFT routines
      ! Initialize variables
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx), INTENT(IN)  :: fhat, ghat     ! Vorticity and streamfunction in Fourier space
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx), INTENT(OUT) :: Jhat           ! Convection term in Fourier space
      COMPLEX(pr), DIMENSION(:,:), INTENT(IN)       :: fhat, ghat                     ! Vorticity and streamfunction in Fourier space
      COMPLEX(pr), DIMENSION(:,:), INTENT(OUT)      :: Jhat                           ! Convection term in Fourier space
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx) :: fhat_x, fhat_y, ghat_x, ghat_y ! Temporary matrices to store complex Fourier space values
      REAL(pr), DIMENSION(1:n_nse(1), 1:local_Ny)   :: J                              ! Temporary matrix to store Jacobian in physical space
      REAL(pr), DIMENSION(1:n_nse(1), 1:local_Ny)   :: fx, fy, gx, gy                 ! Temporary matrices to store derivatives in physical space
      INTEGER                                       :: i1, i2                         ! Temporary integers for loops

      ! Compute x derivative of vorticity
      DO i2=1,local_Nx
        DO i1=1,n_nse(2)
          fhat_x(i1,i2) = CMPLX(0.0_pr,K1(i2+local_x_offset),pr)*fhat(i1,i2) ! Derivative wrt x of vorticity
          fhat_y(i1,i2) = CMPLX(0.0_pr,K2(i1),pr)*fhat(i1,i2)                ! Derivative wrt y of vorticity
          ghat_x(i1,i2) = CMPLX(0.0_pr,K1(i2+local_x_offset),pr)*ghat(i1,i2) ! Derivative wrt x of streamfunction
          ghat_y(i1,i2) = CMPLX(0.0_pr,K2(i1),pr)*ghat(i1,i2)                ! Derivative wrt y of streamfunction
        END DO
      END DO
      ! Transform to physical space
      ! *** This can be improved upon, using the FFTW Advanced interface to group all transforms into one (rather difficult task, as of 2020)
      CALL fftbwd(fhat_x, fx) ! Derivative wrt x of vorticity transform
      CALL fftbwd(fhat_y, fy) ! Derivative wrt y of vorticity transform
      CALL fftbwd(ghat_x, gx) ! Derivative wrt x of streamfunction transform
      CALL fftbwd(ghat_y, gy) ! Derivative wrt y of streamfunction transform

      ! Compute Jacobian, convective derivative (u,v).grad(w) and obtain correct sign for the Jacobian
      !J = -(fx*gy - fy*gx) ! 2DNS
      J = -(1.0_pr/2.0_pr)*(fx*fx + fy*fy) ! 2DKS
      !J = 0 ! Linearized equation
      ! Compute Fourier transform of Jacobian
      CALL fftfwd(J, Jhat)
      ! Dealias
      CALL dealiasing(Jhat)
    END SUBROUTINE J_bilinear_form

    !==========================================================================
    ! *** Perform dealiasing ***
    ! Input and Output:  ff - vorticity in Fouier space
    !==========================================================================
    SUBROUTINE dealiasing(ff)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Nx, ksq, Kcut, MACH_EPSILON
      ! Initialize variables
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx), INTENT(INOUT) :: ff     ! Vorticity in Fourier space
      COMPLEX(pr), DIMENSION(:,:), INTENT(INOUT) :: ff     ! Vorticity in Fourier space
      REAL(pr)                                                     :: mode   ! Temporary value for Fourier mode
      INTEGER                                                      :: i1, i2 ! Temporary integers for loops

      ! Loop through wavenumbers to dealias
      DO i2 = 1,local_Nx
        DO i1 = 1, n_nse(2)
          ! Ensure that we have mean zero
          IF (ksq(i1, i2) > MACH_EPSILON) THEN
            ! Determine the values of the Fourier modes
            mode = SQRT( ksq(i1, i2) )
            ! Use a Gaussian spectral filter with cutoff for dealiasing
            ff(i1,i2) = ff(i1,i2)*EXP(-36.0_pr*(mode/Kcut)**36) ! 2DNS
            !ff(i1,i2) = ff(i1,i2) ! 2DKS
          ELSE
            ff(i1,i2) = 0.0_pr
          END IF
        END DO
      END DO
    END SUBROUTINE dealiasing

    !==========================================================================
    ! *** Calculate the inner product between two functions (note recursive)***
    ! Input:      f    - function in physical space
    !             g    - function in physical space
    !          mytype  - function space to perform inner product
    ! Output: inn_prod - inner product result
    !==========================================================================
    RECURSIVE FUNCTION inner_product(f,g,mytype) RESULT (inn_prod)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Nx, local_Ny, local_x_offset, K1, K2, dV,ksq, ell, Statinfo
      ! Load subroutines
      USE fftwfunction, ONLY: fftfwd, fftbwd  ! FFT routines
      USE mpi                                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
!      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny), INTENT(IN) :: f,g     ! Functions for inner product in physical space
      REAL(pr), DIMENSION(:,:), INTENT(IN)           :: f,g            ! Functions for inner product in physical space
      CHARACTER(len=*),         INTENT(IN)           :: mytype         ! String for function space
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)  :: fhat           ! Temporary matrix to store complex Fourier space values
      REAL(pr),    DIMENSION(1:n_nse(1), 1:local_Ny) :: ftmp           ! Temporary matrix to store derivatives in physical space
      REAL(pr)                                       :: inn_prod       ! Scalar value of inner product
      REAL(pr)                                       :: local_inn_prod ! Temporary value of local results
      INTEGER                                        :: i1, i2         ! Temporary integers for loops

      SELECT CASE (mytype) ! Function space for inner product
      CASE ("L2") ! L^2 inner product
        local_inn_prod = 0.0_pr ! Start with zero value
        ! Periodic domain, so guassian quadrature is equivalent to taking a double sum over the domain
        DO i2=1,local_Ny
          DO i1=1,n_nse(1)
            local_inn_prod = local_inn_prod + f(i1,i2)*g(i1,i2)*dV ! Integration of product of functions over entire domain
          END DO
        END DO
        
        CALL MPI_ALLREDUCE(local_inn_prod, inn_prod, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, Statinfo) ! Sum over all processor

      CASE ("H1") ! H^1 inner product
        CALL fftfwd(f, fhat) ! Transform to Fourier space to compute derivatives
        ! Periodic domain, so we can integrate by parts and compute by one minus Laplacian on one function
        DO i2=1,local_Nx
          DO i1=1,n_nse(2)
            fhat(i1,i2) = (1.0_pr + (ell**2)*ksq(i1, i2))*fhat(i1,i2)
            !fhat(i1,i2) = (1.0_pr + (ell**2)*ksq(i1, i2))*ABS(fhat(i1,i2))**2
          END DO
        END DO
        
        CALL fftbwd(fhat, ftmp) ! Transform back to physical space, to compute the L2 inner product
        inn_prod = inner_product(ftmp, ftmp, "L2") ! Now, compute over L2 inner product    
      
      CASE ("H2") ! H^2 inner product
        CALL fftfwd(f, fhat) ! Transform to Fourier space to compute derivatives
        ! Periodic domain, so we can integrate by parts and compute by one minus Laplacian on one function
        DO i2=1,local_Nx
          DO i1=1,n_nse(2)
            fhat(i1,i2) = ((1.0_pr + (ell**2)*ksq(i1, i2))**2)*fhat(i1,i2) 
            !fhat(i1,i2) = ((1.0_pr + (ell**2)*ksq(i1, i2))**2)*ABS(fhat(i1,i2))**2
          END DO
        END DO
        
        CALL fftbwd(fhat, ftmp)  ! Transform back to physical space, to compute the L2 inner product
        inn_prod = inner_product(ftmp, ftmp, "L2") ! Now, compute over L2 inner product

      CASE ("Hn1") ! H^(-1) inner product
        CALL fftfwd(f, fhat) ! Transform to Fourier space to compute derivatives
        ! Periodic domain, so we can integrate by parts and compute by one minus Laplacian on one function
        DO i2=1,local_Nx
          DO i1=1,n_nse(2)
            fhat(i1,i2) = fhat(i1,i2)/(1.0_pr + (ell**2)*ksq(i1, i2))
            !fhat(i1,i2) = (ABS(fhat(i1,i2))**2)/(1.0_pr + (ell**2)*ksq(i1, i2))
          END DO
        END DO
        
        CALL fftbwd(fhat, ftmp) ! Transform back to physical space, to compute the L2 inner product
        inn_prod = inner_product(ftmp, ftmp, "L2") ! Now, compute over L2 inner product

      CASE ("H1semi") ! H^1 semi-norm inner product
        CALL fftfwd(f, fhat) ! Transform to Fourier space to compute derivatives

        ! Periodic domain, so we can integrate by parts and compute by negative Laplacian on one function
        DO i2=1,local_Nx
          DO i1=1,n_nse(2)
            fhat(i1,i2) = ksq(i1, i2)*fhat(i1,i2)
          END DO
        END DO

        CALL fftbwd(fhat, ftmp) ! Transform back to physical space, to compute the L2 inner product
        inn_prod = inner_product(ftmp, ftmp, "L2") ! Now, compute over L2 inner product
      END SELECT

    END FUNCTION inner_product

    !==========================================================================
    ! *** Calculate the mean of the function ***
    ! Input:      f    - function in physical space
    ! Output: mean_val - mean result
    !==========================================================================
    FUNCTION mean_vort(f) RESULT (mean_val)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Nx, local_Ny, local_x_offset, dV, Statinfo
      ! Load subroutines
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
      REAL(pr), DIMENSION(:,:), INTENT(IN) :: f            ! Function for to calculate the mean of in physical space
      REAL(pr)                             :: mean_val     ! Scalar value of mean
      REAL(pr)                             :: local_mean   ! Temporary value of local results
      INTEGER                              :: i1, i2       ! Temporary integers for loops

      ! Start with zero value
      local_mean = 0.0_pr
      ! Periodic domain, so guassian quadrature is equivalent to taking a double sum over the domain
      DO i2=1,local_Ny
        DO i1=1,n_nse(1)
          ! Integration of product of functions over entire domain
          local_mean = local_mean + f(i1,i2)*dV
        END DO
      END DO
      ! Sum over all processor (only on rank 0 needs a copy)
      CALL MPI_REDUCE(local_mean, mean_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, Statinfo)
    END FUNCTION mean_vort

    !==========================================================================
    ! *** Integration over time domain **
    ! Input:  f - Intergrand, over time
    ! Output: r - Integral result
    !==========================================================================
    FUNCTION trapzT(f) RESULT(r)
      ! Load variables
      USE global_variables, ONLY: pr, t_vec, numel_T, dt, dtnp, np, rank, Statinfo
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
      REAL(pr), DIMENSION(:), INTENT(IN) :: f          ! Intergrand, to integrate over time
      REAL(pr)                           :: r          ! Resulting integral
      REAL(pr)                           :: r_local    ! Local integral value
      INTEGER                            :: sval, tval ! Integer position for local processors

      ! Starting value for local processor
      sval = 1+rank*dtnp
      ! Terminal value for local processor
      IF (rank == np-1) THEN
        tval = numel_T ! Last processor will pickup the slack, if the number doesn't perfectly divide
      ELSE
        tval = (rank+1)*dtnp
      END IF

      ! Integrate over time, using the trapezoidal rule using each processor compute a portion
      ! Factor of a half comes from trapezoidal rule, not cost functional
      r_local = 0.5_pr*SUM( ( dt )*( f(sval:tval) + f(sval+1:tval+1) ) )

      CALL MPI_ALLREDUCE(r_local, r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, Statinfo)

    END FUNCTION trapzT

    !==========================================================================
    ! *** Calculate the streamfunction for the adjoint system ***
    ! Input:     z - Adjoint vorticity field, in physical space
    !           fx - Derivative of vorticity, wrt x, in physical space
    !           fy - Derivative of vorticity, wrt y, in physical space
    ! Output: ghat - Adjoint streamfunction, in Fourier space
    !==========================================================================
    SUBROUTINE adj_stream(z, fx, fy, ghat)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Nx, local_Ny, MACH_EPSILON, ksq, local_x_offset, K1, K2
      ! Load subroutines
      USE fftwfunction, ONLY: fftfwd, fftbwd  ! FFT routines
      ! Initialize variables
      REAL(pr),    DIMENSION(:, :), INTENT(IN)  :: z                             ! Adjoint vorticity field in physical space
      REAL(pr),    DIMENSION(:, :), INTENT(IN)  :: fx, fy                        ! Derivative of vorticity in physical space
      COMPLEX(pr), DIMENSION(:, :), INTENT(OUT) :: ghat                          ! Adjoint streamfunction field in Fourier space
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)  :: nonlinx_hat, nonliny_hat ! Temporary matrices to store Fourier space values
      REAL(pr),    DIMENSION(1:n_nse(1), 1:local_Ny) :: nonlinx, nonliny         ! Temporary matrices for storage in physical space
      COMPLEX(pr)                                    :: RHS                      ! Temporary for storing RHS of equation
      INTEGER                                        :: i1, i2                   ! Temporary integers for loops

      ! Compute x component of nonlinear term
      nonlinx = fx*z
      ! Compute Fourier transform of nonlinear component
      CALL fftfwd(nonlinx, nonlinx_hat)
      ! Dealias (better to do it individually)
      CALL dealiasing(nonlinx_hat)

      ! Compute y component of nonlinear term
      nonliny = fy*z
      ! Compute Fourier transform of nonlinear component
      CALL fftfwd(nonliny, nonliny_hat)
      ! Dealias (better to do it individually)
      CALL dealiasing(nonliny_hat)

      ! Loop through wavenumbers to determine adjoint streamfunction
      DO i2 = 1,local_Nx
        DO i1 = 1,n_nse(2)
          ! RHS for Poisson equation
          RHS = CMPLX(0.0_pr,K2(i1),pr)*nonlinx_hat(i1,i2) - CMPLX(0.0_pr,K1(i2+local_x_offset),pr)*nonliny_hat(i1,i2)
          ! Determine adjoint stream function from Poisson equation
          IF (ksq(i1, i2) > MACH_EPSILON) THEN
            ghat(i1,i2) = -RHS/ksq(i1, i2)
          ELSE
            ghat(i1,i2) = 0.0_pr
          END IF
        END DO
      END DO

    END SUBROUTINE adj_stream

    !==========================================================================
    ! *** Convection term for the adjoint system ***
    ! Input:     z - Adjoint vorticity field, in physical space
    !            u - Derivative of vorticity, wrt y, in physical space
    !            v - Derivative of vorticity, wrt x, in physical space
    ! Output: ghat - Adjoint streamfunction, in Fourier space
    !==========================================================================
    SUBROUTINE adj_conv(z_hat, u, v, ghat)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Ny
      ! Load subroutines
      USE fftwfunction, ONLY: fftfwd, fftbwd  ! FFT routines
      ! Initialize variables
      COMPLEX(pr), DIMENSION(:, :), INTENT(IN)    :: z_hat  ! Adjoint vorticity field
      REAL(pr),    DIMENSION(:, :), INTENT(IN)    :: u, v   ! Velocity components
      COMPLEX(pr), DIMENSION(:, :), INTENT(OUT)   :: ghat   ! Adjoint convection field
      REAL(pr), DIMENSION(1:n_nse(1), 1:local_Ny) :: nonlin ! Nonlinear term
      REAL(pr), DIMENSION(1:n_nse(1), 1:local_Ny) :: zx, zy ! Derivatives of the adjoint field

      ! Compute the derivatives of the adjoint field
      CALL cal_deriv(z_hat, zx, zy)
      ! Calculate the nonlinear term
      nonlin = u*zx - v*zy
      ! Compute Fourier transform of nonlinear term
      CALL fftfwd(nonlin, ghat)
      ! Dealias
      CALL dealiasing(ghat)
    END SUBROUTINE adj_conv

    !==========================================================================
    ! *** Calculate derivatives in x and y ***
    ! Input:   fhat    - 2D field, in Fourier space
    ! Output:    fx    - Derivative of function wrt x component, in physical space
    !            fy    - Derivative of function wrt y component, in physical space
    !==========================================================================
    SUBROUTINE cal_deriv(fhat, fx, fy)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Nx, local_x_offset, K1, K2
      ! Load subroutines
      USE fftwfunction, ONLY: fftfwd, fftbwd  ! FFT routines
      ! Initialize variables
      COMPLEX(pr), DIMENSION(:,:), INTENT(IN)        :: fhat           ! Input field in Fourier space
      REAL(pr),    DIMENSION(:,:), INTENT(OUT)       :: fx, fy         ! Output derivatives in x and y in physical space
      COMPLEX(pr), DIMENSION(1:n_nse(2), 1:local_Nx) :: fhat_x, fhat_y ! Temporary matrices to store Fourier space values
      INTEGER                                        :: i1, i2         ! Temporary integers for loops

      ! Compute x and y derivative of vorticity
      DO i2=1,local_Nx
        DO i1=1,n_nse(2)
          fhat_x(i1,i2) = CMPLX(0.0_pr,K1(i2+local_x_offset),pr)*fhat(i1,i2)
          fhat_y(i1,i2) = CMPLX(0.0_pr,K2(i1),pr)*fhat(i1,i2)
        END DO
      END DO
      ! Transform to physical space
      CALL fftbwd(fhat_x, fx)
      CALL fftbwd(fhat_y, fy)
    END SUBROUTINE cal_deriv

    !==========================================================================
    ! *** Calculate Sobolev Gradient from L2 Gradient ***
    ! Input:   fhat    - L2 gradient, in Fourier space
    ! Output:     g    - Sobolev gradient, in physical space
    !==========================================================================
    SUBROUTINE SobolevGrad(fhat, g)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Nx, local_x_offset, ksq, ell
      ! Load subroutines
      USE fftwfunction, ONLY: fftfwd, fftbwd  ! FFT routines
      ! Initialize variables
      COMPLEX(pr), DIMENSION(:,:), INTENT(IN)        :: fhat           ! Input field in Fourier space
      REAL(pr),    DIMENSION(:,:), INTENT(OUT)       :: g              ! Output derivatives in x and y in physical space
      COMPLEX(pr), DIMENSION(1:n_nse(2), 1:local_Nx) :: ghat           ! Temporary matrix to store Fourier space values
      INTEGER                                        :: i1, i2         ! Temporary integers for loops

      ! Compute Sobolev gradient
      DO i2=1,local_Nx
        DO i1=1,n_nse(2)
          ghat(i1,i2) = fhat(i1,i2)/(1.0_pr + (ell**2)*ksq(i1, i2))
        END DO
      END DO
      ! Transform to physical space
      CALL fftbwd(ghat, g)
    END SUBROUTINE SobolevGrad

    !==========================================================================
    ! *** Projection operator to enforce norm constraint ***
    ! Input and Output:  ff - vorticity in physical space
    !==========================================================================
    SUBROUTINE ProjConstr(ff)
      ! Load variables
      USE global_variables, ONLY: pr, normconstr, PalinIV
      ! Initialize variables
      REAL(pr), DIMENSION(:,:), INTENT(INOUT) :: ff     ! Vorticity in Fourier space


      ! Function space
      SELECT CASE (normconstr)
      ! H1 semi-norm
      CASE ("H1semi")
        ! Normalize, so norm is equal to PalinIV to enforce constraint on initial Palinstrophy
        ff = SQRT(PalinIV)*ff/SQRT(inner_product(ff, ff, normconstr))
      CASE DEFAULT
        ! Normalize value to 1
        ff = ff/SQRT(inner_product(ff, ff, normconstr))
      END SELECT


    END SUBROUTINE ProjConstr


END MODULE
