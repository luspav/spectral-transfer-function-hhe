!***********************************************************************
! Spectral Radiative Transfer Function for H and He in the IGM
!***********************************************************************
!
! Computes the frequency-dependent spectral transfer function F(nu, z)
! for radiation propagating through the intergalactic medium (IGM)
! containing hydrogen (HI) and helium (HeI, HeII) during cosmic
! reionization. The transfer function accounts for:
!   - Photoionization continuum opacity (HI, HeI, HeII)
!   - Lyman-series line opacity (HI, HeI, HeII)
!   - Thomson scattering by free electrons
!
! References:
!   - Bosman, Becker & Inoue (2022) — neutral hydrogen fraction fits
!   - Kageura & Bosman (2025)       — extended xHI approximation
!   - Glazer & Bosman (2018)        — alternative xHI approximation
!
! Compilation:
!   macOS:   gfortran -isysroot $(xcrun --sdk macosx --show-sdk-path) -O3 -o rtf ./rtf_H-He.f90
!   Linux:   gfortran -O3 -o rtf rtf_H-He.f90
!   Windows: gfortran -O3 -o rtf.exe rtf_H-He.f90
!
! Usage:
!   ./rtf          (or: time ./rtf)
!
! Input files (must be in the working directory):
!   HOS.dat   — HI Lyman-series oscillator strengths (index, wavelength [A], f)
!   HeIOS.dat — HeI line oscillator strengths          (index, wavelength [A], f)
!
! Output files:
!   rtf_H-He_1_z_<z>.dat — two-column tables: frequency [THz] vs F(nu)
!***********************************************************************

PROGRAM main
  IMPLICIT NONE

  ! Double precision kind parameter
  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  ! ---- Grid dimensions ----
  INTEGER, PARAMETER :: n_xf     = 32000   ! Number of frequency grid points
  INTEGER, PARAMETER :: n_z      = 9       ! Number of redshift nodes
  INTEGER, PARAMETER :: n_ly     = 39      ! Number of HI Lyman-series lines
  INTEGER, PARAMETER :: n_lyHeI  = 10      ! Number of HeI lines

  ! ---- Model switch ----
  ! klxHI selects the xHI(z) approximation for z > 5.15:
  !   1 — Bosman+22 + Kageura+25  (KBBI)
  !   2 — Bosman+22 + Glazer+18   (GBI)
  INTEGER, PARAMETER :: klxHI = 1

  ! =====================================================================
  ! Physical constants (CGS units)
  ! =====================================================================
  REAL(dp), PARAMETER :: hp   = 6.626D-27          ! Planck constant [erg s]
  REAL(dp), PARAMETER :: kb   = 1.38D-16           ! Boltzmann constant [erg/K]
  REAL(dp), PARAMETER :: mH   = 1.67262192369D-24  ! Hydrogen atom mass [g]
  REAL(dp), PARAMETER :: muH  = 1.2D0              ! Mean atomic mass in units of mH
  REAL(dp), PARAMETER :: me   = 9.10938356D-28     ! Electron mass [g]
  REAL(dp), PARAMETER :: sigT = 6.6524587D-25      ! Thomson scattering cross-section [cm^2]
  REAL(dp), PARAMETER :: h2k  = hp / kb * 1.0D12   ! h/k scaled to THz [K/THz]
  REAL(dp), PARAMETER :: c    = 2.99792D10          ! Speed of light [cm/s]
  REAL(dp), PARAMETER :: pi   = 3.1415926D0

  ! =====================================================================
  ! Cosmological parameters (Planck 2018)
  ! =====================================================================
  REAL(dp), PARAMETER :: H0  = 6.736D6 / 3.086D24  ! Hubble constant [1/s]
  REAL(dp), PARAMETER :: hc  = 0.6736D0             ! Dimensionless Hubble parameter h
  REAL(dp), PARAMETER :: Omm = 0.3153D0             ! Total matter density parameter
  REAL(dp), PARAMETER :: OmL = 0.6847D0             ! Dark energy (Lambda) density parameter
  REAL(dp), PARAMETER :: Omb = 0.0493D0             ! Baryon density parameter
  REAL(dp), PARAMETER :: Yp  = 0.2446D0             ! Primordial helium mass fraction

  ! =====================================================================
  ! HeI reionization model: xHeI(z) = xHI(z) * [eta + (1-eta)*W(z)]
  !   where W(z) = 0.5 * [1 + tanh((z - zt) / dzt)]
  ! =====================================================================
  REAL(dp), PARAMETER :: zt  = 5.5D0   ! Transition redshift for HeI
  REAL(dp), PARAMETER :: dzt = 0.5D0   ! Transition width
  REAL(dp), PARAMETER :: eta = 0.1D0   ! Minimum HeI/HI ratio

  ! =====================================================================
  ! HeIII reionization model: xHeIII(z) = 0.5*[1 + tanh((zHeIII - z)/dzHeIII)]
  ! =====================================================================
  REAL(dp), PARAMETER :: zHeIII  = 3.1D0  ! HeII reionization redshift
  REAL(dp), PARAMETER :: dzHeIII = 0.4D0  ! HeII reionization width

  ! =====================================================================
  ! Observation / output redshift
  ! =====================================================================
  REAL(dp), PARAMETER :: zo = 0.0D0   ! Observer redshift (z = 0)
  REAL(dp), PARAMETER :: Dv = 1.78D02 ! Halo density contrast at virialization

  ! =====================================================================
  ! Present-day number densities [cm^-3]
  ! =====================================================================
  REAL(dp), PARAMETER :: nH0  = 0.10693D31 * H0**2 * Omb * (1.0D0 - Yp)  ! Hydrogen
  REAL(dp), PARAMETER :: nHe0 = 0.26732D30 * H0**2 * Omb * Yp            ! Helium

  ! =====================================================================
  ! HI photoionization cross-section: Verner et al. fitting parameters
  !   sigma(E) = sigma0 * (xpi - 1)^2 * xpi^(P/2 - 5.5) / (1 + sqrt(xpi/ya))^P
  !   where xpi = E / E0H (photon energy in units of ionization threshold)
  ! =====================================================================
  REAL(dp), PARAMETER :: sigma0 = 0.11083D-13  ! Reference cross-section [cm^2]
  REAL(dp), PARAMETER :: ya     = 23.424D0     ! Fitting parameter
  REAL(dp), PARAMETER :: P      = 2.3745D0     ! Fitting parameter
  REAL(dp), PARAMETER :: E0H    = 0.1637D-11   ! HI ionization energy [erg] (13.6 eV)

  ! =====================================================================
  ! Ionization thresholds [THz]
  ! =====================================================================
  REAL(dp), PARAMETER :: xfpi      = 0.3284D16    ! HI  threshold frequency [Hz] (= E0H / hp)
  REAL(dp), PARAMETER :: xfpiHeII  = 4.0D0 * 0.3284D16  ! HeII threshold (= 4 * xfpi for Z=2)

  ! =====================================================================
  ! HeI photoionization cross-section fitting parameters
  !   sigma_HeI(E) = 7.4e-18 * (a1/x^sp + (1-a1)/x^(1+sp))
  !                + 7.33e-22 * (E/keV)^{-3.5} * (1 + a2/sqrt(x) * exp(a3/sqrt(x)))
  !   where x = E / EpiHeI
  ! =====================================================================
  REAL(dp), PARAMETER :: a1      = 7.3861D00   ! Fitting parameter
  REAL(dp), PARAMETER :: a2      = -3.2491D00  ! Fitting parameter
  REAL(dp), PARAMETER :: a3      = 1.1783D00   ! Fitting parameter
  REAL(dp), PARAMETER :: sp      = 3.9119D00   ! Power-law index
  REAL(dp), PARAMETER :: EpiHeI  = 24.58D00 * 1.602D-12  ! HeI ionization energy [erg] (24.58 eV)
  REAL(dp), PARAMETER :: xfpiHeI = EpiHeI / hp            ! HeI threshold frequency [Hz]

  ! =====================================================================
  ! Work arrays
  ! =====================================================================
  REAL(dp), DIMENSION(n_xf) :: xf_vals       ! Frequency grid [THz]
  REAL(dp), DIMENSION(n_xf) :: F_vals        ! Transfer function values

  ! Auxiliary arrays (available for diagnostics if uncommented in the output)
  REAL(dp), DIMENSION(n_xf) :: y_vals, sigma_vals
  REAL(dp), DIMENSION(n_xf) :: taucHI_vals, tauLS_vals, tauLyHI_vals
  REAL(dp), DIMENSION(n_xf) :: taucHeI_vals, sigmaHeI_vals, sigma_Ly_a_vals

  REAL(dp), DIMENSION(n_z)  :: z_values, tauThom_vals, xHI_vals

  ! Lyman-series atomic data arrays
  REAL(dp), DIMENSION(n_ly)     :: lambda_lyi, f_lyi           ! HI line wavelengths [A] and oscillator strengths
  REAL(dp), DIMENSION(n_lyHeI)  :: lambdaHeI_lyi, fHeI_lyi     ! HeI line wavelengths [A] and oscillator strengths
  COMMON lambda_lyi, f_lyi, lambdaHeI_lyi, fHeI_lyi

  REAL(dp) :: xf_min, xf_max, xf_step, z, t_start, t_end, t_z_start
  INTEGER  :: i, j

  CALL CPU_TIME(t_start)
  WRITE(*, '(A)') '=== Spectral Transfer Function: H + He ==='
  WRITE(*, '(A, I5, A)') 'Frequency grid: ', n_xf, ' points (100 - 30000 THz)'
  WRITE(*, '(A, I2, A)') 'Redshift nodes: ', n_z, ' values'
  WRITE(*, *)

  ! -------------------------------------------------------------------
  ! Read atomic data: HI Lyman-series oscillator strengths
  ! -------------------------------------------------------------------
  WRITE(*, '(A)') 'Reading atomic data...'
  OPEN(UNIT=10, FILE='HOS.dat', STATUS='UNKNOWN', ACTION='READ')
  DO i = 1, n_ly
    READ(10, *) j, lambda_lyi(i), f_lyi(i)
  END DO
  CLOSE(10)

  ! -------------------------------------------------------------------
  ! Read atomic data: HeI line oscillator strengths
  ! -------------------------------------------------------------------
  OPEN(UNIT=11, FILE='HeIOS.dat', STATUS='UNKNOWN', ACTION='READ')
  DO i = 1, n_lyHeI
    READ(11, *) j, lambdaHeI_lyi(i), fHeI_lyi(i)
  END DO
  CLOSE(11)

  ! -------------------------------------------------------------------
  ! Build logarithmic frequency grid from xf_min to xf_max [THz]
  ! -------------------------------------------------------------------
  xf_min  = 100.0D0
  xf_max  = 30000.0D0
  xf_step = (DLOG10(xf_max) - DLOG10(xf_min)) / REAL(n_xf - 1, dp)
  DO i = 1, n_xf
    xf_vals(i) = 10**(DLOG10(xf_min) + REAL(i - 1, dp) * xf_step)
  END DO

  ! -------------------------------------------------------------------
  ! Define the source-redshift nodes at which the transfer function
  ! is evaluated (observer at zo = 0)
  ! -------------------------------------------------------------------
  z_values = (/ zo + 5.0D0,  zo + 5.5D0,  zo + 6.0D0,  zo + 6.5D0, &
                zo + 7.0D0,  zo + 8.0D0,  zo + 10.0D0, zo + 12.0D0, &
                zo + 15.0D0 /)

  ! -------------------------------------------------------------------
  ! Main computation loop: for each source redshift, compute F(nu, z)
  ! and write frequency–transmission table to file
  ! -------------------------------------------------------------------
  WRITE(*, *)
  DO i = 1, n_z
    z = z_values(i)
    CALL CPU_TIME(t_z_start)
    WRITE(*, '(A, I2, A, I2, A, F5.1, A)', ADVANCE='NO') &
      '[', i, '/', n_z, '] z = ', z, ' ... '

    OPEN(UNIT=12, &
         FILE='rtf_H-He_1_z_'//TRIM(ADJUSTL(WRITEF(z)))//'.dat', &
         STATUS='UNKNOWN', ACTION='WRITE')

    DO j = 1, n_xf
      F_vals(j) = F(xf_vals(j), z, zo)
      WRITE(12, '(E16.8, 1X, E16.8)') xf_vals(j), F_vals(j)
    END DO

    CLOSE(12)
    CALL CPU_TIME(t_end)
    WRITE(*, '(A, F7.1, A)') 'done (', t_end - t_z_start, 's)'
  END DO

  WRITE(*, *)
  WRITE(*, '(A, F8.1, A)') 'Total time: ', t_end - t_start, 's'
  WRITE(*, '(A)') 'All output files written.'

CONTAINS

  ! ===================================================================
  ! WRITEF — convert a real number to a trimmed string for file names
  ! ===================================================================
  FUNCTION WRITEF(x) RESULT(str)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: x
    CHARACTER(LEN=20)    :: str
    WRITE(str, '(F0.1)') x
  END FUNCTION WRITEF

  ! ===================================================================
  ! xHI — neutral hydrogen fraction as a function of redshift
  !
  ! Uses piecewise analytic fits:
  !   z <= 5.15 : Bosman, Becker & Inoue (2022) rational-exponential fit
  !   z >  5.15 : logistic model with coefficients from either
  !               Kageura & Bosman (2025) [klxHI=1] or
  !               Glazer & Bosman (2018)  [klxHI=2]
  !
  ! The pivot point is xHI(5.15) = 2.77125e-5 (Bosman average).
  ! ===================================================================
  FUNCTION xHI(z) RESULT(xHI_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: z
    REAL(dp) :: xHI_val, ax1, ax2, ax3, ax4, b1, b2, b3, tz, s1, s2
    REAL(dp) :: Az1, Bz1, Cz1, zC1, Az2, Bz2, zC2, ch1, zn1, Lz, xHI5

    ! BBI (Bosman-Becker-Inoue) fitting coefficients for z <= 5.15
    ax1 = 159.1068D0
    ax2 = 36.4547D0
    ax3 = -3.29829D0
    ax4 = -1.07885D0
    b1  = -100.000D0
    b2  = -42.4907D0
    b3  = -4.86781D0

    ! KB (Kageura-Bosman) coefficients for z > 5.15 (klxHI = 1)
    Az1 = 11.7083D0
    Bz1 = -4.06938D0
    Cz1 = 0.523253D0
    zC1 = 0.343817D0

    ! GB (Glazer-Bosman) coefficients for z > 5.15 (klxHI = 2)
    Az2 = 10.7053D0
    Bz2 = -1.52525D0
    zC2 = 1.14447D0

    ! Pivot value: Bosman averaged xHI at z = 5.15
    xHI5 = 2.77125D-05
    tz   = z - 5.15D0

    IF (z .LE. 5.15D0) THEN
      ! Low-redshift regime: rational-exponential fit (BBI)
      ch1  = tz * (ax1 + ax2*tz + ax3*tz**2 + ax4*tz**3)
      zn1  = 1.0 + b1*tz + b2*tz**2 + b3*tz**3
      xHI_val = xHI5 * DEXP(ch1 / zn1)
    ELSE
      ! High-redshift regime: logistic model
      IF (klxHI == 1) THEN
        ! KBBI approximation (Kageura & Bosman 2025)
        s1 = DLOG(1.0 + tz / zC1)
        Lz = DLOG(xHI5 / (1.0 - xHI5)) + Az1*s1 + Bz1*s1**2 + Cz1*s1**3
      ELSE
        ! GBI approximation (Glazer & Bosman 2018)
        s2 = DLOG(1.0 + tz / zC2)
        Lz = DLOG(xHI5 / (1.0 - xHI5)) + Az2*s2 + Bz2*s2**2
      END IF
      xHI_val = 1.0 / (1.0 + EXP(-1.0 * Lz))
    END IF
  END FUNCTION xHI

  ! ===================================================================
  ! xpi — dimensionless photon energy in units of HI ionization threshold
  !
  ! xpi = h * nu_obs / E0H, where nu_obs = xf * (1+z) / (1+zo) * 1e12
  ! ===================================================================
  FUNCTION xpi(xf, z, zo) RESULT(xpi_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: xpi_val
    xpi_val = 0.10D13 * hp * xf * (z + 1.0D0) / ((1.0D0 + zo) * E0H)
  END FUNCTION xpi

  ! ===================================================================
  ! sigma — HI photoionization cross-section [cm^2]
  !
  ! Verner et al. analytic fit. Returns zero below the HI threshold.
  ! ===================================================================
  FUNCTION sigma(xf, z, zo) RESULT(sigma_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: sigma_val, xpi_val

    xpi_val = xpi(xf, z, zo)

    IF (xf * (1.0D0 + z) / (1.0D0 + zo) < xfpi / 1.0D12) THEN
      sigma_val = 0.0D0
    ELSE
      sigma_val = sigma0 * (xpi_val - 1.0D0)**2 * xpi_val**(P/2.0D0 - 5.5D0) / &
                  (1.0D0 + DSQRT(xpi_val / ya))**P
    END IF
  END FUNCTION sigma

  ! ===================================================================
  ! sigmaHeII — HeII photoionization cross-section [cm^2]
  !
  ! Hydrogen-like with Z=2: threshold is 4 x HI, cross-section is 1/4.
  ! ===================================================================
  FUNCTION sigmaHeII(xf, z, zo) RESULT(sigmaHeII_v)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: sigmaHeII_v, xpi_val

    xpi_val = xpi(xf, z, zo) / 4.0D0

    IF (xf * (1.0D0 + z) / (1.0D0 + zo) < xfpiHeII / 1.0D12) THEN
      sigmaHeII_v = 0.0D0
    ELSE
      sigmaHeII_v = 0.25D0 * sigma0 * (xpi_val - 1.0D0)**2 * &
                    xpi_val**(P/2.0D0 - 5.5D0) / &
                    (1.0D0 + DSQRT(xpi_val / ya))**P
    END IF
  END FUNCTION sigmaHeII

  ! ===================================================================
  ! sigmaHeI — HeI photoionization cross-section [cm^2]
  !
  ! Two-component fit: power-law near threshold + high-energy tail.
  ! ===================================================================
  FUNCTION sigmaHeI(xf, z, zo) RESULT(sigmaHeI_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: sigmaHeI_val, Eval, x

    ! Photon energy at the source redshift [erg]
    Eval = hp * xf * 1.0D12 * (1.0D0 + z) / (1.0D0 + zo)
    ! Energy in units of HeI ionization threshold
    x = Eval / EpiHeI

    IF (xf * (1.0D0 + z) / (1.0D0 + zo) < xfpiHeI / 1.0D12) THEN
      sigmaHeI_val = 0.0D0
    ELSE
      sigmaHeI_val = 7.4D-18 * (a1 / x**sp + (1.0D0 - a1) / x**(1.0D0 + sp)) + &
                     7.33D-22 / (Eval / 1.602D-09)**3.5 * &
                     (1.0D0 + a2 / DSQRT(x) / DEXP(a3 / DSQRT(x)))
    END IF
  END FUNCTION sigmaHeI

  ! ===================================================================
  ! tauThom — Thomson scattering optical depth from zo to z
  !
  ! tau_T = (c * nH0 / H0) * integral_{zo}^{z} x_e(z') * (1+z')^2
  !         / sqrt(Omm*(1+z')^3 + OmL) * sigT dz'
  !
  ! The free-electron fraction x_e accounts for ionized H, HeII, HeIII.
  ! Integration uses Simpson's rule with n_steps = 8000.
  ! ===================================================================
  FUNCTION tauThom(z, zo) RESULT(tauThom_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: z, zo
    REAL(dp) :: tauThom_val, integrand, h, x
    REAL(dp) :: xHeI, xHeII, xHeIII, xHIv, xev
    REAL(dp) :: xzHeIII, arHeIII, dyzHeIII, W
    INTEGER  :: n_steps, i

    ! Precompute HeIII reionization tanh parameters
    xzHeIII  = (1.0D0 + zHeIII)**1.5
    dyzHeIII = 1.5D0 * DSQRT(1.0D0 + zHeIII) * dzHeIII

    n_steps     = 8000
    h           = (z - zo) / REAL(n_steps, dp)
    tauThom_val = 0.0D0

    DO i = 0, n_steps - 1
      x = zo + REAL(i, dp) * h

      ! HeI neutral fraction via tanh model linked to xHI
      W    = 0.5D0 * (1.0D0 + DTANH((x - zt) / dzt))
      xHeI = xHI(x) * (eta + (1.0D0 - eta) * W)

      ! HeIII ionized fraction via tanh reionization model
      arHeIII = (xzHeIII - (1.0D0 + x)**1.5) / dyzHeIII
      xHeIII  = 0.5D0 * (1.0D0 + DTANH(arHeIII))

      ! HeII = complement of HeI + HeIII
      xHeII = 1.0D0 - xHeI - xHeIII

      ! Free-electron fraction: contributions from ionized H and He species
      xHIv = xHI(x)
      xev  = 1.0D0 - xHIv + (xHeII + 2.0D0 * xHeIII) * nHe0 / nH0

      integrand = xev * ((1.0D0 + x)**2) * &
                  sigT / DSQRT(Omm * ((1.0D0 + x)**3) + OmL)

      ! Simpson's rule weights
      IF (i == 0 .OR. i == n_steps - 1) THEN
        tauThom_val = tauThom_val + integrand
      ELSE IF (MOD(i, 2) == 0) THEN
        tauThom_val = tauThom_val + 2.0D0 * integrand
      ELSE
        tauThom_val = tauThom_val + 4.0D0 * integrand
      END IF
    END DO

    tauThom_val = c * nH0 * h / 3.0D0 * tauThom_val / H0
  END FUNCTION tauThom

  ! ===================================================================
  ! taucHI — HI continuum (photoionization) optical depth from zo to z
  !
  ! tau_c^HI = (c * nH0 / H0) * integral xHI(z') * sigma_HI(nu, z')
  !            * (1+z')^2 / sqrt(...) dz'
  ! ===================================================================
  FUNCTION taucHI(xf, z, zo) RESULT(taucHI_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: taucHI_val, integrand, h, x, xHI_val
    INTEGER  :: n_steps, i

    n_steps    = 8000
    h          = (z - zo) / REAL(n_steps, dp)
    taucHI_val = 0.0D0

    DO i = 0, n_steps - 1
      x       = zo + REAL(i, dp) * h
      xHI_val = xHI(x)

      integrand = xHI_val * ((1.0D0 + x)**2) * &
                  sigma(xf, x, zo) / DSQRT(Omm * ((1.0D0 + x)**3) + OmL)

      ! Simpson's rule weights
      IF (i == 0 .OR. i == n_steps - 1) THEN
        taucHI_val = taucHI_val + integrand
      ELSE IF (MOD(i, 2) == 0) THEN
        taucHI_val = taucHI_val + 2.0D0 * integrand
      ELSE
        taucHI_val = taucHI_val + 4.0D0 * integrand
      END IF
    END DO

    taucHI_val = c * nH0 * h / 3.0D0 * taucHI_val / H0
  END FUNCTION taucHI

  ! ===================================================================
  ! taucHeII — HeII continuum (photoionization) optical depth from zo to z
  ! ===================================================================
  FUNCTION taucHeII(xf, z, zo) RESULT(taucHeII_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: taucHeII_val, integrand, h, x
    REAL(dp) :: xHeI, xHeII, xHeIII
    REAL(dp) :: xzHeIII, arHeIII, dyzHeIII, W
    INTEGER  :: n_steps, i

    xzHeIII  = (1.0D0 + zHeIII)**1.5
    dyzHeIII = 1.5D0 * DSQRT(1.0D0 + zHeIII) * dzHeIII

    n_steps       = 8000
    h             = (z - zo) / REAL(n_steps, dp)
    taucHeII_val  = 0.0D0

    DO i = 0, n_steps - 1
      x = zo + REAL(i, dp) * h

      ! Compute HeII fraction at this redshift
      W       = 0.5D0 * (1.0D0 + DTANH((x - zt) / dzt))
      xHeI    = xHI(x) * (eta + (1.0D0 - eta) * W)
      arHeIII = (xzHeIII - (1.0D0 + x)**1.5) / dyzHeIII
      xHeIII  = 0.5D0 * (1.0D0 + DTANH(arHeIII))
      xHeII   = 1.0D0 - xHeI - xHeIII

      integrand = xHeII * ((1.0D0 + x)**2) * &
                  sigmaHeII(xf, x, zo) / DSQRT(Omm * ((1.0D0 + x)**3) + OmL)

      ! Simpson's rule weights
      IF (i == 0 .OR. i == n_steps - 1) THEN
        taucHeII_val = taucHeII_val + integrand
      ELSE IF (MOD(i, 2) == 0) THEN
        taucHeII_val = taucHeII_val + 2.0D0 * integrand
      ELSE
        taucHeII_val = taucHeII_val + 4.0D0 * integrand
      END IF
    END DO

    taucHeII_val = c * nHe0 * h / 3.0D0 * taucHeII_val / H0
  END FUNCTION taucHeII

  ! ===================================================================
  ! taucHeI — HeI continuum (photoionization) optical depth from zo to z
  ! ===================================================================
  FUNCTION taucHeI(xf, z, zo) RESULT(taucHeI_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: taucHeI_val, integrand, h, x, xHeI, W
    INTEGER  :: n_steps, i

    n_steps      = 8000
    h            = (z - zo) / REAL(n_steps, dp)
    taucHeI_val  = 0.0D0

    DO i = 0, n_steps - 1
      x = zo + REAL(i, dp) * h

      ! HeI neutral fraction
      W    = 0.5D0 * (1.0D0 + DTANH((x - zt) / dzt))
      xHeI = xHI(x) * (eta + (1.0D0 - eta) * W)

      integrand = xHeI * ((1.0D0 + x)**2) * &
                  sigmaHeI(xf, x, zo) / DSQRT(Omm * ((1.0D0 + x)**3) + OmL)

      ! Simpson's rule weights
      IF (i == 0 .OR. i == n_steps - 1) THEN
        taucHeI_val = taucHeI_val + integrand
      ELSE IF (MOD(i, 2) == 0) THEN
        taucHeI_val = taucHeI_val + 2.0D0 * integrand
      ELSE
        taucHeI_val = taucHeI_val + 4.0D0 * integrand
      END IF
    END DO

    taucHeI_val = c * nHe0 * h / 3.0D0 * taucHeI_val / H0
  END FUNCTION taucHeI

  ! ===================================================================
  ! tauLyHI — HI Lyman-series line optical depth
  !
  ! Sums Gunn-Peterson optical depth for each Lyman line (Ly-alpha
  ! through Ly-39). Each line contributes only when the observed
  ! frequency falls within the resonance redshift window.
  ! ===================================================================
  FUNCTION tauLyHI(xf, z, zo) RESULT(tauLyHI_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: tauLyHI_val, xHI_val, nu_i, la_i, f_i, xz, tau_i
    INTEGER  :: i

    tauLyHI_val = 0.0D0

    DO i = 1, n_ly
      f_i  = f_lyi(i)
      la_i = lambda_lyi(i) * 1.0D-08     ! Convert Angstroms to cm
      nu_i = c / la_i / 1.0D12           ! Line centre frequency [THz]
      tau_i = 0.0D0

      ! Check if observed frequency lies within this line's redshift window
      IF (xf >= nu_i * (1.0D0 + zo) / (1.0D0 + z) .AND. xf < nu_i) THEN
        xz      = nu_i * (1.0D0 + zo) / xf - 1.0D0   ! Absorption redshift
        xHI_val = xHI(xz)
        ! Gunn-Peterson optical depth for this line
        tau_i = 2.8378D28 * la_i * f_i * H0 * (1.0D0 - Yp) * Omb * &
                (1.0D0 + xz)**3 * xHI_val / DSQRT(Omm * (1.0D0 + xz)**3 + OmL)
      END IF

      tauLyHI_val = tauLyHI_val + tau_i
    END DO
  END FUNCTION tauLyHI

  ! ===================================================================
  ! tauLyHeII — HeII Lyman-series line optical depth
  !
  ! HeII is hydrogen-like with Z=2, so its Lyman wavelengths are 1/4
  ! of the HI values. The HeII fraction is computed from the helium
  ! ionization model.
  ! ===================================================================
  FUNCTION tauLyHeII(xf, z, zo) RESULT(tauLyHeII_v)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: tauLyHeII_v, nu_i, la_i, f_i, xz, tau_i
    REAL(dp) :: xHeI, xHeII, xHeIII
    REAL(dp) :: xzHeIII, arHeIII, dyzHeIII, W
    INTEGER  :: i

    xzHeIII  = (1.0D0 + zHeIII)**1.5
    dyzHeIII = 1.5D0 * DSQRT(1.0D0 + zHeIII) * dzHeIII

    tauLyHeII_v = 0.0D0

    DO i = 1, n_ly
      f_i  = f_lyi(i)
      la_i = lambda_lyi(i) * 0.25D-08    ! HeII wavelengths = HI/4 [cm]
      nu_i = c / la_i / 1.0D12           ! Line centre frequency [THz]
      tau_i = 0.0D0

      IF (xf >= nu_i * (1.0D0 + zo) / (1.0D0 + z) .AND. xf < nu_i) THEN
        xz = nu_i * (1.0D0 + zo) / xf - 1.0D0

        ! Compute HeII fraction at absorption redshift
        W       = 0.5D0 * (1.0D0 + DTANH((xz - zt) / dzt))
        xHeI    = xHI(xz) * (eta + (1.0D0 - eta) * W)
        arHeIII = (xzHeIII - (1.0D0 + xz)**1.5) / dyzHeIII
        xHeIII  = 0.5D0 * (1.0D0 + DTANH(arHeIII))
        xHeII   = 1.0D0 - xHeI - xHeIII

        tau_i = 7.0945D27 * la_i * f_i * H0 * Yp * Omb * &
                (1.0D0 + xz)**3 * xHeII / DSQRT(Omm * (1.0D0 + xz)**3 + OmL)
      END IF

      tauLyHeII_v = tauLyHeII_v + tau_i
    END DO
  END FUNCTION tauLyHeII

  ! ===================================================================
  ! tauLyHeI — HeI line optical depth (using HeI-specific line data)
  ! ===================================================================
  FUNCTION tauLyHeI(xf, z, zo) RESULT(tauLyHeI_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: tauLyHeI_val, nu_i, la_i, f_i, xz, tau_i, xHeI, W
    INTEGER  :: i

    tauLyHeI_val = 0.0D0

    DO i = 1, n_lyHeI
      f_i  = fHeI_lyi(i)
      la_i = lambdaHeI_lyi(i) * 1.0D-08  ! Convert Angstroms to cm
      nu_i = c / la_i / 1.0D12           ! Line centre frequency [THz]
      tau_i = 0.0D0

      IF (xf >= nu_i * (1.0D0 + zo) / (1.0D0 + z) .AND. xf < nu_i) THEN
        xz = nu_i * (1.0D0 + zo) / xf - 1.0D0

        ! HeI neutral fraction at absorption redshift
        W    = 0.5D0 * (1.0D0 + DTANH((xz - zt) / dzt))
        xHeI = xHI(xz) * (eta + (1.0D0 - eta) * W)

        tau_i = 7.0945D27 * la_i * f_i * H0 * Yp * Omb * &
                (1.0D0 + xz)**3 * xHeI / DSQRT(Omm * (1.0D0 + xz)**3 + OmL)
      END IF

      tauLyHeI_val = tauLyHeI_val + tau_i
    END DO
  END FUNCTION tauLyHeI

  ! ===================================================================
  ! F — total spectral transfer function
  !
  ! Combines all opacity sources:
  !   tau_total = tau_c^HI + tau_Ly^HI
  !            + tau_c^HeI + tau_Ly^HeI
  !            + tau_c^HeII + tau_Ly^HeII
  !
  ! Returns F = exp(-tau_total), clamped to 0 for tau > 70 to avoid
  ! floating-point underflow.
  !
  ! Note: Thomson scattering (tauThom) is computed but currently
  ! excluded from the total — uncomment to include.
  ! ===================================================================
  FUNCTION F(xf, zp, zo) RESULT(F_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, zp, zo
    REAL(dp) :: F_val, tauz, tauzHI, tauzHeI, tauzHeII, tauTh

    ! HI contributions: Lyman-series lines + continuum
    tauzHI   = tauLyHI(xf, zp, zo)   + taucHI(xf, zp, zo)
    ! HeI contributions: lines + continuum
    tauzHeI  = tauLyHeI(xf, zp, zo)  + taucHeI(xf, zp, zo)
    ! HeII contributions: Lyman-series lines + continuum
    tauzHeII = taucHeII(xf, zp, zo)  + tauLyHeII(xf, zp, zo)
    ! Thomson scattering (computed but not added to total)
    tauTh    = tauThom(zp, zo)

    tauz = tauzHI + tauzHeI + tauzHeII  ! + tauTh

    F_val = 0.0D0
    IF (tauz <= 70.D0) THEN
      F_val = DEXP(-tauz)
    END IF
  END FUNCTION F

END PROGRAM main
