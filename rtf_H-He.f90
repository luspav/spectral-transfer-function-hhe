!=======================================================================
! Program: Spectral Transfer Function Calculator (rtf_H-He.f90)
! Description: Computes the hydrogen/helium spectral transfer function,
!   the primary data product discussed in the accompanying research article.
! Authors:
!   - Bohdan Novosyadlyj <bnovos@gmail.com>
!   - Pavlo Kopach <luspav@gmail.com>
! License: MIT License (see README.md for full text).
!=======================================================================

!***********************************************************************
! Fortran90 Code
!***********************************************************************

! # How to use:
! # MacOS
! gfortran -isysroot $(xcrun --sdk macosx --show-sdk-path) -o ss_sed ./ss_sed.f90
! # Linux
! gfortran -o ss_sed ss_sed.f90 or ifort -o ss_sed ss_sed.f90
! ./ss_sed or time ./ss_sed
! # Windows
! gfortran -o ss_sed ss_sed.f90 or ifort -o ss_sed ss_sed.f90
! activate ss_sed.exe by Enter


PROGRAM main
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: n_xf = 32000
  INTEGER, PARAMETER :: n_z = 9
  INTEGER, PARAMETER :: n_ly = 39
  INTEGER, PARAMETER :: n_lyHeI = 10 
  INTEGER, PARAMETER :: klxHI = 2 !1- approximation xHI for Bosman et al. (2022) + Kageura et al. (2025)
                                   !2- approximation xHI for Bosman et al. (2022) + Glazer et al. (2018)
  ! Constants
!  REAL(dp), PARAMETER :: Rsun = 6.96342D10        ! Star radius in cm
  REAL(dp), PARAMETER :: hp = 6.626D-27           ! Planck constant in erg*s
  REAL(dp), PARAMETER :: kb = 1.38D-16            ! Boltzmann constant in erg/K
  REAL(dp), PARAMETER :: mH = 1.67262192369D-24   ! Mass of hydrogen atom in g 
  REAL(dp), PARAMETER :: muH = 1.2D0              ! середня атомна одиниця маси?
  REAL(dp), PARAMETER :: me = 9.10938356D-28      ! Mass of electron in g
  REAL(dp), PARAMETER :: sigT = 6.6524587D-25     ! effective crosssection of Thompson scattering (cm^2)
  REAL(dp), PARAMETER :: h2k = hp / kb * 1.0D12   ! Planck constant divided by Boltzmann constant times 1e12
  REAL(dp), PARAMETER :: c = 2.99792D10           ! Speed of light in cm/s
  REAL(dp), PARAMETER :: H0 = 6.736D6 / 3.086D24  ! Hubble constant in 1/s
  REAL(dp), PARAMETER :: hc = 0.6736D0            ! dimensionless Hubble constant
  REAL(dp), PARAMETER :: Omm = 0.3153D0           ! Matter density parameter
  REAL(dp), PARAMETER :: OmL = 0.6847D0           ! Dark energy density parameter
  REAL(dp), PARAMETER :: Omb = 0.0493D0           ! Baryonic matter density parameter
  REAL(dp), PARAMETER :: Yp = 0.2446D0            ! Primordial mass fraction of helium
  REAL(dp), PARAMETER :: zt = 5.5D0               ! xHeI(z) = xHI(z)(η + (1 − η) W(z),  W(z) =  0.5[1 + tanh(z-zt)/ ∆z]
  REAL(dp), PARAMETER :: dzt = 0.5D0
  REAL(dp), PARAMETER :: eta = 1.0D0
  REAL(dp), PARAMETER :: zHeIII=3.1D0             ! xHeIII(z) =  0.5[1 + tanh(zHeII-z)/dzHeIII]
  REAL(dp), PARAMETER :: dzHeIII=0.4D0 
!  REAL(dp), PARAMETER :: xHII0 = 0.999999D0          ! Fraction of ionized hydrogen at z = 6
!  REAL(dp), PARAMETER :: zre = 8.0D0              ! Redshift of reionization
!  REAL(dp), PARAMETER :: Dy = 2.0D0               ! Transition width in redshift space
  REAL(dp), PARAMETER :: zo = 0.0D0               ! target redshift
  REAL(dp), PARAMETER :: Dv = 1.78D02             ! Density contrast of halos at the moment of virializatio

  !REAL(dp), PARAMETER :: pi = 3.1415926535897932384626433832795D0
  REAL(dp), PARAMETER :: pi = 3.1415926D0
  ! Total hydrogen number density at z = 0
  REAL(dp), PARAMETER :: nH0 = 0.10693D31 * H0**2 * Omb * (1.0D0 - Yp)
 ! Total hydrogen number density at z = 0
  REAL(dp), PARAMETER :: nHe0 = 0.26732D30 * H0**2 * Omb * Yp 
  ! Pre-computed value for the reionization redshift
!  REAL(dp), PARAMETER :: yzre = (1.0D0 + zre)**1.5D0
  ! Reference cross-section for photoionization in cm^2
  REAL(dp), PARAMETER :: sigma0 = 0.11083D-13
  ! Fitting parameters
  REAL(dp), PARAMETER :: ya = 23.424D0
  REAL(dp), PARAMETER :: P = 2.3745D0
  ! Hydrogen ionization energy in erg
  REAL(dp), PARAMETER :: E0H = 0.1637D-11
  !Ionization potential of HI in Hz
  REAL(dp), PARAMETER :: xfpi = 0.3284D16
  !Ionization potential of HeII in Hz
  REAL(dp), PARAMETER :: xfpiHeII = 4.0D0*0.3284D16  
  !Fitting parameter in HeI photoinization cross-section
  REAL(dp), PARAMETER :: a1=7.3861D00
  !Fitting parameter in HeI photoinization cross-section
  REAL(dp), PARAMETER :: a2=-3.2491D00
  !Fitting parameter in HeI photoinization cross-section
  REAL(dp), PARAMETER :: a3=1.1783D00
  !Fitting parameter in HeI photoinization cross-section
  REAL(dp), PARAMETER :: sp=3.9119D00
  !Potential of HeI photoinization
  REAL(dp), PARAMETER ::    EpiHeI = 24.58D00*1.602D-12
  !Threshold frequency of HeI photoionization
  REAL(dp), PARAMETER :: xfpiHeI=EpiHeI/hP

  ! Variables
  REAL(dp), DIMENSION(n_xf) :: xf_vals, sigma_Ly_a_vals
  REAL(dp), DIMENSION(n_xf) :: y_vals, sigma_vals, tauc_vals,  tauLS_vals,  tauGP_vals, taucHeI_vals
  REAL(dp), DIMENSION(n_xf) :: F_vals, sigmaHeI_vals 
 
  REAL(dp), DIMENSION(n_z) :: z_values, tauThom_vals, xHI_vals 
  REAL(dp), DIMENSION(n_ly) :: lambda_lyi,f_lyi
  REAL(dp), DIMENSION(n_lyHeI) :: lambdaHeI_lyi,fHeI_lyi
  COMMON lambda_lyi,f_lyi,lambdaHeI_lyi,fHeI_lyi

  REAL(dp) :: xf_min, xf_max, xf_step, z 
  INTEGER :: i, j
 
  
  
  OPEN(UNIT=10, FILE='../AtomicData/HOS.dat', STATUS='UNKNOWN', ACTION='READ')
  DO i=1,n_ly
  READ(10,*) j,lambda_lyi(i),f_lyi(i)
  END DO
  CLOSE(10)
  
  OPEN(UNIT=11, FILE='../AtomicData/HeIOS.dat', STATUS='UNKNOWN', ACTION='READ')
  DO i=1,n_lyHeI
  READ(11,*) j,lambdaHeI_lyi(i),fHeI_lyi(i)
  END DO 
  CLOSE(11)
    OPEN(UNIT=13, FILE='gff.dat', STATUS='UNKNOWN', ACTION='WRITE')
  ! Initialize xf_valssigma_Ly_a_vals(j)
  xf_min = 100.0D0
  xf_max = 30000.0D0
  xf_step = (dlog10(xf_max) - dlog10(xf_min)) / REAL(n_xf - 1, dp)
  DO i = 1, n_xf
    xf_vals(i) = 10**(dlog10(xf_min) + REAL(i - 1, dp) * xf_step)
  END DO
 
  z_values = (/zo+5.0D0, zo+5.5D0, zo+6.0D0, zo+6.5D0, zo+7.0D0, zo+8.0D0, zo+10.0D0, zo+12.0D0, zo+15.0D0/)
!    z_values = (/zo+15.0D0/)
  ! Loop over z_values to compute other functions
  DO i = 1,n_z
    z = z_values(i)
!    tauThom_vals=tauThom(z, zo)
    ! Open files for each z
    OPEN(UNIT=12, FILE='rtf_H-He_2_z_'//TRIM(ADJUSTL(WRITEF(z)))//'.dat', STATUS='UNKNOWN', ACTION='WRITE')
 
    DO j = 1, n_xf
!    tauLS_vals(j)=tauLyHeII(xf_vals(j), z, zo)
!    tauc_vals(j)=taucHeII(xf_vals(j), z, zo)
    F_vals(j) = F(xf_vals(j), z, zo)
    WRITE(12, '(E16.8, 1X, 9E16.8)') xf_vals(j), F_vals(j)!, tauLS_vals(j), tauc_vals(j) 
    END DO

    ! Close files for this z
    CLOSE(12)
 
  END DO

CONTAINS

  ! Function to convert real to string
  FUNCTION WRITEF(x) RESULT(str)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: x
    CHARACTER(LEN=20) :: str
    WRITE(str, '(F0.1)') x
  END FUNCTION WRITEF
  
  ! Function xHI
  FUNCTION xHI(z) RESULT(xHI_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: z
    REAL(dp) :: xHI_val,AH,BH,CH,zre,dzre,k 
    REAL(dp) :: xHI_12,xHIk,zk,ay,by,xHI_15
    IF (klxHI == 1) then
    AH=0.25066D-04                     
    BH=0.90177D0-AH
    zre=0.60596D+01  
    dzre=0.33974D0  
    k=0.49822D+01
    ELSE
    AH=0.24115D-04     
    BH=0.99978D0-AH
    zre=0.67042D+01
    dzre=0.55666D+00
    k=0.42456D+01
    ENDIF
      CH = 1.0D0/DEXP((z-zre)/dzre)
      xHI_val = AH+BH/(1.0D0+CH)**k
    END FUNCTION xHI 

  ! Function xpi
  FUNCTION xpi(xf, z, zo) RESULT(xpi_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: xpi_val
    xpi_val = 0.10D13 * hp * xf * (z + 1.0D0) / ((1.0D0 + zo) * E0H)
  END FUNCTION xpi
 
! Function sigma
  FUNCTION sigma(xf, z, zo) RESULT(sigma_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: sigma_val, xpi_val
    xpi_val = xpi(xf, z, zo)
    IF (xf * (1.0D0 + z)/(1.0D0 + zo) < xfpi / 1.0D12) THEN
      sigma_val = 0.0D0
    ELSE
      sigma_val = sigma0*(xpi_val-1.0D0)**2*xpi_val**(P/2.0D0-5.5D0)/&
                  (1.0D0+DSQRT(xpi_val/ya))**P
    END IF
  END FUNCTION sigma
  
! Function sigmaHeII
  FUNCTION sigmaHeII(xf, z, zo) RESULT(sigmaHeII_v)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: sigmaHeII_v, xpi_val
    xpi_val = xpi(xf, z, zo)/4.0D0
    IF (xf * (1.0D0 + z)/(1.0D0 + zo) < xfpiHeII / 1.0D12) THEN
      sigmaHeII_v = 0.0D0
    ELSE
      sigmaHeII_v = 0.25D0*sigma0*(xpi_val-1.0D0)**2*xpi_val**(P/2.0D0-5.5D0)/&
                  (1.0D0+DSQRT(xpi_val/ya))**P
    END IF
  END FUNCTION sigmaHeII    
  
  ! Function sigmaHeI
  FUNCTION sigmaHeI(xf, z, zo) RESULT(sigmaHeI_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: sigmaHeI_val, Eval, x
     
    Eval= hP*xf*1.0D12*(1.0D0 + z)/(1.0D0 + zo)
    x=Eval/EpiHeI
    IF (xf * (1.0D0 + z)/(1.0D0 + zo) < xfpiHeI / 1.0D12) THEN
      sigmaHeI_val = 0.0D0
    ELSE
      sigmaHeI_val = 7.4D-18*(a1/x**sp+(1.0D0-a1)/x**(1.0D0+sp))+&
    7.33D-22/(Eval/1.602D-09)**3.5*(1.0D0+a2/DSQRT(x)/DEXP(a3/DSQRT(x)))
    END IF
  END FUNCTION sigmaHeI 
  
! Function tauThompson
  FUNCTION tauThom(z, zo) RESULT(tauThom_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: z, zo
    REAL(dp) :: tauThom_val, integrand, h, x, ye, xHeI, xHeII, xHeIII, xHIv, xev
    REAL(dp) :: xzHeIII, arHeIII, thHeIII, dyzHeIII, W     
    INTEGER :: n_steps, i
    
    xzHeIII=(1.0D0+zHeIII)**1.5
    dyzHeIII=1.5D0*DSQRT(1.0D0+zHeIII)*dzHeIII 

    n_steps = 8000
    h = (z - zo) / REAL(n_steps, dp)
    tauThom_val = 0.0D0

    DO i = 0, n_steps - 1
      x = zo + REAL(i, dp) * h
      W=0.5D0*(1.0D0+DTANH((x-zt)/dzt))    
      xHeI=xHI(x)*(eta+(1.0D0-eta)*W)
     arHeIII=(xzHeIII-(1.0D0+x)**1.5)/dyzHeIII
     xHeIII=0.5D0*(1.0D0+DTANH(arHeIII))
     xHeII=1.0D0-xHeI-xHeIII      
    xHIv=xHI(x)
    xev=1.0D0-xHIv+(xHeII+2.0D0*xHeII)*nHe0/nH0
  
      integrand = xev * ((1.0D0 + x)**2) * &
                  sigT / DSQRT(Omm * ((1.0D0 + x)**3) + OmL)
      IF (i == 0 .OR. i == n_steps - 1) THEN
        tauThom_val = tauThom_val + integrand
      ELSE IF (MOD(i, 2) == 0) THEN
        tauThom_val = tauThom_val + 2.0D0 * integrand
      ELSE
        tauThom_val = tauThom_val + 4.0D0 * integrand
      END IF
!      IF (x == 15.0D0) THEN
!      write (*,'(6E12.4)') x,xHeI,xHeII,xHeIII
!      END IF
    END DO

    tauThom_val = c * nH0 * h / 3.0D0 * tauThom_val / H0
    
!        write (*,'(6E12.4)') z,tauThom_val,xHIv,xHeI,xHeII,xHeIII

  END FUNCTION tauThom    
  
  ! Function tauc
  FUNCTION tauc(xf, z, zo) RESULT(tauc_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: tauc_val, integrand, h, x,xHI_val
    INTEGER :: n_steps, i

    n_steps = 8000
    h = (z - zo) / REAL(n_steps, dp)
    tauc_val = 0.0D0

    DO i = 0, n_steps - 1
      x = zo + REAL(i, dp) * h
      xHI_val=xHI(x)
      integrand = xHI_val * ((1.0D0 + x)**2) * &
                  sigma(xf, x, zo) / DSQRT(Omm * ((1.0D0 + x)**3) + OmL)
      IF (i == 0 .OR. i == n_steps - 1) THEN
        tauc_val = tauc_val + integrand
      ELSE IF (MOD(i, 2) == 0) THEN
        tauc_val = tauc_val + 2.0D0 * integrand
      ELSE
        tauc_val = tauc_val + 4.0D0 * integrand
      END IF
    END DO

    tauc_val = c * nH0 * h / 3.0D0 * tauc_val / H0
  END FUNCTION tauc
  
! Function taucHeII
  FUNCTION taucHeII(xf, z, zo) RESULT(taucHeII_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: taucHeII_val, integrand, h, x, ye, xHeI, xHeII, xHeIII 
    REAL(dp) ::  xzHeIII, arHeIII, thHeIII, dyzHeIII, W   
    INTEGER :: n_steps, i

    xzHeIII=(1.0D0+zHeIII)**1.5
    dyzHeIII=1.5D0*DSQRT(1.0D0+zHeIII)*dzHeIII 
    n_steps = 8000
    h = (z - zo) / REAL(n_steps, dp)
    taucHeII_val = 0.0D0

    DO i = 0, n_steps - 1
      x = zo + REAL(i, dp) * h
      W=0.5D0*(1.0D0+DTANH((x-zt)/dzt))    
      xHeI=xHI(x)*(eta+(1.0D0-eta)*W)
     arHeIII=(xzHeIII-(1.0D0+x)**1.5)/dyzHeIII
     xHeIII=0.5D0*(1.0D0+DTANH(arHeIII))
     xHeII=1.0D0-xHeI-xHeIII      
      
      integrand = xHeII * ((1.0D0 + x)**2) * &
                  sigmaHeII(xf, x, zo) / DSQRT(Omm * ((1.0D0 + x)**3) + OmL)
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

  ! Function taucHeI
  FUNCTION taucHeI(xf, z, zo) RESULT(taucHeI_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: taucHeI_val, integrand, h, x, ye, xHeI, W    
    INTEGER :: n_steps, i
! 
    n_steps = 8000
    h = (z - zo) / REAL(n_steps, dp)
    taucHeI_val = 0.0D0

    DO i = 0, n_steps - 1
      x = zo + REAL(i, dp) * h
      W=0.5D0*(1.0D0+DTANH((x-zt)/dzt))    
      xHeI=xHI(x)*(eta+(1.0D0-eta)*W)
      integrand = xHeI * ((1.0D0 + x)**2) * &
                  sigmaHeI(xf, x, zo) / DSQRT(Omm * ((1.0D0 + x)**3) + OmL)
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
  
    FUNCTION tauGP(xf, z, zo) RESULT(tauGP_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: tauGP_val, xHI_val, nu_i, la_i, f_i, xz, tau_i
    INTEGER :: i
    tauGP_val=0.0D0
    DO i = 1,n_ly
    f_i = f_lyi(i)
    la_i = lambda_lyi(i)*1.0D-08
    nu_i = c/la_i/1.0D12
    tau_i = 0.0D0
    IF (xf.ge.nu_i*(1.0D0+zo)/(1.0D0+z).AND.xf.lt.nu_i) THEN
    xz = nu_i*(1.0D0+zo)/xf-1.0D0
    xHI_val = xHI(xz)    
    tau_i =  2.8378D28*la_i*f_i*H0*(1.0D0-Yp)*Omb*(1.0D0+xz)**3*xHI_val/DSQRT(Omm*(1.0D0+xz)**3+OmL)   
    END IF
    tauGP_val = tauGP_val + tau_i
    END DO
    
  END FUNCTION tauGP
  
    FUNCTION tauLyHeII(xf, z, zo) RESULT(tauLyHeII_v)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: tauLyHeII_v, nu_i, la_i, f_i, xz, ye, tau_i, xHeI, xHeII, xHeIII
    REAL(dp) :: xzHeIII, arHeIII, thHeIII, dyzHeIII, W     
    INTEGER :: i

    xzHeIII=(1.0D0+zHeIII)**1.5
    dyzHeIII=1.5D0*DSQRT(1.0D0+zHeIII)*dzHeIII 
    tauLyHeII_v=0.0D0
    DO i = 1,n_ly
    f_i = f_lyi(i)
    la_i = lambda_lyi(i)*0.25D-08
    nu_i = c/la_i/1.0D12
    tau_i = 0.0D0
    IF (xf.ge.nu_i*(1.0D0+zo)/(1.0D0+z).AND.xf.lt.nu_i) THEN
    xz = nu_i*(1.0D0+zo)/xf-1.0D0
      W=0.5D0*(1.0D0+DTANH((xz-zt)/dzt))    
      xHeI=xHI(xz)*(eta+(1.0D0-eta)*W)
     arHeIII=(xzHeIII-(1.0D0+xz)**1.5)/dyzHeIII
     xHeIII=0.5D0*(1.0D0+DTANH(arHeIII))
    xHeII=1.0D0-xHeI-xHeIII              
    tau_i = 7.0945D27*la_i*f_i*H0*Yp*Omb*(1.0D0+xz)**3*xHeII/&
    DSQRT(Omm*(1.0D0+xz)**3+OmL)   
    END IF
    tauLyHeII_v = tauLyHeII_v + tau_i
    END DO
    
  END FUNCTION tauLyHeII  
  
    FUNCTION tauLyHeI(xf, z, zo) RESULT(tauLyHeI_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, z, zo
    REAL(dp) :: tauLyHeI_val, nu_i, la_i, f_i, xz, ye, tau_i, xHeI, W 
    INTEGER :: i
    
    tauLyHeI_val=0.0D0
    DO i = 1,n_lyHeI
    f_i = fHeI_lyi(i)
    la_i = lambdaHeI_lyi(i)*1.0D-08
    nu_i = c/la_i/1.0D12
    tau_i = 0.0D0
    IF (xf.ge.nu_i*(1.0D0+zo)/(1.0D0+z).AND.xf.lt.nu_i) THEN
    xz = nu_i*(1.0D0+zo)/xf-1.0D0
    W=0.5D0*(1.0D0+DTANH((xz-zt)/dzt))    
    xHeI=xHI(xz)*(eta+(1.0D0-eta)*W)   
    tau_i = 7.0945D27*la_i*f_i*H0*Yp*Omb*(1.0D0+xz)**3*xHeI/DSQRT(Omm*(1.0D0+xz)**3+OmL)  
    END IF
    tauLyHeI_val = tauLyHeI_val + tau_i
    END DO
    
  END FUNCTION tauLyHeI    

  ! Function F
  FUNCTION F(xf, zp, zo) RESULT(F_val)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xf, zp, zo
    REAL(dp) :: F_val, tauz 
      tauz = tauLyHeI(xf, zp, zo)  + taucHeI(xf, zp, zo) + tauc(xf, zp, zo) + &
      tauGP(xf, zp, zo) + taucHeII(xf, zp, zo) + tauLyHeII(xf, zp, zo) + tauThom(zp, zo) 
!      tauz = tauGP(xf, zp, zo) + tauc(xf, zp, zo)
!      tauz = tauLyHeI(xf, zp, zo)  + taucHeI(xf, zp, zo)
!      tauz = taucHeII(xf, zp, zo) + tauLyHeII(xf, zp, zo)
      F_val = 0.0D0      
      IF (tauz.le.70.D0) THEN
      F_val = DEXP(-tauz) 
      END IF
   
  END FUNCTION F
 

  
END PROGRAM main
