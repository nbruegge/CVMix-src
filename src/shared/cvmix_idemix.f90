
module cvmix_idemix
! This module contains the main computations of the IDEMIX 1 parameterization (described in "A Global Model for the Diapycnal
! Diffusivity Induced by Internal Gravity Waves", Olbers&Eden 2013) of Internal wave energy and its dissipation
! 
! @see
!
! @see
!
!
!! @author Nils Brueggemann, University of Hamburg / MPIMET
!!
!! @par Copyright
!! 2002-2013 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.

use cvmix_kinds_and_types,    only : cvmix_r8,                     &
                                      CVMIX_OVERWRITE_OLD_VAL,     &
                                      CVMIX_SUM_OLD_AND_NEW_VALS,  &
                                      CVMIX_MAX_OLD_AND_NEW_VALS,  &
                                      cvmix_data_type,             &
                                      cvmix_PI,                    & 
                                      cvmix_global_params_type

use cvmix_utils,              only : solve_tridiag


!USE mo_exception,          ONLY: finish

implicit none
private 
save


!public member functions

public :: cvmix_init_idemix
public :: cvmix_coeffs_idemix
public :: cvmix_put_idemix
!public :: cvmix_coeffs_idemix_low
public :: calc_idemix_v0
public :: gofx2  ! fixme: only used by IDEMIX public?
public :: hofx1  ! fixme: public?
public :: hofx2  ! fixme: public?

!=================================================================================
!---------------------------------------------------------------------------------
! Interface to call the IDEMIX parameterization
!---------------------------------------------------------------------------------

interface cvmix_coeffs_idemix
    module procedure cvmix_coeffs_idemix_low
    module procedure cvmix_coeffs_idemix_wrap
end interface cvmix_coeffs_idemix

interface cvmix_put_idemix
    module procedure cvmix_put_idemix_int
    module procedure cvmix_put_idemix_real
end interface cvmix_put_idemix

!=================================================================================

  ! cvmix_idemix_params_type contains the necessary parameters for IDEMIX mixing
  type, public :: cvmix_idemix_params_type
    private

       real(cvmix_r8) :: tau_v       ! time scale for vertical symmetrisation (sec)
       real(cvmix_r8) :: tau_h       ! time scale for horizontal symmetrisation, only necessary for lateral diffusion (sec)
       real(cvmix_r8) :: gamma        ! constant of order one derived from the shape of the spectrum in m space (dimensionless)
       real(cvmix_r8) :: jstar        ! spectral bandwidth in modes (dimensionless)
       real(cvmix_r8) :: mu0          ! dissipation parameter (dimensionless)
      
       ! Flag for what to do with old values of CVmix_vars%[MTS]diff  
       integer :: handle_old_vals
  end type cvmix_idemix_params_type

  type(cvmix_idemix_params_type), target :: CVmix_idemix_params_saved

!CHARACTER(LEN=*), PARAMETER :: module_name = 'cvmix_idemix'

contains

!=================================================================================

subroutine cvmix_init_idemix(         &
  tau_v, &
  tau_h, &
  gamma, &
  jstar,  &
  mu0,    &
  old_vals,                           &
  CVmix_idemix_params_user            &
  )

    !character(len=*), optional, intent(in) :: mix_scheme,                     &
    !                                          old_vals
    character(len=*), optional, intent(in) :: old_vals

! !INPUT PARAMETERS:
    real(cvmix_r8),optional, intent(in) ::      &
      tau_v,                                    & 
      tau_h,                                    & 
      gamma,                                    &
      jstar,                                    &
      mu0

! !OUTPUT PARAMETERS:
    type(cvmix_idemix_params_type), optional, target, intent(inout) ::        &
                                               CVmix_idemix_params_user

!EOP
!BOC

    type(cvmix_idemix_params_type), pointer :: CVmix_idemix_params_out

    if (present(CVmix_idemix_params_user)) then
      CVmix_idemix_params_out => CVmix_idemix_params_user
    else
      CVmix_idemix_params_out => CVmix_idemix_params_saved
    end if


    if (present(tau_v)) then
      call cvmix_put_idemix('tau_v', tau_v,  CVmix_idemix_params_user)
    else
      call cvmix_put_idemix('tau_v',1.d0*86400.0 ,  CVmix_idemix_params_user)
    end if
    
    if (present(tau_h)) then
      call cvmix_put_idemix('tau_h', tau_h,  CVmix_idemix_params_user)
    else
      call cvmix_put_idemix('tau_h', 15.d0*86400.0,  CVmix_idemix_params_user)
    end if
    
    if (present(gamma)) then
      call cvmix_put_idemix('gamma', gamma,  CVmix_idemix_params_user)
    else
      call cvmix_put_idemix('gamma', 1.57d0,  CVmix_idemix_params_user)
    end if
    
    if (present(jstar)) then
      call cvmix_put_idemix('jstar', jstar,  CVmix_idemix_params_user)
    else
      call cvmix_put_idemix('jstar', 10.d0 ,  CVmix_idemix_params_user)
    end if

    if (present(mu0)) then
      call cvmix_put_idemix('mu0', mu0,  CVmix_idemix_params_user)
    else
      call cvmix_put_idemix('mu0', 4.d0/3.0 ,  CVmix_idemix_params_user)
    end if


    if (present(old_vals)) then
      select case (trim(old_vals))
        case ("overwrite")
          call cvmix_put_idemix('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,    &
                               cvmix_idemix_params_user)
        case ("sum")
          call cvmix_put_idemix('handle_old_vals', CVMIX_SUM_OLD_AND_NEW_VALS, &
                               cvmix_idemix_params_user)
        case ("max")
          call cvmix_put_idemix('handle_old_vals', CVMIX_MAX_OLD_AND_NEW_VALS, &
                               cvmix_idemix_params_user)
        case DEFAULT
          print*, "ERROR: ", trim(old_vals), " is not a valid option for ",   &
                  "handling old values of diff and visc."
          stop 1
      end select
    else
      call cvmix_put_idemix('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,        &
                               cvmix_idemix_params_user)
    end if

end subroutine cvmix_init_idemix

!=================================================================================

! !IROUTINE: cvmix_coeffs_idemix_wrap
! !INTERFACE:

  subroutine cvmix_coeffs_idemix_wrap(CVmix_vars,  CVmix_idemix_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the KPP boundary layer mixing
!  parameterization.
!
! FIXME: nils: this subroutine is never called but it is CVMix coding standard. 
! Here cvmix_coeffs_idemix is
! called directly. If this is used first debug it.

! FIXME: nils: Do these comments make sence for anyone?
! This subroutine is necessary to handle old/new values and to hand over the IDEMIX parameter set in previous subroutine
! This subroutine should be called from calling ocean model or driver
!
!\\
!\\
!
! !USES:
!  only those used by entire module.


    type(cvmix_idemix_params_type), intent(in), optional, target ::           &
                                              CVmix_idemix_params_user

! !INPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    real(cvmix_r8), dimension(CVmix_vars%max_nlev+1) ::       & 
      new_E_iw                                          ,& 
      cvmix_int_1                                       ,&
      cvmix_int_2                                       ,&
      cvmix_int_3                                       ,&
      iwe_Ttot                                          ,&
      iwe_Tdif                                          ,&
      !iwe_Thdi                                          ,&
      iwe_Tdis                                          ,&
      iwe_Tsur                                          ,&
      iwe_Tbot                                          ,&
      new_KappaM                                        ,&
      new_KappaH                                        ,&
      c0                                                ,&
      v0!                                                ,&
      !new_iw_diss                                       ! 
    
    integer ::                                           &
      nlev                                              ,&
      max_nlev                                          !
    ! FIXME: nils: for debugging delet later
    !integer :: i,j, tstep_count
    !logical :: debug
    
    type(cvmix_idemix_params_type), pointer :: idemix_constants_in
    
    !CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_coeffs_idemix_wrap'
    
    idemix_constants_in => CVmix_idemix_params_saved
    if (present( CVmix_idemix_params_user)) then
      idemix_constants_in =>  CVmix_idemix_params_user
    end if
    nlev = CVmix_vars%nlev
    max_nlev = CVmix_vars%max_nlev
    
    
    ! FIXME: nils: put a security stop in case this routine is called in the current state
    ! call to actual computation of IDEMIX parameterization
    !write(*,*) 'I am wrapping'
    !stop
    !CALL finish(method_name,'I am wrapping')
    
    call cvmix_coeffs_idemix( &
                             ! parameter
                             dtime           = CVmix_vars%dtime,           &
                             dzw             = CVmix_vars%dzw,             &
                             dzt             = CVmix_vars%dzt,             &
                             nlev            = nlev,                      &
                             max_nlev        = max_nlev,                  &
                             coriolis        = CVmix_vars%coriolis,        &
                             ! essentials
                             iwe_old         = CVmix_vars%E_iw,            &
                             iwe_new         = new_E_iw,                  &
                             forc_iw_surface = CVmix_vars%forc_iw_surface, &
                             forc_iw_bottom  = CVmix_vars%forc_iw_bottom,  &
                             ! FIXME: nils: better output IDEMIX Ri directly
                             alpha_c         = CVmix_vars%alpha_c,         &
                             ! only for Osborn shortcut 
                             ! FIXME: nils: put this to cvmix_tke
                             KappaM_out      = new_KappaM,                &
                             KappaH_out      = new_KappaH,                &
                             Nsqr            = CVmix_vars%Nsqr_iface,      &
                             ! diagnostics
                             iwe_Ttot        = iwe_Ttot,                  &
                             iwe_Tdif        = iwe_Tdif,                  &
                             !iwe_Thdi        = iwe_Thdi,                  &
                             iwe_Tdis        = iwe_Tdis,                  &
                             iwe_Tsur        = iwe_Tsur,                  &
                             iwe_Tbot        = iwe_Tbot,                  &
                             c0              = c0,                        &
                             v0              = v0,                        &
                             ! debugging
                             !debug = debug,             & ! FIXME: nils: for debuging
                             !i = i,                     & ! FIXME: nils: for debuging
                             !j = j,                     & ! FIXME: nils: for debuging
                             !tstep_count = tstep_count, & ! FIXME: nils: for debuging
                             cvmix_int_1     = cvmix_int_1,               &
                             cvmix_int_2     = cvmix_int_2,               &
                             cvmix_int_3     = cvmix_int_3,               &
                             CVmix_idemix_params_user =  CVmix_idemix_params_user)
    
    ! FIXME: nils: This should probably be cvmix_update_idemix. However, it is not used
    ! anyway.
    ! update CVmix_vars to new values
    !call cvmix_update_tke(idemix_constants_in%handle_old_vals,            &
    !                      nlev,                                           &
    !                      iw_diss_out = CVmix_vars%iw_diss,                &
    !                      new_iw_diss = new_iw_diss,                      &
    !                      iwe_new     = CVmix_vars%E_iw,                   &
    !                      new_E_iw    = new_E_iw)

  end subroutine cvmix_coeffs_idemix_wrap

subroutine calc_idemix_v0(nlev, max_nlev, Nsqr, dzw, coriolis, &
                          v0,  CVmix_idemix_params_user)
  integer, intent(in) ::                                          &
    nlev, max_nlev                                                !,&

  real(cvmix_r8), intent(in)                              ::      & 
    coriolis                                                        !

  !logical, intent(in) :: debug

  real(cvmix_r8), dimension(max_nlev+1), intent(in) ::                &
    dzw
  
  real(cvmix_r8), dimension(max_nlev+1), intent(in)           ::      &
    Nsqr                                                 !,&

  real(cvmix_r8), dimension(max_nlev+1), intent(out) ::               &
    v0                                                           !,&

  ! IDEMIX namelist parameters
  real(cvmix_r8)                                          ::      & 
    cstar                                                        ,& ! 
    tau_h                                                        ,& !
    gamma                                                        ,& !
    jstar                                                        ,& !
    mu0                                                          ,& !
    bN0                                                          !,&

  real(cvmix_r8)                                          ::      & 
    fxa 

  integer                                                 ::      &
    k


  type(cvmix_idemix_params_type), intent(in), optional, target ::  CVmix_idemix_params_user
  type(cvmix_idemix_params_type), pointer :: idemix_constants_in

  ! FIXME: nils: Is this necessary?
  idemix_constants_in => CVmix_idemix_params_saved
  if (present( CVmix_idemix_params_user)) then
    idemix_constants_in =>  CVmix_idemix_params_user
  end if

  ! set idemix_constants locally
  tau_h = idemix_constants_in%tau_h
  gamma = idemix_constants_in%gamma
  mu0   = idemix_constants_in%mu0
  jstar = idemix_constants_in%jstar
 
  ! calculate cstar from OE13 Eq. (13)
  bN0=0.0
  do k=2,nlev
    bN0 = bN0 + max(0.0_cvmix_r8,Nsqr(k))**0.5*dzw(k) 
  enddo
  cstar = max(1e-2_cvmix_r8,bN0/(cvmix_PI*jstar) )
     
  ! calculate horizontal representative group velocity v0
  ! v0: OE13 Eq. (A9)
  do k=1,nlev+1
    fxa = max(0.0_cvmix_r8,Nsqr(k))**0.5/(1d-22 + abs(coriolis) )
    v0(k)=max(0.0_cvmix_r8, gamma*cstar*hofx2(fxa))

    ! set v0 to zero to prevent horizontal iwe propagation in mixed layer
    if ( fxa<1.0_cvmix_r8 ) then
      v0(k) = 0.0_cvmix_r8
    endif

    !! for debugging:
    !if (debug .eqv. .true.) then
    !  write(*,*) "Nsqr = ", Nsqr(k) 
    !  write(*,*) "gamma = ", gamma
    !  write(*,*) "fxa = ", fxa
    !  write(*,*) "cstar = ", cstar 
    !  write(*,*) "hofx2(fxa) = ", hofx2(fxa)
    !  write(*,*) 'v0 = ', v0(k)
    !end if
  enddo
  !v0 = min(3d-1, v0)
end subroutine calc_idemix_v0

!=================================================================================
! This subroutine contains the actual computation of IDEMIX
subroutine cvmix_coeffs_idemix_low( &
                            ! parameter
                            dtime,                 &
                            dzw,                   &
                            dzt,                   &
                            nlev,                  &
                            max_nlev,              &
                            coriolis,              &
                            ! essentials
                            iwe_old,               & ! in
                            iwe_new,               & ! out
                            forc_iw_surface,       & ! in
                            forc_iw_bottom,        & ! in
                            ! FIXME: nils: better output IDEMIX Ri directly
                            alpha_c,               & ! out
                            ! only for Osborn shortcut
                            ! FIXME: nils: put this to cvmix_tke
                            KappaM_out,            & ! FIXME: nils: put to tke?
                            KappaH_out,            & ! FIXME: nils: put to tke?
                            Nsqr,                  & ! FIXME: nils: put to tke?
                            ! diagnostics
                            iwe_Ttot,              & ! diagnostic
                            iwe_Tdif,              & ! diagnostic
                            !iwe_Thdi,              & ! diagnostic
                            iwe_Tdis,              & ! diagnostic
                            iwe_Tsur,              & ! diagnostic
                            iwe_Tbot,              & ! diagnostic
                            c0,                    &
                            v0,                    &
                            ! debugging
                            !debug,                 & ! FIXME: nils: for debuging
                            !i,                     & ! FIXME: nils: for debuging
                            !j,                     & ! FIXME: nils: for debuging
                            !tstep_count,           & ! FIXME: nils: for debuging
                            cvmix_int_1,           & ! FIXME: nils: for debuging
                            cvmix_int_2,           & ! FIXME: nils: for debuging
                            cvmix_int_3,           & ! FIXME: nils: for debuging
                             CVmix_idemix_params_user &
                            )

  
   type(cvmix_idemix_params_type), intent(in), optional, target ::  CVmix_idemix_params_user
  
   integer, intent(in) ::                                          &
     nlev                                                         ,&
     max_nlev                                                         
  
   real(cvmix_r8), dimension(max_nlev+1), intent(inout) ::             &
      KappaM_out                                                  ,&
      KappaH_out
  
   ! FIXME: nils: for debuging
   !integer, intent(in) :: i, j, tstep_count
   !logical, intent(in) :: debug
  
   real(cvmix_r8), dimension(max_nlev), intent(in) ::                &
     dzw
  
   real(cvmix_r8), dimension(max_nlev+1), intent(in)           ::      &
     Nsqr                                                         ,&
     iwe_old                                                      ,&
     !old_iw_diss                                                  ,& 
     dzt                                                             !
  
   ! diagnostics
   real(cvmix_r8), dimension(max_nlev+1), intent(out) ::               &
     !iw_diss_out                                                  ,& 
     iwe_new                                                      ,&
     cvmix_int_1                                                  ,&
     cvmix_int_2                                                  ,&
     cvmix_int_3                                                  ,&
     iwe_Ttot                                                     ,&
     iwe_Tdif                                                     ,&
     iwe_Tdis                                                     ,&
     iwe_Tsur                                                     ,&
     iwe_Tbot                                                     ,&
     c0                                                           ,&
     v0                                                           ,&
     alpha_c

   !real(cvmix_r8), dimension(max_nlev+1), intent(in) ::               &
   !  iwe_Thdi
  
  real(cvmix_r8), intent(in)                              ::      & 
    forc_iw_bottom                                               ,& !
    forc_iw_surface                                              ,& !
    dtime                                                        ,& !
    coriolis                                                        !
 
  integer                                                 ::      &
    k
 
  ! coefficients for the tri-diagonal solver
  real(cvmix_r8), dimension(max_nlev+1)                       ::      &
    a_dif                                                        ,& !
    b_dif                                                        ,& !
    c_dif                                                        ,& !
    a_tri                                                        ,& !
    b_tri                                                        ,& !
    c_tri                                                        ,& !
    d_tri
 
  real(cvmix_r8), dimension(max_nlev+1)                       ::      &
    delta                                                        ,& !
    iwe_max                                                      ,& ! 
    forc                                                            ! 
 
  ! IDEMIX namelist parameters
  real(cvmix_r8)                                          ::      & 
    cstar                                                        ,& ! 
    tau_v                                                        ,& !
    tau_h                                                        ,& !
    gamma                                                        ,& !
    jstar                                                        ,& !
    mu0                                                          ,& !
    bN0                                                             !
 
  real(cvmix_r8)                                          ::      & 
    fxa 
 
  type(cvmix_idemix_params_type), pointer ::idemix_constants_in

  ! initialize variables
  iwe_new     = 0.0_cvmix_r8
  cvmix_int_1 = 0.0_cvmix_r8
  cvmix_int_2 = 0.0_cvmix_r8
  cvmix_int_3 = 0.0_cvmix_r8
  iwe_Ttot    = 0.0_cvmix_r8
  iwe_Tdif    = 0.0_cvmix_r8
  iwe_Tdis    = 0.0_cvmix_r8
  iwe_Tsur    = 0.0_cvmix_r8
  iwe_Tbot    = 0.0_cvmix_r8
  c0          = 0.0_cvmix_r8
  v0          = 0.0_cvmix_r8
  alpha_c     = 0.0_cvmix_r8
  a_dif       = 0.0_cvmix_r8
  b_dif       = 0.0_cvmix_r8
  c_dif       = 0.0_cvmix_r8
  a_tri       = 0.0_cvmix_r8
  b_tri       = 0.0_cvmix_r8
  c_tri       = 0.0_cvmix_r8
  d_tri       = 0.0_cvmix_r8
  delta       = 0.0_cvmix_r8
  iwe_max     = 0.0_cvmix_r8
  forc        = 0.0_cvmix_r8
 
  ! FIXME: nils: Is this necessary?
  idemix_constants_in => CVmix_idemix_params_saved
  if (present( CVmix_idemix_params_user)) then
    idemix_constants_in =>  CVmix_idemix_params_user
  end if
 
  ! set idemix_constants locally
  tau_v = idemix_constants_in%tau_v
  tau_h = idemix_constants_in%tau_h
  gamma = idemix_constants_in%gamma
  mu0   = idemix_constants_in%mu0
  jstar = idemix_constants_in%jstar
 
  ! calculate cstar from OE13 Eq. (13)
  bN0=0.0_cvmix_r8
  do k=2,nlev
    bN0 = bN0 + max(0.0_cvmix_r8,Nsqr(k))**0.5*dzw(k) 
  enddo
  cstar = max(1e-2_cvmix_r8,bN0/(cvmix_PI*jstar) )
     
  ! calculate vertical and horizontal representative group velocities c0 and v0
  ! c0: OE13 Eq. (13) 
  ! alpha_c iwe**2: dissipation of internal wave energy (OE13 Eq. (15))
  do k=1,nlev+1
    fxa = max(0.0_cvmix_r8,Nsqr(k))**0.5/(1e-22_cvmix_r8 + abs(coriolis) )
    c0(k)=max(0.0_cvmix_r8, gamma*cstar*gofx2(fxa) )
    v0(k)=max(0.0_cvmix_r8, gamma*cstar*hofx2(fxa))
    !v0(k)=0.5
    alpha_c(k) = max( 1e-4_cvmix_r8, mu0*acosh(max(1.0_cvmix_r8,fxa))*abs(coriolis)/cstar**2 )

    ! set v0 to zero to prevent horizontal iwe propagation in mixed layer
    if ( fxa<1.0_cvmix_r8 ) then
      v0(k) = 0.0_cvmix_r8
    endif
  enddo
 
  !---------------------------------------------------------------------------------
  ! initialize forcing
  forc(:)=0.0_cvmix_r8
 
  ! add tendency of horizontal diffusion (is calculated externally)
  !forc(:) = forc(:) + iwe_Thdi(:)
 
  !---------------------------------------------------------------------------------
  ! prevent negative dissipation of IW energy
  ! FIXME: Carsten thinks we don't need this
  iwe_max = max(0.0_cvmix_r8, iwe_old)
 
 
  ! vertical diffusion and dissipation is solved implicitely 
  !---------------------------------------------------------------------------------
  ! assignment of tridiagonal matrix
  !---------------------------------------------------------------------------------
  ! |b1 c1 0  0  0  | (E1) = (d1)
  ! |a2 b2 c2 0  0  | (E2) = (d2)
  ! |0  a3 b3 c3 0  | (E3) = (d3)
  ! |0  0  a4 b4 c4 | (E4) = (d4)
  ! |0  0  0  an bn | (En) = (dn)
  !
  ! d1 = diss_1 + surf_forc 
  ! dn = diss_n + bott_forc 
  ! 
    
  ! vertical flux
  do k=1,nlev
   delta(k) = tau_v/dzw(k) * 0.5*(c0(k)+c0(k+1))
  enddo
  delta(nlev+1) = 0.0          ! delta(nlev+1) is never used
 
  ! -- a -- 
  do k=2,nlev+1
    a_dif(k) = delta(k-1)*c0(k-1)/dzt(k)
  enddo
  a_dif(1) = 0.0 ! not part of the diffusion matrix, thus value is arbitrary
 
  ! -- b -- 
  do k=2,nlev
    b_dif(k) = (delta(k-1)*c0(k)+delta(k)*c0(k))/dzt(k)
  enddo
 
  ! Neumann boundary conditions
  k = 1
  b_dif(k) = delta(k)*c0(k)/dzt(k)
  k = nlev+1
  b_dif(k) = delta(k-1)*c0(k)/dzt(k)
 
  ! -- c-- 
  do k=1,nlev
    c_dif(k) = delta(k)*c0(k+1)/dzt(k)
  enddo
  c_dif(nlev+1) = 0.0 ! not part of the diffusion matrix, thus value is arbitrary
 
  !--- construct tridiagonal matrix to solve diffusion and dissipation implicitely
  a_tri = -dtime*a_dif
  b_tri = 1+dtime*b_dif
  ! FIXME: nils: Should dissipation also be in first and last layer?
  b_tri(2:nlev) = b_tri(2:nlev) + dtime*alpha_c(2:nlev)*iwe_max(2:nlev)
  c_tri = -dtime*c_dif
   
  ! -- d -- 
  d_tri(1:nlev+1) = iwe_old(1:nlev+1) + dtime*forc(1:nlev+1)
  d_tri(nlev+1)   = d_tri(nlev+1)      + dtime*forc_iw_bottom/dzt(nlev+1) 
  d_tri(1)        = d_tri(1)           + dtime*forc_iw_surface/dzt(1)
 
  ! solve the tri-diag matrix 
  call solve_tridiag(a_tri, b_tri, c_tri, d_tri, iwe_new, nlev+1)
 
  ! --- diagnose implicite tendencies (only for diagnostics)
  ! vertical diffusion of E_iw
  do k=2,nlev
    iwe_Tdif(k) = a_dif(k)*iwe_new(k-1) - b_dif(k)*iwe_new(k) + c_dif(k)*iwe_new(k+1)
  enddo
  k = 1
  iwe_Tdif(k) = - b_dif(k)*iwe_new(k) + c_dif(k)*iwe_new(k+1)
  k = nlev+1
  iwe_Tdif(k) = a_dif(k)*iwe_new(k-1) - b_dif(k)*iwe_new(k)
 
  ! dissipation of E_iw
  iwe_Tdis = 0.0
  ! FIXME: nils: dissipation also in first or last layer?
  !iwe_Tdis(1:nlev+1) =  -alpha_c(1:nlev+1) * iwe_max(1:nlev+1) * iwe_new(1:nlev+1)
  iwe_Tdis(2:nlev) =  -alpha_c(2:nlev) * iwe_max(2:nlev) * iwe_new(2:nlev)
 
  iwe_Tsur(1)      = forc_iw_surface/dzt(1) 
  iwe_Tbot(nlev+1) = forc_iw_bottom/dzt(nlev+1)
 
  iwe_Ttot = (iwe_new-iwe_old)/dtime
 
  ! IDEMIX only shortcut: derive diffusivity and viscosity using Osbourne relation
  KappaH_out = 0.0
  KappaM_out = 0.0
  do k=2,nlev
    KappaH_out(k) =  0.2/(1.0+0.2) * (-1.0*iwe_Tdis(k)) / max(1e-12_cvmix_r8, Nsqr(k))
    KappaH_out(k) = max(1e-9_cvmix_r8, KappaH_out(k))
    KappaH_out(k) = min(1.0_cvmix_r8, KappaH_out(k))
    KappaM_out(k) =  10.0 * KappaH_out(k)
  enddo
 
  !---------------------------------------------------------------------------------
  ! rest is for debuggin only
  !---------------------------------------------------------------------------------
  cvmix_int_1 = Nsqr
  cvmix_int_2 = alpha_c 
  cvmix_int_3 = c0
 
  ! debugging: 
  !if (debug .eqv. debug) then
!  if (.false.) then
  !if (i==45 .and. j==10) then
  !if (i==45 .and. j==45) then
!     write(*,*) ' ===================== '
 
!     write(*,*) 'dtime = ', dtime
!     write(*,*) 'delta = ', delta
!     write(*,*) 'dzw = ', dzw
!     write(*,*) 'c0 = ', c0
!     write(*,*) 'a_tri = ', a_tri
!     write(*,*) 'b_tri = ', b_tri
!     write(*,*) 'c_tri = ', c_tri
!     write(*,*) 'd_tri = ', d_tri
!     write(*,*) 'forc_iw_surface = ', forc_iw_surface
!     write(*,*) 'iwe_new = ', iwe_new
 
!     write(*,*) 'iwe_Ttot = ', iwe_Ttot
!     write(*,*) 'iwe_Tdif = ', iwe_Tdif
!     write(*,*) 'iwe_Thdi = ', iwe_Thdi
!     write(*,*) 'iwe_Tdis = ', iwe_Tdis
!     write(*,*) 'iwe_Tsur = ', iwe_Tsur
!     write(*,*) 'iwe_Tbot = ', iwe_Tbot
!     write(*,*) 'iwe_Tres = ', iwe_Ttot-(iwe_Tdif+iwe_Thdi+iwe_Tdis+iwe_Tsur+iwe_Tbot)
 
!    write(*,*) 'tau_v = ', tau_v
!    write(*,*) 'tau_h = ', tau_h
!    write(*,*) 'gamma = ', gamma
!    write(*,*) 'jstar = ', jstar
!    write(*,*) 'mu0 = ', mu0
 
    !stop
  !endif
  !endif
!  endif
 
end subroutine cvmix_coeffs_idemix_low

!=================================================================================


function gofx2(x1)
!=======================================================================
! a function g(x)	! adapted from pyOM 
!=======================================================================
 implicit none
 real(cvmix_r8) :: gofx2,x1,x2,c
 real(cvmix_r8), parameter :: pi = 3.14159265358979323846264338327950588
 x2=max(3.0_cvmix_r8,x1)
 c= 1.-(2./pi)*asin(1./x2)
 gofx2 = 2/pi/c*0.9*x2**(-2./3.)*(1-exp(-x2/4.3))
end function gofx2

function hofx2(x1)
!=======================================================================
! a function h(x) 	! adapted from pyOM
!=======================================================================
 implicit none
 real(cvmix_r8) :: hofx2,x1,x2
 real(cvmix_r8), parameter :: pi = 3.14159265358979323846264338327950588
 x2 = max(10.0_cvmix_r8, x1) ! by_nils: it has to be x2>1
 hofx2 = (2./pi)/(1.-(2./pi)*asin(1./x2)) * (x2-1.)/(x2+1.)
end function hofx2

function hofx1(x)
!=======================================================================
! a function h(x) 	!from pyOM
!=======================================================================
 implicit none
 real(cvmix_r8) :: hofx1,x
 real(cvmix_r8), parameter :: pi = 3.14159265358979323846264338327950588
 hofx1 = (2./pi)/(1.-(2./pi)*asin(1./x)) * (x-1.)/(x+1.)
end function hofx1

!=================================================================================


!=================================================================================

subroutine cvmix_put_idemix_int(varname,val, CVmix_idemix_params_user)

! This subroutine puts integer values to IDEMIX variables
!IN
  character(len=*),           intent(in) :: varname
  integer,                    intent(in) :: val
!OUT   
  type(cvmix_idemix_params_type), intent(inout), target, optional::  CVmix_idemix_params_user
  type(cvmix_idemix_params_type), pointer :: idemix_constants_out

  idemix_constants_out=>CVmix_idemix_params_saved
  if (present( CVmix_idemix_params_user)) then
    idemix_constants_out=>  CVmix_idemix_params_user
  end if

  select case(trim(varname))

    ! FIXME: can be deleted
    case('handle_old_vals')
      idemix_constants_out%handle_old_vals=val
    
  end select
    
end subroutine cvmix_put_idemix_int

subroutine cvmix_put_idemix_real(varname, val, CVmix_idemix_params_user)

! !DESCRIPTION:
!  Write a real value into a cvmix\_idemix\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_idemix_params_type), intent(inout), target, optional ::           &
                                              CVmix_idemix_params_user

!EOP
!BOC

    type(cvmix_idemix_params_type), pointer :: CVmix_idemix_params_out

    CVmix_idemix_params_out => CVmix_idemix_params_saved
    if (present(CVmix_idemix_params_user)) then
      CVmix_idemix_params_out => CVmix_idemix_params_user
    end if

  select case(trim(varname))
      case ('tau_v')
        CVmix_idemix_params_out%tau_v = val
      case ('tau_h')
        CVmix_idemix_params_out%tau_h = val
      case ('gamma')
        CVmix_idemix_params_out%gamma = val
      case ('jstar')
        CVmix_idemix_params_out%jstar = val
      case ('mu0')
        CVmix_idemix_params_out%mu0 = val
      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
  end select
    
end subroutine cvmix_put_idemix_real

!=================================================================================

end module cvmix_idemix 
