!BOP
!\newpage
! !ROUTINE: cvmix_tke_driver

! FIXME: Update text
! !DESCRIPTION: A routine to test the Large, et al., implementation of shear
!  mixing. Inputs are the coefficients used in Equation (28) of the paper.
!  The diffusivity coefficient is output from a single column to allow
!  recreation of the paper's Figure 3. Note that here each "level" of the
!  column denotes a different local gradient Richardson number rather than a
!  physical ocean level. All memory is declared in the driver, and the CVMix
!  data type points to the local variables.
!\\
!\\

! !INTERFACE:

Subroutine cvmix_tke_driver(nlev, max_nlev)

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_zero,               &
                                    cvmix_one,                &
                                    cvmix_data_type
  ! FIXME delete shear
  use cvmix_shear,           only : cvmix_init_shear,         &
                                    cvmix_coeffs_shear
  use cvmix_tke,             only : init_tke,           &
                                    cvmix_coeffs_tke, &
                                    integrate_tke
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : save_variable!, write_simulation_info
  use cvmix_io,              only : cvmix_io_open,            &
                                    cvmix_output_write,       &
#ifdef _NETCDF
                                    cvmix_output_write_att,   &
#endif
                                    cvmix_io_close

  implicit none

! !INPUT PARAMETERS:
  integer, intent(in) :: nlev,               &! number of levels for column
                         max_nlev             ! number of columns in memory

!EOP
!BOC

  ! Global parameter
  integer, parameter :: ncol = 6

  ! CVMix datatypes
  !type(cvmix_data_type)                  :: CVmix_vars_tke_1D
  ! FIXME: delete
  !type(cvmix_data_type)                  :: CVmix_vars_LMD_1D, CVmix_vars_PP_1D
  !type(cvmix_data_type), dimension(ncol) :: CVmix_vars_PP_2D

  ! All regression tests look at Richardson numbers in [0,1]
  !real(cvmix_r8), dimension(:), allocatable, target :: Ri_g

  ! "1D" variables will be 2 x nlev+1 (1 for LMD, 1 for PP)
  ! "2D" variables will be ncol x nlev+1 (for PP)
  real(cvmix_r8), dimension(:), allocatable, target :: Mdiff_1D, Tdiff_1D
  !real(cvmix_r8), dimension(:,:), allocatable, target :: Mdiff_2D, Tdiff_2D

  ! Hard-code in parameters for each c in Table 1
  !real(cvmix_r8), dimension(ncol) :: PP_nu_zero_2D, PP_exp_2D, PP_alpha_2D
  !real(cvmix_r8)                  :: PP_nu_zero, PP_exp, PP_alpha

  !integer :: icol, kw, fid
  !integer :: fid
  integer :: kw

  integer :: i, j, tstep_count, ll, nt
  real(cvmix_r8), dimension(:), allocatable, target :: dzw, dzt
  real(cvmix_r8), dimension(:), allocatable, target :: Nsqr, Ssqr
  real(cvmix_r8), dimension(:), allocatable, target :: cvmix_dummy_1, &
    cvmix_dummy_2, cvmix_dummy_3
  real(cvmix_r8), dimension(:), allocatable, target :: tke_Tbpr, tke_Tspr, &
    tke_Tdif ,tke_Tdis, tke_Twin, tke_Tiwf, tke_Tbck, tke_Ttot
  real(cvmix_r8), dimension(:), allocatable, target :: tke, tke_Lmix, tke_Pr, &
    tke_iw_alpha_c, tke_plc
  real(cvmix_r8), dimension(:), allocatable, target :: Av_old, kv_old
  real(cvmix_r8), dimension(:), allocatable, target :: tke_iwe
  !real(cvmix_r8) :: forc_tke_surf_2D, bottom_fric_2D, forc_rho_surf_2D
  real(cvmix_r8) :: forc_tke_surf_2D, forc_rho_surf_2D
  real(cvmix_r8) :: dtime
  real(cvmix_r8) :: OceanReferenceDensity, grav

  ! Namelist variables
  ! KPP mixing parameters for column
  real(cvmix_r8) :: LMD_nu_zero, LMD_Ri_zero, LMD_exp

  namelist/tke_nml/LMD_nu_zero, LMD_Ri_zero, LMD_exp

  print*, "Active levels: ", nlev
  print*, "Levels allocated in memory: ", max_nlev

  tstep_count = 0
  nt = 10
  i = 1
  j = 1

  ! Allocate memory to store viscosity and diffusivity values
  allocate(Mdiff_1D(max_nlev+1), Tdiff_1D(max_nlev+1))
  !allocate(Mdiff_2D(ncol,max_nlev+1), Tdiff_2D(ncol,max_nlev+1))
  !allocate(Nsqr(1:max_nlev+1))
  !allocate(Ssqr(1:max_nlev+1))

  allocate(dzw(max_nlev), dzt(max_nlev+1))
  allocate(Nsqr(max_nlev+1), Ssqr(max_nlev+1))
  allocate(cvmix_dummy_1(max_nlev+1), &
    cvmix_dummy_2(max_nlev+1), cvmix_dummy_3(max_nlev+1))
  allocate(tke_Tbpr(max_nlev+1), tke_Tspr(max_nlev+1), &
    tke_Tdif(max_nlev+1) ,tke_Tdis(max_nlev+1), tke_Twin(max_nlev+1), &
    tke_Tiwf(max_nlev+1), tke_Tbck(max_nlev+1), tke_Ttot(max_nlev+1))
  allocate(tke(max_nlev+1), tke_Lmix(max_nlev+1), &
    tke_Pr(max_nlev+1), tke_iw_alpha_c(max_nlev+1), tke_plc(max_nlev+1))
  !allocate(tke_Twin(max_nlev+1), tke_Tiwf(max_nlev+1))
  allocate(Av_old(max_nlev+1), kv_old(max_nlev+1))
  allocate(tke_iwe(max_nlev+1))

  dzw(:) = 10.
  dzt(:) = 10.
  tke_iwe = 0.

  Nsqr(1) = 1e-6
  do kw=2,max_nlev+1
    Nsqr(kw) = Nsqr(kw-1) * 100./1000.
  enddo
  Ssqr(:) = 1e-4
  !Nsqr(:) = 0.0
  !Ssqr(:) = 0.0
  forc_tke_surf_2D = 1.0 

  CALL init_tke(tke_min=1d-6)
  print*, "After init_tke()"

  do ll = 1,nt
    tstep_count = tstep_count + 1
    !CALL cvmix_coeffs_tke( &
    CALL integrate_tke( &
        ! parameter
        !i = i,                             &
        !j = j,                             &
        !tstep_count  = tstep_count,        &
        dtime        = dtime,              &
        dzw          = dzw(:),             &
        dzt          = dzt(:),             &
        nlev         = max_nlev,           &
        max_nlev     = max_nlev,           &
        tke_old      = tke(:),             & ! in 
        tke_new      = tke(:),             & ! out
        !
        Ssqr         = Ssqr(:),            & ! in
        Nsqr         = Nsqr(:),            & ! in
        KappaM_out   = Mdiff_1D(:),        & ! out
        KappaH_out   = Tdiff_1D(:),        & ! out
        tke_Tbpr     = tke_Tbpr(:),        &
        tke_Tspr     = tke_Tspr(:),        &
        tke_Tdif     = tke_Tdif(:),        &
        tke_Tdis     = tke_Tdis(:),        &
        tke_Twin     = tke_Twin(:),        &
        tke_Tiwf     = tke_Tiwf(:),        &
        tke_Tbck     = tke_Tbck(:),        &
        tke_Ttot     = tke_Ttot(:),        &
        cvmix_int_1  = cvmix_dummy_1(:),   & !
        cvmix_int_2  = cvmix_dummy_2(:),   & !
        cvmix_int_3  = cvmix_dummy_3(:),   & !
        !
        tke_Lmix     = tke_Lmix(:),        &
        tke_Pr       = tke_Pr(:),          &
        forc_tke_surf= forc_tke_surf_2D,   &
        E_iw         = tke_iwe(:),            & ! for IDEMIX Ri
        !bottom_fric  = bottom_fric_2D,     &
        !old_kappaM   = Av_old(:),          & ! in
        !old_KappaH   = kv_old(:),          & ! in
        iw_diss      = tke_Tiwf(:),        & 
        forc_rho_surf= forc_rho_surf_2D,   &
        rho_ref      = OceanReferenceDensity, &
        grav         = grav,               &
        alpha_c      = tke_iw_alpha_c(:)  & ! for IDEMIX Ri
        )

    write(*,*) ' ------------- '
    write(*,*) 'll = ', ll
    write(*,*) 'tke = ', tke
    write(*,*) 'Nsqr = ', Nsqr
    write(*,*) 'Ssqr = ', Ssqr
    CALL save_variable('./out/', tke, 'tke', 'turbulent kinetic energy', ll, max_nlev+1, 'm2/s2')
    CALL save_variable('./out/', tke_Lmix, 'tke_Lmix', 'TKE mixing length', ll, max_nlev+1, 'm')
    CALL save_variable('./out/', tke_Pr, 'tke_Pr', 'TKE Prandtl number', ll, max_nlev+1, '')
    CALL save_variable('./out/', Mdiff_1D, 'Av', 'vertical viscosity', ll, max_nlev+1, 'm2/s')
    CALL save_variable('./out/', Tdiff_1D, 'kv', 'vertical diffusivity', ll, max_nlev+1, 'm2/s')
  enddo
  print*, "After integrate_tke()"
!
!  ! Initialization for 1D TKE test
!  call cvmix_put(CVmix_vars_TKE_1D, 'nlev',     nlev)
!  call cvmix_put(CVmix_vars_TKE_1D, 'max_nlev', max_nlev)
!  CVmix_vars_TKE_1D%tke => tke(:)
!  CVmix_vars_TKE_1D%tke_Lmix => tke_Lmix(:)
!  CVmix_vars_TKE_1D%tke_Pr => tke_Pr(:)
!  write(*,*) CVmix_vars_TKE_1D%tke(:)

!  ! Point CVmix_vars values to memory allocated above
!  CVmix_vars_TKE%Mdiff_iface           => Mdiff_1D(2,:)
!  CVmix_vars_TKE%Tdiff_iface           => Tdiff_1D(2,:)
!  CVmix_vars_TKE%ShearRichardson_iface => Ri_g

!   Set TKE single column parameters
!  !call cvmix_init_shear(mix_scheme='PP', PP_nu_zero=PP_nu_zero,               &
!  !                      PP_alpha=PP_alpha, PP_exp=PP_exp)
!  call cvmix_coeffs_shear(CVmix_vars_TKE)
!

  write(*,*) "All done!"
End Subroutine cvmix_tke_driver
