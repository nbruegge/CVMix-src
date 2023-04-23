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
                                    cvmix_data_type,          &
                                    cvmix_strlen
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
  real(cvmix_r8), dimension(:), allocatable, target :: Mdiff_1D, Tdiff_1D

  integer :: i, j, tstep_count, ll, llout, kw
  real(cvmix_r8), dimension(:), allocatable, target :: dzw, dzt
  real(cvmix_r8), dimension(:), allocatable, target :: Nsqr, Ssqr
  real(cvmix_r8), dimension(:), allocatable, target :: cvmix_dummy_1, &
    cvmix_dummy_2, cvmix_dummy_3
  real(cvmix_r8), dimension(:), allocatable, target :: tke_Tbpr, tke_Tspr, &
    tke_Tdif ,tke_Tdis, tke_Twin, tke_Tiwf, tke_Tbck, tke_Ttot
  real(cvmix_r8), dimension(:), allocatable, target :: tke, tke_Lmix, tke_Pr, &
    tke_iw_alpha_c, tke_plc
  real(cvmix_r8), dimension(:), allocatable, target :: tke_in, tke_out
  real(cvmix_r8), dimension(:), allocatable, target :: Av_old, kv_old
  real(cvmix_r8), dimension(:), allocatable, target :: tke_iwe
  !real(cvmix_r8) :: forc_tke_surf_2D, bottom_fric_2D, forc_rho_surf_2D
  real(cvmix_r8) :: forc_tke_surf_2D, forc_rho_surf_2D
  real(cvmix_r8) :: OceanReferenceDensity, grav

  ! Namelist variables
  ! TKE mixing parameters
  real(cvmix_r8) :: dtime
  integer :: nt, nt_output
  character(len=cvmix_strlen) :: path_out

  ! read namelist
  namelist/tke_nml/dtime, nt, nt_output, path_out
  read(*, nml=tke_nml)

  print*, "Active levels: ", nlev
  print*, "Levels allocated in memory: ", max_nlev
  print*, "time step: ", dtime
  print*, "number of time steps: ", nt
  print*, "path_out: ", path_out

  ! Allocate memory to store viscosity and diffusivity values
  allocate(Mdiff_1D(max_nlev+1), Tdiff_1D(max_nlev+1))
  allocate(dzw(max_nlev), dzt(max_nlev+1))
  allocate(Nsqr(max_nlev+1), Ssqr(max_nlev+1))
  allocate(cvmix_dummy_1(max_nlev+1), &
    cvmix_dummy_2(max_nlev+1), cvmix_dummy_3(max_nlev+1))
  allocate(tke_Tbpr(max_nlev+1), tke_Tspr(max_nlev+1), &
    tke_Tdif(max_nlev+1) ,tke_Tdis(max_nlev+1), tke_Twin(max_nlev+1), &
    tke_Tiwf(max_nlev+1), tke_Tbck(max_nlev+1), tke_Ttot(max_nlev+1))
  allocate(tke(max_nlev+1), tke_Lmix(max_nlev+1), &
    tke_Pr(max_nlev+1), tke_iw_alpha_c(max_nlev+1), tke_plc(max_nlev+1))
  allocate(tke_in(max_nlev+1), tke_out(max_nlev+1))
  allocate(Av_old(max_nlev+1), kv_old(max_nlev+1))
  allocate(tke_iwe(max_nlev+1))

  ! Initialize variables
  llout = 0
  tstep_count = 0
  i = 1
  j = 1

  dzw(:) = 10.
  dzt(:) = 10.
  tke_iwe = 0.

  ! stratification
  Nsqr(1) = 1e-6
  do kw=2,max_nlev+1
    Nsqr(kw) = Nsqr(kw-1) * 0.92
  enddo
  Nsqr(1:10) = -1e-8

  ! shear
  Ssqr(:) = 1e-6
  do kw=2,max_nlev+1
    Ssqr(kw) = Ssqr(kw-1) * 0.92
  enddo

  ! surface focing
  forc_tke_surf_2D = 1.e-4

  CALL init_tke(tke_min=1d-6)
  print*, "After init_tke()"

  do ll = 1,nt
    write(*,*) "ll = ", ll
    tstep_count = tstep_count + 1
    llout = llout + 1
    CALL integrate_tke( &
        dtime        = dtime,              &
        dzw          = dzw(:),             &
        dzt          = dzt(:),             &
        nlev         = max_nlev,           &
        max_nlev     = max_nlev,           &
        tke_old      = tke_in(:),             & ! in 
        tke_new      = tke_out(:),             & ! out
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

    tke(:) = tke_out(:)
    tke_in(:) = tke_out(:)

    !write(*,*) ' ------------- '
    !write(*,*) 'll = ', ll
    !write(*,*) 'tke = ', tke

    if (llout == nt_output) then
      write(*,*) 'Writing out at ll = ', ll
      CALL save_variable(path_out, Ssqr, 'Ssqr', 'shear squared', ll, max_nlev+1, '1/s2')
      CALL save_variable(path_out, Nsqr, 'Nsqr', 'shear squared', ll, max_nlev+1, '1/s2')
      CALL save_variable(path_out, tke, 'tke', 'turbulent kinetic energy', ll, max_nlev+1, 'm2/s2')
      CALL save_variable(path_out, tke_Lmix, 'tke_Lmix', 'TKE mixing length', ll, max_nlev+1, 'm')
      CALL save_variable(path_out, tke_Pr, 'tke_Pr', 'TKE Prandtl number', ll, max_nlev+1, '')
      CALL save_variable(path_out, Mdiff_1D, 'Av', 'vertical viscosity', ll, max_nlev+1, 'm2/s')
      CALL save_variable(path_out, Tdiff_1D, 'kv', 'vertical diffusivity', ll, max_nlev+1, 'm2/s')
      CALL save_variable(path_out, tke_Ttot, 'tke_Ttot', 'tke tendency tot', ll, max_nlev+1, 'm2/s3')
      CALL save_variable(path_out, tke_Tbpr, 'tke_Tbpr', 'tke tendency bpr', ll, max_nlev+1, 'm2/s3')
      CALL save_variable(path_out, tke_Tspr, 'tke_Tspr', 'tke tendency spr', ll, max_nlev+1, 'm2/s3')
      CALL save_variable(path_out, tke_Tdif, 'tke_Tdif', 'tke tendency dif', ll, max_nlev+1, 'm2/s3')
      CALL save_variable(path_out, tke_Tdis, 'tke_Tdis', 'tke tendency dif', ll, max_nlev+1, 'm2/s3')
      CALL save_variable(path_out, tke_Tbck, 'tke_Tbck', 'tke tendency bck', ll, max_nlev+1, 'm2/s3')
      CALL save_variable(path_out, tke_Twin, 'tke_Twin', 'tke tendency win', ll, max_nlev+1, 'm2/s3')
      CALL save_variable(path_out, tke_Tiwf, 'tke_Tiwf', 'tke tendency iwf', ll, max_nlev+1, 'm2/s3')
      llout = 0
    endif
  enddo
  print*, "After integrate_tke()"

  write(*,*) "All done!"
End Subroutine cvmix_tke_driver
