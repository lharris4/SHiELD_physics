module topo_drag_driver_mod

!=======================================================================
! TOPOGRAPHIC DRAG CLOSURE -- Garner (2005)
! Driver for compatibility with the blocking in SHiELD_physics for input
!  (output will be handled on original grid)
! This should be called from the FV3GFS layer.
!=======================================================================

!-----------------------------------------------------------------------
!  Calculates horizontal velocity tendency due to topographic drag
!-----------------------------------------------------------------------

use          mpp_mod, only: input_nml_file
use  mpp_domains_mod, only: domain2d
use          fms_mod, only: error_mesg, FATAL, NOTE, &
                            mpp_pe, mpp_root_pe, stdout, stdlog,       &
                            check_nml_error, write_version_number
use      fms2_io_mod, only: read_data, get_variable_size, variable_exists, file_exists, &
                            FmsNetcdfDomainFile_t, register_variable_attribute, &
                            register_restart_field, register_axis, unlimited, &
                            open_file, read_restart, write_restart, close_file, &
                            register_field, write_data, get_global_io_domain_indices, &
                            FmsNetcdfFile_t
use    constants_mod, only: Radian
use horiz_interp_mod, only: horiz_interp_type, horiz_interp_init, &
                            horiz_interp_new, horiz_interp, horiz_interp_del
use    topo_drag_mod, only: frcrit, alin, anonlin, beta, gamma, epsi,               & !namelist items
                            h_frac, zref_fac, tboost, pcut, samp, max_udt,          &
                            no_drag_frac, max_pbl_frac,                             &
                            do_conserve_energy, keep_residual_flux, do_pbl_average, &
                            use_mg_scaling, use_mask_for_pbl, ricrit,               &
                            use_pbl_from_vert_turb, use_uref_4stable
use    GFS_typedefs,  only: GFS_sfcprop_type
use block_control_mod,only: block_control_type
use physcons,         only: Grav => con_g,  Cp_Air => con_cp,  Pi => con_pi, &
                            Rdgas=> con_rd, cvap  => con_cvap
use IPD_typedefs,     only: IPD_control_type
use atmosphere_mod,   only: atmosphere_grid_bdry

implicit none

private

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

logical :: module_is_initialized = .false.

! horizontal array size

integer :: nlon, nlat
integer :: kd=0

! arrays defined by topo_drag_init:

real, allocatable, dimension(:,:) :: t11, t21, t12, t22    ! drag tensor
real, allocatable, dimension(:,:) :: hmin, hmax

! parameters:

real, parameter :: u0=1.0       ! arbitrary velocity scale for diagnostics
real, parameter :: xl=80.0e3    ! arbitrary horiz length scale for diagnostics
real, parameter :: ro=1.2       ! arbitrary density scale for diagnostics
real, parameter :: lapse=Grav/Cp_Air ! adiabatic temperature lapse rate
real, parameter :: tiny=1.0e-20

real, parameter :: resolution=60.0 ! # of points per degree in topo datasets
real, parameter :: frint=0.5

integer, parameter :: ipts=360*resolution
integer, parameter :: jpts=180*resolution

type(domain2d), pointer     :: topo_domain


NAMELIST /topo_drag_nml/                                               &
  frcrit, alin, anonlin, beta, gamma, epsi,                            &
  h_frac, zref_fac, tboost, pcut, samp, max_udt,                       &
  no_drag_frac, max_pbl_frac,                                          &
  do_conserve_energy, keep_residual_flux, do_pbl_average,              &
  use_mg_scaling, use_mask_for_pbl, use_pbl_from_vert_turb,            &    !stg
  use_uref_4stable, ricrit

public topo_drag_init, topo_drag_end
public topo_drag_restart

contains

!#######################################################################

!=======================================================================
 subroutine topo_drag_register_tile_restart(restart)

     type(FmsNetcdfDomainFile_t), intent(inout) :: restart
     character(len=8), dimension(3)             :: dim_names

     dim_names(1) = "xaxis_1"
     dim_names(2) = "yaxis_1"
     dim_names(3) = "Time"
     call register_axis(restart, dim_names(1), "x")
     call register_axis(restart, dim_names(2), "y")
     call register_axis(restart, dim_names(3), unlimited)

     !< Register the domain decomposed dimensions as variables so that the combiner can work
     !! correctly
     call register_field(restart, dim_names(1), "double", (/dim_names(1)/))
     call register_field(restart, dim_names(2), "double", (/dim_names(2)/))

     call register_restart_field(restart, "t11", t11, dim_names)
     call register_restart_field(restart, "t12", t12, dim_names)
     call register_restart_field(restart, "t21", t21, dim_names)
     call register_restart_field(restart, "t22", t22, dim_names)
     call register_restart_field(restart, "hmin", hmin, dim_names)
     call register_restart_field(restart, "hmax", hmax, dim_names)

end subroutine

subroutine topo_drag_init (domain, Sfcprop, Atm_block, Model, enforce_rst_cksum)

type(domain2D), target, intent(in) :: domain
type(GFS_sfcprop_type), intent(in) :: Sfcprop(:)
type(block_control_type), intent(in)    :: Atm_block
type(IPD_control_type),   intent(inout) :: Model
logical, intent(in) :: enforce_rst_cksum

character(len=128) :: msg
character(len=64)  :: restart_fname='INPUT/topo_drag.res.nc'
character(len=64)  :: topography_file='INPUT/poztopog.nc'
character(len=64)  :: dragtensor_file='INPUT/dragelements.nc'
character(len=3)   :: tensornames(4) = (/ 't11', 't21', 't12', 't22' /)

logical :: found_field(4)

real, parameter :: bfscale=1.0e-2      ! buoyancy frequency scale [1/s]

real, allocatable, dimension(:)   :: xdatb, ydatb
real, allocatable, dimension(:,:) :: zdat, zout
real, allocatable, dimension(:,:) :: lonb_loc, latb_loc
type (horiz_interp_type) :: Interp
real :: exponent, hmod

integer :: isc, jsc
integer :: n
integer :: io, ierr, unit_nml, logunit
integer :: i, j
integer :: nb, ix
integer :: siz(4)
type(FmsNetcdfDomainFile_t) :: Topo_restart !< Fms2io domain decomposed fileobj
type(FmsNetcdfFile_t) :: topography_fileobj, dragtensor_fileobj !< Fms2io fileobj

  if (module_is_initialized) return

  isc = Model%isc
  jsc = Model%jsc

  nlon = Model%nx
  nlat = Model%ny

! read namelist

   read (input_nml_file, nml=topo_drag_nml, iostat=io)
   ierr = check_nml_error(io,'topo_drag_nml')

! write version number and namelist to logfile

  call write_version_number (version, tagname)
  logunit = stdlog()
  if (mpp_pe() == mpp_root_pe())                                       &
                                    write (logunit, nml=topo_drag_nml)

  topo_domain => domain
  allocate (t11(nlon,nlat))
  allocate (t21(nlon,nlat))
  allocate (t12(nlon,nlat))
  allocate (t22(nlon,nlat))
  allocate (hmin(nlon,nlat))
  allocate (hmax(nlon,nlat))

  if (gamma == beta + epsi) gamma = gamma + tiny

! read restart file

  if ( open_file(Topo_restart, restart_fname, "read", topo_domain, is_restart = .true.) ) then

     if (mpp_pe() == mpp_root_pe()) then
        write ( msg, '("Reading restart file: ",a40)' ) restart_fname
        call error_mesg('topo_drag_mod', msg, NOTE)
     endif

     call topo_drag_register_tile_restart(Topo_restart)
     call read_restart(Topo_restart, ignore_checksum=enforce_rst_cksum)
     call close_file(Topo_restart)

  else if (file_exists(topography_file) .and.                           &
           file_exists(dragtensor_file)) then

!    read and interpolate topography datasets

     if (mpp_pe() == mpp_root_pe()) then
        write ( msg, '("Reading topography file: ",a)')                &
                                                        trim(topography_file)
        call error_mesg('topo_drag_mod', msg, NOTE)
     endif

     if( .not. open_file(topography_fileobj, topography_file, "read")) then
         call error_mesg('topo_drag_mod', "Error opening topography file", FATAL)
     endif

     if( .not. open_file(dragtensor_fileobj, dragtensor_file, "read")) then
         call error_mesg('topo_drag_mod', "Error opening dragtensor file", FATAL)
     endif

     ! check for correct field size in topography
     call get_variable_size(topography_fileobj, trim('hpoz'), siz(1:2))
     if (siz(1) /= ipts .or. siz(2) /= jpts) then
         call error_mesg('topo_drag_mod', 'Field \"hpoz\" in file '//  &
                   trim(topography_file)//' has the wrong size', FATAL)
     endif
     
     allocate (xdatb(ipts+1))
     allocate (ydatb(jpts+1))
     allocate (zdat(ipts,jpts))
     allocate (zout(nlon,nlat))

     do i=1,ipts+1
        xdatb(i) = (i-1)/resolution / Radian
     enddo
     do j=1,jpts+1
        ydatb(j) = (-90.0 + (j-1)/resolution) / Radian
     enddo

     allocate (lonb_loc(nlon+1,nlat+1))
     allocate (latb_loc(nlon+1,nlat+1))
     call atmosphere_grid_bdry (lonb_loc, latb_loc, global=.false.)

     ! initialize horizontal interpolation

     call horiz_interp_init
     call horiz_interp_new ( Interp, xdatb, ydatb, lonb_loc, latb_loc, interp_method="conservative" )
     deallocate(lonb_loc, latb_loc)

     call read_data (topography_fileobj, 'hpoz', zdat)

     exponent = 2. - gamma
     zdat = max(0., zdat)**exponent
     call horiz_interp ( Interp, zdat, zout )

     hmax = (abs(zout)*(gamma + 2.)/(2.*gamma) * &
          (1. - h_frac**(2.*gamma))/(1. - h_frac**(gamma + 2.)))**(1.0/exponent)
     hmin = hmax*h_frac

     if (mpp_pe() == mpp_root_pe()) then
        write ( msg, '("Reading drag tensor file: ",a)')             &
                                                trim(dragtensor_file)
        call error_mesg('topo_drag_mod', msg, NOTE)
     endif

     ! check for correct field size in tensor file

     call get_variable_size(dragtensor_fileobj, tensornames(1), siz(1:2))
     if (siz(1) /= ipts .or. siz(2) /= jpts) then
         call error_mesg('topo_drag_mod', 'Fields in file ' &
         //trim(dragtensor_file)//' have the wrong size', FATAL)
     endif

     do n=1,4
        found_field(n) = variable_exists(dragtensor_fileobj, tensornames(n))
        if (.not. found_field(n)) cycle
        call read_data (dragtensor_fileobj, tensornames(n), zdat)
        call horiz_interp ( Interp, zdat, zout )
        if ( tensornames(n) == 't11' ) then
           t11 = zout/bfscale
        else if ( tensornames(n) == 't21' ) then
           t21 = zout/bfscale
        else if ( tensornames(n) == 't12' ) then
           t12 = zout/bfscale
        else if ( tensornames(n) == 't22' ) then
           t22 = zout/bfscale
        endif
     enddo

     if (.not. found_field(3)) t12 = t21

     deallocate (zdat, zout)
     call horiz_interp_del ( Interp )

     call close_file(topography_fileobj)
     call close_file(dragtensor_fileobj)

  else

     call ERROR_MESG ('topo_drag_init',                                &
                'No sub-grid orography available for topo_drag', FATAL)

  endif

  !SHIELD: perform blocking of topography arrays.
  do nb = 1, Atm_block%nblks
     do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        Sfcprop(nb)%t11(ix) = t11(i,j)
        Sfcprop(nb)%t12(ix) = t12(i,j)
        Sfcprop(nb)%t21(ix) = t21(i,j)
        Sfcprop(nb)%t22(ix) = t22(i,j)
        Sfcprop(nb)%hmin(ix) = hmin(i,j)
        Sfcprop(nb)%hmax(ix) = hmax(i,j)
     enddo
  enddo

  module_is_initialized = .true.

end subroutine topo_drag_init

!=======================================================================

subroutine topo_drag_end

! writes static arrays to restart file

  if (mpp_pe() == mpp_root_pe() ) then
     call error_mesg('topo_drag_mod', 'Writing netCDF formatted restart file: RESTART/topo_drag.res.nc', NOTE)
  endif

  call topo_drag_restart

  module_is_initialized = .false.

end subroutine topo_drag_end

!#######################################################################
! <SUBROUTINE NAME="topo_drag_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine topo_drag_restart(timestamp)
      character(len=*), intent(in), optional :: timestamp
      type(FmsNetcdfDomainFile_t) :: Topo_restart
      character(len=128)  :: restart_fname

      if (present(timestamp)) then
          restart_fname='RESTART/'//trim(timestamp)//'.topo_drag.res.nc'
      else
          restart_fname='RESTART/topo_drag.res.nc'
      endif

      if (.not. open_file(Topo_restart, restart_fname, "overwrite", topo_domain, is_restart = .true.)) then
         call error_mesg("topo_drag_mod", "The topo_drag tiled restart file does not exist", fatal)
      endif

      call topo_drag_register_tile_restart(Topo_restart)
      call write_restart(Topo_restart)
      call add_domain_dimension_data(Topo_restart)
      call close_file(Topo_restart)
end subroutine topo_drag_restart
! </SUBROUTINE>

!#######################################################################

!< Add_dimension_data: Adds dummy data for the domain decomposed axis
subroutine add_domain_dimension_data(fileobj)
  type(FmsNetcdfDomainFile_t) :: fileobj !< Fms2io domain decomposed fileobj
  integer, dimension(:), allocatable :: buffer !< Buffer with axis data
  integer :: is, ie !< Starting and Ending indices for data

    call get_global_io_domain_indices(fileobj, "xaxis_1", is, ie, indices=buffer)
    call write_data(fileobj, "xaxis_1", buffer)
    deallocate(buffer)

    call get_global_io_domain_indices(fileobj, "yaxis_1", is, ie, indices=buffer)
    call write_data(fileobj, "yaxis_1", buffer)
    deallocate(buffer)

end subroutine add_domain_dimension_data

endmodule topo_drag_driver_mod
