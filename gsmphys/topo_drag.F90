module topo_drag_mod

!=======================================================================
! TOPOGRAPHIC DRAG CLOSURE -- Garner (2005)
! Re-written to be compatible with the blocking in SHiELD_physics
!=======================================================================

!-----------------------------------------------------------------------
!  Calculates horizontal velocity tendency due to topographic drag
!-----------------------------------------------------------------------

!use    constants_mod, only: Grav, Cp_Air, Rdgas, Pi
      use physcons, Grav => con_g,  Cp_Air => con_cp,  Pi => con_pi, &
                    Rdgas=> con_rd, cvap  => con_cvap

implicit none

private

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

logical :: module_is_initialized = .false.

! horizontal array size

integer :: kd=0

! parameters:

real, parameter :: u0=1.0       ! arbitrary velocity scale for diagnostics
real, parameter :: xl=80.0e3    ! arbitrary horiz length scale for diagnostics
real, parameter :: ro=1.2       ! arbitrary density scale for diagnostics
real, parameter :: lapse=Grav/Cp_Air ! adiabatic temperature lapse rate
real, parameter :: tiny=1.0e-20

real, parameter :: frint=0.5

! parameters in namelist (topo_drag_nml):

real, public :: &
   frcrit=0.7   &      ! critical value of Froude number for nonlinear flow
  ,alin=1.0     &      ! amplitude of propagating drag
  ,anonlin=5.0  &      ! amplitude of nonpropagating drag
  ,gamma=0.4    &      ! exponent in aspect ratio power law
  ,epsi=0.0     &      ! exponent in distribution power law
  ,beta=0.5     &      ! bluntness of topographic features
  ,h_frac=0.0   &      ! ratio of min to max subgrid mountain height
  ,zref_fac=1.0 &      ! adjusts level separating breaking/laminar flow
  ,tboost=1.0   &      ! surface T boost to improve PBL height estimate
  ,pcut=0.0     &      ! high-level cutoff pressure for momentum forcing
  ,samp=1.0     &      ! correction for coarse sampling of d2v/dz2
  ,max_udt=3.e-3     & ! upper bound on acceleration [m/s2]
  ,no_drag_frac=0.05 & ! fraction of lower atmosphere with no breaking
  ,max_pbl_frac=0.50  ! max fraction of lower atmosphere in PBL
logical, public :: &
   do_conserve_energy=.true. & ! conserve total energy?
  ,keep_residual_flux=.true. & ! redistribute residual pseudomomentum?
  ,do_pbl_average=.false.    & ! average u,rho,N over PBL for baseflux?
  ,use_mg_scaling=.false.    & ! base flux saturates with value 'usat'?
  ,use_mask_for_pbl=.false.  & ! use bottom no_drag_layer as pbl?
  ,use_pbl_from_lock=.false. &  ! use pbl height from Lock boundary scheme  
  ,use_uref_4stable=.false.  

public topo_drag
contains

!#######################################################################

subroutine topo_drag (                                                 &
                                           is, delt, uwnd, vwnd, atmp, &
                                           pfull, phalf, zfull, zhalf, & 
                                                  u_ref, v_ref, z_pbl, & !bqx+ z_pbl
                                       t11, t21, t12, t22, hmin, hmax, & !topo arrays
                              dtaux, dtauy, dtaux_np, dtauy_np, dtemp, &
                              taux, tauy, taus, me )

integer, intent(in) :: is, me
real,    intent(in) :: delt

! INPUT
! -----

! UWND     Zonal wind (dimensioned IX x KDIM)
! VWND     Meridional wind (dimensioned IX x KDIM)
! ATMP     Temperature at full levels (IX x KDIM)
! PFULL    Pressure at full levels (IX x KDIM)
! PHALF    Pressure at half levels (IX x KDIM+1)
! ZFULL    Height at full levels (IX x KDIM)
! ZHALF    Height at half levels (IX x KDIM+1)

real, intent(in), dimension(:,:) :: uwnd, vwnd, atmp
real, intent(in), dimension(:)   :: u_ref, v_ref, z_pbl  !bqx+
real, intent(in), dimension(:,:) :: pfull, phalf, zfull, zhalf

! In SHiELD, pass in the topography arrays as arguments
!   and let the driver handle the memory
real, intent(in), dimension(:) :: t11, t21, t12, t22    ! drag tensor
real, intent(in), dimension(:) :: hmin, hmax

! OUTPUT
! ------

! DTAUX,DTAUY  Tendency of the vector wind in m/s^2 (IX x KDIM)
! DTEMP        Tendency of the temperature in K/s (IX x KDIM)
! TAUX,TAUY    Base momentum flux in kg/m/s^2 (IX) for diagnostics
! TAUS         clipped saturation momentum flux (IX x KDIM) for diagnostics

real, intent(out), dimension(:)   :: taux, tauy
real, intent(out), dimension(:,:) :: dtaux, dtauy, dtaux_np, dtauy_np, dtemp, taus

integer, dimension(size(zfull,1)) :: kpbl, knod, kcut
real,    dimension(size(zhalf,1),size(zhalf,2)) :: tausat

! work arrays

real, dimension(size(zfull,1)) :: taub, taul, taup, taun
real, dimension(size(zfull,1)) :: frulo, fruhi, frunl, rnorm

integer :: idim
integer :: i, k, kdim, km
real    :: dz

  idim = size(uwnd,1)
  kdim = size(uwnd,2)

! estimate height of pbl

  call get_pbl ( atmp, zfull, pfull, phalf, kpbl, knod, kcut )
!
  if (use_pbl_from_lock) then
     kpbl = kdim
     do k = kdim, 2, -1
       where ( zfull(:,k) < zhalf(:,kdim+1) + z_pbl(:) )
           kpbl(:) = k - 1           ! the first full model level above PBL
       endwhere
    enddo
  endif

! calculate base flux

!!$!! Check if inputs are right
!!$  do i=1,idim
!!$        write(100+me,*) "First check:"
!!$        write(100+me,*) i, u_ref(i), v_ref(i), z_pbl(i)
!!$        write(100+me,*) t11(i), t21(i), t12(i), t22(i)
!!$        write(100+me,*) hmin(i), hmax(i)
!!$        do k=KDIM-2, KDIM
!!$           write(100+me,'(I, 5G)') k, zfull(i,k), pfull(i,k), uwnd(i,k), vwnd(i,k), atmp(i,k)
!!$        enddo
!!$        !write(100+me,*) KDIM+1, zhalf(i,KDIM+1), phalf(i,KDIM+1)
!!$        write(100+me,*) ' '
!!$     exit !loop
!!$  enddo

  call base_flux (                                                     &
                                                 is, uwnd, vwnd, atmp, &
                                                  u_ref, v_ref, z_pbl, &  !bqx+
                                                   t11, t12, t21, t22, &
                                                           hmin, hmax, &                                             
                                             taux, tauy, dtaux, dtauy, &
                                               taub, taul, taup, taun, &
                                           frulo, fruhi, frunl, rnorm, &
                                     zfull, zhalf, pfull, phalf, kpbl, me )


!!$!! DEBUG CODE
!!$  do i=1,idim
!!$     !if (hmin(i) > 50.) then
!!$        write(100+me,*) "Second check:"
!!$        write(100+me,*) i, u_ref(i), v_ref(i), z_pbl(i)
!!$        write(100+me,*) hmin(i), hmax(i)
!!$        write(100+me,*) taux(i), tauy(i)
!!$        write(100+me,*) taub(i), frulo(i), fruhi(i), frunl(:)
!!$        write(100+me,*) ' '
!!$     !endif
!!$     exit !loop
!!$  enddo

! calculate saturation flux profile

  call satur_flux (                                                    &
                                                     uwnd, vwnd, atmp, &
                                                   taup, taub, tausat, &
                                                  frulo, fruhi, frunl, &
                        dtaux, dtauy, zfull, pfull, phalf, kpbl, kcut )

!!$  do i=1,idim
!!$     !if (hmin(i) > 50.) then
!!$        write(100+me,*) "Third check:"
!!$        write(100+me,*) i, u_ref(i), v_ref(i), z_pbl(i)
!!$        write(100+me,*) hmin(i), hmax(i)
!!$        write(100+me,*) taux(i), tauy(i)
!!$        write(100+me,*) ' '
!!$     !endif
!!$     exit !loop
!!$  enddo

! calculate momentum tendency

  call topo_drag_tend (                                                &
                                               delt, uwnd, vwnd, atmp, &
                                       taux, tauy, taul, taun, tausat, &
  dtaux, dtauy, dtaux_np, dtauy_np, dtemp, zfull, zhalf, pfull, phalf, kpbl ) !bqx+ dtaux_np, dtauy_np

! put saturation flux profile into 'taus' for diagnostics

  do k=1,kdim
!     taus(:,k) = 0.5*rnorm(:)*(tausat(:,k) + tausat(:,k+1))
     taus(:,k) = tausat(:,k+1)*taub(:)/taul(:)
  enddo

! put total drag into 'taux,tauy' for diagnostics

  taup = taup - tausat(:,1)
  taub = (taup + taun)/taul
  taux = taux*taub
  tauy = tauy*taub


!!! Choose a debug cell
!!! DEBUG CODE
 
!!$  do i=1,idim
!!$!     if (hmin(i) > 50.) then
!!$     write(100+me,*) "Final check"
!!$        write(100+me,*) i, u_ref(i), v_ref(i), z_pbl(i)
!!$        write(100+me,*) 
!!$        do k=1, KDIM
!!$           write(100+me,'(I, 9G)') k, dtaux(i,k), dtauy(i,k), tausat(i,k)
!!$        enddo
!!$        write(100+me,*) KDIM+1, taux(i), tauy(i), taul(i), taun(i)
!!$        write(100+me,*) ' '
!!$!     endif
!!$     exit !loop
!!$  enddo

!!! END DEBUG CODE

end subroutine topo_drag

!=======================================================================
                                  
subroutine base_flux (                                                 &
                                                 is, uwnd, vwnd, atmp, &
                                                  u_ref, v_ref, z_pbl, & !bqx
                                                   t11, t12, t21, t22, &
                                                           hmin, hmax, &
                                             taux, tauy, dtaux, dtauy, &
                                               taub, taul, taup, taun, &
                                           frulo, fruhi, frunl, rnorm, &
                                     zfull, zhalf, pfull, phalf, kpbl, me )

integer, intent(in) :: is, me
real, intent(in),  dimension(:,:) :: uwnd, vwnd, atmp
real, intent(in),  dimension(:)   :: u_ref, v_ref, z_pbl
real, intent(in),  dimension(:,:) :: zfull, zhalf, pfull, phalf
real, intent(in),  dimension(:)   :: t11, t12, t21, t22
real, intent(in),  dimension(:)   :: hmin, hmax
real, intent(out), dimension(:)   :: taux, tauy
real, intent(out), dimension(:,:) :: dtaux, dtauy
real, intent(out), dimension(:)   :: taub, taul, taup, taun
real, intent(out), dimension(:)   :: frulo, fruhi, frunl, rnorm
integer, intent(in), dimension(:) :: kpbl

real, dimension(size(uwnd,1)) :: ubar, vbar

integer :: i, idim, id
integer :: k, kdim, kb, kbp, kt, km

real :: usat, bfreq2, bfreq, dphdz, vtau, d2udz2, d2vdz2
real :: dzfull, dzhalf, dzhalf1, dzhalf2, density
real :: frmin, frmax, frmed, frumin, frumax, frumed, fruclp, fruclm
real :: rnormal, gterm, hterm, fru0, frusat
real :: usum, vsum, n2sum, delp

logical :: do_prt = .true.

  idim = size(uwnd,1)
  kdim = size(uwnd,2)

! compute base flux

     do i=1,idim
        usum = 0.
        vsum = 0.
        kt = kpbl(i)
        kb = max(kd,kt)
        do k=kt,kb
           delp = phalf(i,k+1) - phalf(i,k)
           usum = usum + uwnd(i,k)*delp
           vsum = vsum + vwnd(i,k)*delp
        enddo
        ubar(i) = usum/(phalf(i,kb+1) - phalf(i,kt))
        vbar(i) = vsum/(phalf(i,kb+1) - phalf(i,kt))
     enddo

  if (use_uref_4stable) then
   where( z_pbl (:) == 0. )
    ubar(:) = u_ref(:)
    vbar(:) = v_ref(:)
   endwhere
  endif

     do i=1,idim
        id = is+i-1
        kt = kpbl(i)
        kb = max(kd,kt)
        kbp = min(kdim, kb+1) !bqx+
        dzfull = zhalf(i,kt) - zhalf(i,kb+1)
        density = (phalf(i,kb+1) - phalf(i,kt))/(Grav*dzfull)
        dzfull = zfull(i,kt-1) - zfull(i,kbp)
        bfreq2 = Grav*((atmp(i,kt-1) - atmp(i,kbp))/dzfull+lapse)/&
                  (0.5*(atmp(i,kt-1) + atmp(i,kbp)))
!
        bfreq = sqrt(max(tiny, bfreq2))

!       included 'alin' 4/2015

        taux(i) = (ubar(i)*t11(id) + vbar(i)*t21(id))      &
                                                   *bfreq*density
        tauy(i) = (ubar(i)*t12(id) + vbar(i)*t22(id))      &
                                                   *bfreq*density

        taub(i) = max(tiny, sqrt(taux(i)**2 + tauy(i)**2))

!!$!!! DEBUG CODE 
!!$        if (do_prt) then
!!$           write(100+me,*) 'base_flux: ', i, id, taux(i), tauy(i)
!!$           write(100+me,*) i, id, u_ref(i), v_ref(i), z_pbl(i)
!!$           write(100+me,*) t11(i), t21(i), t12(i), t22(i)
!!$           write(100+me,*) hmin(i), hmax(i)
!!$           write(100+me,*) ubar(i), t11(i), vbar(i), t21(i)
!!$           write(100+me,*) density, dzfull, phalf(i,kb+1), phalf(i,kt)
!!$           do_prt = .false.
!!$        endif
!!$           
!!$!!! END DEBUG CODE

!       min/max Froude numbers based on low-level flow

        vtau = max(tiny, -(ubar(i)*taux(i)                         &
                         + vbar(i)*tauy(i))/taub(i))
        frmax = hmax(id)*bfreq / vtau
        frmin = hmin(id)*bfreq / vtau
        frmed = frcrit + frint

!       linear momentum flux associated with min/max Froude numbers

        dphdz = bfreq / vtau
        usat = sqrt(density/ro) * vtau / sqrt(dphdz*xl)
        frusat = frcrit*usat

        frumin = frmin*usat
        frumax = frmax*usat
        frumed = frmed*usat

        frumax = max(frumax,frumin + tiny)
        fruclp = min(frumax,max(frumin,frusat))
        fruclm = min(frumax,max(frumin,frumed))
        fru0 = (u0/vtau)*usat

!       total drag in linear limit

        rnormal =                                                      &
                   (frumax**(2.0*gamma - epsi)                         &
                  - frumin**(2.0*gamma - epsi))/(2.0*gamma - epsi)
        rnormal = fru0**gamma * ro/rnormal  

        taul(i) =                                                    &
                   (frumax**(2.0 + gamma - epsi)                       &
                  - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi)

!       separate propagating and nonpropagating parts of total drag

        gterm = frusat**(beta + 1.0)*                                  &
                (frumax**(gamma - epsi - beta)                         &
               - fruclp**(gamma - epsi - beta))/(gamma - epsi - beta)

        taup(i) =  alin *                                             &
                  ( (fruclp**(2.0 + gamma - epsi)                       &
                   - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi) &
                                                       + frusat*gterm )

        taun(i) = anonlin*usat/(1.0 + beta) *                        &
                 ( (frumax**(1.0 + gamma - epsi)                       &
                  - fruclp**(1.0 + gamma - epsi))/(1.0 + gamma - epsi) &
                                                      - gterm )

!       5/2015 mg option: depth of blocking ~ U/N, not h

        if (use_mg_scaling) taun(i) = taun(i)/max(frmax,frcrit)

        fruhi(i) = frumax
        frulo(i) = frumin
        frunl(i) = frusat
        rnorm(i) = rnormal

     enddo
! wind component opposite the drag at full levels (stored as 'dtaux')

  do k=1,kdim
     do i=1,idim
        dtaux(i,k) =                                              &
             -(uwnd(i,k)*taux(i) + vwnd(i,k)*tauy(i))/taub(i)
     enddo
  enddo

! curvature of wind at full levels (stored as 'dtauy')

  dtauy = 0.

     do i=1,idim
!        kt = kpbl(i)
        kt = min(kdim-1, kpbl(i)) !bqx+
        do k=2,kt
           dzfull = zhalf(i,k) - zhalf(i,k+1)
           dzhalf1 = zfull(i,k-1) - zfull(i,k)
           dzhalf2 = zfull(i,k) - zfull(i,k+1)
           d2udz2 = ((uwnd(i,k-1) - uwnd(i,k  ))/dzhalf1           &
                   - (uwnd(i,k  ) - uwnd(i,k+1))/dzhalf2)/dzfull
           d2vdz2 = ((vwnd(i,k-1) - vwnd(i,k  ))/dzhalf1           &
                   - (vwnd(i,k  ) - vwnd(i,k+1))/dzhalf2)/dzfull
           dtauy(i,k) = -(d2udz2*taux(i) + d2vdz2*tauy(i))/      &
                                                              taub(i)
        enddo

     enddo

end subroutine base_flux

!=======================================================================

subroutine satur_flux (                                                &
                                                     uwnd, vwnd, atmp, &
                                                   taup, taub, tausat, &
                                                  frulo, fruhi, frunl, &
                        dtaux, dtauy, zfull, pfull, phalf, kpbl, kcut )

real, intent(in),  dimension (:,:) :: uwnd, vwnd, atmp
real, intent(in),  dimension (:,:) :: dtaux, dtauy
real, intent(in),  dimension (:,:) :: zfull, pfull, phalf
real, intent(in),  dimension (:)   :: taup
real, intent(out), dimension (:)   :: taub
real, intent(out), dimension (:,:) :: tausat
real, intent(in),  dimension (:)   :: frulo, fruhi, frunl
integer, intent(in), dimension (:) :: kpbl, kcut

real, dimension(size(zfull,1)) :: usat

real :: dzhalf, gterm, gterm0, density
real :: bfreq2, bfreq, vtau, d2vtau, dphdz, xl1
real :: frumin, frumax, fruclp, frusat, frusat0, fruclp0

integer :: i, idim
integer :: k, kdim, k1

  idim = size(uwnd,1)
  kdim = size(uwnd,2)

! get vertical profile of propagating part of momentum flux

  usat = frunl/frcrit

  do k=kdim,2,-1
        do i=1,idim

!          buoyancy frequency, velocity and density at half levels

           dzhalf = zfull(i,k-1) - zfull(i,k)
           density = (pfull(i,k) - pfull(i,k-1))/(Grav*dzhalf)
           bfreq2 = Grav*                                              &
                       ((atmp(i,k-1) - atmp(i,k))/dzhalf + lapse)/ & 
                   (0.5*(atmp(i,k-1) + atmp(i,k)))
           bfreq = sqrt(max(tiny, bfreq2))
           
           vtau = max(tiny, 0.5*(dtaux(i,k-1) + dtaux(i,k)))

!          WKB correction of vertical wavelength

           d2vtau = 0.5*(dtauy(i,k-1) + dtauy(i,k))
           xl1 = xl*max(0.5, min(2.0, 1.0 - samp*vtau*d2vtau/(bfreq*bfreq)))

!          min/max and critical momentum flux values at half levels

           dphdz = bfreq / vtau
           usat(i) = min(usat(i),sqrt(density/ro) * vtau/sqrt(dphdz*xl1))
           frusat = frcrit*usat(i)

           frumin = frulo(i)
           frumax = fruhi(i)
           fruclp = min(frumax,max(frumin,frusat))
           frusat0 = frunl(i)
           fruclp0 = min(frumax,max(frumin,frusat0))

!          propagating part of momentum flux (from WKB or EP)

           gterm0 = (frumax**(gamma - epsi - beta)                     &
                - fruclp0**(gamma - epsi - beta))/(gamma - epsi - beta)
           gterm = (fruclp0**(gamma - epsi)                            &
                               - fruclp**(gamma - epsi))/(gamma - epsi)

           tausat(i,k) = alin *                                      &
                 ( (fruclp**(2.0 + gamma - epsi)                       &
                  - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi) &
                         + frusat**2.0*(gterm0*frusat0**beta + gterm) )
        enddo
  enddo

! make propagating flux constant with height in zero-drag top layer
! changed 5/2014

  k1 = maxval(kcut)
  do k=k1,1,-1
     where (k <= kcut)
        tausat(:,k) = tausat(:,k+1)
     endwhere
  enddo

! make propagating flux constant with height in zero-drag surface layer

  k1 = minval(kpbl)
  do k=kdim+1,k1+1,-1
     where (k > kpbl)
        tausat(:,k) = taup
     endwhere
  enddo

! redistribute residual forcing

  if ( keep_residual_flux ) then
     taub(:) = tausat(:,1)/(phalf(:,kdim+1) - phalf(:,1))
     do k=1,kdim
        tausat(:,k) = tausat(:,k)                                  &
                        - taub(:)*(phalf(:,kdim+1) - phalf(:,k))
     enddo
  endif

endsubroutine satur_flux

!=======================================================================

subroutine topo_drag_tend (                                            &
                                               delt, uwnd, vwnd, atmp, &
                                       taux, tauy, taul, taun, tausat, &
                dtaux, dtauy, dtaux_np, dtauy_np, dtemp, zfull, zhalf, &
                                                    pfull, phalf, kpbl )

real, intent(in) :: delt
real, intent(in), dimension(:,:)    :: uwnd, vwnd, atmp
real, intent(in), dimension(:,:)    :: zfull, zhalf, pfull, phalf
real, intent(in), dimension(:)      :: taux, tauy, taul, taun
real, intent(inout), dimension(:,:) :: tausat
real, intent(inout), dimension(:,:) :: dtaux, dtauy, dtemp
real, intent(out), dimension(:,:)   :: dtaux_np, dtauy_np !bqx
integer, intent(in), dimension (:)  :: kpbl

real, parameter :: bfmin=0.7e-2, bfmax=1.7e-2  ! min/max buoyancy freq [1/s]
real, parameter :: vvmin=1.0                   ! minimum surface wind [m/s]

integer,dimension(size(zfull,1)) :: kref
real :: dzhalf, zlast, rscale, phase, bfreq, bfreq2, vtau
real :: gfac, gfac1, dp, weight, wtsum, taunon, taunon1

integer :: i, idim
integer :: k, kdim, kr, kt

real,dimension(size(zfull,1)) :: dx, dy

  idim = size(uwnd,1)
  kdim = size(uwnd,2)

! find reference level for non-propagating drag (z ~ pi U/N)

!the following re-orients the drag to align with low-level wind
!do i=1,idim
!k = knod(i)
!gfac = sqrt( (taux(i)**2 + tauy(i)**2) / max (tiny, uwnd(i,k)**2 + vwnd(i,k)**2) )
!dx(i) = gfac * uwnd(i,k)
!dy(i) = gfac * vwnd(i,k)
!enddo

     do i=1,idim
        k = kpbl(i)
!stg        k = kdim
        phase = 0.0
        zlast = zhalf(i,k)
        do while (phase <= Pi*zref_fac .and. k > 1)
           k = k-1
           vtau = 0.5*(dtaux(i,k-1) + dtaux(i,k))
           dzhalf = zfull(i,k-1) - zfull(i,k)
           bfreq2 = Grav*                                              &
                       ((atmp(i,k-1) - atmp(i,k))/dzhalf + lapse)/ &
                   (0.5*(atmp(i,k-1) + atmp(i,k)))
           bfreq = sqrt(max(tiny, bfreq2))
           rscale = max(bfmin, min(bfmax, bfreq))/max(vvmin, vtau)
           dzhalf = zfull(i,k-1) - zlast
           phase = phase + dzhalf*rscale
           zlast = zfull(i,k-1)
        enddo
        kref(i) = k
     enddo

! CALCULATE DECELERATION DUE TO PROPAGATING DRAG (~-rho^-1 dtau/dz)

  do k=1,kdim
        do i=1,idim
          dp = phalf(i,k+1) - phalf(i,k)
          gfac = tausat(i,k+1) - tausat(i,k)
          gfac1 = gfac*Grav/(dp*taul(i))
          dtaux(i,k) = gfac1*taux(i)
          dtauy(i,k) = gfac1*tauy(i)
!dtaux(i,k) = -gfac1*dx(i)
!dtauy(i,k) = -gfac1*dy(i)
        enddo
  enddo

! CALCULATE DECELERATION DUE TO NON-PROPAGATING DRAG
  dtaux_np = 0. !bqx
  dtauy_np = 0. !bqx

     do i=1,idim
        kr = kref(i)
        kt = kpbl(i)
!stg        kt = kdim
        gfac = taun(i)/taul(i) * Grav
        wtsum = 0.0
        do k=kr,kt
           dp = phalf(i,k+1) - phalf(i,k)
           weight = pfull(i,k) - phalf(i,kr)
           wtsum = wtsum + dp*weight
        enddo
        taunon = 0.0
        do k=kr,kt
           weight = pfull(i,k) - phalf(i,kr)
           gfac1 = gfac*weight/wtsum
           dtaux(i,k) = dtaux(i,k) + gfac1*taux(i)
           dtauy(i,k) = dtauy(i,k) + gfac1*tauy(i)
!bqx
           dtaux_np(i,k) = gfac1*taux(i)
           dtauy_np(i,k) = gfac1*tauy(i)
!dtaux(i,k) = dtaux(i,k) - gfac1*dx(i)
!dtauy(i,k) = dtauy(i,k) - gfac1*dy(i)
           dp = phalf(i,k+1) - phalf(i,k)
           taunon = taunon + gfac1*dp
           taunon1 = taunon*taul(i)/Grav
           tausat(i,k) = tausat(i,k) + taunon1
        enddo
        do k=kt+1,kdim
           tausat(i,k) = tausat(i,k) + taunon1
        enddo
     enddo

  dtaux = max(-max_udt, min(max_udt, dtaux))   !stg
  dtauy = max(-max_udt, min(max_udt, dtauy))

! CALCULATE HEATING TO CONSERVE TOTAL ENERGY

  if (do_conserve_energy) then
     dtemp = -((uwnd + 0.5*delt*dtaux)*dtaux                           &
             + (vwnd + 0.5*delt*dtauy)*dtauy)/Cp_Air
  else
     dtemp = 0.0
  endif

end subroutine topo_drag_tend

!=======================================================================

subroutine get_pbl ( atmp, zfull, pfull, phalf, kpbl, knod, kcut )

integer, intent(out), dimension(:) :: kpbl, knod, kcut
real, intent(in), dimension(:,:)   :: atmp
real, intent(in), dimension(:,:)   :: zfull, pfull, phalf

real, dimension(size(pfull,1)) :: ppbl, pbot
real, dimension(size(pfull,1)) :: tbot, zbot

integer :: i, idim
integer :: k, kdim

  idim = size(atmp,1)
  kdim = size(atmp,2)

     do i=1,idim
        ppbl(i) = (1.0 - max_pbl_frac)*phalf(i,kdim+1)
        pbot(i) = (1.0 - no_drag_frac)*phalf(i,kdim+1)
        tbot(i) = atmp(i,kdim) + tboost
        zbot(i) = zfull(i,kdim)
     enddo

! find highest model level in no-drag surface layer

  knod = kdim-1

  do k=kdim-2,2,-1
     where ( pfull(:,k) >= pbot(:) )
        knod = k
     endwhere
  enddo

! find lowest model level in no-drag top layer

  kcut = 1

  do k=2,kdim
     where ( pfull(:,k) <= pcut )
        kcut = k
     endwhere
  enddo

  if (kd == 0 .and. do_pbl_average) kd = kdim-1

  if ( use_mask_for_pbl ) then
     kpbl = knod
     return
  endif

! find the first layer above PBL

  kpbl = kdim-1

  do k=kdim-2,2,-1
     where ( pfull(:,k) >= ppbl(:) .and.                          &
        tbot(:) - atmp(:,k) > lapse*(zfull(:,k) - zbot(:)) )
        kpbl = k -1
     endwhere
  enddo

end subroutine get_pbl

endmodule topo_drag_mod
