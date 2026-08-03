! SNAP: Servere Nuclear Accident Programme
! Copyright (C) 1992-2017   Norwegian Meteorological Institute

! This file is part of SNAP. SNAP is free software: you can
! redistribute it and/or modify it under the terms of the
! GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module rwalkML
  USE iso_fortran_env, only: real64, int32

  implicit none
  private

  real(real64), save :: vrdbla ! l-eta above mixing height
  real(real64), save :: tfactor_v ! tfactor=tstep/tmix
  real(real64), save :: tsqrtfactor_v ! tsqrtfactor_v=1/sqrt(tmix/tstep)
  real(real64), save :: tfactor_h ! tfactor=tstep/thour
  real(real64), save :: tsqrtfactor_h ! tsqrtfactor_h=1/sqrt(thour/tstep)
  real(real64), save :: tstep

  real(real64), parameter :: hmax = 2500.0 ! maximum mixing height
  real(real64), parameter :: tmix_v = 15.0*60.0 ! Characteristic mixing time = 15 min (to reach full bl-height)
  real(real64), parameter :: tmix_h = 15.0*60.0 ! Horizontal base-time time = 15 min (to reach ax^b width)
  real(real64), parameter :: lmax = 0.28 ! Maximum l-eta in the mixing layer
  real(real64), parameter :: labove = 0.03 ! Standard l-eta above the mixing layer
  real(real64), parameter :: entrainment = 0.10 ! Entrainment zone = 10%*h
  real, parameter :: dt_min = 1.0 ! minimum timestep size for adaptive timestepping

  ! Values for random number generation
  integer(int32), parameter :: max_rands=10000000
  integer(int32) :: nrand=1
  real :: rands(max_rands) ! initialise array to store random numbers

  real(real64), save, public :: a_in_bl = 0.5
  real(real64), save, public :: a_above_bl = 0.25
  real(real64), save, public :: b = 0.875
  logical, save, public :: turb_homogeneous = .FALSE.
  logical, save, public :: well_mixed_test = .FALSE.
  logical, save, public :: blfullmix = .FALSE.

  ! Turbulence flags
  logical, save, public :: diffusion_in_metres             = .FALSE.
  logical, save, public :: turbulence_fields_required      = .FALSE.
  logical, save, public :: density_correction              = .FALSE.
  logical, save, public :: scheme_is_tke                   = .FALSE.
  logical, save, public :: scheme_is_random_walk           = .FALSE.
  logical, save, public :: scheme_uses_adaptive_above_bl   = .FALSE.
  integer, save, public :: diffusion_scheme_id = 0
  integer, save, public :: bl_id = 0

  character(len=64), save, public :: diffusion_scheme = ''
  character(len=64), save, public :: bl_definition = ''
  character(len=64), save, public :: meteo_type = ''
  character(len=64), save, public :: entrainment_scheme = ''

  public rwalk_init, diffusion_fields, air_density, turbulence_master, eta_to_metres, &
         metres_to_eta

  contains

!> Initialise constants needed for rwalk
subroutine rwalk_init(timestep)
  use init_random_seedML, only: generate_normal_randoms
!> time step in seconds (trajectory calculations)
  real, intent(in) :: timestep

  tfactor_v = timestep/tmix_v
  tsqrtfactor_v=sqrt(tfactor_v)
  tfactor_h = timestep/tmix_h
  tsqrtfactor_h=sqrt(tfactor_h)
  tstep = timestep

  ! l-eta above mixing height
  vrdbla = labove*tsqrtfactor_v

  ! If diffusion scheme not default SNAP, then generate normally distributed random nums
  if (diffusion_scheme_id /= 0) then
    call generate_normal_randoms(rands, max_rands)
  endif

end subroutine

subroutine turbulence_master(part,pextra, dt_remaining, adaptive)
  USE particleML, only: extraParticle, Particle

    !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle (u, v, rmx, rmy)
  type(extraParticle), intent(inout) :: pextra
  !> full mixing in boundarylayer (true=old,false=new)

  ! Variables for determining adaptive time step size
  real, intent(in) :: dt_remaining
  logical, intent(in) :: adaptive

  if (diffusion_scheme_id == 2) then
    ! Check if particle within abl
    if (part%zmetres.lt.part%hbl) then
      call flexpart_diffusion_within_abl(part,pextra, dt_remaining, adaptive)
    else
      call flexpart_diffusion_above_abl(part, pextra, dt_remaining)
    endif
  elseif (diffusion_scheme_id == 5) then
    ! Check if particle within abl
    if (part%zmetres.lt.part%hbl) then
      call variable_k_inhomog(part, pextra)
    else
      call variable_k_inhomog_above(part, pextra)
    endif
  elseif (diffusion_scheme_id == 3) then
    ! BL vs FT handled by minimum values in flexpart function
    if (adaptive) then
      call random_walk_name(part, pextra, dt_remaining, adaptive)
    else
      if (part%zmetres.lt.part%hbl) then
        call random_walk_name_fixed_below2(part, pextra, dt_remaining, adaptive)
      else
        call random_walk_name_fixed_above(part, pextra, dt_remaining, adaptive)
      endif
    endif
  elseif (diffusion_scheme_id == 4) then
    ! Check if particle within abl
    if (part%zmetres.lt.part%hbl) then
      call fixed_k_within_bl(part)
    else
      call fixed_k_above_bl(part)
    endif
  elseif (diffusion_scheme_id == 6) then
    call tke_diffusion(part, pextra, dt_remaining, adaptive)
  elseif (diffusion_scheme_id == 7) then
    ! Hybrid near-field far-field scheme
    if (adaptive) then
      call random_walk_name(part, pextra, dt_remaining, adaptive)
    else
      if (part%zmetres.lt.part%hbl) then
        call variable_k_within_bl(part, pextra)
      else
        call variable_k_above_bl(part, pextra)
      endif
    endif
  else 
    call rwalk(blfullmix, part, pextra)
  endif

end subroutine

!>  Purpose:  Diffusion, in and above boudary layer.
!>
!>  Method:   Random walk.
!>
!> ::rwalk_init must be run before rwalk
subroutine rwalk(blfullmix,part,pextra)
!   24.04.2009 Jerzy Bartnicki: Model particle which goes below the
!   ground or above the top boundary in the random walk is reflected
!   26.03.2011 Jerzy Bartnicki: New parameterization of vertical diffusion in the
!   mixing layer. l-eta proportional to mixing height and the time step.
!   For mixing height = 2500 m and time step = 15 min:
!   In ABL: l-eta=0.28
!   Above ABL: l-eta=0.003
!   For 200< mixing height<2500 and arbitrary time step:
!   In ABL: l-eta=0.28*(mh/2500m)*(tstep/tstep-mix)
!   Above ABL: l-eta=0.003*(tstep/tstep-mix)
!   Entrainment zone = 10%*h
  USE particleML, only: extraParticle, Particle
!> full mixing in boundarylayer (true=old,false=new)
  logical, intent(in) :: blfullmix
!> particle with information
  type(Particle), intent(inout)  :: part
!> extra information regarding the particle (u, v, rmx, rmy)
  type(extraParticle), intent(in) :: pextra

  real(real64) :: rnd(3), rl, vabs
  real(real64) :: rv, top_entrainment, bl_entrainment_thickness

  real(real64) :: a

! the random_number function returns 3 (x,y,z) random real numbers between 0.0 and 1.0
  call random_number(rnd)
  rnd = rnd - 0.5

! horizontal diffusion
  if (part%z > part%tbl) then ! in boundary layer
    a = a_in_bl
  else ! above boundary layer
    a = a_above_bl
  endif

  if (.NOT.well_mixed_test) then
    vabs = hypot(pextra%u, pextra%v)
    rl = 2*a*((vabs*tmix_h)**b) * tsqrtfactor_h ! sqrt error/sigma propagation
    part%x = part%x + rl*rnd(1)*pextra%rmx
    part%y = part%y + rl*rnd(2)*pextra%rmy
  endif


! vertical diffusion
  if (part%z <= part%tbl) then ! Above boundary layer

    ! if ((part%z + vrdbla*rnd(3)).gt.part%tbl) then ! reflect off BL top
    !   part%z = 2 * part%tbl - (part%z + vrdbla*rnd(3))
    ! else
    !   part%z = part%z + vrdbla*rnd(3)
    ! endif
    part%z = part%z + vrdbla*rnd(3)

  else ! In boundary layer
    bl_entrainment_thickness = (1.0 - part%tbl)*(1.+entrainment)
    if (blfullmix .or. (tsqrtfactor_v .gt. 1.0)) then ! full mixing
      part%z = 1.0 - bl_entrainment_thickness*(rnd(3)+0.5)
    else ! vertical mixing split in smaller time-steps   
      rv  = (1-part%tbl)*tsqrtfactor_v

      part%z = part%z + rv*rnd(3)

      !... reflection from the ABL top
      !... but allow for entrainment
      ! top_entrainment 10% higher than tbl
      top_entrainment = max(0., 1.0 - bl_entrainment_thickness)
      if(part%z < top_entrainment) then
        part%z = 2.0*part%tbl - part%z 
      endif

    !... reflection from the bottom
      if(part%z > 1.0) then
        part%z = 2.0 - part%z
      endif

    !..vertical limits
      part%z = min(part%z, 1.0d0)
      part%z = max(part%z, real(top_entrainment, kind=kind(part%z)))
    end if
  end if
end subroutine rwalk

subroutine flexpart_diffusion_within_abl(part, pextra, dt_remaining, adaptive)
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  real, intent(in) :: dt_remaining
  logical, intent(in) :: adaptive

  real :: sigu, sigv, sigw ! Turbulent velocity standard deviations
  real :: tlu, tlv, tlw ! Lagrangian timescales
  real :: ru, rv, rw ! Turbulence correlation terms
  real :: delz ! Turbulent vertical displacement (m)
  real :: dsigwdz
  real :: dttlw
  real :: density_corr ! Density correction term
  real :: scaled_height
  real :: wst ! convective scale velocity
  real :: ol ! obukhov length
  real :: ust ! friction velocity

  ! Dimensionless height 
  if (turb_homogeneous) then ! Take middle BL value if homogeneous
    scaled_height = 0.5
  else
    scaled_height = max(part%zmetres/part%hbl, 1e-3)
  endif

  ust = pextra%ust
  wst = pextra%wst
  ol = pextra%ol

  if (well_mixed_test) then
    wst = 1.5
    ol = -100
    ust = 0.5
  endif

  ! Case 1, Neutral Conditions
  if (part%hbl/ABS(ol).lt.1.) then
    ust = max(1.e-4, ust)

    ! Eq. 7.25 Hanna 1982: 
    sigu = 2.0 * ust * EXP(-3.e-4*part%zmetres/ust)
    sigu = MAX(sigu, 1.e-5)

    ! Eq. 7.26 Hanna 1982:
    sigv = 1.3 * ust * EXP(-2.e-4*part%zmetres/ust)
    sigv=max(sigv,1.e-5)
    sigw=sigv

    ! Vertical gradient of sigw
    dsigwdz=(-2.e-4*sigw)/ust

    ! Lagrangian timescales
    tlu=0.5*part%zmetres/sigw/(1.+1.5e-3*part%zmetres/ust)
    tlv=tlu
    tlw=tlu

  ! Case 2 , Unstable Conditions
  elseif (ol.lt.0.) then

    sigu = ust * (12 - 0.5 * part%hbl/ol)**0.33333
    sigv = sigu

    sigw = (1.2 * wst**2 *(1-0.9*scaled_height)*(scaled_height)**0.6666 + (1.8 - 1.4*scaled_height)*ust**2)**0.5
    dsigwdz=0.5/sigw/part%hbl*(-1.4*ust**2+wst**2*(0.8*max(scaled_height,1.e-3)**(-.33333)-1.8*scaled_height**0.66666))

    ! Lagrangian timescales
    tlu = 0.15 * part%hbl/sigu
    tlv = tlu
    if (part%zmetres.lt.abs(ol)) then
      tlw=0.1*part%zmetres/(sigw*(0.55-0.38*abs(part%zmetres/ol)))
    else if (scaled_height.lt.0.1) then
      tlw=0.59*part%zmetres/sigw
    else
      tlw=0.15*part%hbl/sigw*(1.-exp(-5*scaled_height))
    endif
    
  ! Case 3, Stable Conditions 
  else
    sigu=2.*ust*(1.-scaled_height) !. 7.20, Hanna
    sigv=1.3*ust*(1.-scaled_height) !. 7.19, Hanna
    sigu=max(sigu,1.e-6)
    sigv=max(sigv,1.e-6)
    sigw=sigv !. 7.19, Hanna
    dsigwdz=-1.3*ust/part%hbl

    ! Lagrangian timescales
    tlu=0.15*part%hbl/sigu*(sqrt(scaled_height)) !. 7.22, Hanna
    tlv=0.07*part%hbl/sigv*(sqrt(scaled_height)) !. 7.23, Hanna
    tlw=0.1*part%hbl/sigw*scaled_height**0.8

  endif

  ! Clamp lagrangian timescales
  tlu=max(10.,tlu)
  tlv=max(10.,tlv)
  tlw=max(30.,tlw)

  part%tlw = tlw
  part%dsigwdz = dsigwdz
  part%sigw = sigw

  ! Based on turbulence parameters, calculate size of integration time step
  call set_turbulence_timestep(part, dt_remaining, adaptive)

  ! Calculate turbulent horizontal velocities
  if (nrand+1.gt.max_rands) nrand=1
  if (part%ptstep/tlu.lt..5) then
    part%turbvelu=(1.-part%ptstep/tlu)*part%turbvelu+rands(nrand)*sigu*sqrt(2.*part%ptstep/tlu)
  else
    ru=exp(-part%ptstep/tlu)
    part%turbvelu=ru*part%turbvelu+rands(nrand)*sigu*sqrt(1.-ru**2)
  endif
  if (part%ptstep/tlv.lt..5) then
    part%turbvelv=(1.-part%ptstep/tlv)*part%turbvelv+rands(nrand+1)*sigv*sqrt(2.*part%ptstep/tlv)
  else
    rv=exp(-part%ptstep/tlv)
    part%turbvelv=rv*part%turbvelv+rands(nrand+1)*sigv*sqrt(1.-rv**2)
  endif
  nrand=nrand+2

  ! Transform turbulent velocities from along/crosswind to x/y
  call align_turbvels(part, pextra)

  density_corr=pextra%rhograd/pextra%rho

  ! ratio of time step to lagrangian timescale for autocorrelation
  dttlw = part%ptstep/tlw

  if (nrand+1.gt.max_rands) nrand=1
  ! Calculate turbulent vertical velocity
  if (dttlw.lt..5) then ! Small adaptive time steps
    part%turbvelw=((1.-dttlw)*part%turbvelw + part%ptstep*(dsigwdz+density_corr*sigw)) * part%icbt &
       + sqrt(2.*dttlw) * rands(nrand)
    delz=part%turbvelw*sigw*part%ptstep
  else ! larger time steps
    rw=exp(-dttlw)
    part%turbvelw=(rw*part%turbvelw +tlw*(1.-rw)*(dsigwdz+density_corr*sigw)) * part%icbt &
         + sqrt(1.-rw**2) * rands(nrand)
    delz=part%turbvelw*sigw*part%ptstep 
  endif
  nrand=nrand+2

  ! Calculate new vertical position
  if (abs(delz).gt.part%hbl) then
    delz=mod(delz,part%hbl)
  endif

  ! Reflection and position updates
  if (delz.lt.-part%zmetres) then         ! reflection at ground
    part%zmetres = -part%zmetres - delz
    part%icbt = -1
  else if (delz.gt.(part%hbl-part%zmetres)) then ! reflection at top
    part%zmetres = -part%zmetres-delz+2.*part%hbl
    part%icbt = -1
  else                         ! no reflection
    part%zmetres = part%zmetres+delz
    part%icbt = 1
  endif

end subroutine flexpart_diffusion_within_abl

subroutine flexpart_diffusion_above_abl(part, pextra, dt_remaining)

  USE particleML, only: extraParticle, Particle

  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle 
  type(extraParticle), intent(inout) :: pextra

  real, intent(in) :: dt_remaining

  ! turbulence factors for the troposphere
  real, parameter :: d_trop = 50.0
  real :: uxscale

  call set_turbulence_timestep(part, dt_remaining, .false.)

  ! assume within troposphere
  uxscale=sqrt(2.*d_trop/part%ptstep)
  if (nrand+1.gt.max_rands) nrand=1
  part%turbvelu=rands(nrand)*uxscale
  part%turbvelv=rands(nrand+1)*uxscale 
  nrand=nrand+2
  part%turbvelw=0

end subroutine flexpart_diffusion_above_abl

subroutine random_walk_name(part, pextra, dt_remaining, adaptive)
  USE particleML, only: extraParticle, Particle

  !> particle with information
  type(Particle), intent(inout) :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  real, intent(in) :: dt_remaining
  logical, intent(in) :: adaptive

  real :: kappa ! Von Karmans constant
  real :: sigu, sigv, sigw ! Turbulent velocity standard deviations
  real :: tlu, tlv, tlw ! Lagrangian timescales
  real :: delz ! Turbulent vertical displacement (m)
  real :: zeta
  real :: wst, ust, ol, h, z
  real :: eps, c
  real :: dttlu, dttlv, dttlw
  real :: dsigwdz, dsigw2dz
  real :: ru, rv, rw
  real :: sigw2
  real :: density_corr
  real :: zc

  ust = pextra%ust
  wst = pextra%wst
  ol = pextra%ol
  z = part%zmetres
  h = part%hbl
  zeta = max(z/h, 1e-3)
  zc = zeta * h ! clamped particle height
  c = 2.0
  kappa = 0.4

  if (density_correction) then
    density_corr = pextra%rhograd/pextra%rho
  else
    density_corr = 0
  endif
  density_corr = 0

  if (well_mixed_test) then
    ol = -100
    ust = 0.5
    wst = ust * (h / (kappa * abs(ol)))**(1.0/3.0)
  endif

  if (zeta < 1.0) then
    if (ol <= 0.0) then 
      ! Unstable BL
      sigu = sqrt(0.4 * wst**2 + 4*ust**2 * (1- zeta)**1.5)
      sigu = max(sigu, 0.25)
      sigv = sigu

      sigw2 = 1.2 * wst**2 * zeta**0.666 * (1 - zeta) + 1.69 * ust**2 * (1 - zeta)**1.5
      sigw = sqrt(sigw2)

      dsigw2dz = (1/h) * (1.2 * wst**2 * (2.0/3.0 * zeta**(-1.0/3.0) * (1 - zeta) - zeta**(2.0/3.0)) &
            - 2.535 * ust**2 * (1 - zeta)**(0.5))
      dsigwdz = 0.5 * dsigw2dz / sigw

      eps = (1.5 - 1.2*zeta**0.3333)*wst**3/h &
            + ust**3/(kappa*zc)*(1.0 - zeta)

      eps = max(1.0e-5, eps)

      tlu = 2.0 * sigu**2 / (c * eps)
      tlv = tlu
      tlw = 2.0 * sigw**2 / (c * eps)

      if (sigw < 0.1) then
        sigw = 0.1
        dsigwdz = 0
      endif

    else 
      ! Stable BL
      sigu = 2.0*ust*(1.0-zeta)**0.75
      sigu = max(sigu,0.25)
      sigv = sigu

      eps = (ust**3 / (kappa * zc) + (4 * ust**3) / (kappa * ol)) * (1 - zeta)
      eps = max(1.0e-5, eps)

      sigw = 1.3*ust*(1.0-zeta)**0.75
      dsigwdz = -0.975 * ust / h  * (1 - zeta)**(-0.25)

      tlu = 2.0 * sigu**2 / (c * eps)
      tlv = tlu
      tlw = 2.0 * sigw**2 / (c * eps)

      if (sigw < 0.1) then
        sigw = 0.1
        dsigwdz = 0
      endif

    endif

    ! Clamp lagrangian timescales
    tlu=max(10.,tlu)
    tlv=max(10.,tlv)
    tlw=max(20.,tlw)

  else
    sigu = 0.25
    sigv = sigu
    sigw = 0.1
    tlu = 300.0
    tlv = tlu
    tlw = 100.0
    dsigwdz = 0.0
  endif

  if (turb_homogeneous) dsigwdz=0

  part%tlw = tlw
  part%sigw = sigw
  part%dsigwdz = dsigwdz

  ! Based on turbulence parameters, calculate size of integration time step
  call set_turbulence_timestep(part, dt_remaining, adaptive)

  dttlu = part%ptstep/tlu
  dttlv = part%ptstep/tlv

  ! Calculate turbulent horizontal velocities
  if (nrand+1.gt.max_rands) nrand=1
  if (dttlu.lt.0.5) then
    part%turbvelu = (1.0-dttlu) * part%turbvelu + rands(nrand) * sigu * sqrt(2.0*dttlu)
  else
    ru = exp(-dttlu)
    part%turbvelu = ru * part%turbvelu + rands(nrand) * sigu * sqrt(1.0-ru**2)
  endif

  if (dttlv.lt.0.5) then
    part%turbvelv = (1.0-dttlv) * part%turbvelv + rands(nrand+1) * sigv * sqrt(2.0*dttlv)
  else
    rv = exp(-dttlv)
    part%turbvelv = rv * part%turbvelv + rands(nrand+1) * sigv * sqrt(1.0-rv**2)
  endif
  nrand = nrand + 2

  dttlw = part%ptstep/tlw

  if (nrand.gt.max_rands) nrand=1
  ! Calculate turbulent vertical velocity
  if (dttlw.lt.0.5) then
    part%turbvelw = (1.0-dttlw)*part%turbvelw + part%ptstep*(dsigwdz+density_corr*sigw) &
                    + sqrt(2.0*dttlw)*rands(nrand)
  else
    rw = exp(-dttlw)
    part%turbvelw = rw*part%turbvelw + tlw*(1.0-rw)*(dsigwdz+density_corr*sigw) &
                    + sqrt(1.0-rw**2)*rands(nrand)
  endif
  nrand = nrand + 1

  ! Calculate new vertical position
  delz=part%turbvelw*sigw*part%ptstep

  ! Reflection and position updates
  if (delz.lt.-part%zmetres) then ! reflection at ground
    part%zmetres = -part%zmetres - delz
    part%turbvelw = -part%turbvelw
  else ! no reflection
    part%zmetres = part%zmetres+delz
  endif

end subroutine random_walk_name

subroutine random_walk_name_fixed_below(part, pextra, dt_remaining, adaptive)
  use particleML, only: extraParticle, Particle
  implicit none

  type(Particle), intent(inout) :: part
  type(extraParticle), intent(inout) :: pextra

  real, intent(in) :: dt_remaining
  logical, intent(in) :: adaptive

  real :: kappa, c
  real :: sigu, sigv, sigw, sigw2
  real :: tlu, tlv, tlw
  real :: delz
  real :: zeta, z, h
  real :: wst, ust, ol
  real :: eps
  real :: dttlu, dttlv, dttlw
  real :: dsigw2dz
  real :: ru, rv, rw
  real :: density_corr
  real :: znew, zc

  ust = pextra%ust
  wst = pextra%wst
  ol = pextra%ol
  z = part%zmetres
  h = part%hbl
  zeta = min(max(z/h, 1.0e-3), 1.0 - 1.0e-3)
  zc = zeta * h ! clamped particle height
  c = 2.0
  kappa = 0.4

  if (well_mixed_test) then
    ol  = -500.0
    ust = 0.4
    wst = ust * (h / (kappa * abs(ol)))**(1.0/3.0)
  endif

  if (density_correction) then
    density_corr = pextra%rhograd / pextra%rho
  else
    density_corr = 0.0
  endif

  density_corr = 0.0

  if (ol <= 0.0) then
    ! Unstable BL
    sigu = sqrt(0.4 * wst**2 + 4*ust**2 * (1- zeta)**1.5)
    sigu = max(sigu, 0.25)
    sigv = sigu

    sigw2 = 1.2 * wst**2 * zeta**0.666 * (1 - zeta) + 1.69 * ust**2 * (1 - zeta)**1.5
    sigw = sqrt(sigw2)

    dsigw2dz = (1/h) * (1.2 * wst**2 * (2.0/3.0 * zeta**(-1.0/3.0) * (1 - zeta) - zeta**(2.0/3.0)) &
           - 2.535 * ust**2 * (1 - zeta)**(0.5))

    eps = (1.5 - 1.2*zeta**0.3333)*wst**3/h &
          + ust**3/(kappa*zc)*(1.0 - zeta)
    eps = max(1.0e-5, eps)

    tlu = 2.0 * sigu**2 / (c * eps)
    tlv = tlu
    tlw = 2.0 * sigw**2 / (c * eps)

    if (sigw < 0.1) then
      sigw = 0.1
      dsigw2dz = 0
    endif

  else
    ! Stable BL
    sigu = 2.0*ust*(1.0-zeta)**0.75
    sigu = max(sigu,0.25)
    sigv = sigu

    sigw = 1.3*ust*(1.0-zeta)**0.75

    eps = (ust**3 / (kappa * zc) + (4 * ust**3) / (kappa * ol)) * (1 - zeta)
    eps = max(1.0e-5, eps)

    dsigw2dz = -2.535 * ust**2 / h * (1 - zeta)**0.5

    tlu = 2.0 * sigu**2 / (c * eps)
    tlv = tlu
    tlw = 2.0 * sigw**2 / (c * eps)

    if (sigw < 0.1) then
      sigw = 0.1
      dsigw2dz = 0
    endif

  endif

  ! Clamp lagrangian timescales
  tlu=max(10.,tlu)
  tlv=max(10.,tlv)
  tlw=max(20.,tlw)

  part%tlw = tlw
  part%sigw = sigw

  call set_turbulence_timestep(part, dt_remaining, adaptive)

  dttlu = part%ptstep / tlu
  dttlv = part%ptstep / tlv
  dttlw = part%ptstep / tlw

  if (nrand+1 .gt. max_rands) nrand = 1
  if (dttlu < 0.5) then
    part%turbvelu = (1.0-dttlu) * part%turbvelu + rands(nrand) * sigu * sqrt(2.0*dttlu)
  else
    ru = exp(-dttlu)
    part%turbvelu = ru * part%turbvelu + rands(nrand) * sigu * sqrt(1.0-ru**2)
  endif
  if (dttlv < 0.5) then
    part%turbvelv = (1.0-dttlv) * part%turbvelv + rands(nrand+1) * sigv * sqrt(2.0*dttlv)
  else
    rv = exp(-dttlv)
    part%turbvelv = rv * part%turbvelv + rands(nrand+1) * sigv * sqrt(1.0-rv**2)
  endif
  nrand = nrand + 2

  if (nrand .gt. max_rands) nrand = 1
  rw = exp(-dttlw)
  part%turbvelw = rw * part%turbvelw &
                + tlw * (1.0 - rw) * (dsigw2dz + density_corr * sigw**2) &
                + sqrt(1.0 - rw**2) * sigw * rands(nrand)
  nrand = nrand + 1

  delz = part%turbvelw * part%ptstep
  znew = part%zmetres + delz

  do
    if (znew < 0.0) then
      znew = -znew
      part%turbvelw = -part%turbvelw
    else if (znew > h) then
      znew = 2.0 * h - znew
      part%turbvelw = -part%turbvelw
    else
      exit
    end if
  end do

  part%zmetres = znew

end subroutine random_walk_name_fixed_below

subroutine random_walk_name_fixed_below2(part, pextra, dt_remaining, adaptive)
  use particleML, only: extraParticle, Particle
  implicit none

  type(Particle), intent(inout) :: part
  type(extraParticle), intent(inout) :: pextra

  real, intent(in) :: dt_remaining
  logical, intent(in) :: adaptive

  real :: kappa, c
  real :: sigu, sigv, sigw, sigw2
  real :: tlu, tlv, tlw
  real :: delz
  real :: zeta, z, h
  real :: wst, ust, ol
  real :: eps
  real :: dttlu, dttlv, dttlw
  real :: dsigwdz, dsigw2dz
  real :: ru, rv, rw
  real :: density_corr
  real :: znew, zc

  ust = pextra%ust
  wst = pextra%wst
  ol = pextra%ol
  z = part%zmetres
  h = part%hbl
  zeta = min(max(z/h, 1.0e-3), 1.0 - 1.0e-3)
  zc = zeta * h ! clamped particle height
  c = 2.0
  kappa = 0.4

  if (well_mixed_test) then
    ol  = -500.0
    ust = 0.4
    wst = ust * (h / (kappa * abs(ol)))**(1.0/3.0)
  endif

  if (density_correction) then
    density_corr = pextra%rhograd / pextra%rho
  else
    density_corr = 0.0
  endif

  density_corr = 0.0

  if (ol <= 0.0) then
    ! Unstable BL

    sigu = sqrt(0.4 * wst**2 + 4*ust**2 * (1- zeta)**1.5)
    sigu = max(sigu, 0.25)
    sigv = sigu

    sigw2 = 1.2 * wst**2 * zeta**0.666 * (1 - zeta) + 1.69 * ust**2 * (1 - zeta)**1.5
    sigw = sqrt(sigw2)

    dsigw2dz = (1/h) * (1.2 * wst**2 * (2.0/3.0 * zeta**(-1.0/3.0) * (1 - zeta) - zeta**(2.0/3.0)) &
           - 2.535 * ust**2 * (1 - zeta)**(0.5))
    dsigwdz = 0.5 * dsigw2dz / sigw

    eps = (1.5 - 1.2*zeta**0.3333)*wst**3/h &
          + ust**3/(kappa*zc)*(1.0 - zeta)
    eps = max(1.0e-5, eps)

    tlu = 2.0 * sigu**2 / (c * eps)
    tlv = tlu
    tlw = 2.0 * sigw**2 / (c * eps)

    if (sigw < 0.1) then
      sigw = 0.1
      dsigwdz = 0
    endif

  else
    ! Stable BL
    sigu = 2.0*ust*(1.0-zeta)**0.75
    sigu = max(sigu,0.25)
    sigv = sigu

    sigw = 1.3*ust*(1.0-zeta)**0.75

    eps = (ust**3 / (kappa * zc) + (4 * ust**3) / (kappa * ol)) * (1 - zeta)
    eps = max(1.0e-5, eps)

    dsigwdz = -0.975 * ust / h  * (1 - zeta)**(-0.25)

    tlu = 2.0 * sigu**2 / (c * eps)
    tlv = tlu
    tlw = 2.0 * sigw**2 / (c * eps)

    if (sigw < 0.1) then
      sigw = 0.1
      dsigwdz = 0
    endif

  endif

  ! Clamp lagrangian timescales
  tlu=max(10.,tlu)
  tlv=max(10.,tlv)
  tlw=max(20.,tlw)

  part%tlw = tlw
  part%sigw = sigw

  call set_turbulence_timestep(part, dt_remaining, adaptive)

  dttlu = part%ptstep / tlu
  dttlv = part%ptstep / tlv
  dttlw = part%ptstep / tlw

  if (nrand+1 .gt. max_rands) nrand = 1
  if (dttlu < 0.5) then
    part%turbvelu = (1.0-dttlu) * part%turbvelu + rands(nrand) * sigu * sqrt(2.0*dttlu)
  else
    ru = exp(-dttlu)
    part%turbvelu = ru * part%turbvelu + rands(nrand) * sigu * sqrt(1.0-ru**2)
  endif
  if (dttlv < 0.5) then
    part%turbvelv = (1.0-dttlv) * part%turbvelv + rands(nrand+1) * sigv * sqrt(2.0*dttlv)
  else
    rv = exp(-dttlv)
    part%turbvelv = rv * part%turbvelv + rands(nrand+1) * sigv * sqrt(1.0-rv**2)
  endif
  nrand = nrand + 2

  if (nrand .gt. max_rands) nrand = 1
  rw = exp(-dttlw)
  part%turbvelw = rw * part%turbvelw &
                + tlw * (1.0 - rw) * (dsigwdz + density_corr * sigw) &
                + sqrt(1.0 - rw**2) * rands(nrand)
  nrand = nrand + 1

  delz = part%turbvelw * sigw * part%ptstep
  znew = part%zmetres + delz

  do
    if (znew < 0.0) then
      znew = -znew
      part%turbvelw = -part%turbvelw
    else if (znew > h) then
      znew = 2.0 * h - znew
      part%turbvelw = -part%turbvelw
    else
      exit
    end if
  end do

  part%zmetres = znew

end subroutine random_walk_name_fixed_below2

subroutine random_walk_name_fixed_above(part, pextra, dt_remaining, adaptive)
  use particleML, only: extraParticle, Particle
  implicit none

  type(Particle), intent(inout) :: part
  type(extraParticle), intent(inout) :: pextra

  real, intent(in) :: dt_remaining
  logical, intent(in) :: adaptive

  ! turbulence factors for the troposphere
  real, parameter :: ku_trop = 18.75
  real, parameter :: kz_trop = 1.0

  if (nrand+2.gt.max_rands) nrand=1
  part%turbvelu=rands(nrand)*sqrt(2.*ku_trop/part%ptstep)
  part%turbvelv=rands(nrand+1)*sqrt(2.*ku_trop/part%ptstep) 
  part%turbvelw=rands(nrand+2)*sqrt(2.*kz_trop/part%ptstep)
  nrand=nrand+3

end subroutine random_walk_name_fixed_above

subroutine variable_k_within_bl(part, pextra) 
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  real :: kappa ! Von Karmans constant
  real :: sigu, sigv, sigw ! Turbulent velocity standard deviations
  real :: tlu, tlv, tlw ! Lagrangian timescales
  real :: delz ! Turbulent vertical displacement (m)
  real :: scaled_height
  real :: wst, ust, ol, h
  real :: eps, c

  ! Dimensionless height, take middle value for homogeneous profiles
  scaled_height = 0.5

  kappa = 0.4
  c = 2 ! constant, values for this disagree. 3 from Sawford
  wst = pextra%wst
  ust = pextra%ust
  ol = pextra%ol
  h = part%hbl

  if (well_mixed_test) then
    ol = -100
    ust = 0.5
    wst = ust * (h / (kappa * abs(ol)))**(1.0/3.0)
  endif

  ! Case 1, Stable Conditions
  if (ol.gt.0.) then
    sigu=2.*ust*(1.-scaled_height) !. 7.20, Hanna
    sigu=max(sigu,0.25)
    sigv=sigu
    sigw=1.3*ust*(1.-scaled_height)
    sigw=max(sigw,0.1)

    ! Lagrangian timescales
    tlu=0.07*(part%hbl/sigv)*(scaled_height)**0.5 !. 7.22, Hanna
    tlv=tlu
    tlw=0.1*part%hbl/sigw*(scaled_height)**0.8
  
  ! Case 2 , Unstable Conditions
  else
    sigu = (0.4*wst**2 + (5 - 4*scaled_height)*ust**2)**0.5
    sigu=max(sigu,0.25)
    sigv = sigu
    sigw = (1.2 * wst**2  *(1-0.9*scaled_height)*(scaled_height)**0.666 + (1.8 - 1.4*scaled_height)*ust**2)**0.5
    sigw=max(sigw,0.1) 

    eps = (1.5 - 1.2 * (scaled_height)**0.333)*(wst**3/part%hbl) + (ust**3 * (1-0.8*scaled_height)/(kappa*(part%hbl*0.5)))

    tlu = 2*sigu**2 / (c*eps)
    tlv = tlu
    tlw = 2*sigw**2 / (c*eps)
  endif

  ! Clamp lagrangian timescales
  tlu=max(10.,tlu)
  tlv=max(10.,tlv)
  tlw=max(30.,tlw)

  ! Calculate Turbulent Velocities, Ryall and Maryon 1998
  if (nrand+2.gt.max_rands) nrand=1
  part%turbvelu = ((2*(sigu**2 * tlu))/tstep)**0.5 * rands(nrand)
  part%turbvelv = ((2*(sigv**2 * tlv))/tstep)**0.5 * rands(nrand+1)
  part%turbvelw = ((2*(sigw**2 * tlw))/tstep)**0.5 * rands(nrand+2)
  nrand=nrand+3

  delz=part%turbvelw*tstep 

  part%zmetres = part%zmetres + delz
  do while (part%zmetres < 0.0 .or. part%zmetres > part%hbl)
    if (part%zmetres < 0.0) part%zmetres = -part%zmetres
    if (part%zmetres > part%hbl) part%zmetres = 2.0*part%hbl - part%zmetres
  end do

end subroutine variable_k_within_bl

subroutine variable_k_above_bl(part, pextra) 
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  real :: sigu, sigv, sigw ! Turbulent velocity standard deviations
  real :: tlu, tlv, tlw ! Lagrangian timescales
  real :: delz ! Turbulent vertical displacement (m)

  sigu = 0.25
  sigv = sigu
  sigw = 0.1
  tlu = 300
  tlv = tlu
  tlw = 100

  if (nrand+2.gt.max_rands) nrand=1
  part%turbvelu = ((2*(sigu**2 * tlu))/tstep)**0.5 * rands(nrand)
  part%turbvelv = ((2*(sigv**2 * tlv))/tstep)**0.5 * rands(nrand+1)
  part%turbvelw = ((2*(sigw**2 * tlw))/tstep)**0.5 * rands(nrand+2)
  nrand=nrand+3

  ! Calculate new vertical position
  delz=part%turbvelw*tstep 

  ! Reflect at inversion (from above) so it is symmetric with BL routine
  if (part%zmetres+delz .lt. part%hbl) then
    part%zmetres = 2*part%hbl - (part%zmetres + delz)
  else
    part%zmetres = part%zmetres + delz
  endif

end subroutine variable_k_above_bl

subroutine variable_k_inhomog(part, pextra)
  use particleML, only: extraParticle, Particle
  implicit none

  type(Particle), intent(inout)      :: part
  type(extraParticle), intent(inout) :: pextra

  real :: sigu, sigv, sigw
  real :: tlu, tlv, tlw
  real :: Kz, Ku, Kv, dKdz
  real :: delz
  real :: ust, wst, ol
  real :: c, kappa
  real :: eps, depsdz
  real :: zeta, S, dSdz, h, z
  real :: B

  ! turbulence factors for the troposphere
  real, parameter :: ku_trop = 18.75
  real, parameter :: kz_trop = 1.0

  ust = pextra%ust
  wst = pextra%wst
  ol = pextra%ol
  z = part%zmetres
  h = part%hbl
  zeta = z/h 
  c = 2.0
  kappa = 0.4

  if (well_mixed_test) then
    ol = -100.0
    ust = 0.6
    wst = ust * (h / (kappa * abs(ol)))**(1.0/3.0)
  endif

  if (ol <= 0.0) then

    sigu = sqrt(0.4 * wst**2 + 4*ust**2 * (1- zeta)**1.5)
    sigu = max(sigu, 0.25)
    sigv = sigu

    ! Let S = sigw**2
    S = 1.2*wst**2 * zeta**0.6666 * (1.0 - zeta) + 1.69*ust**2 * (1.0 - zeta)**1.5
    dSdz = (1/h) * (1.2 * wst**2 * (2.0/3.0 * zeta**(-1.0/3.0) * (1 - zeta) - zeta**(2.0/3.0)) &
          - 2.535 * ust**2 * (1 - zeta)**(0.5))

    ! Turbulence dissipation rate
    eps = (1.5 - 1.2*zeta**0.3333)*wst**3/h &
        + ust**3/(kappa*z)*(1.0 - zeta)

    depsdz = -1/(h**2) * (0.4 * wst**3 * zeta**(-2.0/3.0) + ust**3/kappa * zeta**(-2))

    tlu = (2 * sigu**2) / (c * eps)
    tlv = tlu

    Ku = sigu**2 * tlu
    Ku = max(Ku, ku_trop)
    Kv = Ku

    Kz = 2.0*S*S/(c*eps)
    dKdz = Kz * (2*dSdz / S - depsdz / eps)
    
    if (Kz < 1.0) then
      Kz = 1.0
      dKdz = 0.0
    endif

  else
    ! Stable Conditions
    sigu = 2.0*ust*(1.0-zeta)**0.75
    sigu = max(sigu,0.25)
    sigv = sigu
    sigw = 1.3*ust*(1.0-zeta)**0.75

    eps = (ust**3 / (kappa * z) + (4 * ust**3) / (kappa * ol)) * (1 - zeta)
  
    tlu = (2 * sigu**2) / (c * eps)
    tlv = tlu

    B = 2.0*(1.3**4)*ust*kappa*ol/c     

    Ku = sigu**2 * tlu
    Ku = max(Ku, ku_trop)
    Kv = Ku

    ! Kz = sigw**4 / (C * eps), expanded and simplified
    Kz = B * (1 - zeta)**2 * (z/(ol + 4*z))
    dKdz = Kz * ((1/z) - 2/(h-z) - 4/(ol + 4*z))

    if (Kz < 1.0) then
      Kz = 1.0
      dKdz = 0.0
    endif
  endif
  
  if (turb_homogeneous) dKdz = 0.0

  if (nrand+2 .gt. max_rands) nrand = 1
  part%turbvelu = sqrt(2.0*(Ku)/tstep) * rands(nrand)
  part%turbvelv = sqrt(2.0*(Kv)/tstep) * rands(nrand+1)
  part%turbvelw = dKdz + sqrt(2.0*(Kz)/tstep) * rands(nrand+2)
  nrand = nrand + 3

  delz = part%turbvelw * tstep

  part%zmetres = part%zmetres + delz

  do while (part%zmetres < 0.0 .or. part%zmetres > part%hbl)
    if (part%zmetres < 0.0) then
      part%zmetres = -part%zmetres
    else if (part%zmetres > part%hbl) then
      part%zmetres = 2.0 * part%hbl - part%zmetres
    endif
  end do

end subroutine variable_k_inhomog

subroutine variable_k_inhomog_above(part, pextra)

  USE particleML, only: extraParticle, Particle

  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle 
  type(extraParticle), intent(inout) :: pextra

  ! turbulence factors for the troposphere
  real, parameter :: ku_trop = 18.75
  real, parameter :: kz_trop = 1.0

  if (nrand+1.gt.max_rands) nrand=1
  part%turbvelu=rands(nrand)*sqrt(2.*ku_trop/part%ptstep)
  part%turbvelv=rands(nrand+1)*sqrt(2.*ku_trop/part%ptstep) 
  part%turbvelw=rands(nrand+2)*sqrt(2.*kz_trop/part%ptstep)
  nrand=nrand+3
   
end subroutine variable_k_inhomog_above

subroutine fixed_k_within_bl(part)
  USE particleML, only: extraParticle, Particle

  !> particle with information
  type(Particle), intent(inout) :: part

  ! Locals
  real :: rnd(1)
  integer, parameter :: hor_diffu = 5300 ! m^2 s^-1 (BL horizontal diffusion)
  real, parameter :: K_ft = 1.5 ! m^2 s^-1 (FT vertical diffusion)
  real :: dt, delz

  dt = tstep

  if (nrand+2.gt.max_rands) nrand=1
  part%turbvelu = ((2*hor_diffu)/tstep)**0.5 * rands(nrand)
  part%turbvelv = ((2*hor_diffu)/tstep)**0.5 * rands(nrand+1)
  part%turbvelw = sqrt((2.0 * K_ft) / dt) * rands(nrand + 2)
  nrand=nrand+3

  call random_number(rnd)
  part%zmetres = rnd(1) * part%hbl

  ! Apply a vertical step equal in size to FT diffusion
  delz = part%turbvelw * dt

  ! Reflection and position updates
  if (delz.lt.-part%zmetres) then         ! reflection at ground
    part%zmetres = -part%zmetres - delz
  else if (delz.gt.(part%hbl-part%zmetres)) then ! reflection at top
    part%zmetres = -part%zmetres-delz+2.*part%hbl
  else                         ! no reflection
    part%zmetres = part%zmetres+delz
  endif

end subroutine fixed_k_within_bl

subroutine fixed_k_above_bl(part)
  use particleML, only: extraParticle, Particle
  implicit none

  type(Particle),      intent(inout) :: part

  real, parameter :: hor_diffu_ft = 5300.0/4.0 ! m^2 s^-1
  real, parameter :: K_ft = 1.5 ! m^2 s^-1
  real :: dt, delz

  dt = tstep

  ! Turbulent velocities
  if (nrand + 2 .gt. max_rands) nrand = 1
  part%turbvelu = sqrt((2.0 * hor_diffu_ft) / dt) * rands(nrand)
  part%turbvelv = sqrt((2.0 * hor_diffu_ft) / dt) * rands(nrand+1)
  part%turbvelw = sqrt((2.0 * K_ft) / dt) * rands(nrand+2)
  nrand = nrand + 3

  ! Apply a vertical step
  delz = part%turbvelw * dt

  if (part%zmetres+delz .lt. part%hbl) then
    part%zmetres = 2*part%hbl - (part%zmetres + delz)
  else
    part%zmetres = part%zmetres + delz
  endif

end subroutine fixed_k_above_bl

subroutine tke_diffusion(part, pextra, dt_remaining, adaptive)
  USE particleML, only: extraParticle, Particle
  USE snapfldML, only: hinterf
  USE snapdimML, only: nk

  type(Particle), intent(inout)  :: part
  type(extraParticle), intent(inout) :: pextra

  real, intent(in) :: dt_remaining
  logical, intent(in) :: adaptive

  integer :: i, j, k
  real :: sigu, sigv, sigw
  real :: part_z
  real :: tlu, tlv, tlw
  real :: ru, rv, rw
  real :: delz, rhoaux
  real :: dttlw

  i = part%x
  j = part%y

  ! Identify model layer of particle
  do k = 1, nk-1
    if (part%zmetres <= hinterf(i, j, k)) exit
  end do
  part_z = part%zmetres
  
  call calc_turb_params_tke(i, j, k, sigu, sigv, sigw, tlu, tlv, tlw)

  part%tlw = tlw

  rhoaux = pextra%rhograd/pextra%rho

  call set_turbulence_timestep(part, dt_remaining, adaptive)

  dttlw = part%ptstep/tlw

  if (nrand+2.gt.max_rands) nrand=1

  ! Vertical turbulence
  if (dttlw.lt..5) then
    part%turbvelw = (1.-dttlw)*part%turbvelw+rands(nrand)*sqrt(2.*dttlw)+part%ptstep*rhoaux*sigw
  else
      rw = exp(-dttlw)
    part%turbvelw = rw*part%turbvelw+rands(nrand)*sqrt(1.-rw**2)+tlw*(1.-rw)*rhoaux*sigw
  end if

  delz = part%turbvelw * sigw * part%ptstep

  call vertical_reflection_step(i, j, k, part_z, tlu, tlv, tlw, sigu, sigv, sigw, part%turbvelw, delz)

  part%zmetres = part_z

  ! Horizontal turbulence
  ru = exp(-part%ptstep / tlu)
  rv = exp(-part%ptstep / tlv)
  part%turbvelu = ru * part%turbvelu + rands(nrand+1) * sqrt(1.0 - ru**2) * sigu
  part%turbvelv = rv * part%turbvelv + rands(nrand+2) * sqrt(1.0 - rv**2) * sigv

  nrand = nrand + 3

end subroutine tke_diffusion

subroutine vertical_reflection_step(i, j, k, part_z, tlu, tlv, tlw, sigu, sigv, sigw, part_turbvelw, delz)

  USE snapfldML, only: hlevel2, hinterf
  USE snapdimML, only: nk

  integer, intent(in) :: i, j
  integer, intent(inout) :: k
  real, intent(inout) :: tlu, tlv, sigu, sigv, sigw, tlw
  real(4), intent(inout) :: part_z, delz
  real(8), intent(inout) :: part_turbvelw

  integer :: dir, k_c
  real :: z_bot, z_top
  real :: ts, ratio, sigw_c, tlw_c
  real :: rnd(1)
  logical :: reflect

  do
    call random_number(rnd)
    ! Calculate layer boundaries
    z_bot = hinterf(i,j,k-1)
    z_top = hinterf(i,j,k)

    ! Calculate time scale and adjust displacement
    ts = delz / (part_turbvelw * sigw)

    ! Determine if crossing will occur
    if (part_z + delz < z_bot) then
      k_c = k - 1
      dir = -1
    else if (part_z + delz > z_top) then
      if (k == nk) then
        part_z = hlevel2(i,j,nk-1)
        exit
      end if
      k_c = k + 1
      dir = 1

    else
      part_z = part_z + delz
      exit
    end if

    reflect = .false.

    ! Find time until particle reaches boundary
    ! Also set particle to interface height
    if (dir == 1) then
      ts = ts * (1 - ((z_top - part_z) / delz))
      part_z = z_top
    else 
      ts = ts * (1 - ((z_bot - part_z) / delz))
      part_z = z_bot
    endif

    ! Handle model boundary reflection
    if (k == 1 .AND. dir == -1) then
      reflect = .true.
    else
      ! Update turbulence parameters for new layer
      call calc_turb_params_tke(i, j, k_c, sigu, sigv, sigw_c, tlu, tlv, tlw_c)

      ! Compute transmission probability
      ! If next layer is more or equally turbulent then transfer
      ! If next layer is less turbulent, have probability of transfer
      ratio = sigw_c / sigw

      ! Draw random number for transmission
      if (rnd(1) >= ratio) reflect = .true.

    end if

    ! Apply reflection or transmission
    if (reflect) then
      part_turbvelw = -part_turbvelw
    else
      sigw = sigw_c
      k = k_c
      tlw = tlw_c
    end if

    delz = part_turbvelw * sigw * ts

  end do

end subroutine vertical_reflection_step

subroutine calc_turb_params_tke(i, j, k, sigu, sigv, sigw, tlu , tlv, tlw)

  USE particleML, only: extraParticle, Particle
  USE snapfldML, only: tke, hlevel2, dudxprof, dvdyprof, dwdzprof, pttprof, pttrefprof
  USE snapdimML, only: nk

  integer, intent(in) :: i, j, k
  real, intent(out) :: tlu, tlv, tlw, sigu, sigv, sigw

  integer :: kp
  real :: tke_z, yl, yl_up, yl_down, sum, e1
  real :: fu2, fv2, fw2
  real, parameter :: g = 9.80665

  tke_z = max(tke(i, j, k), 1e-6)

  ! Compute BL89 mixing length
  e1 = -g / pttrefprof(k) * (pttprof(k) - pttprof(k+1)) * (hlevel2(i, j, k+1) - hlevel2(i, j, k))
  if (e1 >= tke_z) then
    yl = hlevel2(i, j, k+1) - hlevel2(i, j, k)
  else
    ! Upward
    sum = 0.0
    kp = k+1
    do while (kp < nk)
      sum = sum - (g / pttrefprof(kp) * (pttprof(k) - pttprof(kp)) * (hlevel2(i, j, kp) - hlevel2(i, j, kp-1)))
      if (sum >= tke_z) exit
      kp = kp + 1
    end do
    yl_up = hlevel2(i, j, kp) - hlevel2(i, j, k)
    ! Downward
    sum = 0.0
    kp = k
    do while (kp > 2)
      sum = sum - (g / pttrefprof(kp) * (pttprof(kp) - pttprof(k+1)) * (hlevel2(i, j, kp+1) - hlevel2(i, j, kp)))
      if (sum >= tke_z) exit
      kp = kp - 1
    end do
    yl_down = hlevel2(i, j, k+1) - hlevel2(i, j, kp)
    yl = ((yl_up**(-2.0/3.0) + yl_down**(-2.0/3.0))/2.0)**(-3.0/2.0)
  end if

  ! Calculate Anisotropy fractions
  fu2 = max(0.,1.-(yl/5./sqrt(tke_z)*dudxprof(k)))/3.
  fv2 = max(0.,1.-(yl/5./sqrt(tke_z)*dvdyprof(k)))/3.
  fw2 = max(0.,(1.-yl/5./sqrt(tke_z)*dwdzprof(k)))/3.

  ! Variances
  sigu = sqrt(2.0 * tke_z * fu2)
  sigv = sqrt(2.0 * tke_z * fv2)
  sigw = sqrt(2.0 * tke_z * fw2)

  ! Timescales
  tlu = 2.0 * yl / max(sigu, 1e-6)
  tlv = 2.0 * yl / max(sigv, 1e-6)
  tlw = 2.0 * yl / max(sigw, 1e-6)
  tlu = max(10.0, tlu)
  tlv = max(10.0, tlv)
  tlw = max(30.0, tlw)

end subroutine calc_turb_params_tke

subroutine diffusion_fields
  use snapfldML, only: ps_io, hbl_io, surface_stress, hflux, tv, obukhov_l_io, u_star_io, w_star_io, tv
  use snapdimML, only: nx, ny
  use, intrinsic :: ieee_arithmetic

  real, parameter :: r=287, g=9.81, k=0.4, cpa=1004.6

  real :: rho_a(nx, ny)
  real :: fhsfc(nx, ny) ! Surface kinematic heat flux 

  ! Calculate surface air density
  rho_a = (ps_io*100) / (tv(:,:,2) * r)

  ! Calculate friction velocity
  u_star_io = sqrt(surface_stress/(rho_a))

  ! surface kinematic heat flux = H0/(rho*cpa)
  fhsfc = hflux/(rho_a*cpa)

  ! Calculate the obukhov length
  obukhov_l_io = - (tv(:,:,2) * u_star_io**3)/(k*g*fhsfc)

  if (bl_id == 1) then
    hbl_io=600
  endif

  ! Calculate the convective velocity scale
  w_star_io = u_star_io*(hbl_io/(k * abs(obukhov_l_io)))**(1.0/3.0)

  where (ieee_is_nan(w_star_io))
    w_star_io = 0.0
  end where
  
end subroutine diffusion_fields

subroutine air_density
  use snapfldML, only: spec_humid, ps_io, hlevel_io, t2, rho_io, rhograd_io, pressures, tv
  use snapdimML, only: nx, ny, nk
  use snapgrdML, only: alevel, blevel

  real, parameter :: r = 287.0

  integer :: i, j, k

  ! Compute virtual temperature where defined (k >= 2)
  tv = t2 * (1.0 + 0.608 * spec_humid)

  ! Surface virtual temperature is undefined so copy from level 2
  do j = 1, ny
    do i = 1, nx
      tv(i,j,1) = tv(i,j,2)
    end do
  end do

  ! Compute pressure at full levels (including the surface level)
  do j = 1, ny
    do i = 1, nx
      do k = 1, nk
        pressures(i,j,k) = alevel(k) * 100 + blevel(k) * ps_io(i,j) * 100.0
      end do
    end do
  end do

  ! Density 
  do j = 1, ny
    do i = 1, nx
      do k = 2, nk
        rho_io(i,j,k) = pressures(i,j,k) / (r * tv(i,j,k))
      end do
      ! Fill synthetic surface 
      rho_io(i,j,1) = rho_io(i,j,2)
    end do
  end do

  ! Interior points
  do k = 3, nk-1
    do j = 1, ny
      do i = 1, nx
        rhograd_io(i,j,k) = (rho_io(i,j,k+1) - rho_io(i,j,k-1)) / &
                        (hlevel_io(i,j,k+1) - hlevel_io(i,j,k-1))
      end do
    end do
  end do

  ! Bottom boundary, forward difference
  do j = 1, ny
    do i = 1, nx
      rhograd_io(i,j,2) = (rho_io(i,j,3) - rho_io(i,j,2)) &
                       / (hlevel_io(i,j,3) - hlevel_io(i,j,2))
    end do
  end do

  ! Top boundary, backward difference
  do j = 1, ny
    do i = 1, nx
      rhograd_io(i,j,nk) = (rho_io(i,j,nk) - rho_io(i,j,nk-1)) &
                        / (hlevel_io(i,j,nk) - hlevel_io(i,j,nk-1))
    end do
  end do

end subroutine

subroutine align_turbvels(part, pextra)
  use particleML, only: extraParticle, Particle
  implicit none

  type(Particle),      intent(inout) :: part
  type(extraParticle), intent(in)    :: pextra

  real :: umean, vmean, mag_inv, cosphi, sinphi
  real :: du_old, dv_old
  real, parameter :: eps = 1.e-30

  umean = pextra%u
  vmean = pextra%v
  mag_inv = 1.0 / max(sqrt(umean*umean + vmean*vmean), eps)

  cosphi = umean * mag_inv
  sinphi = vmean * mag_inv

  du_old = part%turbvelu
  dv_old = part%turbvelv

  part%turbvelu = du_old * cosphi - dv_old * sinphi
  part%turbvelv = du_old * sinphi + dv_old * cosphi
end subroutine align_turbvels

subroutine eta_to_metres(part)

  use particleML, only: extraParticle, Particle
  use snapfldML,  only: hlevel2, ps2
  use snapgrdML,  only: vlevel, alevel, blevel
  use snapdimML,  only: nk

  implicit none

  type(Particle), intent(inout) :: part

  integer :: i, j, k
  real :: eta
  real :: z1, z2
  real :: p1, p2, px, ps
  real :: frac

  i = part%x
  j = part%y
  eta = part%z
  ps = ps2(i,j) * 100.0

  ! Find vertical layer in eta 
  do k = 1, nk-1
    if (eta > vlevel(k+1)) exit
  end do
  k = max(1, min(k, nk-1))

  ! Fraction in eta
  frac = (eta - vlevel(k)) / (vlevel(k+1) - vlevel(k))
  frac = max(0.0, min(1.0, frac))

  ! Heights
  z1 = hlevel2(i,j,k)
  z2 = hlevel2(i,j,k+1)

  ! Pressures
  p1 = alevel(k)*100 + blevel(k) * ps
  p2 = alevel(k+1)*100 + blevel(k+1) * ps

  ! Pressure at particle
  px = p1 * (1.0 - frac) + p2 * frac

  ! Log-pressure height interpolation
  if (p1 > 0.0 .and. p2 > 0.0) then
    part%zmetres = z1 + (z2 - z1) / log(p2/p1) * log(px/p1)
  else
    part%zmetres = z1 * (1.0 - frac) + z2 * frac
  end if

end subroutine eta_to_metres

subroutine metres_to_eta(part)

  use particleML, only: extraParticle, Particle
  use snapfldML,  only: hlevel2, ps2
  use snapgrdML,  only: vlevel, alevel, blevel
  use snapdimML,  only: nk

  implicit none

  type(Particle), intent(inout) :: part

  integer :: i, j, k
  real :: z
  real :: z1, z2
  real :: p1, p2, px, ps
  real :: frac

  i = part%x
  j = part%y
  z = part%zmetres
  ps = ps2(i,j) * 100.0

  do k = 1, nk-1
    if (z < hlevel2(i,j,k+1)) exit
  end do
  k = max(1, min(k, nk-1))

  z1 = hlevel2(i,j,k)
  z2 = hlevel2(i,j,k+1)

  p1 = alevel(k)*100 + blevel(k) * ps
  p2 = alevel(k+1)*100 + blevel(k+1) * ps

  if (p1 > 0.0 .and. p2 > 0.0) then
    px = p1 * exp( log(p2/p1) * (z - z1) / (z2 - z1) )
    frac = (px - p1) / (p2 - p1)
  else
    frac = (z - z1) / (z2 - z1)
  end if

  frac = max(0.0, min(1.0, frac))

  part%z = vlevel(k) * (1.0 - frac) + vlevel(k+1) * frac

end subroutine metres_to_eta

subroutine set_turbulence_timestep(part, dt_remaining, adaptive)
  use particleML, only: Particle
  implicit none

  type(Particle), intent(inout) :: part
  real, intent(in) :: dt_remaining
  logical, intent(in) :: adaptive

  if (adaptive) then
     part%ptstep = min(part%tlw, &
                      part%hbl / max(2.0 * abs(part%turbvelw * part%sigw), 1.e-5), &
                      0.5 / max(abs(part%dsigwdz), 1.e-6)) * 0.10
     part%ptstep = max(part%ptstep, dt_min)
     part%ptstep = min(part%ptstep, dt_remaining)
  else
     part%ptstep = dt_remaining
  end if
end subroutine

end module rwalkML
