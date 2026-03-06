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
  USE iso_fortran_env, only: real64

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
  real(real64), parameter :: entrainment = 0.1 ! Entrainment zone = 10%*h

  ! Values for random number generation
  integer, parameter :: max_rands=6000000
  integer :: nrand=1
  real :: rands(max_rands) ! initialise array to store random numbers

  real(real64), save, public :: a_in_bl = 0.5
  real(real64), save, public :: a_above_bl = 0.25
  real(real64), save, public :: b = 0.875
  logical, save, public :: turb_homogeneous = .TRUE.
  character(len=64), save, public :: diffusion_scheme = ''
  character(len=64), save, public :: bl_definition = ''
  character(len=64), save, public :: meteo_type = ''
  character(len=64), save, public :: record_stats = ''
  character(len=64), save, public :: entrainment_scheme = ''

  public rwalk_init, diffusion_fields, air_density, turbulence_master, eta_to_metres

  contains

!> Initialise constants needed for rwalk
subroutine rwalk_init(timestep)
  use init_random_seedML, only: generate_normal_randoms
!> time step in seconds (trajectory calculations)
  real, intent(in) :: timestep
  integer :: i, unit

  tfactor_v = timestep/tmix_v
  tsqrtfactor_v=sqrt(tfactor_v)
  tfactor_h = timestep/tmix_h
  tsqrtfactor_h=sqrt(tfactor_h)
  tstep = timestep

  ! l-eta above mixing height
  vrdbla = labove*tsqrtfactor_v

  if (diffusion_scheme=='random_walk_flexpart' .OR. diffusion_scheme=='random_walk_name' .OR. &
      diffusion_scheme=='variable_k' .OR. diffusion_scheme=='constant_k') then
    call generate_normal_randoms(rands, max_rands)
  endif

end subroutine

subroutine turbulence_master(blfullmix,part,pextra)
  USE particleML, only: extraParticle, Particle

    !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle (u, v, rmx, rmy)
  type(extraParticle), intent(inout) :: pextra
  !> full mixing in boundarylayer (true=old,false=new)
  logical, intent(in) :: blfullmix

  if (bl_definition == 'constant') then
    part%hbl = 600
  endif

  if (diffusion_scheme == 'random_walk_flexpart') then
    ! Convert particle height from eta to metres
    call eta_to_metres(part, pextra)

    ! Check if particle within abl
    if (part%z.gt.part%tbl) then
      call flexpart_diffusion_within_abl(part,pextra)
      call metres_to_eta(part, pextra)
    else
      call flexpart_diffusion_above_abl(part, pextra)
      !call metres_to_eta(part, pextra)
    endif
  elseif (diffusion_scheme == 'variable_k') then
    ! Convert particle height from eta to metres
    call eta_to_metres(part, pextra)
    ! Check if particle within abl
    if (part%z.gt.part%tbl) then
      call variable_k_name_within_bl(part, pextra)
      call metres_to_eta(part, pextra)
    else
      call variable_k_name_above_bl(part, pextra)
      call metres_to_eta(part, pextra)
    endif
  elseif (diffusion_scheme == 'random_walk_name') then
    ! Convert particle height from eta to metres
    call eta_to_metres(part, pextra)
    ! Check if particle within abl
    if (part%z.gt.part%tbl) then
      call name_random_walk_profile_within_bl(part, pextra)
      call metres_to_eta(part, pextra)
    else
      call name_random_walk_profile_above_bl(part, pextra)
      call metres_to_eta(part, pextra)
    endif
  elseif (diffusion_scheme == 'constant_k') then
    ! Convert particle height from eta to metres
    call eta_to_metres(part, pextra)
    ! Check if particle within abl
    if (part%z.gt.part%tbl) then
      call constant_k_name_within_bl(part, pextra)
      call metres_to_eta(part, pextra)
    else
      call constant_k_name_above_bl(part, pextra)
      call metres_to_eta(part, pextra)
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

  if (bl_definition == 'constant') then
    part%hbl = 600
    part%tbl = 0.929817438
  endif

! horizontal diffusion
  if (part%z > part%tbl) then ! in boundary layer
    a = a_in_bl
  else ! above boundary layer
    a = a_above_bl
  endif

  vabs = hypot(pextra%u, pextra%v)
  rl = 2*a*((vabs*tmix_h)**b) * tsqrtfactor_h ! sqrt error/sigma propagation
  part%x = part%x + rl*rnd(1)*pextra%rmx
  part%y = part%y + rl*rnd(2)*pextra%rmy


! vertical diffusion
  if (part%z <= part%tbl) then ! Above boundary layer
      part%z = part%z + vrdbla*rnd(3)

  else ! In boundary layer
    bl_entrainment_thickness = (1.0 - part%tbl)*(1.+entrainment)
    if (blfullmix .or. (tsqrtfactor_v .gt. 1.0)) then ! full mixing
      part%z = 1.0 - bl_entrainment_thickness*(rnd(3)+0.5)
    else ! vertical mixing splittet in smaller time-steps   
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

subroutine flexpart_diffusion_within_abl(part, pextra) ! Based on Hanna1 from FLEXPART code, turbswitch FALSE
  USE particleML, only: extraParticle, Particle
  use snapfldML, only: rho, rhograd
  USE snapgrdML, only: ivlayer
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  integer :: i, j, k
  real :: sigu, sigv, sigw ! Turbulent velocity standard deviations
  real :: tlu, tlv, tlw ! Lagrangian timescales
  real :: ru, rv, rw ! Lagrangian timescales
  real :: delz ! Turbulent vertical displacement (m)
  real :: dsigw2dz 
  real :: dttlw
  real :: rhoaux ! Density correction factor
  real :: s1, s2
  integer :: part_vert_index
  real :: scaled_height
  real :: top_entrainment
  real :: wst ! convective scale velocity
  real :: ol ! obukhov length

  ! Dimensionless height 
  if (turb_homogeneous) then ! Take middle BL value if homogeneous
    scaled_height = 0.5
  else
    scaled_height = part%zmetres/part%hbl
  endif

  ! Get particle position
  i = part%x
  j = part%y
  part_vert_index = part%z*10000
  k = ivlayer(part_vert_index) ! Vertical layer of particle, want to use this instead of ivlevel

  wst = pextra%wst
  ol = pextra%ol

  ! Case 1, Neutral Conditions
  if (part%hbl/ABS(ol).lt.1.) then
    pextra%ust = max(1.e-4, pextra%ust)

    ! Eq. 7.25 Hanna 1982: 
    sigu = 2.0 * pextra%ust * EXP(-3.e-4*part%zmetres/pextra%ust)
    sigu = MAX(sigu, 1.e-5)

    ! Eq. 7.26 Hanna 1982:
    sigv = 1.3 * pextra%ust * EXP(-2.e-4*part%zmetres/pextra%ust)
    sigv=max(sigv,1.e-5)
    sigw=sigv

    ! Vertical gradient of sigw
    dsigw2dz=-6.76e-4*pextra%ust*exp(-4.e-4*part%zmetres/pextra%ust)

    ! Lagrangian timescales
    tlu=0.5*part%zmetres/sigw/(1.+1.5e-3*part%zmetres/pextra%ust)
    tlv=tlu
    tlw=tlu

  ! Case 2 , Unstable Conditions
  elseif (ol.lt.0.) then
    ! Eq. 4.15 Caughey 1982
    sigu=pextra%ust*(12.-0.5*part%hbl/ol)**0.33333
    sigu=MAX(sigu,1.e-6)
    sigv=sigu

    ! Eq. 7.15 Hanna 1982
    if (scaled_height.lt.0.03) then
      sigw=0.96*wst*(3*scaled_height-ol/part%hbl)**0.33333
      dsigw2dz=1.8432*wst*wst/part%hbl*(3*scaled_height-ol/part%hbl)**(-0.33333)
    else if (scaled_height.lt.0.4) then
      s1=0.96*(3*scaled_height-ol/part%hbl)**0.33333
      s2=0.763*scaled_height**0.175
      if (s1.lt.s2) then
        sigw=wst*s1
        dsigw2dz=1.8432*wst*wst/part%hbl*(3*scaled_height-ol/part%hbl)**(-0.33333)
      else
        sigw=wst*s2
        dsigw2dz=0.203759*wst*wst/part%hbl*scaled_height**(-0.65)
      endif
    else if (scaled_height.lt.0.96) then
      sigw=0.722*wst*(1-scaled_height)**0.207
      dsigw2dz=-.215812*wst*wst/part%hbl*(1-scaled_height)**(-0.586)
    else if (scaled_height.lt.1.00) then 
      sigw=0.37*wst
      dsigw2dz=0.
    endif
    sigw=max(sigw,1.e-6) 

    ! Determine average Lagrangian time scale
    ! Eq. 7.17 Hanna  1982
    tlu=0.15*part%hbl/sigu
    tlv=tlu
    if (part%zmetres.lt.ABS(ol)) then
      tlw=0.1*part%zmetres/(sigw*(0.55-0.38*ABS(part%zmetres/ol)))
    else if (scaled_height.lt.0.1) then
      tlw=0.59*part%zmetres/sigw
    else
      tlw=0.15*part%hbl/sigw*(1.-EXP(-5*scaled_height))
    endif

  ! Case 3, Stable Conditions 
  else
    sigu=2.*pextra%ust*(1.-scaled_height) !. 7.20, Hanna
    sigv=1.3*pextra%ust*(1.-scaled_height) !. 7.19, Hanna
    sigu=max(sigu,1.e-6)
    sigv=max(sigv,1.e-6)
    sigw=sigv !. 7.19, Hanna
    dsigw2dz=3.38*pextra%ust*pextra%ust*(scaled_height-1.)/part%hbl

    ! Lagrangian timescales
    tlu=0.15*part%hbl/sigu*(sqrt(scaled_height)) !. 7.22, Hanna
    tlv=0.467*tlu
    tlw=0.1*part%hbl/sigw*scaled_height**0.8

  endif

  ! Clamp lagrangian timescales
  tlu=max(10.,tlu)
  tlv=max(10.,tlv)
  tlw=max(30.,tlw)

  part%tlw = tlw

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

  ! Calculate new horizontal positions. Maybe should only update at the end?
  part%x = part%x + part%turbvelu * part%ptstep*pextra%rmx
  part%y = part%y + part%turbvelv * part%ptstep*pextra%rmy

  ! Factor for density correction, k+1 to avoid artificial first layer
  rhoaux=rhograd(i,j,k+1)/rho(i,j,k+1)

  ! ratio of time step to lagrangian timescale for autocorrelation
  dttlw = part%ptstep/tlw

  if (nrand+1.gt.max_rands) nrand=1
  ! Calculate turbulent vertical velocity
  rw=exp(-dttlw)
  part%turbvelw=(rw*part%turbvelw+rands(nrand)*sqrt(1.-rw**2)*sigw &
        +tlw*(1.-rw)*(dsigw2dz+rhoaux*sigw**2)) * part%icbt
  delz=part%turbvelw*part%ptstep 
  nrand=nrand+1

  ! Calculate new vertical position
  if (abs(delz).gt.part%hbl) then
    delz=mod(delz,part%hbl)
  endif

  if (entrainment_scheme=='snap') then
    ! SNAP ENTRAINMENT SCHEME --------------------------------------------------
    part%zmetres = part%zmetres + delz
    ! ... reflection from the ABL top
    ! ... but allow for entrainment
    ! top_entrainment 10% higher than tbl
    top_entrainment = part%hbl + part%hbl * 0.1
    if (part%zmetres > top_entrainment) then
      part%zmetres = 2.0*part%hbl - part%zmetres
    endif

    if (part%zmetres < 0.0) then
        ! reflect across ground (z=0)
        part%zmetres = -part%zmetres
    endif

    ! cap at entrainment height
    part%zmetres = min(part%zmetres, top_entrainment)

    ! enforce ground
    part%zmetres = max(part%zmetres, 0.0)
    ! --------------------------------------------------
  else
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
  endif

end subroutine flexpart_diffusion_within_abl

subroutine flexpart_diffusion_above_abl(part, pextra)

  USE particleML, only: extraParticle, Particle

  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle 
  type(extraParticle), intent(inout) :: pextra

  ! turbulence factors for the troposphere
  real :: d_trop=50.
  
  real :: uxscale

  ! assume within troposphere
  uxscale=sqrt(2.*d_trop/part%ptstep)
  if (nrand+1.gt.max_rands) nrand=1
  part%turbvelu=rands(nrand)*uxscale
  part%turbvelv=rands(nrand+1)*uxscale 
  nrand=nrand+2
  part%turbvelw=0

  part%x = part%x + part%turbvelu * part%ptstep*pextra%rmx
  part%y = part%y + part%turbvelv * part%ptstep*pextra%rmy

end subroutine flexpart_diffusion_above_abl

subroutine name_random_walk_profile_within_bl(part, pextra) 
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  integer :: i, j
  real :: k ! Von Karmans constant
  real :: sigu, sigv, sigw ! Turbulent velocity standard deviations
  real :: tlu, tlv, tlw ! Lagrangian timescales
  real :: delz ! Turbulent vertical displacement (m)
  real :: scaled_height
  real :: wst, ust
  real :: eps, c
  real :: dttlu, dttlv, dttlw
  real :: dsigwdz

  ! Dimensionless height 
  if (turb_homogeneous) then ! Take middle BL value if homogeneous
    scaled_height = 0.5
  else
    scaled_height = part%zmetres/part%hbl
  endif

  k = 0.4
  c = 2 ! constant, values for this disagree. 3 from Sawford
  wst = pextra%wst
  ust = pextra%ust
  
  ! Case 1, Stable Conditions
  if (pextra%ol.gt.0.) then
    sigu=2.*ust*(1.-scaled_height) 
    sigu=max(sigu,0.25)
    sigv=sigu
    sigw=1.3*ust*(1.-scaled_height)
    sigw=max(sigw,0.1)

    ! Lagrangian timescales
    tlu=0.07*(part%hbl/sigv)*(scaled_height)**0.5 
    tlv=tlu
    tlw=0.1*(part%hbl/sigw)*(scaled_height)**0.8

    dsigwdz = (-1.3 * ust) / part%hbl
  
  ! Case 2 , Unstable Conditions
  else
    sigu = (0.4*wst**2 + (5 - 4*scaled_height)*ust**2)**0.5
    sigu=max(sigu,0.25)
    sigv = sigu
    sigw = (1.2 * wst**2  *(1-0.9*scaled_height)*(scaled_height)**0.666 + (1.8 - 1.4*scaled_height)*ust**2)**0.5
    sigw=max(sigw,0.1) 

    eps = (1.5 - 1.2 * (scaled_height)**0.333)*((wst**3)/part%hbl) + (ust**3 * (1-0.8*scaled_height)/(k*part%zmetres))
    eps = max(1e-4, eps)

    tlu = 2*sigu**2 / (c*eps)
    tlv = tlu
    tlw = 2*sigw**2 / (c*eps)

    ! From FLEXPART Code ??
    dsigwdz=0.5/sigw/part%hbl*(-1.4*ust**2+wst**2* &
         (0.8*max(scaled_height,1.e-3)**(-.33333)-1.8*scaled_height**0.66666))
  endif

  if (turb_homogeneous) dsigwdz=0

  ! Clamp lagrangian timescales
  tlu=min(300.,tlu)
  tlv=min(300.,tlv)
  tlw=min(100.,tlw)

  tlu=max(20.,tlu)
  tlv=max(20.,tlv)
  tlw=max(20.,tlw)

  part%tlw = tlw

  dttlu = part%ptstep/tlu
  dttlv = part%ptstep/tlv

  ! Calculate Turbulent Velocities, Ryall and Maryon 1998
  if (nrand+1.gt.max_rands) nrand=1
  part%turbvelu = part%turbvelu*(1-(dttlu)) + (2*sigu**2 * dttlu)**0.5*rands(nrand)
  part%turbvelv = part%turbvelv*(1-(dttlv)) + (2*sigv**2 * dttlv)**0.5*rands(nrand+1)
  nrand=nrand+2

  ! ratio of time step to lagrangian timescale for autocorrelation
  dttlw = part%ptstep/tlw

  ! Calculate new horizontal positions
  part%x = part%x + part%turbvelu * part%ptstep*pextra%rmx
  part%y = part%y + part%turbvelv * part%ptstep*pextra%rmy

  if (nrand+1.gt.max_rands) nrand=1
  ! Calculate turbulent vertical velocity
  part%turbvelw=part%turbvelw*(1-(dttlw)) + (2*sigw**2 * dttlw)**0.5*rands(nrand) & 
                + (part%ptstep/sigw) * dsigwdz * (sigw**2 + part%turbvelw**2) 
  nrand=nrand+1

  ! Calculate new vertical position
  delz=part%turbvelw*part%ptstep 

  ! Reflection and position updates
  if (delz.lt.-part%zmetres) then         ! reflection at ground
    part%zmetres = -part%zmetres - delz
    part%turbvelw = -part%turbvelw
  else if (delz.gt.(part%hbl-part%zmetres)) then ! reflection at top
    part%zmetres = -part%zmetres-delz+2.*part%hbl
    part%turbvelw = -part%turbvelw
  else                         ! no reflection
    part%zmetres = part%zmetres+delz
  endif

end subroutine name_random_walk_profile_within_bl

subroutine name_random_walk_profile_above_bl(part, pextra) 
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  real :: sigu, sigv, sigw ! Turbulent velocity standard deviations
  real :: tlu, tlv, tlw ! Lagrangian timescales
  real :: delz ! Turbulent vertical displacement (m)
  real :: dttlw, dttlu, dttlv
  
  sigu = 0.25
  sigv = 0.25
  sigw = 0.1
  tlu = 300
  tlv = 300
  tlw = 100

  part%tlw = tlw

  dttlu = part%ptstep/tlu
  dttlv = part%ptstep/tlv

  ! Calculate Turbulent Velocities, Ryall and Maryon 1998
  if (nrand+1.gt.max_rands) nrand=1
  part%turbvelu = part%turbvelu*(1-(dttlu)) + (2*sigu**2 * dttlu)**0.5*rands(nrand)
  part%turbvelv = part%turbvelv*(1-(dttlv)) + (2*sigv**2 * dttlv)**0.5*rands(nrand+1)
  nrand=nrand+2

  ! ratio of time step to lagrangian timescale for autocorrelation
  dttlw = part%ptstep/tlw

  ! Calculate new horizontal positions
  part%x = part%x + part%turbvelu * part%ptstep*pextra%rmx
  part%y = part%y + part%turbvelv * part%ptstep*pextra%rmy

  if (nrand+1.gt.max_rands) nrand=1
  ! Calculate turbulent vertical velocity
  part%turbvelw=part%turbvelw*(1-(dttlw)) + (2*sigw**2 * dttlw)*rands(nrand)
  nrand=nrand+1

  ! Calculate new vertical position
  delz=part%turbvelw*part%ptstep  

  part%zmetres = part%zmetres + delz

end subroutine name_random_walk_profile_above_bl

subroutine variable_k_name_within_bl(part, pextra) 
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  integer :: k ! Von Karmans constant
  real :: sigu, sigv, sigw ! Turbulent velocity standard deviations
  real :: tlu, tlv, tlw ! Lagrangian timescales
  real :: delz ! Turbulent vertical displacement (m)
  real :: scaled_height
  real :: top_entrainment
  real :: wst, ust
  real :: eps, c

  ! Dimensionless height, take middle value for homogeneous profiles
  scaled_height = 0.5

  k = 0.4
  c = 2 ! constant, values for this disagree. 3 from Sawford
  wst = pextra%wst
  ust = pextra%ust
  
  ! Case 1, Stable Conditions
  if (pextra%ol.gt.0.) then
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

    eps = (1.5 - 1.2 * (scaled_height)**0.333)*(wst/part%hbl) + (ust**3 * (1-0.8*scaled_height)/(k*part%zmetres))

    tlu = 2*sigu**2 / (c*eps)
    tlv = tlu
    tlw = 2*sigw**2 / (c*eps)
  endif

  ! Clamp lagrangian timescales
  tlu=max(300.,tlu)
  tlv=max(300.,tlv)
  tlw=max(100.,tlw)

  ! Calculate Turbulent Velocities, Ryall and Maryon 1998
  if (nrand+3.gt.max_rands) nrand=1
  part%turbvelu = ((2*(sigu**2 * tlu))/tstep)**0.5 * rands(nrand)
  part%turbvelv = ((2*(sigv**2 * tlv))/tstep)**0.5 * rands(nrand+1)
  part%turbvelw = ((2*(sigw**2 * tlw))/tstep)**0.5 * rands(nrand+2)
  nrand=nrand+3

  ! Calculate new horizontal positions
  part%x = part%x + part%turbvelu * tstep*pextra%rmx
  part%y = part%y + part%turbvelv * tstep*pextra%rmy

  ! Calculate new vertical position
  delz=part%turbvelw*tstep 
  if (abs(delz).gt.part%hbl) then
    delz=mod(delz,part%hbl)
  endif

  if (entrainment_scheme=='snap') then
    ! SNAP ENTRAINMENT SCHEME --------------------------------------------------
    part%zmetres = part%zmetres + delz
    ! ... reflection from the ABL top
    ! ... but allow for entrainment
    ! top_entrainment 10% higher than tbl
    top_entrainment = part%hbl + part%hbl * 0.1
    if (part%zmetres > top_entrainment) then
      part%zmetres = 2.0*part%hbl - part%zmetres
    endif

    if (part%zmetres < 0.0) then
        ! reflect across ground (z=0)
        part%zmetres = -part%zmetres
    endif

    ! cap at entrainment height
    part%zmetres = min(part%zmetres, top_entrainment)

    ! enforce ground
    part%zmetres = max(part%zmetres, 0.0)
    ! --------------------------------------------------
  else
    ! Reflection and position updates
    if (delz.lt.-part%zmetres) then         ! reflection at ground
      part%zmetres = -part%zmetres - delz
    else if (delz.gt.(part%hbl-part%zmetres)) then ! reflection at top
      part%zmetres = -part%zmetres-delz+2.*part%hbl
    else                         ! no reflection
      part%zmetres = part%zmetres+delz
    endif
  endif
end subroutine variable_k_name_within_bl

subroutine variable_k_name_above_bl(part, pextra) 
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  real :: sigu, sigv, sigw ! Turbulent velocity standard deviations
  real :: tlu, tlv, tlw ! Lagrangian timescales
  real :: delz ! Turbulent vertical displacement (m)
  real :: top_entrainment

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

  ! Calculate new horizontal positions
  part%x = part%x + part%turbvelu * tstep * pextra%rmx
  part%y = part%y + part%turbvelv * tstep * pextra%rmy

  ! Calculate new vertical position
  delz=part%turbvelw*tstep 
  if (abs(delz).gt.part%hbl) then
    delz=mod(delz,part%hbl)
  endif

  ! Calculate new vertical position
  part%zmetres = part%zmetres+delz

end subroutine variable_k_name_above_bl

subroutine constant_k_name_within_bl(part, pextra) 
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  integer :: hor_diffu=5300
  real :: top_entrainment=0.03
  real(real64) :: rnd(1)
  real :: mixing_top

  mixing_top = part%hbl + part%hbl*top_entrainment

  call random_number(rnd)
  rnd = rnd - 0.5

  ! Calculate Turbulent Velocities, Ryall and Maryon 1998
  if (nrand+1.gt.max_rands) nrand=1
  part%turbvelu = ((2*hor_diffu)/tstep)**0.5 * rands(nrand)
  part%turbvelv = ((2*hor_diffu)/tstep)**0.5 * rands(nrand+1)
  nrand=nrand+2

  ! Calculate new horizontal positions
  part%x = part%x + part%turbvelu * tstep*pextra%rmx
  part%y = part%y + part%turbvelv * tstep*pextra%rmy

  part%zmetres = (rnd(1) + 0.5) * mixing_top

end subroutine constant_k_name_within_bl

subroutine constant_k_name_above_bl(part, pextra) 
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  real :: delz ! Turbulent vertical displacement (m)
  integer :: hor_diffu=5300/4
  real :: vert_diffu=1.5

  ! Calculate Turbulent Velocities, Ryall and Maryon 1998
  if (nrand+2.gt.max_rands) nrand=1
  part%turbvelu = ((2*hor_diffu)/tstep)**0.5 * rands(nrand)
  part%turbvelv = ((2*hor_diffu)/tstep)**0.5 * rands(nrand+1)
  part%turbvelw = ((2*vert_diffu)/tstep)**0.5 * rands(nrand+2)
  nrand=nrand+3

  ! Calculate new horizontal positions
  part%x = part%x + part%turbvelu * tstep*pextra%rmx
  part%y = part%y + part%turbvelv * tstep*pextra%rmy

  delz = part%turbvelw*tstep

  part%zmetres = part%zmetres+delz
  
end subroutine constant_k_name_above_bl

subroutine diffusion_fields(u_star, w_star, obukhov_length)
  use snapfldML, only: ps2, t2m, hbl2, surface_stress, hflux, tv
  use snapdimML, only: nx, ny
  use, intrinsic :: ieee_arithmetic

  real, intent(out) :: u_star(:, :)
  real, intent(out) :: w_star(:, :)
  real, intent(out) :: obukhov_length(:, :)

  real, parameter :: r=287, g=9.81, k=0.4, cpa=1004.6

  real :: rho_a(nx, ny)

  ! Calculate surface air density
  rho_a = (ps2*100) / (t2m * r)

  ! Calculate friction velocity
  u_star = sqrt(surface_stress/(rho_a))

  ! Calculate the obukhov length. Negative sign removed as ECMWF convention is positive for downward flux
  obukhov_length = rho_a * cpa * t2m * (u_star**3)/(k*g*hflux)

  ! Calculate the convective velocity scale, p.622/118 stull
  ! surface kinematic heat flux = H0/(rho*cpa)
  if (bl_definition == 'constant') then
    hbl2=600
  endif
  w_star = ((g*hbl2*-hflux)/(tv(:,:,2)*rho_a*cpa))**0.333

  where (ieee_is_nan(w_star))
    w_star = 0.0
  end where

end subroutine diffusion_fields

subroutine air_density(rho, rhograd, pressures, tv)
  use snapfldML, only: spec_humid, ps2, hlevel2, t2
  use snapdimML, only: nx, ny, nk
  use snapgrdML, only: alevel, blevel, ivlayer

  real, parameter :: r = 287.0

  real, intent(out) :: rho(:, :, :)
  real, intent(out) :: rhograd(:, :, :)
  real, intent(out) :: pressures(:, :, :)
  real, intent(out) :: tv(:,:,:)

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
        pressures(i,j,k) = alevel(k) * 100 + blevel(k) * ps2(i,j) * 100.0
      end do
      pressures(i,j,2) = ps2(i,j) * 100.0
    end do
  end do

  ! Density 
  do j = 1, ny
    do i = 1, nx
      do k = 2, nk
        rho(i,j,k) = pressures(i,j,k) / (r * tv(i,j,k))
      end do
      ! Fill synthetic surface 
      rho(i,j,1) = rho(i,j,2)
    end do
  end do

  ! Interior points
  do k = 3, nk-1
    do j = 1, ny
      do i = 1, nx
        rhograd(i,j,k) = (rho(i,j,k+1) - rho(i,j,k-1)) / &
                        (hlevel2(i,j,k+1) - hlevel2(i,j,k-1))
      end do
    end do
  end do

  ! Bottom boundary, forward difference
  do j = 1, ny
    do i = 1, nx
      rhograd(i,j,2) = ( rho(i,j,3) - rho(i,j,2) ) &
                       / ( hlevel2(i,j,3) - hlevel2(i,j,2) )
    end do
  end do

  ! Top boundary, backward difference
  do j = 1, ny
    do i = 1, nx
      rhograd(i,j,nk) = ( rho(i,j,nk) - rho(i,j,nk-1) ) &
                        / ( hlevel2(i,j,nk) - hlevel2(i,j,nk-1) )
    end do
  end do

end subroutine

subroutine eta_to_metres(part, pextra)

  use particleML, only: extraParticle, Particle
  use snapfldML,  only: hlevel2, ps2
  use snapgrdML,  only: vlevel, alevel, blevel
  use snapdimML,  only: nk

  implicit none

  type(Particle), intent(inout) :: part
  type(extraParticle), intent(inout) :: pextra

  integer :: i, j, k
  real :: eta
  real :: z1, z2
  real :: p1, p2, px
  real :: frac

  i = part%x
  j = part%y
  eta = part%z

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
  p1 = alevel(k)*100 + blevel(k) * ps2(i,j) * 100.0
  p2 = alevel(k+1)*100 + blevel(k+1) * ps2(i,j) * 100.0

  ! Pressure at particle (linear in frac)
  px = p1 * (1.0 - frac) + p2 * frac

  ! Log-pressure height interpolation
  if (p1 > 0.0 .and. p2 > 0.0) then
    part%zmetres = z1 + (z2 - z1) / log(p2/p1) * log(px/p1)
  else
    part%zmetres = z1 * (1.0 - frac) + z2 * frac
  end if

end subroutine eta_to_metres

subroutine metres_to_eta(part, pextra)

  use particleML, only: extraParticle, Particle
  use snapfldML,  only: hlevel2, ps2
  use snapgrdML,  only: vlevel, alevel, blevel
  use snapdimML,  only: nk

  implicit none

  type(Particle), intent(inout) :: part
  type(extraParticle), intent(inout) :: pextra

  integer :: i, j, k
  real :: z
  real :: z1, z2
  real :: p1, p2, px
  real :: frac

  i = part%x
  j = part%y
  z = part%zmetres

  do k = 2, nk-1
    if (z < hlevel2(i,j,k+1)) exit
  end do
  k = max(1, min(k, nk-1))

  z1 = hlevel2(i,j,k)
  z2 = hlevel2(i,j,k+1)

  p1 = alevel(k)*100 + blevel(k) * ps2(i,j) * 100.0
  p2 = alevel(k+1)*100 + blevel(k+1) * ps2(i,j) * 100.0

  if (p1 > 0.0 .and. p2 > 0.0) then
    px = p1 * exp( log(p2/p1) * (z - z1) / (z2 - z1) )
    frac = (px - p1) / (p2 - p1)
  else
    frac = (z - z1) / (z2 - z1)
  end if

  frac = max(0.0, min(1.0, frac))

  part%z = vlevel(k) * (1.0 - frac) + vlevel(k+1) * frac

end subroutine metres_to_eta

end module rwalkML
