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

  ! Values for random number generation
  integer(int32), parameter :: max_rands=6000000
  integer(int32) :: nrand=1
  real :: rands(max_rands) ! initialise array to store random numbers

  real(real64), save, public :: a_in_bl = 0.5
  real(real64), save, public :: a_above_bl = 0.25
  real(real64), save, public :: b = 0.875
  logical, save, public :: turb_homogeneous = .TRUE.
  logical, save, public :: well_mixed_test = .FALSE.
  character(len=64), save, public :: diffusion_scheme = ''
  character(len=64), save, public :: bl_definition = ''
  character(len=64), save, public :: meteo_type = ''
  character(len=64), save, public :: entrainment_scheme = ''

  public rwalk_init, diffusion_fields, air_density, turbulence_master, eta_to_metres, interp_tke_to_hybrid_field

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

  ! If diffusion scheme not default SNAP, then generate normally distributed random nums
  if (diffusion_scheme /= '') then
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
  elseif (diffusion_scheme == 'TKE') then
    call eta_to_metres(part, pextra)
    call tke_diffusion(part, pextra)
    call metres_to_eta(part, pextra)
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

subroutine flexpart_diffusion_within_abl(part, pextra) ! Based on Hanna1 from FLEXPART code, turbswitch FALSE
  USE particleML, only: extraParticle, Particle
  use snapfldML, only: rho, rhograd, ps2, t2m, hbl2, surface_stress, hflux, tv
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
  real :: turb_delu, turb_delv ! turbulent displacements in u and v directions

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

  if (well_mixed_test) then
    wst = 1.5
    ol = -100
    pextra%ust = 0.5
  endif

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

  if (.not.well_mixed_test) then
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

    turb_delu = part%turbvelu * part%ptstep
    turb_delv = part%turbvelv * part%ptstep

    ! Transform turbulent displacements from along/crosswind to x/y
    call align_turbvels(part, pextra, turb_delu, turb_delv)

    ! Calculate new horizontal positions
    part%x = part%x + turb_delu*pextra%rmx
    part%y = part%y + turb_delv*pextra%rmy
  endif

  ! Vertical interpolation of rho and rhograd
  call interp_rho(part, pextra)
  rhoaux=pextra%rhograd/pextra%rho

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

  if (.not.well_mixed_test) then
    part%x = part%x + part%turbvelu * part%ptstep*pextra%rmx
    part%y = part%y + part%turbvelv * part%ptstep*pextra%rmy
  endif

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

  if (well_mixed_test) then
    wst = 1.5
    pextra%ol = -100
    ust = 0.5
  endif

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

    dsigwdz = (1.0/(2.0 * sigw * part%hbl)) *(wst**2*(0.8*scaled_height**(-1.0/3.0) - 1.8*scaled_height**(2.0/3.0)) - 1.4*ust**2)
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

  if (.not.well_mixed_test) then
    part%x = part%x + part%turbvelu * part%ptstep*pextra%rmx
    part%y = part%y + part%turbvelv * part%ptstep*pextra%rmy
  endif

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
  if (.not.well_mixed_test) then
    part%x = part%x + part%turbvelu * part%ptstep*pextra%rmx
    part%y = part%y + part%turbvelv * part%ptstep*pextra%rmy
  endif

  if (nrand+1.gt.max_rands) nrand=1
  ! Calculate turbulent vertical velocity
  part%turbvelw=part%turbvelw*(1-(dttlw)) + (2*sigw**2 * dttlw)**0.5*rands(nrand)
  nrand=nrand+1

  ! Calculate new vertical position
  delz=part%turbvelw*part%ptstep  

  !Reflect at inversion (from above) so it is symmetric with BL routine
  if (part%zmetres+delz .lt. part%hbl) then
    part%zmetres = part%hbl + (part%hbl - part%zmetres - delz)
    part%turbvelw = -part%turbvelw
  else
    part%zmetres = part%zmetres + delz
  endif

end subroutine name_random_walk_profile_above_bl

subroutine variable_k_name_within_bl(part, pextra) 
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  real :: k ! Von Karmans constant
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

  if (well_mixed_test) then
    wst = 1.5
    pextra%ol = -100
    ust = 0.5
  endif
  
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

    eps = (1.5 - 1.2 * (scaled_height)**0.333)*(wst**3/part%hbl) + (ust**3 * (1-0.8*scaled_height)/(k*(part%hbl*0.5)))

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

  if (.not.well_mixed_test) then
    ! Calculate new horizontal positions
    part%x = part%x + part%turbvelu * tstep*pextra%rmx
    part%y = part%y + part%turbvelv * tstep*pextra%rmy
  endif

  delz=part%turbvelw*tstep 

  ! Calculate new vertical position
  if (abs(delz).gt.part%hbl) then
    delz=mod(delz,part%hbl)
  endif

  ! Reflection and position updates
  if (delz.lt.-part%zmetres) then         ! reflection at ground
    part%zmetres = -part%zmetres - delz
  else if (delz.gt.(part%hbl-part%zmetres)) then ! reflection at top
    part%zmetres = -part%zmetres-delz+2.*part%hbl
  else                         ! no reflection
    part%zmetres = part%zmetres+delz
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
  real :: znew, ztop

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

  if (.not.well_mixed_test) then
    ! Calculate new horizontal positions
    part%x = part%x + part%turbvelu * tstep * pextra%rmx
    part%y = part%y + part%turbvelv * tstep * pextra%rmy
  endif

  ! Calculate new vertical position
  delz=part%turbvelw*tstep 

  ! Reflect at inversion (from above) so it is symmetric with BL routine
  if (part%zmetres+delz .lt. part%hbl) then
    part%zmetres = 2*part%hbl - (part%zmetres + delz)
  else
    part%zmetres = part%zmetres + delz
  endif

end subroutine variable_k_name_above_bl

subroutine constant_k_name_within_bl(part, pextra) 
  USE particleML, only: extraParticle, Particle
  
  !> particle with information
  type(Particle), intent(inout)  :: part
  !> extra information regarding the particle
  type(extraParticle), intent(inout) :: pextra

  ! Locals
  real :: rnd(1)
  real :: hbl, sigma_ft, delta_ext, mixing_top

  integer, parameter :: hor_diffu = 5300 ! m^2 s^-1 (BL horizontal diffusion)
  real, parameter :: K_ft = 1.5 ! m^2 s^-1 (FT vertical diffusion)
  real, parameter :: alpha = 1.7 ! entrainment zone extension factor

  ! Compute free troposphere displacement length scale
  sigma_ft = sqrt(2.0 * K_ft * tstep) 
  delta_ext = alpha * sigma_ft

  mixing_top = part%hbl + delta_ext

  call random_number(rnd)

  if (.not.well_mixed_test) then
    ! Calculate Turbulent Velocities, Ryall and Maryon 1998
    if (nrand+1.gt.max_rands) nrand=1
    part%turbvelu = ((2*hor_diffu)/tstep)**0.5 * rands(nrand)
    part%turbvelv = ((2*hor_diffu)/tstep)**0.5 * rands(nrand+1)
    nrand=nrand+2

    ! Calculate new horizontal positions
    part%x = part%x + part%turbvelu * tstep*pextra%rmx
    part%y = part%y + part%turbvelv * tstep*pextra%rmy
  endif

  part%zmetres = rnd(1) * mixing_top

end subroutine constant_k_name_within_bl

subroutine constant_k_name_above_bl(part, pextra)
  use particleML, only: extraParticle, Particle
  implicit none

  type(Particle),      intent(inout) :: part
  type(extraParticle), intent(inout) :: pextra

  real, parameter :: hor_diffu_ft = 5300.0/4.0 ! m^2 s^-1
  real, parameter :: K_ft = 1.5 ! m^2 s^-1
  real, parameter :: alpha = 1.7 ! entrainment zone extension factor
  real :: dt, delz, z0, z1, h, a, b, p_cross
  real :: rnd(2)
  real :: hbl, sigma_ft, delta_ext, mixing_top

  dt = tstep
  h = part%hbl ! BL height at the particle location/time

  ! Turbulent velocities
  if (nrand + 2 .gt. max_rands) nrand = 1
  part%turbvelu = sqrt((2.0 * hor_diffu_ft) / dt) * rands(nrand)
  part%turbvelv = sqrt((2.0 * hor_diffu_ft) / dt) * rands(nrand + 1)
  part%turbvelw = sqrt((2.0 * K_ft) / dt) * rands(nrand + 2)
  nrand = nrand + 3

  if (.not. well_mixed_test) then
    part%x = part%x + part%turbvelu * dt * pextra%rmx
    part%y = part%y + part%turbvelv * dt * pextra%rmy
  endif

  ! Vertical step
  z0 = part%zmetres
  delz = part%turbvelw * dt
  z1 = z0 + delz

  ! Compute free troposphere displacement length scale
  sigma_ft = sqrt(2.0 * K_ft * tstep) 
  delta_ext = alpha * sigma_ft

  mixing_top = part%hbl + delta_ext

  if (z1 <= mixing_top) then
    ! End-point crosses into BL, then instant-mix within BL
    call random_number(rnd)
    part%zmetres = rnd(1) * mixing_top
    return
  end if

  ! No crossing, then accept the FT step
  part%zmetres = z1

end subroutine constant_k_name_above_bl

subroutine tke_diffusion(part, pextra)
  USE particleML, only: extraParticle, Particle
  USE snapfldML, only: tke_hyb, rho, rhograd, t2, pressures, hlevel2
  USE snapgrdML, only: ivlayer
  USE snapdimML, only: nk
  use, intrinsic :: ieee_arithmetic

  type(Particle), intent(inout)  :: part
  type(extraParticle), intent(inout) :: pextra

  integer :: i, j, k, kp, part_vert_index
  real :: sigu, sigv, sigw
  real :: tke_z, yl, yl_up, yl_down, sum, e1, part_z
  real :: tlu, tlv, tlw
  real :: ru, rv, rw
  real :: delz, dt, rhoaux
  real :: dttlw
  real :: pttprof(nk), pttrefprof(nk)
  integer :: indz, indzp
  real :: r, cp, g

  ! Physical constants
  r = 287.0
  cp = 1004.0
  g = 9.81

  i = part%x
  j = part%y
  part_vert_index = int(part%z * 10000)
  k = ivlayer(part_vert_index)
  dt = part%ptstep
  part_z = part%zmetres

  call calc_turb_params_tke(i, j, k, sigu, sigv, sigw, tlu, tlv, tlw)

  part%tlw = tlw

  ! Interpolate air density and gradient to part vert position
  call interp_rho(part, pextra)
  rhoaux = pextra%rhograd/pextra%rho

  dttlw = dt/tlw

  if (nrand+2.gt.max_rands) nrand=1

  ! Vertical turbulence
  if (dttlw.lt..5) then
      part%turbvelw = (1.-dttlw)*part%turbvelw+rands(nrand)*sqrt(2.*dttlw)+dt*rhoaux*sigw
  else
      rw = exp(-dttlw)
      part%turbvelw = rw*part%turbvelw+rands(nrand)*sqrt(1.-rw**2)+tlw*(1.-rw)*rhoaux*sigw
  end if

  delz = part%turbvelw * sigw * dt

  call vertical_reflection_step(i, j, k, part_z, tlu, tlv, tlw, sigu, sigv, sigw, part%turbvelw, delz)

  part%zmetres = part_z

  ! Horizontal turbulence
  ru = exp(-dt / tlu)
  rv = exp(-dt / tlv)
  part%turbvelu = ru * part%turbvelu + rands(nrand+1) * sqrt(1.0 - ru**2) * sigu
  part%turbvelv = rv * part%turbvelv + rands(nrand+2) * sqrt(1.0 - rv**2) * sigv

  if (.not.well_mixed_test) then
    part%x = part%x + part%turbvelu * dt * pextra%rmx
    part%y = part%y + part%turbvelv * dt * pextra%rmy
  endif

  nrand = nrand + 3

end subroutine tke_diffusion

subroutine vertical_reflection_step(i, j, k, part_z, tlu, tlv, tlw, sigu, sigv, sigw, part_turbvelw, delz)

  USE snapfldML, only: hlevel2, hlayer2, hinterf, tke_hyb
  USE snapgrdML, only: ivlayer
  USE snapdimML, only: nk
  use, intrinsic :: ieee_arithmetic

  ! Arguments:
  integer, intent(in) :: i, j
  integer, intent(inout) :: k
  real, intent(inout) :: tlu, tlv, sigu, sigv, sigw, tlw
  real(4), intent(inout) :: part_z, delz
  real(8), intent(inout) :: part_turbvelw

  integer :: indz_c, indzp_c, dir, ind, part_vert_index, k_c, interface_ind, num_loops=1
  real :: z_bot, z_top, z_bot_next, z_top_next
  real :: ts, ratio, sigw_c, tlw_c
  real :: rnd(1)
  logical :: reflect
  real :: zt_tmp

  do
    call random_number(rnd)
    ! Calculate layer boundaries
    z_bot = hinterf(i,j,k)
    z_top = hinterf(i,j,k+1)

    ! Temp fix, ivlevel (eta) indexing doesn't exactly line up with hlevel (m)
    part_z = max(z_bot, part_z)
    part_z = min(z_top, part_z)

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
  USE snapfldML, only: tke_hyb, rho, rhograd, t2, pressures, hlevel2
  USE snapgrdML, only: ivlayer
  USE snapdimML, only: nk

  integer, intent(in) :: i, j, k
  real, intent(inout) :: tlu, tlv, tlw, sigu, sigv, sigw

  integer :: kp, part_vert_index
  real :: tke_z, yl, yl_up, yl_down, sum, e1
  real :: ru, rv, rw
  real :: delz, dt, rhoaux
  real :: pttprof(nk), pttrefprof(nk)
  integer :: indz, indzp
  real :: r, cp, g

  do kp = 1, nk
     pttprof(kp) = t2(i, j, kp)
     ! Reference profile: use local value, should try to use domain average?
     pttrefprof(kp) = pttprof(kp)
  end do

  indz = k
  indzp = min(k+1, nk)

  tke_z = tke_hyb(i, j, k)

  ! Compute BL89 mixing length
  e1 = -g / pttrefprof(indz) * (pttprof(indz) - pttprof(indzp)) * (hlevel2(i, j, indzp) - hlevel2(i, j, indz))
  if (e1 >= tke_z) then
    yl = hlevel2(i, j, indzp) - hlevel2(i, j, indz)
  else
    ! Upward
    sum = 0.0
    kp = indzp
    do while (kp < nk)
      sum = sum - (g / pttrefprof(kp) * (pttprof(indz) - pttprof(kp)) * (hlevel2(i, j, kp) - hlevel2(i, j, kp-1)))
      if (sum >= tke_z) exit
      kp = kp + 1
    end do
    yl_up = hlevel2(i, j, kp) - hlevel2(i, j, indz)
    ! Downward
    sum = 0.0
    kp = indz
    do while (kp > 2)
      sum = sum - (g / pttrefprof(kp) * (pttprof(kp) - pttprof(indzp)) * (hlevel2(i, j, kp+1) - hlevel2(i, j, kp)))
      if (sum >= tke_z) exit
      kp = kp - 1
    end do
    yl_down = hlevel2(i, j, indzp) - hlevel2(i, j, kp)
    yl = ((yl_up**(-2.0/3.0) + yl_down**(-2.0/3.0))/2.0)**(-3.0/2.0)
  end if

  ! 3D Turbulent velocity scales (assuming isotropic for now)
  sigu = sqrt(2.0 * tke_z / 3.0)
  sigv = sqrt(2.0 * tke_z / 3.0)
  sigw = sqrt(2.0 * tke_z / 3.0)

  ! Timescales
  tlu = 2.0 * yl / max(sigu, 1e-6)
  tlv = 2.0 * yl / max(sigv, 1e-6)
  tlw = 2.0 * yl / max(sigw, 1e-6)
  tlu = max(10.0, tlu)
  tlv = max(10.0, tlv)
  tlw = max(30.0, tlw)

end subroutine calc_turb_params_tke

subroutine diffusion_fields
  use snapfldML, only: ps2, t2m, hbl2, surface_stress, hflux, tv, obukhov_l_io, u_star_io, w_star_io
  use snapdimML, only: nx, ny
  use, intrinsic :: ieee_arithmetic

  real, parameter :: r=287, g=9.81, k=0.4, cpa=1004.6

  real :: rho_a(nx, ny)
  real :: fhsfc(nx, ny) ! Surface kinematic heat flux 

  ! Calculate surface air density
  rho_a = (ps2*100) / (t2m * r)

  ! Calculate friction velocity
  u_star_io = sqrt(surface_stress/(rho_a))

  ! surface kinematic heat flux = H0/(rho*cpa)
  fhsfc = hflux/(rho_a*cpa)

  ! Calculate the obukhov length
  obukhov_l_io = - (tv(:,:,2) * u_star_io**3)/(k*g*fhsfc)

  if (bl_definition == 'constant') then
    hbl2=600
  endif

  ! Calculate the convective velocity scale
  w_star_io = ((g/tv(:,:,2))*hbl2*fhsfc)**0.333

  where (ieee_is_nan(w_star_io))
    w_star_io = 0.0
  end where
  
end subroutine diffusion_fields

subroutine air_density
  use snapfldML, only: spec_humid, ps2, hlevel2, t2, rho, rhograd, pressures, tv
  use snapdimML, only: nx, ny, nk
  use snapgrdML, only: alevel, blevel, ivlayer

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

subroutine interp_rho(part, pextra)
  use particleML, only: Particle, extraParticle
  use snapfldML,  only: hlevel2, ps2, rho, rhograd
  use snapgrdML,  only: vlevel, alevel, blevel
  use snapdimML,  only: nk
  implicit none

  type(Particle),      intent(inout) :: part
  type(extraParticle), intent(inout) :: pextra

  integer :: i, j, k
  real :: eta
  real :: z1, z2
  real :: p1, p2, px
  real :: frac, w, denom
  real :: y1, y2, m1, m2
  real, parameter :: tiny_denom = 1.0e-12

  ! Grid indices and particle eta
  i = part%x
  j = part%y
  eta = part%z

  ! 1) Find bracketing layer in eta
  do k = 1, nk-1
    if (eta > vlevel(k+1)) exit
  end do
  k = max(1, min(k, nk-1))

  ! 2) Fraction in eta
  frac = (eta - vlevel(k)) / (vlevel(k+1) - vlevel(k))
  frac = max(0.0, min(1.0, frac))

  ! 3) Heights at bracketing levels
  z1 = hlevel2(i,j,k)
  z2 = hlevel2(i,j,k+1)

  ! 4) Pressures at bracketing levels (Pa)
  p1 = alevel(k)*100.0 + blevel(k) * ps2(i,j) * 100.0
  p2 = alevel(k + 1)*100.0 + blevel(k + 1) * ps2(i,j) * 100.0

  ! 5) Pressure at particle (linear in frac)
  px = p1 * (1.0 - frac) + p2 * frac

  ! 6) Particle geometric height using log-pressure interpolation 
  if (p1 > 0.0 .and. p2 > 0.0) then
    part%zmetres = z1 + (z2 - z1) / log(p2/p1) * log(px/p1)
  else
    part%zmetres = z1 * (1.0 - frac) + z2 * frac
  end if

  ! 7) Values at bracketing levels
  y1 = rho(i,j,k)
  y2 = rho(i,j,k+1)
  m1 = rhograd(i,j,k)     ! assumed dρ/dz at level k
  m2 = rhograd(i,j,k+1)   ! assumed dρ/dz at level k+1

  ! 8) Log-pressure interpolation weight
  if (p1 > 0.0 .and. p2 > 0.0) then
    denom = log(p2/p1)
    if (abs(denom) > tiny_denom .and. px > 0.0) then
      w = log(px/p1) / denom
      w = max(0.0, min(1.0, w))
    else
      ! Fallback if pressures are nearly identical
      w = frac
    end if
  else
    ! Fallback if invalid pressures
    w = frac
  end if

  ! 9) Interpolate rho and rhograd using log-pressure weight
  pextra%rho = y1*(1.0 - w) + y2*w
  pextra%rhograd = m1*(1.0 - w) + m2*w

end subroutine interp_rho

subroutine align_turbvels(part, pextra, turb_delu, turb_delv)
  use particleML, only: extraParticle, Particle
  implicit none

  type(Particle),      intent(in)    :: part
  type(extraParticle), intent(in)    :: pextra
  real,                intent(inout) :: turb_delu, turb_delv

  real :: umean, vmean, mag_inv, cosphi, sinphi
  real :: du_old, dv_old
  real, parameter :: eps = 1.e-30

  umean = pextra%u
  vmean = pextra%v
  mag_inv = 1.0 / max(sqrt(umean*umean + vmean*vmean), eps)

  cosphi = umean * mag_inv
  sinphi = vmean * mag_inv

  du_old = turb_delu
  dv_old = turb_delv

  turb_delu = du_old * cosphi - dv_old * sinphi
  turb_delv = du_old * sinphi + dv_old * cosphi
end subroutine align_turbvels


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

subroutine interp_tke_profile_to_hybrid(p_tke, logp_tke, tke_prof, ntke, &
                                        p_prof, tke_hyb_prof, nlev)
  implicit none
  integer, intent(in) :: ntke, nlev
  real(kind=8), intent(in)  :: p_tke(ntke), logp_tke(ntke)
  real(kind=8), intent(in)  :: tke_prof(ntke)       ! TKE on fixed p-levels
  real(kind=8), intent(in)  :: p_prof(nlev)         ! pressures at model levels
  real(kind=8), intent(out) :: tke_hyb_prof(nlev)   ! TKE at model levels

  integer :: k, k1, k2, klo, khi, kmid
  real(kind=8) :: logp_h, w
  real(kind=8) :: pmin, pmax

  pmin = p_tke(1)
  pmax = p_tke(ntke)

  do k = 1, nlev
     logp_h = log(p_prof(k))

     ! Clip outside TKE pressure range
     if (p_prof(k) <= pmin) then
        tke_hyb_prof(k) = tke_prof(1)
     else if (p_prof(k) >= pmax) then
        tke_hyb_prof(k) = tke_prof(ntke)
     else
        ! Find bracketing TKE levels in log-pressure space
        ! klo = 1
        ! khi = ntke
        do k1 = 1, ntke-1
          if (logp_tke(k1) <= logp_h .and. logp_tke(k1+1) >= logp_h) then
              k2 = k1 + 1
              exit
          end if
        end do

        ! Linear interpolation in log-pressure
        w = (logp_h - logp_tke(k1)) / (logp_tke(k2) - logp_tke(k1))
        tke_hyb_prof(k) = (1.d0 - w) * tke_prof(k1) + w * tke_prof(k2)
     end if
  end do

end subroutine interp_tke_profile_to_hybrid

subroutine interp_tke_to_hybrid_field
  use snapdimML, only: nx, ny, nk
  use snapfldML, only: tke, pressures, tke_hyb
  implicit none

  integer, parameter :: ntke = 15

  ! TKE fixed pressure levels in hPa
  real(kind=8), dimension(ntke), parameter :: p_tke_hpa = &
       (/ 200.d0, 250.d0, 500.d0, 600.d0, &
          700.d0, 750.d0, 800.d0, 825.d0, 850.d0, 875.d0, &
          900.d0, 925.d0, 950.d0, 975.d0, 1000.d0 /)

  real(kind=8), dimension(ntke), parameter :: p_tke_pa  = p_tke_hpa * 100.d0
  real(kind=8), dimension(ntke), parameter :: logp_tke  = log(p_tke_pa)

  integer :: i, j, k
  real(kind=8), dimension(ntke) :: tke_prof
  real(kind=8), dimension(nk)   :: p_prof, tke_hyb_prof

  do j = 1, ny
    do i = 1, nx

      ! Extract TKE profile at this (i,j) over fixed pressure levels
      do k = 1, ntke
        tke_prof(k) = tke(i,j,k)
      end do

      ! Strange values of TKE at TOA, set to zero
      tke_prof(1) = 0

      ! Extract pressure profile at this (i,j) on model levels
      do k = 1, nk
        p_prof(k) = pressures(i,j,k)   ! Pa
      end do

      ! Interpolate this column
      call interp_tke_profile_to_hybrid(p_tke_pa, logp_tke, tke_prof, ntke, &
                                        p_prof, tke_hyb_prof, nk)

      ! Store result
      do k = 1, nk
        tke_hyb(i,j,k) = tke_hyb_prof(k)
      end do

    end do
  end do

end subroutine interp_tke_to_hybrid_field

end module rwalkML
