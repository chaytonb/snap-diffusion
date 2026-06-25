module advance_particleML
  implicit none
  private

  public advance_particle_position

  contains

subroutine advance_particle_position(part, pextra, tnow, tstep, rt1, rt2, tf1, tf2, adaptive)

  USE particleML, only: Particle, extraParticle
  USE bldpML, only: set_constant_bl
  USE rwalkML, only: bl_definition

  type(Particle), intent(inout)      :: part
  type(extraParticle), intent(inout) :: pextra
  real, intent(in)                   :: tnow, tstep, rt1, rt2, tf1, tf2 ! Time interpolation variables
  logical, intent(in)                :: adaptive

  if (bl_definition == 'constant') call set_constant_bl(part)

  if (adaptive) then
     call step_adaptive_loop(part, pextra, tnow, tstep, rt1, rt2, tf1, tf2)
  else
     call step_standard_single(part, pextra, tnow, tstep, rt1, rt2, tf1, tf2)
  endif

end subroutine advance_particle_position

subroutine step_adaptive_loop(part, pextra, tnow, tstep, rt1, rt2, tf1, tf2)

  USE particleML, only: extraParticle, Particle
  USE posintML, only: posint_newlevel, posint_vert, calculate_gradient_profiles
  USE snapfldML, only: hlevel2
  USE rwalkML, only: diffusion_scheme, turbulence_master, well_mixed_test
  USE forwrdML, only: forwrd
  USE snapdimML, only: nk

  type(Particle), intent(inout)  :: part
  type(extraParticle), intent(inout) :: pextra
  real, intent(in) :: tnow, tstep, rt1, rt2, tf1, tf2

  real :: ux, vx ! Accumumlated displacements in horizontal direction
  logical :: interpol_exists(nk) ! Tracking which vertical levels have been interpolated to already
  real:: uprof(nk), vprof(nk), wprof(nk), rhoprof(nk), rhogradprof(nk)
  integer :: k0 
  real :: t_local, dt_remaining
  integer :: i, j

  ux = 0.0
  vx = 0.0
  t_local = 0.0
  interpol_exists = .false.

  ! Calculate advective velocities
  if (.not. well_mixed_test) call forwrd(tf1, tf2, tnow, tstep, part, pextra)

  if (diffusion_scheme == 'TKE') call calculate_gradient_profiles(part, rt1, rt2)

  ! Apply adaptive timesteps in the ABL
  if (part%zmetres < part%hbl .OR. diffusion_scheme == 'TKE' .OR. diffusion_scheme == 'random_walk_name') then
    do while (t_local < tstep)

      i = part%x
      j = part%y
      ! Identify bracketing vertical model levels
      do k0 = 1, nk-1
        if (part%zmetres <= hlevel2(i, j, k0+1)) exit
      end do
      k0 = max(1, min(k0, nk-1))

      dt_remaining = tstep - t_local

      ! Horizontally interpolate to particle position on bracketing levels of particle
      if (.not.interpol_exists(k0)) then
        call posint_newlevel(part, pextra, k0, uprof, vprof, wprof, rhoprof, rhogradprof, rt1, rt2)
        interpol_exists(k0) = .true.
      endif  
      if (.not.interpol_exists(k0+1)) then
        call posint_newlevel(part, pextra, k0+1, uprof, vprof, wprof, rhoprof, rhogradprof, rt1, rt2)
        interpol_exists(k0+1) = .true.
      endif

      ! Interpolate vertically to particle position from bracketing levels
      call posint_vert(part, pextra, uprof, vprof, wprof, rhoprof, rhogradprof, k0)

      ! Calculate turbulent velocities and apply vertical turbulent displacement
      call turbulence_master(part, pextra, dt_remaining, .TRUE.)

      ! Accumulate horizontal velocities
      if (.not. well_mixed_test) then
        ux = ux + part%ptstep * (pextra%u + part%turbvelu)
        vx = vx + part%ptstep * (pextra%v + part%turbvelv)
        ! Apply vertical advection
        part%zmetres = max(part%zmetres + pextra%w * part%ptstep, 0.5)
      endif

      t_local = t_local + part%ptstep

      ! ABL top exit condition, complete step above ABL
      if (part%zmetres > part%hbl .and. diffusion_scheme /= 'TKE' .and. diffusion_scheme /= 'random_walk_name') exit

    end do
  endif

  ! Catch for residual time if the particle left the BL early or if particle already above ABL
  if (t_local < tstep) then

     dt_remaining = tstep - t_local

     if (.not. well_mixed_test) call forwrd(tf1, tf2, tnow+t_local, part%ptstep, part, pextra)

     call turbulence_master(part, pextra, dt_remaining, .FALSE.)
     
     if (.not. well_mixed_test) then
        ux = ux + part%ptstep * (pextra%u + part%turbvelu)
        vx = vx + part%ptstep * (pextra%v + part%turbvelv)
        part%zmetres = part%zmetres + pextra%w * part%ptstep
     endif
  endif

  ! Apply accumulated horizontal displacement
  part%x = part%x + ux * pextra%rmx
  part%y = part%y + vx * pextra%rmy

end subroutine step_adaptive_loop


subroutine step_standard_single(part, pextra, tnow, tstep, rt1, rt2, tf1, tf2)

  USE forwrdML, only: forwrd
  USE rwalkML, only: turbulence_master, well_mixed_test, diffusion_scheme, diffusion_in_metres
  USE particleML, only: extraParticle, Particle
  USE posintML, only: vert_interpol_rho_only

  type(Particle), intent(inout) :: part
  type(extraParticle), intent(inout) :: pextra
  real, intent(in) :: tnow, tstep, rt1, rt2, tf1, tf2
  real :: dt_remaining

  dt_remaining = tstep
  part%ptstep = tstep

  ! Calculate advective velocities
  if (.not. well_mixed_test) call forwrd(tf1, tf2, tnow, tstep, part, pextra)

  ! Interpolated density values required for Langevin scheme
  if (diffusion_scheme == 'random_walk_flexpart') call vert_interpol_rho_only(part, pextra, rt1, rt2)

  ! Calculate turbulent velocities and apply vertical turbulent displacement
  call turbulence_master(part, pextra, dt_remaining, .FALSE.)

  ! Apply horizontal advection and diffusion
  if (.not.well_mixed_test .and. diffusion_in_metres) then
    part%x = part%x + (part%turbvelu + pextra%u) * part%ptstep * pextra%rmx
    part%y = part%y + (part%turbvelv + pextra%v) * part%ptstep * pextra%rmy
    ! Apply vertical advection
    part%zmetres = part%zmetres + pextra%w * part%ptstep
  endif

end subroutine step_standard_single

end module advance_particleML