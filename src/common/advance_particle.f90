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

  ! PULL INTERPOLATION: Only calculate tbl if the scheme/definition requires it
  if (bl_definition == 'constant') call set_constant_bl(part)

  if (adaptive) then
     call step_adaptive_loop(part, pextra, tnow, tstep, rt1, rt2, tf1, tf2)
  else
     call step_standard_single(part, pextra, tnow, tstep, tf1, tf2)
  endif

end subroutine advance_particle_position


subroutine step_adaptive_loop(part, pextra, tnow, tstep, rt1, rt2, tf1, tf2)

  USE particleML, only: extraParticle, Particle
  USE posintML, only: posint_newlevel, posint_vert
  USE snapfldML, only: hlevel2
  USE rwalkML, only: diffusion_scheme, turbulence_master, well_mixed_test
  USE forwrdML, only: forwrd
  USE snapdimML, only: nk

  type(Particle), intent(inout)  :: part
  type(extraParticle), intent(inout) :: pextra
  real, intent(in)                   :: tnow, tstep, rt1, rt2, tf1, tf2

  real :: dt_min = 1 ! minimum timestep size for adaptive timestepping
  real :: ux, vx ! Accumumlated displacements in horizontal direction
  logical, allocatable :: interpol_exists(:)
  real, allocatable :: uprof(:), vprof(:), wprof(:), rhoprof(:), rhogradprof(:)
  integer :: k0 
  real :: t_local
  integer :: i, j

  allocate(interpol_exists(nk), uprof(nk), vprof(nk), wprof(nk), rhoprof(nk), rhogradprof(nk))

  ux = 0.0
  vx = 0.0
  t_local = 0.0
  interpol_exists = .false.

  if (.not. well_mixed_test) then
    call forwrd(tf1, tf2, tnow, tstep, part, pextra)
  endif

  do while (t_local < tstep)

    i = part%x
    j = part%y
    ! Identify bracketing vertical model levels
    do k0 = 1, nk-1
      if (part%zmetres <= hlevel2(i, j, k0+1)) exit
    end do
    k0 = max(1, min(k0, nk-1))

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

    ! Calculate physics scales and timesteps
    call turbulence_master(part, pextra)

    ! Enforce time step constrains 
    part%ptstep = min(part%tlw, &
                            part%hbl / max(2.0 * abs(part%turbvelw), 1.e-4), &
                            0.5 / max(abs(part%dsigwdz), 1.e-6)) * 0.1

    part%ptstep = max(part%ptstep, dt_min) ! enforce minimum size for steps

    ! Ensure it does not exceed the global timestep window
    part%ptstep = min(part%ptstep, (tstep - t_local))

    ! Accumulate horizontal velocities
    if (.not. well_mixed_test) then
        ux = ux + part%ptstep * (pextra%u + part%turbvelu)
        vx = vx + part%ptstep * (pextra%v + part%turbvelv)
    endif

    ! Advance vertically
    part%zmetres = max(part%zmetres + pextra%w * part%ptstep, 0.5)
    t_local = t_local + part%ptstep

    ! ABL top exit condition
    if (part%zmetres > part%hbl .and. &
        diffusion_scheme /= 'random_walk_name' .and. &
        diffusion_scheme /= 'TKE') exit

  end do

  ! Catch for residual time if the particle left the BL early
  if (t_local < tstep) then

     part%ptstep = tstep - t_local

     if (.not. well_mixed_test) call forwrd(tf1, tf2, tnow+t_local, part%ptstep, part, pextra)

     call turbulence_master(part, pextra)
     
     if (.not. well_mixed_test) then
        ux = ux + part%ptstep * (pextra%u + part%turbvelu)
        vx = vx + part%ptstep * (pextra%v + part%turbvelv)
     endif
     part%zmetres = part%zmetres + pextra%w * part%ptstep
  endif

  ! Apply final horizontal displacement
  part%x = part%x + ux * pextra%rmx
  part%y = part%y + vx * pextra%rmy

end subroutine step_adaptive_loop


subroutine step_standard_single(part, pextra, tnow, tstep, tf1, tf2)

  USE forwrdML, only: forwrd
  USE rwalkML, only: turbulence_master, well_mixed_test
  USE particleML, only: extraParticle, Particle

  type(Particle), intent(inout) :: part
  type(extraParticle), intent(inout) :: pextra
  real, intent(in) :: tnow, tstep, tf1, tf2

  part%ptstep = tstep

  if (.not. well_mixed_test) then
     call forwrd(tf1, tf2, tnow, tstep, part, pextra)
  endif

  call turbulence_master(part, pextra)

  part%x = part%x + (part%turbvelu + pextra%u) * tstep * pextra%rmx
  part%y = part%y + (part%turbvelv + pextra%v) * tstep * pextra%rmy

  part%zmetres = part%zmetres + pextra%w * tstep

end subroutine step_standard_single

end module advance_particleML