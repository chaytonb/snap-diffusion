! SNAP: Servere Nuclear Accident Programme
! Copyright (C) 1992-2026   Norwegian Meteorological Institute
! License: GNU GPL v3 or later

module posintML
  implicit none
  private

  public :: posint, posint_vert, posint_newlevel

  contains

!> Interpolation in 2D
!> private, inlined only for better readability
  pure real function interp(a00, a10, a01, a11, c1, c2, c3, c4)
    real, intent(in) :: a00, a10, a01, a11, c1, c2, c3, c4
    interp = c1*a00 + c2*a10 + c3*a01 + c4*a11
  end function interp


!> Purpose:  Interpolation of boundary layer top and height
!>           and precipitation intensity to particle positions
subroutine posint(part,rt1,rt2,pextra)
  USE particleML, only: Particle, extraParticle
  USE snapgrdML, only: gparam
  USE rwalkML, only: turbulence_fields_required
  USE snapfldML, only: xm, ym, bl1, bl2, hbl1, hbl2, precip, &
    u_star1, u_star2, obukhov_l1, obukhov_l2, w_star1, w_star2

!> particle
!> (with all particles at the same horizontal position)
  type(Particle), intent(inout) :: part
!> fractions from last timestep
  real, intent(in) ::    rt1
!> fractions to next timestep
  real, intent(in) :: rt2
!> extra information interpolated to the particle
!> position, mainly rmx/rmy/prc
  type(extraParticle), intent(out) :: pextra

  integer :: i,j
  real :: dxgrid,dygrid,dx,dy,c1,c2,c3,c4,bl,hbl,rmx,rmy,ol,ust,wst
  real :: pr

  if (.not.part%is_active()) then
    ! No need to compute any properties on the particle
    return
  endif

  dxgrid=gparam(7)
  dygrid=gparam(8)

  !..for horizontal interpolations
  ! i,j = center of the (lower left) grid cell in which the particle is located
  ! dx,dy = distance of the particle to the center of the grid cell in grid units
  i=int(part%x)
  j=int(part%y)
  dx=part%x-i
  dy=part%y-j
  c1=(1.-dy)*(1.-dx)
  c2=(1.-dy)*dx
  c3=dy*(1.-dx)
  c4=dy*dx

 !..interpolation

  !..top of boundary layer
  bl= rt1*interp(bl1(i,j), bl1(i+1,j), bl1(i,j+1), bl1(i+1,j+1), c1, c2, c3, c4) &
      +rt2*interp(bl2(i,j), bl2(i+1,j), bl2(i,j+1), bl2(i+1,j+1), c1, c2, c3, c4)
  !..height of boundary layer
  hbl= rt1*interp(hbl1(i,j), hbl1(i+1,j), hbl1(i,j+1), hbl1(i+1,j+1), c1, c2, c3, c4) &
      +rt2*interp(hbl2(i,j), hbl2(i+1,j), hbl2(i,j+1), hbl2(i+1,j+1), c1, c2, c3, c4)


  if (turbulence_fields_required) then
    !..friction velocity
    ust= rt1*interp(u_star1(i,j), u_star1(i+1,j), u_star1(i,j+1), u_star1(i+1,j+1), c1, c2, c3, c4) &
      +rt2*interp(u_star2(i,j), u_star2(i+1,j), u_star2(i,j+1), u_star2(i+1,j+1), c1, c2, c3, c4)

    !..obukhov length
    ol= rt1*interp(obukhov_l1(i,j), obukhov_l1(i+1,j), obukhov_l1(i,j+1), obukhov_l1(i+1,j+1), c1, c2, c3, c4) &
      +rt2*interp(obukhov_l2(i,j), obukhov_l2(i+1,j), obukhov_l2(i,j+1), obukhov_l2(i+1,j+1), c1, c2, c3, c4)

    !..convective scale velocity
    wst= rt1*interp(w_star1(i,j), w_star1(i+1,j), w_star1(i,j+1), w_star1(i+1,j+1), c1, c2, c3, c4) &
      +rt2*interp(w_star2(i,j), w_star2(i+1,j), w_star2(i,j+1), w_star2(i+1,j+1), c1, c2, c3, c4)

    pextra%ol=ol
    pextra%wst=wst
    pextra%ust=ust

  endif

  !..map ratio
  rmx= interp(xm(i,j), xm(i+1,j), xm(i,j+1), xm(i+1,j+1), c1, c2, c3, c4)
  rmy= interp(ym(i,j), ym(i+1,j), ym(i,j+1), ym(i+1,j+1), c1, c2, c3, c4)

  !..precipitation intensity (mm/hour)
  pr= interp(precip(i,j), precip(i+1,j), precip(i,j+1), precip(i+1,j+1), c1, c2, c3, c4)  

  !..update boundary layer top and height, map ratio and precipitation
  part%tbl=bl
  part%hbl=hbl
  pextra%rmx=rmx/dxgrid
  pextra%rmy=rmy/dygrid
  pextra%prc=pr


end subroutine posint

subroutine posint_vert(part, pextra, uprof, vprof, wprof, rhoprof, rhogradprof, k)
  
  USE particleML, only: Particle, extraParticle
  USE snapfldML, only: hlevel2

  type(Particle), intent(in) :: part
  type(extraParticle), intent(inout) :: pextra
  real, intent(in) :: uprof(:), vprof(:), wprof(:), rhoprof(:), rhogradprof(:)
  integer, intent(in) :: k
  
  real :: particle_z, below_level, above_level, frac
  integer :: i,j

  particle_z = part%zmetres
  i = part%x
  j = part%y

  below_level = hlevel2(i, j, k)
  above_level = hlevel2(i, j, k+1)

  if (above_level /= below_level) then
    frac = (particle_z - below_level) / (below_level)
  else
    frac = 0.0
  end if

  frac = max(0.0, min(1.0, frac))

  ! interpolate fields linearly at vertical position
  pextra%u = uprof(k) * (1.0-frac) + uprof(k+1) * frac
  pextra%v = vprof(k) * (1.0-frac) + vprof(k+1) * frac
  pextra%w = wprof(k) * (1.0-frac) + wprof(k+1) * frac
  pextra%rho = rhoprof(k) * (1.0-frac) + rhoprof(k+1) * frac
  pextra%rhograd = rhogradprof(k) * (1.0-frac) + rhogradprof(k+1) * frac

end subroutine posint_vert

subroutine posint_newlevel(part, pextra, k, uprof, vprof, wprof, rhoprof, rhogradprof, rt1, rt2)
  use particleML, only: Particle, extraParticle
  use snapdimML, only: nk
  USE snapfldML, only: u1, u2, v1, v2, w_z1, w_z2, rho1, rho2, rhograd1, rhograd2

  type(Particle), intent(in)    :: part
  type(extraParticle), intent(inout) :: pextra
  !> fractions from last timestep
  real, intent(in) :: rt1
  !> fractions to next timestep
  real, intent(in) :: rt2
  integer, intent(in) :: k
  real, intent(inout) :: uprof(:), vprof(:), wprof(:), rhoprof(:), rhogradprof(:)
  
  integer :: i,j
  real :: dx,dy,c1,c2,c3,c4

  if (k < 1 .or. k > nk) return

  !..for horizontal interpolations
  ! i,j = center of the (lower left) grid cell in which the particle is located
  ! dx,dy = distance of the particle to the center of the grid cell in grid units
  i=int(part%x)
  j=int(part%y)
  dx=part%x-i
  dy=part%y-j
  c1=(1.-dy)*(1.-dx)
  c2=(1.-dy)*dx
  c3=dy*(1.-dx)
  c4=dy*dx

  !..interpolation

  !.. u velocity
  uprof(k) = rt1*interp(u1(i,j,k), u1(i+1,j,k), u1(i,j+1,k), u1(i+1,j+1,k), c1, c2, c3, c4) &
            +rt2*interp(u2(i,j,k), u2(i+1,j,k), u2(i,j+1,k), u2(i+1,j+1,k), c1, c2, c3, c4)

  ! .. v velocity
  vprof(k) = rt1*interp(v1(i,j,k), v1(i+1,j,k), v1(i,j+1,k), v1(i+1,j+1,k), c1, c2, c3, c4) &
            +rt2*interp(v2(i,j,k), v2(i+1,j,k), v2(i,j+1,k), v2(i+1,j+1,k), c1, c2, c3, c4)

  ! .. w velocity
  wprof(k) = rt1*interp(w_z1(i,j,k), w_z1(i+1,j,k), w_z1(i,j+1,k), w_z1(i+1,j+1,k), c1, c2, c3, c4) &
            +rt2*interp(w_z2(i,j,k), w_z2(i+1,j,k), w_z2(i,j+1,k), w_z2(i+1,j+1,k), c1, c2, c3, c4)

  ! .. air density
  rhoprof(k) = rt1*interp(rho1(i,j,k), rho1(i+1,j,k), rho1(i,j+1,k), rho1(i+1,j+1,k), c1, c2, c3, c4) &
              +rt2*interp(rho2(i,j,k), rho2(i+1,j,k), rho2(i,j+1,k), rho2(i+1,j+1,k), c1, c2, c3, c4)

  ! .. air density gradient
  rhogradprof(k) = rt1*interp(rhograd1(i,j,k), rhograd1(i+1,j,k), rhograd1(i,j+1,k), rhograd1(i+1,j+1,k), c1, c2, c3, c4) &
              +rt2*interp(rhograd2(i,j,k), rhograd2(i+1,j,k), rhograd2(i,j+1,k), rhograd2(i+1,j+1,k), c1, c2, c3, c4)

end subroutine posint_newlevel

end module posintML
