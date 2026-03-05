! SNAP: Servere Nuclear Accident Programme
! Copyright (C) 1992-2026   Norwegian Meteorological Institute
! License: GNU GPL v3 or later

module posintML
  implicit none
  private

  public :: posint

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
  ! i,j = lower left corner
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

  !..friction velocity
  ust= rt1*interp(u_star1(i,j), u_star1(i+1,j), u_star1(i,j+1), u_star1(i+1,j+1), c1, c2, c3, c4) &
    +rt2*interp(u_star2(i,j), u_star2(i+1,j), u_star2(i,j+1), u_star2(i+1,j+1), c1, c2, c3, c4)

  !..obukhov length
  ol= rt1*interp(obukhov_l1(i,j), obukhov_l1(i+1,j), obukhov_l1(i,j+1), obukhov_l1(i+1,j+1), c1, c2, c3, c4) &
    +rt2*interp(obukhov_l2(i,j), obukhov_l2(i+1,j), obukhov_l2(i,j+1), obukhov_l2(i+1,j+1), c1, c2, c3, c4)

  !..convective scale velocity
  wst= rt1*interp(w_star1(i,j), w_star1(i+1,j), w_star1(i,j+1), w_star1(i+1,j+1), c1, c2, c3, c4) &
    +rt2*interp(w_star2(i,j), w_star2(i+1,j), w_star2(i,j+1), w_star2(i+1,j+1), c1, c2, c3, c4)

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
  pextra%ol=ol
  pextra%wst=wst
  pextra%ust=ust

end subroutine posint
end module posintML
