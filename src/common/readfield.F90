! SNAP: Servere Nuclear Accident Programme
! Copyright (C) 1992-2021   Norwegian Meteorological Institute
! License: GNU General Public License v3.0 or later

!> Module with utilities to read files with fimex or netcdf
module readfieldML
#if defined(FIMEX)
  use readfield_fiML, only: readfield_fi
#endif
  use readfield_ncML, only: readfield_nc
  implicit none
  private
  public readfield, readfield_and_compute
  contains
  subroutine readfield(ftype, istep, backward, itimei, ihr1, ihr2, itimefi, ierror)
    USE datetime, only: datetime_t
    character(len=*), intent(in) :: ftype
!> current timestep (always positive), negative istep means reset
    integer, intent(in) :: istep
!> whether meteorology should be read backwards
    logical, intent(in) :: backward
!> minimal time-offset after itimei
    integer, value :: ihr1
!> maximal time-offset after itimei
    integer, value :: ihr2
!> time since last file input
    type(datetime_t), intent(in) :: itimei
!> final time (output)
    type(datetime_t), intent(out) :: itimefi
!> error (output)
    integer, intent(out) :: ierror

    if (ftype == 'netcdf') then
      call readfield_nc(istep, backward, itimei, ihr1, ihr2, itimefi, ierror)
#if defined(FIMEX)
    else if (ftype == 'fimex') then
      call readfield_fi(istep, backward, itimei, ihr1, ihr2, itimefi, ierror)
#endif
    else
      write(*,'(A)') "Error: unknown file type in readfield"
      ierror = 1
    end if
  end subroutine readfield

  subroutine readfield_and_compute(ftype, istep, backward, itimei, ihr1, ihr2, time_file, ierror)
    USE bldpML, only: bldp
    USE compheightML, only: compheight
    USE datetime, only: datetime_t
    USE snapdebug, only: iulog, idebug
    USE snapgrdML, only: gparam, igtype
    USE forwrdML, only: w_eta_to_m
    USE snapfldML, only: tke
    USE rwalkML, only: bl_id, air_density, diffusion_fields, turbulence_fields_required, &
                       diffusion_in_metres, scheme_is_tke
    USE, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    USE iso_fortran_env, only: error_unit
!> file type (netcdf or fimex)
    character(len=*), intent(in) :: ftype
!> current timestep (always positive), negative istep means reset
    integer, intent(in) :: istep
!> whether meteorology should be read backwards
    logical, intent(in) :: backward
!> minimal time-offset after itimei
    integer, value :: ihr1
!> maximal time-offset after itimei
    integer, value :: ihr2
!> time since last file input
    type(datetime_t), intent(in) :: itimei
!> final time (output)
    type(datetime_t), intent(out) :: time_file
!> error (output)
    integer, intent(out) :: ierror


    call readfield(ftype, istep, backward, itimei, ihr1, ihr2, time_file, ierror)
    if (idebug >= 1) then
      write (iulog, *) "igtype, gparam(8): ", igtype, gparam
    end if
    if (ierror /= 0) then
      write (iulog, *) "Error in readfield, exiting"
      error stop "Error in readfield, exiting"
    end if

    !..compute model level heights
    call compheight

    !..calculate boundary layer (top and height)
    if (bl_id == 0) then
        call bldp
    endif

    ! Compute required fields for turbulence parametrisations
    if (turbulence_fields_required) then
      call air_density 
      call diffusion_fields
    endif

    ! Clean TKE values
    if (scheme_is_tke) then
      tke(:, :, 1) = tke(:, :, 2)
      where (ieee_is_nan(tke))
        tke = 0.0
      end where
    endif

    ! If diffusion scheme operates in metre space, convert vertical velocity to m/s
    if (diffusion_in_metres) call w_eta_to_m

  end subroutine readfield_and_compute

end module readfieldML
