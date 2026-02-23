subroutine fl_thermchron
!For Thermochronology
!
!Chase Shyu <iamhemry@gmail.com>
!Jan. 14th, 2013
!Refactored for geoflac-master Jan 2026

use arrays
use marker_data
use params
implicit none

character :: screen*200
integer :: kk, i, i_closure_if
double precision :: t_temp, t_rate, t_closure, delta_temp, sterr
double precision :: t_temp_last

sterr = 0.d-4
!12 format(1I10,5F15.7)

!if (ireset.eq.1) then
!  iwait = 1
!  return
!end if

if (mod(nloop, 10).ne.0) return
if (i_prestress .eq. 1 ) return

! Update present temperature on marker
!$OMP Parallel do default(none) shared(nmarkers,mark_dead) private(kk)
do kk = 1, nmarkers
  if (mark_dead(kk).eq.0) cycle
  call temp2marker(kk)
end do
!$OMP end parallel do

! do not calculate the closure temperature just after remeshing
if (iwait .eq. 1) then
  iwait = 0
  return
end if

do i = 1, nchron

!$OMP Parallel do default(none) private(kk,t_temp,t_rate,t_temp_last, &
!$OMP                           i_closure_if,t_closure,screen,delta_temp) &
!$OMP                           shared(mark_dead,mark_chron_if,mark_temp,mark_cooling_rate, &
!$OMP                           mark_chron_temp,mark_chron_time, &
!$OMP                           chron_ref,time,dt,i,nloop, &
!$OMP                           sterr,nmarkers,nchron_fpair)
! calculate if closure and assign age
    do kk = 1, nmarkers
      if (mark_dead(kk).eq.0) cycle

      i_closure_if = mark_chron_if(i, kk)

!      if (iupdate_temp_rate .le. 0) cycle
      t_temp = mark_temp(kk)
      t_rate = mark_cooling_rate(kk)
      ! find the closure temperature of t_rate
      if(t_rate .lt. 0.d0)then
        call interpolate_closure(-1.d0*t_rate,t_closure,i)
      else
        t_closure = chron_ref(i,1,2)
      end if
      ! for the marker which already closure
      if (t_temp - t_closure .le. sterr .and. i_closure_if .eq. 1) then
        mark_chron_temp(i, kk) = t_closure
        cycle
      else if (t_closure - t_temp .le. sterr .and. i_closure_if .eq. 0)then
        mark_chron_temp(i, kk) = t_closure
        mark_chron_time(i, kk) = time
        cycle
      end if

      ! for the rock that still too hot to preserve the track
      if (t_temp - t_closure .gt. sterr) then
        i_closure_if = 0
        mark_chron_time(i, kk) = time
      ! for the rock that have been heated and cooled down,
      else if (t_closure - t_temp .gt. sterr) then
        i_closure_if = 1
        mark_chron_time(i, kk) = time
      end if

      mark_chron_if(i, kk) = i_closure_if
      mark_chron_temp(i, kk) = t_closure
    end do
!$OMP end parallel do

end do

return
end subroutine fl_thermchron


subroutine interpolate_closure(t_rate,t_closure,i_ref)
use params
implicit none

double precision :: t_rate, t_closure
integer :: i_ref

integer :: low, high, mid
double precision :: r_low, r_high, c_low, c_high

low = 1
high = nchron_fpair(i_ref)

! Boundary checks (Assuming rates are sorted ascendingly)
if (t_rate <= chron_ref(i_ref, low, 1)) then
    t_closure = chron_ref(i_ref, low, 2)
    return
elseif (t_rate >= chron_ref(i_ref, high, 1)) then
    t_closure = chron_ref(i_ref, high, 2)
    return
endif

! Binary search
do while (high - low > 1)
    mid = (low + high) / 2
    if (chron_ref(i_ref, mid, 1) > t_rate) then
        high = mid
    else
        low = mid
    endif
end do

! Linear Interpolation
r_low = chron_ref(i_ref, low, 1)
r_high = chron_ref(i_ref, high, 1)
c_low = chron_ref(i_ref, low, 2)
c_high = chron_ref(i_ref, high, 2)

t_closure = c_low + (t_rate - r_low) * (c_high - c_low) / (r_high - r_low)

return
end subroutine interpolate_closure
