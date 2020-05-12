subroutine fl_thermchron
!For Thermochronology
!
!Chase Shyu <iamhemry@gmail.com>
!Jan. 14th, 2013

use arrays
use marker_data
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

character :: screen*200

sterr = 0.d-4
!12 format(1I10,5F15.7)

if (ireset.eq.1) then
  iwait = 1
  return
end if

if (mod(nloop, 10).ne.0) return
if (i_prestress .eq. 1 ) return

! Update present temperature on marker
do kk = 1, nmarkers
  if (mark(kk)%dead.eq.0) cycle
  call temp2marker(kk)
end do

! do not calculate the closure temperature just after remeshing
if (iwait .eq. 1) then
  iwait = 0
  return
end if

do i = 1, nchron

!$OMP Parallel do default(none) private(kk,t_temp,t_rate,t_temp_last, &
!$OMP                           i_closure_if,t_closure,screen,delta_temp) &
!$OMP                           shared(mark, chron_ref,time,dt,sec_year,i,nloop, &
!$OMP                           sterr)
! calculate if closure and assign age
    do kk = 1, nmarkers
      if (mark(kk)%dead.eq.0) cycle

      i_closure_if = mark(kk)%chron_if(i)

!      if (iupdate_temp_rate .le. 0) cycle
      t_temp = mark(kk)%temp
      t_rate = mark(kk)%cooling_rate
      ! find the closure temperature of t_rate
      if(t_rate .lt. 0.d0)then
        call interpolate_closure(-1.d0*t_rate,t_closure,i)
      else
        t_closure = chron_ref(i,1,2)
      end if
      ! for the marker which already closure
      if (t_temp - t_closure .le. sterr .and. i_closure_if .eq. 1) then
        mark(kk)%chron_temp(i) = t_closure
        cycle
      else if (t_closure - t_temp .le. sterr .and. i_closure_if .eq. 0)then
        mark(kk)%chron_temp(i) = t_closure
        mark(kk)%chron_time(i) = time
        cycle
      end if

      ! for the rock that still too hot to preserve the track
      if (t_temp - t_closure .gt. sterr) then
        i_closure_if = 0
        mark(kk)%chron_time(i) = time
      ! for the rock that have been heated and cooled down,
      else if (t_closure - t_temp .gt. sterr) then
        i_closure_if = 1
        mark(kk)%chron_time(i) = time
      end if

      mark(kk)%chron_if(i) = i_closure_if
      mark(kk)%chron_temp(i) = t_closure
    end do
!$OMP end parallel do

end do

return
end subroutine fl_thermchron


subroutine interpolate_closure(t_rate,t_closure,i_ref)

include 'precision.inc'
include 'params.inc'
integer*4 :: i_rate(3)

i_rate(1) = 1
i_rate(2) = nchron_fpair(i_ref)
i_rate(3) = int(nchron_fpair(i_ref)/2.d0,4)

!Judge if t_rate is lower or hugher than reference
if(chron_ref(i_ref,i_rate(1),1) > t_rate) then
  t_closure = chron_ref(i_ref,i_rate(1),2)
  return
else if(chron_ref(i_ref,i_rate(2),1) < t_rate) then
  t_closure = chron_ref(i_ref,i_rate(2),2)
  return
end if
do
  if (chron_ref(i_ref,i_rate(3),1) > t_rate) then
    i_rate(1) = i_rate(1)
    i_rate(2) = i_rate(3)
    i_rate(3) = idnint(sum(i_rate(1:2))/2.d0)
  else
    i_rate(1) = i_rate(3)
    i_rate(2) = i_rate(2)
    i_rate(3) = idnint(sum(i_rate(2:3))/2.d0)
  end if
  if ((i_rate(3) == i_rate(1)) .or. (i_rate(3) == i_rate(2))) exit
end do

delta_rate = chron_ref(i_ref,i_rate(2),1) - chron_ref(i_ref,i_rate(1),1)
delta_closure = chron_ref(i_ref,i_rate(2),2) - chron_ref(i_ref,i_rate(1),2)
rate_interp_ratio = (t_rate - chron_ref(i_ref,i_rate(1),1))/delta_rate
t_closure = chron_ref(i_ref,i_rate(1),2) + rate_interp_ratio * delta_closure

return
end subroutine interpolate_closure

