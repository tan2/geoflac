subroutine temp2marker(kk)
use marker_data
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

real*8 :: tt(3)
real*8 :: bar(3)

! from ntriag, get element number
  n = mark(kk)%ntriag
  k = MOD(n - 1, 2) + 1
  j = MOD((n - k) / 2, nz-1) + 1
  i = (n - k) / 2 / (nz - 1) + 1

  bar(1) = mark(kk)%a1
  bar(2) = mark(kk)%a2
  bar(3) = 1.d0 - bar(1) - bar(2)

  if (MOD(n,2).eq.1) then
    tt(1) = temp(j  ,i  )
    tt(2) = temp(j+1,i  )
    tt(3) = temp(j  ,i+1)
  else
    tt(1) = temp(j  ,i+1)
    tt(2) = temp(j+1,i  )
    tt(3) = temp(j+1,i+1)
  endif

  t_temp = sum(tt*bar)

  ! record the max temperature the rock have been
  if (t_temp .gt. mark(kk)%tempmax) mark(kk)%tempmax = t_temp

  t_temp_last = mark(kk)%temp
  delta_temp = t_temp - t_temp_last
  delta_time = time - mark(kk)%update_time

  mark(kk)%temp = t_temp
  mark(kk)%update_time = time

!  if ( mark(kk)%cooling_rate .lt. -1000. .and. mark(kk)%cooling_rate .gt. 0.) then
!    print*, nloop,i,j,mark(kk)%id,mark(kk)%temp,t_temp_last,delta_temp, delta_time, mark(kk)%cooling_rate
!  end if

  if (iwait.eq.1) return

  ! Calculate the cooling rate
  mark(kk)%cooling_rate = delta_temp*sec_year* 1.d6/delta_time

return
end subroutine temp2marker