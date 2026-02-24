subroutine temp2marker(kk)
use marker_data
use arrays
use params
implicit none

integer :: kk
real*8 :: tt(3)
real*8 :: bar(3)
real*8 :: t_temp, t_temp_last, delta_temp, delta_time
integer :: n, k, j, i

! from ntriag, get element number
  n = mark_ntriag(kk)
  k = MOD(n - 1, 2) + 1
  j = MOD((n - k) / 2, nz-1) + 1
  i = (n - k) / 2 / (nz - 1) + 1

  bar(1) = mark_a1(kk)
  bar(2) = mark_a2(kk)
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
  if (t_temp .gt. mark_tempmax(kk)) mark_tempmax(kk) = t_temp

  t_temp_last = mark_temp(kk)
  delta_temp = t_temp - t_temp_last
  delta_time = time - mark_update_time(kk)

  mark_temp(kk) = t_temp
  mark_update_time(kk) = time

!  if ( mark_cooling_rate(kk) .lt. -1000. .and. mark_cooling_rate(kk) .gt. 0.) then
!    print*, nloop,i,j,mark_id(kk),mark_temp(kk),t_temp_last,delta_temp, delta_time, mark_cooling_rate(kk)
!  end if

  if (iwait.eq.1) return

  ! Calculate the cooling rate
  if (abs(delta_time) > 1.d-10) then
      mark_cooling_rate(kk) = delta_temp*sec_year* 1.d6/delta_time
  else
      mark_cooling_rate(kk) = 0.0d0
  endif

return
end subroutine temp2marker