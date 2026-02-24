subroutine init_thermochron
use arrays
use marker_data
use params
implicit none

integer :: kk, i, j, k, n

! Read the reference of thernochronology & set initial condition
  if (ithermochron .gt. 0) call read_thermochron_reference

  unreset_time = -100.*sec_year*1e6

  iwait = 0

  do kk = 1, nmarkers
    if (mark_dead(kk).eq.0) cycle
    n = mark_ntriag(kk)
    k = mod(n - 1, 2) + 1
    j = mod((n - k) / 2, nz-1) + 1
    i = (n - k) / 2 / (nz - 1) + 1

    mark_chron_time(:, kk) = unreset_time
    mark_chron_if(:, kk) = 1
    
    ! Initialize variables correctly
    mark_temp(kk) = 0.0d0
    mark_tempmax(kk) = -1000.0d0
    mark_update_time(kk) = time
    
    call temp2marker(kk)
  end do

return
end subroutine init_thermochron


subroutine read_thermochron_reference
use arrays
use params
implicit none

integer :: i, j, num, line, ios

open(12, file=chron_file)

  do i = 1, nchron
      do
        call AdvanceToNextInputLine( 12, line )
        read(12,*,iostat=ios) num
        if (ios .ne. 0) exit 
        if (num == 0) exit
      end do
      
      
      do j = 1, nchron_fpair(i)
          call AdvanceToNextInputLine( 12, line )
          read(12,*,iostat=ios) chron_ref(i,j,:)
          if (ios .ne. 0) then
              print *, 'Error reading thermochron data for chron ', i, ' point ', j
              stop
          endif
      end do
  end do

close(12)
return
end subroutine read_thermochron_reference
