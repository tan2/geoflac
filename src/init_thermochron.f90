subroutine init_thermochron
use arrays
use marker_data
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

! Read the reference of thernochronology & set initial condition
  if (ithermochron .gt. 0) call read_thermochron_reference

  unreset_time = -1.*sec_year*1e6

  iwait = 0

  do kk = 1, nmarkers
    if (mark(kk)%dead.eq.0) cycle
    n = mark(kk)%ntriag
    k = mod(n - 1, 2) + 1
    j = mod((n - k) / 2, nz-1) + 1
    i = (n - k) / 2 / (nz - 1) + 1

    mark(kk)%chron_time(:) = unreset_time
    mark(kk)%chron_if(:) = 1
    call temp2marker(kk)
  end do

return
end subroutine init_thermochron


subroutine read_thermochron_reference
use arrays
include 'precision.inc'
include 'params.inc'

open(12, file=chron_file)

  do i = 1, nchron
      do
        call AdvanceToNextInputLine( 12 )
        read(12,*) num
        if (num /= 0) then
          do j = 1, nchron_fpair(i-1)
            call AdvanceToNextInputLine( 12 )
            read(12,*)
          end do
          cycle
        end if
        do
          call AdvanceToNextInputLine( 12 )
          read(12,*) num
          if (num /= 1) then
            do j = 1, nchron_fpair(i)
              call AdvanceToNextInputLine( 12 )
              read(12,*)
            end do
            cycle
          else
            call AdvanceToNextInputLine( 12 )
            do j = 1, nchron_fpair(i)
              read(12,*) chron_ref(i,j,:)
            end do
            exit
          end if
        end do
        exit
      end do
  end do

close(12)
return
end subroutine read_thermochron_reference

