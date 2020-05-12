MODULE marker_data
  integer, parameter :: max_markers = 10000000
  SAVE
  type marker
     sequence
     real(8) :: a1,a2         ! baricentric coordinates
     real(8) :: x,y           ! Euler coordinates
     real(8) :: age           ! creation time
     integer(4) :: dead
     integer(4) :: ntriag     ! number of FE-triangle
     integer(4) :: phase
     integer(4) :: ID         ! unique ID-number

     real(8) :: temp          ! temperature
     real(8) :: tempmax       ! the max temperature have been
     real(8) :: cooling_rate  ! cooling rate
     real(8) :: update_time   ! the time at last thermochron calculate
     real(8) :: chron_time(3)  ! closure time of thermochron
     real(8) :: chron_temp(3) ! the closure temperatures of thermochron
     integer(4) :: chron_if(3) ! if closure of thermochron
  end type marker

  type(marker) :: mark(max_markers)

END MODULE marker_data
