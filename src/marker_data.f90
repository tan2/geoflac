MODULE marker_data
  integer, parameter :: max_markers = 10000000
  SAVE
  type marker
     sequence
     real(8) :: a1(max_markers), a2(max_markers) ! baricentric coordinates
     real(8) :: x(max_markers), y(max_markers)   ! Euler coordinates
     real(8) :: age(max_markers)           ! creation time
     integer(4) :: dead(max_markers)
     integer(4) :: ntriag(max_markers)     ! number of FE-triangle
     integer(4) :: phase(max_markers)
     integer(4) :: ID(max_markers)         ! unique ID-number
  end type marker

  type(marker) :: mark

END MODULE marker_data
