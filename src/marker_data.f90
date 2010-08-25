MODULE marker_data
  SAVE
  type marker
     sequence
     integer(4) :: dead
     integer(4) :: ntriag     ! number of FE-triangle
     integer(4) :: phase
     integer(4) :: ID         ! unique ID-number
     real(8) :: a1,a2         ! baricentric coordinates
     real(8) :: x,y           ! Euler coordinates
  end type marker

  type(marker) :: mark(10000000)

END MODULE marker_data
