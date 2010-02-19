MODULE marker_data
  SAVE
  type marker
     sequence
     real(8) :: a1,a2         ! baricentric coordinates
     real(8) :: x,y           ! Euler coordinates
     real(8) :: maps
     real(8) :: meII
     real(8) :: mpres
     real(8) :: mtemp
     integer(4) :: ID         ! unique ID-number
     integer(4) :: ntriag     ! number of FE-triangle
     integer(4) :: phase
     integer(4) :: dead
  end type marker

  type(marker) :: mark(10000000)

END MODULE marker_data
