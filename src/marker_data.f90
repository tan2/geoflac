MODULE marker_data
  integer, parameter :: max_markers = 10000000
  SAVE
  real(8) :: mark_a1(max_markers), mark_a2(max_markers) ! baricentric coordinates
  real(8) :: mark_x(max_markers), mark_y(max_markers)   ! Euler coordinates
  real(8) :: mark_age(max_markers)           ! creation time
  integer(4) :: mark_dead(max_markers)
  integer(4) :: mark_ntriag(max_markers)     ! number of FE-triangle
  integer(4) :: mark_phase(max_markers)
  integer(4) :: mark_ID(max_markers)         ! unique ID-number

END MODULE marker_data
