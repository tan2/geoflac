MODULE marker_data
  integer, parameter :: max_markers = 10000000
  SAVE
  double precision, allocatable :: mark_a1(:), mark_a2(:) ! baricentric coordinates
  double precision, allocatable :: mark_x(:), mark_y(:)   ! Euler coordinates
  double precision, allocatable :: mark_age(:)            ! creation time
  integer, allocatable :: mark_dead(:)
  integer, allocatable :: mark_ntriag(:)      ! number of FE-triangle
  integer, allocatable :: mark_phase(:)
  integer, allocatable :: mark_ID(:)          ! unique ID-number

  !$ACC declare create(mark_a1, mark_a2, mark_x, mark_y, mark_age, &
  !$ACC                mark_dead, mark_ntriag, mark_phase, mark_ID)
  contains

  subroutine allocate_markers(nz, nx)
    implicit none

    integer, intent(in) :: nz, nx
    integer :: max_markers
    max_markers = nz * nx * 32

    allocate(mark_a1(max_markers), mark_a2(max_markers), &
             mark_x(max_markers), mark_y(max_markers), &
             mark_age(max_markers), &
             mark_dead(max_markers), &
             mark_ntriag(max_markers), &
             mark_phase(max_markers), &
             mark_ID(max_markers))

  end subroutine
END MODULE marker_data
