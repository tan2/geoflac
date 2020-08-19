MODULE marker_data
  integer, parameter :: max_markers = 10000000

  !!! maximum number of ELEMENTS !!!
  integer, parameter :: max_markers_per_elem=32

  SAVE
  double precision, allocatable :: mark_a1(:), mark_a2(:) ! baricentric coordinates
  double precision, allocatable :: mark_x(:), mark_y(:)   ! Euler coordinates
  double precision, allocatable :: mark_age(:)            ! creation time
  integer, allocatable :: mark_dead(:)
  integer, allocatable :: mark_ntriag(:)      ! number of FE-triangle
  integer, allocatable :: mark_phase(:)
  integer, allocatable :: mark_ID(:)          ! unique ID-number

  integer, allocatable :: mark_id_elem(:,:,:), nmark_elem(:,:)

  contains

  subroutine allocate_markers(nz, nx)
    implicit none

    integer, intent(in) :: nz, nx
    integer :: max_markers
    max_markers = nz * nx * max_markers_per_elem

    allocate(mark_a1(max_markers), mark_a2(max_markers), &
             mark_x(max_markers), mark_y(max_markers), &
             mark_age(max_markers), &
             mark_dead(max_markers), &
             mark_ntriag(max_markers), &
             mark_phase(max_markers), &
             mark_ID(max_markers))

    allocate(mark_id_elem(max_markers_per_elem, nz-1, nx-1))
    allocate(nmark_elem(nz-1, nx-1))

  end subroutine

  subroutine add_marker(x, y, iph, age, kk, j, i, inc)
    ! Add a marker at physical coordinate (x, y), with phase iph and age, to
    ! element (j, i). The current (before adding thsi marker) marker size
    ! is kk. If (x, y) is not within the element, inc is set to 0 and
    ! marker not added. Otherwise, marker is added to "mark" array and kk
    ! incremented by 1.
    !
    ! *** This subroutine is not thread-safe. DON'T CALL IT WITHIN
    ! *** OPENMP/OMP SECTION.
  
    use arrays
    use params
    implicit none
    integer :: iph, kk, j, i, inc
    double precision :: x, y, age
    integer :: ntr
    double precision :: bar1, bar2
  
    character*200 msg
  
    call check_inside(x , y, bar1, bar2, ntr, i, j, inc)
    if(inc.eq.0) return
  
    if(nmark_elem(j,i) == max_markers_per_elem) then
        !write(msg*) 'Too many markers at element:', i, j, nmark_elem(j,i)
        !call SysMsg(msg)
        !call SysMsg('Marker skipped, not added!')
        inc = -1
        return
    endif
    ! recording the id of markers belonging to the element
    nmark_elem(j,i) = nmark_elem(j,i) + 1
    mark_id_elem(nmark_elem(j,i),j,i) = kk + 1

    kk = kk + 1
  
    mark_x(kk) = x
    mark_y(kk) = y
    mark_dead(kk) = 1
    mark_ID(kk) = kk
    mark_a1(kk) = bar1
    mark_a2(kk) = bar2
    mark_age(kk) = age
    mark_ntriag(kk) = ntr
    mark_phase(kk) = iph

  end subroutine add_marker


  subroutine newphase2marker (j1, j2, i1, i2, iph)
    use arrays
    use params
    implicit none
    
    integer :: j1, j2, i1, i2, iph, &
               kk, n, j, i

    !$OMP parallel do private(i,j,n,kk)
    ! reset the markers within elements in the rectangular region
    do i = i1, i2
      do j = j1, j2

        do n = 1 , nmark_elem(j,i)
          kk = mark_id_elem(n,j,i)
          mark_phase(kk) = iph
        enddo

      enddo
    enddo
    !$OMP end parallel do

    iphase(j1:j2,i1:i2) = iph
    phase_ratio(:,j1:j2,i1:i2) = 0.d0
    phase_ratio(iph,j1:j2,i1:i2) = 1.d0
    
    return
  end subroutine newphase2marker


  subroutine count_phase_ratio(j, i, iph)
    use params
    use arrays

    integer, intent(in) :: j, i
    integer, intent(out) :: iph
    integer :: n, kk, ncounters(maxph), kph, nm

    ncounters = 0
    do n = 1, nmark_elem(j,i)
      kk = mark_id_elem(n,j,i)
      iph = mark_phase(kk)
      ncounters(iph) = ncounters(iph) + 1
    enddo

    phase_ratio(1:nphase,j,i) = ncounters(1:nphase) / (nmark_elem(j,i) * 1.0d0)

    ! the phase of this element is the most abundant marker phase
    !kph = maxloc(ncounters)
    iph = 1
    nm = 0
    do kk = 1, nphase
       if (ncounters(kk) > nm) then
           nm = ncounters(kk)
           iph = kk
       end if
    end do

  end subroutine count_phase_ratio
END MODULE marker_data
