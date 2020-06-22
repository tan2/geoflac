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


  !$ACC declare create(mark_a1, mark_a2, mark_x, mark_y, mark_age, &
  !$ACC                mark_dead, mark_ntriag, mark_phase, mark_ID, &
  !$ACC                mark_id_elem, nmark_elem)
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

    !$ACC update device(mark_a1, mark_a2, mark_x, mark_y, mark_age, &
    !$ACC               mark_dead, mark_ntriag, mark_phase, mark_ID, &
    !$ACC               mark_id_elem, nmark_elem)

  end subroutine

  subroutine add_marker(x, y, iph, age, kk, j, i, inc)
    ! Add a marker at physical coordinate (x, y), with phase iph and age, to
    ! element (j, i). The current (before adding thsi marker) marker size
    ! is kk. If (x, y) is not within the element, inc is set to 0 and
    ! marker not added. Otherwise, marker is added to "mark" array and kk
    ! incremented by 1.
    !

    !$ACC routine seq
    !$ACC routine(check_inside) seq
    use arrays
    use params
    implicit none
    integer, intent(in) :: iph, j, i
    integer, intent(out) :: inc
    integer, intent(inout) :: kk
    double precision, intent(in) :: x, y, age
    integer :: ntr, nmark_local, kk_local
    double precision :: bar1, bar2
  
    !character*200 msg
  
    call check_inside(x , y, bar1, bar2, ntr, i, j, inc)
    if(inc.eq.0) return
  
    if(nmark_elem(j,i) == max_markers_per_elem) then
        !write(msg*) 'Too many markers at element:', i, j, nmark_elem(j,i)
        !call SysMsg(msg)
        !call SysMsg('Marker skipped, not added!')
        inc = 0
        return
    endif

    !$ACC atomic read
    kk_local = kk
    kk_local = kk_local + 1
    !$ACC atomic write
    kk = kk_local

    ! recording the id of markers belonging to the element
    !$ACC atomic read
    nmark_local = nmark_elem(j,i)
    nmark_local = nmark_local + 1
    !$ACC atomic write
    nmark_elem(j,i) = nmark_local
    mark_id_elem(nmark_local,j,i) = kk_local

    if(kk_local > max_markers) then
      !call SysMsg('ADD_MARKER: # of markers exceeds max. value. Please increase mark array size.')
      stop 15
    endif

    mark_x(kk_local) = x
    mark_y(kk_local) = y
    mark_dead(kk_local) = 1
    mark_ID(kk_local) = kk_local
    mark_a1(kk_local) = bar1
    mark_a2(kk_local) = bar2
    mark_age(kk_local) = age
    mark_ntriag(kk_local) = ntr
    mark_phase(kk_local) = iph

  end subroutine add_marker


  subroutine marker2elem
    use myrandom_mod
    use arrays
    use params
    implicit none
  
    integer :: kph(1), i, j, kinc, inc, iseed, nm
    double precision :: x1, x2, y1, y2, xx, yy, rx, ry
  
    !character*200 msg
  
    ! Interpolate marker properties into elements
    ! Find the triangle in which each marker belongs
  
    !$ACC parallel private(iseed, kph, nm)
    iseed = 0
    !$ACC loop collapse(2)
    do i = 1 , nx-1
        do j = 1 , nz-1
            kinc = nmark_elem(j,i)
  
            !  if there are too few markers in the element, create a new one
            !  with age 0 (similar to initial marker)
            !if(kinc.le.4) then
            !    write(msg,*) 'marker2elem: , create a new marker in the element (i,j))', i, j
            !    call SysMsg(msg)
            !endif
            do while (kinc.le.4)
                x1 = min(cord(j  ,i  ,1), cord(j+1,i  ,1))
                y1 = min(cord(j  ,i  ,2), cord(j  ,i+1,2))
                x2 = max(cord(j+1,i+1,1), cord(j  ,i+1,1))
                y2 = max(cord(j+1,i+1,2), cord(j+1,i  ,2))
  
                call myrandom(iseed, rx)
                call myrandom(iseed, ry)
  
                xx = x1 + rx*(x2-x1)
                yy = y1 + ry*(y2-y1)
  
                call add_marker(xx, yy, iphase(j,i), 0.d0, nmarkers, j, i, inc)
                if(inc.eq.0) cycle
  
                kinc = kinc + 1
            enddo
  
            call count_phase_ratio(j,i)
  
        enddo
    enddo
  
    !$ACC end parallel
    return
  end subroutine marker2elem
  

  subroutine newphase2marker (j1, j2, i1, i2, iph)
    !$ACC routine gang
    use arrays
    use params
    implicit none
    
    integer :: j1, j2, i1, i2, iph, &
               kk, n, j, i

    !$ACC loop collapse(2) private(i,j,n,kk)
    !$OMP parallel do private(i,j,n,kk)
    ! reset the markers within elements in the rectangular region
    do i = i1, i2
      do j = j1, j2
        !$ACC loop
        do n = 1 , nmark_elem(j,i)
          kk = mark_id_elem(n,j,i)
          mark_phase(kk) = iph
        enddo
        !$ACC end loop
      enddo
    enddo
    !$OMP end parallel do
    !$ACC end loop

    iphase(j1:j2,i1:i2) = iph
    phase_ratio(:,j1:j2,i1:i2) = 0.d0
    phase_ratio(iph,j1:j2,i1:i2) = 1.d0
    
    return
  end subroutine newphase2marker


  subroutine count_phase_ratio(j, i)
    !$ACC routine vector
    use params
    use arrays

    integer, intent(in) :: j, i
    integer :: n, kk, ncounters(maxph), iph, kph, nm

    ncounters = 0
    !$ACC loop
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
    iphase(j,i) = iph
  end subroutine count_phase_ratio
END MODULE marker_data
