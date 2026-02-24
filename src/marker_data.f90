MODULE marker_data
  integer :: max_markers

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
  
  ! Thermochronology
  double precision, allocatable :: mark_temp(:)          ! temperature
  double precision, allocatable :: mark_tempmax(:)       ! the max temperature have been
  double precision, allocatable :: mark_cooling_rate(:)  ! cooling rate
  double precision, allocatable :: mark_update_time(:)   ! the time at last thermochron calculate
  double precision, allocatable :: mark_chron_time(:,:)  ! closure time of thermochron (3, max)
  double precision, allocatable :: mark_chron_temp(:,:)  ! the closure temperatures of thermochron
  integer, allocatable :: mark_chron_if(:,:)             ! if closure of thermochron

  !$ACC declare create(max_markers)

  contains

  subroutine allocate_markers
    use params
    implicit none

    max_markers = nz * nx * max_markers_per_elem
    !$ACC update device(max_markers) async(1)

    allocate(mark_a1(max_markers), mark_a2(max_markers), &
             mark_x(max_markers), mark_y(max_markers), &
             mark_age(max_markers), &
             mark_dead(max_markers), &
             mark_ntriag(max_markers), &
             mark_phase(max_markers), &
             mark_ID(max_markers))

    allocate(mark_temp(max_markers), &
             mark_tempmax(max_markers), &
             mark_cooling_rate(max_markers), &
             mark_update_time(max_markers), &
             mark_chron_time(nchron, max_markers), &
             mark_chron_temp(nchron, max_markers), &
             mark_chron_if(nchron, max_markers))

    allocate(mark_id_elem(max_markers_per_elem, nz-1, nx-1))
    allocate(nmark_elem(nz-1, nx-1))

  end subroutine


  subroutine add_marker(x, y, iph, age, j, i, inc)
    !$ACC routine seq
    !$ACC routine(check_inside) seq
    
    use arrays
    use params
    implicit none
    integer :: iph, j, i, inc
    double precision :: x, y, age
    integer :: ntr, kk, nm
    double precision :: bar1, bar2
    integer :: n
    double precision :: t_closure

    !character*200 msg

    call check_inside(x , y, bar1, bar2, ntr, i, j, inc)
    if(inc.eq.0) return

    !$OMP atomic capture
    !$ACC atomic capture
    nmarkers = nmarkers + 1
    kk = nmarkers
    !$ACC end atomic
    !$OMP end atomic

    !$OMP atomic capture
    !$ACC atomic capture
    nmark_elem(j,i) = nmark_elem(j,i) + 1
    nm = nmark_elem(j,i)
    !$ACC end atomic
    !$OMP end atomic

    if(nm > max_markers_per_elem .or. kk > max_markers) then
        !write(msg*) 'Too many markers at element:', i, j, nm
        !call SysMsg(msg)
        !call SysMsg('Marker skipped, not added!')
        inc = -1
        return
    endif

    ! recording the id of markers belonging to the element
    mark_id_elem(nm,j,i) = kk

    mark_x(kk) = x
    mark_y(kk) = y
    mark_dead(kk) = 1
    mark_ID(kk) = kk
    mark_a1(kk) = bar1
    mark_a2(kk) = bar2
    mark_age(kk) = age
    mark_ntriag(kk) = ntr
    mark_phase(kk) = iph
   
    if (ithermochron > 0) then
        call temp2marker(kk)
        mark_update_time(kk) = time
        mark_cooling_rate(kk) = 0.0d0
        mark_tempmax(kk) = -1000.0d0

        ! Inherit thermochron state directly from element grid arrays
        if (i > 1 .and. i < nx-1) then
            mark_chron_time(:, kk) = chron_time(:, j, i)
            mark_chron_temp(:, kk) = chron_temp(:, j, i)
            mark_chron_if(:, kk)   = chron_if(:, j, i)
        else
            mark_chron_time(:, kk) = unreset_time
            mark_chron_temp(:, kk) = 0.0d0
            mark_chron_if(:, kk)   = 1
        endif

        do n = 1, nchron
            t_closure = chron_ref(n, 1, 2)
            if (mark_temp(kk) .gt. t_closure) then
                mark_chron_if(n, kk) = 0
                mark_chron_time(n, kk) = time
            endif
        enddo
    endif

  end subroutine add_marker


  subroutine newphase2marker (j1, j2, i1, i2, iph)
    !$ACC routine worker
    use arrays
    use params
    implicit none
    
    integer :: j1, j2, i1, i2, iph, &
               kk, n, j, i

    !$ACC loop collapse(2)
    !$OMP parallel do private(i,j,n,kk)
    ! reset the markers within elements in the rectangular region
    do i = i1, i2
      do j = j1, j2

        do n = 1 , nmark_elem(j,i)
          kk = mark_id_elem(n,j,i)
          mark_phase(kk) = iph
        enddo

        iphase(j,i) = iph
        phase_ratio(:,j,i) = 0.d0
        phase_ratio(iph,j,i) = 1.d0

      enddo
    enddo
    !$OMP end parallel do

    return
  end subroutine newphase2marker


  subroutine count_phase_ratio(j, i)
    !$ACC routine seq
    use params
    use arrays

    integer, intent(in) :: j, i
    integer :: n, kk, ncounters(maxph), iph, kph, nm

    if (nmark_elem(j,i) == 0) stop 'No markers in element'

    ncounters = 0
    do n = 1, nmark_elem(j,i)
      kk = mark_id_elem(n,j,i)
      iph = mark_phase(kk)
      ncounters(iph) = ncounters(iph) + 1
    enddo

    do kk = 1, nphase
      phase_ratio(kk,j,i) = ncounters(kk) / (nmark_elem(j,i) * 1.0d0)
    enddo

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
