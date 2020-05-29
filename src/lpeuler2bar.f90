subroutine lpeuler2bar
USE marker_data
use arrays
use params
use euler2bar
implicit none

integer :: n, k, i, j, ntr, inc
double precision :: xx, yy, bar1, bar2
character*200 msg

!$ACC parallel
nphase_counter(:,:,:) = 0
ntopmarker(:) = 0
mark_id_elem(:,:,:) = 0
nmark_elem(:,:) = 0

!$ACC loop
!$OMP parallel do private(n,k,j,i,xx,yy,bar1,bar2,ntr,inc)
do n = 1 , nmarkers
    if (mark_dead(n).eq.0) cycle

    ! from ntriag, get element number
    k = mod(mark_ntriag(n) - 1, 2) + 1
    j = mod((mark_ntriag(n) - k) / 2, nz-1) + 1
    i = (mark_ntriag(n) - k) / 2 / (nz - 1) + 1

    xx = mark_x(n)
    yy = mark_y(n)
    call euler2bar(xx,yy,bar1,bar2,ntr,i,j,inc)
    if (inc.eq.0) then
        !write(*,*) i,j,ntr,xx,yy
        mark_dead(n) = 0
        cycle
    endif

    mark_a1(n) = bar1
    mark_a2(n) = bar2
    mark_ntriag(n) = ntr

    !$OMP critical (lpeulerbar1)
    !$ACC atomic update
    nphase_counter(mark_phase(n),j,i) = nphase_counter(mark_phase(n),j,i) + 1
    !$OMP end critical (lpeulerbar1)

    if(nmark_elem(j, i) == max_markers_per_elem) then
        !write(msg,*) 'Too many markers at elements in lpeuler2bar:', i, j, n
        !call SysMsg(msg)
        mark_dead(n) = 0
        cycle
    endif
    !$OMP critical (lpeulerbar2)
    ! recording the id of markers belonging to surface elements
    !$ACC atomic update
    nmark_elem(j, i) = nmark_elem(j, i) + 1
    !$ACC atomic write
    mark_id_elem(nmark_elem(j, i), j, i) = n
    !$OMP end critical (lpeulerbar2)
enddo
!$OMP end parallel do
!$ACC end parallel
return

end subroutine lpeuler2bar
