subroutine lpeuler2bar
USE marker_data
use arrays
use params

include 'precision.inc'

character*200 msg

nphase_counter(:,:,:) = 0
ntopmarker(:) = 0

!$OMP parallel do private(n,k,j,i,xx,yy,bar1,bar2,ntr,inc)
do n = 1 , nmarkers
    if (mark%dead(n).eq.0) cycle

    ! from ntriag, get element number
    k = mod(mark%ntriag(n) - 1, 2) + 1
    j = mod((mark%ntriag(n) - k) / 2, nz-1) + 1
    i = (mark%ntriag(n) - k) / 2 / (nz - 1) + 1

    xx = mark%x(n)
    yy = mark%y(n)
    call euler2bar(xx,yy,bar1,bar2,ntr,i,j,inc)
    if (inc.eq.0) then
        !write(*,*) i,j,ntr,xx,yy
        mark%dead(n) = 0
        cycle
    endif

    mark%a1(n) = bar1
    mark%a2(n) = bar2
    mark%ntriag(n) = ntr

    !$OMP critical (lpeulerbar1)
    nphase_counter(mark%phase(n),j,i) = nphase_counter(mark%phase(n),j,i) + 1
    !$OMP end critical (lpeulerbar1)

    if(j == 1) then
        if(ntopmarker(i) == max_markers_per_elem) then
            write(msg,*) 'Too many markers at surface elements in lpeuler2bar:', i, n
            call SysMsg(msg)
            mark%dead(n) = 0
            cycle
        endif
        !$OMP critical (lpeulerbar2)
        ! recording the id of markers belonging to surface elements
        ntopmarker(i) = ntopmarker(i) + 1
        itopmarker(ntopmarker(i), i) = n
        !$OMP end critical (lpeulerbar2)
    end if
enddo
!$OMP end parallel do
return

end subroutine lpeuler2bar
