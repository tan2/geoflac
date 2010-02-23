subroutine lpeuler2bar
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
parameter( min_elmarkers = 0, max_elmarkers = 12 )
!type(marker):: mark(nmarkers)
!!!$DIR LOOP_PARALLEL
!!!$DIR LOOP_PRIVATE(k,xx,yy,bar1,bar2,ntr)
do k = 1 , nmarkers
    if (mark(k)%dead.eq.0) cycle
    nmtriag = mark(k)%ntriag
    ntr = nmtriag
    nmelemt = int((nmtriag-1)/2+1)
    if (mod(nmelemt,nz-1).eq.0) then
        ii = nmelemt/(nz-1)
    else
        ii = int((nmelemt/(nz-1))+1)
    endif
    jj = nmelemt-(nz-1)*(ii-1)
    xx = mark(k)%x
    yy = mark(k)%y
    call euler2bar(xx,yy,bar1,bar2,ntr,ii,jj,inc)
    if (inc.eq.0) then
        !write(*,*) ii,jj,ntr,xx,yy
        mark(k)%dead = 0
        bar1 = 1.e27
        bar2 = 1.e27
        ntr = 0
    endif
    mark(k)%a1 = bar1
    mark(k)%a2 = bar2
    mark(k)%ntriag = ntr
enddo
return

end subroutine lpeuler2bar
