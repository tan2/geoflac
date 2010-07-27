subroutine newphase2marker (ik,j,ntriang) 
USE marker_data
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

! reset the markers in this element

!!$DIR PREFER_PARALLEL
!!!$DIR LOOP_PARALLEL
!!!$DIR LOOP_PRIVATE(i,ntest,xntest,yntest,xxik,ik,xxi)

do i = 1 , nmarkers
    if (mark(i)%dead.eq.0) cycle
    n = mark(i)%ntriag
    if (n.ne.ntriang) cycle
    mark(i)%maps = aps (j,ik)
    mark(i)%meII = strainII(j,ik)
    mark(i)%mpres = stressI(j,ik)
    tmpr = 0.25*(temp(j,ik)+temp(j+1,ik)+temp(j,ik+1)+temp(j+1,ik+1))
    mark(i)%mtemp = tmpr
    nphase_counter(j,ik,mark(i)%phase) = nphase_counter(j,ik,mark(i)%phase) - 1
    mark(i)%phase = iphase(j,ik)
    nphase_counter(j,ik,mark(i)%phase) = nphase_counter(j,ik,mark(i)%phase) + 1
enddo

phase_ratio(j,ik,:) = 0.0
phase_ratio(j,ik,iphase(j,ik)) = 1.0

return
end
