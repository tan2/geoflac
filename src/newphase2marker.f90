subroutine newphase2marker (j,i,iph)
USE marker_data
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

! reset the markers in this element

iphase(j,i) = iph
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

! Calculate triangle number in which the markers belong
ntriang1 = 2 * ( (nz-1)*(i-1)+j-1) + 1
ntriang2 = 2 * ( (nz-1)*(i-1)+j-1) + 2

!XXX: add omp directive
do kk = 1 , nmarkers
    if (mark(kk)%dead.eq.0) cycle
    n = mark(kk)%ntriag
    if (n.ne.ntriang1 .or. n.ne.ntriang2) cycle
    mark(kk)%maps = aps (j,i)
    mark(kk)%meII = strainII(j,i)
    mark(kk)%mpres = stressI(j,i)
    mark(kk)%mtemp = tmpr
    nphase_counter(j,i,mark(kk)%phase) = nphase_counter(j,i,mark(kk)%phase) - 1
    mark(kk)%phase = iph
    nphase_counter(j,i,iph) = nphase_counter(j,i,iph) + 1
enddo

phase_ratio(j,i,:) = 0.0
phase_ratio(j,i,iph) = 1.0

return
end
