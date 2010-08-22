subroutine elem2marker 
USE marker_data
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

! Interpolate element properties into markers 
! Find the element in which each marker belongs

!$OMP Parallel do private(i,j,k,n,tmpr)

do n = 1 , nmarkers
    if (mark(n)%dead.eq.0) cycle

    ! from ntriag, get element number
    k = mod(mark(n)%ntriag - 1, 2) + 1
    j = mod((mark(n)%ntriag - k) / 2, nz-1) + 1
    i = (mark(n)%ntriag - k) / 2 / (nz - 1) + 1

    mark(n)%maps = aps(j,i)
    mark(n)%meII = strainII(j,i)
    mark(n)%mpres = stressI(j,i)
    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
    mark(n)%mtemp = tmpr
enddo
!$OMP end parallel do
return
end









 
