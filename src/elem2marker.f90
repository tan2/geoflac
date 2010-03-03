subroutine elem2marker 
USE marker_data
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
parameter( min_elmarkers = 0, max_elmarkers = 12 )
!type(marker) :: mark(nmarkers)

! Interpolate element properties into markers 
! Find the element in which each marker belongs

!$OMP Parallel do private(i,j,n,ntest,xntest,yntest,k,xxik,ik,xxi,tmpr)

do i = 1 , nmarkers
    if (mark(i)%dead.eq.0) cycle
    n = mark(i)%ntriag
! if triangle number is pair or impair
     ntest = n/2
     xntest = float(n)/2
     yntest = float(ntest)
     if ((yntest-xntest).eq.0) then
         k = 2
     else
         k = 1
     endif
     do j = 1, nz - 1
     xxik = (1. + ((float(n)-float(k))/2-float(j) + 1)/(float(nz)-1))
     ik   = int(xxik)
     xxi  = float(ik)
        if ((abs(xxik)-abs(xxi)).eq.0) goto 33
     enddo
33            mark(i)%maps = aps (j,ik)
            mark(i)%meII = strainII(j,ik)
            mark(i)%mpres = stressI(j,ik)
            tmpr = 0.25*(temp(j,ik)+temp(j+1,ik)+temp(j,ik+1)+temp(j+1,ik+1))
            mark(i)%mtemp = tmpr
enddo
!$OMP end parallel do
return
end









 
