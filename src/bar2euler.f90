subroutine bar2euler
use arrays
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
!type(marker) :: mark(nmarkers)
! calculate the new paramters for the triangles

!$OMP parallel private(x1,x2,x3,y1,y2,y3,n,det,i,k,j,ba1,ba2,xnum,xdem)
!$OMP do
do i = 1 , nx-1
do j = 1 , nz-1
do k = 1 , 2

  if (k.eq.1) then
      x1 = cord(j  ,i  ,1)
      x2 = cord(j+1,i  ,1)
      x3 = cord(j  ,i+1,1)
      y1 = cord(j  ,i  ,2)
      y2 = cord(j+1,i  ,2)
      y3 = cord(j  ,i+1,2)
   else  !if (k.eq.2) then
      x1 = cord(j  ,i+1,1)
      x2 = cord(j+1,i  ,1)
      x3 = cord(j+1,i+1,1)
      y1 = cord(j  ,i+1,2)
      y2 = cord(j+1,i  ,2)
      y3 = cord(j+1,i+1,2)
   endif

! Calculate triangle number in which marker belongs
         n = 2 * ( (nz-1)*(i-1)+j-1) + k

! Calculate triangle properties
        det=( (x2*y3-y2*x3) - (x1*y3-y1*x3) + (x1*y2-y1*x2) )

!Find the parameters ONLY for 2 vertices
        xmpt(n,1,1) = (x2*y3-y2*x3)/det
        xmpt(n,1,2) = (y2-y3)/det
        xmpt(n,1,3) = (x3-x2)/det
        xmpt(n,2,1) = (x3*y1-y3*x1)/det
        xmpt(n,2,2) = (y3-y1)/det
        xmpt(n,2,3) = (x1-x3)/det
enddo
enddo
enddo
!$OMP end do

!$OMP do
do i = 1 , nmarkers
     if (mark(i)%dead.eq.0) cycle
     n = mark(i)%ntriag
     ba1 = mark(i)%a1 
     ba2 = mark(i)%a2 
! Calculate eulerian from barycentic coordinates
     xnum = ba2*xmpt(n,1,3)-xmpt(n,2,1)*xmpt(n,1,3)-xmpt(n,2,3)*ba1+xmpt(n,2,3)*xmpt(n,1,1)
     xdem = xmpt(n,1,3)*xmpt(n,2,2)-xmpt(n,2,3)*xmpt(n,1,2)
     mark(i)%x   = xnum/xdem
     mark(i)%y   = (ba1 - xmpt(n,1,1) - xmpt(n,1,2)*(xnum/xdem))/xmpt(n,1,3) 
enddo
!$OMP end do
!$OMP end parallel
return
end









 
