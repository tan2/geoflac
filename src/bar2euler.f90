subroutine bar2euler
use arrays
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ xmpt(2,3,mnz*mnx*2)

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
        xmpt(1,1,n) = (x2*y3-y2*x3)/det
        xmpt(1,2,n) = (y2-y3)/det
        xmpt(1,3,n) = (x3-x2)/det
        xmpt(2,1,n) = (x3*y1-y3*x1)/det
        xmpt(2,2,n) = (y3-y1)/det
        xmpt(2,3,n) = (x1-x3)/det
enddo
enddo
enddo
!$OMP end do

! For brevity xmp1(1,2,n) --> m12 etc
!
! a1 = m11 + m12*x + m13*y
! a2 = m21 + m22*x + m23*y
!
! solve for x and y

!$OMP do
do i = 1 , nmarkers
     if (mark(i)%dead.eq.0) cycle
     n = mark(i)%ntriag
     ba1 = mark(i)%a1 
     ba2 = mark(i)%a2 
! Calculate eulerian from barycentic coordinates
     xnum = ba2*xmpt(1,3,n)-xmpt(2,1,n)*xmpt(1,3,n)-xmpt(2,3,n)*ba1+xmpt(2,3,n)*xmpt(1,1,n)
     xdem = xmpt(1,3,n)*xmpt(2,2,n)-xmpt(2,3,n)*xmpt(1,2,n)
     mark(i)%x   = xnum/xdem
     mark(i)%y   = (ba1 - xmpt(1,1,n) - xmpt(1,2,n)*(xnum/xdem))/xmpt(1,3,n)
enddo
!$OMP end do
!$OMP end parallel
return
end









 
