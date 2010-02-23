subroutine euler2bar(x,y,bar1,bar2,ntr,ii,jj,inc)

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
parameter( min_elmarkers = 0, max_elmarkers = 12 )

! find the triangle in which the marker belongs
inc = 0 
kkk = 0

ibeg = ii-10
iend = ii+10
jbeg = jj-10
jend = jj+10
if (ibeg.le.0) ibeg = 1
if (jbeg.le.0) jbeg = 1
if (iend.ge.nx) iend = nx-1
if (jend.ge.nz) jend = nz-1
do j = jbeg, jend
do i = ibeg ,iend
do k = 1 , 2
ninters = 0
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
ymax1 = max(y1,y2)
ymax = max(ymax1,y3)
ymin1 = min(y1,y2)
ymin = min(ymin1,y3)
xmax1 = max(x1,x2)
xmax = max(xmax1,x3)
xmin1 = min(x1,x2)
xmin = min(xmin1,x3)

! check if point is above or under triangle
      if (y.gt.ymax) cycle
      if (y.lt.ymin) cycle

! check if point is on the left or the right of the triangle
      if (x.gt.xmax) cycle
      if (x.lt.xmin) cycle

! calculate function for each vertices
      if (x1.eq.x2) then
         xinter1 = x1
         goto 30
      endif
      xa1 = (y2 - y1)/(x2 - x1)
      xb1 = y1 - xa1*x1
      if (xa1.eq.0.) then
      xinter1 = 0.
      else
      xinter1 = (y-xb1)/xa1
      endif
30    continue
      if (x1.eq.x3) then
         xinter2 = x1
         goto 40
      endif
      xa2 = (y3 - y1)/(x3 - x1)
      xb2 = y1 - xa2*x1
      if(xa2.eq.0.) then
      xinter2 = 0.
      else
      xinter2 = (y-xb2)/xa2
      endif
40 continue
      if (x2.eq.x3) then
         xinter3 = x2
         goto 50
      endif
      xa3 = (y3 - y2)/(x3 - x2)
      xb3 = y2 - xa3*x2
      if (xa3.eq.0.) then
      xinter3 = 0.
      else
      xinter3 = (y-xb3)/xa3
      endif
50 continue
! calculate intersection of y with vertice for xi > x
      if (xinter1.lt.x) ninters = 0
      if (xinter1.ge.x.and.xinter1.le.xmax) ninters = ninters + 1 
      if (xinter2.lt.x) ninters = ninters 
      if (xinter2.ge.x.and.xinter2.le.xmax) ninters = ninters + 1 
      if (xinter3.lt.x) ninters = ninters 
      if (xinter3.ge.x.and.xinter3.le.xmax) ninters = ninters + 1 
! check for intersectons with nodes or vertices kill marker if it is the case
      if (y.eq.y1) then
         if (x.eq.x1) then
            ninters = 1
         endif
      endif 
      if (y.eq.y2) then
         if (x.eq.x2) then
            ninters = 1
         endif
      endif 
      if (y.eq.y3) then
         if (x.eq.x3) then
            ninters = 1
         endif
      endif 
! continue or stop if ninters is eq to 1
      if (ninters.eq.0.or.ninters.eq.2) cycle
! do nothing if marker is dead
!      if (mar
      if (ninters.eq.1) then
!         write(*,*) x1,x2,x3,y1,y2,y3,cord(j,i,1),cord(j,i+1,1),cord(j+1,i,1),cord(j,i,2),cord(j,i+1,2) 
         iiiiin = iiiiin +1
         inc = inc + 1
! Calculate triangle number in which marker belongs
         ntr = 2 * ( (nz-1)*(i-1)+j-1) + k
! Calculate triangle properties
        det=( (x2*y3-y2*x3) - (x1*y3-y1*x3) + (x1*y2-y1*x2) )
          n=ntr
!Find the parameters ONLY for 2 vertices
        xmpt(n,1,1) = (x2*y3-y2*x3)/det
        xmpt(n,1,2) = (y2-y3)/det
        xmpt(n,1,3) = (x3-x2)/det
        xmpt(n,2,1) = (x3*y1-y3*x1)/det
        xmpt(n,2,2) = (y3-y1)/det
        xmpt(n,2,3) = (x1-x3)/det
! Calculate barycentic coordinates
        bar1 = xmpt(n,1,1) + xmpt(n,1,2)*x + y*xmpt(n,1,3) 
        bar2 = xmpt(n,2,1) + xmpt(n,2,2)*x + y*xmpt(n,2,3) 
!        if(i.gt.92.and.j.lt.12.and.mark) write(*,*) nmarkers,det,x1,x2,x3,y1,y2,y3
!              goto 110
endif
enddo
enddo
enddo
!110 continue
!if (inc.eq.0) then
!x = x + 3.
!y = y - 3.
!if (kkk.gt.500) goto 113
!kkk = kkk + 1
!write(*,*) kkk,x,y
!goto 111
!endif
!113 continue
return
end
