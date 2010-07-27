subroutine init_marker

USE marker_data

use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
real(8), dimension(9) :: x_tr,y_tr
parameter( min_elmarkers = 0, max_elmarkers = 12 )

nphase_counter = 0
max_markers = 9*(nz-1)*(nx-1)

nelemts  = (nx-1)*(nz-1)
! define euler coordinate of the markers
! Distribute evenly first then randomize the distribution
! to start 9 markers per elements
nmarkers = 0
kk = 0
! zones with 9 markers per elements
! calculate the id (element number) of the zones of high res

!call random_seed
!write(333,*) 'Call to random_seed(), result may be stochastic'

do j = 1 , nz-1
    do i = 1 , nx-1
        dx = cord(j,i+1,1)-cord(j,i,1)
        dy = cord(j+1,i,2)-cord(j,i,2)

        x_tr(1) = cord(j,i,1) + dx*0.25
        y_tr(1) = cord(j,i,2) + dy*0.25
        x_tr(2) = cord(j,i,1) + dx*0.5
        y_tr(2) = cord(j,i,2) + dy*0.25
        x_tr(3) = cord(j,i,1) + dx*0.75
        y_tr(3) = cord(j,i,2) + dy*0.25

        x_tr(4) = cord(j,i,1) + dx*0.25
        y_tr(4) = cord(j,i,2) + dy*0.5
        x_tr(5) = cord(j,i,1) + dx*0.5
        y_tr(5) = cord(j,i,2) + dy*0.5
        x_tr(6) = cord(j,i,1) + dx*0.75
        y_tr(6) = cord(j,i,2) + dy*0.5

        x_tr(7) = cord(j,i,1) + dx*0.25
        y_tr(7) = cord(j,i,2) + dy*0.75
        x_tr(8) = cord(j,i,1) + dx*0.5
        y_tr(8) = cord(j,i,2) + dy*0.75
        x_tr(9) = cord(j,i,1) + dx*0.75
        y_tr(9) = cord(j,i,2) + dy*0.75


! randomize the new coordinates inside the element
        do l = 1 , 9 
            nmarkers = nmarkers + 1
            call random_number(rx)
            rx = 0.5 - rx
            ddx = dx*rx/5.
            ddy = dy*rx/5.
            x_tr(l) = x_tr(l)+ddx
            y_tr(l) = y_tr(l)+ddy
      
! define the markers for each elements in columns (9*nelemts=nmarkers) 
            kk = kk +1 
            mark(kk)%x = x_tr(l)
            mark(kk)%y = y_tr(l)
            mark(kk)%dead = 1 
            mark(kk)%ID = kk 
            xx = mark(kk)%x 
            yy = mark(kk)%y 
            ii = i
            jj = j
            call euler2bar(xx,yy,bar1,bar2,ntr,ii,jj,inc)
            if(ntr.eq.0) write(*,*) ii,jj,xx,yy
            mark(kk)%a1 = bar1
            mark(kk)%a2 = bar2
!            write(*,*) i,j,k,kk,xx,yy,mark(kk)%a1,mark(kk)%a2 
            mark(kk)%ntriag = ntr 
            mark(kk)%phase = iphase(j,i)

! XXX: hard-coded marker phase
! Special case of TAiwan thin layer of sediment accross the box until marge
            if (j.eq.1.and.l.le.3) mark(kk)%phase = 10.

            mark(kk)%maps   = aps(j,i)
            mark(kk)%meII   = 0. 
            mark(kk)%mpres  = 0.
            mark(kk)%mtemp  = 0.

            nphase_counter(j,i,mark(kk)%phase) = nphase_counter(j,i,mark(kk)%phase) + 1
        enddo
    enddo
enddo
write(333,*) '# of markers', nmarkers
return
end subroutine init_marker

