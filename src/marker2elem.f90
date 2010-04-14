subroutine marker2elem 
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
parameter( min_elmarkers = 0, max_elmarkers = 12 )
!type(marker):: mark(nmarkers)

kkk = nmarkers
xiph = 0.
! Interpolate marker properties into elements 
! Find the triangle in which each marker belongs
!old_aps = aps
o_phasez = phasez
ophase_ratio = phase_ratio
!aps = 0.
phasez = 0.
phase_ratio = 0.
do j = 1 , nz-1
do i = 1 , nx-1
        kinc = 0
        meml= 0
        phase_counter = 0.0
         do k = 1 , 2
            do l = 1 , kkk 
            if(mark(l)%dead.eq.0) cycle
            n = 2*( (nz-1)*(i-1) +j-1) + k 
            ntra = mark(l)%ntriag
            if (n.ne.ntra) then
               cycle 
            else
               kinc = kinc + 1

               phasez(j,i) = phasez(j,i)+ mark(l)%phase
               meml = l
! To linearly interpolate for density, conductivities, viscosities calculate phases ratios
              do m = 1,nphasl
         if(mark(l)%phase.eq.lphase(m)) phase_counter(m) = phase_counter(m) + 1
              enddo
             endif
          enddo ! m
          enddo ! k

!  if there are too few markers in the elmenent, create a new one
           if(kinc.le.2) then
           nmarkers = nmarkers+1
          kinc = kinc +1
               x1 = cord(j  ,i+1,1)
               x2 = cord(j+1,i  ,1)
               x3 = cord(j+1,i+1,1)
               y1 = cord(j  ,i+1,2)
               y2 = cord(j+1,i  ,2)
               y3 = cord(j+1,i+1,2)
               xx = x2 +0.5*(x3-x2)  
               mark(nmarkers)%x = xx
               yy = y1 +0.5*(y3-y1)  
               mark(nmarkers)%y = yy
 ! Calculate barycentic coordinates
               ii = i
               jj = j
               call euler2bar(xx,yy,bar1,bar2,ntr,ii,jj,inc) 
           mark(nmarkers)%a1 = bar1 
           mark(nmarkers)%a2 = bar2 
           mark(nmarkers)%maps = aps(j,i) 
           mark(nmarkers)%ntriag = ntr 
           mark(nmarkers)%dead = 1
           if(inc.eq.0) mark(nmarkers)%dead = 0
! Calculate strain
           mark(nmarkers)%meII = strainII(j,i)
! Calculate pressure and temperature
           tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
           mark(nmarkers)%mpres = stressI(j,i)
           mark(nmarkers)%mtemp = tmpr 
           mark(nmarkers)%ID = nmarkers 
! assign phase to the new marker
           if(kinc.gt.1) then
           mark(nmarkers)%phase = mark(meml)%phase 
           phasez(j,i) = phasez(j,i)+mark(nmarkers)%phase
!           write(*,*) i,j,xiph,n,ntr,xiph,kinc,mark(nmarkers)%phase
! if there were no more markers in the new element
           else 
!          xiph = 0.20*(o_phasez(j,i)+o_phasez(j,i-1)+o_phasez(j,i+1)+o_phasez(j-1,i)+o_phasez(j+1,i))
           xiph = o_phasez(j,i)
!	if (i.eq.1) xiph = (o_phasez(j-1,i)+o_phasez(j+1,i))/2.
!	if (i.eq.nx-1) xiph = (o_phasez(j-1,i)+o_phasez(j+1,i))/2.
!	if (j.eq.1) xiph = (o_phasez(j,i+1)+o_phasez(j,i-1))/2.
!	if (j.eq.nz-1) xiph = (o_phasez(j,i-1)+o_phasez(j,i+1))/2.
!	if (j.eq.1.and.i.eq.1) xiph = (o_phasez(j,i+1)+o_phasez(j+1,i)+o_phasez(j+1,i+1))/3.
!	if (j.eq.nz-1.and.i.eq.nx-1) xiph = (o_phasez(j,i-1)+o_phasez(j-1,i)+o_phasez(j-1,i-1))/3.
!	if (j.eq.1.and.i.eq.nx-1) xiph = (o_phasez(j,i-1)+o_phasez(j+1,i)+o_phasez(j+1,i-1))/3.
!	if (j.eq.nz-1.and.i.eq.1) xiph = (o_phasez(j,i+1)+o_phasez(j-1,i)+o_phasez(j-1,i+1))/3.
!	if (i.eq.1) xiph = o_phasez(j,5)
           mark(nmarkers)%phase = iphase(xiph)
           phasez(j,i) = phasez(j,i)+mark(nmarkers)%phase
!write(*,*) xiph,o_phasez(j,i-1),o_phasez(j,i+1),o_phasez(j-1,i),o_phasez(j+1,i)
!write(333,*) xiph,o_phasez(j,i-1),o_phasez(j,i+1),o_phasez(j-1,i),o_phasez(j+1,i)
          endif
         do m = 1,nphasl
         if(mark(nmarkers)%phase.eq.lphase(m)) phase_counter(m) = phase_counter(m) + 1
         enddo
         endif
do kn = 1 , nphasl
phase_ratio (j,i,kn) = phase_counter(kn)/kinc     
enddo
         phasez(j,i) = phasez(j,i)/kinc

enddo
enddo

return
end
