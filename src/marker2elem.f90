subroutine marker2elem 
  USE marker_data

  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
  common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
  parameter( min_elmarkers = 0, max_elmarkers = 12 )
  !type(marker):: mark(nmarkers)

  ! Interpolate marker properties into elements
  ! Find the triangle in which each marker belongs


  do j = 1 , nz-1
      do i = 1 , nx-1
          kinc = sum(nphase_counter(j,i,1:nphase))

          !  if there are too few markers in the elmenent, create a new one
          do while (kinc.le.2)
              call SysMsg('marker2elem: too few markers in the elmenent, create a new one')
              nmarkers = nmarkers + 1
              kinc = kinc + 1
              x1 = cord(j  ,i+1,1)
              x2 = cord(j+1,i  ,1)
              x3 = cord(j+1,i+1,1)
              y1 = cord(j  ,i+1,2)
              y2 = cord(j+1,i  ,2)
              y3 = cord(j+1,i+1,2)

              xx = x2 + 0.5*(x3-x2)
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
              mark(nmarkers)%phase = iphase(j,i)
              nphase_counter(j,i,mark(nmarkers)%phase) = nphase_counter(j,i,mark(nmarkers)%phase) + 1
          enddo

          phase_ratio(j,i,1:nphase) = nphase_counter(j,i,1:nphase) / float(kinc)

      enddo
  enddo

  return
end subroutine marker2elem
