subroutine newphase2marker (ik,j,ntriang) 
USE marker_data
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
parameter( min_elmarkers = 0, max_elmarkers = 12 )
!type(marker) :: mark(nmarkers)

! Interpolate element properties into markers 
! Find the element in which each marker belongs
!!$DIR PREFER_PARALLEL
!!!$DIR LOOP_PARALLEL
!!!$DIR LOOP_PRIVATE(i,ntest,xntest,yntest,xxik,ik,xxi)

do i = 1 , nmarkers
      if (mark(i)%dead.eq.0) cycle
    n = mark(i)%ntriag
         if (n.ne.ntriang) cycle
            mark(i)%maps = aps (j,ik)
            mark(i)%meII = strainII(ik,j)
            mark(i)%mpres = stressI(ik,j)
            tmpr = 0.25*(temp(j,ik)+temp(j+1,ik)+temp(j,ik+1)+temp(j+1,ik+1))
            mark(i)%mtemp = tmpr
     if (phasez(j,ik).eq.5..or.phasez(j,ik).eq.8..or.phasez(j,ik).eq.9. &
          .or.phasez(j,ik).eq.1..or.phasez(j,ik).eq.7. &
           .or.phasez(j,ik).eq.12.or.phasez(j,ik).eq.15. ) then
                mark(i)%phase = phasez(j,ik)
            endif
enddo

return
end
