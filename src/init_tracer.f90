subroutine init_tracer
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ nfreemarkers,ndeadmarkers,mpt(mnz*mnx*2,2,3)
parameter( min_elmarkers = 0, max_elmarkers = 12 )
!type(marker) :: mark(nmarkers)

common /tracers/ tracerid(mnz*mnx*2)
kint = 0
nmtracers = 0
iik = 0 
jjk = 0
nzz = 0
do m = 1, nzone_tracer 
   nzz = nzz + (ity2(m)-ity1(m)+1)
enddo
do n = 1, nzone_tracer
   do i = itx1(n),itx2(n) 
         iik = iik + 1
      do j = ity1(n),ity2(n)
         jjk = jjk + 1
         do k = 1,2
         do l = 1,nmarkers
               ntr = 2*( (nz-1)*(i-1) +j-1) + k
               mm = 2*( (nzz-1)*(iik-1)+jjk - 1)
            if (ntr.ne.mark(l)%ntriag) then
               cycle
            else
               kint = kint + 1
               if(kint.eq.1) then
               nmtracers = nmtracers + 1
               tracerid(nmtracers) = mark(l)%ID
!               write(*,*) nmtracers, tracerid(mm+k),mm+k               
               endif
               endif
          enddo
          enddo
          kint = 0
    enddo
enddo
enddo
return
end
