subroutine outmarker
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
parameter( min_elmarkers = 0, max_elmarkers = 12 )
parameter( kindr=4, kindi=4 )
!type(marker) :: mark (nmarkers)
real(kindr) D1d(nmarkers)
integer(kindi) D1i(nmarkers)

nrec = 0
D1d = 0.
! define record number and write it to contents
if( lastout .eq. 1 ) then
    nrec = 1
    open (1,file='_outmarkers.0')
else
    open (1,file='_outmarkers.0',status='old',err=5)

    do while (.TRUE.)
        read( 1, *, end=10 ) nrec
    end do
    5 continue
    open (1,file='_outmarkers.0')
    nrec = 0
    10 continue
    nrec = nrec + 1
    backspace(1)
endif
write( 1, '(i4,1x,i8,1x,i8,1x,f6.2)' ) nrec, nloop,nmarkers, time/sec_year/1.e6
close(1)

! Coordinates  [km]
nwords = nmarkers 
!write(*,*) nmarkers,nwords
do i = 1, nmarkers
 D1d(i)= mark(i)%x
! if (i.lt.1000) write(*,*) mark(i)%x, D1d(i)
enddo
open (1,file='outmarkxx.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) D1d 
close (1)

D1d = 0.

do i = 1,nmarkers
 D1d(i)= mark(i)%y
enddo
open (1,file='outmarkyy.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec)  D1d 
close (1)


D1i = 0
do l = 1,nmarkers
 D1i(l)= mark(l)%phase
enddo
open (1,file='outmarkphase.0',access='direct',recl=nwords*kindi)
write (1,rec=nrec) D1i
close (1)

return
end
