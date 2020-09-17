subroutine outtracer
USE marker_data
use arrays
use params

include 'precision.inc'
parameter( kindr=4 )
real xtrak(nmtracers),ytrak(nmtracers),temptrak(nmtracers),phtrak(nmtracers)
real prestrak(nmtracers),straintrak(nmtracers)
real(kindr) D1d(nmtracers)

integer :: i, id, j, k, kk, n, nn, nrec, nwords

nrec = 0
! define record number and write it to contents
if( lastout .eq. 1 ) then
    nrec = 1
    open (1,file='_tracers.0')
else
    open (1,file='_tracers.0',status='old',err=5)
    do while (.TRUE.)
        read( 1, *, end=10 ) nrec,nmtracers
    end do
    5 continue
    open (1,file='_tracers.0',position='append')
    nrec = 0
    10 continue
    nrec = nrec + 1
    backspace(1)
endif
write( 1, '(i6,1x,i8,1x,i8,1x,f10.6)' ) nrec, nmtracers,nloop,  time/sec_year/1.d6
close(1)

! Coordinates  [km]
nwords = nmtracers
call bar2euler

do kk = 1,nmtracers
    id = idtracer(kk)

    if(mark_dead(id) .eq. 0) then
        xtrak(kk) = 0.d0
        ytrak(kk) = 0.d0
        temptrak(kk) = 0.d0
        prestrak(kk) = 0.d0
        straintrak(kk) = 0.d0
        phtrak(kk) = 0.d0
    else
        n = mark_ntriag(id)
        nn = (n-1)/2
        k = mod(n-1, 2) + 1
        j = mod(nn, nz-1) + 1
        i = nn/(nz-1) + 1

        xtrak(kk) = real(mark_x(id))
        ytrak(kk) = real(mark_y(id))
        tmpr = real(0.25d0*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1)))
        temptrak(kk) = real(tmpr)
        prestrak(kk) = real(stressI(j,i))
        straintrak(kk) = real(strainII(j,i))
        phtrak(kk) = int(mark_phase(id))
    endif
enddo

D1d = 0.d0
do i = 1, nmtracers
D1d(i) = xtrak(i)
enddo
open (1,file='outtrackxx.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) D1d 
close (1)
do i = 1, nmtracers
D1d(i) = ytrak(i)
enddo
open (1,file='outtrackyy.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) D1d 
close (1)
do i = 1, nmtracers
D1d(i) = temptrak(i)
enddo
open (1,file='outtracktemp.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) D1d 
close (1)
do i = 1, nmtracers
D1d(i) = prestrak(i)
enddo
open (1,file='outtrackpres.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) D1d 
close (1)
do i = 1, nmtracers
D1d(i) = straintrak(i)
enddo
open (1,file='outtrackstrain.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) D1d 
close (1)
do i = 1, nmtracers
D1d(i) = phtrak(i)
enddo
open (1,file='outtrackphase.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) D1d
close (1)

return
end
