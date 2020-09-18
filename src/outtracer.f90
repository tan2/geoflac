subroutine outtracer
USE marker_data
use arrays
use params

include 'precision.inc'
integer, parameter :: kindr=4
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

do kk = 1,nmtracers
    id = idtracer(kk)

    if(mark_dead(id) .eq. 0) then
        tracerx(kk) = 0.d0
        tracerz(kk) = 0.d0
        tracert(kk) = 0.d0
        tracerp(kk) = 0.d0
        tracere(kk) = 0.d0
        tracerph(kk) = 0
    else
        n = mark_ntriag(id)
        nn = (n-1)/2
        k = mod(n-1, 2) + 1
        j = mod(nn, nz-1) + 1
        i = nn/(nz-1) + 1

        ba1 = mark_a1(kk)
        ba2 = mark_a2(kk)
        ba3 = 1.0d0 - ba1 - ba2

        if (k .eq. 1) then
            i1 = i
            i2 = i
            i3 = i + 1
            j1 = j
            j2 = j + 1
            j3 = j
        else
            i1 = i + 1
            i2 = i
            i3 = i + 1
            j1 = j
            j2 = j + 1
            j3 = j + 1
        endif

        ! interpolate nodal values to the marker
        x = cord(j1,i1,1)*ba1 + cord(j2,i2,1)*ba2 + cord(j3,i3,1)*ba3
        z = cord(j1,i1,2)*ba1 + cord(j2,i2,2)*ba2 + cord(j3,i3,2)*ba3
        tmpr = temp(j1,i1)*ba1 + temp(j2,i2)*ba2 + temp(j3,i3)*ba3

        tracerx(kk) = real(x) * 1.e-3
        tracerz(kk) = real(z) * 1.e-3
        tracert(kk) = real(tmpr)
        tracerp(kk) = real(stressI(j,i)) * (-1e-8)
        tracere(kk) = real(strainII(j,i))
        tracerph(kk) = mark_phase(id)
    endif
enddo

nwords = nmtracers

! Coordinates  [km]
open (1,file='tracerx.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracerx
close (1)
open (1,file='tracerz.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracerz
close (1)
! Temperature [Celsius]
open (1,file='tracert.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracert
close (1)
! Pressure [kbar]
open (1,file='tracerp.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracerp
close (1)
! Strain II
open (1,file='tracere.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracere
close (1)
! Phase
open (1,file='tracerph.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracerph
close (1)

return
end
