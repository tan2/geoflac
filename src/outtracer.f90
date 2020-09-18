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
        vx = vel(j1,i1,1)*ba1 + vel(j2,i2,1)*ba2 + vel(j3,i3,1)*ba3
        vz = vel(j1,i1,2)*ba1 + vel(j2,i2,2)*ba2 + vel(j3,i3,2)*ba3
        tmpr = temp(j1,i1)*ba1 + temp(j2,i2)*ba2 + temp(j3,i3)*ba3
        ! no interpolation for elemental values
        si = stressI(j,i)
        sxx = 0.25d0*(stress0(j,i,1,1)+stress0(j,i,1,2)+stress0(j,i,1,3)+stress0(j,i,1,4)) - si
        szz = 0.25d0*(stress0(j,i,2,1)+stress0(j,i,2,2)+stress0(j,i,2,3)+stress0(j,i,2,4)) - si
        sxz = 0.25d0*(stress0(j,i,3,1)+stress0(j,i,3,2)+stress0(j,i,3,3)+stress0(j,i,3,4)) - si

        tracerx(kk) = real(x) * 1.e-3
        tracerz(kk) = real(z) * 1.e-3
        tracervx(kk) = real(vx) * 100 * sec_year
        tracervz(kk) = real(vz) * 100 * sec_year
        tracert(kk) = real(tmpr)
        tracerp(kk) = real(si) * (-1e-8)
        tracersxx(kk) = real(sxx) * 1e-8
        tracerszz(kk) = real(szz) * 1e-8
        tracersxz(kk) = real(sxz) * 1e-8
        tracere(kk) = real(strainII(j,i))
        traceredot(kk) = real(srateII(j,i))
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
! Velocity [cm/yr]
open (1,file='tracervx.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracervx
close (1)
open (1,file='tracervz.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracervz
close (1)
! Temperature [Celsius]
open (1,file='tracert.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracert
close (1)
! Pressure [kbar]
open (1,file='tracerp.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracerp
close (1)
! Stress [kbar]
open (1,file='tracersxx.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracersxx
close (1)
open (1,file='tracerszz.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracerszz
close (1)
open (1,file='tracersxz.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracersxz
close (1)
! Strain II
open (1,file='tracere.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracere
close (1)
! Strain rate II
open (1,file='traceredot.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) traceredot
close (1)
! Phase
open (1,file='tracerph.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) tracerph
close (1)

return
end
