subroutine outtracer
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
parameter( min_elmarkers = 0, max_elmarkers = 12 )
parameter( kindr=4 )
!type(marker) :: mark(nmarkers)
real xik(nmtracers),timtrk(nmtracers),xtrak(nmtracers),ytrak(nmtracers),temptrak(nmtracers)
real prestrak(nmtracers),straintrak(nmtracers)
real(kindr) D1d(nmtracers)
common /tracers/ tracerid(mnz*mnx*2)

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
endif
write(*,*) nmtracers
write( 1, '(i4,1x,i8,1x,i8,1x,f6.2)' ) nrec, nmtracers,nloop,  time/sec_year/1.e6
close(1)
! Coordinates  [km]
nwords = nmtracers
  call bar2euler
write(*,*) nmtracers,nmarkers
do i = 1,nmtracers
do j = 1,nmarkers
if (tracerid(i).eq.mark(j)%ID) then
     n = mark(j)%ntriag
     nn = (n-1)/2
! if triangle number is pair or impair
     k = mod(n-1, 2) + 1
     l = mod(nn, nz-1) + 1
     ik = nn/(nz-1) + 1

     mark(j)%meII = strainII(l,ik)
     mark(j)%mpres = stressI(l,ik)
     tmpr = 0.25*(temp(l,ik)+temp(l+1,ik)+temp(l,ik+1)+temp(l+1,ik+1))
     mark(j)%mtemp = tmpr
!write(*,*) i,k,n,ik,l, mark(j)%meII,mark(j)%mpres,mark(j)%mtemp
xik(i) = float (i)
timtrk(i) = time/sec_year/1.e6
xtrak(i) = mark(j)%x
ytrak(i) = mark(j)%y
temptrak(i) = mark(j)%mtemp
prestrak(i) = mark(j)%mpres
straintrak(i) = mark(j)%meII
endif
enddo
enddo
D1d = 0.
do i = 1, nmtracers
D1d(i) = xik(i)
enddo
open (1,file='outtrackID.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) D1d 
close (1)
do i = 1, nmtracers
D1d(i) = timtrk(i)
enddo
open (1,file='outtracktime.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) D1d 
close (1)
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

return
end
