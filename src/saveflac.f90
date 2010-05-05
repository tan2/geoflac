!---------------------------------------------------------------
!      Saving state
!---------------------------------------------------------------

subroutine saveflac
use arrays
USE marker_data
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

parameter( kindr=4, kindi=4 )
real(kindr), allocatable :: dum1(:),dum2(:,:),dum3(:,:,:),dum4(:,:,:,:)
integer(kindi), allocatable :: dum11(:)
real*8 rtime, rdt

! define record number and write it to contents

if( lastsave .eq. 1 ) then
    nrec = 1
    open (1,file='_contents.rs')
else
    open (1,file='_contents.rs',status='old',err=5)
    do while (.TRUE.)
        read( 1, *, end=10 ) nrec
    end do
    5 continue
    open (1,file='_contents.rs',position='append')
    nrec = 0
    10 continue
    nrec = nrec + 1
endif
write( 1, '(i4,1x,i8,1x,f6.2,1x,i9,1x,i9)' ) nrec, nloop, time/sec_year/1.e6, &
     nmarkers,nmtracers
close(1)

! Time and dt
open (1,file='time.rs',access='direct',recl=2*8) 
rtime = time
rdt = dt
write (1,rec=nrec) rtime, rdt
close (1) 


! Coordinates and velocities
allocate( dum3(nz,nx,2) )

nwords = nz*nx*2

dum3(1:nz,1:nx,1:2) = cord(1:nz,1:nx,1:2)
open (1,file='cord.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum3
close (1)

dum3(1:nz,1:nx,1:2) = vel(1:nz,1:nx,1:2)
open (1,file='vel.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum3
close (1)

deallocate( dum3 )


! Strain
allocate( dum3(3,nz-1,nx-1) )

nwords = 3*(nz-1)*(nx-1)

dum3(1:3,1:nz-1,1:nx-1) = strain(1:nz-1,1:nx-1,1:3)
open (1,file='strain.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum3
close (1)

deallocate( dum3 )


! Stress
allocate( dum4(4,4,nz-1,nx-1) )

nwords = 4*4*(nz-1)*(nx-1)
dum4(1:4,1:4,1:nz-1,1:nx-1) = stress0(1:nz-1,1:nx-1,1:4,1:4)
open (1,file='stress.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum4
close (1)

deallocate( dum4 )


! 2-D (nx*nz) arrays - nodes defined
allocate( dum2(nz,nx) )

nwords = nz*nx

! Temperature
dum2(1:nz,1:nx) = temp(1:nz,1:nx)
open (1,file='temp.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum2
close (1)

deallocate( dum2 )


! 2-D (nx-1)*(nz-1) arrays - elements defined
allocate( dum2(nz-1,nx-1) )

nwords = (nz-1)*(nx-1)

! Phases
dum2(1:nz-1,1:nx-1) = phasez(1:nz-1,1:nx-1)
open (1,file='phasez.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum2
close (1)

! Plastic strain
dum2(1:nz-1,1:nx-1) = aps(1:nz-1,1:nx-1)
open (1,file='aps.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum2
close (1)


! Heat sources
dum2(1:nz-1,1:nx-1) = source(1:nz-1,1:nx-1)
open (1,file='source.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum2
close (1)
deallocate(dum2)

if(iint_marker.eq.1) then

allocate(dum1(nmarkers))
nwords= nmarkers
! Markers
do i = 1,nmarkers
dum1(i) = mark(i)%x
enddo
open (1,file='xmarker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%y
enddo
open (1,file='ymarker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%a1
enddo
open (1,file='xa1marker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%a2
enddo
open (1,file='xa2marker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%maps
enddo
open (1,file='xmapsmarker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%meII
enddo
open (1,file='xmeIImarker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%mpres
enddo
open (1,file='xmpresmarker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%mtemp
enddo
open (1,file='xmtempmarker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
deallocate(dum1)

allocate(dum11(nmarkers))
do i = 1,nmarkers
dum11(i) = mark(i)%ID
enddo
open (1,file='xIDmarker.rs',access='direct',recl=nwords*kindi) 
write (1,rec=nrec) dum11 
close (1)
do i = 1,nmarkers
dum11(i) = mark(i)%ntriag
!write(*,*) mark(i)%ntriag,dum11(i)
enddo
open (1,file='xntriagmarker.rs',access='direct',recl=nwords*kindi) 
write (1,rec=nrec) dum11 
close (1)
do i = 1,nmarkers
dum11(i) = mark(i)%phase
enddo
open (1,file='xphasemarker.rs',access='direct',recl=nwords*kindi) 
write (1,rec=nrec) dum11 
close (1)
deallocate(dum11)

endif

return 
end
