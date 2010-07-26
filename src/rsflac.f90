
! Re-starting FLAC

subroutine rsflac
use arrays
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


parameter( kindr=8, kindi=4 )

real(kindr), allocatable :: dum1(:),dum2(:,:)
integer(kindi), allocatable :: dum11(:), idum2(:,:)
real*8 rtime, rdt

! Try to open 'restart.rec'. 
! If file does not exist - restart from last record in '_contents.rs'
! If exists - read nrec from it.
open(1,file='restart.rec', status='old', err=10)
read(1,*) nrec
close(1)
goto 20

10 continue
open( 1, file='_contents.rs', status='old' )
do while (.TRUE.)
    read( 1, *, end=30 ) nrec, nloop, time_my,nmarkers,nmtracers
!    write(*,*) nrec, nloop,time_my,nmarkers,nmtracers
end do
30 close(1)
goto 40

20 continue
open( 1, file='_contents.rs', status='old' )
do while (.TRUE.)
    read( 1, *, end=60 ) nrecf, nloop, time_my
    if( nrecf .eq. nrec ) goto 50
end do
60 call SysMsg('RESTART: could not find record number given in RESTART.REC in _CONTENTS.RS')
stop

50 endfile(1)
goto 30


40 continue

! Read time and dt
open (1,file='time.rs',access='direct',recl=2*8) 
read (1,rec=nrec) rtime, rdt
close (1)
time = rtime
dt = rdt
time_t = time

! Coordinates and velocities
nwords = nz*nx*2

open (1,file='cord.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) cord
close (1)
open (1,file='vel.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) vel
close (1)

! Strain
nwords = 3*(nz-1)*(nx-1)

open (1,file='strain.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) strain
close (1)


! Stress
nwords = 4*4*(nx-1)*(nz-1)

open (1,file='stress.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) stress0
close (1)


! 2-D (nx*nz) arrays - nodes defined
nwords = nz*nx

! Temperature
open (1,file='temp.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) temp
close (1)


! 2-D (nx-1)*(nz-1) arrays - elements defined
allocate( dum2(nz-1,nx-1) )

nwords = (nz-1)*(nx-1)

! Phases
allocate( idum2(nz-1,nx-1) )
open (1,file='phase.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) idum2
close (1)
iphase(1:nz-1,1:nx-1) = idum2(1:nz-1,1:nx-1)
deallocate( idum2 )

! Check if viscous rheology present
ivis_present = 0
do i = 1,nx-1
    do j = 1, nz-1
        iph = iphase(j,i)
        if( irheol(iph).eq.3 .or. irheol(iph).ge.11 ) ivis_present = 1
    end do
end do

! Plastic strain
open (1,file='aps.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
aps(1:nz-1,1:nx-1) = dum2(1:nz-1,1:nx-1)

! Heat sources
open (1,file='source.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
source(1:nz-1,1:nx-1) = dum2(1:nz-1,1:nx-1)

deallocate( dum2 )
nphasl =11
lphase(1) = 2
lphase(2) = 3
lphase(3) = 4
lphase(4) = 6
lphase(5) = 7
lphase(6) = 8
lphase(7) = 9
lphase(8) = 10 
lphase(9) = 11 
lphase(10) = 12 
lphase(11) = 14 



if (iint_marker.eq.1) then


! Markers
nwords= nmarkers
allocate (dum1(nmarkers))
! Markers
open (1,file='xmarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%x = dum1(i)
enddo


open (1,file='ymarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%y = dum1(i)
enddo


open (1,file='xa1marker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%a1 = dum1(i)
enddo


open (1,file='xa2marker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%a2 = dum1(i)
enddo


open (1,file='xmapsmarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%maps = dum1(i)
enddo


open (1,file='xmeIImarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%meII = dum1(i)
enddo


open (1,file='xmpresmarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%mpres = dum1(i)
enddo


open (1,file='xmtempmarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%mtemp = dum1(i)
enddo
deallocate(dum1)


allocate(dum11(nmarkers))

open (1,file='xIDmarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark(i)%ID = dum11(i)
enddo


open (1,file='xntriagmarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark(i)%ntriag = dum11(i)
!write(*,*) dum11(i),mark(i)%ntriag
enddo

open (1,file='xphasemarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark(i)%phase = dum11(i)
enddo

deallocate(dum11)

call marker2elem

endif

! Pressure at the bottom: pisos 
if( nyhydro .eq. 2 ) then
    open(1,file='pisos.rs')
    read(1,*) pisos
    close (1)
endif

! Calculate AREAS (Important: iphase is needed to calculate area!)
call init_areas

! Distribution of REAL masses to nodes
call rmasses

! Boundary conditions
call init_bc

! Inertial masses and time steps (elastic, maxwell and max_thermal)
call dt_mass

return
end
