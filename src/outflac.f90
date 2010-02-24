!---------------------------------------------------------------
!      Output
!---------------------------------------------------------------

subroutine outflac
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

parameter( kindr=4, kindi=2 )

real(kindr), allocatable :: D1d(:),De(:,:),Dn2(:,:,:)
real(kindr) rtime

! define record number and write it to contents
if( lastout .eq. 1 ) then
    nrec = 1
    open (1,file='_contents.0')
else
    open (1,file='_contents.0',status='old',err=5)
    do while (.TRUE.)
        read( 1, *, end=10 ) nrec
    end do
    5 continue
    open (1,file='_contents.0')
    nrec = 0
    10 continue
    nrec = nrec + 1
    backspace(1)
endif
write( 1, '(i4,1x,i8,1x,f6.2)' ) nrec, nloop, time/sec_year/1.e6
close(1)

! Time
open (1,file='time.0',access='direct',recl=kindr)
rtime = time
write (1,rec=nrec) rtime
close (1) 


! Coordinates in [km]
allocate( Dn2(nz,nx,2) )

nwords = nz*nx*2
Dn2(1:nz,1:nx,1:2) = cord(1:nz,1:nx,1:2) / 1000
open (1,file='mesh.0',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) Dn2
close (1)

deallocate( Dn2 )


! 2-D (nx-1)*(nz-1) arrays - elements defined
allocate( De(nz-1,nx-1) )

nwords = (nz-1)*(nx-1)

! Velocities in [cm/year]
if( io_vel.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = 0.25*( vel(j,i,1)+vel(j+1,i,1)+vel(j,i+1,1)+vel(j+1,i+1,1) ) * sec_year*100
!            De(j,i) = 0.25*( force(j,i,1)+force(j+1,i,1)+force(j,i+1,1)+force(j+1,i+1,1) )
        end do
    end do
    open (1,file='vx.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)

    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = 0.25*( vel(j,i,2)+vel(j+1,i,2)+vel(j,i+1,2)+vel(j+1,i+1,2) ) * sec_year*100
!            De(j,i) = 0.25*( force(j,i,2)+force(j+1,i,2)+force(j,i+1,2)+force(j+1,i+1,2) )
        end do
    end do
    open (1,file='vz.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Strain rate II
if( io_srII.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            if( e2sr(j,i).ne.0. ) then
                De(j,i) = dlog10( e2sr(j,i) )
            else
                De(j,i) = 0
            endif
        enddo
    enddo
    open (1,file='srII.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Strain
if( io_eII.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = strainII(i,j)
        end do
    end do
    open (1,file='eII.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif

    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = strainI(i,j)
        end do
    end do
    open (1,file='eI.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)

! Density 
if( io_mark.eq.1 ) then
!    De(1:nz-1,1:nx-1) = rmarker(1:nz-1,1:nx-1)
   do j = 1, nz-1
   do i = 1, nx-1
          De(j,i) = Eff_dens(i,j)
   enddo
   enddo
    open (1,file='phas.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)
endif


! APS
if( io_aps.eq.1 ) then
    De(1:nz-1,1:nx-1) = aps(1:nz-1,1:nx-1)
    open (1,file='aps.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Stress II in [kb]
if( io_sII.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = stressII(i,j) * 1.e-8
        end do
    end do
    open (1,file='sII.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Sxx in [kb]
if( io_sxx.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            sxx = 0.25 * (stress0(1,1,j,i)+stress0(1,2,j,i)+stress0(1,3,j,i)+stress0(1,4,j,i) )
            De(j,i) = ( sxx-stressI(i,j) ) * 1.e-8
        end do
    end do
    open (1,file='sxx.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Szz in [kb]
if( io_szz.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            szz = 0.25 * (stress0(2,1,j,i)+stress0(2,2,j,i)+stress0(2,3,j,i)+stress0(2,4,j,i) )
            De(j,i) = ( szz-stressI(i,j) ) * 1.e-8
        end do
    end do
    open (1,file='szz.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Sxz in [kb]
if( io_sxz.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            sxz = 0.25 * (stress0(3,1,j,i)+stress0(3,2,j,i)+stress0(3,3,j,i)+stress0(3,4,j,i)) * 1.e-8
            De(j,i) = sxz * 1.e-8
        end do
    end do
    open (1,file='sxz.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Pressure in [kb]
if( io_pres.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = stressI(i,j) * 1.e-8
        end do
    end do
    open (1,file='pres.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Temperature
if( io_temp.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = 0.25*( temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1) )
        end do
    end do
    open (1,file='temp.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Melt 
if( io_melt.eq.1 ) then
    De(1:nz-1,1:nx-1) = phasez(1:nz-1,1:nx-1)
    open (1,file='melt.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Viscosities (log)
if( io_visc.eq.1 ) then
    De(1:nz-1,1:nx-1) = dlog10( visn(1:nz-1,1:nx-1) )
    open (1,file='visc.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif



! Heat sources
if( io_src.eq.1 ) then
    De(1:nz-1,1:nx-1) = source(1:nz-1,1:nx-1)
    open (1,file='src.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Energy dissipation
if( io_diss.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            iph = iphase(i,j,phasez(j,i))
            De(j,i) = shrheat(j,i)/den(iph)/hs
        enddo
    enddo
    open (1,file='diss.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif

deallocate( De )


! 1-D nx array - nodes defined
allocate( D1d(nx) )
nwords = nx

! Surface heat flow
if( io_hfl.eq.1 ) then
    do i = 1,nx
        ii = min(i,nx-1)
        dtmpr = temp(2,i) - temp(1,i)
        dl = -(cord(2,i,2)-cord(1,i,2))/1000
        D1d(i) = Eff_conduct(ii,1) * dtmpr/dl
    end do
    open (1,file='hfl.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) D1d
    close (1)
endif


! Topo
if( io_topo.eq.1 ) then
    do i = 1,nx
        D1d(i) = cord(1,i,2)/1000
    end do
    open (1,file='topo.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) D1d
    close (1)
endif

deallocate( D1d )


return 
end
