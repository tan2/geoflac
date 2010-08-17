!-*- F90 -*-

! --------- Flac ------------------------- 

subroutine flac

use arrays
include 'precision.inc' 
include 'params.inc'
include 'arrays.inc'

#ifndef USE_CUDA

! Update Thermal State
! Skip the therm calculations if itherm = 3
call fl_therm

if (itherm .eq.2) goto 500  ! Thermal calculation only

! Calculation of strain rates from velocity
call fl_srate

! Changing marker phases
call change_phase

! Update stresses by constitutive law (and mix isotropic stresses)
call fl_rheol

! update stress boundary conditions
if (ynstressbc.eq.1.) call bc_update

! Calculations in a node: forces, balance, velocities, new coordinates
call fl_node

! New coordinates
call fl_move

! Adjust real masses due to temperature
if( mod(nloop,ifreq_rmasses).eq.0 ) call rmasses

! Adjust inertial masses or time step due to deformations
if( mod(nloop,ifreq_imasses) .eq. 0 ) call dt_mass

! Adjust time Step 
call dt_adjust

500 continue

#else

interface
   subroutine cu_flac(pforce, pbalance, pvel, &
     pcord, pstress0, ptemp, &
     prmass, pamass, &
     parea, pdvol, pstrain, &
     boff, &
     pbc, pncod, &
     time, time_t, &
     dtmax_therm, dt, &
     nloop, itherm, movegrid, &
     ifreq_rmasses, ifreq_imasses, &
     nx, nz) bind(c)
     use iso_c_binding
     implicit none
     type(c_ptr), target, value :: pforce, pbalance, pvel, pcord, pstress0, ptemp, &
          prmass, pamass, parea, pdvol, pstrain, pbc, pncod

     real*8 boff, time, time_t, dtmax_therm, dt
     integer nloop, itherm, movegrid, ifreq_rmasses, ifreq_imasses, nx, nz
   end subroutine cu_flac
end interface

call cu_flac(pforce, pbalance, pvel, &
     pcord, pstress0, ptemp, &
     prmass, pamass, &
     parea, pdvol, pstrain, &
     boff, &
     pbc, pncod, &
     time, time_t, &
     dtmax_therm, dt, &
     nloop, itherm, movegrid, &
     ifreq_rmasses, ifreq_imasses, &
     nx, nz)

#endif


#if 0

! write the first three iterations of arrays to file
! the file is used to compared f90 and cuda results
irec = nloop - nloop_restarted + 1
open(111, file='force.rr',action='write',access='direct',recl=8*nz*nx*2)
write(111, rec=irec) force
close(111)
open(111, file='balance.rr',action='write',access='direct',recl=8*nz*nx*2)
write(111, rec=irec) balance
close(111)
open(111, file='vel.rr',action='write',access='direct',recl=8*nz*nx*2)
write(111, rec=irec) vel
close(111)
open(111, file='stress0.rr',action='write',access='direct',recl=8*(nz-1)*(nx-1)*4*4)
write(111, rec=irec) stress0
close(111)
open(111, file='cord.rr',action='write',access='direct',recl=8*nz*nx*2)
write(111, rec=irec) cord
close(111)
open(111, file='rmass.rr',action='write',access='direct',recl=8*nz*nx)
write(111, rec=irec) rmass
close(111)
open(111, file='area.rr',action='write',access='direct',recl=8*(nz-1)*(nx-1)*4)
write(111, rec=irec) area
close(111)
open(111, file='dvol.rr',action='write',access='direct',recl=8*(nz-1)*(nx-1)*4)
write(111, rec=irec) dvol
close(111)
open(111, file='strain.rr',action='write',access='direct',recl=8*(nz-1)*(nx-1)*3)
write(111, rec=irec) strain
close(111)
write(*,*) 'nloop:', nloop, irec

if(irec .ge. 1) stop 111

#endif


return
end subroutine flac
