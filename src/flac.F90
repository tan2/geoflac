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
if( time-time_t .gt. dtmax_therm/10) call fl_therm

if (itherm .eq.2) goto 500  ! Thermal calculation only

! Calculation of strain rates from velocity
call fl_srate

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
     pbc, pncod, &
     time, time_t, &
     dtmax_therm, dt, &
     nloop, itherm, &
     ifreq_rmasses, ifreq_imasses, &
     nx, nz) bind(c)
     use iso_c_binding
     implicit none
     type(c_ptr), value :: pforce, pbalance, pvel, pcord, pstress0, ptemp, &
          prmass, pamass, pbc, pncod

     real*8 time, time_t, dtmax_therm, dt
     integer nloop, itherm, ifreq_rmasses, ifreq_imasses, nx, nz
   end subroutine cu_flac
end interface

call cu_flac(pforce, pbalance, pvel, &
     pcord, pstress0, ptemp, &
     prmass, pamass, &
     pbc, pncod, &
     time, time_t, &
     dtmax_therm, dt, &
     nloop, itherm, &
     ifreq_rmasses, ifreq_imasses, &
     nx, nz)

#endif

return
end subroutine flac
