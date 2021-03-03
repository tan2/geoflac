!-*- F90 -*-

! --------- Flac ------------------------- 

subroutine flac

use arrays
use params
include 'precision.inc' 

! Update Thermal State
! Skip the therm calculations if itherm = 3
call fl_therm

if (itherm .eq.2) goto 500  ! Thermal calculation only

! Calculation of strain rates from velocity
call fl_srate

! Changing marker phases
! XXX: change_phase is slow, don't call it every loop
if( mod(nloop, ifreq_rmasses).eq.0 ) call change_phase

! Update stresses by constitutive law (and mix isotropic stresses)
call fl_rheol

! update stress boundary conditions
if (nystressbc.eq.1) call bc_update

! Calculations in a node: forces, balance, velocities, new coordinates
call fl_node

! New coordinates
call fl_move

! Adjust real masses due to temperature
if( mod(nloop,ifreq_rmasses).eq.0 ) call rmasses

! Adjust inertial masses or time step due to deformations
if( mod(nloop,ifreq_imasses) .eq. 0 ) call dt_mass

500 continue


return
end subroutine flac
