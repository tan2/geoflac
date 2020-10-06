!-*- F90 -*-

! --------- Flac ------------------------- 

subroutine flac

use arrays
use params
use nvtx_mod
include 'precision.inc' 

! Update Thermal State
! Skip the therm calculations if itherm = 3
call nvtxStartRange('fl_therm')
call fl_therm
call nvtxEndRange()

if (itherm .eq.2) goto 500  ! Thermal calculation only

! Calculation of strain rates from velocity
call nvtxStartRange('fl_srate')
call fl_srate
!$ACC update device(dtavg,nsrate)
call nvtxEndRange()

! Changing marker phases
! XXX: change_phase is slow, don't call it every loop
call nvtxStartRange('change_phase')
if( mod(nloop, ifreq_rmasses).eq.0 ) call change_phase
call nvtxEndRange()

! Update stresses by constitutive law (and mix isotropic stresses)
call nvtxStartRange('fl_rheol')
call fl_rheol
call nvtxEndRange()

! update stress boundary conditions
call nvtxStartRange('bc_update')
if (nystressbc.eq.1) call bc_update
call nvtxEndRange()

! Calculations in a node: forces, balance, velocities, new coordinates
call nvtxStartRange('fl_node')
call fl_node
call nvtxEndRange()

! New coordinates
call nvtxStartRange('fl_move')
call fl_move
call nvtxEndRange()

! Adjust real masses due to temperature
call nvtxStartRange('rmasses')
if( mod(nloop,ifreq_rmasses).eq.0 ) call rmasses
call nvtxEndRange()

! Adjust inertial masses or time step due to deformations
call nvtxStartRange('dt_mass')
if( mod(nloop,ifreq_imasses) .eq. 0 ) call dt_mass
call nvtxEndRange()

500 continue


return
end subroutine flac
