
! Setup some parameters (rmass,amass,initial stress,vel,viscosity)

subroutine setflac
use arrays
use params
implicit none

nloop = 0
time = 0.
!$ACC update device(nloop,time)

! Mesh generator
call init_cord

! Initial accumulated plastic strain
aps = 0

! Initial velocity
vel = 0

dvol = 0
strain = 0

! Phases in the mesh
call init_phase

! Setup markers
if (iint_marker.eq.1) then
call init_marker
endif
! Setup tracers
if (iint_tracer.eq.1) call init_tracer

! Inverse Areas of triangles
call init_areas

! Initiate temperature field
call init_temp

! Calculation of the initial STRESSES (as hydrostatic)
call init_stress

! Setup boundary conditions
call init_bc

! Distribution of REAL masses to nodes
call rmasses

! Initialization of viscosity
if( ivis_present.eq.1 ) call init_visc

! Inertial masses and time steps (elastic and maxwell)
call dt_mass
dt = min( dt_elastic, dt_maxwell )
!$ACC update device(dt)

! Initiate parameters for stress averaging
dtavg=0
nsrate=-1

!Initialization
!$ACC update device(temp, vel, stress0, force, balance, amass, rmass, &
!$ACC               area, dvol, strain, bc, ncod, junk2, xmpt, tkappa, &
!$ACC               iphase, nphase_counter, ntopmarker, itopmarker, irheol_fl, &
!$ACC               nopbou, ncodbou, idtracer, phase_ratio, dtopo, dhacc, extrusion, &
!$ACC               andesitic_melt_vol, extr_acc, strainr, aps, visn, e2sr, &
!$ACC               temp0, source, shrheat, bcstress, &
!$ACC               pt, barcord, cold, cnew, numtr, &
!$ACC               se2sr, sshrheat, dtavg, nsrate)


return
end
