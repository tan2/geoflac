
! Setup some parameters (rmass,amass,initial stress,vel,viscosity)

subroutine setflac
use arrays
use params
use marker_data
implicit none

nloop = 0
time = 0.d0
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
    ! Setup tracers
    if (iint_tracer.eq.1) call init_tracer
endif

! Inverse Areas of triangles
call init_areas

! Initiate temperature field
call init_temp

! Calculation of the initial STRESSES (as hydrostatic)
call init_stress

! Setup boundary conditions
call init_bc

temp0 = temp
shrheat = 0
sshrheat = 0
dtopo = 0
extrusion = 0
andesitic_melt_vol = 0
e2sr = 1d-16
se2sr = 1d-16

!! These arrays are not used below and can be uploaded to GPU early
!$ACC update device(force, &
!$ACC               dvol, strain, bc, ncod, junk2, xmpt, tkappa, &
!$ACC               nopbou, ncodbou, dtopo, dhacc, extrusion, &
!$ACC               andesitic_melt_vol, extr_acc, strainr, aps, &
!$ACC               temp0, source, shrheat, bcstress, &
!$ACC               pt, barcord, cold, cnew, numtr, &
!$ACC               se2sr, sshrheat) async(1)

! Distribution of REAL masses to nodes
!$ACC update device(area, cord, temp, stress0, iphase, phase_ratio)
call rmasses

! Initialization of viscosity
!$ACC update device(e2sr)
if( ivis_present.eq.1 ) call init_visc

! Inertial masses and time steps (elastic and maxwell)
!$ACC update device(vel)
call dt_mass

! Initiate parameters for stress averaging
dtavg=0
nsrate=-1
!$ACC update device(dtavg, nsrate) async(1)

return
end
