! -*- F90 -*-

module arrays
  implicit none

  save

  ! fortran array allocatable
  real*8, allocatable:: cord(:,:,:), temp(:,:), vel(:,:,:), stress0(:,:,:,:), &
       force(:,:,:), balance(:,:,:), amass(:,:), rmass(:,:), &
       area(:,:,:), dvol(:,:,:), strain(:,:,:), bc(:,:,:)

  integer, allocatable:: ncod(:,:,:)

  ! temporary array
  real*8, allocatable :: junk2(:,:)

  !!! maximum number of ELEMENTS !!!
  integer, parameter :: mnz=200, mnx=700, max_markers_per_elem=32

  integer, allocatable :: iphase(:,:), nphase_counter(:,:,:), &
      ntopmarker(:), itopmarker(:,:), &
      irheol_fl(:,:), &
      nopbou(:,:), ncodbou(:,:)

  real*8, allocatable :: phase_ratio(:,:,:), &
      dtopo(:), dhacc(:), extrusion(:), &
      andesitic_melt_vol(:), extr_acc(:), &
      strainr(:,:,:,:), &
      aps(:,:),visn(:,:),e2sr(:,:), &
      temp0(:,:),source(:,:),shrheat(:,:), &
      bcstress(:,:)

  !$acc declare create(cord, temp, vel, stress0, force, balance, amass, rmass, area, dvol, strain, bc, ncod, junk2)
contains

  subroutine allocate_arrays(nz, nx)
    implicit none

    integer, intent(in) :: nz, nx

    allocate(cord(nz, nx, 2))
    allocate(temp(nz, nx))
    allocate(vel(nz, nx, 2))
    allocate(stress0(nz-1, nx-1, 4, 4))
    allocate(force(nz, nx, 2))
    allocate(balance(nz, nx, 2))
    allocate(amass(nz, nx))
    allocate(rmass(nz, nx))
    allocate(area(nz-1, nx-1, 4))
    allocate(dvol(nz-1, nx-1, 4))
    allocate(strain(nz-1, nx-1, 3))
    allocate(bc(nz, nx, 2))

    allocate(ncod(nz, nx, 2))
    allocate(iphase(nz, nx))
    allocate(nphase_counter(20, nz, nx))
    allocate(ntopmarker(nx))
    allocate(itopmarker(max_markers_per_elem, nx))
    allocate(irheol_fl(nz, nx))
    allocate(nopbou((nz+nx)*2, 4))
    allocate(ncodbou((nz+nx)*2, 3))

    allocate(phase_ratio(20, nz, nx))
    allocate(dtopo(nx+1))
    allocate(dhacc(nx+1))
    allocate(extrusion(nx))
    allocate(andesitic_melt_vol(nx))
    allocate(extr_acc(nx))
    allocate(strainr(3, 4, nz, nx))
    allocate(aps(nz, nx))
    allocate(visn(nz, nx))
    allocate(e2sr(nz, nx))
    allocate(temp0(nz+1, nx+1))
    allocate(source(nz, nx))
    allocate(shrheat(nz, nx))
    allocate(bcstress((nz+nx)*2, 3))

    allocate(junk2(nz, nx))

  end subroutine allocate_arrays

end module arrays
