! -*- F90 -*-

module arrays
  implicit none

  ! fortran array pointer
  real*8, pointer, save :: cord(:,:,:), temp(:,:), vel(:,:,:), stress0(:,:,:,:), &
       force(:,:,:), balance(:,:,:), amass(:,:), rmass(:,:), &
       area(:,:,:), dvol(:,:,:), strain(:,:,:), bc(:,:,:)

  integer, pointer, save :: ncod(:,:,:)

  ! temporary array
  real*8, pointer, save :: junk2(:,:)

  !!! maximum number of ELEMENTS !!!
  integer, parameter :: mnz=200, mnx=700, max_markers_per_elem=32

  integer :: iphase(mnz,mnx), nphase_counter(20,mnz,mnx), &
      ntopmarker(mnx), itopmarker(max_markers_per_elem,mnx), &
      irheol_fl(mnz,mnx), &
      nopbou((mnz+mnx)*2,4), ncodbou((mnz+mnx)*2,3)

  real*8 :: phase_ratio(20,mnz,mnx), &
      dtopo(mnx+1), dhacc(mnx+1), extrusion(mnx), &
      andesitic_melt_vol(mnx), extr_acc(mnx), &
      strainr(3,4,mnz,mnx), &
      aps(mnz,mnx),visn(mnz,mnx),e2sr(mnz,mnx), &
      temp0(mnz+1,mnx+1),source(mnz,mnx),shrheat(mnz,mnx), &
      bcstress((mnz+mnx)*2,3)

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

    allocate(junk2(nz, nx))

  end subroutine allocate_arrays

end module arrays
