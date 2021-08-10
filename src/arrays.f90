! -*- F90 -*-

module arrays
  implicit none

  save

  ! fortran array allocatable
  double precision, allocatable:: cord(:,:,:), temp(:,:), vel(:,:,:), stress0(:,:,:,:), &
       force(:,:,:), amass(:,:), rmass(:,:), &
       area(:,:,:), dvol(:,:,:), strain(:,:,:), bc(:,:,:)

  integer, allocatable :: ncod(:,:,:), iphase(:,:), &
      nopbou(:,:), ncodbou(:,:), jmoho(:)

  double precision, allocatable :: phase_ratio(:,:,:), &
      dtopo(:), dhacc(:), extrusion(:), &
      fmelt(:,:), fmagma(:,:), extr_acc(:), &
      strainr(:,:,:,:), flux(:,:,:,:), &
      aps(:,:),visn(:,:),e2sr(:,:), &
      temp0(:,:),source(:,:),shrheat(:,:), &
      bcstress(:,:)

  double precision, allocatable :: se2sr(:,:), sshrheat(:,:)

  ! remeshing arrays
  double precision, allocatable :: pt(:,:,:), barcord(:,:,:), &
            cold(:,:,:), cnew(:,:,:), &
            cordo(:,:,:), dhnew(:), extnew(:)
  integer, allocatable :: numtr(:,:)

  ! temporary array
  double precision, allocatable :: dummyn(:,:), dummye(:,:), stmpn(:)
  integer, allocatable :: itmp(:,:)

contains

  subroutine allocate_arrays(nz, nx, nphase)
    implicit none

    integer, intent(in) :: nz, nx, nphase
    integer, parameter :: ntri = 4    ! triangles per element

    allocate(cord(nz, nx, 2))
    allocate(temp(nz, nx))
    allocate(vel(nz, nx, 2))
    allocate(stress0(nz-1, nx-1, 4, ntri))
    allocate(force(nz, nx, 2))
    allocate(amass(nz, nx))
    allocate(rmass(nz, nx))
    allocate(area(nz-1, nx-1, ntri))
    allocate(dvol(nz-1, nx-1, ntri))
    allocate(strain(nz-1, nx-1, 3))
    allocate(bc(nz, nx, 2))

    allocate(ncod(nz, nx, 2))
    allocate(iphase(nz-1, nx-1))
    allocate(nopbou(((nz-1)+(nx-1))*2, 4))
    allocate(ncodbou(((nz-1)+(nx-1))*2, 3))
    allocate(jmoho(nx-1))

    allocate(phase_ratio(nphase, nz-1, nx-1))
    allocate(dtopo(nx))
    allocate(dhacc(nx-1))
    allocate(extrusion(nx-1))
    allocate(fmelt(nz-1, nx-1))
    allocate(fmagma(nz-1, nx-1))
    allocate(extr_acc(nx-1))
    allocate(strainr(3, ntri, nz-1, nx-1))
    allocate(flux(2, 2, nz-1, nx-1))
    allocate(aps(nz-1, nx-1))
    allocate(visn(nz-1, nx-1))
    allocate(e2sr(nz-1, nx-1))
    allocate(temp0(nz, nx))
    allocate(source(nz-1, nx-1))
    allocate(shrheat(nz-1, nx-1))
    allocate(bcstress(((nz-1)+(nx-1))*2, 3))

    allocate(se2sr(nz-1,nx-1))
    allocate(sshrheat(nz-1,nx-1))

    allocate(pt((nz-1)*(nx-1)*2, 2, 3), barcord(nz, nx, 3), &
             cold(nz, nx, 2), cnew(nz, nx, 2))
    allocate(numtr(nz, nx))
    allocate(cordo(nz,nx,2))
    allocate(dhnew(nx-1), extnew(nx-1))

    ! tmp arrays used in subroutines
    allocate(dummyn(nz, nx))
    allocate(dummye(nz-1, nx-1))
    allocate(stmpn(max(nx,nz)))
    allocate(itmp(nz, nx))

  end subroutine allocate_arrays

end module arrays
