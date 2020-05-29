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
  real*8, allocatable :: junk2(:,:), xmpt(:,:,:), tkappa(:)

  integer, allocatable :: iphase(:,:), &
      irheol_fl(:,:), &
      nopbou(:,:), ncodbou(:,:), idtracer(:)

  real*8, allocatable :: phase_ratio(:,:,:), &
      dtopo(:), dhacc(:), extrusion(:), &
      andesitic_melt_vol(:), extr_acc(:), &
      strainr(:,:,:,:), &
      aps(:,:),visn(:,:),e2sr(:,:), &
      temp0(:,:),source(:,:),shrheat(:,:), &
      bcstress(:,:)

  double precision, allocatable :: se2sr(:,:), sshrheat(:,:)
  
  ! remeshing arrays
  real*8, allocatable :: pt(:,:,:), barcord(:,:,:), &
            cold(:,:,:), cnew(:,:,:)
  integer, allocatable :: numtr(:,:)


  !$ACC declare create(cord, temp, vel, stress0, force, balance, amass, rmass, &
  !$ACC                area, dvol, strain, bc, ncod, junk2, xmpt, tkappa, &
  !$ACC                iphase, irheol_fl, &
  !$ACC                nopbou, ncodbou, idtracer, phase_ratio, dtopo, dhacc, extrusion, &
  !$ACC                andesitic_melt_vol, extr_acc, strainr, aps, visn, e2sr, &
  !$ACC                temp0, source, shrheat, bcstress, se2sr, sshrheat, &
  !$ACC                pt, barcord, cold, cnew, numtr)
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
    allocate(iphase(nz-1, nx-1))
    allocate(irheol_fl(nz-1, nx-1))
    allocate(nopbou(((nz-1)+(nx-1))*2, 4))
    allocate(ncodbou(((nz-1)+(nx-1))*2, 3))
    allocate(idtracer((nz-1)*(nx-1)*2))

    allocate(phase_ratio(20, nz-1, nx-1))
    allocate(dtopo(nx))
    allocate(dhacc(nx))
    allocate(extrusion(nx-1))
    allocate(andesitic_melt_vol(nx-1))
    allocate(extr_acc(nx-1))
    allocate(strainr(3, 4, nz-1, nx-1))
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

    ! tmp arrays used in subroutines
    allocate(junk2(nz, nx))
    allocate(xmpt(2,3,(nz-1)*(nx-1)*2))
    allocate(tkappa(nx))

    !$ACC update device(cord, temp, vel, stress0, force, balance, amass, rmass, &
    !$ACC               area, dvol, strain, bc, ncod, junk2, xmpt, tkappa, &
    !$ACC               iphase, irheol_fl, &
    !$ACC               nopbou, ncodbou, idtracer, phase_ratio, dtopo, dhacc, extrusion, &
    !$ACC               andesitic_melt_vol, extr_acc, strainr, aps, visn, e2sr, &
    !$ACC               temp0, source, shrheat, bcstress, se2sr, sshrheat, &
    !$ACC               pt, barcord, cold, cnew, numtr)

  end subroutine allocate_arrays

end module arrays
