! -*- F90 -*-

module arrays
#ifdef USE_CUDA
  use iso_c_binding
#endif
  implicit none

  ! fortran array pointer
  real*8, pointer, save :: cord(:,:,:), temp(:,:), vel(:,:,:), stress0(:,:,:,:), &
       force(:,:,:), balance(:,:,:), amass(:,:), rmass(:,:), &
       area(:,:,:), dvol(:,:,:), strain(:,:,:), bc(:,:,:)

  integer, pointer, save :: ncod(:,:,:)

  ! temporary array
  real*8, pointer, save :: junk2(:,:)

#ifdef USE_CUDA

  ! mirror pointer for C
  type(c_ptr), target, save :: pcord, ptemp, pvel, pstress0, &
       pforce, pbalance, pamass, prmass, &
       parea, pdvol, pstrain, pbc, &
       pncod


  !-----------------------------------
  ! functions from external library
  !-----------------------------------
  interface
     subroutine cudaMallocHost(p, n) bind(c, name='cudaMallocHost')
       use iso_c_binding
       implicit none
       type(c_ptr) :: p
       integer(c_size_t), value :: n
     end subroutine cudaMallocHost
  end interface

  interface
     subroutine cudaFreeHost(p) bind(c, name='cudaFreeHost')
       use iso_c_binding
       implicit none
       type(c_ptr), value :: p
     end subroutine cudaFreeHost
  end interface

#endif


contains

  subroutine allocate_arrays(nz, nx)
#ifdef USE_CUDA
    use iso_c_binding
#endif
    implicit none

    integer, intent(in) :: nz, nx

#ifndef USE_CUDA

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

#else


    real*8, pointer :: tmp1d(:)
    integer, pointer :: itmp1d(:)


    ! allocate as 1D array, then reshape the dimension
    call allocate_double(pcord, tmp1d, nz*nx*2)
    call c_f_pointer(pcord, cord, [nz, nx, 2])

    call allocate_double(ptemp, tmp1d, nz*nx)
    call c_f_pointer(ptemp, temp, [nz, nx])

    call allocate_double(pvel, tmp1d, nz*nx*2)
    call c_f_pointer(pvel, vel, [nz, nx, 2])

    call allocate_double(pstress0, tmp1d, (nz-1)*(nx-1)*4*4)
    call c_f_pointer(pstress0, stress0, [nz-1, nx-1, 4, 4])

    call allocate_double(pforce, tmp1d, nz*nx*2)
    call c_f_pointer(pforce, force, [nz, nx, 2])

    call allocate_double(pbalance, tmp1d, nz*nx*2)
    call c_f_pointer(pbalance, balance, [nz, nx, 2])

    call allocate_double(pamass, tmp1d, nz*nx)
    call c_f_pointer(pamass, amass, [nz, nx])

    call allocate_double(prmass, tmp1d, nz*nx)
    call c_f_pointer(prmass, rmass, [nz, nx])

    call allocate_double(parea, tmp1d, (nz-1)*(nx-1)*4)
    call c_f_pointer(parea, area, [nz-1, nx-1, 4])

    call allocate_double(pdvol, tmp1d, (nz-1)*(nx-1)*4)
    call c_f_pointer(pdvol, dvol, [nz-1, nx-1, 4])

    call allocate_double(pstrain, tmp1d, (nz-1)*(nx-1)*3)
    call c_f_pointer(pstrain, strain, [nz-1, nx-1, 3])

    call allocate_double(pbc, tmp1d, nz*nx*2)
    call c_f_pointer(pbc, bc, [nz, nx, 2])


    call allocate_int(pncod, itmp1d, nz*nx*2)
    call c_f_pointer(pncod, ncod, [nz, nx, 2])

    allocate(junk2(nz, nx))

  contains

    subroutine allocate_double(p, x, n)
      use iso_c_binding
      implicit none
      type(c_ptr), target, intent(out) :: p
      real*8, pointer, intent(out) :: x(:)
      integer, intent(in) :: n

      integer(c_size_t) :: m
      real*8 dummy

      m = n * sizeof(dummy)

      call cudaMallocHost(p, m)
      call c_f_pointer(p, x, [n])
      !write(*,*) m, n, p
    end subroutine allocate_double


    subroutine allocate_int(p, x, n)
      use iso_c_binding
      implicit none
      type(c_ptr), target, intent(out) :: p
      integer, pointer, intent(out) :: x(:)
      integer, intent(in) :: n

      integer(c_size_t) :: m
      integer dummy

      m = n * sizeof(dummy)

      call cudaMallocHost(p, m)
      call c_f_pointer(p, x, [n])
      !write(*,*) m, n, p
    end subroutine allocate_int

#endif

  end subroutine allocate_arrays

end module arrays
