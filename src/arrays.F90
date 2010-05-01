! -*- F90 -*-

module arrays
  use iso_c_binding
  implicit none

  ! fortran array pointer
  real*8, pointer, save :: cord(:,:,:), vel(:,:,:), stress0(:,:,:,:), &
       force(:,:,:), balance(:,:,:)


#ifdef USE_CUDA

  ! mirror pointer for C
  type(c_ptr), target, save :: pcord, pvel, pstress0, &
       pforce, pbalance

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
    use iso_c_binding
    implicit none

    integer, intent(in) :: nz, nx

#ifndef USE_CUDA

    allocate(cord(nz, nx, 2))
    allocate(vel(nz, nx, 2))
    allocate(stress0(nz, nx, 4, 4))
    allocate(force(nz, nx, 2))
    allocate(balance(nz, nx, 2))

#else


    real*8, pointer :: tmp1d(:)

    ! allocate as 1D array, then reshape the dimension
    call allocate_double(pcord, tmp1d, nz*nx*2)
    call c_f_pointer(pcord, cord, [nz, nx, 2])

    call allocate_double(pvel, tmp1d, nz*nx*2)
    call c_f_pointer(pvel, vel, [nz, nx, 2])

    call allocate_double(pstress0, tmp1d, nz*nx*4*4)
    call c_f_pointer(pstress0, stress0, [nz, nx, 4, 4])

    call allocate_double(pforce, tmp1d, nz*nx*2)
    call c_f_pointer(pforce, force, [nz, nx, 2])

    call allocate_double(pbalance, tmp1d, nz*nx*2)
    call c_f_pointer(pbalance, balance, [nz, nx, 2])

  contains

    subroutine allocate_double(p, x, n)
      use iso_c_binding
      implicit none
      type(c_ptr), target, intent(out) :: p
      real*8, pointer, intent(out) :: x(:)
      integer, intent(in) :: n

      integer(c_size_t) :: m
      real*8 dummy

      m = n * c_sizeof(dummy)
      write(*,*) m, n

      call cudaMallocHost(p, m)
      call c_f_pointer(p, x, [n])
    end subroutine allocate_double


    subroutine allocate_int(p, x, n)
      use iso_c_binding
      implicit none
      type(c_ptr), target, intent(out) :: p
      integer, pointer, intent(out) :: x(:)
      integer, intent(in) :: n

      integer(c_size_t) :: m
      integer dummy

      m = n * c_sizeof(dummy)
      write(*,*) m, n

      call cudaMallocHost(p, m)
      call c_f_pointer(p, x, [n])
    end subroutine allocate_int

#endif

  end subroutine allocate_arrays

end module arrays
