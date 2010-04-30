! -*- F90 -*-

module arrays
  use iso_c_binding
  implicit none

  ! maximum number of ELEMENTS
  integer, parameter :: mnz=100, mnx=800

  ! fortran array pointer
  real*8, pointer, save :: cord(:,:,:)


#ifdef USE_CUDA

  ! mirror pointer for C
  type(c_ptr), target, save :: pcord

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

  subroutine allocate_arrays
    use iso_c_binding
    implicit none

#ifndef USE_CUDA

    allocate(cord(mnz+1, mnx+1, 2))

#else


    real*8, pointer :: tmp1d(:)
    real*8, pointer :: tmp3d(:,:,:)

    ! allocate as 1D array, then reshape the dimension
    call allocate_double(pcord, tmp1d, (mnz+1)*(mnx+1)*2)
    call c_f_pointer(pcord, tmp3d, [mnz+1, mnx+1, 2])


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
