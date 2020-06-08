
subroutine init_visc
  use arrays
  use params
  use matprops
  implicit none

  integer :: i, j
  !$ACC parallel loop collapse(2)
  !$OMP parallel do
  do i = 1,nx-1  
      do j = 1,nz-1
          visn(j,i) = Eff_visc(j,i)
      end do
  end do
  !$OMP end parallel do
  !$ACC end parallel
  return
end subroutine init_visc
