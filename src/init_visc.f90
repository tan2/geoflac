
subroutine init_visc
  use arrays
  use params
  include 'precision.inc'

  do i = 1,nx-1  
      do j = 1,nz-1
          visn(j,i) = Eff_visc(j,i)
      end do
  end do
  return
end subroutine init_visc
