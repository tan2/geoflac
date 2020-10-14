
!  Distribution of real masses in nodes

subroutine rmasses
!$ACC routine(Eff_dens) seq
use arrays
use params
implicit none

double precision, parameter :: c1d12 = 1.d0/12.d0
integer :: i, j, iblk, jblk
double precision :: dens, Eff_dens

!   Calcualtion of the TRUE GRAVITATIONAL ZONE MASSES
!-----------------------------------
! THE area(n,it) is inverse of "real" DOUBLE area (=1./det) =>
! area (n,it) ( in program) = 1./(2 real_area)
! real_area = 0.5* (1./area(n,t))
!-----------------------------------

!$ACC kernels async(1)
rmass = 0
!$ACC end kernels

do iblk = 1, 2
    do jblk = 1, 2
        !$OMP parallel do private(dens)
        !$ACC parallel loop collapse(2) async(1)
        do i = iblk, nx-1, 2
            do j = jblk, nz-1, 2

                !  Area and densities of zones
                dens = Eff_dens( j,i )

                ! Distribution 1/3 of the mass of each element to the nodes
                ! *0.5 - becuase 1/area is double of real area; *0.5 - 2 grids

                ! (1) Element A:
                rmass(j  ,i  )=rmass(j  ,i  )+c1d12/area(j,i,1)*dens
                rmass(j+1,i  )=rmass(j+1,i  )+c1d12/area(j,i,1)*dens
                rmass(j  ,i+1)=rmass(j  ,i+1)+c1d12/area(j,i,1)*dens
                ! (2) Element B:
                rmass(j+1,i+1)=rmass(j+1,i+1)+c1d12/area(j,i,2)*dens
                rmass(j+1,i  )=rmass(j+1,i  )+c1d12/area(j,i,2)*dens
                rmass(j  ,i+1)=rmass(j  ,i+1)+c1d12/area(j,i,2)*dens

                ! (3) Element C:
                rmass(j  ,i  )=rmass(j  ,i  )+c1d12/area(j,i,3)*dens
                rmass(j+1,i  )=rmass(j+1,i  )+c1d12/area(j,i,3)*dens
                rmass(j+1,i+1)=rmass(j+1,i+1)+c1d12/area(j,i,3)*dens

                ! (4) Element D:
                rmass(j  ,i  )=rmass(j  ,i  )+c1d12/area(j,i,4)*dens
                rmass(j+1,i+1)=rmass(j+1,i+1)+c1d12/area(j,i,4)*dens
                rmass(j  ,i+1)=rmass(j  ,i+1)+c1d12/area(j,i,4)*dens
            enddo
        enddo
        !$OMP end parallel do
    enddo
enddo

return
end
