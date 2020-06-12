
!  Calculations of forces from stresses
subroutine fl_node
use arrays
use params
implicit none

!  1 - 3
!  |   |
!  2 - 4
!
!  diagonal / :
!
!   A:        B:
!
!  1---3         1
!  | /         / |
!  2         2---3
!
!  diagonal \ :
!
!   C:        D:
!
!  1          1---3
!  | \         \  |
!  2---3          2
!
!    assemblage of forces is COUNTRE CLOCK-WISE !
!

integer :: i, j, k
double precision :: drat, fx, fy, &
                  p_est, rosubg, &
                  press_norm_l, dlx_l, dly_l, &
                  press_norm_r, dlx_r, dly_r, &
                  rho_water_g, water_depth

drat = dt / dt_elastic
if (drat .lt. 1.d0) drat = 1.d0

!!**$ACC parallel private(i, j, fx, fy, &
!!**$ACC                  p_est, rosubg, &
!!**$ACC                  press_norm_l, dlx_l, dly_l, &
!!**$ACC                  press_norm_r, dlx_r, dly_r, &
!!**$ACC                  rho_water_g, water_depth)
!!**$ACC loop collapse(2)
!$OMP parallel private(i, j, fx, fy, &
!$OMP                  p_est, rosubg, &
!$OMP                  press_norm_l, dlx_l, dly_l, &
!$OMP                  press_norm_r, dlx_r, dly_r, &
!$OMP                  rho_water_g, water_depth)
!$OMP do
!$ACC parallel loop collapse(2) private(fx, fy)
do i = 1,nx
    do j = 1,nz
        if(nystressbc.eq.0) then
           force(j,i,1) = 0
           force(j,i,2) = 0
           balance(j,i,1)=0
           balance(j,i,2)=0
        endif
        ! REGULAR PART - forces from stresses
        
        ! Element (j-1,i-1). Triangles B,C,D
        if ( j.ne.1 .and. i.ne.1 ) then
            ! triangle B
            ! side 2-3
            fx = stress0(j-1,i-1,1,2) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,3,2) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            fy = stress0(j-1,i-1,3,2) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,2,2) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 3-1
            fx = stress0(j-1,i-1,1,2) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,3,2) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i-1,3,2) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,2,2) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

            ! triangle C
            ! side 2-3
            fx = stress0(j-1,i-1,1,3) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,3,3) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            fy = stress0(j-1,i-1,3,3) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,2,3) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 3-1
            fx = stress0(j-1,i-1,1,3) * (cord(j-1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,3,3) * (cord(j-1,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i-1,3,3) * (cord(j-1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,2,3) * (cord(j-1,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

            ! triangle D
            ! side 1-2
            fx = stress0(j-1,i-1,1,4) * (cord(j  ,i  ,2)-cord(j-1,i-1,2)) - &
                 stress0(j-1,i-1,3,4) * (cord(j  ,i  ,1)-cord(j-1,i-1,1))
            fy = stress0(j-1,i-1,3,4) * (cord(j  ,i  ,2)-cord(j-1,i-1,2)) - &
                 stress0(j-1,i-1,2,4) * (cord(j  ,i  ,1)-cord(j-1,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i-1,1,4) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,3,4) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i-1,3,4) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,2,4) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

        endif

        ! Element (j-1,i). Triangles A,B,C.
        if ( j.ne.1 .and. i.ne.nx ) then
            ! triangle A
            ! side 1-2
            fx = stress0(j-1,i  ,1,1) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,3,1) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            fy = stress0(j-1,i  ,3,1) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,2,1) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i  ,1,1) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,3,1) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i  ,3,1) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,2,1) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

            ! triangle B
            ! side 1-2
            fx = stress0(j-1,i  ,1,2) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 stress0(j-1,i  ,3,2) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            fy = stress0(j-1,i  ,3,2) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 stress0(j-1,i  ,2,2) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i  ,1,2) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,3,2) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i  ,3,2) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,2,2) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

            ! triangle C
            ! side 1-2
            fx = stress0(j-1,i  ,1,3) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,3,3) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            fy = stress0(j-1,i  ,3,3) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,2,3) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i  ,1,3) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,3,3) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i  ,3,3) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,2,3) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

        endif
        
        ! Element (j,i-1). Triangles A,B,D
        if ( j.ne.nz .and. i.ne.1 ) then
            ! triangle A
            ! side 2-3
            fx = stress0(j  ,i-1,1,1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 stress0(j  ,i-1,3,1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            fy = stress0(j  ,i-1,3,1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 stress0(j  ,i-1,2,1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i-1,1,1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,3,1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i-1,3,1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,2,1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

            ! triangle B
            ! side 1-2
            fx = stress0(j  ,i-1,1,2) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,3,2) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i-1,3,2) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,2,2) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i-1,1,2) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,3,2) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            fy = stress0(j  ,i-1,3,2) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,2,2) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

            ! triangle D
            ! side 2-3
            fx = stress0(j  ,i-1,1,4) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,3,4) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            fy = stress0(j  ,i-1,3,4) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,2,4) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i-1,1,4) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,3,4) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i-1,3,4) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,2,4) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

        endif

        ! Element (j,i). Triangles A,C,D
        if ( j.ne.nz .and. i.ne.nx ) then
            ! triangle A
            ! side 1-2
            fx = stress0(j  ,i  ,1,1) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,3,1) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i  ,3,1) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,2,1) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i  ,1,1) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,3,1) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            fy = stress0(j  ,i  ,3,1) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,2,1) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

            ! triangle C
            ! side 1-2
            fx = stress0(j  ,i  ,1,3) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,3,3) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i  ,3,3) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,2,3) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i  ,1,3) * (cord(j  ,i  ,2)-cord(j+1,i+1,2)) - &
                 stress0(j  ,i  ,3,3) * (cord(j  ,i  ,1)-cord(j+1,i+1,1))
            fy = stress0(j  ,i  ,3,3) * (cord(j  ,i  ,2)-cord(j+1,i+1,2)) - &
                 stress0(j  ,i  ,2,3) * (cord(j  ,i  ,1)-cord(j+1,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

            ! triangle D
            ! side 1-2
            fx = stress0(j  ,i  ,1,4) * (cord(j+1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,3,4) * (cord(j+1,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i  ,3,4) * (cord(j+1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,2,4) * (cord(j+1,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i  ,1,4) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,3,4) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            fy = stress0(j  ,i  ,3,4) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,2,4) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25d0*fx
            force(j,i,2) = force(j,i,2) - 0.25d0*fy
            balance(j,i,1) = balance(j,i,1) + 0.25d0*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25d0*abs(fy)

        endif

        ! GRAVITY FORCE
        force(j,i,2) = force(j,i,2) - rmass(j,i)*g
        balance(j,i,2) = balance(j,i,2) + abs(rmass(j,i)*g)

  enddo
enddo
!$OMP end do
!$ACC end parallel

! BOUNDARY CONDITIONS
if(nyhydro.gt.0) then
    !$ACC parallel loop private(rho_water_g, water_depth, press_norm_l, &
    !$ACC                       dlx_l, dly_l, dlx_r, dly_r, press_norm_l, press_norm_r)
    !$OMP do
    do i=1,nx

        ! pressure from water sea on top
        rho_water_g = 1030.d0 * g
        if(i.lt.nx) then
            water_depth = 0.5d0*(cord(1,i+1,2)+cord(1,i,2))
        else
            water_depth = 0.5d0*(cord(1,i-1,2)+cord(1,i,2))
        endif

        if (water_depth.lt.0.d0) then ! No water (above sea level)
            if(i.eq.1) then
                press_norm_l = 0
                dlx_l = 0
                dly_l = 0
                press_norm_r = rho_water_g*((cord(1,i+1,2)+cord(1,i,2))/2.d0)
                dlx_r = cord(1,i+1,1)-cord(1,i  ,1)
                dly_r = cord(1,i+1,2)-cord(1,i  ,2)
            elseif(i.eq.nx) then
                press_norm_l = rho_water_g*((cord(1,i-1,2)+cord(1,i,2))/2.d0)
                dlx_l = cord(1,i  ,1)-cord(1,i-1,1)
                dly_l = cord(1,i  ,2)-cord(1,i-1,2)
                press_norm_r = 0
                dlx_r = 0
                dly_r = 0
            else
                press_norm_l = rho_water_g*((cord(1,i-1,2)+cord(1,i,2))/2.d0)
                dlx_l = cord(1,i  ,1)-cord(1,i-1,1)
                dly_l = cord(1,i  ,2)-cord(1,i-1,2)
                press_norm_r = rho_water_g*((cord(1,i+1,2)+cord(1,i,2))/2.d0)
                dlx_r = cord(1,i+1,1)-cord(1,i  ,1)
                dly_r = cord(1,i+1,2)-cord(1,i  ,2)
            endif
            force(1,i,1) = force(1,i,1)-0.5d0*press_norm_l*dly_l-0.5d0*press_norm_r*dly_r
            force(1,i,2) = force(1,i,2)+0.5d0*press_norm_l*dlx_l+0.5d0*press_norm_r*dlx_r
            balance(1,i,1) = 1.0d+17
        endif
    enddo
    !$OMP end do
    !$ACC end parallel

    !$ACC parallel loop private(p_est, rosubg, &
    !$ACC                       dlx_l, dly_l, dlx_r, dly_r, press_norm_l, press_norm_r)
    !$OMP do
    do i=1,nx

        ! bottom support - Archimed force (normal to the surface, shear component = 0)
        p_est = pisos + 0.5d0*(den(iphsub)+drosub)*g*(cord(nz,i,2)-rzbo)
        rosubg = g * (den(iphsub)+drosub) * (1-alfa(iphsub)*temp(nz,i)+beta(iphsub)*p_est)

        if(i.eq.1) then
            press_norm_l = 0
            dlx_l = 0
            dly_l = 0

            press_norm_r = pisos-rosubg*((cord(nz,i+1,2)+cord(nz,i,2))/2.d0-rzbo)
            dlx_r = cord(nz,i+1,1)-cord(nz,i  ,1)
            dly_r = cord(nz,i+1,2)-cord(nz,i  ,2)
        elseif(i.eq.nx) then
            press_norm_l = pisos-rosubg*((cord(nz,i-1,2)+cord(nz,i,2))/2.d0-rzbo)
            dlx_l = cord(nz,i  ,1)-cord(nz,i-1,1)
            dly_l = cord(nz,i  ,2)-cord(nz,i-1,2)

            press_norm_r = 0
            dlx_r = 0
            dly_r = 0
        else
            press_norm_l = pisos-rosubg*((cord(nz,i-1,2)+cord(nz,i,2))/2.d0-rzbo)
            dlx_l = cord(nz,i  ,1)-cord(nz,i-1,1)
            dly_l = cord(nz,i  ,2)-cord(nz,i-1,2)

            press_norm_r = pisos-rosubg*((cord(nz,i+1,2)+cord(nz,i,2))/2.d0-rzbo)
            dlx_r = cord(nz,i+1,1)-cord(nz,i  ,1)
            dly_r = cord(nz,i+1,2)-cord(nz,i  ,2)
        endif
            
        force(nz,i,1) = force(nz,i,1)-0.5d0*press_norm_l*dly_l-0.5d0*press_norm_r*dly_r
        force(nz,i,2) = force(nz,i,2)+0.5d0*press_norm_l*dlx_l+0.5d0*press_norm_r*dlx_r

        balance(nz,i,1) = 1.0d+17
        !write(*,*) i,pisos,force(nz,i,1),force(nz,i,2),press_norm_l,press_norm_r,dlx_l,dlx_r,dly_l,dly_r

    enddo
    !$OMP end do
    !$ACC end parallel
endif

boff = 0

!XXX:  reduction(max:boff)
!$ACC data copyin(drat)
!$ACC parallel loop collapse(2)
!$OMP do reduction(max:boff)
do i=1,nx
    do j=1,nz

        ! BALANCE-OFF
        if( iand(ncod(j,i,1),1).eq.1 .or. j.le.n_boff_cutoff ) then
            balance(j,i,1) = 0
        else
            balance(j,i,1) = abs(force(j,i,1)) / (balance(j,i,1) + 1.0d-9)
        endif

        if( iand(ncod(j,i,2),2).eq.2 .or. j.le.n_boff_cutoff ) then
            balance(j,i,2) = 0
        else
            balance(j,i,2) = abs(force(j,i,2)) / (balance(j,i,2) + 1.0d-9)
        endif

        ! DAMPING
        if( iand(ncod(j,i,1),1).ne.1 .and. abs(vel(j,i,1)).gt.1.0d-13 ) then
            force(j,i,1) = force(j,i,1) - demf*sign(force(j,i,1),vel(j,i,1))
        endif

        if( iand(ncod(j,i,2),2).ne.2 .and. abs(vel(j,i,2)).gt.1.0d-13 ) then
            force(j,i,2) = force(j,i,2) - demf*sign(force(j,i,2),vel(j,i,2))
        endif

        ! VELOCITIES FROM FORCES
        if( ncod(j,i,1) .eq. 1 ) then
            vel(j,i,1) = bc(j,i,1) 
        else
            vel(j,i,1) = vel(j,i,1) + dt*force(j,i,1)/(amass(j,i)*drat*drat)
        endif
        if( ncod(j,i,2) .eq. 1 ) then
            vel(j,i,2) = bc(j,i,2)
        else
            vel(j,i,2) = vel(j,i,2) + dt*force(j,i,2)/(amass(j,i)*drat*drat)
        endif
        ! MAX balance-off
        boff = max(boff,balance(j,i,1))
        boff = max(boff,balance(j,i,2))

    end do
end do
!$OMP end do
!$OMP end parallel
!$ACC end parallel
!$ACC end data

! Prestress to form the topo when density differences are present WITHOUT PUSHING OR PULLING!
if (i_prestress.eq.1.and.time.lt.600.d3*sec_year) then
     !$ACC parallel loop
     do k = 1,2
        do i = 1, nx
            vel(nz,i,k) = 0
        enddo
        do j = 1, nz
            vel(j,1,k) = 0
            vel(j,nx,k) = 0
        enddo
    enddo
    !$ACC end parallel loop
endif
return
end subroutine fl_node
