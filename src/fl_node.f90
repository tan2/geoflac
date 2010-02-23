
!  Calculations of forces from stresses
subroutine fl_node
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

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

boff = 0

drat = dt / dt_elastic
if (drat .lt. 1.) drat = 1.

!$OMP parallel private(i, j, fx, fy, &
!$OMP                  p_est, rosubg, &
!$OMP                  press_norm_l, dlx_l, dly_l, &
!$OMP                  press_norm_r, dlx_r, dly_r, &
!$OMP                  iunknown, rho_water, water_depth)
!
!$OMP do
do i = 1,nx
    do j = 1,nz
        if(ynstressbc.eq.0.) then
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
            fx = stress0(1,2,j-1,i-1) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(3,2,j-1,i-1) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            fy = stress0(3,2,j-1,i-1) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(2,2,j-1,i-1) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(1,2,j-1,i-1) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(3,2,j-1,i-1) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(3,2,j-1,i-1) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(2,2,j-1,i-1) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle C
            ! side 2-3
            fx = stress0(1,3,j-1,i-1) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(3,3,j-1,i-1) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            fy = stress0(3,3,j-1,i-1) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(2,3,j-1,i-1) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(1,3,j-1,i-1) * (cord(j-1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(3,3,j-1,i-1) * (cord(j-1,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(3,3,j-1,i-1) * (cord(j-1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(2,3,j-1,i-1) * (cord(j-1,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle D
            ! side 1-2
            fx = stress0(1,4,j-1,i-1) * (cord(j  ,i  ,2)-cord(j-1,i-1,2)) - &
                 stress0(3,4,j-1,i-1) * (cord(j  ,i  ,1)-cord(j-1,i-1,1))
            fy = stress0(3,4,j-1,i-1) * (cord(j  ,i  ,2)-cord(j-1,i-1,2)) - &
                 stress0(2,4,j-1,i-1) * (cord(j  ,i  ,1)-cord(j-1,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(1,4,j-1,i-1) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(3,4,j-1,i-1) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(3,4,j-1,i-1) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(2,4,j-1,i-1) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif

        ! Element (j-1,i). Triangles A,B,C.
        if ( j.ne.1 .and. i.ne.nx ) then
            ! triangle A
            ! side 1-2
            fx = stress0(1,1,j-1,i  ) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(3,1,j-1,i  ) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            fy = stress0(3,1,j-1,i  ) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(2,1,j-1,i  ) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(1,1,j-1,i  ) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(3,1,j-1,i  ) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(3,1,j-1,i  ) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(2,1,j-1,i  ) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle B
            ! side 1-2
            fx = stress0(1,2,j-1,i  ) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 stress0(3,2,j-1,i  ) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            fy = stress0(3,2,j-1,i  ) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 stress0(2,2,j-1,i  ) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(1,2,j-1,i  ) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(3,2,j-1,i  ) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(3,2,j-1,i  ) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(2,2,j-1,i  ) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle C
            ! side 1-2
            fx = stress0(1,3,j-1,i  ) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(3,3,j-1,i  ) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            fy = stress0(3,3,j-1,i  ) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(2,3,j-1,i  ) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(1,3,j-1,i  ) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(3,3,j-1,i  ) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(3,3,j-1,i  ) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(2,3,j-1,i  ) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif
        
        ! Element (j,i-1). Triangles A,B,D
        if ( j.ne.nz .and. i.ne.1 ) then
            ! triangle A
            ! side 2-3
            fx = stress0(1,1,j  ,i-1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 stress0(3,1,j  ,i-1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            fy = stress0(3,1,j  ,i-1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 stress0(2,1,j  ,i-1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(1,1,j  ,i-1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(3,1,j  ,i-1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(3,1,j  ,i-1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(2,1,j  ,i-1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle B
            ! side 1-2
            fx = stress0(1,2,j  ,i-1) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(3,2,j  ,i-1) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(3,2,j  ,i-1) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(2,2,j  ,i-1) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(1,2,j  ,i-1) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(3,2,j  ,i-1) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            fy = stress0(3,2,j  ,i-1) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(2,2,j  ,i-1) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle D
            ! side 2-3
            fx = stress0(1,4,j  ,i-1) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(3,4,j  ,i-1) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            fy = stress0(3,4,j  ,i-1) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(2,4,j  ,i-1) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(1,4,j  ,i-1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(3,4,j  ,i-1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(3,4,j  ,i-1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(2,4,j  ,i-1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif

        ! Element (j,i). Triangles A,C,D
        if ( j.ne.nz .and. i.ne.nx ) then
            ! triangle A
            ! side 1-2
            fx = stress0(1,1,j  ,i  ) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(3,1,j  ,i  ) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(3,1,j  ,i  ) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(2,1,j  ,i  ) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(1,1,j  ,i  ) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(3,1,j  ,i  ) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            fy = stress0(3,1,j  ,i  ) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(2,1,j  ,i  ) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle C
            ! side 1-2
            fx = stress0(1,3,j  ,i  ) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(3,3,j  ,i  ) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(3,3,j  ,i  ) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(2,3,j  ,i  ) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(1,3,j  ,i  ) * (cord(j  ,i  ,2)-cord(j+1,i+1,2)) - &
                 stress0(3,3,j  ,i  ) * (cord(j  ,i  ,1)-cord(j+1,i+1,1))
            fy = stress0(3,3,j  ,i  ) * (cord(j  ,i  ,2)-cord(j+1,i+1,2)) - &
                 stress0(2,3,j  ,i  ) * (cord(j  ,i  ,1)-cord(j+1,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle D
            ! side 1-2
            fx = stress0(1,4,j  ,i  ) * (cord(j+1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(3,4,j  ,i  ) * (cord(j+1,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(3,4,j  ,i  ) * (cord(j+1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(2,4,j  ,i  ) * (cord(j+1,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(1,4,j  ,i  ) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(3,4,j  ,i  ) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            fy = stress0(3,4,j  ,i  ) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(2,4,j  ,i  ) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif

        ! GRAVITY FORCE
        force(j,i,2) = force(j,i,2) - rmass(j,i)*g
        balance(j,i,2) = balance(j,i,2) + abs(rmass(j,i)*g)

        ! BOUNDARY CONDITIONS
        ! pressure from water sea on top
        if (nyhydro.gt.0.and. j.eq.1)  then
            rho_water = -10300.
            water_depth = 0.5*(cord(1,i+1,2)+cord(1,i,2))
            if (water_depth.lt.0.) then ! No water (above sea level)
                if(i.eq.1) then
                  press_norm_l = 0
                  dlx_l = 0
                  dly_l = 0
                  press_norm_r = rho_water*((cord(1,i+1,2)+cord(1,i,2))/2.)
                  dlx_r = cord(1,i+1,1)-cord(1,i  ,1)
                  dly_r = cord(1,i+1,2)-cord(1,i  ,2)
                elseif(i.eq.nx) then
                  press_norm_l = rho_water*((cord(1,i-1,2)+cord(1,i,2))/2.)
                  dlx_l = cord(1,i  ,1)-cord(1,i-1,1)
                  dly_l = cord(1,i  ,2)-cord(1,i-1,2)
                  press_norm_r = 0
                  dlx_r = 0
                  dly_r = 0
                else
                  press_norm_l = rho_water*((cord(1,i-1,2)+cord(1,i,2))/2.)
                  dlx_l = cord(1,i  ,1)-cord(1,i-1,1)
                  dly_l = cord(1,i  ,2)-cord(1,i-1,2)
                  press_norm_r = rho_water*((cord(1,i+1,2)+cord(1,i,2))/2.)
                  dlx_r = cord(1,i+1,1)-cord(1,i  ,1)
                  dly_r = cord(1,i+1,2)-cord(1,i  ,2)
                endif
                force(1,i,1) = force(1,i,1)-0.5*press_norm_l*dly_l-0.5*press_norm_r*dly_r
                force(1,i,2) = force(1,i,2)+0.5*press_norm_l*dlx_l+0.5*press_norm_r*dlx_r
                balance(1,i,1) = 1.e+17
            endif
        endif
            
  enddo
enddo
!$OMP end do
if(nyhydro.gt.0) then 
!$OMP do
    do i=1,nx

        ! BOUNDARY CONDITIONS

        ! bottom support - Archimed force (normal to the surface, shear component = 0)
        !write(*,*) force(nz,i,1),force(nz,i,2)            
        p_est = pisos + 0.5*(den(iphsub)+drosub)*g*(cord(nz,i,2)-rzbo)
        rosubg = g * (den(iphsub)+drosub) * (1-alfa(iphsub)*temp(nz,i)+beta(iphsub)*p_est)

        if(i.eq.1) then
            press_norm_l = 0
            dlx_l = 0
            dly_l = 0

            press_norm_r = pisos-rosubg*((cord(nz,i+1,2)+cord(nz,i,2))/2-rzbo)
            dlx_r = cord(nz,i+1,1)-cord(nz,i  ,1)
            dly_r = cord(nz,i+1,2)-cord(nz,i  ,2)
        elseif(i.eq.nx) then
            press_norm_l = pisos-rosubg*((cord(nz,i-1,2)+cord(nz,i,2))/2-rzbo)
            dlx_l = cord(nz,i  ,1)-cord(nz,i-1,1)
            dly_l = cord(nz,i  ,2)-cord(nz,i-1,2)

            press_norm_r = 0
            dlx_r = 0
            dly_r = 0
        else
            press_norm_l = pisos-rosubg*((cord(nz,i-1,2)+cord(nz,i,2))/2-rzbo)
            dlx_l = cord(nz,i  ,1)-cord(nz,i-1,1)
            dly_l = cord(nz,i  ,2)-cord(nz,i-1,2)

            press_norm_r = pisos-rosubg*((cord(nz,i+1,2)+cord(nz,i,2))/2-rzbo)
            dlx_r = cord(nz,i+1,1)-cord(nz,i  ,1)
            dly_r = cord(nz,i+1,2)-cord(nz,i  ,2)
        endif
            
        force(nz,i,1) = force(nz,i,1)-0.5*press_norm_l*dly_l-0.5*press_norm_r*dly_r
        force(nz,i,2) = force(nz,i,2)+0.5*press_norm_l*dlx_l+0.5*press_norm_r*dlx_r

        balance(nz,i,1) = 1.e+17
        !write(*,*) i,pisos,force(nz,i,1),force(nz,i,2),press_norm_l,press_norm_r,dlx_l,dlx_r,dly_l,dly_r

    enddo
!$OMP end do
endif
!$OMP do reduction(max:boff)
do i=1,nx
    do j=1,nz
        
        ! BALANCE-OFF
        if( iand(ncod(j,i,1),1).eq.1 .or. j.le.n_boff_cutoff ) then
            balance(j,i,1) = 0
        else
            balance(j,i,1) = abs(force(j,i,1)) / (balance(j,i,1) + 1.e-9)
        endif

        if( iand(ncod(j,i,2),2).eq.2 .or. j.le.n_boff_cutoff ) then
            balance(j,i,2) = 0
        else
            balance(j,i,2) = abs(force(j,i,2)) / (balance(j,i,2) + 1.e-9)
        endif

        ! DAMPING
        if( iand(ncod(j,i,1),1).ne.1 .and. abs(vel(j,i,1)).gt.1.e-13 ) then
            force(j,i,1) = force(j,i,1) - demf*sign(force(j,i,1),vel(j,i,1))
        endif

        if( iand(ncod(j,i,2),2).ne.2 .and. abs(vel(j,i,2)).gt.1.e-13 ) then
            force(j,i,2) = force(j,i,2) - demf*sign(force(j,i,2),vel(j,i,2))
        endif

        ! VELOCITIES FROM FORCES
        iunknown = 0
        if( ncod(j,i,1) .eq. 1 ) then
            vel(j,i,1) = bc(j,i,1) 
!            vel(j,i,1) = 0.0

!        write(*,*) i,j,vel(j,i,1)
        else
            vel(j,i,1) = vel(j,i,1) + dt*force(j,i,1)/(amass(j,i)*drat*drat)
        endif
        if( ncod(j,i,2) .eq. 1 ) then
            vel(j,i,2) = bc(j,i,2)
            if(iunknown.eq.1) then
                vel(j,i,2) = bc(j,i,2)* sin(time*3.14159/(2*sec_year))
                write(*,*) bc(j,i,2), sin(time*3.14159/(2*sec_year))
            endif
!        write(*,*) i,j,vel(j,i,2)
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
! Prestress to form the topo when density differences are present WITHOUT PUSHING OR PULLING!
if (i_prestress.eq.1.and.(time/sec_year).lt.200000.) then
    do k = 1,2
        do i = 1, nx
            vel(nz,i,k) = 0.
        enddo
        do j = 1, nz
            vel(j,1,k) = 0.
            vel(j,nx,k) = 0.
        enddo
    enddo
endif
return
end subroutine fl_node
