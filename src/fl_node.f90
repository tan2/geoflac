
!  Calculations of forces from stresses
subroutine fl_node
use arrays
use params
include 'precision.inc'

integer, parameter :: kk(3, 4) = reshape((/ 2,3,4, 1,2,3, 1,2,4, 1,3,4 /), (/3, 4/))
integer, parameter :: mm(3, 4) = reshape((/ 3,3,2, 2,2,2, 3,1,3, 1,1,1 /), (/3, 4/))

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

drat = dt / dt_elastic
if (drat .lt. 1.d0) drat = 1.d0

!$OMP parallel private(i, j, k, fx, fy, factor, &
!$OMP                  p_est, rosubg, &
!$OMP                  press_norm_l, dlx_l, dly_l, &
!$OMP                  press_norm_r, dlx_r, dly_r, &
!$OMP                  iunknown, rho_water_g, water_depth)
!
!$OMP do
!$ACC parallel loop collapse(2) async(1)
do i = 1,nx
    do j = 1,nz
        if(nystressbc.eq.0) then
           force(j,i,1) = 0
           force(j,i,2) = 0
        endif
        ! REGULAR PART - forces from stresses
        
        ! Element (j-1,i-1). Triangles B,C,D
        if ( j.ne.1 .and. i.ne.1 ) then
            do k = 1, 3
                factor = 0.25d0 / area(j-1,i-1,kk(k,1))
                fx = factor * (stress0(j-1,i-1,1,kk(k,1)) * shpdx(j-1,i-1,mm(k,1),kk(k,1)) + &
                               stress0(j-1,i-1,3,kk(k,1)) * shpdz(j-1,i-1,mm(k,1),kk(k,1)))
                fy = factor * (stress0(j-1,i-1,3,kk(k,1)) * shpdx(j-1,i-1,mm(k,1),kk(k,1)) + &
                               stress0(j-1,i-1,2,kk(k,1)) * shpdz(j-1,i-1,mm(k,1),kk(k,1)))
                force(j,i,1) = force(j,i,1) - fx
                force(j,i,2) = force(j,i,2) - fy
            enddo
        endif

        ! Element (j-1,i). Triangles A,B,C.
        if ( j.ne.1 .and. i.ne.nx ) then
            do k = 1, 3
                factor = 0.25d0 / area(j-1,i,kk(k,2))
                fx = factor * (stress0(j-1,i,1,kk(k,2)) * shpdx(j-1,i,mm(k,2),kk(k,2)) + &
                               stress0(j-1,i,3,kk(k,2)) * shpdz(j-1,i,mm(k,2),kk(k,2)))
                fy = factor * (stress0(j-1,i,3,kk(k,2)) * shpdx(j-1,i,mm(k,2),kk(k,2)) + &
                               stress0(j-1,i,2,kk(k,2)) * shpdz(j-1,i,mm(k,2),kk(k,2)))
                force(j,i,1) = force(j,i,1) - fx
                force(j,i,2) = force(j,i,2) - fy
            enddo
        endif

        ! Element (j,i-1). Triangles A,B,D
        if ( j.ne.nz .and. i.ne.1 ) then
            do k = 1, 3
                factor = 0.25d0 / area(j,i-1,kk(k,3))
                fx = factor * (stress0(j,i-1,1,kk(k,3)) * shpdx(j,i-1,mm(k,3),kk(k,3)) + &
                               stress0(j,i-1,3,kk(k,3)) * shpdz(j,i-1,mm(k,3),kk(k,3)))
                fy = factor * (stress0(j,i-1,3,kk(k,3)) * shpdx(j,i-1,mm(k,3),kk(k,3)) + &
                               stress0(j,i-1,2,kk(k,3)) * shpdz(j,i-1,mm(k,3),kk(k,3)))
                force(j,i,1) = force(j,i,1) - fx
                force(j,i,2) = force(j,i,2) - fy
            enddo
        endif

        ! Element (j,i). Triangles A,C,D
        if ( j.ne.nz .and. i.ne.nx ) then
            do k = 1, 3
                factor = 0.25d0 / area(j,i,kk(k,4))
                fx = factor * (stress0(j,i,1,kk(k,4)) * shpdx(j,i,mm(k,4),kk(k,4)) + &
                               stress0(j,i,3,kk(k,4)) * shpdz(j,i,mm(k,4),kk(k,4)))
                fy = factor * (stress0(j,i,3,kk(k,4)) * shpdx(j,i,mm(k,4),kk(k,4)) + &
                               stress0(j,i,2,kk(k,4)) * shpdz(j,i,mm(k,4),kk(k,4)))
                force(j,i,1) = force(j,i,1) - fx
                force(j,i,2) = force(j,i,2) - fy
            enddo
        endif

        ! GRAVITY FORCE
        force(j,i,2) = force(j,i,2) - rmass(j,i)*g

  enddo
enddo
!$OMP end do

! BOUNDARY CONDITIONS
if(nyhydro.gt.0) then
    !$OMP do
    !$ACC parallel loop async(1)
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
        endif
    enddo
    !$OMP end do

    !$OMP do
    !$ACC parallel loop async(1)
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

        !write(*,*) i,pisos,force(nz,i,1),force(nz,i,2),press_norm_l,press_norm_r,dlx_l,dlx_r,dly_l,dly_r

    enddo
    !$OMP end do
endif

!$OMP do
!$ACC parallel loop collapse(2) async(1)
do i=1,nx
    do j=1,nz

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
    end do
end do
!$OMP end do
!$OMP end parallel
! Prestress to form the topo when density differences are present WITHOUT PUSHING OR PULLING!
if (i_prestress.eq.1.and.time.lt.600.d3*sec_year) then
     !$ACC parallel loop collapse(2) async(1)
     do k = 1,2
        do i = 1, nx
            vel(nz,i,k) = 0
        enddo
        do j = 1, nz
            vel(j,1,k) = 0
            vel(j,nx,k) = 0
        enddo
    enddo
endif
return
end subroutine fl_node
