
! Calculate thermal field in explict form

subroutine fl_therm
!$ACC routine(Eff_cp) seq
!$ACC routine(Eff_conduct) seq
use arrays
use params
include 'precision.inc'

double precision :: D(3,3)  ! diffusion operator
double precision :: diff_elem(nz-1, nx-1)
double precision :: inv_area

tan_mzone = tan(0.5d0 * angle_mzone * 3.14159265358979323846d0 / 180.d0)
! max. width of the magma zone @ moho (as if melting occurs at 200 km)
ihalfwidth_mzone = ceiling(tan_mzone * 200e3 / dxmin)

! real_area = 0.5* (1./area(n,t))
! Calculate Fluxes in every triangle
!      flux (num_triangle, direction(x,y), j, i)
!
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

! saving old temperature
if (istress_therm > 0 .or. itype_melting == 1) then
    !$ACC kernels async(1)
    temp0(:,:) = temp(:,:)
    !$ACC end kernels
endif

!$OMP Parallel private(i,j,iph,cp_eff,cond_eff,dissip,quad_area, &
!$OMP                  fr_lambda, &
!$OMP                  delta_fmagma,deltaT,real_area13,area_n,rhs, &
!$OMP                  jm,area_ratio,z_moho,z_melt,x_melt,h,x,z,inv_area)

if (itype_melting == 1) then
    ! M: fmegma, magma fraction in the element
    ! dM/dt = P - M * fr_lambda
    ! Rearrange after forward Euler for dM/dt
    ! M(t+dt) = M(t) + P * dt - M(t) * fr_lambda * dt
    !   where P is the production rate and the second term is the freezing rate
    !   fr_lambda is temperature dependent: lower T freezes faster
    ! fr_lambda = lambda_freeze * exp(-lambda_freeze_tdep * (T - Ttop))
    ! P = fmelt * prod_magma * R_zone *  (A_zone / A_melt)
    !   fmelt: melt fraction, calculated in chage_phase()
    !   prod_magma: magma production rate, how rapid a melt migrate away
    !               from the active melting region.
    !   A_melt: the area of the melting element
    !   R_zone, A_zone: the magma distribution ratio and the area of the two
    !                   zones (crust and mantle). R_crust + R_mantle == 1.0
    !   area_ratio = A_zone / A_melt
    !
    ! Magma freezing will release latent heat to the element
    ! rho * cp * deltaT = rho * latent_heat * delta_fmagma
    ! Rearrage to: deltaT = delta_fmagma * latent_heat / cp
    ! This heat is distributed to the 4 corners evenly
    !
    !$OMP do
    !$ACC parallel loop collapse(2) async(1)
    do i = 1,nx-1
        do j = 1,nz-1
            cp_eff = Eff_cp( j,i )
            tmpr = 0.25d0*(temp0(j,i)+temp0(j+1,i)+temp0(j,i+1)+temp0(j+1,i+1))
            fr_lambda = lambda_freeze * exp(-lambda_freeze_tdep * (tmpr-t_top))
            delta_fmagma = min(fmagma(j,i), fmagma(j,i) * dt * fr_lambda)
            fmagma(j,i) = fmagma(j,i) - delta_fmagma

            ! latent heat released by freezing magma
            deltaT = delta_fmagma * latent_heat_magma / cp_eff / 4
            !$OMP atomic update
            !$ACC atomic update
            temp(j  ,i  ) = temp(j  ,i  ) + deltaT
            !$OMP atomic update
            !$ACC atomic update
            temp(j  ,i+1) = temp(j  ,i+1) + deltaT
            !$OMP atomic update
            !$ACC atomic update
            temp(j+1,i  ) = temp(j+1,i  ) + deltaT
            !$OMP atomic update
            !$ACC atomic update
            temp(j+1,i+1) = temp(j+1,i+1) + deltaT
        end do
    enddo
    !$OMP end do

    !$OMP do
    !$ACC parallel loop collapse(2) async(1)
    do i = 1,nx-1
        ! XXX: Assume melting cannot happen above the moho. (j > jmoho(i)) is always true
        ! starting from j=1, not jmoho(i), so that the loop can be collapsed and faster in computation
        do j = 1,nz-1
            jm = jmoho(i)
            quad_area = 0.5d0/area(j,i,1) + 0.5d0/area(j,i,2) ! area of this element
            if (j>jm .and. fmelt(j,i) > 0) then
                ! This element is under melting. The melt will migrate to the
                ! mantle and crust above. Treat the migration as instantaneous.

                ! Within crust, melts migrate by diking, propagate upward vertically
                ! area_ratio: the area of the crust column / the area of the melting element
                !     ~ the thickness of the crust column / the thickness of the melting element
                area_ratio = (cord(1,i,2)+cord(1,i+1,2)-cord(jm,i,2)-cord(jm,i+1,2)) / &
                    (cord(j,i,2)+cord(j,i+1,2)-cord(j+1,i,2)-cord(j+1,i+1,2))
                do jj = 1, jm
                    !$OMP atomic update
                    !$ACC atomic update
                    fmagma(jj,i) = fmagma(jj,i) + fmelt(j,i) * (1d0 - ratio_mantle_mzone) * area_ratio * prod_magma * dt
                enddo

                ! Within mantle, melts migrate by percolation, propagate upward slantly
                z_moho = 0.5d0 * (cord(jm,i,2) + cord(jm,i+1,2))
                z_melt = 0.5d0 * (cord(j+1,i,2) + cord(j+1,i+1,2))
                x_melt = 0.5d0 * (cord(j+1,i,1) + cord(j+1,i+1,1))
                h = z_moho - z_melt
                ! area_ratio: the area of the mantle triangle / the area of the melting element
                area_ratio = h * h * tan_mzone / quad_area
                ! ii: the potential region of magma distribution zone at moho
                do ii = max(1,i-ihalfwidth_mzone), min(nx-1,i+ihalfwidth_mzone)
                    do jj = jmoho(ii)+1, j
                        x = 0.5d0 * (cord(jj,ii,1) + cord(jj,ii+1,1))
                        z = 0.5d0 * (cord(jj,ii,2) + cord(jj,ii+1,2))
                        if (abs(x - x_melt) <= tan_mzone * (z - z_melt)) then
                            !$OMP atomic update
                            !$ACC atomic update
                            fmagma(jj,ii) = fmagma(jj,ii) + fmelt(j,i) * ratio_mantle_mzone * area_ratio * prod_magma * dt
                        endif
                    enddo
                enddo
            endif
            fmagma(j,i) = min(fmagma(j,i), fmagma_max)
        enddo
    enddo
endif

!$ACC parallel loop collapse(2) async(1)
!$OMP do
do i = 1,nx-1
    do j = 1,nz-1

        iph = iphase(j,i)

        ! Calculating effective material properties
        cp_eff = Eff_cp( j,i )
        cond_eff = Eff_conduct( j,i )

        ! if shearh-heating flag is true
        if( ishearh.eq.1 .and. itherm.ne.2 ) then
            dissip = shrheat(j,i)
        else
            dissip = 0
        endif

        ! diffusivity
        diff_elem(j,i) = cond_eff/den(iph)/cp_eff

        ! Additional sources - radiogenic and shear heating
        dummye(j,i) = ( source(j,i) + dissip/den(iph) ) / cp_eff

    end do
end do    
!$OMP end do

!$ACC parallel loop collapse(2) async(1)
!$OMP do
do i = 1,nx
    do j = 1,nz

        rhs = 0.d0
        area_n = 0.d0

        ! Element (j-1,i-1). Triangle B
        if ( j.ne.1 .and. i.ne.1 ) then
            inv_area = 1.d0/area(j-1,i-1,2)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j-1,i-1) * real_area13 - 0.5d0 * diff_elem(j-1,i-1) * inv_area * ( &
                temp(j-1,i  ) * (shpdx(j-1,i-1,1,2)*shpdx(j-1,i-1,3,2) + shpdz(j-1,i-1,1,2)*shpdz(j-1,i-1,3,2)) + &
                temp(j  ,i-1) * (shpdx(j-1,i-1,2,2)*shpdx(j-1,i-1,3,2) + shpdz(j-1,i-1,2,2)*shpdz(j-1,i-1,3,2)) + &
                temp(j  ,i  ) * (shpdx(j-1,i-1,3,2)*shpdx(j-1,i-1,3,2) + shpdz(j-1,i-1,3,2)*shpdz(j-1,i-1,3,2)) )
        endif

        ! Element (j-1,i). Triangles A,B
        if ( j.ne.1 .and. i.ne.nx ) then
            ! triangle A
            inv_area = 1.d0/area(j-1,i,1)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j-1,i) * real_area13 - 0.5d0 * diff_elem(j-1,i) * inv_area * ( &
                temp(j-1,i  ) * (shpdx(j-1,i,1,1)*shpdx(j-1,i,2,1) + shpdz(j-1,i,1,1)*shpdz(j-1,i,2,1)) + &
                temp(j  ,i  ) * (shpdx(j-1,i,2,1)*shpdx(j-1,i,2,1) + shpdz(j-1,i,2,1)*shpdz(j-1,i,2,1)) + &
                temp(j-1,i+1) * (shpdx(j-1,i,3,1)*shpdx(j-1,i,2,1) + shpdz(j-1,i,3,1)*shpdz(j-1,i,2,1)) )

            ! triangle B
            inv_area = 1.d0/area(j-1,i,2)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j-1,i) * real_area13 - 0.5d0 * diff_elem(j-1,i) * inv_area * ( &
                temp(j-1,i+1) * (shpdx(j-1,i,1,2)*shpdx(j-1,i,2,2) + shpdz(j-1,i,1,2)*shpdz(j-1,i,2,2)) + &
                temp(j  ,i  ) * (shpdx(j-1,i,2,2)*shpdx(j-1,i,2,2) + shpdz(j-1,i,2,2)*shpdz(j-1,i,2,2)) + &
                temp(j  ,i+1) * (shpdx(j-1,i,3,2)*shpdx(j-1,i,2,2) + shpdz(j-1,i,3,2)*shpdz(j-1,i,2,2)) )
        endif
        
        ! Element (j,i-1). Triangles A,B
        if ( j.ne.nz .and. i.ne.1 ) then
            ! triangle A
            inv_area = 1.d0/area(j,i-1,1)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j,i-1) * real_area13 - 0.5d0 * diff_elem(j,i-1) * inv_area * ( &
                temp(j  ,i-1) * (shpdx(j,i-1,1,1)*shpdx(j,i-1,3,1) + shpdz(j,i-1,1,1)*shpdz(j,i-1,3,1)) + &
                temp(j+1,i-1) * (shpdx(j,i-1,2,1)*shpdx(j,i-1,3,1) + shpdz(j,i-1,2,1)*shpdz(j,i-1,3,1)) + &
                temp(j  ,i  ) * (shpdx(j,i-1,3,1)*shpdx(j,i-1,3,1) + shpdz(j,i-1,3,1)*shpdz(j,i-1,3,1)) )

            ! triangle B
            inv_area = 1.d0/area(j,i-1,2)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j,i-1) * real_area13 - 0.5d0 * diff_elem(j,i-1) * inv_area * ( &
                temp(j  ,i  ) * (shpdx(j,i-1,1,2)*shpdx(j,i-1,1,2) + shpdz(j,i-1,1,2)*shpdz(j,i-1,1,2)) + &
                temp(j+1,i-1) * (shpdx(j,i-1,2,2)*shpdx(j,i-1,1,2) + shpdz(j,i-1,2,2)*shpdz(j,i-1,1,2)) + &
                temp(j+1,i  ) * (shpdx(j,i-1,3,2)*shpdx(j,i-1,1,2) + shpdz(j,i-1,3,2)*shpdz(j,i-1,1,2)) )
        endif

        ! Element (j,i). Triangle A
        if ( j.ne.nz .and. i.ne.nx ) then
            inv_area = 1.d0/area(j,i,1)
            real_area13 = inv_area / 6.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j,i) * real_area13 - 0.5d0 * diff_elem(j,i) * inv_area * ( &
                temp(j  ,i  ) * (shpdx(j,i,1,1)*shpdx(j,i,1,1) + shpdz(j,i,1,1)*shpdz(j,i,1,1)) + &
                temp(j+1,i  ) * (shpdx(j,i,2,1)*shpdx(j,i,1,1) + shpdz(j,i,2,1)*shpdz(j,i,1,1)) + &
                temp(j  ,i+1) * (shpdx(j,i,3,1)*shpdx(j,i,1,1) + shpdz(j,i,3,1)*shpdz(j,i,1,1)) )
        endif

        ! Update Temperature by Eulerian method 
        temp(j,i) = temp(j,i)+rhs*dt/area_n
    end do
end do
!$OMP end do

! Boundary conditions (top and bottom)
!$ACC parallel loop async(1)
!$OMP do
do i = 1,nx

    temp(1,i) = t_top

    if( itemp_bc.eq.1 ) then
        temp(nz,i) = bot_bc
    elseif( itemp_bc.eq.2 ) then
        cond_eff = Eff_conduct( nz-1, min(i,nx-1) )
        temp(nz,i) = temp(nz-1,i)  +  bot_bc * ( cord(nz-1,i,2)-cord(nz,i,2) ) / cond_eff
    endif

end do
!$OMP end do

! Boundary conditions: dt/dx =0 on left and right  
!$ACC parallel loop async(1)
!$OMP do
do j = 1,nz
    temp(j ,1)  = temp(j,2)
    temp(j, nx) = temp(j,nx-1)
end do
!$OMP end do
!$OMP end parallel

return
end 
