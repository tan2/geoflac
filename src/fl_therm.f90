
! Calculate thermal field in explict form

subroutine fl_therm
!$ACC routine(Eff_cp) seq
!$ACC routine(Eff_conduct) seq
use arrays
use params
include 'precision.inc'

double precision, parameter :: heat_latent_magma = 4.2d5  ! J/kg, latent heat of freezing magma
double precision :: D(3,3)  ! diffusion operator

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
if (istress_therm.gt.0) then
    !$ACC kernels async(1)
    temp0(:,:) = temp(:,:)
    !$ACC end kernels
endif

!$OMP Parallel private(i,j,iph,cp_eff,cond_eff,dissip,diff,quad_area, &
!$OMP                  x1,x2,x3,x4,y1,y2,y3,y4,t1,t2,t3,t4,tmpr, &
!$OMP                  ihw,delta_fmagma,qs,real_area13,area_n,rhs)

if (itype_melting == 1) then
    ! M: fmegma, magma fraction in the element
    ! dM/dt = P - M lambda
    !   where P is the production rate and the second term is the freezing rate
    !
    ! Magma production rate:
    ! M(t+dt) = M(t) + P dt
    !   where P is the magma production rate.
    ! P = fmelt * prod_magma * R_zone *  A_zone / A_melt
    !
    !$OMP do
    !$ACC parallel loop collapse(2) async(1)
    do i = 1,nx-1
        ! XXX: Assume melting cannot happen above the moho. (j > jmoho(i)) is always true
        do j = 1,nz-1
            jm = jmoho(i)
            quad_area = 0.5d0/area(j,i,1) + 0.5d0/area(j,i,2) ! area of melting element
            if (j>jm+1 .and. fmelt(j,i) > 0) then
                ! Within crust, melts migrate by diking, propagate upward vertically
                ! area_ratio: the area of the crust column / the area of the melting element
                !     ~ the thickness of the crust column / the thickness of the melting element
                area_ratio = (cord(1,i,2)+cord(1,i+1,2)-cord(jm+1,i,2)-cord(jm+1,i+1,2)) / &
                    (cord(j,i,2)+cord(j,i+1,2)-cord(j+1,i,2)-cord(j+1,i+1,2))
                do jj = 1, jm
                    !$ACC atomic update
                    fmagma(jj,i) = fmagma(jj,i) + fmelt(j,i) * ratio_crust_mzone * area_ratio * prod_magma * dt
                enddo
                ! Within mantle, melts migrate by percolation, propagate upward slantly
                depth_moho = 0.5d0 * (cord(jm+1,i,2) + cord(jm+1,i+1,2))
                ! area_ratio: the area of the mantle triangle / the area of the melting element
                area_ratio = 0.5d0 * (depth_moho - 0.5d0 *(cord(j,i,2)+cord(j,i+1,2))) * width_mzone / quad_area
                do ii = max(1,i-ihalfwidth_mzone), min(nx-1,i+ihalfwidth_mzone)
                    do jj = jmoho(ii)+1, j
                        ihw = ihalfwidth_mzone * (j - jj + 1) / (j - jmoho(ii) + 1)
                        if (abs(ii-i) <= ihw) then
                            !$ACC atomic update
                            fmagma(jj,ii) = fmagma(jj,ii) + fmelt(j,i) * ratio_mantle_mzone * area_ratio * prod_magma * dt
                        endif
                    enddo
                enddo
            endif
            fmagma(j,i) = min(fmagma(j,i), fmagma_max)
        enddo
    enddo

    ! Magma freezing rate:
    ! [M(t+dt) - M(t)] / dt  = - M(t) lambda
    ! M(t+dt) = M(t) (1 - dt * lambda)
    !   where lambda is temperature dependent: low T has faster decay
    !   lambda = lambda_freeze * exp(-lambda_freeze_tdep * (T - Tsolidus))
    !
    !$OMP do
    !$ACC parallel loop collapse(2) async(1)
    do i = 1,nx-1
        do j = 1,nz-1
            !iph = iphase(j,i)
            cp_eff = Eff_cp( j,i )
            tmpr = 0.25d0*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

            fmagma_old = fmagma(j,i)
            delta_fmagma = fmagma(j,i) * dt * lambda_freeze * exp(-lambda_freeze_tdep * (tmpr-t_top))
            fmagma(j,i) = max(fmagma(j,i) - delta_fmagma, 0d0)

            ! latent heat released by freezing magma
            !$ACC atomic update
            temp(j  ,i  ) = temp(j  ,i  ) + delta_fmagma * heat_latent_magma / cp_eff / 4
            !$ACC atomic update
            temp(j  ,i+1) = temp(j  ,i+1) + delta_fmagma * heat_latent_magma / cp_eff / 4
            !$ACC atomic update
            temp(j+1,i  ) = temp(j+1,i  ) + delta_fmagma * heat_latent_magma / cp_eff / 4
            !$ACC atomic update
            temp(j+1,i+1) = temp(j+1,i+1) + delta_fmagma * heat_latent_magma / cp_eff / 4
        end do
    enddo
    !$OMP end do
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
        diff = cond_eff/den(iph)/cp_eff

        ! Calculate fluxes in two triangles
        x1 = cord (j  ,i  ,1)
        y1 = cord (j  ,i  ,2)
        x2 = cord (j+1,i  ,1)
        y2 = cord (j+1,i  ,2)
        x3 = cord (j  ,i+1,1)
        y3 = cord (j  ,i+1,2)
        x4 = cord (j+1,i+1,1)
        y4 = cord (j+1,i+1,2)
        t1 = temp (j   ,i  )
        t2 = temp (j+1 ,i  )
        t3 = temp (j   ,i+1)
        t4 = temp (j+1 ,i+1)

        ! Additional sources - radiogenic and shear heating
        tmpr = 0.25d0*(t1 + t2 + t3 + t4)
        !dummye(j,i) = ( source(j,i) + dissip/den(iph) - 600.d0*cp_eff*Eff_melt(iph,tmpr)) / cp_eff
        dummye(j,i) = ( source(j,i) + dissip/den(iph) ) / cp_eff

        ! (1) A element:
        flux(1,1,j,i) = -diff * ( t1*(y2-y3)+t2*(y3-y1)+t3*(y1-y2) ) * area(j,i,1)
        flux(1,2,j,i) = -diff * ( t1*(x3-x2)+t2*(x1-x3)+t3*(x2-x1) ) * area(j,i,1)
 
        ! (2) B element: Interchange of numeration: (1 -> 3,  3 -> 4)
        flux(2,1,j,i) = -diff * ( t3*(y2-y4)+t2*(y4-y3)+t4*(y3-y2) ) * area(j,i,2)
        flux(2,2,j,i) = -diff * ( t3*(x4-x2)+t2*(x3-x4)+t4*(x2-x3) ) * area(j,i,2)

    end do
end do    
!$OMP end do

!$ACC parallel loop collapse(2) async(1)
!$OMP do
do i = 1,nx
    do j = 1,nz

        rhs = 0
        area_n = 0

        ! Element (j-1,i-1). Triangle B
        if ( j.ne.1 .and. i.ne.1 ) then

            ! side 2-3
            qs = flux(2,1,j-1,i-1) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 flux(2,2,j-1,i-1) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            rhs = rhs + 0.5d0*qs

            ! side 3-1
            qs = flux(2,1,j-1,i-1) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 flux(2,2,j-1,i-1) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j-1,i-1,2)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j-1,i-1)*real_area13

        endif

        ! Element (j-1,i). Triangles A,B
        if ( j.ne.1 .and. i.ne.nx ) then

            ! triangle A
            ! side 1-2
            qs = flux(1,1,j-1,i  ) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 flux(1,2,j-1,i  ) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            rhs = rhs + 0.5d0*qs

            ! side 2-3
            qs = flux(1,1,j-1,i  ) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 flux(1,2,j-1,i  ) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j-1,i  ,1)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j-1,i  )*real_area13

            ! triangle B
            ! side 1-2
            qs = flux(2,1,j-1,i  ) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 flux(2,2,j-1,i  ) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            rhs = rhs + 0.5d0*qs

            ! side 2-3
            qs = flux(2,1,j-1,i  ) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 flux(2,2,j-1,i  ) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j-1,i  ,2)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j-1,i  )*real_area13

        endif
        
        ! Element (j,i-1). Triangles A,B
        if ( j.ne.nz .and. i.ne.1 ) then

            ! triangle A
            ! side 2-3
            qs = flux(1,1,j  ,i-1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 flux(1,2,j  ,i-1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            rhs = rhs + 0.5d0*qs

            ! side 3-1
            qs = flux(1,1,j  ,i-1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 flux(1,2,j  ,i-1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j  ,i-1,1)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j  ,i-1)*real_area13

            ! triangle B
            ! side 1-2
            qs = flux(2,1,j  ,i-1) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 flux(2,2,j  ,i-1) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            ! side 3-1
            qs = flux(2,1,j  ,i-1) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 flux(2,2,j  ,i-1) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j  ,i-1,2)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j  ,i-1)*real_area13

        endif

        ! Element (j,i). Triangle A
        if ( j.ne.nz .and. i.ne.nx ) then

            ! side 1-2
            qs = flux(1,1,j  ,i  ) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 flux(1,2,j  ,i  ) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            rhs = rhs + 0.5d0*qs

            ! side 3-1
            qs = flux(1,1,j  ,i  ) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 flux(1,2,j  ,i  ) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            rhs = rhs + 0.5d0*qs

            real_area13 = 0.5d0/area(j  ,i  ,1)/3.d0
            area_n = area_n + real_area13
            rhs = rhs + dummye(j  ,i  )*real_area13

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
