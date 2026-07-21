
! Calculation of strain rates from velocities 
! Mixed of volumetric strain_rates (for incompressible flow)
subroutine fl_srate
use arrays
use params
implicit none
integer :: i,j,k
double precision :: vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
         em,eda,edb,s11,s22,s12, &
         srII,srI,srs2,stII
double precision :: x1,y1,x2,y2,x3,y3,x4,y4
double precision :: shpdx_loc(3, 4), shpdz_loc(3, 4)

dtavg = dtavg + dt
!$ACC update device(dtavg) async(2)

!$OMP parallel do &
!$OMP private(i,j,k, &
!$OMP         vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP         em,eda,edb,s11,s22,s12, &
!$OMP         srII,srI,srs2,stII, &
!$OMP         x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP         shpdx_loc,shpdz_loc)
!$ACC parallel loop collapse(2) async(1)
do i = 1,nx-1
    do j = 1,nz-1
 
        ! Velocities:
        vx1 = vel(j  ,i  ,1)
        vy1 = vel(j  ,i  ,2)
        vx2 = vel(j+1,i  ,1)
        vy2 = vel(j+1,i  ,2)
        vx3 = vel(j  ,i+1,1)
        vy3 = vel(j  ,i+1,2)
        vx4 = vel(j+1,i+1,1)
        vy4 = vel(j+1,i+1,2)



        ! 4 different elements (combinations):
        ! b: [d/dx],   c: [d/dy]
 
        ! Coordinates:
        x1 = cord(j  ,i  ,1)
        y1 = cord(j  ,i  ,2)
        x2 = cord(j+1,i  ,1)
        y2 = cord(j+1,i  ,2)
        x3 = cord(j  ,i+1,1)
        y3 = cord(j  ,i+1,2)
        x4 = cord(j+1,i+1,1)
        y4 = cord(j+1,i+1,2)

        ! Triangle A (k=1)
        shpdx_loc(1, 1) = (y2 - y3) * area(j, i, 1)
        shpdx_loc(2, 1) = (y3 - y1) * area(j, i, 1)
        shpdx_loc(3, 1) = (y1 - y2) * area(j, i, 1)

        shpdz_loc(1, 1) = (x3 - x2) * area(j, i, 1)
        shpdz_loc(2, 1) = (x1 - x3) * area(j, i, 1)
        shpdz_loc(3, 1) = (x2 - x1) * area(j, i, 1)

        ! Triangle B (k=2)
        shpdx_loc(1, 2) = (y2 - y4) * area(j, i, 2)
        shpdx_loc(2, 2) = (y4 - y3) * area(j, i, 2)
        shpdx_loc(3, 2) = (y3 - y2) * area(j, i, 2)

        shpdz_loc(1, 2) = (x4 - x2) * area(j, i, 2)
        shpdz_loc(2, 2) = (x3 - x4) * area(j, i, 2)
        shpdz_loc(3, 2) = (x2 - x3) * area(j, i, 2)

        ! Triangle C (k=3)
        shpdx_loc(1, 3) = (y2 - y4) * area(j, i, 3)
        shpdx_loc(2, 3) = (y4 - y1) * area(j, i, 3)
        shpdx_loc(3, 3) = (y1 - y2) * area(j, i, 3)

        shpdz_loc(1, 3) = (x4 - x2) * area(j, i, 3)
        shpdz_loc(2, 3) = (x1 - x4) * area(j, i, 3)
        shpdz_loc(3, 3) = (x2 - x1) * area(j, i, 3)

        ! Triangle D (k=4)
        shpdx_loc(1, 4) = (y4 - y3) * area(j, i, 4)
        shpdx_loc(2, 4) = (y3 - y1) * area(j, i, 4)
        shpdx_loc(3, 4) = (y1 - y4) * area(j, i, 4)

        shpdz_loc(1, 4) = (x3 - x4) * area(j, i, 4)
        shpdz_loc(2, 4) = (x1 - x3) * area(j, i, 4)
        shpdz_loc(3, 4) = (x4 - x1) * area(j, i, 4)

        ! (1) A element:
        strainr(1,1,j,i) = vx1 * shpdx_loc(1, 1) + vx2 * shpdx_loc(2, 1) + vx3 * shpdx_loc(3, 1)
        strainr(2,1,j,i) = vy1 * shpdz_loc(1, 1) + vy2 * shpdz_loc(2, 1) + vy3 * shpdz_loc(3, 1)
        strainr(3,1,j,i) = 0.5d0*(vx1 * shpdz_loc(1, 1) + vx2 * shpdz_loc(2, 1) + vx3 * shpdz_loc(3, 1) + &
                                  vy1 * shpdx_loc(1, 1) + vy2 * shpdx_loc(2, 1) + vy3 * shpdx_loc(3, 1))
 
        !(2) B element: Interchange of numeration: (1 -> 3,  3 -> 4)
        strainr(1,2,j,i) = vx3 * shpdx_loc(1, 2) + vx2 * shpdx_loc(2, 2) + vx4 * shpdx_loc(3, 2)
        strainr(2,2,j,i) = vy3 * shpdz_loc(1, 2) + vy2 * shpdz_loc(2, 2) + vy4 * shpdz_loc(3, 2)
        strainr(3,2,j,i) = 0.5d0*(vx3 * shpdz_loc(1, 2) + vx2 * shpdz_loc(2, 2) + vx4 * shpdz_loc(3, 2) + &
                                  vy3 * shpdx_loc(1, 2) + vy2 * shpdx_loc(2, 2) + vy4 * shpdx_loc(3, 2))
 
        ! (3) C element: ( 3 -> 4 )
        strainr(1,3,j,i) = vx1 * shpdx_loc(1, 3) + vx2 * shpdx_loc(2, 3) + vx4 * shpdx_loc(3, 3)
        strainr(2,3,j,i) = vy1 * shpdz_loc(1, 3) + vy2 * shpdz_loc(2, 3) + vy4 * shpdz_loc(3, 3)
        strainr(3,3,j,i) = 0.5d0*(vx1 * shpdz_loc(1, 3) + vx2 * shpdz_loc(2, 3) + vx4 * shpdz_loc(3, 3) + &
                                  vy1 * shpdx_loc(1, 3) + vy2 * shpdx_loc(2, 3) + vy4 * shpdx_loc(3, 3))

        ! (4) D element: (2 -> 4 )
        strainr(1,4,j,i) = vx1 * shpdx_loc(1, 4) + vx4 * shpdx_loc(2, 4) + vx3 * shpdx_loc(3, 4)
        strainr(2,4,j,i) = vy1 * shpdz_loc(1, 4) + vy4 * shpdz_loc(2, 4) + vy3 * shpdz_loc(3, 4)
        strainr(3,4,j,i) = 0.5d0*(vx1 * shpdz_loc(1, 4) + vx4 * shpdz_loc(2, 4) + vx3 * shpdz_loc(3, 4) + &
                                  vy1 * shpdx_loc(1, 4) + vy4 * shpdx_loc(2, 4) + vy3 * shpdx_loc(3, 4))

        
        !  Mixed discretization of Cundall
        if ( mix_strain .eq. 1 ) then
            do k = 1, 3, 2
                em = 0.5d0*(strainr(1,k,j,i)+strainr(2,k,j,i)+strainr(1,k+1,j,i)+strainr(2,k+1,j,i))
                eda = strainr(1,k,j,i)-strainr(2,k,j,i)
                edb = strainr(1,k+1,j,i)-strainr(2,k+1,j,i)
                strainr(1,k,j,i) = 0.5d0*(em+eda)
                strainr(2,k,j,i) = 0.5d0*(em-eda)
                strainr(1,k+1,j,i) = 0.5d0*(em+edb)
                strainr(2,k+1,j,i) = 0.5d0*(em-edb)
            enddo
        endif


        ! integration for averaging of strain rate and dissipation function
        s11 = 0.25d0 * sum(strainr(1,:,j,i))
        s22 = 0.25d0 * sum(strainr(2,:,j,i))
        s12 = 0.25d0 * sum(strainr(3,:,j,i))
        srII = 0.5d0 * sqrt((s11-s22)**2 + 4*s12*s12)
        srI = (s11+s22)/2
        srs2 = (s11-srI)*(s11-srI) + (s22-srI)*(s22-srI) + 2*s12*s12
        se2sr(j,i,1) = se2sr(j,i,1) + s11*dt
        se2sr(j,i,2) = se2sr(j,i,2) + s22*dt
        se2sr(j,i,3) = se2sr(j,i,3) + s12*dt

        s11 = 0.25d0 * sum(stress0(j,i,1,:))
        s22 = 0.25d0 * sum(stress0(j,i,2,:))
        s12 = 0.25d0 * sum(stress0(j,i,3,:))
        stII = 0.5d0 * sqrt((s11-s22)**2 + 4*s12*s12)
        if( srII.ne.0.d0 ) sshrheat(j,i) = sshrheat(j,i) + stII/srII*srs2*dt
    enddo
enddo
!$OMP end parallel do
! following block is needed for averaging

!$ACC wait(2)
! re-initialisation after navgsr steps
if( nsrate .eq. ifreq_avgsr ) then
    !$OMP parallel do
    !$ACC parallel loop collapse(2) async(1)
    do i = 1,nx-1
        do j = 1, nz-1
            e2sr(j,i) = (0.5d0 * sqrt((se2sr(j,i,1)-se2sr(j,i,2))**2 + 4*se2sr(j,i,3)*se2sr(j,i,3))) / dtavg
            shrheat(j,i) = sshrheat(j,i) / dtavg
            se2sr(j,i,:) = 0.d0
            sshrheat(j,i) = 0.d0
        end do
    end do
    !$OMP end parallel do
    dtavg = 0
    nsrate = 0
elseif( nsrate .eq. -1 ) then
    !$OMP parallel do
    !$ACC parallel loop collapse(2) async(1)
    do i = 1,nx-1
        do j = 1, nz-1
            e2sr(j,i) = (0.5d0 * sqrt((se2sr(j,i,1)-se2sr(j,i,2))**2 + 4*se2sr(j,i,3)*se2sr(j,i,3))) / dtavg
            shrheat(j,i) = sshrheat(j,i) / dtavg
        end do
    end do
    !$OMP end parallel do
endif

nsrate = nsrate + 1
!--------------

return
end 
