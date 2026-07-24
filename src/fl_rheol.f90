
!  Rheology (Update stresses depending on rheology)
!  Calculate total finite strain and plastic strain  
    
subroutine fl_rheol
!$ACC routine(pre_plast) seq
!$ACC routine(elastic) seq
!$ACC routine(maxwell) seq
!$ACC routine(plastic) seq
!$ACC routine(Eff_visc) seq
use arrays
use params
implicit none

double precision :: depl(4)
double precision :: s11p(4),s22p(4),s12p(4),s33p(4),s11v(4),s22v(4),s12v(4),s33v(4)
double precision :: bulkm,rmu,coh,phi,psi, &
                    stherm,hardn,vis, &
                    de11,de22,de12,de33,dv, &
                    diss, poiss, &
                    quad_area, s0, s0a,s0b, &
                    sII_plas, sII_visc, young
double precision :: vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
                    x1,y1,x2,y2,x3,y3,x4,y4, &
                    em,eda,edb, &
                    s11,s22,s12,srII,srI,srs2,stII
double precision :: shpdx_loc(3, 4), shpdz_loc(3, 4), sr_loc(3, 4)
double precision :: Eff_visc
integer :: i, j, k, iph, irh, &
           ipls

! Step 1: Strain rate calculation and dissipation integration
dtavg = dtavg + dt
!$ACC update device(dtavg) async(2)

!$OMP parallel do &
!$OMP private(i,j,k, &
!$OMP         vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP         em,eda,edb,s11,s22,s12, &
!$OMP         srII,srI,srs2,stII, &
!$OMP         x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP         shpdx_loc,shpdz_loc,sr_loc)
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
        sr_loc(1,1) = vx1 * shpdx_loc(1, 1) + vx2 * shpdx_loc(2, 1) + vx3 * shpdx_loc(3, 1)
        sr_loc(2,1) = vy1 * shpdz_loc(1, 1) + vy2 * shpdz_loc(2, 1) + vy3 * shpdz_loc(3, 1)
        sr_loc(3,1) = 0.5d0*(vx1 * shpdz_loc(1, 1) + vx2 * shpdz_loc(2, 1) + vx3 * shpdz_loc(3, 1) + &
                             vy1 * shpdx_loc(1, 1) + vy2 * shpdx_loc(2, 1) + vy3 * shpdx_loc(3, 1))

        !(2) B element: Interchange of numeration: (1 -> 3,  3 -> 4)
        sr_loc(1,2) = vx3 * shpdx_loc(1, 2) + vx2 * shpdx_loc(2, 2) + vx4 * shpdx_loc(3, 2)
        sr_loc(2,2) = vy3 * shpdz_loc(1, 2) + vy2 * shpdz_loc(2, 2) + vy4 * shpdz_loc(3, 2)
        sr_loc(3,2) = 0.5d0*(vx3 * shpdz_loc(1, 2) + vx2 * shpdz_loc(2, 2) + vx4 * shpdz_loc(3, 2) + &
                             vy3 * shpdx_loc(1, 2) + vy2 * shpdx_loc(2, 2) + vy4 * shpdx_loc(3, 2))

        ! (3) C element: ( 3 -> 4 )
        sr_loc(1,3) = vx1 * shpdx_loc(1, 3) + vx2 * shpdx_loc(2, 3) + vx4 * shpdx_loc(3, 3)
        sr_loc(2,3) = vy1 * shpdz_loc(1, 3) + vy2 * shpdz_loc(2, 3) + vy4 * shpdz_loc(3, 3)
        sr_loc(3,3) = 0.5d0*(vx1 * shpdz_loc(1, 3) + vx2 * shpdz_loc(2, 3) + vx4 * shpdz_loc(3, 3) + &
                             vy1 * shpdx_loc(1, 3) + vy2 * shpdx_loc(2, 3) + vy4 * shpdx_loc(3, 3))

        ! (4) D element: (2 -> 4 )
        sr_loc(1,4) = vx1 * shpdx_loc(1, 4) + vx4 * shpdx_loc(2, 4) + vx3 * shpdx_loc(3, 4)
        sr_loc(2,4) = vy1 * shpdz_loc(1, 4) + vy4 * shpdz_loc(2, 4) + vy3 * shpdz_loc(3, 4)
        sr_loc(3,4) = 0.5d0*(vx1 * shpdz_loc(1, 4) + vx4 * shpdz_loc(2, 4) + vx3 * shpdz_loc(3, 4) + &
                             vy1 * shpdx_loc(1, 4) + vy4 * shpdx_loc(2, 4) + vy3 * shpdx_loc(3, 4))

        if ( mix_strain .eq. 1 ) then
            do k = 1, 3, 2
                em = 0.5d0*(sr_loc(1,k)+sr_loc(2,k)+sr_loc(1,k+1)+sr_loc(2,k+1))
                eda = sr_loc(1,k)-sr_loc(2,k)
                edb = sr_loc(1,k+1)-sr_loc(2,k+1)
                sr_loc(1,k) = 0.5d0*(em+eda)
                sr_loc(2,k) = 0.5d0*(em-eda)
                sr_loc(1,k+1) = 0.5d0*(em+edb)
                sr_loc(2,k+1) = 0.5d0*(em-edb)
            enddo
        endif

        ! integration for averaging of strain rate and dissipation function
        s11 = 0.25d0 * sum(sr_loc(1,:))
        s22 = 0.25d0 * sum(sr_loc(2,:))
        s12 = 0.25d0 * sum(sr_loc(3,:))
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

! Step 2: Constitutive stress and strain updates
!$OMP Parallel Private(i,j,k,iph,irh,bulkm,rmu,coh,phi,psi, &
!$OMP                  stherm,hardn,vis, &
!$OMP                  de11,de22,de12,de33,dv, &
!$OMP                  s11p,s22p,s12p,s33p, &
!$OMP                  s11v,s22v,s12v,s33v, &
!$OMP                  depl,ipls,diss, &
!$OMP                  sII_plas,sII_visc, &
!$OMP                  quad_area,s0a,s0b,s0, &
!$OMP                  vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP                  x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP                  em,eda,edb,shpdx_loc,shpdz_loc,sr_loc)
!$OMP do schedule(guided)
!$ACC parallel loop gang vector collapse(2) private(depl,s11p,s22p,s12p,s33p,s11v,s22v,s12v,s33v,shpdx_loc,shpdz_loc,sr_loc) async(1)
do i = 1,nx-1
    do j = 1,nz-1
        ! iphase (j,i) is number of a phase NOT a rheology
        iph = iphase(j,i)
        irh = irheol(iph)

        ! Calculate element strain rates on the fly
        vx1 = vel(j  ,i  ,1)
        vy1 = vel(j  ,i  ,2)
        vx2 = vel(j+1,i  ,1)
        vy2 = vel(j+1,i  ,2)
        vx3 = vel(j  ,i+1,1)
        vy3 = vel(j  ,i+1,2)
        vx4 = vel(j+1,i+1,1)
        vy4 = vel(j+1,i+1,2)

        x1 = cord(j  ,i  ,1)
        y1 = cord(j  ,i  ,2)
        x2 = cord(j+1,i  ,1)
        y2 = cord(j+1,i  ,2)
        x3 = cord(j  ,i+1,1)
        y3 = cord(j  ,i+1,2)
        x4 = cord(j+1,i+1,1)
        y4 = cord(j+1,i+1,2)

        shpdx_loc(1, 1) = (y2 - y3) * area(j, i, 1)
        shpdx_loc(2, 1) = (y3 - y1) * area(j, i, 1)
        shpdx_loc(3, 1) = (y1 - y2) * area(j, i, 1)

        shpdz_loc(1, 1) = (x3 - x2) * area(j, i, 1)
        shpdz_loc(2, 1) = (x1 - x3) * area(j, i, 1)
        shpdz_loc(3, 1) = (x2 - x1) * area(j, i, 1)

        shpdx_loc(1, 2) = (y2 - y4) * area(j, i, 2)
        shpdx_loc(2, 2) = (y4 - y3) * area(j, i, 2)
        shpdx_loc(3, 2) = (y3 - y2) * area(j, i, 2)

        shpdz_loc(1, 2) = (x4 - x2) * area(j, i, 2)
        shpdz_loc(2, 2) = (x3 - x4) * area(j, i, 2)
        shpdz_loc(3, 2) = (x2 - x3) * area(j, i, 2)

        shpdx_loc(1, 3) = (y2 - y4) * area(j, i, 3)
        shpdx_loc(2, 3) = (y4 - y1) * area(j, i, 3)
        shpdx_loc(3, 3) = (y1 - y2) * area(j, i, 3)

        shpdz_loc(1, 3) = (x4 - x2) * area(j, i, 3)
        shpdz_loc(2, 3) = (x1 - x4) * area(j, i, 3)
        shpdz_loc(3, 3) = (x2 - x1) * area(j, i, 3)

        shpdx_loc(1, 4) = (y4 - y3) * area(j, i, 4)
        shpdx_loc(2, 4) = (y3 - y1) * area(j, i, 4)
        shpdx_loc(3, 4) = (y1 - y4) * area(j, i, 4)

        shpdz_loc(1, 4) = (x3 - x4) * area(j, i, 4)
        shpdz_loc(2, 4) = (x1 - x3) * area(j, i, 4)
        shpdz_loc(3, 4) = (x4 - x1) * area(j, i, 4)

        sr_loc(1,1) = vx1 * shpdx_loc(1, 1) + vx2 * shpdx_loc(2, 1) + vx3 * shpdx_loc(3, 1)
        sr_loc(2,1) = vy1 * shpdz_loc(1, 1) + vy2 * shpdz_loc(2, 1) + vy3 * shpdz_loc(3, 1)
        sr_loc(3,1) = 0.5d0*(vx1 * shpdz_loc(1, 1) + vx2 * shpdz_loc(2, 1) + vx3 * shpdz_loc(3, 1) + &
                             vy1 * shpdx_loc(1, 1) + vy2 * shpdx_loc(2, 1) + vy3 * shpdx_loc(3, 1))

        sr_loc(1,2) = vx3 * shpdx_loc(1, 2) + vx2 * shpdx_loc(2, 2) + vx4 * shpdx_loc(3, 2)
        sr_loc(2,2) = vy3 * shpdz_loc(1, 2) + vy2 * shpdz_loc(2, 2) + vy4 * shpdz_loc(3, 2)
        sr_loc(3,2) = 0.5d0*(vx3 * shpdz_loc(1, 2) + vx2 * shpdz_loc(2, 2) + vx4 * shpdz_loc(3, 2) + &
                             vy3 * shpdx_loc(1, 2) + vy2 * shpdx_loc(2, 2) + vy4 * shpdx_loc(3, 2))

        sr_loc(1,3) = vx1 * shpdx_loc(1, 3) + vx2 * shpdx_loc(2, 3) + vx4 * shpdx_loc(3, 3)
        sr_loc(2,3) = vy1 * shpdz_loc(1, 3) + vy2 * shpdz_loc(2, 3) + vy4 * shpdz_loc(3, 3)
        sr_loc(3,3) = 0.5d0*(vx1 * shpdz_loc(1, 3) + vx2 * shpdz_loc(2, 3) + vx4 * shpdz_loc(3, 3) + &
                             vy1 * shpdx_loc(1, 3) + vy2 * shpdx_loc(2, 3) + vy4 * shpdx_loc(3, 3))

        sr_loc(1,4) = vx1 * shpdx_loc(1, 4) + vx4 * shpdx_loc(2, 4) + vx3 * shpdx_loc(3, 4)
        sr_loc(2,4) = vy1 * shpdz_loc(1, 4) + vy4 * shpdz_loc(2, 4) + vy3 * shpdz_loc(3, 4)
        sr_loc(3,4) = 0.5d0*(vx1 * shpdz_loc(1, 4) + vx4 * shpdz_loc(2, 4) + vx3 * shpdz_loc(3, 4) + &
                             vy1 * shpdx_loc(1, 4) + vy4 * shpdx_loc(2, 4) + vy3 * shpdx_loc(3, 4))

        if ( mix_strain .eq. 1 ) then
            do k = 1, 3, 2
                em = 0.5d0*(sr_loc(1,k)+sr_loc(2,k)+sr_loc(1,k+1)+sr_loc(2,k+1))
                eda = sr_loc(1,k)-sr_loc(2,k)
                edb = sr_loc(1,k+1)-sr_loc(2,k+1)
                sr_loc(1,k) = 0.5d0*(em+eda)
                sr_loc(2,k) = 0.5d0*(em-eda)
                sr_loc(1,k+1) = 0.5d0*(em+edb)
                sr_loc(2,k+1) = 0.5d0*(em-edb)
            enddo
        endif

        ! Elastic modules & viscosity & plastic properties
        bulkm = rl(iph) + 2.d0*rm(iph)/3.d0
        rmu   = rm(iph)

        ! Thermal stresses (alfa_v = 3.e-5 1/K)
        stherm = 0.d0
        if (istress_therm.gt.0) stherm = -alfa(iph)*bulkm*(temp(j,i)-temp0(j,i))


        ! Preparation of plastic properties
        if (irh.eq.6 .or. irh.ge.11) call pre_plast(i,j,coh,phi,psi,hardn)
              
        ! Re-evaluate viscosity
        if (irh.eq.3 .or. irh.ge.11) then 
            if( mod(nloop,ifreq_visc).eq.0) visn(j,i) = Eff_visc(j,i)
        endif
        vis = visn(j,i)

        ! Cycle by triangles
        do k = 1,4

            ! Incremental strains
            de11 = sr_loc(1,k)*dt
            de22 = sr_loc(2,k)*dt
            de12 = sr_loc(3,k)*dt
            de33 = 0.d0
            dv = dvol(j,i,k)
            s11p(k) = stress0(j,i,1,k) + stherm 
            s22p(k) = stress0(j,i,2,k) + stherm 
            s12p(k) = stress0(j,i,3,k) 
            s33p(k) = stress0(j,i,4,k) + stherm
            s11v(k) = s11p(k)
            s22v(k) = s22p(k)
            s12v(k) = s12p(k)
            s33v(k) = s33p(k)

            if (irh.eq.1) then
                ! elastic
                call elastic(bulkm,rmu,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de12)
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)

            elseif (irh.eq.3) then
                ! viscous
                call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),de11,de22,de33,de12,dv,&
                             ndim,dt)
                stress0(j,i,1,k) = s11v(k)
                stress0(j,i,2,k) = s22v(k)
                stress0(j,i,3,k) = s12v(k)
                stress0(j,i,4,k) = s33v(k)

            elseif (irh.eq.6) then
                ! plastic
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                             ten_off,ndim)
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)

            elseif (irh.ge.11) then 
                ! Mixed rheology (Maxwell or plastic)
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                    s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                    ten_off,ndim)
                call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                    de11,de22,de33,de12,dv,&
                    ndim,dt)
            endif
        enddo

        if( irh.ge.11 ) then
            ! deside - elasto-plastic or viscous deformation
            sII_plas = (s11p(1)+s11p(2)+s11p(3)+s11p(4)-s22p(1)-s22p(2)-s22p(3)-s22p(4))**2 &
                     + 4*(s12p(1)+s12p(2)+s12p(3)+s12p(4))**2

            sII_visc = (s11v(1)+s11v(2)+s11v(3)+s11v(4)-s22v(1)-s22v(2)-s22v(3)-s22v(4))**2 &
                     + 4*(s12v(1)+s12v(2)+s12v(3)+s12v(4))**2

            if (sII_plas .lt. sII_visc) then
                do k = 1, 4
                    stress0(j,i,1,k) = s11p(k)
                    stress0(j,i,2,k) = s22p(k)
                    stress0(j,i,3,k) = s12p(k)
                    stress0(j,i,4,k) = s33p(k)
                end do
            else 
                do k = 1, 4
                    stress0(j,i,1,k) = s11v(k)
                    stress0(j,i,2,k) = s22v(k)
                    stress0(j,i,3,k) = s12v(k)
                    stress0(j,i,4,k) = s33v(k)
                end do
            endif
        endif


        ! Averaging of isotropic stresses for pair of elements
        if (mix_stress .eq. 1 ) then
        
            ! For A and B couple:
            ! area(n,it) is INVERSE of "real" DOUBLE area (=1./det)
            quad_area = 1.d0/area(j,i,1) + 1.d0/area(j,i,2)
            s0a=0.5d0*(stress0(j,i,1,1)+stress0(j,i,2,1))
            s0b=0.5d0*(stress0(j,i,1,2)+stress0(j,i,2,2))
            s0=(s0a/area(j,i,2)+s0b/area(j,i,1))/quad_area
            stress0(j,i,1,1) = stress0(j,i,1,1) - s0a + s0
            stress0(j,i,2,1) = stress0(j,i,2,1) - s0a + s0
            stress0(j,i,1,2) = stress0(j,i,1,2) - s0b + s0
            stress0(j,i,2,2) = stress0(j,i,2,2) - s0b + s0

            ! For C and D couple:
            quad_area = 1.d0/area(j,i,3) + 1.d0/area(j,i,4)
            s0a=0.5d0*(stress0(j,i,1,3)+stress0(j,i,2,3))
            s0b=0.5d0*(stress0(j,i,1,4)+stress0(j,i,2,4))
            s0=(s0a/area(j,i,4)+s0b/area(j,i,3))/quad_area
            stress0(j,i,1,3) = stress0(j,i,1,3) - s0a + s0
            stress0(j,i,2,3) = stress0(j,i,2,3) - s0a + s0
            stress0(j,i,1,4) = stress0(j,i,1,4) - s0b + s0
            stress0(j,i,2,4) = stress0(j,i,2,4) - s0b + s0
        endif

        if (irh.eq.6 .or. irh.ge.11) then
            !  ACCUMULATED PLASTIC STRAIN
            ! Average the strain for pair of the triangles
            ! Note that area (n,it) is inverse of double area !!!!!
            aps(j,i) = aps(j,i) &
                + 0.5d0*( depl(1)/area(j,i,2)+depl(2)/area(j,i,1) ) / (1.d0/area(j,i,1)+1.d0/area(j,i,2)) &
                + 0.5d0*( depl(3)/area(j,i,4)+depl(4)/area(j,i,3) ) / (1.d0/area(j,i,3)+1.d0/area(j,i,4))
            if( aps(j,i) .lt. 0.d0 ) aps(j,i) = 0.d0

            !	write(*,*) depl(1),depl(2),depl(3),depl(4),area(j,i,1),area(j,i,2),area(j,i,3),area(j,i,4)

            ! LINEAR HEALING OF THE PLASTIC STRAIN
            if (tau_heal .ne. 0.d0) &
                 aps (j,i) = aps (j,i)/(1.d0+dt/tau_heal)
        end if

        ! TOTAL FINITE STRAIN
        strain(j,i,1) = strain(j,i,1) + 0.25d0*dt*(sr_loc(1,1)+sr_loc(1,2)+sr_loc(1,3)+sr_loc(1,4))
        strain(j,i,2) = strain(j,i,2) + 0.25d0*dt*(sr_loc(2,1)+sr_loc(2,2)+sr_loc(2,3)+sr_loc(2,4))
        strain(j,i,3) = strain(j,i,3) + 0.25d0*dt*(sr_loc(3,1)+sr_loc(3,2)+sr_loc(3,3)+sr_loc(3,4))

    enddo
enddo
!$OMP end do
!$OMP end parallel

return
end
