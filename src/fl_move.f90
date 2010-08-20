
! Move grid and adjust stresses due to rotation

subroutine fl_move
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'omp_lib.h'


! Move Grid
if (movegrid .eq. 0) return

! UPDATING COORDINATES

!$OMP parallel
!$OMP do
do i = 1,nx
!    write(*,*) cord(j,i,1),cord(j,i,2),vel(j,i,1),vel(j,i,2),dt
    cord(:,i,1) = cord(:,i,1) + vel(:,i,1)*dt
    cord(:,i,2) = cord(:,i,2) + vel(:,i,2)*dt
!    write(*,*) cord(j,i,1),cord(j,i,2)
enddo
!$OMP end do
!$OMP end parallel

! Diffuse topography
if( topo_kappa.gt.0. .OR. bottom_kappa.gt.0. ) call diff_topo


!$OMP parallel private(i,j,x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP                  vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP                  det,dw12,s11,s22,s12)
!$OMP do
!--- Adjusting Stresses And Updating Areas Of Elements
do  i = 1,nx-1
    do  j = 1,nz-1

        ! Coordinates
        x1 = cord (j  ,i  ,1)
        y1 = cord (j  ,i  ,2)
        x2 = cord (j+1,i  ,1)
        y2 = cord (j+1,i  ,2)
        x3 = cord (j  ,i+1,1)
        y3 = cord (j  ,i+1,2)
        x4 = cord (j+1,i+1,1)
        y4 = cord (j+1,i+1,2)

        ! Velocities
        vx1 = vel (j  ,i  ,1)
        vy1 = vel (j  ,i  ,2)
        vx2 = vel (j+1,i  ,1)
        vy2 = vel (j+1,i  ,2)
        vx3 = vel (j  ,i+1,1)
        vy3 = vel (j  ,i+1,2)
        vx4 = vel (j+1,i+1,1)
        vy4 = vel (j+1,i+1,2)

        ! (1) Element A:
        det=((x2*y3-y2*x3)-(x1*y3-y1*x3)+(x1*y2-y1*x2))
        dvol(j,i,1) = det*area(j,i,1) - 1
        area(j,i,1) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x3-x2)+vx2*(x1-x3)+vx3*(x2-x1) - &
            vy1*(y2-y3)-vy2*(y3-y1)-vy3*(y1-y2))/det*dt
        s11 = stress0(j,i,1,1)
        s22 = stress0(j,i,2,1)
        s12 = stress0(j,i,3,1)
        stress0(j,i,1,1) = s11 + s12*2.*dw12
        stress0(j,i,2,1) = s22 - s12*2.*dw12
        stress0(j,i,3,1) = s12 + dw12*(s22-s11)

        ! rotate strains 
        s11 = strain(j,i,1)
        s22 = strain(j,i,2)
        s12 = strain(j,i,3)
        strain(j,i,1) = s11 + s12*2.*dw12
        strain(j,i,2) = s22 - s12*2.*dw12
        strain(j,i,3) = s12 + dw12*(s22-s11)

        ! (2) Element B:
        det=((x2*y4-y2*x4)-(x3*y4-y3*x4)+(x3*y2-y3*x2))
        dvol(j,i,2) = det*area(j,i,2) - 1
        area(j,i,2) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx3*(x4-x2)+vx2*(x3-x4)+vx4*(x2-x3) - &
           vy3*(y2-y4)-vy2*(y4-y3)-vy4*(y3-y2))/det*dt
        s11 = stress0(j,i,1,2)
        s22 = stress0(j,i,2,2)
        s12 = stress0(j,i,3,2)
        stress0(j,i,1,2) = s11 + s12*2.*dw12
        stress0(j,i,2,2) = s22 - s12*2.*dw12
        stress0(j,i,3,2) = s12 + dw12*(s22-s11)

        ! (3) Element C:
        det=((x2*y4-y2*x4)-(x1*y4-y1*x4)+(x1*y2-y1*x2))
        dvol(j,i,3) = det*area(j,i,3) - 1
        area(j,i,3) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x4-x2)+vx2*(x1-x4)+vx4*(x2-x1) - &
           vy1*(y2-y4)-vy2*(y4-y1)-vy4*(y1-y2))/det*dt
        s11 = stress0(j,i,1,3)
        s22 = stress0(j,i,2,3)
        s12 = stress0(j,i,3,3)
        stress0(j,i,1,3) = s11 + s12*2.*dw12
        stress0(j,i,2,3) = s22 - s12*2.*dw12
        stress0(j,i,3,3) = s12 + dw12*(s22-s11)

        ! (4) Element D:
        det=((x4*y3-y4*x3)-(x1*y3-y1*x3)+(x1*y4-y1*x4))
        dvol(j,i,4) = det*area(j,i,4) - 1
        area(j,i,4) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x3-x4)+vx4*(x1-x3)+vx3*(x4-x1) - &
            vy1*(y4-y3)-vy4*(y3-y1)-vy3*(y1-y4))/det*dt
        s11 = stress0(j,i,1,4)
        s22 = stress0(j,i,2,4)
        s12 = stress0(j,i,3,4)
        stress0(j,i,1,4) = s11 + s12*2.*dw12
        stress0(j,i,2,4) = s22 - s12*2.*dw12
        stress0(j,i,3,4) = s12 + dw12*(s22-s11)
    enddo
enddo
!$OMP end do
!$OMP end parallel
return
end subroutine fl_move


!============================================================
! Diffuse topography
!============================================================
subroutine diff_topo
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
parameter( min_elmarkers = 0, max_elmarkers = 12 )

dimension dh(mnx+1),xtgt(mnx)
! Make sure the top layer is made of mobile sediments

if (mod(nloop,10000).eq.0) then
  do ii = 1,nx-1
      call newphase2marker(1,ii,11)
  enddo
endif
!EROSION PROCESSES
if( topo_kappa .gt. 0. ) then             
    xtgtmax = 0.
    dzmax = 0.
    do i = 2, nx-1
        water_depth = 0.5*(cord(1,i+1,2)+cord(1,i,2))
        if (water_depth.lt.0) then
          topo_kappa2 = topo_kappa/10
        else
          topo_kappa2 = topo_kappa
        endif
        xtgt(i) = abs((cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1))) 
        xtgtmax = max(xtgt(i),xtgtmax)
        snder = ( (cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1)) - &
            (cord(1,i  ,2)-cord(1,i-1,2))/(cord(1,i  ,1)-cord(1,i-1,1)) ) / &
            (cord(1,i+1,1)-cord(1,i-1,1))
        dh(i) = topo_kappa2 * dt * snder
!        alphaE = log(10.0)/20000/20000
!        dh(i)= -1.0 * topo_kappa * dt * exp(-1.0*alphaE*(cord(1,i,1)-100000)**2)
        if( dabs(dh(i)).gt.dzmax) dzmax = dabs(dh(i))
!        write(*,*) i,dh(i),snder,dt,topo_kappa
    end do

    if( dzmax*(nz-1) .gt. dabs(rzbo) ) then
         write(*,*) dzmax,nz,rzbo,dt,topo_kappa2,snder
        call SysMsg('DIFF_TOPO: divergence!')
        stop
    endif
!    if(xtgtmax.gt.0.4) then
    do i = 2, nx-1
        cord(1,i,2) = cord(1,i,2) + dh(i)
!        write(*,*) dh(i)
    end do
!   endif
    cord(1,1 ,2) = cord(1,2   ,2)
    cord(1,nx,2) = cord(1,nx-1,2)

endif


! diffuse also bottom boundary
if( bottom_kappa .gt. 0. ) then
    
    dzmax = 0.
    dhsum = 0.
    do i = 2, nx-1
        snder = ( (cord(nz,i+1,2)-cord(nz,i  ,2))/(cord(nz,i+1,1)-cord(nz,i  ,1)) - &
            (cord(nz,i  ,2)-cord(nz,i-1,2))/(cord(nz,i  ,1)-cord(nz,i-1,1)) ) / &
            (cord(nz,i+1,1)-cord(nz,i-1,1))
        dh(i) = bottom_kappa * dt * snder
        dhsum = dhsum + dh(i)
        if( dabs(dh(i)).gt.dzmax) dzmax = dabs(dh(i))
    end do

    new_cord = cord(nz,2,2) + dh(2)
    dh(1) = new_cord - cord(nz,1,2)
    dhsum = dhsum + dh(1)
    if( dabs(dh(1)).gt.dzmax) dzmax = dabs(dh(1))

    new_cord = cord(nz,nx-1,2) + dh(nx-1)
    dh(nx) = new_cord - cord(nz,nx,2)
    dhsum = dhsum + dh(nx)
    if( dabs(dh(nx)).gt.dzmax) dzmax = dabs(dh(nx))

    if( dzmax*(nz-1) .gt. dabs(rzbo) ) then
        call SysMsg('DIFF_TOPO: divergency at the bottom!')
        stop
    endif

    deltah = dhsum/nx
    do i = 1, nx
        cord(nz,i,2) = cord(nz,i,2) + dh(i) - deltah
    end do

endif
       
return
end subroutine diff_topo


subroutine diff_topo_old
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension dh(mnx+1)

             
dzmax = 0.

do i = 2, nx-1
    dx = cord(1,i,1)-cord(1,i-1,1)
    snder = ( (cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1)) - &
              (cord(1,i  ,2)-cord(1,i-1,2))/(cord(1,i  ,1)-cord(1,i-1,1)) ) / &
            (cord(1,i+1,1)-cord(1,i-1,1))
    dh(i) = topo_kappa * dt * snder
    if( dabs(dh(i)).gt.dzmax) dzmax = dabs(dh(i))
end do

if( dzmax*(nz-1) .gt. dabs(rzbo) ) then
    call SysMsg('DIFF_TOPO: divergency!')
    stop
endif

do i = 2, nx-1
    cord(1,i,2) = cord(1,i,2) + dh(i)
end do


! diffuse also bottom boundary
if( bottom_kappa .ne. 0. ) then
    do i = 2, nx-1
        dx = cord(nz,i,1)-cord(nz,i-1,1)
        snder = (cord(nz,i+1,2)-2*cord(nz,i,2)+cord(nz,i-1,2))/4/dx/dx
        dh(i) = bottom_kappa * dt * snder
        if( dabs(dh(i)).gt.dzmax) dzmax = dabs(dh(i))
    end do

    if( dzmax*(nz-1) .gt. dabs(rzbo) ) then
        call SysMsg('DIFF_TOPO: divergency!')
        stop
    endif

    do i = 2, nx-1
        cord(nz,i,2) = cord(nz,i,2) + dh(i)
    end do

endif
       
return
end subroutine diff_topo_old
