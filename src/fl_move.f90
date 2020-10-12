
! Move grid and adjust stresses due to rotation

subroutine fl_move
use arrays
use params
include 'precision.inc'


! Move Grid
if (movegrid .eq. 0) return

! UPDATING COORDINATES

!$OMP parallel private(i)
!$OMP do
!$ACC kernels
do i = 1,nx
!    write(*,*) cord(j,i,1),cord(j,i,2),vel(j,i,1),vel(j,i,2),dt
    cord(:,i,1) = cord(:,i,1) + vel(:,i,1)*dt
    cord(:,i,2) = cord(:,i,2) + vel(:,i,2)*dt
!    write(*,*) cord(j,i,1),cord(j,i,2)
enddo
!$ACC end kernels
!$OMP end do
!$OMP end parallel

! Diffuse topography
if( topo_kappa.gt.0.d0) call diff_topo


!$OMP parallel private(i,j,x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP                  vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP                  det,dw12,s11,s22,s12)
!$OMP do
!--- Adjusting Stresses And Updating Areas Of Elements
!$ACC parallel loop collapse(2)
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
        area(j,i,1) = 1.d0/det

        ! Adjusting stresses due to rotation
        dw12 = 0.5d0*(vx1*(x3-x2)+vx2*(x1-x3)+vx3*(x2-x1) - &
            vy1*(y2-y3)-vy2*(y3-y1)-vy3*(y1-y2))/det*dt
        s11 = stress0(j,i,1,1)
        s22 = stress0(j,i,2,1)
        s12 = stress0(j,i,3,1)
        stress0(j,i,1,1) = s11 + s12*2*dw12
        stress0(j,i,2,1) = s22 - s12*2*dw12
        stress0(j,i,3,1) = s12 + dw12*(s22-s11)

        ! rotate strains 
        s11 = strain(j,i,1)
        s22 = strain(j,i,2)
        s12 = strain(j,i,3)
        strain(j,i,1) = s11 + s12*2*dw12
        strain(j,i,2) = s22 - s12*2*dw12
        strain(j,i,3) = s12 + dw12*(s22-s11)

        ! (2) Element B:
        det=((x2*y4-y2*x4)-(x3*y4-y3*x4)+(x3*y2-y3*x2))
        dvol(j,i,2) = det*area(j,i,2) - 1
        area(j,i,2) = 1.d0/det

        ! Adjusting stresses due to rotation
        dw12 = 0.5d0*(vx3*(x4-x2)+vx2*(x3-x4)+vx4*(x2-x3) - &
           vy3*(y2-y4)-vy2*(y4-y3)-vy4*(y3-y2))/det*dt
        s11 = stress0(j,i,1,2)
        s22 = stress0(j,i,2,2)
        s12 = stress0(j,i,3,2)
        stress0(j,i,1,2) = s11 + s12*2*dw12
        stress0(j,i,2,2) = s22 - s12*2*dw12
        stress0(j,i,3,2) = s12 + dw12*(s22-s11)

        ! (3) Element C:
        det=((x2*y4-y2*x4)-(x1*y4-y1*x4)+(x1*y2-y1*x2))
        dvol(j,i,3) = det*area(j,i,3) - 1
        area(j,i,3) = 1.d0/det

        ! Adjusting stresses due to rotation
        dw12 = 0.5d0*(vx1*(x4-x2)+vx2*(x1-x4)+vx4*(x2-x1) - &
           vy1*(y2-y4)-vy2*(y4-y1)-vy4*(y1-y2))/det*dt
        s11 = stress0(j,i,1,3)
        s22 = stress0(j,i,2,3)
        s12 = stress0(j,i,3,3)
        stress0(j,i,1,3) = s11 + s12*2*dw12
        stress0(j,i,2,3) = s22 - s12*2*dw12
        stress0(j,i,3,3) = s12 + dw12*(s22-s11)

        ! (4) Element D:
        det=((x4*y3-y4*x3)-(x1*y3-y1*x3)+(x1*y4-y1*x4))
        dvol(j,i,4) = det*area(j,i,4) - 1
        area(j,i,4) = 1.d0/det

        ! Adjusting stresses due to rotation
        dw12 = 0.5d0*(vx1*(x3-x4)+vx4*(x1-x3)+vx3*(x4-x1) - &
            vy1*(y4-y3)-vy4*(y3-y1)-vy3*(y1-y4))/det*dt
        s11 = stress0(j,i,1,4)
        s22 = stress0(j,i,2,4)
        s12 = stress0(j,i,3,4)
        stress0(j,i,1,4) = s11 + s12*2*dw12
        stress0(j,i,2,4) = s22 - s12*2*dw12
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
use params
include 'precision.inc'

!EROSION PROCESSES
if( topo_kappa .gt. 0.d0 ) then

    !$ACC kernels
    topomean = sum(cord(1,:,2)) / nx
    stmpn = topo_kappa ! elevation-dep. topo diffusivity
    ! higher elevation has higher erosion rate
    do i = 1, nx
        if (cord(1,i,1) > topomean) stmpn(i) = topo_kappa * (1 + (cord(1,i,1) - topomean) * fac_kappa)
    enddo

    do i = 2, nx-1

        snder = ( stmpn(i+1)*(cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1)) - &
            stmpn(i-1)*(cord(1,i  ,2)-cord(1,i-1,2))/(cord(1,i  ,1)-cord(1,i-1,1)) ) / &
            (cord(1,i+1,1)-cord(1,i-1,1))
        dtopo(i) = dt * snder
    end do

    dtopo(1) = dtopo(2)
    dtopo(nx) = dtopo(nx-1)
    cord(1,1:nx,2) = cord(1,1:nx,2) + dtopo(1:nx)

    ! accumulated topo change since last resurface
    dhacc(1:nx-1) = dhacc(1:nx-1) + 0.5d0 * (dtopo(1:nx-1) + dtopo(2:nx))
    !$ACC end kernels

    ! adjust markers
    if(mod(nloop, ifreq_avgsr) .eq. 0) then
!!$        print *, 'max sed/erosion rate (m/yr):' &
!!$             , maxval(dtopo(1:nx)) * 3.16d7 / dt &
!!$             , minval(dtopo(1:nx)) * 3.16d7 / dt
        call resurface
    end if
endif

! magma extrusion
!$ACC kernels
if ( .true. ) then
    ! grid spacing in x
    extrusion(1:nx-1) = cord(1,2:nx,1) - cord(1,1:nx-1,1)
    ! height of extrusion = volume of extrusion / grid spacing in x
    extrusion(1:nx-1) = andesitic_melt_vol(1:nx-1) / extrusion(1:nx-1)
    extr_acc(1:nx-1) = extr_acc(1:nx-1) + extrusion(1:nx-1)

    cord(1,1:nx-1,2) = cord(1,1:nx-1,2) + extrusion(1:nx-1)
    cord(1,2:nx  ,2) = cord(1,2:nx  ,2) + extrusion(1:nx-1)
endif
!$ACC end kernels

return
end subroutine diff_topo




subroutine resurface
  !$ACC routine(bar2xy) seq
  !$ACC routine(shape_functions) seq
  !$ACC routine(add_marker_at_top) seq
  use marker_data
  use arrays
  use params
  use phases
  include 'precision.inc'

  dimension shp2(2,3,2)

  !$ACC kernels
  do i = 1, nx-1
      ! averge thickness of this element
      elz = 0.5d0 * (cord(1,i,2) - cord(2,i,2) + cord(1,i+1,2) - cord(2,i+1,2))
      ! change in topo
      chgtopo = dhacc(i)
      ! # of markers in this element
      kinc = nmark_elem(1,i)

      if (abs(chgtopo*kinc) >= elz) then
          ! add/remove markers if topo changed too much
          if (chgtopo > 0.d0) then
              ! sedimentation, add a sediment marker
              !print *, 'add sediment', i, chgtopo, elz
              call add_marker_at_top(i, 0.1d0, time, ksed2, nmarkers)
          else
              ! erosion, remove the top marker
              !print *, 'erosion', i, chgtopo, elz
              ymax = -1d30
              nmax = 0
              kmax = 0
              call shape_functions(1,i,shp2)
              do k = 1, nmark_elem(1,i)
                  n = mark_id_elem(k, 1, i)
                  ntriag = mark_ntriag(n)
                  m = mod(ntriag,2) + 1
                  call bar2xy(mark_a1(n), mark_a2(n), shp2(:,:,m), x, y)
                  if(ymax < y) then
                      ymax = y
                      nmax = n
                      kmax = k
                  endif
              end do
              if (nmax .ne. 0) then
                  mark_dead(nmax) = 0
                  ! replace marker kmax with last marker
                  mark_id_elem(kmax, 1, i) = mark_id_elem(nmark_elem(1, i), 1, i)
                  nmark_elem(1, i) = nmark_elem(1, i) - 1
              endif
          endif

          dhacc(i) = 0

          ! recalculate phase ratio
          call count_phase_ratio(1,i)

      else
          ! nothing to do
      end if

      ! change in topo
      chgtopo2 = extr_acc(i)

      if (chgtopo2*kinc >= elz) then
          ! add/remove markers if topo changed too much
          ! extrusion, add an arc marker
          n_to_add = ceiling(chgtopo2 / elz * kinc)
          dz_ratio = min(chgtopo2 / elz, 1.0d0)
          !print *, 'add arc', i, chgtopo2, elz, n_to_add, dz_ratio
          do ii = 1, n_to_add
              call add_marker_at_top(i, dz_ratio, time, karc1, nmarkers)
          enddo

          extr_acc(i) = 0

          ! recalculate phase ratio
          call count_phase_ratio(1,i)

      end if
  end do
  !$ACC end kernels
  !$ACC update self(nmarkers)

end subroutine resurface


subroutine add_marker_at_top(i, dz_ratio, time, kph, nmarkers)
  !$ACC routine seq
  !$ACC routine(add_marker) seq
  use myrandom_mod
  use marker_data
  use arrays
  include 'precision.inc'

  iseed = nloop
  do while(.true.)
     call myrandom(iseed, r1)
     call myrandom(iseed, r2)
     j = 1

     ! (x1, y1) and (x2, y2)
     x1 = cord(j  ,i,1)*(1-r1) + cord(j  ,i+1,1)*r1
     y1 = cord(j  ,i,2)*(1-r1) + cord(j  ,i+1,2)*r1
     x2 = cord(j+1,i,1)*(1-r1) + cord(j+1,i+1,1)*r1
     y2 = cord(j+1,i,2)*(1-r1) + cord(j+1,i+1,2)*r1

     ! connect the above two points
     ! (this point is not uniformly distributed within the element area
     ! and is biased against the thicker side of the element, but this
     ! point is almost gauranteed to be inside the element)
     r2 = r2 * dz_ratio
     xx = x1*(1-r2) + x2*r2
     yy = y1*(1-r2) + y2*r2

     call add_marker(xx, yy, kph, time, nmarkers, 1, i, inc)
     if(inc==1 .or. inc==-1) exit
     !write(333,*) 'add_marker_at_top failed: ', xx, yy, rx, elz, kph
     !write(333,*) '  ', cord(1,i,:)
     !write(333,*) '  ', cord(1,i+1,:)
     !write(333,*) '  ', cord(2,i,:)
     !write(333,*) '  ', cord(2,i+1,:)
     !call SysMsg('Cannot add marker.')
  end do
end subroutine add_marker_at_top
