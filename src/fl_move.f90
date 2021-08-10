
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
!$ACC parallel loop async(1)
do i = 1,nx
!    write(*,*) cord(j,i,1),cord(j,i,2),vel(j,i,1),vel(j,i,2),dt
    cord(:,i,1) = cord(:,i,1) + vel(:,i,1)*dt
    cord(:,i,2) = cord(:,i,2) + vel(:,i,2)*dt
!    write(*,*) cord(j,i,1),cord(j,i,2)
enddo
!$OMP end do
!$OMP end parallel

! Diffuse topography
if( topo_kappa.gt.0.d0) call diff_topo


!$OMP parallel private(i,j,x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP                  vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP                  det,dw12,s11,s22,s12)
!$OMP do
!--- Adjusting Stresses And Updating Areas Of Elements
!$ACC parallel loop collapse(2) async(1)
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

        if (any(area(j,i,:) <= 0)) then
            ! write(333, *) 'area', j, i, nloop
            ! write(333, *) area(j,i,:)
            ! write(333, *) 'cord:'
            ! write(333, *) cord(j  ,i  ,:)
            ! write(333, *) cord(j  ,i+1,:)
            ! write(333, *) cord(j+1,i  ,:)
            ! write(333, *) cord(j+1,i+1,:)
            ! write(333, *) 'vel:'
            ! write(333, *) vel(j  ,i  ,:)
            ! write(333, *) vel(j  ,i+1,:)
            ! write(333, *) vel(j+1,i  ,:)
            ! write(333, *) vel(j+1,i+1,:)
            ! flush(333)
            ! call SysMsg('Negative area!')
            stop 40
        endif
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
use phases
include 'precision.inc'

!EROSION PROCESSES
if( topo_kappa .gt. 0.d0 ) then

    !$ACC parallel loop async(2)
    do i = 1, nx
        stmpn(i) = topo_kappa ! elevation-dep. topo diffusivity
    enddo

    topomean = 0
    !$ACC parallel loop reduction(+:topomean) async(1)
    do i = 1, nx
        topomean = topomean + cord(1,i,2) / nx
    enddo

    !$ACC wait(2)

    ! !$ACC parallel loop async(1)
    ! do i = 1, nx-1
    !     ! higher erosion for sediments above mean topo
    !     if (iphase(1,i) == ksed2 .and. cord(1,i,2) > topomean) stmpn(i) = stmpn(i) * 10
    ! enddo


    !$ACC parallel loop async(1)
    do i = 2, nx-1

        snder = ( stmpn(i+1)*(cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1)) - &
            stmpn(i-1)*(cord(1,i  ,2)-cord(1,i-1,2))/(cord(1,i  ,1)-cord(1,i-1,1)) ) / &
            (cord(1,i+1,1)-cord(1,i-1,1))
        dtopo(i) = dt * snder
    end do

    !$ACC serial async(1)
    dtopo(1) = dtopo(2)
    dtopo(nx) = dtopo(nx-1)
    !$ACC end serial

    !$ACC parallel loop async(1)
    do i = 1, nx
        ! erosion cannot erode over 0.5x element height
        dtopo(i) = min(dtopo(i), 0.5*(cord(1,i,2)-cord(2,i,2)))
    enddo

    ! accumulated topo change since last resurface
    !$ACC wait(1)
    !$ACC parallel loop async(2)
    do i = 1, nx-1
        dhacc(i) = dhacc(i) + 0.5d0 * (dtopo(i) + dtopo(i + 1))
    enddo

    !$ACC parallel loop async(1)
    do i = 1, nx
        cord(1,i,2) = cord(1,i,2) + dtopo(i)
    enddo
endif

! magma extrusion
arc_extrusion_rate = 1.d0 - ratio_crust_mzone - ratio_mantle_mzone
if (arc_extrusion_rate > 0) then
    !$ACC parallel loop async(1)
    do i = 2, nx-2  ! avoid edge elements, which should not contain arc magma
        totalmelt = 0
        !$ACC loop reduction(+:totalmelt)
        do j = 1, nz-1
            quad_area = 0.5d0/area(j,i,1) + 0.5d0/area(j,i,2)
            ! volume of the melt in this column
            totalmelt = totalmelt + quad_area * fmelt(j,i)
        enddo
        ! height of extrusion in this column
        extrusion(i) = arc_extrusion_rate * dt * totalmelt * prod_magma &
            / (cord(1,i+1,1) - cord(1,i,1) + 0.5d0 * (cord(1,i,1) - cord(1,i-1,1) + cord(1,i+2,1) - cord(1,i+1,1)))
        !print *, i, extrusion(i), totalmelt
        extr_acc(i) = extr_acc(i) + extrusion(i)
        !$ACC atomic update
        cord(1,i,2) = cord(1,i,2) + extrusion(i)
        !$ACC atomic update
        cord(1,i+1,2) = cord(1,i+1,2) + extrusion(i)
    enddo
endif

! adjust markers
if (topo_kappa > 0 .or. arc_extrusion_rate > 0) then
    if(mod(nloop, ifreq_avgsr) .eq. 0) then
!!$        print *, 'max sed/erosion rate (m/yr):' &
!!$             , maxval(dtopo(1:nx)) * 3.16d7 / dt &
!!$             , minval(dtopo(1:nx)) * 3.16d7 / dt
        call resurface
    end if
endif

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

  !$ACC serial private(dz_ratio) async(1)
  do i = 1, nx-1
      ! averge thickness of this element
      elz = 0.5d0 * (cord(1,i,2) - cord(2,i,2) + cord(1,i+1,2) - cord(2,i+1,2))
      ! change in topo
      chgtopo = dhacc(i)
      ! # of markers in this element
      kinc = nmark_elem(1,i)

      ichanged = 0
      if (-chgtopo * kinc >= elz .and. kinc > 1) then
            ! erosion, remove the topmost marker
            ymax = -1d30
            kmax = 0
            ! find the topmost marker in this element
            do k = 1, kinc
                n = mark_id_elem(k, 1, i)
                ntriag = mark_ntriag(n)
                ! get physical coordinate (x, y) of marker n
                m = mark_ntriag(i)
                kk = mod(m-1, 2) + 1
                jj = mod((m - kk) / 2, nz-1) + 1
                ii = (m - kk) / 2 / (nz - 1) + 1
                ba1 = mark_a1(n)
                ba2 = mark_a2(n)
                ba3 = 1.0d0 - ba1 - ba2

                if (kk .eq. 1) then
                  i1 = ii
                  i2 = ii
                  i3 = ii + 1
                  j1 = jj
                  j2 = jj + 1
                  j3 = jj
                else
                  i1 = ii + 1
                  i2 = ii
                  i3 = ii + 1
                  j1 = jj
                  j2 = jj + 1
                  j3 = jj + 1
                endif
                y = cord(j1,i1,2)*ba1 + cord(j2,i2,2)*ba2 + cord(j3,i3,2)*ba3
                if(ymax < y) then
                    ymax = y
                    kmax = k
                endif
            end do
            ! delete (mark it as dead)
            if (kmax .ne. 0) then
                nmax = mark_id_elem(kmax, 1, i)
                !write(6,*) 'erosion', i, nmax, chgtopo
                ! replace marker kmax with the last marker
                mark_id_elem(kmax, 1, i) = mark_id_elem(kinc, 1, i)
                mark_id_elem(kinc, 1, i) = 0
                mark_dead(nmax) = 0
                nmark_elem(1, i) = nmark_elem(1, i) - 1
            endif

            dhacc(i) = 0
            ichanged = 1
      endif

      if (chgtopo * kinc >= elz .and. kinc .ne. max_markers_per_elem) then
            ! sedimentation, add a sediment marker
            !write(6,*) 'sediment', i, chgtopo, elz
            call add_marker_at_top(i, 0.1d0, time, nloop, ksed2)

            dhacc(i) = 0
            ichanged = 1
      endif

      ! change in topo due to volcanism
      chgtopo2 = extr_acc(i)
      if (chgtopo2 * kinc >= elz .and. kinc .ne. max_markers_per_elem) then
            ! extrusion, add an arc marker
            n_to_add = min(ceiling(chgtopo2 / elz * kinc), max_markers_per_elem - kinc)
            dz_ratio = min(chgtopo2 / elz, 1.0d0)
            !write(6,*) 'arc', i, chgtopo2, elz, n_to_add, dz_ratio
            do ii = 1, n_to_add
                call add_marker_at_top(i, dz_ratio, time, nloop+i+ii, karc1)
            enddo

            extr_acc(i) = 0
            ichanged = 1
      endif

      if (ichanged == 1) then
            ! recalculate phase ratio
            call count_phase_ratio(1,i)
      endif
  end do
  !$ACC end serial
  !$ACC update self(nmarkers) async(1)

end subroutine resurface


subroutine add_marker_at_top(i, dz_ratio, time, loop, kph)
  !$ACC routine seq
  !$ACC routine(add_marker) seq
  use myrandom_mod
  use marker_data
  use arrays
  include 'precision.inc'

  iseed = loop + i
  icount = 0
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

     call add_marker(xx, yy, kph, time, 1, i, inc)
     if(inc==1 .or. inc==-1) exit
     icount = icount + 1
     if(icount > 100) stop 134
     !write(333,*) 'add_marker_at_top failed: ', xx, yy, rx, elz, kph
     !write(333,*) '  ', cord(1,i,:)
     !write(333,*) '  ', cord(1,i+1,:)
     !write(333,*) '  ', cord(2,i,:)
     !write(333,*) '  ', cord(2,i+1,:)
     !call SysMsg('Cannot add marker.')
  end do
end subroutine add_marker_at_top
