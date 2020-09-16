subroutine change_phase
USE marker_data
use arrays
use params
use phases

implicit none

integer :: kph(1), jj, j, i, iph, nelem_serp, &
           jbelow, k, kinc, kk, n
double precision :: ratio(20), yy, dep2, depth, press, quad_area, &
                    tmpr, trpres, trpres2, vol_frac_melt

! max. depth (m) of eclogite phase transition
real*8, parameter :: max_basalt_depth = 150.d3
! min. temperature (C) of eclogite phase transition
real*8, parameter :: min_eclogite_temp = 500.d0
real*8, parameter :: min_eclogite_depth = 20d3
real*8, parameter :: mantle_density = 3000.d0

! temperature (C) of serpentine phase transition
real*8, parameter :: serpentine_temp = 550.d0

! temperature (C) and depth (m) of 10% partial melting of upper mantle.
real*8, parameter :: partial_melt_temp = 600.d0
!real*8, parameter :: partial_melt_depth = -70.d3
! thickness of new crust
real*8, parameter :: new_crust_thickness = 7.d3

! search the element for melting
do jj = 1, nz-1
   ! search for crustal depth
   dep2 = 0.25d0*(cord(jj,1,2)+cord(jj+1,1,2)+cord(jj,2,2)+cord(jj+1,2,2))
   if (cord(1,1,2) - dep2 >= new_crust_thickness) exit
end do
j = min(max(2, jj), nz-1)

do i = 1, nx-1
  iph = iphase(j,i)
  if (iph==kmant1 .or. iph==kmant2) then
    tmpr = 0.25d0*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
    if (tmpr > partial_melt_temp) then
      call newphase2marker(1,j-1,i,i,kocean1)
    end if
  end if
end do


! nelem_inject was used for magma injection, reused here for serpentization
nelem_serp = nelem_inject
! rate_inject was used for magma injection, reused here for dehydration melting
vol_frac_melt = rate_inject
andesitic_melt_vol(1:nx-1) = 0

itmp = 0  ! indicates which element has phase-changed markers


!$OMP parallel private(kk,i,j,k,n,tmpr,depth,iph,press,jbelow,trpres,trpres2,kinc,quad_area,yy)
!$OMP do schedule(guided)
do kk = 1 , nmarkers
    if (mark_dead(kk).eq.0) cycle

    ! from ntriag, get element number
    n = mark_ntriag(kk)
    k = mod(n - 1, 2) + 1
    j = mod((n - k) / 2, nz-1) + 1
    i = (n - k) / 2 / (nz - 1) + 1

    if (k .eq. 1) then
       yy = cord(j,i,2)*mark_a1(kk) + cord(j+1,i,2)*mark_a2(kk) + cord(j,i+1,2)*(1-mark_a1(kk)-mark_a2(kk))
       tmpr = temp(j,i)*mark_a1(kk) + temp(j+1,i)*mark_a2(kk) + temp(j,i+1)*(1-mark_a1(kk)-mark_a2(kk))
    else
       yy = cord(j,i+1,2)*mark_a1(kk) + cord(j+1,i,2)*mark_a2(kk) + cord(j+1,i+1,2)*(1-mark_a1(kk)-mark_a2(kk))
       tmpr = temp(j,i+1)*mark_a1(kk) + temp(j+1,i)*mark_a2(kk) + temp(j+1,i+1)*(1-mark_a1(kk)-mark_a2(kk))
    endif

    ! depth below the surface in m
    depth = 0.5d0*(cord(1,i,2)+cord(1,i+1,2)) - yy

    ! # of markers inside quad
    kinc = nmark_elem(j,i)

    !XXX: Some quick checks to skip markers that won't change phase. Might
    !     not be accurate!

    ! If location of this element is too deep, this marker is already
    ! too deep in the mantle, where there is no significant phase change.
    if (depth > 200.d3) cycle

    iph = mark_phase(kk)

    ! Rules of phase changes
    select case(iph)
    case (kcont1, kcont2)
        ! subduction below continent, continent becomes weaker to
        ! facilitate further subduction
        do jbelow = min(j+1,nz-1), min(j+3,nz-1)
            if(phase_ratio(kocean1,jbelow,i) > 0.8d0 .or. &
                 phase_ratio(kocean2,jbelow,i) > 0.8d0 .or. &
                 phase_ratio(karc1,jbelow,i) > 0.8d0 .or. &
                 phase_ratio(ksed1,jbelow,i) > 0.8d0) then
                !$OMP critical (change_phase1)
                itmp(j,i) = 1
                !$OMP end critical (change_phase1)
                mark_phase(kk) = kweak
                exit
            endif
        enddo

        ! XXX: middle crust with high dissipation becomes weaker,
        ! this helps with localization
        !if(tmpr > 300.d0 .and. tmpr < 400.d0 &
        !     .and. stressII(j,i)*strainII(j,i) > 4.d6) then
        !    !$OMP critical (change_phase1)
        !    !itmp(j,i) = 1
        !    !$OMP end critical (change_phase1)
        !    mark_phase(kk) = kweakmc
        !endif

    case (kmant1, kmant2)
        ! subuducted oceanic crust below mantle, mantle is serpentinized
        if(depth > max_basalt_depth) cycle
        ! Phase diagram taken from Ulmer and Trommsdorff, Nature, 1995
        ! Fixed points (730 C, 2.1 GPa) (500 C, 7.5 GPa)
        trpres = 2.1d9 + (7.5d9 - 2.1d9) * (tmpr - 730.d0) / (500.d0 - 730.d0)
        ! Fixed points (730 C, 2.1 GPa) (650 C, 0.2 GPa)
        trpres2 = 2.1d9 + (0.2d9 - 2.1d9) * (tmpr - 730.d0) / (650.d0 - 730.d0)
        press = mantle_density * g * depth
        if (.not. (press < trpres .and. press > trpres2)) cycle
        do jbelow = min(j+1,nz-1), min(j+nelem_serp,nz-1)
            if(phase_ratio(kocean1,jbelow,i) > 0.8d0 .or. &
                phase_ratio(kocean2,jbelow,i) > 0.8d0 .or. &
                phase_ratio(ksed1,jbelow,i) > 0.8d0) then
                !$OMP critical (change_phase1)
                itmp(j,i) = 1
                !$OMP end critical (change_phase1)
                mark_phase(kk) = kserp
                exit
            endif
        enddo
    case (kocean0, kocean1, kocean2)
        ! basalt -> eclogite
        ! phase change pressure
        trpres = -0.3d9 + 2.2d6*tmpr
        press = mantle_density * g * depth
        if (tmpr < min_eclogite_temp .or. depth < min_eclogite_depth .or. press < trpres) cycle
        !$OMP critical (change_phase1)
        itmp(j,i) = 1
        !$OMP end critical (change_phase1)
        mark_phase(kk) = keclg
    case (kserp)
        ! dehydration, serpentinite -> hydrated mantle
        ! Phase diagram taken from Ulmer and Trommsdorff, Nature, 1995
        ! Fixed points (730 C, 2.1 GPa) (500 C, 7.5 GPa)
        trpres = 2.1d9 + (7.5d9 - 2.1d9) * (tmpr - 730.d0) / (500.d0 - 730.d0)
        ! Fixed points (730 C, 2.1 GPa) (650 C, 0.2 GPa)
        trpres2 = 2.1d9 + (0.2d9 - 2.1d9) * (tmpr - 730.d0) / (650.d0 - 730.d0)
        press = mantle_density * g * depth
        if (tmpr < serpentine_temp .or. (press < trpres .and. press > trpres2)) cycle
        !$OMP critical (change_phase1)
        itmp(j,i) = 1
        !$OMP end critical (change_phase1)
        mark_phase(kk) = khydmant
    case (ksed1, ksed2)
        ! dehydration, sediment -> schist/gneiss
        ! from sediment solidus in Nichols et al., Nature, 1994
        if (tmpr < 650d0 .or. depth < 20d3) cycle
        !$OMP critical (change_phase1)
        itmp(j,i) = 1
        !$OMP end critical (change_phase1)
        mark_phase(kk) = kmetased
    case (khydmant)
        if (tmpr > ts(khydmant)) then
            ! area(j,i) is INVERSE of "real" DOUBLE area (=1./det)
            quad_area = 0.5d0/area(j,i,1) + 0.5d0/area(j,i,2)
            andesitic_melt_vol(i) = andesitic_melt_vol(i) + quad_area * vol_frac_melt / kinc

            !$OMP critical (change_phase1)
            itmp(j,i) = 1
            !$OMP end critical (change_phase1)
            mark_phase(kk) = kmant1
        endif
    end select

enddo
!$OMP end do
!$OMP end parallel

! storing plastic strain in temporary array
dummy(1:nz-1,1:nx-1) = aps(1:nz-1,1:nx-1)

!$OMP parallel do private(iph, kinc)
! recompute phase ratio of those changed elements
do i = 1, nx-1
    do j = 1, nz-1

        ! skip unchanged element
        if (itmp(j,i) == 0) cycle

        ! the phase of this element is the most abundant marker phase
        call count_phase_ratio(j,i)

        ! When phase change occurs, the mineral would recrystalize and lost
        ! all plastic strain associated with this marker.
        kinc = nmark_elem(j,i)
        aps(j,i) = max(aps(j,i) - dummy(j,i) / float(kinc), 0d0)
    enddo
enddo
!$OMP end parallel do

return
end subroutine change_phase
