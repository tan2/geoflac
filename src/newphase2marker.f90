module newphase2marker

contains

subroutine newphase2marker (j1, j2, i1, i2, iph)
!$ACC routine seq
USE marker_data
use arrays
use params
implicit none

integer :: j1, j2, i1, i2, iph, &
           kk, n, k, j, i

! reset the markers within elements in the rectangular region

do kk = 1 , nmarkers
    if (mark_dead(kk).eq.0) cycle
    n = mark_ntriag(kk)
    k = mod(n - 1, 2) + 1
    j = mod((n - k) / 2, nz-1) + 1
    i = (n - k) / 2 / (nz - 1) + 1

    if(j>=j1 .and. j<=j2 .and. i>=i1 .and. i<=i2) then
        nphase_counter(mark_phase(kk),j,i) = nphase_counter(mark_phase(kk),j,i) - 1
        mark_phase(kk) = iph
        nphase_counter(iph,j,i) = nphase_counter(iph,j,i) + 1
    endif
enddo

iphase(j1:j2,i1:i2) = iph
phase_ratio(:,j1:j2,i1:i2) = 0.d0
phase_ratio(iph,j1:j2,i1:i2) = 1.d0

return
end subroutine newphase2marker


subroutine change_phase
USE marker_data
use arrays
use params
use phases

implicit none

integer :: ichanged(100*(nx-1)), jchanged(100*(nx-1))
integer :: kph(1), jj, j, i, iph, nelem_serp, nchanged, &
           jbelow, k, kinc, kk, n
double precision :: ratio(20), yy, dep2, depth, press, quad_area, &
                    tmpr, trpres, trpres2, vol_frac_melt

! max. depth (m) of eclogite phase transition
real*8, parameter :: max_basalt_depth = 150.e3
! min. temperature (C) of eclogite phase transition
real*8, parameter :: min_eclogite_temp = 500.
real*8, parameter :: min_eclogite_depth = 20e3
real*8, parameter :: mantle_density = 3000.

! temperature (C) of serpentine phase transition
real*8, parameter :: serpentine_temp = 550.

! temperature (C) and depth (m) of 10% partial melting of upper mantle.
real*8, parameter :: partial_melt_temp = 600.
!real*8, parameter :: partial_melt_depth = -70.e3
! thickness of new crust
real*8, parameter :: new_crust_thickness = 7.e3


! search the element for melting
!$ACC parallel create(ichanged, jchanged, kph, ratio)
!$ACC loop
do jj = 1, nz-1
   ! search for crustal depth
   dep2 = 0.25*(cord(jj,1,2)+cord(jj+1,1,2)+cord(jj,2,2)+cord(jj+1,2,2))
   if (cord(1,1,2) - dep2 >= new_crust_thickness) exit
end do
!$ACC end loop
j = min(max(2, jj), nz-1)

!$ACC loop
do i = 1, nx-1
  iph = iphase(j,i)
  if (iph==kmant1 .or. iph==kmant2) then
    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
    if (tmpr > partial_melt_temp) then
      call newphase2marker(1,j-1,i,i,kocean1)
    end if
  end if
end do
!$ACC end loop

! nelem_inject was used for magma injection, reused here for serpentization
nelem_serp = nelem_inject
! rate_inject was used for magma injection, reused here for dehydration melting
vol_frac_melt = rate_inject
andesitic_melt_vol(1:nx-1) = 0

nchanged = 0

!$ACC loop
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
    depth = 0.5*(cord(1,i,2)+cord(1,i+1,2)) - yy

    ! # of markers inside quad
    kinc = sum(nphase_counter(:,j,i))

    !XXX: Some quick checks to skip markers that won't change phase. Might
    !     not be accurate!

    ! If location of this element is too deep, this marker is already
    ! too deep in the mantle, where there is no significant phase change.
    if (depth > 200.e3) cycle

    iph = mark_phase(kk)

    ! Rules of phase changes
    select case(iph)
    case (kcont1, kcont2)
        ! subduction below continent, continent becomes weaker to
        ! facilitate further subduction
        do jbelow = min(j+1,nz-1), min(j+3,nz-1)
            if(phase_ratio(kocean1,jbelow,i) > 0.8 .or. &
                 phase_ratio(kocean2,jbelow,i) > 0.8 .or. &
                 phase_ratio(karc1,jbelow,i) > 0.8 .or. &
                 phase_ratio(ksed1,jbelow,i) > 0.8) then
                ! XXX: TODO: how to critical section for acc?
                !$OMP critical (change_phase1)
                !$ACC atomic update
                nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
                !$ACC atomic update
                nphase_counter(kweak,j,i) = nphase_counter(kweak,j,i) + 1
                !$ACC atomic update
                nchanged = nchanged + 1
                !$ACC atomic write
                ichanged(nchanged) = i
                !$ACC atomic write
                jchanged(nchanged) = j
                !$OMP end critical (change_phase1)
                mark_phase(kk) = kweak
                exit
            endif
        enddo

        ! XXX: middle crust with high dissipation becomes weaker,
        ! this helps with localization
        !if(tmpr > 300. .and. tmpr < 400. &
        !     .and. stressII(j,i)*strainII(j,i) > 4.e6) then
        !    !$OMP critical (change_phase1)
        !    !$ACC atomic update
        !    nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
        !    !$ACC atomic update
        !    nphase_counter(kweakmc,j,i) = nphase_counter(kweakmc,j,i) + 1
        !    !$ACC atomic update
        !    nchanged = nchanged + 1
        !    !$ACC atomic write
        !    ichanged(nchanged) = i
        !    !$ACC atomic write
        !    jchanged(nchanged) = j
        !    !$OMP end critical (change_phase1)
        !    mark_phase(kk) = kweakmc
        !endif

    case (kmant1, kmant2)
        ! subuducted oceanic crust below mantle, mantle is serpentinized
        if(depth > max_basalt_depth) cycle
        ! Phase diagram taken from Ulmer and Trommsdorff, Nature, 1995
        ! Fixed points (730 C, 2.1 GPa) (500 C, 7.5 GPa)
        trpres = 2.1e9 + (7.5e9 - 2.1e9) * (tmpr - 730.) / (500. - 730.)
        ! Fixed points (730 C, 2.1 GPa) (650 C, 0.2 GPa)
        trpres2 = 2.1e9 + (0.2e9 - 2.1e9) * (tmpr - 730.) / (650. - 730.)
        press = mantle_density * g * depth
        if (.not. (press < trpres .and. press > trpres2)) cycle
        do jbelow = min(j+1,nz-1), min(j+nelem_serp,nz-1)
            if(phase_ratio(kocean1,jbelow,i) > 0.8 .or. &
                phase_ratio(kocean2,jbelow,i) > 0.8 .or. &
                phase_ratio(ksed1,jbelow,i) > 0.8) then
                !$OMP critical (change_phase1)
                !$ACC atomic update
                nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
                !$ACC atomic update
                nphase_counter(kserp,j,i) = nphase_counter(kserp,j,i) + 1
                !$ACC atomic update
                nchanged = nchanged + 1
                !$ACC atomic write
                ichanged(nchanged) = i
                !$ACC atomic write
                jchanged(nchanged) = j
                !$OMP end critical (change_phase1)
                mark_phase(kk) = kserp
                exit
            endif
        enddo
    case (kocean0, kocean1, kocean2)
        ! basalt -> eclogite
        ! phase change pressure
        trpres = -0.3e9 + 2.2e6*tmpr
        press = mantle_density * g * depth
        if (tmpr < min_eclogite_temp .or. depth < min_eclogite_depth .or. press < trpres) cycle
        !$OMP critical (change_phase1)
        !$ACC atomic update
        nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
        !$ACC atomic update
        nphase_counter(keclg,j,i) = nphase_counter(keclg,j,i) + 1
        !$ACC atomic update
        nchanged = nchanged + 1
        !$ACC atomic write
        ichanged(nchanged) = i
        !$ACC atomic write
        jchanged(nchanged) = j
        !$OMP end critical (change_phase1)
        mark_phase(kk) = keclg
    case (kserp)
        ! dehydration, serpentinite -> hydrated mantle
        ! Phase diagram taken from Ulmer and Trommsdorff, Nature, 1995
        ! Fixed points (730 C, 2.1 GPa) (500 C, 7.5 GPa)
        trpres = 2.1e9 + (7.5e9 - 2.1e9) * (tmpr - 730.) / (500. - 730.)
        ! Fixed points (730 C, 2.1 GPa) (650 C, 0.2 GPa)
        trpres2 = 2.1e9 + (0.2e9 - 2.1e9) * (tmpr - 730.) / (650. - 730.)
        press = mantle_density * g * depth
        if (tmpr < serpentine_temp .or. (press < trpres .and. press > trpres2)) cycle
        !$OMP critical (change_phase1)
        !$ACC atomic update
        nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
        !$ACC atomic update
        nphase_counter(khydmant,j,i) = nphase_counter(khydmant,j,i) + 1
        !$ACC atomic update
        nchanged = nchanged + 1
        !$ACC atomic write
        ichanged(nchanged) = i
        !$ACC atomic write
        jchanged(nchanged) = j
        !$OMP end critical (change_phase1)
        mark_phase(kk) = khydmant
    case (ksed1, ksed2)
        ! dehydration, sediment -> schist/gneiss
        ! from sediment solidus in Nichols et al., Nature, 1994
        if (tmpr < 650 .or. depth < 20e3) cycle
        !$OMP critical (change_phase1)
        !$ACC atomic update
        nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
        !$ACC atomic update
        nphase_counter(kmetased,j,i) = nphase_counter(kmetased,j,i) + 1
        !$ACC atomic update
        nchanged = nchanged + 1
        !$ACC atomic write
        ichanged(nchanged) = i
        !$ACC atomic write
        jchanged(nchanged) = j
        !$OMP end critical (change_phase1)
        mark_phase(kk) = kmetased
    case (khydmant)
        if (tmpr > ts(khydmant)) then
            ! area(j,i) is INVERSE of "real" DOUBLE area (=1./det)
            quad_area = 1./(area(j,i,1)+area(j,i,2))
            andesitic_melt_vol(i) = andesitic_melt_vol(i) + quad_area * vol_frac_melt / kinc

            !$OMP critical (change_phase1)
            !$ACC atomic update
            nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
            !$ACC atomic update
            nphase_counter(kmant1,j,i) = nphase_counter(kmant1,j,i) + 1
            !$ACC atomic update
            nchanged = nchanged + 1
            !$ACC atomic write
            ichanged(nchanged) = i
            !$ACC atomic write
            jchanged(nchanged) = j
            !$OMP end critical (change_phase1)
            mark_phase(kk) = kmant1
        endif
    end select

    if(nchanged >= 100*(nx-1)) stop 38
enddo
!$OMP end do
!$OMP end parallel
!$ACC end loop

! storing plastic strain in temporary array
junk2(1:nz-1,1:nx-1) = aps(1:nz-1,1:nx-1)

! recompute phase ratio of those changed elements
!$ACC loop
do k = 1, nchanged
    i = ichanged(k)
    j = jchanged(k)

    !if(minval(nphase_counter(:,j,i)) < 0) then
    !    print *, j, i, nphase_counter(:,j,i)
    !    stop 999
    !endif

    kinc = sum(nphase_counter(:,j,i))
    ratio(1:nphase) = nphase_counter(1:nphase,j,i) / float(kinc)
    kph = maxloc(nphase_counter(:,j,i))

    ! the phase of this element is the most abundant marker phase
    iphase(j,i) = kph(1)
    phase_ratio(1:nphase,j,i) = ratio(1:nphase)

    ! When phase change occurs, the mineral would recrystalize and lost
    ! all plastic strain associated with this marker.
    aps(j,i) = max(aps(j,i) - junk2(j,i) / float(kinc), 0d0)

enddo
!$ACC end loop
!$ACC end parallel

return
end subroutine change_phase

end module newphase2marker
