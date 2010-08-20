subroutine newphase2marker (j,i,iph)
USE marker_data
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

! reset the markers in this element

iphase(j,i) = iph
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

! Calculate triangle number in which the markers belong
ntriang1 = 2 * ( (nz-1)*(i-1)+j-1) + 1
ntriang2 = 2 * ( (nz-1)*(i-1)+j-1) + 2

!XXX: add omp directive
do kk = 1 , nmarkers
    if (mark(kk)%dead.eq.0) cycle
    n = mark(kk)%ntriag
    if (n.ne.ntriang1 .or. n.ne.ntriang2) cycle
    mark(kk)%maps = aps (j,i)
    mark(kk)%meII = strainII(j,i)
    mark(kk)%mpres = stressI(j,i)
    mark(kk)%mtemp = tmpr
    nphase_counter(mark(kk)%phase,j,i) = nphase_counter(mark(kk)%phase,j,i) - 1
    mark(kk)%phase = iph
    nphase_counter(iph,j,i) = nphase_counter(iph,j,i) + 1
enddo

phase_ratio(:,j,i) = 0.0
phase_ratio(iph,j,i) = 1.0

return
end subroutine newphase2marker


subroutine change_phase
USE marker_data
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

integer ichanged(10*mnx), jchanged(10*mnx)
integer kph(1)
integer, parameter :: kocean1 = 3
integer, parameter :: kocean2 = 7
integer, parameter :: kcont1 = 2
integer, parameter :: kcont2 = 6
integer, parameter :: kmant1 = 4
integer, parameter :: kmant2 = 8
integer, parameter :: ksed1 = 10
integer, parameter :: karc1 = 14
integer, parameter :: kweak = 12
integer, parameter :: kserp = 9
integer, parameter :: kweakmc = 15
integer, parameter :: keclg = 13

nchanged = 0


!$OMP parallel private(kk,i,j,k,n,tmpr,iph,jabove,kinc,kph)
!$OMP do
do kk = 1 , nmarkers
    if (mark(kk)%dead.eq.0) cycle

    ! from ntriag, get element number
    n = mark(kk)%ntriag
    k = mod(n - 1, 2) + 1
    j = mod((n - k) / 2, nz-1) + 1
    i = (n - k) / 2 / (nz - 1) + 1

    ! If temperature of this element is too high, this marker is already
    ! too deep in the mantle, where there is no significant phase change.
    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
    if (tmpr.gt.1000.) cycle

    ! depth below the surface in km
    depth = (cord(1,i,2) - 0.5*(cord(j,i,2)+cord(j+1,i,2)))*1e-3
    iph = mark(kk)%phase
    press = den(iph)*g*depth

    select case(iph)
    case (kcont1, kcont2)
        ! subduction below continent, continent becomes weaker to
        ! facilitate further subduction
        do jbelow = j, min(j+2,nz-1)
            if(phase_ratio(kocean1,jbelow,i) > 0.8 .or. &
                 phase_ratio(kocean2,jbelow,i) > 0.8 .or. &
                 phase_ratio(karc1,jbelow,i) > 0.8 .or. &
                 phase_ratio(ksed1,jbelow,i) > 0.8) then
                !$OMP critical (change_phase1)
                nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
                nphase_counter(kweak,j,i) = nphase_counter(kweak,j,i) + 1
                nchanged = nchanged + 1
                !$OMP end critical (change_phase1)
                mark(kk)%phase = kweak
                ichanged(nchanged) = i
                jchanged(nchanged) = j
            endif
        enddo

        ! XXX: middle crust with high dissipation becomes weaker,
        ! this helps with localization
        !if(tmpr > 300. .and. tmpr < 400. &
        !     .and. stressII(j,i)*strainII(j,i) > 4.e6) then
        !    !$OMP critical (change_phase1)
        !    nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
        !    nphase_counter(kweakmc,j,i) = nphase_counter(kweakmc,j,i) + 1
        !    nchanged = nchanged + 1
        !    !$OMP end critical (change_phase1)
        !    mark(kk)%phase = kweakmc
        !    ichanged(nchanged) = i
        !    jchanged(nchanged) = j
        !endif

    case (kmant1, kmant2)
        ! subuducted oceanic crust below mantle, mantle is serpentinized
        if(depth < 50.) then
            do jbelow = j, min(j+2,nz-1)
                if(phase_ratio(kocean1,jbelow,i) > 0.8 .or. &
                     phase_ratio(kocean2,jbelow,i) > 0.8 .or. &
                     phase_ratio(ksed1,jbelow,i) > 0.8) then
                    !$OMP critical (change_phase1)
                    nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
                    nphase_counter(kserp,j,i) = nphase_counter(kserp,j,i) + 1
                    nchanged = nchanged + 1
                    !$OMP end critical (change_phase1)
                    mark(kk)%phase = kserp
                    ichanged(nchanged) = i
                    jchanged(nchanged) = j
                endif
            enddo
        endif
    case (1, kocean1, kocean2)
        ! basalt -> eclogite
        ! phase change pressure
        trpres = -0.3e9 + 2.2e6*tmpr
        if (tmpr.gt.550 .and. (-1.0*press).ge.trpres) then
            mark(kk)%phase = keclg
            !$OMP critical (change_phase1)
            nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
            nphase_counter(keclg,j,i) = nphase_counter(keclg,j,i) + 1
            nchanged = nchanged + 1
            !$OMP end critical (change_phase1)
            ichanged(nchanged) = i
            jchanged(nchanged) = j
        endif
    case (kserp)
        ! remove serpentinite when it goes below 50 km
        if(depth >= 50.) then
            !$OMP critical (change_phase1)
            nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
            nphase_counter(kmant1,j,i) = nphase_counter(kmant1,j,i) + 1
            nchanged = nchanged + 1
            !$OMP end critical (change_phase1)
            ichanged(nchanged) = i
            jchanged(nchanged) = j
        endif
    end select

    if(nchanged >= 10*mnx) stop 38
enddo
!$OMP end do

!$OMP do
! recompute phase ratio of those changed elements
do k = 1, nchanged
    i = ichanged(k)
    j = jchanged(k)

    kinc = sum(nphase_counter(:,j,i))

    !$OMP critical (change_phase2)
    phase_ratio(1:nphase,j,i) = nphase_counter(1:nphase,j,i) / float(kinc)
    !$OMP end critical (change_phase2)

    ! the phase of this element is the most abundant marker phase
    kph = maxloc(nphase_counter(:,j,i))
    !$OMP critical (change_phase3)
    iphase(j,i) = kph(1)
    !$OMP end critical (change_phase3)
enddo
!$OMP end do
!$OMP end parallel
return
end subroutine change_phase
