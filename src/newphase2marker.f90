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
    nphase_counter(j,i,mark(kk)%phase) = nphase_counter(j,i,mark(kk)%phase) - 1
    mark(kk)%phase = iph
    nphase_counter(j,i,iph) = nphase_counter(j,i,iph) + 1
enddo

phase_ratio(j,i,:) = 0.0
phase_ratio(j,i,iph) = 1.0

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

nchanged = 0

kocean1 = 3
kocean2 = 7
kcont1 = 2
kcont2 = 6
kmant1 = 4
kmant2 = 8
ksed1 = 10
karc1 = 14
kweak = 12
kserp = 9
kweakmc = 15

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
    ! too deep in the mantle, don't need to consider it further
    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
    if (tmpr.gt.1000.) cycle

    iph = mark(kk)%phase

    ! subducted crust becomes weaker to facilitate further subduction
    if(iph==kocean1 .or. iph==kocean2 .or. iph==karc1 .or. iph==ksed1) then
        do jabove = max(1,j-3), j-1
            if(phase_ratio(jabove,i,kcont1) > 0.8 .or. &
                 phase_ratio(jabove,i,kcont2) > 0.8 ) then
                ! subducted below continent, weakening
                mark(kk)%phase = kweak
                !$OMP critical (change_phase1)
                nphase_counter(j,i,iph) = nphase_counter(j,i,iph) - 1
                nphase_counter(j,i,kweak) = nphase_counter(j,i,kweak) + 1
                nchanged = nchanged + 1
                !$OMP end critical (change_phase1)
                ichanged(nchanged) = i
                jchanged(nchanged) = j

            else if(phase_ratio(jabove,i,kmant1) > 0.8 .or. &
                 phase_ratio(jabove,i,kmant2) > 0.8 ) then
                ! subuducted below mantle, serpentinization
                mark(kk)%phase = kserp
                !$OMP critical (change_phase1)
                nphase_counter(j,i,iph) = nphase_counter(j,i,iph) - 1
                nphase_counter(j,i,kserp) = nphase_counter(j,i,kserp) + 1
                nchanged = nchanged + 1
                !$OMP end critical (change_phase1)
                ichanged(nchanged) = i
                jchanged(nchanged) = j
            endif
        enddo
    endif

    ! XXX: middle crust with high dissipation becomes weaker,
    ! this helps with localization
    !if(iph==kcont1 .or. iph==kcont2 &
    !     .and. tmpr > 300. .and. tmpr < 400. &
    !     .and. stressII(j,i)*strainII(j,i) > 4.e6) then
    !    mark(kk)%phase = kweakmc
    !    !$OMP critical (change_phase1)
    !    nphase_counter(j,i,iph) = nphase_counter(j,i,iph) - 1
    !    nphase_counter(j,i,kweakmc) = nphase_counter(j,i,kweakmc) + 1
    !    nchanged = nchanged + 1
    !    !$OMP end critical (change_phase1)
    !    ichanged(nchanged) = i
    !    jchanged(nchanged) = j
    !endif

    if(nchanged >= 10*mnx) stop 38
enddo
!$OMP end do

!$OMP do
! recompute phase ratio of those changed elements
do k = 1, nchanged
    i = ichanged(k)
    j = jchanged(k)

    kinc = sum(nphase_counter(j,i,:))

    !$OMP critical (change_phase2)
    phase_ratio(j,i,1:nphase) = nphase_counter(j,i,1:nphase) / float(kinc)
    !$OMP end critical (change_phase2)

    ! the phase of this element is the most abundant marker phase
    kph = maxloc(nphase_counter(j,i,:))
    !$OMP critical (change_phase3)
    iphase(j,i) = kph(1)
    !$OMP end critical (change_phase3)
enddo
!$OMP end do
!$OMP end parallel
return
end subroutine change_phase
