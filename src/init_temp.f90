! Initiate temperature profile

subroutine init_temp
use arrays
use params
include 'precision.inc'

!  Read distribution of temperatures from the dat file
if (irtemp .gt. 0) then
    open( 1, file=tempfile, status='old', err=101 )
    do i = 1, nx
    do j = 1, nz
        read( 1, * ,err=102 ) temp(j,i)
!     if(temp(j,i).ge.1000.d0) temp(j,i) = 1000.d0
    enddo
    enddo
    close(1)

   goto 10
    101 call SysMsg('INIT_TEMP: Cannot open file with temperatures initial distrib!')
    stop 21
    102 call SysMsg('INIT_TEMP: Error reading file with temperatures initial distrib!')
    stop 21
endif

!!  geotherm of a given age accross the box with variable age
do n = 1, nzone_age
    if (n /= 1) then
        if (iph_col_trans(n-1) == 1) cycle
    endif

    do i = ixtb1(n), ixtb2(n)
        if (iph_col_trans(n) == 1) then
            i1 = ixtb1(n)
            i2 = ixtb2(n)
            ratio = (cord(1,i,1) - cord(1,i1,1)) / (cord(1,i2,1) - cord(1,i1,1))
            age_1n = age_1(n) + (age_1(n+1) - age_1(n)) * ratio
            tp1n = tp1(n) + (tp1(n+1) - tp1(n)) * ratio
            tp2n = tp2(n) + (tp2(n+1) - tp2(n)) * ratio
        else
            age_1n = age_1(n)
            tp1n = tp1(n)
            tp2n = tp2(n)
        endif
        call init_geotherm_profile(n, i, age_1n, tp1n, tp2n)
    enddo
enddo

10 continue

! DISTRIBUTE SOURCES in elements
do j = 1,nz-1
    y = -( cord(j+1,1,2)+cord(j,1,2) )/2 / 1000
    source(j,1:nx-1) = hs*exp(-y/hr)
end do

do i = 1, inhom
    ! Initial gaussian temperature perturbation

    ! vertical gaussian
    ! x between (ix1, ix2), y between (iy1, iy2)
    ! gaussian dist. in x direction
    ! linear dist in y direction
    ! xinitaps: amplitude of gaussian
    ! inphase: not used
    if (igeom(i).eq.11) then
        ixc  = (ix1(i)+ix2(i))/2
        iwidth = (ix2(i)-ix1(i))
        amp = xinitaps(i)
        do j = ix1(i),ix2(i)
            pert = amp*exp(-(float(j-ixc)/(0.25d0*float(iwidth)))**2)
            do k = iy1(i),iy2(i)
                pert2 = 1.0d0*(k-iy1(i)) / (iy2(i) - iy1(i))
                temp(k,j) = min(t_bot, temp(k,j)+pert*pert2)
            enddo
        enddo
    endif

    ! slant gaussian
    ! x between (ix1, ix2) at top, shift 1-grid to right for every depth grid
    ! z between (iy1, iy2)
    ! xinitaps: amplitude of gaussian
    ! inphase: not used
    if (igeom(i).eq.13) then
        ixc  = (ix1(i)+ix2(i))/2
        iwidth = (ix2(i)-ix1(i))
        amp = xinitaps(i)
        do k = iy1(i),iy2(i)
            kk = k - iy1(i)
            do j = ix1(i),ix2(i)
                pert = amp*exp(-(float(j-ixc)/(0.25d0*float(iwidth)))**2)
                temp(k,j+kk) = max(t_top, min(t_bot, temp(k,j+kk)+pert))
                !print *, k, j, pert
            enddo
        enddo
    endif

    ! slant gaussian
    ! x between (ix1, ix2) at top, shift 1-grid to left for every depth grid
    ! z between (iy1, iy2)
    ! xinitaps: amplitude of gaussian
    ! inphase: not used
    if (igeom(i).eq.14) then
        ixc  = (ix1(i)+ix2(i))/2
        iwidth = (ix2(i)-ix1(i))
        amp = xinitaps(i)
        do k = iy1(i),iy2(i)
            kk = k - iy1(i)
            do j = ix1(i),ix2(i)
                pert = amp*exp(-(float(j-ixc)/(0.25d0*float(iwidth)))**2)
                temp(k,j-kk) = max(t_top, min(t_bot, temp(k,j-kk)+pert))
                !print *, k, j, pert
            enddo
        enddo
    endif
enddo

!call RedefineTemp

return
end subroutine init_temp


subroutine init_geotherm_profile(n, i, age_1n, tp1n, tp2n)
    !$ACC routine vector
    use arrays
    use params
    include 'precision.inc'

    double precision, parameter :: pi = 3.14159265358979323846d0


    if (ictherm(n)==1) then
        !! Oceanic geotherm (half space cooling model, T&S 3rd ed. Eq(4.113))
        diffusivity = 1.d-6
        age = age_1n
        !$ACC loop
        do j = 1,nz
            f = (cord(1,i,2)-cord(j,i,2)) / sqrt(4 * diffusivity * age * 1.d6 * sec_year)
            temp(j,i) = t_top + (t_bot - t_top) * erf(f)
            !print *, j, age, -cord(j,i,2), temp(j,i)
        enddo
    elseif (ictherm(n)==2) then
        !! Oceanic geotherm (plate cooling cooling model, T&S 3rd ed. Eq(4.130))
        diffusivity = 1.d-6
        age = age_1n
        yL0 = tp1n    ! plate thickness in km
        age_init = age*sec_year*1d6
        tau_d = yL0*yL0*1d6 / (pi*pi*diffusivity)
        !$ACC loop
        do j = 1,nz
            ! depth in km
            y = (cord(1,i,2) - cord(j,i,2)) * 1d-3
            tss = t_top + (t_bot - t_top) * y / yL0
            tt = 0.d0
            do k = 1,100
                tt = tt + 1.d0/k * exp(-k*k*age_init/tau_d) * sin(pi*k*y/yL0)
            enddo
            temp(j,i) = tss + 2.d0/pi*(t_bot-t_top)*tt
            if(temp(j,i)>t_bot .or. y>yL0) temp(j,i) = t_bot
        enddo
    elseif (ictherm(n)==12) then
        !! Continental geotherm (plate cooling model with radiogenic heating)
        !
        ! Starting from the steady state (ss) solution as in T&S 3rd ed. Eq(4.30)
        ! Let the ss moho temperature be tm and ss heatflux be qm.
        ! qm = cond * (t_bot - tm) / (yL0 - ymoho)
        ! Substitute qm to 4.30 to solve for tm
        age = age_1n
        yL0 = tp1n    ! plate thickness in km
        ymoho = tp2n  ! crust thickness in km
        cond = 3.3d0
        dens_c = 2700.d0
        diffusivity = 1.d-6
        !   write(*,*) rzbo, hs, hr
        age_init = age*sec_year*1d6
        tau_d = yL0*yL0*1d6 / (pi*pi*diffusivity)
        rr = ymoho / (yL0 - ymoho)
        tm = (t_top + rr*t_bot + dens_c*hs*hr*hr*1d6/cond*(1d0-exp(-ymoho/hr))) / (1 + rr)
        qm = cond * (t_bot - tm) / (yL0 - ymoho)
        !$ACC loop
        do j = 1,nz
            ! depth in km
            y = (cord(1,i,2) - cord(j,i,2)) * 1d-3

            ! ss part with radiogenic heat
            if (y <= ymoho) then
                tss = t_top + qm/cond*y + (dens_c*hs*hr*hr*1.d+6/cond)*(1d0-exp(-exp(-y/hr)))
            elseif (y <= yL0) then ! below moho, inside lithosphere
                tss = tm + qm/cond*(y-ymoho)
            else
                tss = t_bot
            endif
            ! time-dependent part
            ! see T&S 3rd ed. Eq(4.130)
            tt = 0.d0
            do k = 1,100
                tt = tt + 1.d0/k * exp(-k*k*age_init/tau_d) * sin(pi*k*y/yL0)
            enddo
            temp(j,i) = tss + 2.d0/pi*(t_bot-t_top)*tt
            if(temp(j,i)>t_bot .or. y>yL0) temp(j,i) = t_bot
            !       write(*,*) j,y,tss,yL0,tt
        enddo
    elseif (ictherm(n)==21) then
        !! Constant geotherm gradient at top layer, then T=t_bot all the way to the bottom
        bot_dep = age_1n   ! bottom of the top layer
        do j = 1,nz
            y = (cord(1,i,2)-cord(j,i,2))*1.d-3
            temp(j,i) = y * (t_bot-t_top) / bot_dep
            if(temp(j,i).gt.t_bot) temp(j,i) = t_bot
        enddo
    elseif (ictherm(n)==22) then
        !! Constant geotherm gradient at top two layers, then T=t_bot all the way to the bottom
        temp1 = age_1n
        ylayer1 = tp1n  ! depth in km
        ylayer2 = tp2n  ! depth in km
        !$ACC loop
        do j = 1,nz
            y = (cord(1,i,2)-cord(j,i,2))*1.d-3
            if (y <= ylayer1) then
                temp(j,i) = t_top + (temp1 - t_top) * y / ylayer1
            else
                temp(j,i) = temp1 + (t_bot - temp1) * (y - ylayer1) / (ylayer2 - ylayer1)
            endif
            if(temp(j,i).gt.t_bot) temp(j,i) = t_bot
        enddo
    else
        !call sysmsg('init_temp: ictherm not supported!')
        stop 1
    endif
end subroutine init_geotherm_profile


subroutine sidewalltemp(i1, i2)
  !$ACC routine(init_geotherm_profile) vector

  use arrays, only : temp, cord, source
  use params
  implicit none
  ! This subroutine is intended for remeshing.

  integer, intent(in) :: i1, i2
  double precision :: age_1n, tp1n, tp2n
  integer :: n, i

  if(i1 == 1) then
      ! left sidewall
      n = 1
  else
      ! right sidewall
      n = nzone_age
  endif
  age_1n = age_1(n)
  tp1n = tp1(n)
  tp2n = tp2(n)

  !$ACC parallel loop async(1)
  do i = i1, i2
    call init_geotherm_profile(n, i, age_1n, tp1n, tp2n)
  enddo

  if(i1 == 1) then
      !$ACC parallel loop async(1)
      do i = i1, i2-1
          source(1:nz-1,i) = source(1:nz-1,i2)
      enddo
  else
      !$ACC parallel loop async(1)
      do i = i1, i2-1
          source(1:nz-1,i) = source(1:nz-1,i1-1)
      enddo
  endif
  return
end subroutine sidewalltemp
