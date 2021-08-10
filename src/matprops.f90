!==============================================
! Density
function Eff_dens( j, i)
  !$ACC routine seq
  use arrays
  use params
  use phases
  include 'precision.inc'

  zcord = 0.25d0*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
  tmpr = 0.25d0*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

  ! 660-km phase change with clapeyron slope 600C/50km
  if (zcord < -660d3 + (tmpr - t_bot) * 500.d0/6) then
      Eff_dens = 3450 * (1 - 2d-5 * tmpr)  ! ~6% density jump
      return
  endif


  press = 0
  do ii = 1, 4
      press = press - (stress0(j,i,1,ii)+stress0(j,i,2,ii)+stress0(j,i,4,ii))
  enddo
  press = press / 12

  Eff_dens = 0.d0
  do k = 1, nphase
    ratio = phase_ratio(k,j,i)
    ! when ratio is small, it won't affect the density
    if(ratio .lt. 0.01d0) cycle

    dens = den(k) * ( 1 - alfa(k)*tmpr + beta(k)*press )

    !press = 0
    !press = stressI(j,i)
    press = dens*g*zcord

    dens = den(k) * ( 1 - alfa(k)*tmpr + beta(k)*press )

    Eff_dens = Eff_dens + ratio*dens
  enddo
  return
end function Eff_dens



!=================================================
! Effective Heat Capacity incorporating latent heat
!=================================================
function Eff_cp( j, i )
!$ACC routine seq
use arrays
use params
implicit none

integer :: iph, j, i
double precision :: Eff_cp

iph = iphase(j,i)
Eff_cp = cp(iph)

return
end function Eff_cp


!=================================================
! Effective Thermal Conductivity
!=================================================
function Eff_conduct( j, i )
!$ACC routine seq
use arrays
use params
implicit none

integer :: iph, j, i
double precision :: Eff_conduct

iph = iphase(j,i)
Eff_conduct = conduct(iph)

return
end function Eff_conduct



!=================================================
! Non-Newtonian viscosity
!=================================================

! from Chen and Morgan (1990)
! Uses A value in MPa and but gives viscosity in (Pa * s) 
! Therefore there is a coefficient 1.e+6 

function Eff_visc( j, i )
!$ACC routine seq
use arrays
use params
include 'precision.inc'


do k = 1, inhom
    !Fixed weak zone
    if (igeom(k).eq.100) then
        if (i >= ix1(k) .and. i <= ix2(k) .and. &
            j >= iy1(k) .and. j <= iy2(k)) then
            Eff_visc = v_min
            return
        endif
    endif
enddo

Eff_visc = 0.d0
r=8.31448d0
tmpr = 0.25d0*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

zcord = 0.25d0*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
if (zcord < -660d3) then
    ! T-dep. viscosity
    ! eta=1e22 when T=1330 Celsius, inc. 10x when T is 200 deg lower
    Eff_visc = 1d22 * exp( -0.0115d0 * (tmpr - 1330))
    Eff_visc = min(v_max, max(v_min, Eff_visc))
    return
elseif (zcord < -410d3) then
    Eff_visc = 1d21 * exp( -0.0115d0 * (tmpr - 1330))
    Eff_visc = min(v_max, max(v_min, Eff_visc))
    return
endif

srat = e2sr(j,i)
if( srat .eq. 0 ) srat = vbc/rxbo

do k = 1, nphase
    if(phase_ratio(k,j,i) .lt. 0.01d0) cycle
    pow  =  1.d0/pln(k) - 1.d0
    pow1 = -1.d0/pln(k)

    vis = 0.25d0 * srat**pow*(0.75d0*acoef(k))**pow1* &
          exp(eactiv(k)/(pln(k)*r*(tmpr+273.d0)))*1.d+6

    if (vis .lt. v_min) vis = v_min
    if (vis .gt. v_max) vis = v_max

    ! harmonic mean
    Eff_visc = Eff_visc + phase_ratio(k,j,i) / vis
    !write(*,*) i,j, Eff_visc, vis, tmpr,phase_ratio(k,j,i)
enddo

Eff_visc = 1 / Eff_visc
if (itype_melting == 1) Eff_visc = Eff_visc * exp(weaken_ratio_viscous * fmagma(j,i) / fmagma_max)

! Final cut-off
Eff_visc = min(v_max, max(v_min, Eff_visc))

return
end function Eff_visc
