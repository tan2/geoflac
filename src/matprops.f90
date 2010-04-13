!==============================================
! Density
function Eff_dens( j, i)
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

Eff_dens = 0.

iph = iphase(i,j,phasez(j,i))

zcord = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
 tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
press = 0
do ii = 1, 4
    press = press - (stress0(1,ii,j,i)+stress0(2,ii,j,i)+stress0(4,ii,j,i))/3
enddo
press = press / 4
dens = den(iph) * ( 1 - alfa(iph)*tmpr + beta(iph)*press )

! Effect of melt
fmelt = Eff_melt2(iph, tmpr)
dens = dens * ( 1.-0.1*fmelt )
Eff_dens = dens

if (iint_marker.eq.1) then
Eff_dens = 0.
do k = 1, nphasl
   
iph = lphase(k)
delta_rho= 0.


ratio = phase_ratio (j,i,k) 
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
press = 0
do ii = 1, 4
   press = press - (stress0(1,ii,j,i)+stress0(2,ii,j,i)+stress0(4,ii,j,i))/3
enddo
press = press / 4
dens = den(iph) * ( 1 - alfa(iph)*tmpr + beta(iph)*press )
zcord = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
press = 0
!press = stressI(j,i)
press = dens*g*zcord
!    if (iph.eq.1.or.iph.eq.3.or.iph.eq.7.or.iph.eq.2.or.iph.eq.6) then
!if(tmpr.gt.110.) then
!            trpres = 1.3e9
!           if ((-1.0*press).ge.trpres) then
!              delta_rho = 600.
!           endif
!       endif
!       if (tmpr.gt.550) then
!           trpres = -0.3e9 + 2.2e6*tmpr
!           if ((-1.0*press).ge.trpres) then
!              delta_rho = 400.
!	endif
!        endif
!     endif
dens = (den(iph)+delta_rho) * ( 1 - alfa(iph)*tmpr + beta(iph)*press )
if(iph.eq.11) then
         delta_den = 400.
         zefold = 6000.
!dens = (den(iph) - delta_den*exp(zcord/zefold)) * ( 1 - alfa(iph)*tmpr + beta(iph)*press )
dens = (den(iph) - delta_den*exp((zcord-0.5*(cord(1,i,2)+cord(1,i+1,2)))/zefold)) * ( 1 - alfa(iph)*tmpr + beta(iph)*press )
if (j.eq.1.and.cord(j,i,2).gt.0.) dens = (2750.*(1.-alfa(iph)*tmpr))        
       if (dens.lt.2400.) dens = 2400.
     endif
         
! Effect of melt
!fmelt = Eff_melt2(iph, tmpr)
!dens = dens * ( 1.-0.1*fmelt )

Eff_dens = Eff_dens + ratio*dens

enddo
endif
return
end function Eff_dens


!==============================================
! Melt fraction
!==============================================
function Eff_melt( j, i )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


iph = iphase(i,j,phasez(j,i))
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

Eff_melt = Eff_melt2(iph, tmpr)
return
end function Eff_melt


!==============================================
function Eff_melt2(iph, tmpr)
include 'precision.inc'
include 'params.inc'

if( tmpr .lt. ts(iph) ) then
    ! below solidus
    fm = 0.
elseif( tmpr .lt. tk(iph) ) then
    fm = fk(iph)/(tk(iph)-ts(iph)) * (tmpr-ts(iph))
elseif( tmpr .lt. tl(iph) ) then
    fm = (1.-fk(iph))/(tl(iph)-tk(iph))*(tmpr-tk(iph)) + fk(iph)
else
    fm = 1.
endif

if( fm .lt. 0 ) fm = 0.
if( fm .gt. 1 ) fm = 1.

Eff_melt2 = fm

return
end function Eff_melt2


!=================================================
! Effective Heat Capacity incorporating latent heat
!=================================================
function Eff_cp( j, i )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


HeatLatent = 420000.

iph = iphase(i,j,phasez(j,i))
!Eff_cp = cp(iph)


tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
if( tmpr .lt. ts(iph) ) then
    Eff_cp = cp(iph)
elseif( tmpr .lt. tk(iph) ) then
    Eff_cp = cp(iph) + HeatLatent * fk(iph)/(tk(iph)-ts(iph))
elseif( tmpr .lt. tl(iph) ) then
    Eff_cp = cp(iph) + HeatLatent * (1.-fk(iph))/(tl(iph)-tk(iph))
else
    Eff_cp = cp(iph)
endif


! HOOK
! Intrusions - melting effect - see user_ab.f90
if( if_intrus .eq. 1 ) then
    HeatLatent = 420000.

    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
    if( tmpr .lt. ts(iph) ) then
        Eff_cp = cp(iph)
    elseif( tmpr .lt. tl(iph)+1 ) then
        Eff_cp = cp(iph) + HeatLatent/(tl(iph)-ts(iph))
    else
        Eff_cp = cp(iph)
    endif
endif

return
end function Eff_cp


!=================================================
! Effective Thermal Conductivity
!=================================================
function Eff_conduct( j, i )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
Eff_conduct = 0.

iph = iphase(i,j,phasez(j,i))
cond = conduct(iph)

!if( den(iph) .lt. 3000. ) then  ! for crustal material
!    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
!    if( tmpr.lt.25 ) tmpr = 25.
!    Eff_conduct = -0.38*dlog(tmpr) + 4.06
!endif

! HOOK
! Hydrothermal alteration of thermal diffusivity  - see user_luc.f90
if( if_hydro .eq. 1 ) then
    cond = HydroDiff(i,j)*den(iph)*cp(iph)
endif
 Eff_conduct = cond

if (iint_marker.eq.1) then
Eff_conduct = 0.
do k = 1 , nphasl

iph = lphase(k) 
cond = conduct(iph)

!if( den(iph) .lt. 3000. ) then  ! for crustal material
!    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
!    if( tmpr.lt.25 ) tmpr = 25.
!    Eff_conduct = -0.38*dlog(tmpr) + 4.06
!endif

! HOOK
! Hydrothermal alteration of thermal diffusivity  - see user_luc.f90
if( if_hydro .eq. 1 ) then
    cond = HydroDiff(i,j)*den(iph)*cp(iph)
endif
Eff_conduct = Eff_conduct + phase_ratio(j,i,k)*cond
enddo
endif
!write(*,*) Eff_conduct, cond
return
end function Eff_conduct



!=================================================
! Non-Newtonian viscosity
!=================================================

! from Chen and Morgan (1990)
! Uses A value in MPa and but gives viscosity in (Pa * s) 
! Therefore there is a coefficient 1.e+6 

function Eff_visc( j, i )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
Eff_visc = 0.
r=8.31448e0
iph = iphase(i,j,phasez(j,i))
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
srat = e2sr(j,i)
if( srat .eq. 0 ) srat = vbc/rxbo

pow  =  1./pln(iph) - 1.
pow1 = -1./pln(iph) 

vis = 0.25 * srat**pow*(0.75*acoef(iph))**pow1* &
        exp(eactiv(iph)/(pln(iph)*r*(tmpr+273.)))*1.e+6

! Effect of melt
fmelt_crit = 0.05
fmelt = Eff_melt2(iph, tmpr)
if( fmelt .gt. 0. ) then
    if( fmelt .lt. fmelt_crit ) then
        vislog = fmelt/fmelt_crit*dlog10(v_min/vis) + dlog10(vis)
        vis = 10.**vislog
    else
        vis = v_min
    endif
endif


! limiting from above (quasi-Peierls)
!sIImax = 5.e+8
!vis_peierls = sIImax / srat / 2
!if( vis .gt. vis_peierls ) vis = vis_peierls


! Final cut-off
if (vis .lt. v_min) vis = v_min
if (vis .gt. v_max) vis = v_max

Eff_visc = vis

if (iint_marker.eq.1) then
Eff_visc = 0.
do k = 1, nphasl

iph = lphase(k) 
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
srat = e2sr(j,i)
press = stressI(j,i)

xeactiv = eactiv(iph)
    if (iph.eq.20.or.iph.eq.30) then
       if (tmpr.gt.550) then
           trpres = -0.3e9 + 2.2e6*tmpr
           if ((-1.0*press).ge.trpres) then
               xeactiv  = 3.5e5
           endif
        endif
     endif


if( srat .eq. 0 ) srat = vbc/rxbo

pow  =  1./pln(iph) - 1.
pow1 = -1./pln(iph) 

vis = 0.25 * srat**pow*(0.75*acoef(iph))**pow1* &
        exp(xeactiv/(pln(iph)*r*(tmpr+273.)))*1.e+6

!increase continental mantle viscosity (one order of magnitude)
if (iph.eq.5) vis = 10*vis
!!!!!!!!!!!!!!!!

! Effect of melt
!fmelt_crit = 0.05
!fmelt = Eff_melt2(iph, tmpr)
!if( fmelt .gt. 0. ) then
!    if( fmelt .lt. fmelt_crit ) then
!        vislog = fmelt/fmelt_crit*dlog10(v_min/vis) + dlog10(vis)
!        vis = 10.**vislog
!    else
!        vis = v_min
    !endif
!endif


! limiting from above (quasi-Peierls)
!sIImax = 5.e+8
!vis_peierls = sIImax / srat / 2
!if( vis .gt. vis_peierls ) vis = vis_peierls

!write(*,*) vis,srat,pln(iph),xeactiv
! Final cut-off
if (vis .lt. v_min) vis = v_min
if (vis .gt. v_max) vis = v_max
if(iph.eq.8.and.vis.le.1.e19)  vis = 1.e19     
if(iph.eq.4.and.vis.le.1.e19)  vis = 1.e19     
Eff_visc = Eff_visc + phase_ratio(j,i,k)*vis
!write(*,*) i,j, Eff_visc, vis, tmpr,phase_ratio(j,i,k)
enddo
if(vis.eq.0.) write(*,*) e2sr(j,i),phase_ratio(j,i,k),Eff_visc,vis
endif
return
end function Eff_visc
