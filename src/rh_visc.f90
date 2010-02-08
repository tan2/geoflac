! from Chen and Morgan (1990)
! Uses A value in MPa and but gives viscosity in (Pa * s) 
! Therefore there is a coefficient 1.e+6 

function vis_creep(i,j)
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

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
if( fmelt(j,i) .gt. 0. ) then
    if( fmelt(j,i) .lt. fmelt_crit ) then
        vislog = fmelt(j,i)/fmelt_crit*dlog10(v_min/vis) + dlog10(vis)
        vis = 10.**vislog
    else
        vis = v_min
    endif
endif


! experiments with Peierls law
if (iph .eq. 115 ) then
    PL_sigma = 8.5e+9
    PL_H = 5.4e+5
    PL_A = 5.7e+11
    PL_q = 2.

    s11avg = 0.25*(stress0(1,1,j,i) + stress0(1,2,j,i) + stress0(1,3,j,i) + stress0(1,4,j,i))
    s22avg = 0.25*(stress0(2,1,j,i) + stress0(2,2,j,i) + stress0(2,3,j,i) + stress0(2,4,j,i))
    s12avg = 0.25*(stress0(3,1,j,i) + stress0(3,2,j,i) + stress0(3,3,j,i) + stress0(3,4,j,i))
    sII = sqrt(0.25*(s11avg-s22avg)**2 + s12avg**2)

    if (sII .eq. 0) then
        vs_peierls = 1.e+30
    else
        vs_peierls = sII/PL_A * exp( PL_H/r/(tmpr+273)*(1.-sII/PL_sigma)**PL_q )
    endif

    vis = 1. / (1./vis + 1./vs_peierls)
endif


! Final cut-off
if (vis .lt. v_min) vis = v_min
if (vis .gt. v_max) vis = v_max

vis_creep = vis

return
end
