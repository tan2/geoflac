! ---- Choice of the time step  
 
subroutine dt_adjust
use params
use arrays
implicit none

double precision :: dt_min

!$ACC serial
dt_min = min(dt_elastic, dt_maxwell)
dt = max (dt, dt_min )
!$ACC end serial
return
end
