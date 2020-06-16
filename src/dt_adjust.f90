! ---- Choice of the time step  
 
subroutine dt_adjust
use params
include 'precision.inc'

dt_min = min(dt_elastic, dt_maxwell)
dt = max (dt, dt_min )
return
end
