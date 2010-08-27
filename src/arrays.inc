! -*- F90 -*-

!!! maximum number of ELEMENTS !!!
parameter(mnz=100, mnx=400)

common /phase_arrays/ iphase(mnz,mnx), nphase_counter(20,mnz,mnx), phase_ratio(20,mnz,mnx)

common /all_arrays/ strainr(3,4,mnz,mnx)

common /rheol_arrays/ aps(mnz,mnx),visn(mnz,mnx),e2sr(mnz,mnx),irheol_fl(mnz,mnx)

common /temp_arrays/ temp0(mnz+1,mnx+1),source(mnz,mnx),shrheat(mnz,mnx)
    
common /boundc/ bcstress((mnz+mnx)*2,3), &
     nopbou((mnz+mnx)*2,4), ncodbou((mnz+mnx)*2,3)

