! -*- F90 -*-

!!! maximum number of ELEMENTS !!!
parameter(mnz=80, mnx=400)

common /geom_arrays/ iphase(mnz,mnx), nphase_counter(mnz,mnx,20)

common /all_arrays/ strainr(3,4,mnz,mnx)

common /rheol_arrays/ aps(mnz,mnx),visn(mnz,mnx),e2sr(mnz,mnx),phase_ratio(mnz,mnx,20),irheol_fl(mnz,mnx)

common /temp_arrays/ temp0(mnz+1,mnx+1),source(mnz,mnx),shrheat(mnz,mnx)
    
common /boundc/ bcstress((mnz+mnx)*2,3), &
     nopbou((mnz+mnx)*2,4), ncodbou((mnz+mnx)*2,3)

