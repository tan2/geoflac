! -*- F90 -*-

!!! maximum number of ELEMENTS !!!
parameter(mnz=10, mnx=50 )

common /geom_arrays/ area(4,mnz,mnx),phasez(mnz,mnx), &
    dvol(4,mnz,mnx),o_phasez(mnz,mnx)

common /all_arrays/ vel(mnz+1,mnx+1,2),strainr(3,4,mnz,mnx),strain(3,mnz,mnx), &
    force(mnz+1,mnx+1,2),balance(mnz+1,mnx+1,2),rmass(mnz+1,mnx+1),amass(mnz+1,mnx+1), &
    stress0(4,4,mnz,mnx)

common /rheol_arrays/ aps(mnz,mnx),visn(mnz,mnx),e2sr(mnz,mnx),phase_ratio(mnz,mnx,20), ophase_ratio(mnz,mnx,20),irheol_fl(mnz,mnx)

common /temp_arrays/ temp(mnz+1,mnx+1),temp0(mnz+1,mnx+1),source(mnz,mnx),shrheat(mnz,mnx)
    
common /boundc/ bcstress((mnz+mnx)*2,3),bc(mnz+1,mnx+1,2), &
     ncod(mnz+1,mnx+1,2), &
     nopbou((mnz+mnx)*2,4), ncodbou((mnz+mnx)*2,3)

