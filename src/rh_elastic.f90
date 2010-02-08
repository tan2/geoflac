! Linear Elastic Model   (Plane strain)

subroutine elastic(bulkm,rmu,s11,s22,s33,s12,de11,de22,de12,iph)
include 'precision.inc'
include 'params.inc'

c1d3 = 1./3.
c4d3 = 4./3.
c2d3 = 2./3.

a1 = bulkm + c4d3*rmu  
a2 = bulkm - c2d3*rmu  

s11 = s11 + a1*de11 + a2*de22  
s22 = s22 + a2*de11 + a1*de22 
s12 = s12 + 2.*rmu*de12
s33 = s33 + a2*(de11+de22) 
s0 = c1d3 * (s11 + s22 + s33)

return 
end 

!  In  lame coefficients:
!      s11 = s11 + rlam*(de11+de22) + 2.*rmu*de11 
!      s22 = s22 + rlam*(de11+de22) + 2.*rmu*de22
!      s12 = s12 + 2.*rmu*de12
!      s33 = s33 + rlam*(de11+de22)
