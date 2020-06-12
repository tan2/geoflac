!================================
! First invariant of strain
function strainI(iz,ix)
use arrays
include 'precision.inc'

strainI = 0.5d0 * ( strain(iz,ix,1) + strain(iz,ix,2) )

return
end function strainI


!================================
! Second invariant of strain
function strainII(iz,ix)
use arrays
include 'precision.inc'

strainII = 0.5d0 * sqrt((strain(iz,ix,1)-strain(iz,ix,2))**2 + 4*strain(iz,ix,3)**2)

return
end function strainII


!================================
! Second invariant of strain rate
function srateII(iz,ix)
use arrays
include 'precision.inc'

s11 = 0.25d0 * (strainr(1,1,iz,ix)+strainr(1,2,iz,ix)+strainr(1,3,iz,ix)+strainr(1,4,iz,ix))
s22 = 0.25d0 * (strainr(2,1,iz,ix)+strainr(2,2,iz,ix)+strainr(2,3,iz,ix)+strainr(2,4,iz,ix))
s12 = 0.25d0 * (strainr(3,1,iz,ix)+strainr(3,2,iz,ix)+strainr(3,3,iz,ix)+strainr(3,4,iz,ix))
srateII = 0.5d0 * sqrt((s11-s22)**2 + 4*s12*s12)

return
end function srateII


!================================
! First invariant of stress (pressure)
function stressI(iz,ix)
use arrays
include 'precision.inc'

s11 = 0.25d0 * (stress0(iz,ix,1,1)+stress0(iz,ix,1,2)+stress0(iz,ix,1,3)+stress0(iz,ix,1,4))
s22 = 0.25d0 * (stress0(iz,ix,2,1)+stress0(iz,ix,2,2)+stress0(iz,ix,2,3)+stress0(iz,ix,2,4))
s33 = 0.25d0 * (stress0(iz,ix,4,1)+stress0(iz,ix,4,2)+stress0(iz,ix,4,3)+stress0(iz,ix,4,4))
stressI = (s11+s22+s33)/3d0

return
end function stressI


!================================
! Second invariant of stress
function stressII(iz,ix)
use arrays
include 'precision.inc'

s11 = 0.25d0 * (stress0(iz,ix,1,1)+stress0(iz,ix,1,2)+stress0(iz,ix,1,3)+stress0(iz,ix,1,4))
s22 = 0.25d0 * (stress0(iz,ix,2,1)+stress0(iz,ix,2,2)+stress0(iz,ix,2,3)+stress0(iz,ix,2,4))
s12 = 0.25d0 * (stress0(iz,ix,3,1)+stress0(iz,ix,3,2)+stress0(iz,ix,3,3)+stress0(iz,ix,3,4))
s33 = 0.25d0 * (stress0(iz,ix,4,1)+stress0(iz,ix,4,2)+stress0(iz,ix,4,3)+stress0(iz,ix,4,4))
stressII = 0.5d0 * sqrt((s11-s22)**2 + 4*s12*s12)

return
end function stressII


!==================================================
! Get the largest eigenvalue and its eigenvector (with its x-component
! fixed to 1) of the deviatoric of a symmetric 2x2 matrix
subroutine eigen2x2(a11, a22, a12, eigval1, eigvecy1)
include 'precision.inc'

adif = 0.5d0 * (a11 - a22)
eigval1 = sqrt(adif**2 + a12**2)
eigvecy1 = (eigval1 - adif) / a12

end subroutine eigen2x2


  !==================================================
subroutine SysMsg( message )
use params
include 'precision.inc'
character* (*) message

!$OMP critical (sysmsg1)
open( 13, file='sys.msg', position="append", action="write" )
write(13, * ) "Loops:", nloop, "Time[Ma]:", time/sec_year/1.d+6
write(13, * ) message
close(13)
!$OMP end critical (sysmsg1)

return
end
