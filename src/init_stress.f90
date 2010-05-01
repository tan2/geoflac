!---------------------------------------------------
! Calculate  HYDROSTATIC prestresses 
!---------------------------------------------------

subroutine init_stress
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

pi = 3.14159265358979323846
degrad = pi/180.

n = 0
do 51 i = 1,nx-1
    rogh = 0.
    do 52 j = 1,nz-1
        iph = iphase(phasez (j,i))
        if(iph.eq.0) goto 52
        tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
        densT = den(iph) * ( 1 - alfa(iph)*tmpr )
        dh1 = cord (j,i  ,2) - cord (j+1,i  ,2)
        dh2 = cord (j,i+1,2) - cord (j+1,i+1,2)
        dh  = 0.5 * (dh1+dh2)
        dPT = densT * g * dh

        dP = dPT * ( 1 - beta(iph)*rogh ) / ( 1 + beta(iph)/2*dPT )

        press = rogh + 0.5*dP
        do ii = 1,4
            stress0(j,i,1,ii) = -press
            stress0(j,i,2,ii) = -press
            stress0(j,i,3,ii) = 0.
            stress0(j,i,4,ii) = -press
        end do
        !write (*,*) -press
        rogh = rogh + dP
52  continue

    ! Calculate pisos (isostatic pressure at the bottom boundary in the left corner!!)
    ! for AUTOMATIC Wrinkler b.c
    ! pisos is the pressure at y = rzbo 
    if(i.eq.nx-1 .and. nyhydro.eq.2) pisos = rogh
51 continue


open(1,file='pisos.rs')
write(1,*) pisos
close (1)

return
end
