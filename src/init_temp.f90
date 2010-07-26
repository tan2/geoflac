! Initiate temperature profile

subroutine init_temp
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

!allocatable :: a(:,:),b(:)
!allocate( a(nz,nz),b(nz) )
!  Read distribution of temperatures from the dat file
if (irtemp .gt. 0) then
    open( 1, file=tempfile, status='old', err=101 )
    do i = 1, nx
    do j = 1, nz
        read( 1, * ,err=102 ) temp(j,i)
!     if(temp(j,i).ge.1000.) temp(j,i) = 1000.
    enddo
    enddo
    close(1)

   goto 10
    101 call SysMsg('INIT_TEMP: Cannot open file with temperatures initial distrib!')
    stop 21
    102 call SysMsg('INIT_TEMP: Error reading file with temperatures initial distrib!')
    stop 21
endif

! Temperature structure for ridges
! uses setup for viscosity from Alexei
if (iynts.eq.1) then
do i = 1,nx-1
    do j = 1,nz-1
       xc = 0.25*(cord (j,i  ,1) + cord(j+1,i  ,1) + &
                cord (j,i+1,1) + cord(j+1,i+1,1))
       yc = 0.25*(cord (j,i  ,2) + cord(j+1,i  ,2) + &
                cord (j,i+1,2) + cord(j+1,i+1,2))
!       Line
       if (igeotherm .eq.0) geoth = g_y0c
!       Gauss perturbation
       if (igeotherm .eq.1 ) then
       geoth = g_y0c + g_amplitude*exp(-((xc-g_x0)/g_width)**2.)
       endif
!       Linear perturbation
       if (igeotherm .eq.2) then
       if ( abs(g_x0-xc).lt.g_width) geoth = g_y0c+ &
        g_amplitude*(1.-abs(g_x0-xc)/g_width)
       if ( abs(g_x0-xc).ge.g_width) geoth = g_y0c
       endif

! Temperatures
       te0 = tbos
       efold = efoldc

! E-fold depends on x (correction due to lateral change in geotherm)

       if(yc.ge.geoth) then
           temp(j,i)=t_top+((te0-t_top)/geoth)*yc
       else
           temp(j,i)=te0 + ((te0-t_top)/(0.5*geoth))*(yc-geoth)
       endif
       if(temp(j,i).gt.t_bot) temp(j,i) = t_bot
enddo
enddo
do j = 1, nz
   temp(j,nx) = temp(j,nx-2)
enddo
do i = 1, nx
   temp(nz,i) = temp(nz-1,i)
enddo
goto 10
endif
!!  Initial continental geotherm of a given age accross the box with variable age
if (iynts.eq.2) then
cond_c = 2.2
cond_m = 3.3
dens_c = 2700.
dens_m = 3300.
t_boot = 1330.
pi = 3.14159
do n = 1, nzone_age
   tr= dens_c*hs*hr*hr*1.e+6/cond_c*exp(1.-exp(-hc(n)/hr))
      q_m = (t_boot-t_top-tr)/((hc(n)*1000.)/cond_c+((200.e3-(hc(n))*1000.))/cond_m)  
   tm  = t_top + (q_m/cond_c)*hc(n)*1000. + tr
!   write(*,*) rzbo, tr, hs, hr, hc(n), q_m, tm
   age_init = age_1(n)*3.14*1.e+7*1.e+6
   diff_m = cond_m/1000./dens_m
   tau_d = 200.e3*200.e3/(pi*pi*diff_m) 
   do i = ixtb1(n), ixtb2(n)
      do j = 1,nz
! depth in km
         y = -cord(j,i,2)*1.e-3
!  steady state part
         if (y.le.hc(n)) tss = t_top+(q_m/cond_c)*y*1000.+(dens_c*hs*hr*hr*1.e+6/cond_c)*exp(1.-exp(-y/hr))
         if (y.gt.hc(n)) tss = tm + (q_m/cond_m)*1000.*(y-hc(n))

! time-dependent part
         tt = 0.
         pp =-1.
         do k = 1,100
            an = 1.*k
            pp = -pp
            tt = tt +pp/(an)*exp(-an*an*age_init/tau_d)*dsin(pi*k*(200.e3-y*1000.)/(200.e3))
         enddo
      temp(j,i) = tss +2./pi*(t_boot-t_top)*tt
        if(temp(j,i).gt.1330.or.y.gt.200.) temp(j,i)= 1330.
        if (j.eq.1) temp(j,i) = t_top
        t_bot = temp(nz,i)
!       write(*,*) tss,tm,q_m,cond_m,hc(n),y,tt
      enddo
    enddo
enddo 
goto 20
endif      
if (iynts.eq.20) then
cond_c = 2.2
cond_m = 3.3
dens_c = 2700.
dens_m = 3300.
t_boot = 1330.
pi = 3.14159
n=1
   tr= dens_c*hs*hr*hr*1.e+6/cond_c*exp(1.-exp(-hc(n)/hr))
      q_m = (t_boot-t_top-tr)/((hc(n)*1000.)/cond_c+((200.e3-(hc(n))*1000.))/cond_m)  
   tm  = t_top + (q_m/cond_c)*hc(n)*1000. + tr
!   write(*,*) rzbo, tr, hs, hr, hc(n), q_m, tm
   age_init = age_1(n)*3.14*1.e+7*1.e+6 + time 
   diff_m = cond_m/1000./dens_m
   tau_d = 200.e3*200.e3/(pi*pi*diff_m)
   !XXX: temp(:,6:nx) not initialized
   write(*,*) 'Warning: iynts=20 -- temp(:,6:nx) not initialized'
   do i = 1, 5 
      do j = 1,nz
! depth in km
         y = -cord(j,i,2)*1.e-3
!  steady state part
         if (y.le.hc(n)) tss = t_top+(q_m/cond_c)*y*1000.+(dens_c*hs*hr*hr*1.e+6/cond_c)*exp(1.-exp(-y/hr))
         if (y.gt.hc(n)) tss = tm + (q_m/cond_m)*1000.*(y-hc(n))

! time-dependent part
         tt = 0.
         pp =-1.
         do k = 1,100
            an = 1.*k
            pp = -pp
            tt = tt +pp/(an)*exp(-an*an*age_init/tau_d)*dsin(pi*k*(200.e3-y*1000.)/(200.e3))
         enddo
      temp(j,i) = tss +2./pi*(t_boot-t_top)*tt
        if(temp(j,i).gt.1330.or.y.gt.200.) temp(j,i)= 1330.
        if (j.eq.1) temp(j,i) = t_top
!       write(*,*) tss,tm,q_m,cond_m,hc(n),y,tt
      enddo
    enddo
goto 20
endif      

! estimate initial temperature as linear (for first approx. of conductivities)
do j = 1,nz
    temp(j,1:nx) = (t_bot-t_top)/abs(rzbo)*abs(cord(j,1,2)-z0) + t_top
end do

dz = abs(cord(2,1,2)-cord(1,1,2))/1000

!irep = 0
!do while (.true.)
!    a = 0; b = 0;
!    a(1,1) = 1; b(1) = t_top;
!    do j = 2, nz-1
!        a(j,j-1) = cnd(j-1)+4*cnd(j)-cnd(j+1)
!        a(j,j  ) = -8*cnd(j)
!        a(j,j+1) = cnd(j+1)+4*cnd(j)-cnd(j-1)
!        b(j) = -4*dz*dz*htgen(j)
!    enddo
!    a(nz,nz) = 1; b(nz) = t_bot;
!
!    call Gauss(a,nz,b)
!
!    tdiff = 0
!    do j = 1,nz
!        tdiff = max( abs( b(j)-temp(j,1) ), tdiff )
!        temp(j,1:nx) = b(j)
!    end do

!    if( tdiff .lt. 0.1 ) exit

!    irep = irep+1
!    if( irep .gt. 1000 ) then
!        call SysMsg('INIT_TEMP: No convergence !')
!        stop
!    endif

!end do

!deallocate( a,b )

10 continue
20 continue
open( 1, file='temp0.dat' )
do j = 1,nz
    write(1,'(f5.1,1x,f6.1,1x,f6.1,1x,f6.1)') -cord (j,1,2)*1.e-3, temp(j,1)
end do
close(1)


! DISTRIBUTE SOURCES in elements
do j = 1,nz-1
    y = -( cord(j+1,1,2)+cord(j,1,2) )/2 / 1000
    source(j,1:nx-1) = hs*exp(-y/hr)
end do

! Initial quadrilateral temperature perturbation
if( temp_per.gt.0 ) then
    temp(iy1t:iy2t,ix1t:ix2t) = temp(iy1t:iy2t,ix1t:ix2t) + temp_per
endif              


!call RedefineTemp

return
end


function cnd( j )
include 'precision.inc'

cnd = Eff_conduct(j,1)

return
end


function htgen( j )
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

y = - cord(j,1,2)*1.e-3

iph = iphase(j,1)

htgen = den(iph)*hs*exp(-y/hr) * 1.e+6

return
end


!=========================================================
subroutine RedefineTemp
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

write(*,*) 'ATTENTION! Special form of initial temperature distribution !'

return
end

