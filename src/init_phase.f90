!--------------------------------------------------------------
! Initialization of phases in each zone 
!--------------------------------------------------------------

subroutine init_phase
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


!character  phasedat*20

pi = 3.14159265358979323846
degrad = pi/180.




! Main phase
iphase = mphase

! Other phases in horizontal layers 
do k = 1,nphasl
    do i = 1,nx-1
        do j = ltop(k),lbottom(k)
            iphase(j,i) = lphase(k)
            !XXX: hard-coded weak zone layer? remove?
            if(i.ge.170) iphase(j,i) = 8
!making everything west of the continent oceanic mantle....
!!! THAT IS THE LINES FOR THE CONTINENT ON THE RIGHT!!!
!if (i.le.218) iphase(j,i) = 4
!tapering off the continental crust at teh boundary.....
!if (i.ge.219.and.i.le.223.and.j.ge.6) iphase(j,i) = 4
!if (i.ge.224.and.i.le.229.and.j.ge.6) iphase(j,i) = 4
!if (i.ge.230.and.i.le.234.and.j.ge.6) iphase(j,i) = 4
!if (i.ge.235.and.i.le.220.and.j.ge.6) iphase(j,i) = 4

!if (i.ge.241.and.i.le.220.and.j.ge.6) iphase(j,i) = 4
!if (i.ge.221.and.i.le.230.and.j.ge.9) iphase(j,i) = 4
!if (i.ge.231.and.i.le.240.and.j.ge.11) iphase(j,i) = 4
!if (i.ge.241.and.i.le.250.and.j.ge.13) iphase(j,i) = 4
!if (i.ge.188.and.i.le.193.and.j.ge.10) iphase(j,i) = 4
!if (i.ge.194.and.i.le.199.and.j.ge.12) iphase(j,i) = 4
!if (i.ge.200.and.i.le.205.and.j.ge.14) iphase(j,i) = 4
!if (i.ge.206.and.i.le.210.and.j.ge.16) iphase(j,i) = 4

!if (i.le.168.and.j.le.1) iphase(j,i) = 3

!if (i.le.169.and.j.eq.2) iphase(j,i) = 3
!if (i.ge.168.and.j.le.2) iphase(j,i) = 10
!if (i.ge.168.and.i.le.170.and.j.eq.3) iphase(j,i) = 6
!if (i.ge.171.and.i.le.173.and.j.le.4.and.j.ge.3) iphase(j,i) = 6
!if (i.ge.174.and.i.le.218.and.j.le.5.and.j.ge.3) iphase(j,i) = 6

        end do
    end do
end do

!XXX: hard-coded phase layers
nphasl = 12
lphase(1) = 2
lphase(2) = 3
lphase(3) = 4 
lphase(4) = 6
lphase(5) = 7 
lphase(6) = 8 
lphase(7) = 9 
lphase(8) = 10 
lphase(9) = 11 
lphase(10) = 12 
lphase(11) = 14
lphase(12) = 15

!  Read distribution of the phases from the dat file
if (irphase .gt. 0) then
open(12,file='phasedat.inp')
read(12,*) nphasl
do 333 k=1,nphasl
read(12,*) lphase(k)
333 continue
do 332 i=1,nx-1
do 332 j=1,nz-1
!write(*,*) nx,nz
read(12,*) ii,jj,iphase(j,i)
!XXX: hard-coded phase; remove?
if (j.eq.1.and.i.gt.65) iphase(j,i) = 3
if (j.eq.2.and.i.gt.64) iphase(j,i) = 3
if (j.eq.3.and.i.gt.63) iphase(j,i) = 3
332  continue
!XXX: hard-coded phase number; remove?
nphasl = 3
lphase(3) = 3
close(12)

endif                        

! Case with iynts = 2 for continental and collision
if (iynts.eq.10) then
nphasl = 2 
lphase(nphasl-1) = iph_col1(1)
lphase(nphasl) = iph_col2(1)
    do n = 1, nzone_age
       if(n.gt.1) ixtb1(n) = ixtb1(n)-1
       do i = ixtb1(n),ixtb2(n)-1
       do j = 1, nz-1
          y = -cord(j,i,2)*1.e-3
          if (y.lt.hc(n)) iphase(j,i) = iph_col1(n)   
          if (y.ge.hc(n)) iphase(j,i) = iph_col2(n)   
       enddo
       enddo
!   smooth zone in between
!
!       if(iph_col1(n).ne.iph_col1(n-1)) then
!       nphasl = nphasl +1
!       lphase(nphasl) = iph_col1(n)
!       endif
!       if(iph_col2(n).ne.iph_col2(n-1)) then
!       nphasl = nphasl +1
!       lphase(nphasl) = iph_col2(n)
!       endif
    enddo
endif
do i = 1,nphasl
   write(*,*) i,lphase(i)
   write(333,*) i,lphase(i)
enddo
!   Put different rheologies for inclusions 
do i = 1,inhom
    ! Rectangular shape:
    if (igeom(i) .eq.0) then
        do j = ix1(i),ix2(i)
            do k = iy1(i),iy2(i)
                iphase(k,j) = inphase(i)
                aps(k,j)=xinitaps
            end do
        end do
    endif

    ! Gauss shape:
    if (igeom(i).eq.1.or.igeom(i).eq.2) then
        ! symmetric case:
        if (igeom(i).eq.1) then
            ixc  = (ix1(i)+ix2(i))/2  
            iwidth = (ix2(i)-ix1(i))
        else
            ixc    = ix1(i)
            iwidth = ix2(i)-ix1(i)  
        endif
 
        iamp = iy2(i)-iy1(i)
  
        do j = ix1(i),ix2(i)
            itop = itop_geom(j,ixc,iwidth,iamp) 
            do k = iy1(i),iy2(i)
                if (k .ge. (iy2(i)-itop)) then 
                    iphase(k,j) = inphase(i)
                endif
            end do
        end do
    endif

    ! weak zone at 45 degree
    if (igeom (i) .eq.3) then
        do j = ix1(i),ix2(i)
            k = int(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            iphase(k,j) = inphase(i)
        end do
    endif
    
    ! Weak zone in accumulated plastic strain at 45 degree        
    if (igeom (i).eq.4) then
        do j =ix1(i),ix2(i)
            k = int(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            aps(k,j)=xinitaps
            iphase(k,j) = inphase(i)
            !this do loop allows for a less dense Luzon arc by using phase 4
            do ii=1,k-1
              if(iphase(ii,j).eq.8)then
                iphase(ii,j) = 4
              endif
            enddo
        end do
    endif
end do


! Check if viscous rheology present
ivis_present = 0
do i = 1,nx-1
    do j = 1, nz-1
        iph = iphase(j,i)
        if( irheol(iph).eq.3 .or. irheol(iph).ge.11 ) ivis_present = 1
    end do
end do
    

return
end


!==========================================================
! Gauss perturbation at the top of heterogenity 
function itop_geom(j,ixc,iwidth,iamp) 
    
    itop_geom = iamp*exp(-(float(j-ixc)/(0.5*float(iwidth)))**2.)

return
end 


