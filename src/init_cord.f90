! ---------------------------------------------------------------
! Calculation of coordinates at the start
! ---------------------------------------------------------------

subroutine init_cord
use arrays
use params
include 'precision.inc'

dimension rmesh1(nx+nz)

if (ircoord .gt. 0) then
    ! read the coordinates from file
    open(1, file=coordfile, status='old', err=101)

    do i = 1,nx
        do j = 1,nz
            read(1, *, err=102) x, y
            cord(j,i,1) = x
            cord(j,i,2) = y
        end do
    end do

    close(1)
    go to 200

101 call SysMsg('INIT_CORD: Cannot open file with initial coordinates!')
    stop 21
102 call SysMsg('INIT_CORD: Error reading file with initial coordinates!')
    stop 21
end if

! Check dimensions for the mesh (number of elements and sizes)  
call test_grid
!  X - component 
call mesh1 (x0,x0+rxbo,rmesh1,nzonx,nelz_x,sizez_x) 

do i = 1,nx
    cord (:,i,1) =  rmesh1(i)
end do

!  Z - component
do j = 1,nz
    call mesh1 (z0,z0+rzbo,rmesh1,nzony,nelz_y,sizez_y)
    cord (j,:,2) =  rmesh1(j)
end do

! topo perturbation
do i = 1, inhom
    if (igeom(i).eq.20) then
        ! table mountain
        amp = xinitaps(i)
        cord(1, ix1(i):ix2(i), 2) = cord(1, ix1(i):ix2(i), 2) + amp
    elseif (igeom(i).eq.21) then
        ! trapzoidal mountain
        amp = xinitaps(i)
        do j = ix1(i), ix2(i)-1
            cord(1,j,2) = cord(1,j,2) + amp*real(j-ix1(i))/(ix2(i)-ix1(i))
        enddo
        cord(1, ix2(i):iy1(i), 2) = cord(1, ix2(i):iy1(i), 2) + amp
        do j = iy1(i)+1, iy2(i)
            cord(1,j,2) = cord(1,j,2) + amp*real(iy2(i)-j)/(iy2(i)-iy1(i))
        enddo 
    endif
enddo


200 continue

dhacc(1:nx) = 0.d0
extr_acc(1:nx-1) = 0.d0

return
end



!====================================================================
!------Test for mesh dimensions (check size and num. of elem)

subroutine test_grid
use params
include 'precision.inc'


!- X direction
sizesum=0.d0
nelsum =0
do i = 1,nzonx
    sizesum = sizesum + sizez_x(i)
    nelsum =  nelsum  + nelz_x(i)
end do

if ( abs ((sizesum-1.0d0)) .gt. 1.d-4 ) then
    call SysMsg('INIT_CORD: Sum of zones sizes is not correct.(X-direction)')
    stop 16
endif

if (nelsum .ne. (nx-1)) then
    call SysMsg('INIT_CORD: Nelem in zones is not correct.(X-direction)')
    stop 17
endif


!- Y direction
sizesum=0.d0
nelsum =0
do i = 1,nzony 
    sizesum = sizesum + sizez_y(i)
    nelsum =  nelsum  + nelz_y(i)
end do

if ( abs (sizesum-1.d0) .gt. 1.d-4 ) then
    call SysMsg('INIT_CORD: Sum of zones sizes is not correct.(Y-direction)')
    stop 18
endif
 
if (nelsum .ne. (nz-1)) then
    call SysMsg('INIT_CORD: Nelem in zones is not correct.(Y-direction)')
    stop 19
endif

return
end 


!==========================================================================
! ONE-D mesh
!==========================================================================
subroutine mesh1 (x1,x2,xmesh,nzon,nelz,sizez)
use params
implicit none
double precision :: x1, x2, xmesh(nx+nz), sizez(maxzone)
integer :: nzon, nelz(maxzone), k, i, nel
double precision :: xbeg, xend, elemsize

nel = 1
xmesh(nel) = x1

do k = 1,nzon
    xbeg = xmesh(nel)
    xend = xbeg + sizez(k)*(x2-x1)
    if( k.eq.nzon ) xend = x2
    elemsize = (xend-xbeg)/nelz(k)
    do i = 1,nelz(k)
        nel = nel + 1
        xmesh(nel) = xbeg + i*elemsize
    end do
    continue
end do

return
end
