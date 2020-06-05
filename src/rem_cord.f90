subroutine rem_cord(cordo)
use arrays
use params
implicit none
double precision :: cordo(nz,nx,2), rmesh1(nx)
integer :: i, j, ii
double precision :: xl, xr, xx, zcorr, zl, zr, zz, total_area
logical, parameter :: do_volcorrection = .false.

! X - coordinate
do j = 1, nz

    if( mode_rem.eq.1 ) then
        xl = cord(j,1 ,1)
        xr = cord(j,nx,1)
    elseif( mode_rem.eq.11 .OR. mode_rem.eq.3 ) then
        xl = x0
        xr = x0 + rxbo
    endif

    call mesh1( xl,xr,rmesh1,nzonx,nelz_x,sizez_x )

    do i = 1, nx
        cord(j,i,1) = rmesh1(i)
    end do

end do


! Z-coordinate correction for volume change
if( mode_rem.eq.1 .and. do_volcorrection ) then
    zcorr = -( total_area(0)-rzbo*rxbo ) / abs(xr-xl)
else
    zcorr = 0
endif
 

! Z - coordinate
do i = 1, nx

    !  Top and bottom Z coordinates by interpolation from an old grid
    do j = 1, nz, nz-1
        xx = cord(j,i,1)
        if ( xx .le. cordo(j,1,1) ) then
            cord(j,i,2) = cordo(j,1,2)
        elseif ( xx .ge. cordo(j,nx,1) ) then
            cord(j,i,2) = cordo(j,nx,2)
        else
            do ii = 1,nx-1
                xl = cordo(j,ii,1)
                xr = cordo(j,ii+1,1)
                if (xx.ge.xl .and. xx.le.xr) then
                    zl = cordo(j,ii,2)
                    zr = cordo(j,ii+1,2)
                    zz = zl + (xx-xl)*(zr-zl)/(xr-xl)
                    cord(j,i,2) = zz + zcorr
                    exit
                endif
            end do
        endif
    end do
    
    ! For mode_rem=3 bottom is always fixed
    if( mode_rem .eq. 3.or.mode_rem.eq.1 ) cord(nz,i,2) = z0 + rzbo

    ! Creating Mesh inside of the boundaries
    call mesh1 (cord(1,i,2),cord(nz,i,2),rmesh1,nzony,nelz_y,sizez_y)
    do j = 1, nz
        cord(j,i,2) = rmesh1(j)
    end do

end do


return
end
