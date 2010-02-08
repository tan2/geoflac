!================================
! First invariant of strain
function strainI(i,j)
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

strainI = 0.5 * ( strain(1,j,i) + strain(2,j,i) )

return
end


!================================
! Second invariant of strain
function strainII(i,j)
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

strainII = 0.5 * sqrt((strain(1,j,i)-strain(2,j,i))**2 + 4*strain(3,j,i)**2)

return
end


!================================
! Second invariant of strain rate
function srateII(i,j)
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

s11 = 0.25 * (strainr(1,1,j,i)+strainr(1,2,j,i)+strainr(1,3,j,i)+strainr(1,4,j,i))
s22 = 0.25 * (strainr(2,1,j,i)+strainr(2,2,j,i)+strainr(2,3,j,i)+strainr(2,4,j,i))
s12 = 0.25 * (strainr(3,1,j,i)+strainr(3,2,j,i)+strainr(3,3,j,i)+strainr(3,4,j,i))
srateII = 0.5 * sqrt((s11-s22)**2 + 4*s12*s12)

return
end


!================================
! First invariant of stress (pressure)
function stressI(i,j)
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

s11 = 0.25 * (stress0(1,1,j,i)+stress0(1,2,j,i)+stress0(1,3,j,i)+stress0(1,4,j,i))
s22 = 0.25 * (stress0(2,1,j,i)+stress0(2,2,j,i)+stress0(2,3,j,i)+stress0(2,4,j,i))
s33 = 0.25 * (stress0(4,1,j,i)+stress0(4,2,j,i)+stress0(4,3,j,i)+stress0(4,4,j,i))
stressI = (s11+s22+s33)/3

return
end


!================================
! Second invariant of stress
function stressII(i,j)
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

s11 = 0.25 * (stress0(1,1,j,i)+stress0(1,2,j,i)+stress0(1,3,j,i)+stress0(1,4,j,i))
s22 = 0.25 * (stress0(2,1,j,i)+stress0(2,2,j,i)+stress0(2,3,j,i)+stress0(2,4,j,i))
s12 = 0.25 * (stress0(3,1,j,i)+stress0(3,2,j,i)+stress0(3,3,j,i)+stress0(3,4,j,i))
s33 = 0.25 * (stress0(4,1,j,i)+stress0(4,2,j,i)+stress0(4,3,j,i)+stress0(4,4,j,i))
stressII = 0.5 * sqrt((s11-s22)**2 + 4*s12*s12)

return
end


!==================================================
subroutine SysMsg( message )
include 'precision.inc'
include 'params.inc'
character* (*) message

open( 13, file='sys.msg' )
do while (.TRUE.)
    read( 13, *, end=10 )
end do
10 continue
backspace(1)
write(13, '(i9,1x,i7,5x,9HTime[Ma]:,1x,f5.2)' ) nloop, time/sec_year/1.e+6
write(13, * ) message
close(13)

return
end


!=======================================
FUNCTION ERFCC(X)
include 'precision.inc'

Z=ABS(X)      
T=1./(1.+0.5*Z)
ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+ &
      T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+ &
      T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
IF (X.LT.0.) ERFCC=2.-ERFCC

RETURN
END


!===========================================
! Protected dlog10
function dlog10ab( val )
include 'precision.inc'

if ( val .le. 0. ) then
    dlog10ab = 1.e-30
else
    dlog10ab = dlog10( val )
endif

return
end


!===============================
! Gaussian elemination
!===============================
subroutine Gauss(A,N,B)
include 'precision.inc'
dimension A(N,N),B(N),IPIV(N),INDXR(N),INDXC(N)
   
ipiv = 0      

do i=1,N
    BIG=0.
    do j=1,N
        if( IPIV(j).ne.1 ) then
            do k=1,N
                if( IPIV(k).eq.0 ) then
                    if( abs(A(j,k)).gt.BIG ) then
                        BIG = abs( A(j,k) )
                        IROW=j
                        ICOL=k
                    endif
                else if( IPIV(k).gt.1 ) then
                    pause 'Singular matrix'
                    stop
                endif
            end do
        endif
    end do

    IPIV(ICOL) = IPIV(ICOL) + 1
    if( IROW.ne.ICOL ) then
        do l=1,N
            dum = A(IROW,l)
            A(IROW,l) = A(ICOL,l)
            A(ICOL,l) = dum
        end do
        dum = B(IROW)
        B(IROW)=B(ICOL)
        B(ICOL)=dum
    end if

    INDXR(i) = IROW
    INDXC(i) = ICOL
    if( A(ICOL,ICOL).eq.0. ) then
        pause 'Singular matrix.'
        stop
    endif
    pivinv = 1./A(ICOL,ICOL)
    A(ICOL,ICOL) = 1.
    do l=1,N
        A(ICOL,l) = A(ICOL,l) * pivinv
    end do
    B(ICOL) = B(ICOL) * pivinv
    do ll=1,N
        if( ll.ne.ICOL ) then
            dum = A(ll,ICOL)
            A(ll,ICOL) = 0.
            do l=1,N
                A(ll,l) = A(ll,l) - A(ICOL,l)*dum
            end do
            B(ll) = B(ll) - B(ICOL)*dum
        endif
    end do
end do

do l=N,1,-1
    if( INDXR(l).ne.INDXC(l) ) then
        do k=1,N
            dum = A(k,INDXR(l))
            A(k,INDXR(l)) = A(k,INDXC(l))
            A(k,INDXC(l)) = dum
        end do
    end if
end do


return
end
