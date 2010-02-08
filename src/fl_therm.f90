
! Calculate thermal field in explict form

subroutine fl_therm
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension rhs(mnz+1,mnx+1), flux(mnz,mnx,2,2), area_n(mnz+1,mnx+1), add_source(mnz,mnx)

! real_area = 0.5* (1./area(n,t))
! Calculate Fluxes in every triangle
!      flux (j,i,num_triangle, direction(x,y)
!
!  1 - 3
!  |   |
!  2 - 4
!
!  diagonal / :
!
!   A:        B:
!
!  1---3         1
!  | /         / |
!  2         2---3


ntherm = ntherm+1

! saving old temperature
!$DIR PREFER_PARALLEL
temp0 = temp

!dt_therm = time - time_t
dt_therm = dt
!if( dt_therm.gt.dtmax_therm ) then
!    write(*,*) dt_therm,dtmax_therm
!    call SysMsg('DT_THERM is larger than DTMAX_THERM')
!    stop 37
!    return
!endif

!$OMP Parallel private(i,j,iph,cp_eff,cond_eff,dissip,diff,x1,x2,x3,x4,y1,y2,y3,y4,t1,t2,t3,t4,qs,real_area13)
!$OMP do
do i = 1,nx-1
    do j = 1,nz-1

        iph = iphase(i,j,phasez(j,i))

        ! Calculating effective material properties
        cp_eff = Eff_cp( i,j )
        cond_eff = Eff_conduct( i,j )

        ! if shearh-heating flag is true
        if( ishearh.eq.1 .and. itherm.ne.2 ) then
            dissip = shrheat(j,i)
        else
            dissip = 0
        endif

        ! Additional sources - radiogenic and shear heating
        add_source(j,i) = ( source(j,i) + dissip/den(iph) - 600.*cp_eff*Eff_melt(i,j)) / cp_eff

        ! diffusivity
        diff = cond_eff/den(iph)/cp_eff

        ! Calculate fluxes in two triangles
        x1 = cord (j  ,i  ,1)
        y1 = cord (j  ,i  ,2)
        x2 = cord (j+1,i  ,1)
        y2 = cord (j+1,i  ,2)
        x3 = cord (j  ,i+1,1)
        y3 = cord (j  ,i+1,2)
        x4 = cord (j+1,i+1,1)
        y4 = cord (j+1,i+1,2)
        t1 = temp (j   ,i  )
        t2 = temp (j+1 ,i  )
        t3 = temp (j   ,i+1)
        t4 = temp (j+1 ,i+1)

        ! (1) A element:
        flux(j,i,1,1) = -diff * ( t1*(y2-y3)+t2*(y3-y1)+t3*(y1-y2) ) * area(1,j,i)
        flux(j,i,1,2) = -diff * ( t1*(x3-x2)+t2*(x1-x3)+t3*(x2-x1) ) * area(1,j,i)
 
        ! (2) B element: Interchange of numeration: (1 -> 3,  3 -> 4)
        flux(j,i,2,1) = -diff * ( t3*(y2-y4)+t2*(y4-y3)+t4*(y3-y2) ) * area(2,j,i)
        flux(j,i,2,2) = -diff * ( t3*(x4-x2)+t2*(x3-x4)+t4*(x2-x3) ) * area(2,j,i)

    end do
end do    
!$OMP end do


!$OMP do
do i = 1,nx
    do j = 1,nz

        rhs(j,i) = 0
        area_n(j,i) = 0

        ! Element (j-1,i-1). Triangle B
        if ( j.ne.1 .and. i.ne.1 ) then

            ! side 2-3
            qs = flux(j-1,i-1,2,1) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 flux(j-1,i-1,2,2) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            ! side 3-1
            qs = flux(j-1,i-1,2,1) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 flux(j-1,i-1,2,2) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            real_area13 = 0.5/area(2,j-1,i-1)/3.
            area_n(j,i) = area_n(j,i) + real_area13
            rhs(j,i) = rhs(j,i) + add_source(j-1,i-1)*real_area13

        endif

        ! Element (j-1,i). Triangles A,B
        if ( j.ne.1 .and. i.ne.nx ) then

            ! triangle A
            ! side 1-2
            qs = flux(j-1,i  ,1,1) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 flux(j-1,i  ,1,2) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            ! side 2-3
            qs = flux(j-1,i  ,1,1) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 flux(j-1,i  ,1,2) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            real_area13 = 0.5/area(1,j-1,i  )/3.
            area_n(j,i) = area_n(j,i) + real_area13
            rhs(j,i) = rhs(j,i) + add_source(j-1,i  )*real_area13

            ! triangle B
            ! side 1-2
            qs = flux(j-1,i  ,2,1) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 flux(j-1,i  ,2,2) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            ! side 2-3
            qs = flux(j-1,i  ,2,1) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 flux(j-1,i  ,2,2) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            real_area13 = 0.5/area(2,j-1,i  )/3.
            area_n(j,i) = area_n(j,i) + real_area13
            rhs(j,i) = rhs(j,i) + add_source(j-1,i  )*real_area13

        endif
        
        ! Element (j,i-1). Triangles A,B
        if ( j.ne.nz .and. i.ne.1 ) then

            ! triangle A
            ! side 2-3
            qs = flux(j  ,i-1,1,1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 flux(j  ,i-1,1,2) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            ! side 3-1
            qs = flux(j  ,i-1,1,1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 flux(j  ,i-1,1,2) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            real_area13 = 0.5/area(1,j  ,i-1)/3.
            area_n(j,i) = area_n(j,i) + real_area13
            rhs(j,i) = rhs(j,i) + add_source(j  ,i-1)*real_area13

            ! triangle B
            ! side 1-2
            qs = flux(j  ,i-1,2,1) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 flux(j  ,i-1,2,2) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            ! side 3-1
            qs = flux(j  ,i-1,2,1) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 flux(j  ,i-1,2,2) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            real_area13 = 0.5/area(2,j  ,i-1)/3.
            area_n(j,i) = area_n(j,i) + real_area13
            rhs(j,i) = rhs(j,i) + add_source(j  ,i-1)*real_area13

        endif

        ! Element (j,i). Triangle A
        if ( j.ne.nz .and. i.ne.nx ) then

            ! side 1-2
            qs = flux(j  ,i  ,1,1) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 flux(j  ,i  ,1,2) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            ! side 3-1
            qs = flux(j  ,i  ,1,1) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 flux(j  ,i  ,1,2) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            rhs(j,i) = rhs(j,i) + 0.5*qs

            real_area13 = 0.5/area(1,j  ,i  )/3.
            area_n(j,i) = area_n(j,i) + real_area13
            rhs(j,i) = rhs(j,i) + add_source(j  ,i  )*real_area13

        endif

        ! Update Temperature by Eulerian method 
        temp(j,i) = temp(j,i)+rhs(j,i)*dt_therm/area_n(j,i)
    end do
end do
!$OMP end do

! Boundary conditions (top and bottom)
!$OMP do
do i = 1,nx

    temp(1,i) = t_top

    if( itemp_bc.eq.1 ) then
        temp(nz,i) = bot_bc
    elseif( itemp_bc.eq.2 ) then
        if( i.ne. nx ) then
            cond_eff = Eff_conduct( i,nz-1 )
        else
            cond_eff = Eff_conduct( nx-1,nz-1 )
        endif
        temp(nz,i) = temp(nz-1,i)  +  bot_bc * ( cord(nz-1,i,2)-cord(nz,i,2) ) / cond_eff
    endif

end do
!$OMP end do

! Boundary conditions: dt/dx =0 on left and right  
!$OMP do
do j = 1,nz
    temp(j ,1)  = temp(j,2)
    temp(j, nx) = temp(j,nx-1)
end do
!$OMP end do
!$OMP end parallel

time_t = time


! HOOK
! Intrusions - see user_ab.f90
if( if_intrus.eq.1 ) call MakeIntrusions

return
end 
