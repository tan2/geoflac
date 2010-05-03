! -*- F90 -*-

program DREZINA
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'omp_lib.h'


real*4 secnds,time0

! Area
open( 33, file='area.dat' )
open( 333, file='output.asc' )

time0 = secnds(0.0)

! seconds in a year
sec_year = 3.1558e+7

! Read task parameters
call read_params()
call allocate_arrays(nz, nx)

! Try to read save-file contents. If file exist - restart, othewise - new start
open(1,file='_contents.rs',status='old',err=10)

irestart = 1
close(1)
goto 20

10 irestart = 0

20 continue

if ( irestart .eq. 1 ) then  !file exists - restart
    call rsflac
    if( dtout_screen .ne. 0 ) then
        write(6,*) 'you CONTINUE from  ', nloop+1, ' step'
    else
        call SysMsg('you CONTINUE the execution')
    endif
else ! file does not exist - new start
    if( dtout_screen .ne. 0 ) then
        write(6,*) 'you have NEW start conditions'
        write(333,*) 'you have NEW start conditions'
    else
        call SysMsg('you have NEW start conditions')
    endif
    call setflac
    ! Output of initial configuration
    call outflac
 !   if (iint_marker.eq.1) call outmarker
end if


#ifdef USE_CUDA
! copy variables to C and CUDA
call cu_copy_param( &
     irheol, visc, &
     den, alfa, beta, &
     pln, acoef, eactiv, &
     rl, rm, coha, cohdisp, &
     phimean, phidisp, psia, &
     conduct, cp, &
     ts, tl, tk, fk, &
     g, pisos, drosub, &
     rzbo, demf, &
     sec_year, ynstressbc, &
     dt_scale, frac, fracm, &
     strain_inert, vbc, &
     lphase, &
     nx, nz, nyhydro, iphsub, &
     n_boff_cutoff, i_prestress, &
     iint_marker, nphasl, idt_scale &
     )

if(.false.) then
    ! Copy and paste the arguments to here
    ! Otherwise, the compiler won't be able to detect typo with debug build
    write(*,*) &
         irheol, visc, &
         den, alfa, beta, &
         pln, acoef, eactiv, &
         rl, rm, coha, cohdisp, &
         phimean, phidisp, psia, &
         conduct, cp, &
         ts, tl, tk, fk, &
         g, pisos, drosub, &
         rzbo, demf, &
         sec_year, ynstressbc, &
         dt_scale, frac, fracm, &
         strain_inert, vbc, &
         lphase, &
         nx, nz, nyhydro, iphsub, &
         n_boff_cutoff, i_prestress, &
         iint_marker, nphasl, idt_scale
endif

#endif


!      ****************** running ********************************
ireset = 1
dtacc_screen = 0
dtacc_file = 0
dtacc_save = 0
dtacc_tracer = 0
nrec = 1 
i_search = 0
!do index_nobody_would_use=1,1
do while( time .lt. time_max )
  nloop = nloop + 1
!write(*,*) dt
  if( dtout_screen .ne. 0 ) then
    if( dtacc_screen .gt. dtout_screen ) then
       write(*,'(I7,A,F6.3,A,F6.2,A,F6.1,A)') nloop,'''s step. Time[My]=', time/sec_year/1.e+6, &
                ', dt=', dt/sec_year, ',  elapsed-', secnds(time0)/60, ' min'
       write(333,'(I7,A,F6.3,A,F6.2,A,F6.1,A)') nloop,'''s step. Time[My]=', time/sec_year/1.e+6, &
                ', dt=', dt/sec_year, ',  elapsed-', secnds(time0)/60, ' min'

       ! Forces at the boundaries
       if( io_forc.eq.1 ) then
         force_l=0.
         force_r=0.
         do j = 1,nz-1
           sxx = 0.25 * (stress0(j,1,1,1)+stress0(j,1,1,2)+stress0(j,1,1,3)+stress0(j,1,1,4) )
           sxxd = sxx-stressI(j,1)
           dl = cord(j+1,1,2)-cord(j,1,2)
           force_l = force_l+abs(sxxd)*abs(dl)

           sxx = 0.25 * (stress0(j,nx-1,1,1)+stress0(j,nx-1,1,2)+stress0(j,nx-1,1,3)+stress0(j,nx-1,1,4) )
           sxxd = sxx-stressI(j,nx-1)
           dl = cord(j+1,nx-1,2)-cord(j,nx-1,2)
           force_r = force_r+abs(sxxd)*abs(dl)
         end do
         open (1,file='forc.0',access='direct',form='formatted',recl=28)
         write (1,'(f6.2,1x,e10.2,1x,e10.2)',rec=nrec) time/sec_year/1.e6, force_l, force_r
         nrec = nrec + 1
         close (1)
       endif

       dtacc_screen = 0
     endif
  endif

  ! FLAC
  call flac

  if( ireset.eq.1 ) ireset = 0

  ! Remeshing
  if( ny_rem.eq.1 .and. itherm.ne.2 ) then
    if( itest_mesh() .eq. 1 ) then
      call fl_therm
      if(iynts.eq.1) call init_temp
      ! If there are markers recalculate their x,y global coordinate and assign them aps, eII, press, temp
      if(iint_marker.eq.1) then
        call bar2euler
        call elem2marker
      endif
      call re_mesh
      ! If markers are present recalculate a1,a2 local coordinates and assign elements with phase ratio vector
      if (iint_marker.eq.1) then
        call lpeuler2bar
        call marker2elem
      endif
      ireset = 1
      ! call outmarker
    endif
  endif

    ! OUTPUT  
  if( dtout_file .ne. 0 ) then 
    if( dtacc_file .gt. dtout_file ) then
      call outflac
      dtacc_file = 0
    endif
  endif
  if ( iint_tracer.eq.1) then
    if (dtacc_tracer.gt.dt_outtracer) then
      call outtracer
      dtacc_tracer = 0
    endif
  endif
  ! SAVING
  if( dtsave_file .ne. 0 ) then 
    if( dtacc_save .gt. dtsave_file ) then
      call saveflac
      dtacc_save = 0
    endif
  endif

  ! Area
  if( mod(nloop,1000) .eq. 0 ) then
    area_diff = total_area(0)/abs(rzbo*rxbo) - 1
    !write( *,'(i6,1x,e9.2,1x,e9.2,1x,e9.2)' ) nloop, area_diff, devmax, dvmax
    write(33,'(i6,1x,e9.2,1x,e9.2,1x,e9.2)' ) nloop, area_diff, devmax, dvmax
    devmax = 0; dvmax = 0;
    !call flush(33)
  endif

  time = time + dt

  dtacc_screen = dtacc_screen + dt
  dtacc_file   = dtacc_file   + dt
  dtacc_save   = dtacc_save   + dt
  dtacc_tracer = dtacc_tracer + dt 
end do

! Area
close(33)
close(333)

call SysMsg('Congratulations !')
end program DREZINA
