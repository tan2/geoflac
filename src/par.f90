! -*- F90 -*-

program DREZINA
use arrays
use params
use marker_data

character*200 inputfile
real*4 secnds,time0
integer :: narg, iargc, j, itest_mesh, irestart
double precision :: dtacc_file, dtacc_save, dtacc_screen, dtacc_tracer

narg = iargc()
if(narg /= 1) then
    write(*,*) 'usage: flac inputfile'
    stop 1
endif
call getarg(1, inputfile)

open( 333, file='output.asc' )

time0 = secnds(0.0)

! Read task parameters
call read_params(inputfile)
call allocate_arrays(nz, nx)
call allocate_markers(nz, nx)

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
        write(6,*) 'you CONTINUE from  ', nloop, ' step'
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
    if (iint_marker.eq.1) call outmarker
end if


!      ****************** running ********************************
dtacc_screen = 0
dtacc_file = 0
dtacc_save = 0
dtacc_tracer = 0

do while( time .le. time_max )
  if( dtout_screen .ne. 0 ) then
    if( dtacc_screen .gt. dtout_screen ) then
       write(*,'(I10,A,F7.3,A,F8.1,A)') nloop,'''s step. Time[My]=', time/sec_year/1.d+6, &
                ',  elapsed sec-', secnds(time0)
       write(333,'(I10,A,F7.3,A,F8.1,A)') nloop,'''s step. Time[My]=', time/sec_year/1.d+6, &
                ',  elapsed sec-', secnds(time0)
       dtacc_screen = 0
    endif
  endif

  do j = 1, ntest_rem
    ! FLAC
    call flac

    nloop = nloop + 1
    time = time + dt

    dtacc_screen = dtacc_screen + dt
    dtacc_file   = dtacc_file   + dt
    dtacc_save   = dtacc_save   + dt
    dtacc_tracer = dtacc_tracer + dt
  end do

  ! Remeshing
  if( ny_rem.eq.1 .and. itherm.ne.2 ) then
    if( itest_mesh() .eq. 1 ) then
      ! If there are markers recalculate their x,y global coordinate and assign them aps, eII, press, temp
      if(iint_marker.eq.1) then
        call bar2euler
      endif
      call remesh
      ! If markers are present recalculate a1,a2 local coordinates and assign elements with phase ratio vector
      if (iint_marker.eq.1) then
        call lpeuler2bar
        call marker2elem
      endif
    endif
  endif

    ! OUTPUT  
  if( dtout_file .ne. 0 ) then 
    if( dtacc_file .gt. dtout_file ) then
      call outflac
      if (iint_marker.eq.1) call outmarker
      dtacc_file = 0
    endif
  endif
  if ( iint_tracer.eq.1) then
    if (dtacc_tracer.gt.dt_outtracer) then
      call outtracer
      dtacc_tracer = 0
    endif
  endif
  ! SAVING the restart information
  if( dtsave_file .ne. 0 ) then 
    if( dtacc_save .gt. dtsave_file ) then
      call saveflac
      dtacc_save = 0
    endif
  endif

end do

! SAVING the restart information of last step
call saveflac

close(333)

call SysMsg('Congratulations !')
end program DREZINA
