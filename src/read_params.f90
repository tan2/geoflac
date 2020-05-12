! Reads problem parameters from input file

subroutine read_params(inputfile)
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

character*200 inputfile

iu =4 
open( iu, file=inputfile )

! MESH
call AdvanceToNextInputLine( 4 )
read(4,*) nx, nz
open(11,file='nxnz.0')
write(11,*) nx, nz
close(11)

if((nx.gt.mnx) .or. (nz.gt.mnz)) then
    write(*,*) '# of elements exceed maximum. Increase mnx and mnz in "arrays.inc".'
    stop 1
endif

nq = nx*nz
nx = nx+1
nz = nz+1
nmarkers = 9*(nz-1)*(nx-1)
call AdvanceToNextInputLine( 4 )
read(4,*) x0, z0
call AdvanceToNextInputLine( 4 )
read(4,*) rxbo, rzbo
call AdvanceToNextInputLine( 4 )
read(4,*) ircoord, coordfile
if (ircoord .gt. 0) then
    ! ignore next two lines
    call AdvanceToNextInputLine( 4 )
    read(4,*) nzonx
    if (nzonx .ne. 0) stop 'ircoord is set, but nzonx is not 0!'
    call AdvanceToNextInputLine( 4 )
    read(4,*) nzony
    if (nzony .ne. 0) stop 'ircoord is set, but nzony is not 0!'
    go to 177
endif
call AdvanceToNextInputLine( 4 )
read(4,*) nzonx
if (nzonx .eq.0) then
    nzonx = 1
    nelz_x(1) = nx - 1
    sizez_x(1) = 1.
    go to 166
endif
do i = 1, nzonx
    call AdvanceToNextInputLine( 4 )
    read(4,*) nelz_x(i), sizez_x(i)
end do
166 continue

call AdvanceToNextInputLine( 4 )
read(4,*) nzony
if (nzony .eq.0) then
    nzony = 1
    nelz_y(1)    = nz - 1
    sizez_y(1)   = 1.
    go to 177
endif
do i = 1,nzony
    sizez_y(i) = 1.
    call AdvanceToNextInputLine( 4 )
    read(4,*) nelz_y(i), sizez_y(i)
end do
177 continue
call AdvanceToNextInputLine( 4 )
read(4,*) iint_marker, iint_tracer
call AdvanceToNextInputLine( 4 )
read(4,*) nzone_marker
call AdvanceToNextInputLine( 4 )
  do 253 i = 1, nzone_marker
        read(iu,*) imx1(i),imy1(i),imx2(i),imy2(i)
253 continue 
call AdvanceToNextInputLine( 4 )
read(4,*) nzone_tracer, dt_outtracer
dt_outtracer = dt_outtracer * 1000 * sec_year
call AdvanceToNextInputLine( 4 )
  do 233 i = 1, nzone_tracer
        read(iu,*) itx1(i),ity1(i),itx2(i),ity2(i)
233 continue 

! MECHANICAL CONDITIONS
call AdvanceToNextInputLine( 4 )
read(4,*) ynstressbc,ydrsides
call AdvanceToNextInputLine( 4 )
read(4,*) nofbc
call AdvanceToNextInputLine( 4 )
        do 21 i = 1,nofbc
        read(iu,*) nofside(i),nbc1(i),nbc2(i),nbc(i),  &
        bca(i),bcb(i),bcc(i),  &
        bcd(i),bce(i),bcf(i),bcg(i),bch(i),bci(i)
 21     continue

call AdvanceToNextInputLine( 4 )
! hydrostatic pressure applied at the bottom
call AdvanceToNextInputLine( 4 )
read(4,*) nyhydro,pisos,iphsub,drosub,damp_vis
! gravity
call AdvanceToNextInputLine( 4 )
read(4,*) g



! THERMAL CONDITIONS
call AdvanceToNextInputLine( 4 )
read(4,*) i_prestress
call AdvanceToNextInputLine( 4 )
read(4,*) itherm

! assign kinematic velocity
if (itherm.eq.3) then
  do i=1,nofbc
    if (nbc(i).eq.10 .or. nbc(i).eq.30) then
      ind = (nbc(i)+10) / 20
      vel_init(ind) = bca(i)
      do j=1,nofbc
        if (nbc(j).eq.10 .or. nbc(j).eq.30) then
          ind2 = (nbc(j)+10) / 20
          if (ind.eq.ind2) then
            if (vel_init(ind).ne.bca(j)) then
                stop 'boundary velocities should be the same when using itherm=3'
            end if
          end if
        end if
      end do
    end if
  end do
end if

call AdvanceToNextInputLine( 4 )
read(4,*) istress_therm
call AdvanceToNextInputLine( 4 )             ! thermal stresses
read (4,*) ishearh                           ! shear heating
call AdvanceToNextInputLine( 4 )
read (4,*) t_top
call AdvanceToNextInputLine( 4 )
read (4,*) t_bot  
call AdvanceToNextInputLine( 4 )
read (4,*) hs, hr 
! boundary conditions at the bottom (1-T,2-Flux) 
call AdvanceToNextInputLine( 4 )
read (4,*) itemp_bc, bot_bc
if( itemp_bc.eq.2 ) bot_bc = bot_bc/1000  ! convert in W/m3
! temperature pertrubation (rectangular)
call AdvanceToNextInputLine( 4 )
read (4,*) temp_per, ix1t, ix2t, iy1t, iy2t
! Predefined distributions
call AdvanceToNextInputLine( 4 )
read(4,*) irtemp
if ( irtemp .gt. 0 ) then
    call AdvanceToNextInputLine( 4 )
    read(4,*) tempfile
else
    call AdvanceToNextInputLine( 4 )
    read(4,*)
endif
! time scale
call AdvanceToNextInputLine( 4 )
read (4,*) time_scale
! temp structure
call AdvanceToNextInputLine( 4 )
read (4,*) iynts,tbos
call AdvanceToNextInputLine( 4 )
read (4,*) iax1,iay1,ibx1,iby1,icx1,icy1,idx1,idy1 
call AdvanceToNextInputLine( 4 )
read (4,*) nzone_age 
call AdvanceToNextInputLine( 4 )
do i = 1, nzone_age
      read (4,*) age_1(i),hc1(i),hc2(i),hc3(i),hc4(i),iph_col1(i),iph_col2(i), &
           iph_col3(i),iph_col4(i),iph_col5(i),ixtb1(i),ixtb2(i)
enddo
! check smooth nzone_age
do i = 1, nzone_age-1
    iph_col_trans(i) = 0
    if (ixtb2(i) == -1) then
        if (ixtb1(i+1) /= -1) then
            print *, 'Error: ixtb1 is not -1 at', i+1, 'column!'
        endif
        if (iph_col1(i) /= iph_col1(i+1)) then
            print *, 'Error: iph_col1 at', i, i+1, 'columns are not equal!'
        endif
        if (iph_col2(i) /= iph_col2(i+1)) then
            print *, 'Error: iph_col2 at', i, i+1, 'columns are not equal!'
        endif
        if (iph_col3(i) /= iph_col3(i+1)) then
            print *, 'Error: iph_col3 at', i, i+1, 'columns are not equal!'
        endif
        if (iph_col4(i) /= iph_col4(i+1)) then
            print *, 'Error: iph_col4 at', i, i+1, 'columns are not equal!'
        endif
        if (iph_col5(i) /= iph_col5(i+1)) then
            print *, 'Error: iph_col5 at', i, i+1, 'columns are not equal!'
        endif
        ! flag indicates smotth transition
        iph_col_trans(i) = 1
        ! remove '-1' for
        ixtb2(i) = ixtb2(i+1)
        ixtb1(i+1) = ixtb1(i)
    endif
enddo

! RHEOLOGY
call AdvanceToNextInputLine( 4 )
read(4,*) nphase
do i = 1, nphase 
    call AdvanceToNextInputLine( 4 )
    read(4,*) irheol(i),visc(i),den(i),alfa(i),beta(i),pln(i),acoef(i),eactiv(i),rl(i),rm(i), &
         plstrain1(i),plstrain2(i),fric1(i),fric2(i),cohesion1(i),cohesion2(i), &
         dilat1(i),dilat2(i), &
         conduct(i),cp(i),ts(i),tl(i),tk(i),fk(i)
    print *, i
end do
! Flag to take initial phase distribution from a file
call AdvanceToNextInputLine( 4 )
read(4,*) irphase
!write(*,*) irphase

if ( irphase .gt. 0 ) then
    call AdvanceToNextInputLine( 4 )
    read(4,*) phasefile
else
    call AdvanceToNextInputLine( 4 )
    read(4,*)
endif
! main phase
call AdvanceToNextInputLine( 4 )
read(4,*) mphase
! number of horizontal layers
call AdvanceToNextInputLine( 4 )
read(4,*) nphasl
! layers
do i = 1, nphasl
    call AdvanceToNextInputLine( 4 )
    read(4,*) ltop(i), lbottom(i), lphase(i)
end do
! inclusions
call AdvanceToNextInputLine( 4 )
read(4,*) inhom
if( inhom .gt. 50 ) then
    call SysMsg('Read_params: Increase arrays for inhomogenities')
    stop 26
endif
do i = 1, inhom
    call AdvanceToNextInputLine( 4 )
    read(4,*) ix1(i), ix2(i), iy1(i), iy2(i), inphase(i), igeom(i), xinitaps(i)
end do
! Tension cut-off
call AdvanceToNextInputLine( 4 )
read(4,*) ten_off
!linear healing parameter
call AdvanceToNextInputLine( 4 )
read(4,*) tau_heal
! viscosity limits
call AdvanceToNextInputLine( 4 )
read(4,*) v_min, v_max, ivis_shape,efoldc
call AdvanceToNextInputLine( 4 )
read(4,*)igeotherm,g_x0,g_y0c, g_amplitude,g_width 
call AdvanceToNextInputLine( 4 )
read(4,*) ny_inject, nelem_inject, rate_inject 

! THERMOCHRONOLOGY
call AdvanceToNextInputLine( 4 )
read(4,*) ithermochron
call AdvanceToNextInputLine( 4 )
read(4,*) chron_file
call AdvanceToNextInputLine( 4 )
read(4,*) nchron
do i = 1, nchron
    call AdvanceToNextInputLine( 4 )
    read(4,*) chron_name(i), nchron_fpair(i)
end do
open(11,file='chron.0')
do i=1,nchron
  write(11,*) chron_name(i)
end do
close(11)

! REMESHING
call AdvanceToNextInputLine( 4 )
read(4,*)  ny_rem, mode_rem, ntest_rem, angle_rem
if ( mode_rem.ne.1 .and. mode_rem.ne.11 .and. mode_rem.ne.3 ) then
    call SysMsg('Illegal remeshing mode! Allowable - 1, 3 or 11')
    stop
endif
! dx_rem - remeshing criteria for mode_rem=11
call AdvanceToNextInputLine( 4 )
read(4,*)  dx_rem
! Surface process
call AdvanceToNextInputLine( 4 )
read(4,*)  itopodiff
! diffusion of topography
call AdvanceToNextInputLine( 4 )
read(4,*)  topo_kappa, fac_kappa

! kinematic erosion
call AdvanceToNextInputLine( 4 )
read(4,*) ero_middle, ero_range
do i = 1, 2
  call AdvanceToNextInputLine( 4 )
  read(4,*) nero_rate_sect(i)
  do j = 1, nero_rate_sect(i)
    call AdvanceToNextInputLine( 4 )
    read(4,*) ero_rate(i,1,j), ero_rate(i,2,j), ero_duration(i,j)
  end do
end do

ero_rate = ero_rate / (1.d3*sec_year)
ero_duration = ero_duration * 1.d6*sec_year

d1 = ero_duration(1,1) * ( (ero_rate(1,1,1) + ero_rate(2,1,1))/2. + (ero_rate(1,2,1) + ero_rate(2,2,1))/2. )
d1 = d1 + ero_duration(1,2) * ( ( ero_rate(1,1,2) + ero_rate(2,1,2) )/2. + ( ero_rate(1,2,2) + ero_rate(2,2,2) )/2. )
! the ratio of uplifting rate for  generating alltitude
ratiok = (4400. / d1) + 1.

! PROCESS CONTROL
! inertial mass scaling
call AdvanceToNextInputLine( 4 )
read(4,*)  idt_scale
call AdvanceToNextInputLine( 4 )
read(4,*)  dt_scale, strain_inert
call AdvanceToNextInputLine( 4 )
read(4,*)  i_rey,xReyn 
call AdvanceToNextInputLine( 4 )
read(4,*)  ifreq_rmasses
call AdvanceToNextInputLine( 4 )
read(4,*)  ifreq_imasses
call AdvanceToNextInputLine( 4 )
read(4,*)  ifreq_visc
call AdvanceToNextInputLine( 4 )
read(4,*)  ifreq_avgsr
! acceleration parameters
call AdvanceToNextInputLine( 4 )
read(4,*) amul, ratl, ratu
call AdvanceToNextInputLine( 4 )
read(4,*) frac, fracm
call AdvanceToNextInputLine( 4 )
read(4,*) n_boff_cutoff
call AdvanceToNextInputLine( 4 )
read(4,*) movegrid,ndim
call AdvanceToNextInputLine( 4 )
read(4,*) demf, mix_strain, mix_stress

! OUTPUT PARAMETERS
call AdvanceToNextInputLine( 4 )
read(4,*) time_max
time_max = time_max * 1000 * sec_year
call AdvanceToNextInputLine( 4 )
read (4,*) dtout_screen
dtout_screen = dtout_screen * 1000 * sec_year
call AdvanceToNextInputLine( 4 )
read(4,*) dtout_file
dtout_file = dtout_file * 1000 * sec_year
call AdvanceToNextInputLine( 4 )
!read(4,*) nprofil
!do 135 i=1,nprofil
!read(4,*) hv_out(i:i), iprof_out(i) 
!135 continue
!call AdvanceToNextInputLine( 4 )
read(4,*) io_vel,io_srII,io_eII,io_aps,io_sII,io_sxx,io_szz,io_sxz,io_pres, &
    io_temp,io_melt,io_visc,io_phas,io_mark,io_src,io_diss,io_forc,io_hfl,io_topo,io_thermochron
call AdvanceToNextInputLine( 4 )
read(4,*) lastout
call AdvanceToNextInputLine( 4 )
read(4,*) dtsave_file
dtsave_file = dtsave_file * 1000 * sec_year
call AdvanceToNextInputLine( 4 )
read(4,*) lastsave

close (iu)


! ADDITIONAL PARTICULAR INPUT
!call ReadMoreParams()

return
end


!========================================================
subroutine AdvanceToNextInputLine( iu )
character*1 buf

10    read(iu, '(A1)') buf
      if( buf(1:1).eq.';' ) then
          goto 10
      else
          backspace( iu )
          return
      endif

print *, 'AdvanceToNextInputLine: EOF reached!'
stop

print *, 'AdvanceToNextInputLine: Error reading file!'
stop

return
end


subroutine ReadMoreParams()

call ReadIntrusions()  ! - see user_ab.f90

call ReadHydro()       ! - see user_luc.f90

return
end

