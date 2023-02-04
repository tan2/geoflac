! Reads problem parameters from input file

subroutine read_params(inputfile)
use arrays
use params
include 'precision.inc'

character*200 inputfile

iu =4 
open( iu, file=inputfile )

! MESH
call AdvanceToNextInputLine( 4 )
read(4,*) nx, nz
open(11,file='nxnz.0')
write(11,*) nx, nz
close(11)

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
    sizez_x(1) = 1.d0
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
    sizez_y(1)   = 1.d0
    go to 177
endif
do i = 1,nzony
    sizez_y(i) = 1.d0
    call AdvanceToNextInputLine( 4 )
    read(4,*) nelz_y(i), sizez_y(i)
end do
177 continue

! MECHANICAL CONDITIONS
call AdvanceToNextInputLine( 4 )
read(4,*) nystressbc,nydrsides
call AdvanceToNextInputLine( 4 )
read(4,*) nofbc
call AdvanceToNextInputLine( 4 )
        do i = 1,nofbc
        read(iu,*) nofside(i),nbc1(i),nbc2(i),nbc(i),  &
        bca(i),bcb(i),bcc(i),  &
        bcd(i),bce(i),bcf(i),bcg(i),bch(i),bci(i)
        enddo
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
! temp structure
call AdvanceToNextInputLine( 4 )
read (4,*) nzone_age
call AdvanceToNextInputLine( 4 )
do i = 1, nzone_age
    read (4,*) ictherm(i),age_1(i),tp1(i),tp2(i),ixtb1(i),ixtb2(i)
    call AdvanceToNextInputLine( 4 )
    read(4,*) nph_layer(i), (hc(i,j), j=1,nph_layer(i)-1)
    !print *, (hc(i,j), j=1,nph_layer(i)-1)
    read(4,*) (iph_col(i,j), j=1,nph_layer(i))
    !print *, (iph_col(i,j), j=1,nph_layer(i))
enddo
! check smooth nzone_age
do i = 1, nzone_age-1
    iph_col_trans(i) = 0
    if (ixtb2(i) == -1) then
        if (ixtb1(i+1) /= -1) then
            print *, 'Error: ixtb1 is not -1 at', i+1, 'column!'
        endif
        if (ictherm(i) /= ictherm(i+1)) then
            print *, 'Error: ictherm at', i, i+1, 'columns are not equal!'
        endif
        if (nph_layer(i) /= nph_layer(i+1)) then
            print *, 'Error: nph_layer at', i, i+1, 'columns are not equal!'
        endif
        do j = 1, nph_layer(i)
            ! the phases in the columns must match
            if (iph_col(i,j) /= iph_col(i+1,j)) then
                print *, 'Error: iph_col at', i, i+1, 'columns are not equal!'
            endif
        enddo
        ! flag indicates smooth transition
        iph_col_trans(i) = 1
        ! remove '-1'
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
end do
print *, nphase, 'phases'
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

! inclusions
call AdvanceToNextInputLine( 4 )
read(4,*) inhom
if( inhom .gt. maxinh ) then
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
read(4,*) itype_melting, nelem_serp, prod_magma, rho_magma
call AdvanceToNextInputLine( 4 )
read(4,*) angle_mzone, fmagma_max, ratio_mantle_mzone
call AdvanceToNextInputLine( 4 )
read(4,*) latent_heat_magma, lambda_freeze, lambda_freeze_tdep
call AdvanceToNextInputLine( 4 )
read(4,*) weaken_ratio_plastic, weaken_ratio_viscous
weaken_ratio_viscous = log(weaken_ratio_viscous)

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
! diffusion of topography
call AdvanceToNextInputLine( 4 )
read(4,*)  topo_kappa, fac_kappa


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

if (mod(ntest_rem, ifreq_rmasses) + mod(ntest_rem, ifreq_imasses) &
    + mod(ntest_rem, ifreq_visc) + mod(ntest_rem, ifreq_avgsr) .ne. 0) then
        call SysMsg('ntest_rem must be multiples of process frequency')
        stop
endif

call AdvanceToNextInputLine( 4 )
read(4,*) frac, fracm
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
read(4,*) io_vel,io_srII,io_eII,io_aps,io_sII,io_sxx,io_szz,io_sxz,io_pres, &
    io_temp,io_phase,io_visc,io_unused,io_density,io_src,io_diss,io_forc,io_hfl,io_topo
call AdvanceToNextInputLine( 4 )
read(4,*) lastout
call AdvanceToNextInputLine( 4 )
read(4,*) dtsave_file
dtsave_file = dtsave_file * 1000 * sec_year
call AdvanceToNextInputLine( 4 )
read(4,*) lastsave

close (iu)

return
end


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
