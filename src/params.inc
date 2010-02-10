! -*- F90 -*-

common /mesh/ x0,z0,rxbo,rzbo,sizez_x(10),sizez_y(10), &
    dx_rem,angle_rem,anglemin1,anglemint,topo_kappa,bottom_kappa, &
    velbc_l,velbc_r,nx,nz,nq,nzonx,nzony,nelz_x(10),nelz_y(10), &
    ny_rem,mode_rem,ntest_rem,v_min,v_max,ivis_shape,igeotherm,efoldc, &
    g_x0,g_y0c,g_amplitude,g_width,ydrsides,ynstressbc, &
    ny_inject,nelem_inject,rate_inject,nmass_update,nopbmax,iinj, &
    nofbc,nofside(20),nbc1(20),nbc2(20),nbc(20), &
    bca(20),bcb(20),bcc(20),dx_init,xReyn,i_rey,iac_rem, &
    bcd(20),bce(20),bcf(20),bcg(20),bch(20),bci(20),iph_int, & 
    iynts,tbos, iax1,iay1,ibx1,iby1,icx1,icy1,idx1,idy1
common /flac_ctrl/ dt_scale,strain_inert,vbc,amul,ratl,ratu,frac, &
    dtmax_therm,dt_maxwell,fracm,srate0, &
    dt_elastic,dt_elastic0,demf,boff,mix_strain,mix_stress, &
    ivis_present,n_boff_cutoff,idt_scale,ifreq_imasses,ifreq_rmasses, &
	movegrid,ndim,ifreq_visc,nmtracers
    
character phasefile*20,tempfile*20
common /inout/ dtout_screen,dtout_file,dtsave_file,lastsave,lastout, &
        io_vel,io_srII,io_eII,io_aps,io_sII,io_sxx,io_szz, &
		io_sxz,io_pres,io_temp,io_melt,io_visc,io_phas,io_mark, &
		io_src,io_diss,io_forc,io_hfl,io_topo, &
		irphase,irtemp,phasefile,tempfile,irh_mark

common /rheol/ visc(20),den(20),alfa(20),beta(20),pln(20),acoef(20),eactiv(20),rl(20),rm(20), &
    coha(20),cohdisp(20),phimean(20),phidisp(20),psia(20),conduct(20),cp(20), &
	ts(20),tl(20),tk(20),fk(20),ten_off, &
	plstrain(5),fric(5),dilat(5),cohesion(5),xinitaps,tau_heal, &
        plstrain_oc(5), fric_oc(5), dilat_oc(5), cohesion_oc(5), &
        plstrain_n(5), fric_n(5), dilat_n(5), cohesion_n(5), phase_counter(20), &
	nphase,mphase,irheol(20),nysoft,nsegments,iynocean, &
	inhom,ix1(20),ix2(20),iy1(20),iy2(20),inphase(20),igeom(20),nphasl, &
    ltop(20),lbottom(20),lphase(20),max_markers,nmarkers, &
        iint_marker,iint_tracer,nzone_tracer,dt_outtracer,itx1(20),itx2(20),ity1(20),ity2(20), &
        nzone_marker,imx1(20),imx2(20),imy1(20),imy2(20)

common /therm/ t_top,t_bot,hs,hr,temp_per,bot_bc,time_t,dt_therm, &
	itherm,istress_therm,initial_geoth,itemp_bc,ix1t,ix2t,iy1t,iy2t,ishearh, &
	ntherm,hc(20),age_1(20),ixtb1(20),ixtb2(20),nzone_age,i_prestress, &
       iph_col1(20), iph_col2(20)

common /stress/ g,pisos,drosub,damp_vis,nyhydro,iphsub

common /process/ sec_year,time,dt,nloop,irestart,time_max,ifreq_avgsr,ireset

common /hooks/ if_intrus,if_hydro

common /temporary/ devmax, dvmax
