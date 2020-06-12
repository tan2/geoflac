# INPUT PARAMETER FILE FOR GEOFLAC
### (by David Okaya & Eh Tan - November 2011)

-----------------------------------------------------------------------------

### Notes:

1. This document is in Markdown format and is best viewed online at [this link](https://github.com/tan2/geoflac/doc/input_description.md).

2. Any number of comment lines can be present and are preceded by a semi-colon.

3. If the line contains more than one parameter, the parameters are separated by spaces and/or optional comma.

4. The parameter name is directly taken from the source code. Due to markdown syntax, some parameter names that contain '_' (underscore) are replaced with '-' (hyphen).

5. Parameters in italic are repeated for #lines based on the prior parameter. E.g.,

    **nitems**

    **_value1(i), value2(i),_** .... (this line repeats **nitems** times. If **nitems**=0, then NO italic line should be given.)


----
### MESH PARAMETERS

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
|**nx, nz** | 2 int | #Nodes in X and Z directions.|
|**x0, z0** | 2 dbl | Beginning coordinates (Upper Left).|
|**rxbo, rzbo** | 2 dbl | Size: width from left (positive), vertical thickness (negative), or coordinate value of Lower Right corner (depth is negative).|

* In each direction of the mesh, #nodes = #elems + 1 (one more nodes than elements or cells), i.e. **nx** = **nex**+1 and **nz** = **nez**+1.
* First three lines define element spacing of uniform grid. dx = **(rxbo-x0)/nex**, dz = **(rzbo-z0)/nez**. In general, it is good to have dx = dz (square elements).
* **nx** and **nz** cannot exceed the maximum value (**mnx** and **mnz**) set in `arrays.inc`.

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
|**ircoord, coordfile** | int,str | if **ircoord**=1, read the initial coordinates from the file.|
|**nzonx** (>=0) | int | #Zones in X-direction (0-uniform grid).|
|**_nelz-x(i), sizez-x(i)_** | int,dbl | #Elems/zone, size of zone (non-dimensional; percentage of **xbo**). Sum of all **nelz-x(i)** must == **nex**; sum of all **sizez-x(i)**= 1.0).|
|**nzony** (>=0) | int | #Zones in Z-direction (0-uniform grid).|
|**_nelz-y(i), sizez-y(i)_** | int,dbl | #Elems/zone, size of zone (non-dimensional; percentage of **rzbo**). Sum of all **nelz-y(i)** must = **nez**; sum of all **sizez-y(i)** = 1.0).|
|**iint-marker, iint-tracer** | int, int | Flags denoting use of markers, tracers: 0-No, 1-Yes. For the current version, **iint-marker** MUST be 1.|
|**nzone-marker** (>=0) | int | #Rectangular zones with 9 markers.|
|**_imx1(i),imy1(i),imx2(i),imy2(i)_** | 4 int | Ignored. (UL, LR corners of cells for zone markers in older version.)|
|**nzone-tracer** (>=0), **dt-outtracer** | int, dbl | #Rectangular zones with traces (stores x,z,P,T), frequency of storage in Kyrs.|
|**_itx1(i),ity1(i),itx2(i),ity2(i)_** | 4 int | UL, LR corners of cells for zone tracers.|



#### MECHANICAL CONDITIONS

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
|**nystressbc, nydrsides** | 2 int | Stress boundary conditions, lithostatic stress BC on the sides.       Both are 0-No, 1-Yes.|
|**nofbc** | int | # of boundary conditions (i.e., how many segments).|
|**_nofside(i),nbc1(i),nbc2(i),nbc(i),bca(i),bcb(i),bcc(i),bcd(i),bce(i),bcf(i),bcg(i),bch(i),bci(i)_**|  4 int, 9 dbl | See the following text. |

* There are four sides (**nofside**) for a rectangular domain. 1-left, 2-bottom, 3-right, 4-top.
* Each side is composed of one or more segments.  Each segment with its own BC behavior.  The starting and ending nodes of each segment are **nbc1** and **nbc2**.
* These are the following types of boundary conditions (**nbc**):
    +  1 - velz
    +  2 - shear stress (x,z) plane
    + 10 - velx 
    + 20 - normal stress
    + 30 - vely (strike slip version). (Experimental, used at your own risk.)
    + 40 - shear stress out of the plane. (Experimental, used at your own risk.)
    + 50 - velx with visc-dep profile. (Experimental, used at your own risk.)
* If the BC is a velocity in x or z, it is provided in units of meters/sec, a very small number.  For reference, 1 mm/year = 1 km/myr = 3.17E-11 m/s.
* The BC function **fbc(s) = a + b s + c s^2 + (d cos (2pi e s) + f sin (2pi g s))*exp(-(s - i)^2 h^2 )** is measured along this segment as a percentage of the length. **s** is non-dimensional (normalized between 0~1), i.e. **s** = (node-in-segment - **nbc1**) / (**nbc2** - **nbc1**).

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
|**nyhydro,pisos,iphsub,drosub,damp-vis**| int,dbl,int,2 dbl | Winkler-type boundary condition, ie. the bottom boundary is open for flow passing through in and out. Effectively, the bottom is supported by an inviscid fluid.  **nyhydro**:  0-no, 1-yes (manual), 2-yes (auto). If 2, the reference hydrostatic pressure is computed from the column of elements at the _right_ boundary.  **pisos**: reference hydorstatic pressure at compensation depth (used when **nyhydro**=1).  **iphsub**: substratum phase (see _PHASE_ Section below) of the inviscid fluid (whose density affects the restoration stress at the bottom).  **drosub**:  additional density difference of the fluid.  **damp-vis**:  (not used.)  NOTE: Use nyhdro=2, typically iphsub=8 (mantle phase).|
|**g**| dbl | Acceleration of gravity (m/sec^2).  (typically 9.81 or 10).|



#### THERMAL CONDITIONS

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
|**i-prestress**| int | Let model equilibrate given initial state, by running model for 600,000 years before any BC changes are allowed (0-no, 1-yes).  Temperature and stress will redistribute; topography, and mountain roots will isostatically adjust.|
|**itherm**| int | 1-mechanical+thermal calculation, 2-no mechenical. NOTE: usually itherm=1.|
|**istress-therm**| int | Add thermal stress (0-no, 1-yes). |
|**ishearh**| int | Add shear heating.|
|**t-top**| dbl | Surface temperature in Celsius. Boundary AND initial condition.|
|**t-bot**| dbl | Bottom temperature in Celsius. ONLY initial condition. (The bottom TBC is given by **bot-bc** below.)|
|**hs, hr**| 2 dbl | Radiogenic heating of crust: W/kg (~1.e-9). Radiogenic folding depth: km.|
|**itemp-bc, bot-bc**| int, dbl | Thermal boundary conditions at the bottom (1-T, 2-Flux). If **itemp-bc**=1, **bot-bc** is the temperature in Celsius. If itemp-bc=2, bot-bc is the heat flux in mW/m2.|
|**temp-per, ix1t, ix2t, iy1t, iy2t**| dbl, 4 int | Initial thermal perturbation, amplitude and location between rectangular Xleft, Xright, Ztop, Zbottom nodes.|
|**irtemp**| int | Initial thermal distributions from a file? 0-no, 1-yes.|
|**tempfile**| string | Filename containing temperature distribution. When **irtemp**=0, this file is read.  File format is one column of temperatures written in blocks of nz, looping over nx.  E.g., `do ix=1,nx; do iz=1,nz; read T(iz,ix); enddo; enddo;`|
|**time-scale** | dbl | not used.|


##### Initial temperature structure

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
|**iynts, tbos**| int, dbl | **iynts**: 1-initialize temperature of .
     2 - initialize temperature by nzone-age (below).
     Else - linear temperature profile.
    tbos: maximum or bottom temperature (only when iynts=1).      If (iynts = 2) Initially variable thermal age continental and
     oceanic lithosphere across the box.|

|**iax1,iay1,ibx1,iby1,icx1,icy1,idx1,idy1**| 8 int | not used.|
|**nzone-age**| int | # zones of different age (max 20). Used when iynts=2.|
|**_age-1(i),hc1(i),hc2(i),hc3(i),hc4(i),iph-col1(i),iph-col2(i),iph-col3(i),iph-col4(i),iph-col5(i),ixtb1(i),ixtb2(i)_**| 2 dbl, 4 int | **age-1**: thermal age (Myrs). **hc1, hc2, hc3, hc4**: Depths of layer interfaces (km). **hc3** is also the Moho depth. **iph-col1, iph-col2, iph-col3, iph-col4, iph-col5**: layer phase. **ixtb1, ixtb2**: left and right side of zone in nodes.|



#### PHASES & RHEOLOGY

There are 16 phases (rock types) that are defined for this code.  Users can add 5 more.  Users should not override these original 16 phases, because they are involved in phase changes (see newphase2marker.f90 for detail).

The predefined phases are:
 1.  basalt (anhydrate)
 2.  continental crust,  same as #6.
 3.  basalt (oceanic crust), same as #7.
 4.  olivine (mantle),  same as #8.
 5.  metasediment (schist), transformed from #10 and #11.
 6.  continental crust, same as #2.
 7.  basalt (oceanic crust), same as #3.
 8.  olivine (mantle), same as #4.
 9.  serpentinite (weak mantle), transformed from mantle (#4, #8) if residing above subducted oceanic crust (#3, #7 #10).  Will transform to hydrated mantle (#16) when T-P conditions are high enough.
 10. sediment1 
 11. sediment2, generated by erosion.
 12. weak continental crust, transformed from continental crust (#2, #6) if residing above oceanic crust or arc (#3, #7, #10, #14).
 13. eclogite, transformed from basal (#1, #3, #7) when T-P conditions are high enough.
 14. arc crust
 15. weak middle crust, transformed from continental crust (#3, #7) if stressed and heated. Disabled.
 16. hydrated mantle, will transform (partially melt) to mantle (#4) if warmer than olivine wet solidus and generate arc (#14) at surface.

Phase changes are activated and will take place among certain of the defined phases. 
1. When oceanic crust (#3, #7, #10, #14) subducts under continent (#2, #6), this original continent crust  becomes "weak". cont. crust #2, #6 ----> #12 weak crust 
2. when oceanic crust (#3, #7, #10) is subducted under overlying mantle (#4, #8) (e.g., mantle wedge),  this overlying mantle will serpentinize.  This serpentinite will transform back to #16 when brought to greater depth, and further transform back to #4 if partially melted. ocean mantle #4, #8 ---->  #9 serpentini ----> #16 hydrated mantle -> #4 mantle
3. basalt (oceanic crust) will transform into eclogite. basalt (oceanic crust) #3, #7 ----> #13 eclogite
4. erosion produces sediments at top layer. any phase at topmost layer ----> #11 sediment2


| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
|**nphase**| int | Number of different types of phases (materials), max. 20.|
|**_irheol(i),visc(i), den(i), alfa(i),beta(i), pln(i),acoef(i),eactiv(i), rl(i),rm(i), plstrain1(i),plstrain2(i), fric1(i),fric2(i), cohesion1(i),cohesion2(i), dilat1(i),dilat2(i), conduct(i),cp(i), ts(i),tl(i),tk(i),fk(i)_**| int, 24 dbl | See the following list.|


The 
+  **irheol**: Rheology for this phase (see below).
    + 1- elastic,
    + 3- visco-elastic (Maxwell, Non-Newtonian) 
    + 6- elasto-plastic (Mohr-Coulomb) 
    + 11- visco-elasto-plastic (Mohr-Coulomb + Fixed Maxwell, Newtonian) 
    + 12- visco-elasto-plastic (Mohr-Coulomb + Maxwell, Non-Newtonian)

+  **visc**:  not used.
+  **den**:  density at 0 Pa and 0 Celsius.
+  **alfa**:  alpha - coeff. of thermal expansion.
+  **beta**: beta - compressibility.
+  **pln**:  n for viscosity power law
+  **acoef**: A for viscosity power law
+  **eactiv**: E - activation energy for viscosity
+  **rl**:  Lame parameter: rl - lambda
+  **rm**:  Lame parameter: rm - mu
+  **plstrain1**: Mohr-Coulomb: plastic strain of onset softening 
+  **plstrain2**: Mohr-Coulomb: plastic strain of saturated softening 
+  **fric1**: Mohr-Coulomb: friction angle before softening
+  **fric2**: Mohr-Coulomb: friction angle of saturated softening, in degree
+  **cohesion1**: Mohr-Coulomb: cohesion before softening
+  **cohesion2**: Mohr-Coulomb: cohesion of saturated softening, in Pa
+  **fric1**: Mohr-Coulomb: dilation angle before softening
+  **fric2**: Mohr-Coulomb: dilation angle of saturated softening, in degree
+  **conduct**: thermal conductivity
+  **cp**:  heat capacity at constant pressure, J/(kg C)
+  **ts**:  solidus temperature 
+  **tl**:  liquidus temperature 
+  **tk**:  latent heat of fusion temperature
+  **fk**:  coeff of latent heat of fusion?



#### INITIAL PHASE DISTRIBUTION

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|

; Define model rock type distribution via external file or by horizontal layers 
**irphase**    int Flag to take initial phase distribution from a file: 0-N, 1-Y.
**p8hasefile**   string Filename containing phase distribution. 
     If **irphase**=0, this line is read but not used.
     One column of phase id by cell (not by node). 
     Vertical fast, horz slow dimensions.

**mphase**    int Main phase. The most abundant phase.
**nphasl**    int Number of horizontal layers for phases.
     Note: when **nzone-age**>0, the phase structure is determined by **iph-col1, iph-col2** etc. there.
 **_ltop(i), lbottom(i), lphase(i)_**
    3 int **ltop**  vertical cells of top of this layer.
     **lbottom** vertical cells of bototm of this layer.
     **lphase** phase ID of this layer (between 1 ~ nphase)

; Define rock type anomolies.
  Geometries of inhomogeneities: 
   0 -    rectangular 
   1,2 - Gauss shape (not tested, don’t use)
   3 - diagonal line
   4 - diagonal line plus init. plastic strain
  if rectangular, upper left corner is ix1-iy1, lower right corner is ix2-iy2 (in   cells).
  if diagonal line, ix1 < ix2, but the vertical iy1-iy2 can be any values, in order   to define
  diagonal slope polarity (e.g., iy1>1y2 for dipping towards left).
**inhom**    int Number of rock type inclusions (up to 9).
 **_ix1(i), ix2(i), iy1(i), iy2(i), inphase(i), igeom(i), xinitaps(i)_**
   6 int, dbl **ix1**  Inclusion left horizontal cell.
     **ix2**  Inclusion right horizontal cell.
     **iy1**  Inclusion top vertical cell.
     **iy2**  Inclusion bottom vertical cell. 
     **inphase** phase ID of this inclusion (1~nphase). If negative, the original phase is used.
     **igeom** geometry of inclusion (see code above). 
     **xinitaps**  initial plastic strain. E.g., 15 = very weak. 

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
**ten-off**    dbl Tension cut-off.   E.g., 1.E9
**tau-heal**   dbl Linear healing parameter. Plastic strain decreases with      time;  time value (in sec). 1e15 = 30 myr.  For healing, use      1e13 ~ 0.3 myr.
**v-min, v-max, ivis-shape, efoldc**  Viscosity limits
   2 dbl, int, dbl v-min  Minimum viscosity.   E.g., 1E19
      **v-max** Maximum viscosity.   E.g., 3E27
       **ivis-shape** Initial viscosity shape (0-2; see below)
      = 0  constant: visn = vis0
      = 1  linear:    visn = vis0+(yc-geoth)*efoldc
      = 2  gaussian:  visn = vis0*exp((yc-geoth)/efoldc)2
      efoldc folding depth of viscosity in meter.
**igeotherm, g-x0, g-y0c, g-amplitude, g-width**
    int, 4 dbl igeotherm 0-constant, 1-gaussian, 2-linear
     **g-x0**  see below.
     **g-y0c**  see below.
     **g-amplitude** see below.
     **g-width**   see below.
; Used only in viscosity reset calculations when iynts=1 or irheol=11 (never used because usually irheol=12).  igeotherm used for ridge thermal condition (iynts=1) and crust viscosity condition.
 ;  =0 geoth = g-y0c
 ;  =1 geoth = g-y0c + g-amplitude*exp(-((xc-g-x0)/g-width)**2.)
 ;  =2 if (abs(g-x0-xc).lt.g-width) geoth = g-y0c+ g-amplitude*(1.-abs(g-x0-xc)/g-width)
  if (abs(g-x0-xc).ge.g-width) geoth = g-y0c

**ny-inject, nelem-inject, rate-inject** Magma injection / Accretionary.  This 
        line for mid-ocean ridge problem. There
        are thermal parameters that are 
        hardcoded for this (tcutoff-eclo, etc.).
   2 int, dbl **ny-inject**  0-no, 1-left, 2-center
     **nelem-inject** vertical cell to inject.
     **rate-inject** 


#### REMESHING

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
**ny-rem, mode-rem, ntest-rem, angle-rem**
   3 int, dbl **ny-rem**  Activate remeshing? 0-no, 1-yes.
     **mode-rem**  Mode
         =1 keep 4 sides as is. 
         =3 L,R,B restored to initial sides
          top as is.
         =11 L,R restored to initial sides,
          top,bottom as is.
     **ntest-rem**  #time steps to test if remeshing needed.
     **angle-rem**  critical angle to trigger remeshing.

**dx-rem**    dbl remeshing criteria for mode-rem=3 or 11
     Normalized size = **rxb0/nex**.  dx-rem = 0~1

**topo-kappa, fac-kappa** 2 dbl Erosion diffusion of topography: diffusivity and additional factor for top boundary. Use 0 for no erosion, 1E-5 for strong erosion.



#### PROCESS CONTROL

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
; PROCESS CONTROL
! inertial mass scaling
**idt-scale**   int Inertial mass scaling.  0-no scale for DYNAMICS. 
         1-scale for STATICS. 
         2-Automatic scaling. USE THIS.
**dt-scale, strain-inert**  2 dbl control of time step size; dt-scale is used if idt-scale       =1.
     strain-inert is used if idt-scale=2. Decrease it if the model
     gives NaN or Inf result.
**i-rey, xReyn**   int, dbl Buoyancy + Reynolds (no-0,yes-1), Reynolds number

**ifreq-rmasses** int Frequency of re-evaluation of real masses
**ifreq-imasses** int Frequency of re-evaluation of inertial masses
**ifreq-visc**  int Freq. re-evaluation Non-Newtonian VISC (rheol.eq.12)
**ifreq-avgsr** int navgsr - freq. averaging strain rate and dissipation

; Acceleration Parameters (Cundall, 1982)
**amul, ratl, ratu** 3 dbl Internal adaptive time scaling.
**frac, fracm**  2 dbl calculate time step scales
**n-boff-cutoff**  int force balance cut-off.
**movegrid, ndim**  2 int movegrid  (0-N, 1-Y, 2-move under strain-rate. Use =1.
	**ndim** Dim: =2 normal cases; =3 in/out of plane 							stress.
**demf, mix-strain, mix-stress** dbl 2 int **demf** damping coefficient
						Mixing Procedures:	mix-strain 	0-no, 1-yes.
						**mix-stress**	0-no, 1-yes.

#### OUTPUT PARAMETERS

| Parameters  | Types |  Description  |
|:------------|:-----:|:--------------|
; Time parameters below are in thousands years
**time-max**   dbl Max time of calculations (in Kyrs).
**dtout-screen**  dbl Time interval for screen output of calc. progress.
**dtout-file**   dbl Time interval for file output (frequency of output)

; Variables to print to file:   0-no, 1-yes
**io-vel,io-srII,io-eII,io-aps,io-sII,io-sxx,io-szz,io-sxz,io-pres,io-temp,io-phase,io-visc,io-unused,io-density,io-src,io-diss,io-forc,io-hfl,io-topo**
   19 int 
 io-vel  Velocity (Vx,Vz) in cm/year.       vx.0, vz.0
 io-srII  Deviatoric strain rate - second invariant, log10().  srII.0
 io-eII  Dev. strain - second invariant, first invariant.        eII.0, eI.0
 io-aps Total accumulated plastic strain    aps.0
 io-sII  Dev. stress - sqrt(second invariant), in kbars  sII.0
 io-sxx Dev. stress – xx component, in kbars   sxx.0
 io-szz  Dev. stress – zz component, in kbars   szz.0
 io-sxz Stress - xz (shear) component,  in kbars   sxz.0
 io-pres Pressure, in kbars       pres.0
 io-temp Temperatur, in Celsius      temp.0
 io-phase Phase # (1-n as defined in input file).   phase.0
 io-visc Viscosity, log10(Pa-s)      visc.0
 io-unused unused
 io-density Density, in kg/m       density.0
 io-src  Heat source.       src.0
 io-diss Energy dissipation = shearheat/dens/radiog-heat diss.0
 io-forc Force: at left, right boundaries    forc.0
 io-hfl  Surface heat flow (1-D)      hfl.0
 io-topo Topography (1-D), in km     topo.0



 ; By default, these always produced:
   Time of each output panel in sec    time.0
   Mesh coordinates (z,x)x, (z,x)z in km    mesh.0

**lastout**    int Output for last step only (1) or each nout step (0)
      0=append, 1=override in output file.
**dtsave-file**   dbl Time interval for process saving (in Kyrs). For restart *.rs files.
**lastsave**    int Saving the last step only (1) or each nsave step (0)

     To restart,  (a) copy -contents.save -> _contents.rs
       (b) rerun flac with the same input parameters



#### END OF FILE

