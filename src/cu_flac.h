// -*- C++ -*-

#ifndef __cu_flac_h__
#define __cu_flac_h__


#define mnx     50      // max # of elements in x
#define mnz     10      // max # of elements in z
#define ntriag  4       // # of triangles per element
#define nstr    4       // # of stress components
#define ndim    2       // # of spatial dimension


// Macros to wrap index access to Fortran arrays
// See 'arrays.F90' for the dimension extent of each array

#define amass(j,i) (                    \
        amass_d[ ((i)-1)*nz + (j)-1 ]   \
	)

#define area(j,i,d) (                                                   \
        area_d[ ((d)-1)*(nx-1)*(nz-1) + ((i)-1)*(nz-1) + (j)-1 ]        \
	)

#define balance(j,i,d) (                                \
        balance_d[ ((d)-1)*nx*nz + ((i)-1)*nz + (j)-1 ] \
	)

#define bc(j,i,d) (                                     \
        bc_d[ ((d)-1)*nx*nz + ((i)-1)*nz + (j)-1 ]      \
	)

#define cord(j,i,d) (                                   \
        cord_d[ ((d)-1)*nx*nz + ((i)-1)*nz + (j)-1 ]    \
	)

#define dvol(j,i,d) (                                                   \
        dvol_d[ ((d)-1)*(nx-1)*(nz-1) + ((i)-1)*(nz-1) + (j)-1 ]        \
	)

#define force(j,i,d) (                                  \
        force_d[ ((d)-1)*nx*nz + ((i)-1)*nz + (j)-1 ]   \
	)

#define ncod(j,i,d) (                                   \
        ncod_d[ ((d)-1)*nx*nz + ((i)-1)*nz + (j)-1 ]    \
	)

#define rmass(j,i) (                    \
        rmass_d[ ((i)-1)*nz + (j)-1 ]   \
	)

#define strain(j,i,d) (                                                 \
        strain_d[ ((d)-1)*(nx-1)*(nz-1) + ((i)-1)*(nz-1) + (j)-1 ]      \
	)

#define stress0(j,i,d,k) (                                              \
        stress0_d[ ((k)-1)*nstr*(nx-1)*(nz-1) +                         \
                   ((d)-1)*(nx-1)*(nz-1) + ((i)-1)*(nz-1) + (j)-1 ]     \
	)

#define temp(j,i) (                     \
        temp_d[ ((i)-1)*nz + (j)-1 ]    \
	)

#define vel(j,i,d) (                                    \
        vel_d[ ((d)-1)*nx*nz + ((i)-1)*nz + (j)-1 ]     \
	)


#define irheol(i) (PAR.irheol[(i)-1])
#define visc(i) (PAR.visc[(i)-1])
#define den(i) (PAR.den[(i)-1])
#define alfa(i) (PAR.alfa[(i)-1])
#define beta(i) (PAR.beta[(i)-1])
#define pln(i) (PAR.pln[(i)-1])
#define acoef(i) (PAR.acoef[(i)-1])
#define eactiv(i) (PAR.eactiv[(i)-1])
#define rl(i) (PAR.rl[(i)-1])
#define rm(i) (PAR.rm[(i)-1])
#define coha(i) (PAR.coha[(i)-1])
#define cohdisp(i) (PAR.cohdisp[(i)-1])
#define phimean(i) (PAR.phimean[(i)-1])
#define phidisp(i) (PAR.phidisp[(i)-1])
#define psia(i) (PAR.psia[(i)-1])
#define conduct(i) (PAR.conduct[(i)-1])
#define cp(i) (PAR.cp[(i)-1])
#define ts(i) (PAR.ts[(i)-1])
#define tl(i) (PAR.tl[(i)-1])
#define tk(i) (PAR.tk[(i)-1])
#define fk(i) (PAR.fk[(i)-1])
#define lphase(i) (PAR.lphase[(i)-1])


// struct for various input parameters
typedef struct _flac_param {
    double visc[20];
    double den[20];
    double alfa[20];
    double beta[20];
    double pln[20];
    double acoef[20];
    double eactiv[20];
    double rl[20];
    double rm[20];
    double coha[20];
    double cohdisp[20];
    double phimean[20];
    double phidisp[20];
    double psia[20];
    double conduct[20];
    double cp[20];
    double ts[20];
    double tl[20];
    double tk[20];
    double fk[20];

    double g;
    double pisos;
    double drosub;
    double rzbo;
    double demf;
    double sec_year;
    double ynstressbc;
    double dt_scale;
    double frac;
    double fracm;
    double strain_inert;
    double vbc;

    int irheol[20];
    int lphase[20];

    int nyhydro;
    int iphsub;
    int n_boff_cutoff;
    int i_prestress;
    int iint_marker;
    int nphasl;
    int idt_scale;
} flac_param_t;


__constant__ flac_param_t PAR;
// __constant__ memory has implicit static storage (ie. local to each file)!!!
// Functions using PAR must reside in the same file.

#endif

