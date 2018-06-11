#ifndef IGNHSEEN
#define IGNHSEEN

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <string.h>
#include <ctype.h>
#include <assert.h>

/*  CVODE headers  */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense                */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM        */
#include <sundials/sundials_types.h> /* definition of type realtype          */

#include "TC_interface.h"
#include "TC_params.h"

#define MAX(A,B) ( ((A) > (B)) ? (A) : (B) )
#define MIN(A,B) ( ((A) < (B)) ? (A) : (B) )

#define MAXSPECIN  200
#define MAXSPECOUT 200

typedef struct {
  double Tign;
  int ncalls;
} cdt;

/* Time stepping/CVode info */
typedef struct {

  double t;           /* current time */
  double tsta;        /* start time */
  double tend;        /* end   time */
  double deltat;      /* current time step */
  double deltatMax;   /* maximum time step */
  double deltaTemp;   /* maximum temperature change per time step*/

  double Tini;        /* initial temperature */   
  double Temp_id ;    /* temperature threshold for ignition delay */

  double pressure ;   /* atmospheric pressure */  
  double pfac ;       /* pressure factor */

  int oFreq;          /* Output frequency */
  int NiterMax;       /* Maximum no. of iterations */

  double *scal ;      /* array, Nspec+1 long with T and Y's */
  N_Vector y0 ;       /* CVode vector, Nspec+1 long with T and Y's */
  realtype relT;      /* relative tolerance */
  N_Vector absT;      /* CVode vector, absolute tolerances */

  int getDeltatSeq;   /* use prescribed time steps ? */
  double *deltatSeq;  /* array with prescribed time steps*/

  int getIgnDel ;     /* flags for ignition delay and local SA */
  int getSens   ; 

  double CVsmall, CVrelt;        /* CVode parameters */
  int CVmaxord, CVmaxnumsteps;

} tadvInfo ;

/* Output files */
typedef struct {
  FILE *myfile, *myfile1, *myfile2, *myfile3 ;
  int  *specOutIDs;
  int  specoutno; 
  char *SpecOutName;
} ioInfo;

void initCVODE (void **cvode_mem, tadvInfo *advdata, cdt *udata) ;

double doIgn      (tadvInfo *advdata, ioInfo *iodata, void *cvode_mem) ;
double doIgnReinit(tadvInfo *advdata, ioInfo *iodata, void *cvode_mem)  ;
double doIgnDtSeq (tadvInfo *advdata, ioInfo *iodata, void *cvode_mem) ;

void Output(int iflag, int iter, tadvInfo *advdata, ioInfo *iodata) ;

void setup(tadvInfo *advdata, ioInfo *iodata, unsigned int *withTab,
           char *SpecName, double *SpecMsFr, int *specinno, 
           char *mechfile, char *thermofile) ;


void chemrhs(realtype t, N_Vector y, N_Vector ydot, void *f_data) ;
int chemrhswrapper(realtype t, N_Vector y, N_Vector ydot, void *f_data) ;

#ifdef USEJAC
void chemjac(realtype t, N_Vector y, DlsMat Jac, double *jactmp) ; 
int chemjacwrapper(long int N, realtype t, 
                   N_Vector y, N_Vector fy, DlsMat J, void *udata,
		   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) ;
#endif

int g(realtype t, N_Vector y, realtype *gout, void *user_data) ;
int Check_CVflag(void *flagvalue, char *funcname, int opt) ;

#endif
