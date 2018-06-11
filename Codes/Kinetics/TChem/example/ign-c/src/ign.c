#include "ign.h"

void getmeps(double *rtol, double *atol)
{

  double sml = 1.0 ;
  while ( 1.0+sml != 1.0 ) sml *= 0.5 ;
  *rtol = sqrt(2.0*sml) ;
  *atol = *rtol ;
  return ;

}
int main()
{

  tadvInfo advdata;
  ioInfo   iodata;

  /* Integer parameters */
  int Nspec, Nvars, Nreac ;

  /* Character array parameters */
  int sizecMax = 100, specinno ;
  char *mechfile, *thermofile ;

  int    lsnm = LENGTHOFSPECNAME;
  char   *SpecName ;
  double *SpecMsFr ;

  unsigned int withTab ;  /* tabulation flag */

  /* work parameters */
  int i ;
  double time_id1,time_id2 ;

  /* CVODE-related parameters */
  int iter, ierr ;
  cdt udata ;

  mechfile   = (char *) malloc( sizecMax * sizeof(char) ) ;
  thermofile = (char *) malloc( sizecMax * sizeof(char) );
  memset(mechfile  , 0, sizecMax) ;
  memset(thermofile, 0, sizecMax) ;

  /* Arrays for input/output species info */
  SpecName = (char   *) malloc( MAXSPECIN * lsnm *sizeof(char  ) ) ;
  SpecMsFr = (double *) malloc( MAXSPECIN *       sizeof(double) ) ;
  memset(SpecName, 0, MAXSPECIN * lsnm) ;
  for ( i = 0 ; i < MAXSPECIN; i++ ) SpecMsFr[i] = 0.0 ;

  iodata.SpecOutName = (char *) malloc( MAXSPECOUT* lsnm *sizeof(char)) ;
  memset(iodata.SpecOutName, 0, MAXSPECOUT * lsnm) ;

  setup(&advdata,&iodata,&withTab,SpecName, SpecMsFr, &specinno,
        mechfile, thermofile) ;

  /*
                           _
                  ___  ___| |_ _   _ _ __
                 / __|/ _ \ __| | | | '_ \
                 \__ \  __/ |_| |_| | |_) |
                 |___/\___|\__|\__,_| .__/
		 |_|

  */
  /* Initialize TC library */
  TC_initChem( mechfile, thermofile, (int) withTab, 1.0) ;
  free(mechfile  ) ;
  free(thermofile) ;

  /* Set pressure */
  advdata.pressure *= advdata.pfac ;
  TC_setThermoPres(advdata.pressure) ;

  Nspec = TC_getNspec() ;
  Nreac = TC_getNreac() ;
  Nvars  = Nspec+1 ;       /* T+(Nspec) */
  printf("Nspec/Nreac = %d/%d\n",Nspec,Nreac) ;
  advdata.scal = (double *) malloc( (Nspec+1)*sizeof(double) ) ;

  /* Set initial conditions */
  printf("Set initial conditions \n") ;
  advdata.scal[0] = advdata.Tini ;
  for ( i = 1 ; i<Nspec+1 ; i++) advdata.scal[i] = 0.0;

  for ( i = 0 ; i < specinno ; i++)
  {
    int ispec = TC_getSpos( &SpecName[i*lsnm], strlen(&SpecName[i*lsnm]) ) ;
    if (ispec < 0 )
    {
      printf("Error : Could not find species %s -> Exit !\n",&SpecName[i*lsnm]);
      exit(1) ;
    }
    else
      printf("Index of species %s is %d\n",&SpecName[i*lsnm],ispec) ;
    advdata.scal[ispec+1] = SpecMsFr[i] ;
  }
  free(SpecName  ) ;
  free(SpecMsFr  ) ;

  /* Set ID's for custom species output */
  iodata.specOutIDs = NULL;
  if (iodata.specoutno>0) {
    iodata.specOutIDs = (int *) malloc( Nspec*sizeof(int) ) ;
    for ( i = 0 ; i<Nspec ; i++) iodata.specOutIDs[i] = 0;
    for ( i = 0 ; i < iodata.specoutno ; i++) {
      int ispec = TC_getSpos( &iodata.SpecOutName[i*lsnm], strlen(&iodata.SpecOutName[i*lsnm]) ) ;
      if (ispec < 0 ) {
        printf("Error : Could not find species %s -> Exit !\n",&iodata.SpecOutName[i*lsnm]);
        exit(1) ;
      }
      iodata.specOutIDs[ispec] = 1 ;
    }
  }
  free(iodata.SpecOutName);

  /* convert mole to mass fractions */
  double *msfr = (double *) malloc(Nspec * sizeof(double)) ;
  TC_getMl2Ms( &(advdata.scal[1]), Nspec, msfr ) ;
  for ( i = 1 ; i < Nspec+1 ; i++) advdata.scal[i] = msfr[i-1];
  iodata.myfile  = fopen("msfrini.dat","w" ) ;
  for ( i = 0 ; i < Nspec ; i++) fprintf(iodata.myfile,"%20.12e\n",msfr[i]);
  fclose(iodata.myfile);
  free(msfr) ;

  for ( i = 1 ; i < Nspec+1 ; i++)
    if ( fabs(advdata.scal[i]) > 1.e-15 )
      printf("Mass fraction of species %d is %20.12e \n",i-1,advdata.scal[i]);

  /* double *cps0 = (double *) malloc(Nspec * sizeof(double)) ; */
  /* double *cps1 = (double *) malloc(Nspec * sizeof(double)) ; */
  /* double *cps2 = (double *) malloc(Nspec * sizeof(double)) ; */
  /* double *cps3 = (double *) malloc(Nspec * sizeof(double)) ; */
  /* double *cps4 = (double *) malloc(Nspec * sizeof(double)) ; */
  /* TC_getCpSpecMl( 100.0, Nspec, cps0 ) ; */
  /* TC_getCpSpecMl( 100.0, Nspec, cps1 ) ; */
  /* TC_getCpSpecMl( 1200.0, Nspec, cps2 ) ; */
  /* TC_getCpSpecMl( 120000.0, Nspec, cps3 ) ; */
  /* TC_getCpSpecMl( 100000.0, Nspec, cps4 ) ; */
  /* for ( i = 1 ; i < Nspec ; i++) { */
  /*   printf("%d: %20.12e  %20.12e  %20.12e  %20.12e  %20.12e\n",i+1,cps1[i]-cps0[i],cps1[i],cps2[i],cps3[i],cps4[i]-cps3[i]); */
  /* } */
  /* exit(0); */

  /*
             CVODE setup
  */
  advdata.tsta=0.0;   /* start time */

  /* initial condition array */
  advdata.y0 = N_VNew_Serial( Nvars ) ;

  /* cvode tolerances */
  advdata.relT = advdata.CVrelt ;
  advdata.absT = N_VNew_Serial( Nvars ) ;

  /* set initial conditions and absolute tolerances */
  for( i = 0; i < Nvars ; i++)
  {
    double yscal = MAX(advdata.scal[i],0.0) ;
    NV_Ith_S(advdata.y0  ,i) = yscal ;
    NV_Ith_S(advdata.absT,i) = MAX( advdata.relT*yscal, advdata.CVsmall ) ;
  }

  /* declare work space (mostly for jacobian) */
  //udata = (double *) malloc ( (Nvars*Nvars+1)* sizeof(double))  ;
  //udata[Nvars*Nvars] = Temp_id ;
  udata.Tign   = advdata.Temp_id;
  udata.ncalls = 0;

  /* Create cvode solver */
  void *cvode_mem = NULL;
  initCVODE (&cvode_mem, &advdata, &udata) ;

  printf("Starting time advancement : \n") ;
  //t = advdata.tsta ;
  iter = 0      ;

#ifndef NO_OUTPUT
  iodata.myfile  = fopen("ignsol.hdr","w" ) ;
  iodata.myfile1 = fopen("ys.hdr"    ,"w" ) ;
  iodata.myfile2 = fopen("cs.hdr"    ,"w" ) ;
  iodata.myfile3 = fopen("h.hdr"     ,"w" ) ;
  Output(0, iter, &advdata, &iodata ) ;
  fclose(iodata.myfile);
  fclose(iodata.myfile1);
  fclose(iodata.myfile2);
  fclose(iodata.myfile3);

  iodata.myfile  = fopen("ignsol.dat","w" ) ;
  iodata.myfile1 = fopen("ys.out"    ,"w" ) ;
  iodata.myfile2 = fopen("cs.out"    ,"w" ) ;
  iodata.myfile3 = fopen("h.out"     ,"w" ) ;
  Output(1, iter, &advdata, &iodata ) ;
#endif

  iter = advdata.NiterMax ;
  time_id1 = 0.0  ;
  time_id2 = 0.0  ;

  if (  advdata.getDeltatSeq > 0 )
    time_id1 =doIgnDtSeq( &advdata, &iodata, cvode_mem) ;
  else
    if (  advdata.getSens > 0 ) {

      double rtol, atol;
      int iscal, ir ;
      getmeps(&rtol,&atol) ;

      double *ysave    = (double *) malloc(Nvars * sizeof(double));
      double *absTsave = (double *) malloc(Nvars * sizeof(double));
      double arrh, arrhmod, perturb ;
      for ( iscal = 0; iscal < Nvars; iscal++ ) {
	ysave[iscal]    = NV_Ith_S(advdata.y0,iscal  ) ;
	absTsave[iscal] = NV_Ith_S(advdata.absT,iscal) ;
      }

      for ( ir = 0; ir < Nreac; ir++) {
        /* restore IC */
	for ( iscal = 0; iscal < Nvars; iscal++ ) {
	  NV_Ith_S(advdata.y0,iscal)    = ysave[iscal] ;
	  NV_Ith_S(advdata.absT,iscal)  = absTsave[iscal] ;
	}
        /* get Arrhenius factor and perturb "-" */
	ierr = TC_getArhenFor(ir,advdata.getSens-1,&arrh) ;
	if (ierr != 0) {
          printf("ERROR in %s at line %d\n",__FILE__,__LINE__);
          exit(1);
	}
	perturb = fabs(rtol*arrh) ;
	if ( perturb == 0.0 ) perturb = atol ;
	arrhmod = arrh-perturb;
	ierr = TC_chgArhenFor(ir, advdata.getSens-1, arrhmod) ;
	if (ierr != 0) {
          printf("ERROR in %s at line %d\n",__FILE__,__LINE__);
          exit(1);
	}
	iter = advdata.NiterMax ;
	advdata.t = advdata.tend ;
	time_id1 = doIgn( &advdata, &iodata, cvode_mem) ;

        /* restore IC */
	for ( iscal = 0; iscal < Nvars; iscal++ ) {
	  NV_Ith_S(advdata.y0,iscal)    = ysave[iscal] ;
	  NV_Ith_S(advdata.absT,iscal)  = absTsave[iscal] ;
	}

        /* get perturb "+" */
 	arrhmod = arrh+perturb;
	ierr = TC_chgArhenFor(ir, advdata.getSens-1, arrhmod) ;
	if (ierr != 0) {
          printf("ERROR in %s at line %d\n",__FILE__,__LINE__);
          exit(1);
	}
	iter = advdata.NiterMax ;
	advdata.t = advdata.tend ;
	time_id2 = doIgn( &advdata, &iodata, cvode_mem) ;
	double sval = (time_id2-time_id1)/(2.0*perturb);
	printf("Reaction %d: %20.12e %20.12e, %20.12e %20.12e\n", ir, time_id1,time_id2, sval, sval*arrh);
	ierr = TC_chgArhenFor(ir, advdata.getSens-1, arrh) ;
	if (ierr != 0) {
          printf("ERROR in %s at line %d\n",__FILE__,__LINE__);
          exit(1);
	}
      } /* done loop over all reactions */
    } /* done if sensitivity analysis */
    else
      time_id1 = doIgn( &advdata, &iodata, cvode_mem) ;

#ifndef NO_OUTPUT
  Output(1, iter, &advdata, &iodata) ;
  fclose(iodata.myfile) ;
  fclose(iodata.myfile1) ;
  fclose(iodata.myfile2) ;
  fclose(iodata.myfile3) ;
#endif

  printf("Rootfinding function called  : %d times\n",udata.ncalls) ;

  FILE *myfileID=fopen( "tid.dat", "w" );
  fprintf(myfileID,"%20.12e\n",time_id1) ;
  printf("%20.12e\n",time_id1) ;
  fclose(myfileID) ;

  /* Clean-up */
  N_VDestroy_Serial( advdata.y0   ) ;
  N_VDestroy_Serial( advdata.absT ) ;
  CVodeFree( &cvode_mem ) ;

  TC_reset() ;

  if ( advdata.scal      != NULL ) free(advdata.scal     ) ;
  if ( advdata.deltatSeq != NULL ) free(advdata.deltatSeq) ;

  return ( 0 ) ;

}
