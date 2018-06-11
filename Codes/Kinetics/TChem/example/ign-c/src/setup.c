#include "ign.h"

#define FNAMEMAX 100

void setup(tadvInfo *ad, ioInfo *iod, unsigned int *withTab,
           char *SpecName, double *SpecMsFr, int *specinno,
           char *mechfile, char *thermofile)
{
  int i,lsnm ;
  char filDeltatSeq[FNAMEMAX] ;

  lsnm = LENGTHOFSPECNAME;
  /*
      _       __             _ _
   __| | ___ / _| __ _ _   _| | |_ ___
  / _` |/ _ \ |_ / _` | | | | | __/ __|
 | (_| |  __/  _| (_| | |_| | | |_\__ \
  \__,_|\___|_|  \__,_|\__,_|_|\__|___/

  */
  ad->NiterMax  = 100000    ; /* Maximum no. of iterations */
  ad->oFreq     = 10        ; /* Output frequency          */
  ad->Tini      = 1000.0    ; /* Initial temperature  [K]  */
  ad->Temp_id   = 1500.0    ; /* Threshold temperature for ignition delay [K]  */
  ad->deltat    = 1.e-10    ; /* Initial time step size [s]  */
  ad->deltatMax = 1.e-4     ; /* Maximum time step size [s]  */
  ad->tend      = 2.0e0     ; /* End integration time   [s]  */
  ad->deltaTemp = 1.0       ; /* Maximum temperature change per time step [K]  */
  ad->pressure  = 1.01325e5 ; /* Pressure [Pa] */
  ad->pfac      = 1.0       ; /* Pressure scale factor [ ]  */

  ad->getDeltatSeq = 0      ; /* Do not use a sequence of dt  */
  ad->deltatSeq    = NULL   ;

  ad->CVrelt        = 1.e-8  ;
  ad->CVsmall       = 1.e-17 ;
  ad->CVmaxord      = 5      ;
  ad->CVmaxnumsteps = 10000  ;

  strcpy(mechfile  ,"chem.inp" ) ;
  strcpy(thermofile,"therm.dat") ;

  *withTab = 0 ;  /* no tabulation           */
  ad->getIgnDel = 0 ;/* no ignition delay stop  */
  ad->getSens   = 0 ;/* no sensitivity analysis */

  /* read setup file */
  FILE *fsetup = fopen("input.dat","r") ;
  char key[20] ;
  *specinno  = 0 ;
  iod->specoutno = 0 ;
  while ( fscanf(fsetup, "%s", key) != EOF ) {
    if ( strncmp(key,"END",3) == 0 ) break ;
    if ( strncmp(key,"NiterMax"    , 8) == 0 ) fscanf(fsetup, "%d" , &(ad->NiterMax) ) ;
    if ( strncmp(key,"oFreq"       , 5) == 0 ) fscanf(fsetup, "%d" , &(ad->oFreq)    ) ;
    if ( strncmp(key,"withTab"     , 7) == 0 ) fscanf(fsetup, "%d" , withTab         ) ;
    if ( strncmp(key,"getIgnDel"   , 9) == 0 ) fscanf(fsetup, "%d" , &(ad->getIgnDel)) ;
    if ( strncmp(key,"getSens"     , 7) == 0 ) fscanf(fsetup, "%d" , &(ad->getSens)  ) ;
    if ( strncmp(key,"Tini"        , 4) == 0 ) fscanf(fsetup, "%le", &(ad->Tini)     ) ;
    if ( strncmp(key,"Temp_id"     , 7) == 0 ) fscanf(fsetup, "%le", &(ad->Temp_id)  ) ;
    if ( strncmp(key,"deltatMax"   , 9) == 0 ) fscanf(fsetup, "%le", &(ad->deltatMax)) ;
    if ( strncmp(key,"deltat"      , 6) == 0 ) fscanf(fsetup, "%le", &(ad->deltat)   ) ;
    if ( strncmp(key,"deltaTemp"   , 9) == 0 ) fscanf(fsetup, "%le", &(ad->deltaTemp)) ;
    if ( strncmp(key,"tEnd"        , 4) == 0 ) fscanf(fsetup, "%le", &(ad->tend)     ) ;
    if ( strncmp(key,"pfac"        , 4) == 0 ) fscanf(fsetup, "%le", &(ad->pfac)     ) ;
    if ( strncmp(key,"mech"        , 4) == 0 ) fscanf(fsetup, "%s" , mechfile   ) ;
    if ( strncmp(key,"thermo"      , 6) == 0 ) fscanf(fsetup, "%s" , thermofile ) ;

    if ( strncmp(key,"getDeltatSeq",12) == 0 ) fscanf(fsetup, "%d" , &ad->getDeltatSeq ) ;
    if ( strncmp(key,"filDeltatSeq",12) == 0 ) fscanf(fsetup, "%s" , filDeltatSeq ) ;

    if ( strncmp(key,"spec"      ,4) == 0 )
    {
      fscanf(fsetup, "%s%le" ,&SpecName[(*specinno)*lsnm],&SpecMsFr[*specinno]) ;
      (*specinno)++ ;
    }

    if ( strncmp(key,"sout"      ,4) == 0 )
    {
      fscanf(fsetup, "%s" ,&(iod->SpecOutName[(iod->specoutno)*lsnm])) ;
      (iod->specoutno)++ ;
    }

  }
  fclose(fsetup) ;

  /* Convert species names to upper case */
  for ( i = 0 ; i<(*specinno)*lsnm ; i++)
    SpecName[i] = toupper(SpecName[i]);
  for ( i = 0 ; i<(iod->specoutno)*lsnm ; i++)
    iod->SpecOutName[i] = toupper(iod->SpecOutName[i]);

  printf("------------------------------------------\n") ;
  printf("Run parameters : \n") ;
  printf("    NiterMax  = %d\n"     ,ad->NiterMax  ) ;
  printf("    oFreq     = %d\n"     ,ad->oFreq     ) ;
  printf("    Tini      = %20.12e\n",ad->Tini      ) ;
  printf("    Temp_id   = %20.12e\n",ad->Temp_id   ) ;
  printf("    deltat    = %20.12e\n",ad->deltat    ) ;
  printf("    deltatMax = %20.12e\n",ad->deltatMax ) ;
  printf("    deltaTemp = %20.12e\n",ad->deltaTemp ) ;
  printf("    tEnd      = %20.12e\n",ad->tend      ) ;
  printf("    pfac      = %20.12e\n",ad->pfac      ) ;
  printf("    Kinetic model : %s\n",mechfile   ) ;
  printf("    Thermo props  : %s\n",thermofile ) ;
  printf("    Create tables : %d\n",*withTab    ) ;
  printf("    Do ign del    : %d\n",ad->getIgnDel ) ;
  printf("    Do sensitivity: %d\n",ad->getSens   ) ;
  printf("    Species mole fractions : \n"    ) ;
  for ( i = 0 ; i<*specinno ; i++)
    if ( fabs(SpecMsFr[i]) > 1.e-15 )
      printf("    X_%s = %20.12e\n",&SpecName[i*lsnm],SpecMsFr[i]);
  printf("------------------------------------------\n") ;

  /* Read sequence of time steps from file if necessary */
  if ( (ad->getDeltatSeq) == 1 )
  {
    FILE *fseq = NULL ;
    fseq = fopen(filDeltatSeq,"r");
    if ( fseq != NULL )
    {
      fscanf(fseq,"%d",&ad->getDeltatSeq);
      ad->deltatSeq = (double *) malloc((ad->getDeltatSeq)*sizeof(double)) ;
      for ( i=0; i<ad->getDeltatSeq; i++ ) fscanf(fseq,"%le",&ad->deltatSeq[i]) ;
      fclose(fseq) ;

      printf("There are %d entries in the list of deltat : \n",ad->getDeltatSeq);
      for ( i=0; i<ad->getDeltatSeq; i++ ) printf("%d: %24.16e\n",i,ad->deltatSeq[i]) ;
      printf("------------------------------------------\n") ;

    }
    else
    {
      printf("Error in setup() : Could not open %s -> Abort\n",filDeltatSeq) ;
      exit(-1);
    }
  } /* done if (*getDeltatSeq) == 1 */

  return ;

}
