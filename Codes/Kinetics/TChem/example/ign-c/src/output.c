#include "ign.h"

static double *tempNmsfr=NULL, *rhsvals=NULL ;

void Output(int iflag, int iter, tadvInfo *ad, ioInfo *iod)
{
  int i, ierr, Nvars, Nspec, icount, lsnm ;
  double sumY, *scal;


  lsnm  = LENGTHOFSPECNAME;
  scal  = N_VGetArrayPointer(ad->y0);
  Nvars = NV_LENGTH_S(ad->y0) ;
  Nspec = Nvars-1;

  if ( tempNmsfr == NULL ) 
    tempNmsfr = (double *) malloc( (Nspec+1) * sizeof(double) ) ;
  if ( rhsvals   == NULL ) 
    rhsvals   = (double *) malloc( (Nspec+1) * sizeof(double) ) ;

  if ( iflag == 0 ) 
  {
    char *snameshdr = (char *) malloc( Nspec * lsnm * sizeof(char) ) ; 
    ierr = TC_getSnames( Nspec, snameshdr ) ;
    if (ierr != 0) {
      printf("ERROR in %s at line %d\n",__FILE__,__LINE__); 
      exit(1);
    }

    fprintf(iod->myfile,"# 1: iteration, #2: time, #3: deltat, #4: temperature,\n") ;
    if (iod->specOutIDs==NULL) {
      for ( i = 0; i<Nspec; i++ )
        fprintf(iod->myfile,"#%d: Mass fraction of %s \n",i+5,&snameshdr[i*lsnm]);
    }
    else {
      icount = 0;
      for ( i = 0; i<Nspec; i++ ) 
        if (iod->specOutIDs[i]==1) {
          fprintf(iod->myfile,"#%d: Mass fraction of %s \n",icount+5,&snameshdr[i*lsnm]);
          icount++;
	}
    }

    fprintf(iod->myfile1,"#1: time, #2: temperature,\n") ;
    if (iod->specOutIDs==NULL) {
      for ( i = 0; i<Nspec; i++ )
        fprintf(iod->myfile1,"#%d: Mass fraction of %s \n",
                i+3,&snameshdr[i*lsnm]);
    }
    else {
      icount = 0;
      for ( i = 0; i<Nspec; i++ ) 
        if (iod->specOutIDs[i]==1) {
          fprintf(iod->myfile1,"#%d: Mass fraction of %s \n",
                  icount+3,&snameshdr[i*lsnm]);
          icount++;
	}
    }

    fprintf(iod->myfile2,"#1: time, \n") ;
    if (iod->specOutIDs==NULL) {
      for ( i = 0; i<Nspec; i++ )
        fprintf(iod->myfile2,"#%d: Molar concentration of %s [kmol/m3]\n",
               i+2,&snameshdr[i*lsnm]);
    }
    else {
      icount=0;
      for ( i = 0; i<Nspec; i++ ) 
        if (iod->specOutIDs[i]==1) {
          fprintf(iod->myfile2,"#%d: Molar concentration of %s [kmol/m3]\n",
                 icount+2,&snameshdr[i*lsnm]);
	  icount++;
	}
    }
    
    fprintf(iod->myfile3,"# 1: time, #2: specific enthalpy [J/kg]") ;
    free(snameshdr) ;
    return ;
  }

  /* Integrate T+(Nspec) */
  sumY=0.0;
  for (i = 1; i < Nvars; i++) sumY += NV_Ith_S(ad->y0,i) ;

  /* Output */
  printf(" Iter = %-7d,  t[s] = %14.6e, dt[s] = %14.6e, T[K] = %14.6e, sY-1.0 = %14.6e\n",
         iter, ad->t, ad->deltat, NV_Ith_S(ad->y0,0), sumY-1.0 );

  fprintf(iod->myfile,"%-10d  %20.12e  %20.12e", iter, ad->t, ad->deltat );
  if (iod->specOutIDs==NULL) {
    for ( i = 0 ; i<Nvars ; i++) fprintf(iod->myfile,"  %20.12e", scal[i] );
  } else {
    fprintf(iod->myfile,"  %20.12e", NV_Ith_S(ad->y0,0) );
    for ( i = 0; i<Nspec; i++ )
      if (iod->specOutIDs[i]==1)
	fprintf(iod->myfile,"  %20.12e", scal[i+1] );
  }
  fprintf(iod->myfile,"\n");

  fprintf(iod->myfile1,"%20.12e",ad->t);
  if (iod->specOutIDs==NULL) {
    for ( i = 0 ; i<Nvars ; i++) fprintf(iod->myfile1,"  %20.12e", scal[i] );
  } else {
    fprintf(iod->myfile,"  %20.12e", NV_Ith_S(ad->y0,0) );
    for ( i = 0; i<Nspec; i++ )
      if (iod->specOutIDs[i]==1)
	fprintf(iod->myfile1,"  %20.12e", scal[i+1] );
  }
  fprintf(iod->myfile1,"\n");
  
  TC_getMs2Cc ( scal, Nvars, &tempNmsfr[1] ) ;
  fprintf(iod->myfile2,"%20.12e  %20.12e", ad->t, scal[0]);
  if (iod->specOutIDs==NULL) {
    for ( i = 1 ; i<Nspec+1 ; i++) fprintf(iod->myfile2,"  %20.12e",tempNmsfr[i] );
  } else {
    for ( i = 0; i<Nspec; i++ )
      if (iod->specOutIDs[i]==1)
        fprintf(iod->myfile2,"  %20.12e",tempNmsfr[i+1] );
  }
  fprintf(iod->myfile2,"\n");

  TC_getMs2HmixMs ( scal, Nvars, &tempNmsfr[0] ) ;
  fprintf(iod->myfile3,"%20.12e  %20.12e", ad->t, tempNmsfr[0]);

  return ;

}
