#include "ign.h"

double doIgn(tadvInfo *ad, ioInfo *iod, void *cvode_mem)
{

  int    i, iter, cvflag, Nvars, foundid2;
  double tret, der2, time_id1, time_id2, Temp_m1, Temp_m2, time_m1, time_m2, Temp_save ;

  foundid2 = 0 ;
  time_id1 = time_id2 = 0.0     ;
  time_m1  = time_m2  = 0.0     ;
  Temp_m1  = Temp_m2  = NV_Ith_S(ad->y0,0);
  iter = 0    ;

  Nvars = NV_LENGTH_S(ad->y0) ;

  ad->t = ad->tsta ;

  /* Re-initialize cvode */
  cvflag = CVodeReInit(cvode_mem, ad->tsta, ad->y0 );
  Check_CVflag(&cvflag, "CVodeReInit", 1) ;
  cvflag = CVodeSVtolerances(cvode_mem, ad->relT, ad->absT);
  Check_CVflag(&cvflag, "CVodeSVtolerances", 1) ;
  cvflag = CVodeSetStopTime(cvode_mem, ad->tend);
  Check_CVflag(&cvflag, "CVodeSetStopTime", 1) ;

  //printf("%18.12e\n",ad->relT);
  // for (i=0; i<Nvars;i++)
  //  printf("%d %18.12e\n",i,NV_Ith_S(ad->absT,i));

  //exit(0);

  Temp_save = NV_Ith_S(ad->y0,0) ;

#ifdef DEBUGMSG
  N_Vector ydot = N_VNew_Serial( Nvars );
  chemrhs(0.0, ad->y0, ydot, NULL);
  for( i = 0; i < Nvars ; i++) {
    printf("%d %e %e\n",i,NV_Ith_S(ad->y0,i),NV_Ith_S(ydot,i));
  }
  exit(1);
#endif
  while ( ( ad->t < ad->tend ) && (iter < (ad->NiterMax) ) ) {
    /* call cvode */
    ad->t += (ad->deltat) ;
    iter  += 1         ;
    cvflag = CVode(cvode_mem, ad->t, ad->y0, &tret, CV_NORMAL);

    if (cvflag == CV_ROOT_RETURN)
    {
      int rootsfound[1] ;
      int flagroot = CVodeGetRootInfo(cvode_mem, rootsfound);
      Check_CVflag(&flagroot, "CVodeGetRootInfo", 1) ;
#ifdef VERBOSE
      printf("Found root : %20.12e\n",tret);
#endif
      time_id1 = tret ;
      if ( ad->getIgnDel == 1 ) { ad->tend = tret ; (ad->NiterMax) = iter ; return ( time_id1 ) ; }
    }

    /* reset time based on end time from cvode */
    ad->t = tret ;

    /*
       calculate new time step:
       - reduce/increase time step if \Delta T >< 0.25K
     */
    (ad->deltat) = MIN(
		       ad->deltaTemp/MAX(fabs(NV_Ith_S(ad->y0,0)-Temp_save),TCSMALL)*(ad->deltat),
		       2.0*(ad->deltat)
		       );

     /* reset time step based on set limits */
     /* if ( iter < 100 ) */
     /*   (*deltat) = MIN((*deltat),1.e-4) ; */
     /* else */
     ad->deltat = MIN(ad->deltat,ad->deltatMax) ;
     Temp_save = NV_Ith_S(ad->y0,0) ;

#ifndef NO_OUTPUT
    /* output solution */
     if ( iter % (ad->oFreq) == 0 )
    {
      Output(1, iter, ad, iod) ;
    }
#endif

    /* re-calculate absolute tolerances based on recent solution */
    for( i = 0; i < Nvars ; i++)
    {
      double yscal = MAX( NV_Ith_S(ad->y0,i), 0.0 ) ;
      NV_Ith_S(ad->absT,i) = MAX( ad->relT*yscal, ad->CVsmall ) ;
    }
    /*
  printf("%18.12e\n",ad->relT);
  for (i=0; i<Nvars;i++)
    printf("%d %18.12e\n",i,NV_Ith_S(ad->absT,i));
    */
    /* Ignition delay time - second derivative */
    if ( ( foundid2 == 0 ) && ( NV_Ith_S(ad->y0,0) > 1300.0 ) )
    {
      der2 = (NV_Ith_S(ad->y0,0)-Temp_m1)/(ad->t-time_m1)-(Temp_m1-Temp_m2)/(time_m1-time_m2) ;
      if ( der2 <= 0.0 )
      {
    	foundid2 = 1 ;
    	time_id2 = time_m1 ;
	printf("Found time_id2 : %20.12e\n",time_id2) ;
      }
    }

    /* Cycle temperatures and times */
    Temp_m2 = Temp_m1            ;
    Temp_m1 = NV_Ith_S(ad->y0,0) ;
    time_m2 = time_m1            ;
    time_m1 = ad->t              ;

  } /* */

  ad->NiterMax = iter ;
  ad->tend     = tret ;

  return ( time_id1 ) ;

} /* end of doIgn */
