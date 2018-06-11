#include "ign.h"

double doIgnReinit(tadvInfo *ad, ioInfo *iod, void *cvode_mem)
{

  double tret, Temp_save, time_id1 = 0.0 ; 
  int    i, iter, cvflag, Nvars;

  Nvars = NV_LENGTH_S(ad->y0) ;

  ad->t = ad->tsta ;
  iter = 0    ;
  time_id1 = -100.0 ;

  /* Re-initialize cvode */
  cvflag = CVodeReInit(cvode_mem, ad->tsta, ad->y0 );
  Check_CVflag(&cvflag, "CVodeReInit", 1) ;
  cvflag = CVodeSVtolerances(cvode_mem, ad->relT, ad->absT);
  Check_CVflag(&cvflag, "CVodeSVtolerances", 1) ;
  cvflag = CVodeSetStopTime(cvode_mem, ad->tend);
  Check_CVflag(&cvflag, "CVodeSetStopTime", 1) ;

  Temp_save = NV_Ith_S(ad->y0,0);
  while ( ( ad->t < (ad->tend) ) && (iter < (ad->NiterMax) ) )
  {
    /* call cvode */
    ad->t += (ad->deltat) ;
    iter  += 1         ;
    cvflag = CVode(cvode_mem, ad->t, ad->y0, &tret, CV_NORMAL);

    if (cvflag == CV_ROOT_RETURN) 
    {
      int rootsfound[1] ;
      int flagroot = CVodeGetRootInfo(cvode_mem, rootsfound);
      Check_CVflag(&flagroot, "CVodeGetRootInfo", 1) ;
      /* printf("Found root : %20.12e\n",tret); */
      time_id1 = tret ;
      return ( time_id1 ) ; 
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
    (ad->deltat) = MIN((ad->deltat),ad->deltatMax) ;
    Temp_save = NV_Ith_S(ad->y0,0);

    /* re-calculate absolute tolerances based on recent solution */
    for( i = 0; i < Nvars ; i++)
    {
      double yscal = MAX( NV_Ith_S(ad->y0,i), 0.0 ) ;
      NV_Ith_S(ad->absT,i) = MAX( ad->relT*yscal, ad->CVsmall ) ;
    }

  }

  return ( time_id1 ) ;

} /* end of doIgn */
