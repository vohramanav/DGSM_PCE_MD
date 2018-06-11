#include "ign.h"

double doIgnDtSeq(tadvInfo *ad, ioInfo *iod, void *cvode_mem)
{

  double tret; 
  int    i, iter, cvflag, Nvars;
  double time_id1 = 0.0 ;

  Nvars = NV_LENGTH_S(ad->y0) ;

  ad->t = ad->tsta ;
  iter = 0    ;
  while ( iter < ad->getDeltatSeq ) 
  {

/*#ifdef DONOTDOIT*/
    /* Re-initialize cvode */
    cvflag = CVodeReInit(cvode_mem, ad->t, ad->y0 );
    Check_CVflag(&cvflag, "CVodeReInit", 1) ;

    cvflag = CVodeSVtolerances(cvode_mem, ad->relT, ad->absT);
    Check_CVflag(&cvflag, "CVodeSVtolerances", 1) ;

    cvflag = CVodeSetStopTime(cvode_mem, ad->tend);
    Check_CVflag(&cvflag, "CVodeSetStopTime", 1) ;
/*#endif*/

    /* call cvode */
    ad->deltat = ad->deltatSeq[iter] ;
    ad->t     += ad->deltat ;
    iter += 1 ;
    cvflag = CVode(cvode_mem, ad->t, ad->y0, &tret, CV_NORMAL);

    if (cvflag == CV_ROOT_RETURN) 
    {
      int rootsfound[1] ;
      int flagroot = CVodeGetRootInfo(cvode_mem, rootsfound);
      Check_CVflag(&flagroot, "CVodeGetRootInfo", 1) ;
      printf("Found root : %20.12e\n",tret);
      time_id1 = tret ;
      if ( ad->getIgnDel == 1 ) { ad->tend = tret ; (ad->NiterMax) = iter ; return ( time_id1 ) ; }
    }

    /* reset time based on end time from cvode */
    ad->t = tret ;

#ifndef NO_OUTPUT
    /* output solution */
    if ( iter % (ad->oFreq) == 0 )
      Output(1, iter, ad, iod) ;
#endif

    /* re-calculate absolute tolerances based on recent solution */
    for( i = 0; i < Nvars ; i++)
    {
      double yscal = MAX( NV_Ith_S(ad->y0,i), 0.0 ) ;
      NV_Ith_S(ad->absT,i) = MAX( ad->relT*yscal, ad->CVsmall ) ;
    }

  } /* end of while loop */

  (ad->NiterMax) = iter ; 
  (ad->tend)     = tret ;

  return ( time_id1 ) ;

} /* end of doIgn */
