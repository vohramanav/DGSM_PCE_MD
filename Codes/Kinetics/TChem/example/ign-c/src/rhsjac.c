#include "ign.h"

void chemrhs(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{

  int Nv,i ;

  Nv = NV_LENGTH_S(y) ;

  realtype *TandYs = N_VGetArrayPointer(y   );
  realtype *omg    = N_VGetArrayPointer(ydot);
  TC_getSrc ( TandYs, Nv, omg ) ;

  return ;

}

int chemrhswrapper(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{

  chemrhs(t,y,ydot,f_data) ;
  return ( 0 );

}

#ifdef USEJAC

void chemjac(realtype t, N_Vector y, DlsMat Jac, double *jactmp) 
{

  int Ns,i ;
  unsigned int useJacAnl = 1 ; /* use analytic Jacobian */

  Ns = NV_LENGTH_S(y)-1 ;
  realtype *TandYs = N_VGetArrayPointer(y   );
  TC_getJacTYN ( TandYs, Ns, (Jac->data), useJacAnl ) ;

  //for (i=0; i<1; i++)
  //  printf("%d %18e \n",i,(Jac->data)[i]);

  return ;
  
}

int chemjacwrapper(long int N, realtype t, 
                   N_Vector y, N_Vector fy, DlsMat J, void *udata,
		               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) 
{

  double *jactmp = (double*) udata ;
  chemjac( t, y, J, jactmp) ;

  return ( 0 );

}

#endif

int g(realtype t, N_Vector y, realtype *gout, void *user_data)
{

  cdt *udataloc = (cdt *)user_data;
  realtype temper = NV_Ith_S(y,0) ;

  gout[0] = temper-(udataloc->Tign) ;
  udataloc->ncalls += 1;

  return(0);
}



int Check_CVflag(void *flagvalue, char *funcname, int opt)
{

  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if ( ( opt == 0 ) && ( flagvalue == NULL ) )
  {
    printf("CVODE_ERROR: %s failed - returned NULL pointer\n",funcname) ;
    exit(1) ;
  }

  /* Check if flag < 0 */
  else if ( opt == 1 )
  {
    errflag = (int *) flagvalue;
    if ( *errflag < 0 )
    {
      printf("CVODE_ERROR: %s failed with flag = %d\n", funcname, *errflag );
      exit(1) ;
    }
  }

  return ( 0 ) ;

}
