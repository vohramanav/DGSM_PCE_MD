/**
 * \ingroup srcs
 * \brief Returns dimensional/non-dimensional source term for 
 * \f[\frac{\partial T}{\partial t}=\omega_0,\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 *  based on temperature T and species mass fractions Y's
 */

int TCDND_getSrc(double *scal,int Nvars,double *omega)
{

/**
   \param scal : array of Nspec+1 doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables \f$ N_{vars}= N_{spec}+1\f$ 
   \return omega : array of N<sub>spec</sub>+1 source terms (possibly normalized) for temperature and 
                   species mass fractions equations: 
		   omega[0] : [(K/s)], omega[1...N<sub>spec</sub>] : [(1/s)]
*/

  int ans, i ;

  if ( Nvars != TC_Nvars_ ) 
    TC_errorMSG( 20, "TCDND_getSrc", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getSrc(scal, Nvars, omega) ;

  if ( TC_nonDim_ == 1 ) 
  {
    scal[0] /= TC_Tref_ ;
    omega[0] *= (TC_timref_/TC_Tref_) ;
    for ( i = 1 ; i < TC_Nspec_+1 ; i++ ) omega[i] *= TC_timref_ ;
  }

  return ( ans ) ;
 
}
/*
                _   ____            ____                
      __ _  ___| |_/ ___| _ __ ___ / ___|___  _ __  ___ 
     / _` |/ _ \ __\___ \| '__/ __| |   / _ \| '_ \/ __|
    | (_| |  __/ |_ ___) | | | (__| |__| (_) | | | \__ \
     \__, |\___|\__|____/|_|  \___|\____\___/|_| |_|___/
     |___/                                              

*/
/**
 * \ingroup srcs
 * \brief Returns source term (dimensional/non-dimensional) for 
 * \f[\frac{\partial \rho}{\partial t}=\omega_0,\rho\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 * based on \f$\rho\f$ and Y's
 */
int TCDND_getSrcCons(double *scal,int Nvars,double *omega)
{

/**
   \param scal : array of Nspec+1 doubles \f$(\rho,Y_1,Y_2,...,Y_N)\f$
	        density \f$[kg/m^3]\f$, mass fractions Y [] 
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return omega : array of N<sub>spec</sub>+1 source terms for density and species mf
                   conservative formulation: 
		   omega[0] : \f$[kg/(m^3\cdot s)]\f$, 
		   omega[1...N<sub>spec</sub>] : \f$[kg/(m^3\cdot s)]\f$
*/

  int ans, i ;

  if ( Nvars != TC_Nvars_ ) 
    TC_errorMSG( 20, "TCDND_getSrcCons", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getSrcCons(scal, Nvars, omega) ;

  if ( TC_nonDim_ == 1 ) 
  {

    scal[0] /= TC_Tref_ ;
    omega[0] *= (TC_timref_/TC_rhoref_) ;
    for ( i = 1 ; i < TC_Nspec_+1 ; i++ ) omega[i] *= TC_timref_ ;

  }

  return ( ans ) ;

}
/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian (dimensional/non-dimensional) for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_{N-1})\f$ 
 * based on temperature T and species mass fractions Y's
 */
int TCDND_getJacTYNm1anl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TCDND_getJacTYNm1anl", Nspec, TC_Nspec_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getJacTYNm1anl( scal, Nspec, jac) ;

  if ( TC_nonDim_ == 1 ) 
  {
    scal[0] /= TC_Tref_ ;
    for ( i=0; i<Nspec*Nspec; i++ ) jac[i] *= TC_timref_ ;
    /* first column, except first element */
    for ( i=1 ; i < Nspec ; i++ ) jac[i] *= TC_Tref_ ;
    /* first line, except first element*/
    for ( i=Nspec ; i < Nspec*Nspec ; i+=Nspec ) jac[i] /= TC_Tref_ ;
  }

  return ( ans ) ;

} /* done with TCDND_getJacTYNm1anl */

/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian (dimensional/non-dimensional) for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
*/
int TCDND_getJacTYNanl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TCDND_getJacTYNm1anl", Nspec, TC_Nspec_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getJacTYNanl( scal, Nspec, jac) ;

  if ( TC_nonDim_ == 1 ) 
  {
    scal[0] /= TC_Tref_ ;
    for ( i=0; i<(Nspec+1)*(Nspec+1); i++ ) jac[i] *= TC_timref_ ;
    /* first column, except first element */
    for ( i=1 ; i < Nspec+1 ; i++ ) jac[i] *= TC_Tref_ ;
    /* first line, except first element */
    for ( i=Nspec+1 ; i < (Nspec+1)*(Nspec+1) ; i+=(Nspec+1) ) jac[i] /= TC_Tref_ ;
  }

  return ( ans ) ;

} /* done with TCDND_getJacTYNanl */

/**
 * \ingroup jacs
 * \brief Computes (analytical or numerical) Jacobian (dimensional/non-dimensional) 
 * for the system \f$(T,Y_1,Y_2,\ldots,Y_{N-1})\f$ 
 * based on temperature T and species mass fractions Y's
*/
int TCDND_getJacTYNm1(double *scal, int Nspec, double *jac, unsigned int useJacAnl)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \param useJacAnl : flag for Jacobian type (1-analytical,other values-numerical) 
 *  \return jac  : Jacobian array 
 */

  int ans, i ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TCDND_getJacTYNm1", Nspec, TC_Nspec_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getJacTYNm1( scal, Nspec, jac, useJacAnl) ;

  if ( TC_nonDim_ == 1 ) 
  {
    scal[0] /= TC_Tref_ ;
    for ( i=0; i<Nspec*Nspec; i++ ) jac[i] *= TC_timref_ ;
    /* first column, except first element */
    for ( i=1 ; i < Nspec ; i++ ) jac[i] *= TC_Tref_ ;
    /* first line, except first element */
    for ( i=Nspec ; i < Nspec*Nspec ; i+=Nspec ) jac[i] /= TC_Tref_ ;
  }

  return ( ans ) ;

} /* done TCDND_getJacTYNm1 */

/**
 * \ingroup jacs
 * \brief Computes (analytical or numerical) Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
 */
int TCDND_getJacTYN(double *scal, int Nspec, double *jac, unsigned int useJacAnl)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \param useJacAnl : flag for Jacobian type (1-analytical,other values-numerical) 
 *  \return jac  : Jacobian array 
 */

  int ans, i ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TCDND_getJacTYN", Nspec, TC_Nspec_ ) ;
	     
  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getJacTYN( scal, Nspec, jac, useJacAnl) ;

  if ( TC_nonDim_ == 1 ) 
  {
    scal[0] /= TC_Tref_ ;
    for ( i=0; i<(Nspec+1)*(Nspec+1); i++ ) jac[i] *= TC_timref_ ;
    /* first column, except first element */
    for ( i=1 ; i < Nspec+1 ; i++ ) jac[i] *= TC_Tref_ ;
    /* first line, except first element */
    for ( i=Nspec+1 ; i < (Nspec+1)*(Nspec+1) ; i+=(Nspec+1) ) jac[i] /= TC_Tref_ ;
  }

  return ( ans ) ;

} /* done TCDND_getJacTYN */
