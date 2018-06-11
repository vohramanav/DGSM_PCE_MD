/*! \file TC_src.c
 *    \brief Jacobian functions
 */

/*
           _                 _     _                 
          | | __ _  ___ ___ | |__ (_) __ _ _ __  ___ 
       _  | |/ _` |/ __/ _ \| '_ \| |/ _` | '_ \/ __|
      | |_| | (_| | (_| (_) | |_) | | (_| | | | \__ \
       \___/ \__,_|\___\___/|_.__/|_|\__,_|_| |_|___/

 */
/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_{N-1})\f$ 
 * based on T and Y's
 */
int TC_getJacTYNm1anl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i, k, iJacR,iJacF ;
  double rhomix, Wmix, drhodT, *Yspec ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacTYNm1anl", Nspec, TC_Nspec_ ) ;

  ans   = 0 ;
  Yspec = &scal[1] ;

  /* first, get full Jacobian */
  ans = TC_getJacRPTYNanl(scal, Nspec, TC_jacFull) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reduced Jacobian */
  ans = TC_getRhoMixMs(  scal, TC_Nvars_, &rhomix ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix ( Yspec, TC_Nspec_, &Wmix   ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* first column */
  drhodT = -rhomix/scal[0] ;
  for ( i = 0, iJacR=0, iJacF = TC_Nvjac_*2 ; 
	i < Nspec ; 
	i++,iJacR++,iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+2]+TC_jacFull[iJacF]*drhodT ;

  /* rest of the columns */
  for ( k = 0 ; k < Nspec-1 ; k++ )
  {
    double drhodYk = - rhomix * Wmix * ( 1.0 / TC_sMass_[k] - 1.0 / TC_sMass_[Nspec-1] ) ;
    double dYndYk  = - 1.0 ;
    for ( i = 0,iJacR=Nspec*(k+1),iJacF = TC_Nvjac_*2 ; i < Nspec ; 
	  i++, iJacR++, iJacF+=TC_Nvjac_ )
      jac[iJacR] = TC_jacFull[iJacF+3+k] + TC_jacFull[iJacF]*drhodYk 
                 + TC_jacFull[iJacF+3+Nspec-1]*dYndYk ;
  }

  return ( ans ) ;

} /* done with TC_getJacTYNm1anl */

/*
               _      _           _______   ___   _             _ 
     __ _  ___| |_   | | __ _  __|_   _\ \ / / \ | | __ _ _ __ | |
    / _` |/ _ \ __|  | |/ _` |/ __|| |  \ V /|  \| |/ _` | '_ \| |
   | (_| |  __/ || |_| | (_| | (__ | |   | | | |\  | (_| | | | | |
    \__, |\___|\__\___/ \__,_|\___||_|   |_| |_| \_|\__,_|_| |_|_|
    |___/                                                         
*/

/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
*/
int TC_getJacTYNanl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i, k, iJacR,iJacF ;
  double rhomix, Wmix, drhodT, *Yspec ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacTYNm1anl", Nspec, TC_Nspec_ ) ;

  ans   = 0 ;
  Yspec = &scal[1] ;

  /* first, get full Jacobian */
  ans = TC_getJacRPTYNanl(scal, Nspec, TC_jacFull) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reduced Jacobian */

  ans = TC_getRhoMixMs(  scal, TC_Nvars_, &rhomix ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix ( Yspec, TC_Nspec_, &Wmix   ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* first column */
  drhodT = -rhomix/scal[0] ;
  for ( i = 0, iJacR=0, iJacF = TC_Nvjac_*2 ; 
	i < Nspec+1 ; 
	i++,iJacR++,iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+2]+TC_jacFull[iJacF]*drhodT ;

  /* rest of the columns */
  for ( k = 0 ; k < Nspec; k++ )
  {
    double drhodYk = - rhomix * Wmix / TC_sMass_[k]  ;
    for ( i = 0,iJacR=(Nspec+1)*(k+1),iJacF = TC_Nvjac_*2 ; i < Nspec+1 ; 
	  i++, iJacR++, iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+3+k] + TC_jacFull[iJacF]*drhodYk ;
  }

  return ( ans ) ;

} /* done with TC_getJacTYNanl */

/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
*/
int TC_getJacDiagTYNanl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i, k, iJacR,iJacF ;
  double rhomix, Wmix, drhodT, *Yspec ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacTYNm1anl", Nspec, TC_Nspec_ ) ;

  ans   = 0 ;
  Yspec = &scal[1] ;

  /* first, get full Jacobian */
  ans = TC_getJacRPTYNanl(scal, Nspec, TC_jacFull) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reduced Jacobian */

  ans = TC_getRhoMixMs(  scal, TC_Nvars_, &rhomix ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix ( Yspec, TC_Nspec_, &Wmix   ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* first column */
  drhodT = -rhomix/scal[0] ;
  iJacF = TC_Nvjac_*2 ;
  jac[0] = TC_jacFull[iJacF+2]+TC_jacFull[iJacF]*drhodT ;

  /* rest of the columns */
  for ( k = 0 ; k < Nspec; k++ ) {
    double drhodYk = - rhomix * Wmix / TC_sMass_[k]  ;
    i = k+1;
    iJacR=(Nspec+1)*(k+1)+i;
    iJacF = TC_Nvjac_*(2 + i);
    jac[iJacR] = TC_jacFull[iJacF+3+k] + TC_jacFull[iJacF] * drhodYk ;
  }

  return ( ans ) ;

} /* done with TC_getJacDiagTYNanl */


/**
 * \ingroup jacs
 * \brief Computes (analytical or numerical) Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_{N-1})\f$ 
 * based on T and Y's
*/
int TC_getJacTYNm1(double *scal, int Nspec, double *jac, unsigned int useJacAnl)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \param useJacAnl : flag for Jacobian type (1-analytical,other values-numerical) 
 *  \return jac  : Jacobian array 
 */

  int ans, i, k, iJacR,iJacF ;
  double rhomix, Wmix, drhodT, *Yspec ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacTYNm1", Nspec, TC_Nspec_ ) ;

  ans   = 0 ;
  Yspec = &scal[1] ;

  /* first, get full Jacobian */
  if ( useJacAnl )
    ans = TC_getJacRPTYNanl(scal, Nspec, TC_jacFull) ;
  else
    ans = TC_getJacRPTYNnum(scal, Nspec, TC_jacFull) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reduced Jacobian */
  ans = TC_getRhoMixMs(  scal, TC_Nvars_, &rhomix ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix ( Yspec, TC_Nspec_, &Wmix   ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* first column */
  drhodT = -rhomix/scal[0] ;
  for ( i = 0,iJacR=0,iJacF = TC_Nvjac_*2 ; 
	i < Nspec ; 
	i++,iJacR++,iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+2]+TC_jacFull[iJacF]*drhodT ;

  /* rest of the columns */
  for ( k = 0 ; k<Nspec-1 ; k++)
  {
    double drhodYk = - rhomix * Wmix * ( 1.0 / TC_sMass_[k] - 1.0 / TC_sMass_[Nspec-1] ) ;
    double dYndYk  = - 1.0 ;
    for ( i = 0,iJacR=Nspec*(k+1),iJacF = TC_Nvjac_*2 ; i < Nspec ; 
	  i++,iJacR++,iJacF+=TC_Nvjac_ )
      jac[iJacR] = TC_jacFull[iJacF+3+k] + TC_jacFull[iJacF]*drhodYk 
      + TC_jacFull[iJacF+3+Nspec-1]*dYndYk ;
  }
 
  return ( ans ) ;

} /* done TC_getJacTYNm1 */

/**
 * \ingroup jacs
 * \brief Computes (analytical or numerical) Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
 */
int TC_getJacTYN(double *scal, int Nspec, double *jac, unsigned int useJacAnl)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \param useJacAnl : flag for Jacobian type (1-analytical,other values-numerical) 
 *  \return jac  : Jacobian array 
 */

  int ans, i, k, iJacR,iJacF ;
  double rhomix, Wmix, drhodT, *Yspec ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacTYN", Nspec, TC_Nspec_ ) ;

  ans   = 0 ;
  Yspec = &scal[1] ;

  /* first, get full Jacobian */
  if ( useJacAnl )
    ans = TC_getJacRPTYNanl(scal, Nspec, TC_jacFull) ;
  else
    ans = TC_getJacRPTYNnum(scal, Nspec, TC_jacFull) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reduced Jacobian */
  ans = TC_getRhoMixMs(  scal, TC_Nvars_, &rhomix ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix ( Yspec, TC_Nspec_, &Wmix   ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* first column */
  drhodT = -rhomix/scal[0] ;
  for ( i = 0, iJacR=0, iJacF = TC_Nvjac_*2 ; 
	i < Nspec+1 ; 
	i++,iJacR++,iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+2]+TC_jacFull[iJacF]*drhodT ;

  /* rest of the columns */
  for ( k = 0 ; k < Nspec; k++ )
  {
    double drhodYk = - rhomix * Wmix / TC_sMass_[k] ;
    for ( i = 0,iJacR=(Nspec+1)*(k+1),iJacF = TC_Nvjac_*2 ; i < Nspec+1 ; 
	  i++, iJacR++, iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+3+k] + TC_jacFull[iJacF]*drhodYk ;
  }

  /* /\* skip rho and pressure lines and columns *\/   */
  /* for (iJacR = 0, iJacF=2*TC_Nvjac_+2, k = 0;       */
  /*      k<Nspec+1;                                   */
  /*      k++, iJacF +=2 )                             */
  /*   for ( i = 0; i<Nspec+1; i++, iJacR++, iJacF++ ) */
  /*     jac[ iJacR ] = TC_jacFull[ iJacF ] ;          */
	     
  return ( ans ) ;

} /* done TC_getJacTYN */

/*
                     _      _            ____  ____ _______   ___   _ 
           __ _  ___| |_   | | __ _  ___|  _ \|  _ \_   _\ \ / / \ | |
          / _` |/ _ \ __|  | |/ _` |/ __| |_) | |_) || |  \ V /|  \| |
         | (_| |  __/ || |_| | (_| | (__|  _ <|  __/ | |   | | | |\  |
          \__, |\___|\__\___/ \__,_|\___|_| \_\_|    |_|   |_| |_| \_|
          |___/                                                       

*/
/**
 * \ingroup jacs
 * \brief Computes (analytical) Jacobian for the system 
 * \f$(\rho,P,T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
*/
int TC_getJacRPTYN(double *scal, int Nspec, double *jac, unsigned int useJacAnl)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \param useJacAnl : flag for Jacobian type (1-analytical,other values-numerical) 
 *  \return jac  : Jacobian array 
 */

  int ans ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacRPTYN", Nspec, TC_Nspec_ ) ;

  if ( useJacAnl )
    ans = TC_getJacRPTYNanl(scal, Nspec, jac) ;
  else
    ans = TC_getJacRPTYNnum(scal, Nspec, jac) ;

  return ( ans ) ;

}

/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian for the system 
 * \f$(\rho,P,T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
*/
int TC_getJacRPTYNanl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i, j, k, kspec, ir, indx, indxR, iJac, itbdy, ipfal, iord ;
  double *Yspec, t1, t_1, tln ;
  double rhomix,cpmix,cpmixP ;
  double sum1, sum2, sum3 ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacRPTYNanl", Nspec, TC_Nspec_ ) ;

  ans = 0 ;

  /* scal array contains T and all species mass fractions */
  /* clear data */
  for ( i = 0 ; i<TC_Nspec_; i++) TC_omg[i] = 0.0 ;

  /* Compute molar concentrations (kmol/m3) and transform them to moles/cm3 */
  ans = TC_getMs2Cc(scal,TC_Nspec_+1,TC_Xconc) ;
  if ( ans != 0 ) return ( ans ) ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) TC_Xconc[i] *= 1.e-3 ;

  /* set species array */
  Yspec = &scal[1] ;

  /* compute functions and powers of temperature */
  t1  = scal[0] ;
  t_1 = 1.0/t1  ;
  tln = log(t1) ;

  /* get 3rd-body concentrations */
  ans = TC_get3rdBdyConc(TC_Xconc,TC_Mconc) ;
  if ( ans != 0 ) return ( ans ) ;

  /* compute (-ln(T)+dS/R-dH/RT) for each species */
  ans = TC_getgk (t1,t_1,tln) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getgkp(t1,t_1,tln) ;
  if ( ans != 0 ) return ( ans ) ;

  /* compute forward and reverse rate constants */
  ans = TC_getkForRev (t1,t_1,tln) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getkForRevP(t1,t_1) ;
  if ( ans != 0 ) return ( ans ) ;

  /* compute rate-of-progress */
  ans = TC_getRateofProg(TC_Xconc) ;
  if ( ans != 0 ) return ( ans ) ;

  /* compute pressure dependent factors */
  ans = TC_getCrnd( t1, t_1, tln, TC_Xconc, TC_Mconc ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reaction rates */
  for ( i = 0 ; i<TC_Nreac_; i++)
  {

    TC_rop[i] *= TC_Crnd[i] ;

    indx  = i*TC_maxSpecInReac_ ;
    for ( j = 0; j<TC_reacNreac_[i] ; j++)
    {
      int kspec = TC_reacSidx_[indx] ;
      TC_omg[kspec] += ((double) TC_reacNuki_[indx])*TC_rop[i] ;
      indx++ ;
    }

    indx  = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
    for ( j = 0; j<TC_reacNprod_[i] ; j++)
    {
      int kspec = TC_reacSidx_[indx] ;
      TC_omg[kspec] += ((double) TC_reacNuki_[indx])*TC_rop[i] ;
      indx++ ;
    }

  } /* done loop over the number of reactions */

  /* check for reactions with real stoichiometric coefficients */
  if ( TC_nRealNuReac_ > 0 )
    for ( ir = 0 ; ir < TC_nRealNuReac_; ir++ ) {

      i = TC_reacRnu_[ir] ;

      indx  = i *TC_maxSpecInReac_ ;
      indxR = ir*TC_maxSpecInReac_ ;      
      for ( j = 0 ; j<TC_reacNreac_[i] ; j++) {
        int kspec = TC_reacSidx_[indx] ;
        TC_omg[kspec] += TC_reacRealNuki_[indxR]*TC_rop[i] ;
        indx++  ;
        indxR++ ;
      }

      indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      //indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;   // Bugfix Cosmin 2016/12/28
      for ( j = 0; j<TC_reacNprod_[i] ; j++) {
        int kspec = TC_reacSidx_[indx] ;
        TC_omg[kspec] += TC_reacRealNuki_[indxR]*TC_rop[i] ;
	indx++  ;
	indxR++ ;
      }
    } /* done loop over the number of reactions */

  /* transform from mole/(cm3.s) to kg/(m3.s) */
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) TC_omg[i] *= (1.e3*TC_sMass_[i]) ;

  /* assemble jacobian terms, F[3+i,2...] */
  for ( i = 0 ; i < TC_Nvjac_* TC_Nvjac_ ; i++) jac[i] = 0.0 ;
  itbdy = 0 ; 
  ipfal = 0 ; 
  iord  = 0 ;
  for (j = 0 ; j < TC_Nreac_; j++) {

    unsigned int arbord ;
    /* compute {\partial Crnd}/{\partial T,Y1,Y2,...,YN} */
    TC_getCrndDer( j, &itbdy, &ipfal, t1, t_1, tln, TC_Xconc, TC_Mconc ) ;

    arbord = 0 ;
    if ( iord < TC_nOrdReac_ )
      if (TC_reacAOrd_[iord] == j) arbord = 1 ;

    if ( TC_reacScoef_[j] == -1 ) {

      /* reaction has integer stoichiometric coefficients */
      for ( i = 0 ; i < TC_Nspec_ ; i++) {
      
        /* skip if species i is not involved in reaction j */
        int indxIJ = j*TC_Nspec_+i ;
        if ( TC_NuIJ[indxIJ] == 0 ) continue ;
      
        /* add 1st term for the T,Y1,Y2,...,YN - derivatives */
        iJac = TC_Nvjac_*(i+3)+2 ; 
        for ( k = -1 ; k < TC_Nspec_ ; k++,iJac++)
          jac[iJac] += TC_NuIJ[indxIJ]*TC_CrndDer[k+1]*(TC_ropFor[j]-TC_ropRev[j]) ;

        /* add 2nd term for T-derivative F[3+i,2] */
        iJac = TC_Nvjac_*(i+3)+2 ; 
        jac[iJac] += TC_NuIJ[indxIJ]*TC_Crnd[j]*(TC_ropFor[j]*TC_kforP[j]-TC_ropRev[j]*TC_krevP[j]) ;
      
        /* add 2nd term for species derivatives F[3+i,3+k] */
        if ( !arbord ) {
          indx = j*TC_maxSpecInReac_ ;
          for ( kspec = 0; kspec<TC_reacNreac_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = TC_Nvjac_*(i+3)+3+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] += TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            indx++ ;
          } /* done loop over species k in reactants (2nd term) */

          indx  = j*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
          for ( kspec = 0; kspec<TC_reacNprod_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = TC_Nvjac_*(i+3)+3+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] -= TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            indx++ ;
          } /* done loop over species k in products (2nd term) */
        } /* done section for non-arbitrary order reactions */
        else {
          /* arbitrary order reaction */
          indx = iord*TC_maxOrdPar_ ;
          for ( kspec = 0; kspec<TC_maxOrdPar_ ; kspec++) {
            if (TC_specAOidx_[indx]<0) {
              k = -TC_specAOidx_[indx]-1 ;
              iJac = TC_Nvjac_*(i+3)+3+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] += TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            }
            else if (TC_specAOidx_[indx]>0) {
              k = TC_specAOidx_[indx]-1 ;
              iJac = TC_Nvjac_*(i+3)+3+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] -= TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            }
            indx++ ;
          }
	      } /* done section for arbitrary order reaction */
      } /* done loop over species (counter i) for integer stoichiometric coefficients */
    } /* end if integer coeffs */
    else {

      /* reaction has real stoichiometric coefficients */
      for ( i = 0 ; i < TC_Nspec_ ; i++) {
        int jr = TC_reacScoef_[j] ;
        /* skip if species i is not involved in reaction j */
        int indxIJ = jr*TC_Nspec_+i ;
        if ( TC_RealNuIJ[indxIJ] < TCSMALL ) continue ;
      
        /* add 1st term for the T,Y1,Y2,...,YN - derivatives */
        iJac = TC_Nvjac_*(i+3)+2 ; 
        for ( k = -1 ; k < TC_Nspec_ ; k++,iJac++)
          jac[iJac] += TC_RealNuIJ[indxIJ]*TC_CrndDer[k+1]*(TC_ropFor[j]-TC_ropRev[j]) ;

        /* add 2nd term for T-derivative F[3+i,2] */
        iJac = TC_Nvjac_*(i+3)+2 ; 
        jac[iJac] += TC_RealNuIJ[indxIJ]*TC_Crnd[j]*(TC_ropFor[j]*TC_kforP[j]-TC_ropRev[j]*TC_krevP[j]) ;
      
        /* add 2nd term for species derivatives F[3+i,3+k] */
        if ( !arbord ) {

          indx = j *TC_maxSpecInReac_ ;
          for ( kspec = 0; kspec<TC_reacNreac_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = TC_Nvjac_*(i+3)+3+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] += TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            indx++ ;
          } /* done loop over species k in reactants (2nd term) */

          indx = j *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
          for ( kspec = 0; kspec<TC_reacNprod_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = TC_Nvjac_*(i+3)+3+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] -= TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            indx++ ;	  
          } /* done loop over species k in products (2nd term) */
        }
        else {
          /* arbitrary order reaction */
          indx = iord*TC_maxOrdPar_ ;
          for ( kspec = 0; kspec<TC_maxOrdPar_ ; kspec++) {
            if (TC_specAOidx_[indx]<0) {
              k = -TC_specAOidx_[indx]-1 ;
              iJac = TC_Nvjac_*(i+3)+3+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] += TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            }
            else if (TC_specAOidx_[indx]>0) {
              k = TC_specAOidx_[indx]-1 ;
              iJac = TC_Nvjac_*(i+3)+3+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] -= TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            }
            indx++ ;
          }
        } /* done section for arbitrary order reaction */
      } /* done loop over species i */
    } /* done if real stoichiometric coeffs */

    if ( arbord ) iord++ ;

  } /* done loop over the number of reactions */
  
  /* get density, cpmix, species cp, species enthalpies */
  TC_getRhoMixMs  (scal,TC_Nvars_,&rhomix) ;
  TC_getCpMixMsP  (scal,TC_Nvars_,&cpmixP) ; 
  TC_getMs2CpMixMs(scal,TC_Nvars_,&cpmix ) ;
  TC_getHspecMs   (t1  ,TC_Nspec_,TC_hks ) ;

  /* Multiply by appropriate masses */
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
  {
    /* F[3+i,2] */
    iJac = TC_Nvjac_ * ( i + 3 ) + 2 ; 
    jac[iJac] *= TC_sMass_[i] / rhomix ;

    /* F[3+i,3+k] */
    iJac++ ;
    for ( k = 0 ; k < TC_Nspec_ ; k++, iJac++ )
      jac[iJac] *= TC_sMass_[i] ;
  }

  for ( k = 0 ; k < TC_Nspec_ ; k++ )
  {    
    double oMassk ;
    iJac = TC_Nvjac_*3 + 3 + k ;
    oMassk = 1.0 / TC_sMass_[k] ;
    /* F[3+i,3+k] */
    for ( i = 0 ; i < TC_Nspec_ ; i++, iJac += TC_Nvjac_ )
      jac[iJac] *= oMassk ;
  }

  /* transform T-derivatives to SI */
  for ( i=0, iJac=TC_Nvjac_*3+2 ; i < TC_Nspec_ ; i++, iJac+=TC_Nvjac_)
    jac[iJac] *= 1.e3 ;

  /* compute F[3+i,0] */
  iJac = TC_Nvjac_*3 ;
  for ( i=0, iJac=TC_Nvjac_*3 ; i < TC_Nspec_ ; i++, iJac+=TC_Nvjac_)
  {
    for ( k = 0 ; k<TC_Nspec_ ; k++ )
      jac[iJac] += Yspec[k]*jac[iJac+3+k] ;

    jac[iJac] = ( jac[iJac] - TC_omg[i] / rhomix ) / rhomix ;
  }

  /* compute F[2,0] */
  iJac = TC_Nvjac_*2 ;
  for ( k = 0 ; k < TC_Nspec_ ; k++ )
    jac[iJac] -= TC_hks[k]*jac[iJac+TC_Nvjac_*(k+1)] ;
  jac[iJac] /= cpmix ;


  /* compute F[2,2] */
  sum1 = 0.0; 
  sum2 = 0.0; 
  sum3 = 0.0 ;
  for ( k = 0, iJac=TC_Nvjac_*3+2 ; k < TC_Nspec_ ; k++,iJac+=TC_Nvjac_ )
  {
    sum1 += TC_hks [k]*TC_omg[k] ;
    sum2 += TC_cpks[k]*TC_omg[k] ;
    sum3 += TC_hks [k]*jac[iJac] ;
  }
  iJac = TC_Nvjac_*2+2 ;
  jac[iJac] = (sum1*cpmixP/cpmix-sum2)/(rhomix*cpmix)-sum3/cpmix ;

  /* compute F[2,3+i] */
  sum1 /= (rhomix*cpmix*cpmix) ;
  for ( i = 0,iJac=TC_Nvjac_*2+3 ; i < TC_Nspec_ ; i++,iJac++ )
  {
    for ( k = 0 ; k < Nspec ; k++)
      jac[iJac] -= TC_hks[k]*jac[iJac+TC_Nvjac_*(k+1)] ;

    jac[iJac] = jac[iJac]/cpmix+sum1*TC_cpks[i] ;

  }

  return ( ans ) ;

}

/**
 * \ingroup jacs
 * \brief Computes numerical Jacobian for the constant volume system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ based on T and Y's
*/
int TC_getJacCVTYNnum(double *scal, int Nspec, double *jac)
{

  int ans, i, j, iJac ;
  double perturb, *scalp, *srcCV, *srcCVp ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacCVTYNnum", Nspec, TC_Nspec_ ) ;

  ans = 0 ;

  /* allocate arrays */
  scalp  = (double *) malloc( (TC_Nspec_+1) * sizeof(double) ) ;
  srcCVp = (double *) malloc( (TC_Nspec_+1) * sizeof(double) ) ;
  srcCV  = (double *) malloc( (TC_Nspec_+1) * sizeof(double) ) ;

  /* get reference source terms */
  ans = TC_getSrcCV(scal,Nspec+1,srcCV) ; 

  for ( j=0; j<Nspec+1; j++ ) {

    for ( i=0; i<Nspec+1; i++ ) scalp[i] = scal[i] ;

    perturb = fabs(TC_reltol*scal[j]) ;
    if ( perturb == 0.0 ) perturb = TC_abstol ;
    scalp[j] += perturb ;
    ans = TC_getSrcCV(scalp,Nspec+1,srcCVp) ;

    for ( i = 0, iJac = j*(Nspec+1) ; i < Nspec+1 ; i++, iJac ++ )
      jac[iJac] = (srcCVp[i]-srcCV[i])/perturb ;

  }

  free( scalp  );
  free( srcCVp );
  free( srcCV  );

  return (ans);

}

/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian for the constant volume system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ based on T and Y's
*/
int TC_getJacCVTYNanl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i, j, k, kspec, ir, indx, indxR, indxT, iJac, itbdy, ipfal, iord, nvJac ;
  double *Yspec, t1, t_1, tln ;
  double wmix, rhomix,cpmix,cpmixP, gamma, gammaT, gammaYj, c1g, c2g, dc1g, dc2g ;
  double sum1, sum2, sum3, sumYoW, dtmp ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacCVTYNanl", Nspec, TC_Nspec_ ) ;

  ans = 0 ;

  /* re-scale pressure to get appropriate density */
  ans = TC_getRhoMixMs(scal,Nspec+1,&rhomix) ;
  TC_setThermoPres(TC_pressure_ * TC_rho_/rhomix) ;
  rhomix = TC_rho_;

  /* Local NvJac */
  nvJac = Nspec+1;

  /* scal array contains T and all species mass fractions */
  /* clear data */
  for ( i = 0 ; i<TC_Nspec_; i++) TC_omg[i] = 0.0 ;

  /* Compute molar concentrations (kmol/m3) and transform them to moles/cm3 */
  ans = TC_getMs2Cc(scal,TC_Nspec_+1,TC_Xconc) ;
  if ( ans != 0 ) return ( ans ) ;
  double sumXc = 0.0;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) {
    TC_Xconc[i] *= 1.e-3 ;
    sumXc += TC_Xconc[i];
  }

  /* set species array */
  Yspec = &scal[1] ;

  /* compute functions and powers of temperature */
  t1  = scal[0] ;
  t_1 = 1.0/t1  ;
  tln = log(t1) ;

  /* get 3rd-body concentrations */
  ans = TC_get3rdBdyConc(TC_Xconc,TC_Mconc) ;
  if ( ans != 0 ) return ( ans ) ;

  /* compute (-ln(T)+dS/R-dH/RT) for each species */
  ans = TC_getgk (t1,t_1,tln) ; if ( ans != 0 ) return ( ans ) ;
  ans = TC_getgkp(t1,t_1,tln) ; if ( ans != 0 ) return ( ans ) ;

  /* compute forward and reverse rate constants + derivatives */
  ans = TC_getkForRev    (t1,t_1,tln) ; if ( ans != 0 ) return ( ans ) ;
  ans = TC_getkForRevP   (t1,t_1    ) ; if ( ans != 0 ) return ( ans ) ;
  ans = TC_getkForRevPder(   t_1,tln) ; if ( ans != 0 ) return ( ans ) ;

  /* compute rate-of-progress */
  ans = TC_getRateofProg(TC_Xconc) ; if ( ans != 0 ) return ( ans ) ;

  /* compute pressure dependent factors */
  ans = TC_getCrnd( t1, t_1, tln, TC_Xconc, TC_Mconc ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reaction rates */
  for ( i = 0 ; i<TC_Nreac_; i++)
  {

    TC_rop[i] *= TC_Crnd[i] ;

    indx  = i*TC_maxSpecInReac_ ;
    for ( j = 0; j<TC_reacNreac_[i] ; j++)
    {
      int kspec = TC_reacSidx_[indx] ;
      TC_omg[kspec] += ((double) TC_reacNuki_[indx])*TC_rop[i] ;
      indx++ ;
    }

    indx  = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
    for ( j = 0; j<TC_reacNprod_[i] ; j++)
    {
      int kspec = TC_reacSidx_[indx] ;
      TC_omg[kspec] += ((double) TC_reacNuki_[indx])*TC_rop[i] ;
      indx++ ;
    }

  } /* done loop over the number of reactions */

  /* check for reactions with real stoichiometric coefficients */
  if ( TC_nRealNuReac_ > 0 )
    for ( ir = 0 ; ir < TC_nRealNuReac_; ir++ ) {

      i = TC_reacRnu_[ir] ;

      indx  = i *TC_maxSpecInReac_ ;
      indxR = ir*TC_maxSpecInReac_ ;      
      for ( j = 0 ; j<TC_reacNreac_[i] ; j++) {
        int kspec = TC_reacSidx_[indx] ;
        TC_omg[kspec] += TC_reacRealNuki_[indxR]*TC_rop[i] ;
        indx++  ;
        indxR++ ;
      }

      indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      //indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;  // Bugfix Cosmin 2016/12/28
      for ( j = 0; j<TC_reacNprod_[i] ; j++) {
        int kspec = TC_reacSidx_[indx] ;
        TC_omg[kspec] += TC_reacRealNuki_[indxR]*TC_rop[i] ;
        indx++  ;
        indxR++ ;
      }
    } /* done loop over the number of reactions */

  /* transform from mole/(cm3.s) to kg/(m3.s) */
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) TC_omg[i] *= (1.e3*TC_sMass_[i]) ;

  /* assemble jacobian terms, F[1+i,0...] */
  for ( i = 0 ; i < nvJac* nvJac ; i++) jac[i] = 0.0 ;
  itbdy = 0 ; 
  ipfal = 0 ; 
  iord  = 0 ;
  for (j = 0 ; j < TC_Nreac_; j++) {

    unsigned int arbord ;
    /* compute {\partial Crnd}/{\partial T,Y1,Y2,...,YN} */
    TC_getCrndDer( j, &itbdy, &ipfal, t1, t_1, tln, TC_Xconc, TC_Mconc ) ;

    arbord = 0 ;
    if ( iord < TC_nOrdReac_ )
      if (TC_reacAOrd_[iord] == j) arbord = 1 ;

    if ( TC_reacScoef_[j] == -1 ) {
      
      /* reaction has integer stoichiometric coefficients */
      for ( i = 0 ; i < TC_Nspec_ ; i++) {
      
        /* skip if species i is not involved in reaction j */
        int indxIJ = j*TC_Nspec_+i ;
        if ( TC_NuIJ[indxIJ] == 0 ) continue ;
      
        /* add 1st term for the T,Y1,Y2,...,YN - derivatives */
        iJac = nvJac*(1+i) ; 
        for ( k = -1 ; k < TC_Nspec_ ; k++,iJac++)
          jac[iJac] += TC_NuIJ[indxIJ]*TC_CrndDer[k+1]*(TC_ropFor[j]-TC_ropRev[j]) ;

        /* add 2nd term for T-derivative F[3+i,2] */
        iJac = nvJac*(1+i) ; 
        jac[iJac] +=
        TC_NuIJ[indxIJ]*TC_Crnd[j]*(TC_ropFor[j]*TC_kforP[j]-TC_ropRev[j]*TC_krevP[j]) ;
	
        /* add 3rd term for T-derivative F[3+i,2] */
        iJac = nvJac*(1+i) ; 
        jac[iJac] +=
        TC_NuIJ[indxIJ]*TC_Crnd[j]*(TC_ropFor[j]*TC_kforPder[j]-TC_ropRev[j]*TC_krevPder[j])
        * TC_Rcgs_*sumXc ;

        /* add 2nd term for species derivatives F[3+i,3+k] */
        if ( !arbord ) {
          indx = j*TC_maxSpecInReac_ ;
          for ( kspec = 0; kspec<TC_reacNreac_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = nvJac*(1+i)+1+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] += TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            jac[iJac] += TC_NuIJ[indxIJ]*TC_Crnd[j]*TC_ropFor[j]*TC_kforPder[j]*TC_Rcgs_*scal[0];
            indx++ ;
          } /* done loop over species k in reactants (2nd term) */

          indx  = j*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
          for ( kspec = 0; kspec<TC_reacNprod_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = nvJac*(1+i)+1+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] -= TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            jac[iJac] -= TC_NuIJ[indxIJ]*TC_Crnd[j]*TC_ropRev[j]*TC_krevPder[j]*TC_Rcgs_*scal[0];
            indx++ ;
          } /* done loop over species k in products (2nd term) */
        } /* done section for non-arbitrary order reactions */
        else {
          /* arbitrary order reaction */
          indx = iord*TC_maxOrdPar_ ;
          for ( kspec = 0; kspec<TC_maxOrdPar_ ; kspec++) {
            if (TC_specAOidx_[indx]<0) {
              k = -TC_specAOidx_[indx]-1 ;
              iJac = nvJac*(1+i)+1+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] += TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            }
            else if (TC_specAOidx_[indx]>0) {
              k = TC_specAOidx_[indx]-1 ;
              iJac = nvJac*(1+i)+1+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] -= TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            }
            indx++ ;
          }
        } /* done section for arbitrary order reaction */
      } /* done loop over species (counter i) for integer stoichiometric coefficients */
    } /* end if integer coeffs */
    else {
      /* reaction has real stoichiometric coefficients */
      for ( i = 0 ; i < TC_Nspec_ ; i++) {
        int jr = TC_reacScoef_[j] ;
        /* skip if species i is not involved in reaction j */
        int indxIJ = jr*TC_Nspec_+i ;
        if ( TC_RealNuIJ[indxIJ] < TCSMALL ) continue ;
      
        /* add 1st term for the T,Y1,Y2,...,YN - derivatives */
        iJac = nvJac*(1+i) ; 
        for ( k = -1 ; k < TC_Nspec_ ; k++,iJac++)
          jac[iJac] += TC_RealNuIJ[indxIJ]*TC_CrndDer[k+1]*(TC_ropFor[j]-TC_ropRev[j]) ;

        /* add 2nd term for T-derivative F[3+i,2] */
        iJac = nvJac*(1+i) ; 
        jac[iJac] +=
        TC_RealNuIJ[indxIJ]*TC_Crnd[j]*(TC_ropFor[j]*TC_kforP[j]-TC_ropRev[j]*TC_krevP[j]) ;

        /* add 3rd term for T-derivative F[3+i,2] */
        iJac = nvJac*(1+i) ; 
        jac[iJac] +=
        TC_RealNuIJ[indxIJ]*TC_Crnd[j]*(TC_ropFor[j]*TC_kforPder[j]-TC_ropRev[j]*TC_krevPder[j]) ;
      
        /* add 2nd term for species derivatives F[3+i,3+k] */
        if ( !arbord ) {

          indx = j *TC_maxSpecInReac_ ;
          for ( kspec = 0; kspec<TC_reacNreac_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = nvJac*(1+i)+1+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] += TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            indx++ ;
          } /* done loop over species k in reactants (2nd term) */

          indx = j *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
          for ( kspec = 0; kspec<TC_reacNprod_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = nvJac*(1+i)+1+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] -= TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            indx++ ;	  
          } /* done loop over species k in products (2nd term) */
        }
        else {
          /* arbitrary order reaction */
          indx = iord*TC_maxOrdPar_ ;
          for ( kspec = 0; kspec<TC_maxOrdPar_ ; kspec++) {
            if (TC_specAOidx_[indx]<0) {
              k = -TC_specAOidx_[indx]-1 ;
              iJac = nvJac*(1+i)+1+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] += TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            }
            else if (TC_specAOidx_[indx]>0) {
              k = TC_specAOidx_[indx]-1 ;
              iJac = nvJac*(1+i)+1+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] -= TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            }
            indx++ ;
          }
        } /* done section for arbitrary order reaction */
      } /* done loop over species i */
    } /* done if real stoichiometric coeffs */

    if ( arbord ) iord++ ;

  } /* done loop over the number of reactions */
  
  /* get density, cpmix, species cp, species enthalpies */
  //TC_getRhoMixMs  (scal,TC_Nvars_,&rhomix) ;
  TC_getCpMixMsP  (scal,TC_Nvars_,&cpmixP) ; 
  TC_getMs2CpMixMs(scal,TC_Nvars_,&cpmix ) ;
  TC_getHspecMs   (t1  ,TC_Nspec_,TC_hks ) ;

  /* wmix */
  sumYoW = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) sumYoW += Yspec[i] / TC_sMass_[i] ;
  wmix = 1.0 / sumYoW ;

  /* Multiply by appropriate masses */
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
  {
    /* F[1+i,0] */
    iJac = nvJac * ( 1+i ) ; 
    jac[iJac] *= TC_sMass_[i] / rhomix ;

    /* F[1+i,1+k] */
    iJac++ ;
    for ( k = 0 ; k < TC_Nspec_ ; k++, iJac++ )
      jac[iJac] *= TC_sMass_[i] ;
  }

  for ( k = 0 ; k < TC_Nspec_ ; k++ )
  {    
    double oMassk ;
    iJac = nvJac + 1 + k ;
    oMassk = 1.0 / TC_sMass_[k] ;
    /* F[1+i,1+k] */
    for ( i = 0 ; i < TC_Nspec_ ; i++, iJac += nvJac )
      jac[iJac] *= oMassk ;
  }

  /* transform T-derivatives to SI */
  for ( i=0, iJac=nvJac; i < TC_Nspec_ ; i++, iJac+=nvJac)
    jac[iJac] *= 1.e3 ;

  /* compute F[0,0] */
  sum1 = 0.0; 
  sum2 = 0.0; 
  sum3 = 0.0 ;
  for ( k = 0, iJac=nvJac ; k < TC_Nspec_ ; k++,iJac+=nvJac )
  {
    sum1 += TC_hks [k]*TC_omg[k] ;
    sum2 += TC_cpks[k]*TC_omg[k] ;
    sum3 += TC_hks [k]*jac[iJac] ;
  }
  iJac = 0 ;
  jac[iJac] = (sum1*cpmixP/cpmix-sum2)/(rhomix*cpmix)-sum3/cpmix ;

  /* compute F[0,1+i] */
  sum1 /= (rhomix*cpmix*cpmix) ;
  for ( i = 0, iJac=1 ; i < TC_Nspec_ ; i++,iJac++ )
  {
    for ( k = 0 ; k < Nspec ; k++)
      jac[iJac] -= TC_hks[k]*jac[iJac+nvJac*(k+1)] ;

    jac[iJac] = jac[iJac]/cpmix+sum1*TC_cpks[i] ;

  }

  /* additions for CV equations */
  gamma = cpmix / ( cpmix - TC_Runiv_/wmix );
  c1g = gamma ;
  c2g = gamma-1.0 ;
  dc1g = 1.0 ;
  dc2g = 1.0 ;

  sum1 = 0.0; 
  sum2 = 0.0; 
  sum3 = 0.0 ;
  for ( k = 0, iJac=nvJac ; k < TC_Nspec_ ; k++,iJac+=nvJac )
  {
    sum1 += TC_hks[k]*TC_omg[k] ;
    sum2 += TC_omg[k]/TC_sMass_[k] ;
    sum3 += jac[iJac]/TC_sMass_[k] ;
  }
  
  gammaT = -cpmixP * gamma * (gamma-1) / cpmix  ;
  jac[0] = dc1g * gammaT * (-sum1/(rhomix*cpmix))+c1g * jac[0]
          +(dc2g * gammaT * t1 + c2g) * wmix/rhomix * sum2 + c2g * t1 * wmix * sum3;
  
  for ( k = 0, iJac=1 ; k < TC_Nspec_ ; k++, iJac++ )
  {
    gammaYj = gamma * ( TC_cpks[k]/ cpmix - (TC_cpks[k]-TC_Runiv_/TC_sMass_[k])/(cpmix-TC_Runiv_/wmix) ) ;
    sum3 = 0.0 ;
    for ( i = 0 ; i < TC_Nspec_ ; i++ )
      sum3 += jac[iJac+(1+i)*nvJac]/TC_sMass_[i] ;
    jac[iJac] =  dc1g * gammaYj * (-sum1/(rhomix*cpmix)) + c1g * jac[iJac]
          +(dc2g * gammaYj - c2g * wmix / TC_sMass_[k]) * t1 * wmix/rhomix * sum2 + c2g * t1 * wmix * sum3;
  }

  /* transpose jacobian (temporary bugfix-2015/01/25) */
  for ( i = 0; i < TC_Nspec_+1; i++ )
    for ( j = i+1; j < TC_Nspec_+1; j++ ) {
      indx  = i*(TC_Nspec_+1)+j;
      indxT = j*(TC_Nspec_+1)+i;
      dtmp       = jac[indx];
      jac[indx]  = jac[indxT];
      jac[indxT] = dtmp ;
    }

  return ( ans ) ;

}

/*
               _      _            ____  ____ _______   ___   _                       
     __ _  ___| |_   | | __ _  ___|  _ \|  _ \_   _\ \ / / \ | |_ __  _   _ _ __ ___  
    / _` |/ _ \ __|  | |/ _` |/ __| |_) | |_) || |  \ V /|  \| | '_ \| | | | '_ ` _ \ 
   | (_| |  __/ || |_| | (_| | (__|  _ <|  __/ | |   | | | |\  | | | | |_| | | | | | |
    \__, |\___|\__\___/ \__,_|\___|_| \_\_|    |_|   |_| |_| \_|_| |_|\__,_|_| |_| |_|
    |___/                                                                             

*/
/**
 * \ingroup jacs
 * \brief Computes numerical Jacobian for the system 
 * \f$(\rho,P,T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
*/
int TC_getJacRPTYNnum(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i, j, k, iJac ;
  double *Yspec, rhomix, cpmix, cpmixP, perturb, sum1, sum2, sum3 ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacRPTYNnum", Nspec, TC_Nspec_ ) ;

  ans = 0 ;

  /* scal array contains T and all species mass fractions */
  for ( i = 0 ; i<TC_Nspec_+1; i++) TC_scalIn[i] = scal[i] ;
  Yspec = &scal[1] ;

  /* clear data */
  for ( i = 0 ; i < TC_Nspec_; i++ ) TC_omg[i] = TC_omgP[i] = 0.0 ;
  for ( i = 0 ; i < TC_Nvjac_*TC_Nvjac_ ; i++ ) jac[i] = 0.0 ;
  
  ans = TC_getRhoMixMs   ( TC_scalIn   , Nspec+1,&rhomix ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2CpMixMs ( TC_scalIn   , Nspec+1,&cpmix  ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getCpMixMsP( TC_scalIn   , Nspec+1,&cpmixP ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getHspecMs    ( TC_scalIn[0], Nspec  ,TC_hks  ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getCpSpecMs   ( TC_scalIn[0], Nspec  ,TC_cpks ) ;
  if ( ans != 0 ) return ( ans ) ;

  ans = TC_getMs2Cc ( TC_scalIn, Nspec+1, &TC_scalIn[1]) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getTXC2RRms( TC_scalIn, Nspec+1, TC_omg) ;
  if ( ans != 0 ) return ( ans ) ;

  /* perturb temperature for F[3+i,2] */
  /* perturb = fabs(TC_reltol*TC_scalIn[0]+TC_abstol) ; */
  perturb = fabs(TC_reltol*TC_scalIn[0]) ;
  if ( perturb == 0.0 ) perturb = TC_abstol ;
  TC_scalIn[0] += perturb ;
  TC_getTXC2RRms( TC_scalIn, Nspec+1, TC_omgP) ;

  for ( i = 0,iJac = TC_Nvjac_*3+2 ; i<Nspec ; i++,iJac+=TC_Nvjac_ )
    jac[iJac] = (TC_omgP[i]-TC_omg[i])/(rhomix*perturb) ;
 
  /* perturb molar concentrations for F[3+i,3+j] */
  for ( j = 0 ; j < Nspec ; j++ )
  {
    for ( i = 0 ; i < Nspec+1; i++) TC_scalIn[i] = scal[i] ;

    ans = TC_getMs2Cc ( TC_scalIn, Nspec+1, &TC_scalIn[1]) ;
    /* perturb = fabs(TC_reltol*TC_scalIn[j+1]+TC_abstol) ; */
    perturb = fabs(TC_reltol*TC_scalIn[j+1]) ;
    if ( perturb == 0.0 ) perturb = TC_abstol ;
    TC_scalIn[j+1] += perturb ;
    ans = TC_getTXC2RRms( TC_scalIn, Nspec+1, TC_omgP) ;

    for ( i = 0,iJac = TC_Nvjac_*(i+3)+3+j ; i < Nspec ; i++,iJac+=TC_Nvjac_ )
      jac[iJac] = (TC_omgP[i]-TC_omg[i])/(TC_sMass_[j]*perturb) ;
    
  }

  /* compute F[3+i,0] */
  for ( i=0, iJac=TC_Nvjac_*3 ; i < TC_Nspec_ ; i++, iJac+=TC_Nvjac_)
  {
    for ( k = 0 ; k<TC_Nspec_ ; k++ )
      jac[iJac] += Yspec[k]*jac[iJac+3+k] ;

    jac[iJac] = ( jac[iJac] - TC_omg[i] / rhomix ) / rhomix ;
  }

  /* compute F[2,0] */
  iJac = TC_Nvjac_*2 ;
  for ( k = 0 ; k < TC_Nspec_ ; k++ )
    jac[iJac] -= TC_hks[k]*jac[iJac+TC_Nvjac_*(k+1)] ;
  jac[iJac] /= cpmix ;


  /* compute F[2,2] */
  sum1 = 0.0 ;
  sum2 = 0.0 ;
  sum3 = 0.0 ;
  for ( k = 0, iJac=TC_Nvjac_*3+2 ; k < TC_Nspec_ ; k++,iJac+=TC_Nvjac_ )
  {
    sum1 += TC_hks [k]*TC_omg[k] ;
    sum2 += TC_cpks[k]*TC_omg[k] ;
    sum3 += TC_hks [k]*jac[iJac] ;
  }
  iJac = TC_Nvjac_*2+2 ;
  jac[iJac] = (sum1*cpmixP/cpmix-sum2)/(rhomix*cpmix)-sum3/cpmix ;

  /* compute F[2,3+i] */
  sum1 /= (rhomix*cpmix*cpmix) ;
  for ( i = 0,iJac=TC_Nvjac_*2+3 ; i < TC_Nspec_ ; i++,iJac++ )
  {
    for ( k = 0 ; k < Nspec ; k++)
      jac[iJac] -= TC_hks[k]*jac[iJac+TC_Nvjac_*(k+1)] ;

    jac[iJac] = jac[iJac]/cpmix+sum1*TC_cpks[i] ;

  }

  return ( ans ) ;

}
