/*! \file TC_src.c
 *    \brief Source term and Jacobian functions
 */

/*
      _____              _
     /  ___|            | |
     \ `--. _ __ ___    | |_ ___ _ __ _ __ ___  ___
      `--. \ '__/ __|   | __/ _ \ '__| '_ ` _ \/ __|
     /\__/ / | | (__ _  | ||  __/ |  | | | | | \__ \
     \____/|_|  \___(_)  \__\___|_|  |_| |_| |_|___/

*/
/*
                          _   ____
                __ _  ___| |_/ ___| _ __ ___
               / _` |/ _ \ __\___ \| '__/ __|
              | (_| |  __/ |_ ___) | | | (__
               \__, |\___|\__|____/|_|  \___|
               |___/

*/
/**
 * \ingroup srcs
 * \brief Returns source term for
 * \f[\frac{\partial T}{\partial t}=\omega_0,\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 *  based on temperature T and species mass fractions Y's
 */
//#define DEBUGMSG

int TC_getSrc(double *scal,int Nvars,double *omega)
{

/**
   \param scal : array of Nspec+1 doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
	        temperature T [K], mass fractions Y []
   \param Nvars : no. of variables \f$ N_{vars}= N_{spec}+1\f$
   \return omega : array of N<sub>spec</sub>+1 source terms for temperature and species
                   mass fractions: omega[0] : [(K/s)], omega[1...N<sub>spec</sub>] : [(1/s)]
*/

  int i, ans ;
  double temperature, *Yspec, *omegaspec, rhomix, orho, cpmix ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getSrc", Nvars, TC_Nvars_ ) ;

  ans         = 0 ;
  temperature =  scal[0] ;
  Yspec       = &scal[1] ;
  omegaspec   = &omega[1] ;

  /* get species molar reaction rates */
  ans = TC_getTY2RRms(scal,Nvars,omegaspec) ;
  if ( ans != 0 ) return ( ans ) ;

#ifdef DEBUGMSG
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
    printf("TC_getSrc() : Yspec[%-3d] = %e, omega[%-3d] = %e\n",i+1,Yspec[i],i+1,omega[i+1]) ;
#endif

  /* get density, cpmix */
  ans = TC_getRhoMixMs(scal,Nvars,&rhomix) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2CpMixMs(scal,Nvars,&cpmix ) ;
  if ( ans != 0 ) return ( ans ) ;

#ifdef DEBUGMSG
  printf("TC_getSrc() : (rhomix,cpmix)=(%e,%e)\n",rhomix,cpmix) ;
#endif

  /* get species enthalpies */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;
  if ( ans != 0 ) return ( ans ) ;

  /* transform reaction rate to source term (*Wi/rho) */
  orho = 1.0/rhomix ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omegaspec[i] *= orho ;

#ifdef NONNEG
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) 
    if ( ( Yspec[i]<0 ) && ( omegaspec[i] < 0 ) ) omegaspec[i] = 0.0 ;
#endif

  /* compute source term for temperature */
  omega[0] = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[0] -= omegaspec[i]*TC_hks[i] ;
  omega[0] /= cpmix ;

#ifdef DEBUGMSG
  printf("TC_getSrc() : omega[(%-3d] = %e\n",0,omega[0]) ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) printf("TC_getSrc() : omega[(%-3d] = %e\n",i+1,omega[i+1]) ;
  exit(1) ;
#endif

  return ( ans ) ;

}

/**
 * \ingroup srcs
 * \brief Returns \f$S\f$ term for
 * \f[\frac{\partial T}{\partial t}=\omega_0,\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 *  based on temperature T and species mass fractions Y's
 */

int TC_getSmat(double *scal,int Nvars, double *Smat)
{

/**
   \param scal : array of Nspec+1 doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
	        temperature T [K], mass fractions Y []
   \param Nvars : no. of variables \f$ N_{vars}= N_{spec}+1\f$
   \return Smat : array of \f$(N_{spec}+1)\times 2N_{reac}\f$ holding
                  the S components in column major format
*/

  int i, j, kspec, indx, indxR, ans ;
  double temperature, rhomix, orho, orhocp, cpmix ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getSrc", Nvars, TC_Nvars_ ) ;

  ans         = 0 ;
  temperature = scal[0] ;

  /* clean Smat */
  for ( i = 0 ; i < (TC_Nspec_+1)*2*TC_Nreac_ ; i++ ) Smat[i] = 0.0 ;

  /* get density, cpmix */
  ans = TC_getRhoMixMs  (scal,Nvars,&rhomix) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2CpMixMs(scal,Nvars,&cpmix ) ;
  if ( ans != 0 ) return ( ans ) ;
  orho   = 1.0/rhomix ;
  orhocp = orho/cpmix ;

  /* get species enthalpies and multiply by the weights */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;
  if ( ans != 0 ) return ( ans ) ;
  for ( i = 0 ; i<TC_Nspec_; i++) TC_hks[i] *= TC_sMass_[i] ;

  /* assemble matrix based on integer stoichiometric coefficients */
  for ( i = 0 ; i<TC_Nreac_; i++)
  {
    for ( j = 0; j<TC_reacNreac_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+j ;
      kspec = TC_reacSidx_[indx] ;
      Smat[i*(TC_Nspec_+1)]         -= TC_reacNukiDbl_[indx]*TC_hks[kspec] ;
      Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacNukiDbl_[indx]*TC_sMass_[kspec]*orho ;
    }

    for ( j = 0; j<TC_reacNprod_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j ;
      kspec = TC_reacSidx_[indx] ;
      Smat[i*(TC_Nspec_+1)]         -= TC_reacNukiDbl_[indx]*TC_hks[kspec] ;
      Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacNukiDbl_[indx]*TC_sMass_[kspec]*orho ;
    }

  } /* done loop over reactions */

  /* add contributions from real stoichiometric coefficients if any */
  if ( TC_nRealNuReac_ > 0 )
  {
    int ir ;

    for ( ir = 0 ; ir<TC_nRealNuReac_; ir++)
    {

      int i = TC_reacRnu_[ir] ;

      indx  = i *TC_maxSpecInReac_ ;
      indxR = ir*TC_maxSpecInReac_ ;

      for ( j = 0 ; j<TC_reacNreac_[i] ; j++)
      {
	kspec = TC_reacSidx_[indx] ;
        Smat[i*(TC_Nspec_+1)]         -= TC_reacRealNuki_[indxR]*TC_hks[kspec] ;
        Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacRealNuki_[indxR]*TC_sMass_[kspec]*orho ;
	indx++  ;
	indxR++ ;
      }

      indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++)
      {
	kspec = TC_reacSidx_[indx] ;
        Smat[i*(TC_Nspec_+1)        ] -= TC_reacRealNuki_[indxR]*TC_hks[kspec] ;
        Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacRealNuki_[indxR]*TC_sMass_[kspec]*orho ;
	indx++  ;
	indxR++ ;
      }

    } /* done loop over the number of reactions */

  } /* done if TC_nRealNuReac_ > 0 */

  for ( i = 0 ; i< TC_Nreac_; i++)
    Smat[i*(TC_Nspec_+1)] *= orhocp;

  for ( i = 0 ; i< TC_Nreac_*(TC_Nspec_+1); i++)
    Smat[i+TC_Nreac_*(TC_Nspec_+1)] = -Smat[i];

  return ( ans ) ;

}

/**
 * \ingroup srcs
 * \brief Returns source term for constant volume ignition system
 * \f[\frac{\partial T}{\partial t}=\omega_0,\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 *  based on temperature T and species mass fractions Y's
 */

int TC_getSrcCV(double *scal,int Nvars,double *omega)
{

/**
   \param scal : array of Nspec+1 doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
	        temperature T [K], mass fractions Y []
   \param Nvars : no. of variables \f$ N_{vars}= N_{spec}+1\f$
   \return omega : array of N<sub>spec</sub>+1 source terms for temperature and species
                   mass fractions: omega[0] : [(K/s)], omega[1...N<sub>spec</sub>] : [(1/s)]
*/

  int i, ans ;
  double temperature, *Yspec, *omegaspec, rhomix, orho, cpmix, wmix, sumYoW, gamma ;
  double psum ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getSrcCV", Nvars, TC_Nvars_ ) ;

  if ( TC_rhoset_ == 0 ) {
    printf("TC_getSrcCV() - density was not set -> Abort !\n") ;
    exit(1);
  }

  /* re-scale pressure to get appropriate density */
  ans = TC_getRhoMixMs(scal,Nvars,&rhomix) ;
  TC_setThermoPres(TC_pressure_ * TC_rho_/rhomix) ;
  rhomix = TC_rho_; // bugfix based on MV observation 2015/01/19

  // Checking equivalence in formulations
  /* wmix = 0.0 ; */
  /* for ( i = 0 ; i < TC_Nspec_ ; i++ ) wmix += scal[i+1]/TC_sMass_[i] ; */
  /* TC_pressure_ = TC_rho_*TC_Runiv_*scal[0]*wmix; */
  /* TC_prescgs_  = TC_pressure_*10.0 ; */
  /* rhomix = TC_rho_; */

  ans         = 0 ;
  temperature =  scal[0]  ;
  Yspec       = &scal[1]  ;
  omegaspec   = &omega[1] ;

  /* get species molar reaction rates */
  ans = TC_getTY2RRms(scal,Nvars,omegaspec) ;
  if ( ans != 0 ) return ( ans ) ;

#ifdef DEBUGMSG
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
    printf("TC_getSrc() : Yspec[%-3d] = %e, omega[%-3d] = %e\n",i+1,Yspec[i],i+1,omega[i+1]) ;
#endif

  /* get cpmix, gamma, and gamma factors */
  ans = TC_getMs2CpMixMs(scal,Nvars,&cpmix ) ;
  if ( ans != 0 ) return ( ans ) ;
  sumYoW = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) sumYoW += Yspec[i] / TC_sMass_[i] ;
  wmix = 1.0 / sumYoW ;

  /* gamma */
  gamma = cpmix / ( cpmix - TC_Runiv_/wmix );

#ifdef DEBUGMSG
  printf("TC_getSrc() : (cpmix,gamma)=(%e,%e)\n",cpmix,gamma) ;
#endif

  /* get species enthalpies */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;
  if ( ans != 0 ) return ( ans ) ;

  /* transform reaction rate to source term (*Wi/rho) */
  orho = 1.0/rhomix ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omegaspec[i] *= orho ;

#ifdef NONNEG
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) 
    if ( ( Yspec[i]<0 ) && ( omegaspec[i] < 0 ) ) omegaspec[i] = 0.0 ;
#endif

  /* compute source term for temperature */
  omega[0] = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[0] += omegaspec[i]*TC_hks[i] ;
  omega[0] = gamma * ( -omega[0] / cpmix ) ;

  /* add second term */
  psum = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) psum += omegaspec[i]/TC_sMass_[i] ;
  omega[0] += (gamma-1) * temperature * wmix * psum ;


#ifdef DEBUGMSG
  printf("TC_getSrc() : omega[(%-3d] = %e\n",0,omega[0]) ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) printf("TC_getSrcCV() : omega[(%-3d] = %e\n",i+1,omega[i+1]) ;
  exit(1) ;
#endif

  return ( ans ) ;

}

/**
 * \ingroup srcs
 * \brief Returns \f$S\f$ matrix for constant volume ignition system
 * \f[\frac{\partial T}{\partial t}=\omega_0,\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 *  based on temperature T and species mass fractions Y's
 */

int TC_getSmatCV(double *scal,int Nvars,double *Smat)
{

/**
   \param scal : array of Nspec+1 doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
	        temperature T [K], mass fractions Y []
   \param Nvars : no. of variables \f$ N_{vars}= N_{spec}+1\f$
   \return Smat : array of \f$(N_{spec}+1)\times 2N_{reac}\f$ holding
                  the S components in column major format
*/

  int i, j, kspec, indx, indxR, ans ;
  double temperature, *Yspec, rhomix, orho, cpmix, wmix, sumYoW, gamma, c1g, c2g ;
  double psum ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getSmatCV", Nvars, TC_Nvars_ ) ;

  if ( TC_rhoset_ == 0 ) {
    printf("TC_getSmatCV() - density was not set -> Abort !\n") ;
    exit(1);
  }

  /* clean Smat */
  for ( i = 0 ; i < (TC_Nspec_+1)*2*TC_Nreac_ ; i++ ) Smat[i] = 0.0 ;

  /* re-scale pressure to get appropriate density */
  ans = TC_getRhoMixMs(scal,Nvars,&rhomix) ;
  TC_setThermoPres(TC_pressure_ * TC_rho_/rhomix) ;

  ans         = 0 ;
  temperature =  scal[0]  ;
  Yspec       = &scal[1]  ;

  /* get cpmix, gamma, and gamma factors */
  ans = TC_getMs2CpMixMs(scal,Nvars,&cpmix ) ;
  if ( ans != 0 ) return ( ans ) ;
  sumYoW = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) sumYoW += Yspec[i] / TC_sMass_[i] ;
  wmix = 1.0 / sumYoW ;

  orho   = 1.0/TC_rho_ ;

  /* gamma */
  gamma = cpmix / ( cpmix - TC_Runiv_/wmix );
  c1g = gamma*orho/cpmix ;
  c2g = (gamma-1.0)*temperature*wmix*orho ;

  /* get species enthalpies and multiply by the weights */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;
  if ( ans != 0 ) return ( ans ) ;
  for ( i = 0 ; i<TC_Nspec_; i++) TC_hks[i] *= TC_sMass_[i] ;

  /* assemble matrix based on integer stoichiometric coefficients */
  for ( i = 0 ; i<TC_Nreac_; i++)
  {
    for ( j = 0; j<TC_reacNreac_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+j ;
      kspec = TC_reacSidx_[indx] ;
      Smat[i*(TC_Nspec_+1)]         += TC_reacNukiDbl_[indx]*(-TC_hks[kspec]*c1g+c2g) ;
      Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacNukiDbl_[indx]*TC_sMass_[kspec]*orho ;
    }

    for ( j = 0; j<TC_reacNprod_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j ;
      kspec = TC_reacSidx_[indx] ;
      Smat[i*(TC_Nspec_+1)]         += TC_reacNukiDbl_[indx]*(-TC_hks[kspec]*c1g+c2g) ;
      Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacNukiDbl_[indx]*TC_sMass_[kspec]*orho ;
    }

  } /* done loop over reactions */

  /* add contributions from real stoichiometric coefficients if any */
  if ( TC_nRealNuReac_ > 0 )
  {
    int ir ;

    for ( ir = 0 ; ir<TC_nRealNuReac_; ir++)
    {

      int i = TC_reacRnu_[ir] ;

      indx  = i *TC_maxSpecInReac_ ;
      indxR = ir*TC_maxSpecInReac_ ;

      for ( j = 0 ; j<TC_reacNreac_[i] ; j++)
      {
	kspec = TC_reacSidx_[indx] ;
        Smat[i*(TC_Nspec_+1)]         += TC_reacRealNuki_[indxR]*(-TC_hks[kspec]*c1g+c2g) ;
        Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacRealNuki_[indxR]*TC_sMass_[kspec]*orho ;
	indx++  ;
	indxR++ ;
      }

      indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++)
      {
	kspec = TC_reacSidx_[indx] ;
        Smat[i*(TC_Nspec_+1)        ] += TC_reacRealNuki_[indxR]*(-TC_hks[kspec]*c1g+c2g) ;
        Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacRealNuki_[indxR]*TC_sMass_[kspec]*orho ;
	indx++  ;
	indxR++ ;
      }

    } /* done loop over the number of reactions */

  } /* done if TC_nRealNuReac_ > 0 */

  for ( i = 0 ; i< TC_Nreac_*(TC_Nspec_+1); i++)
    Smat[i+TC_Nreac_*(TC_Nspec_+1)] = -Smat[i];

  return ( ans ) ;

}

/**
 * \ingroup srcs
 * \brief Returns source term for
 * \f[\frac{\partial \rho}{\partial t}=\omega_0,\rho\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 * based on \f$\rho\f$ and Y's
 */
int TC_getSrcCons(double *scal,int Nvars,double *omega)
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
  double temperature, rhomix, *Yspec, *omegaspec, cpmix, Wmix, sumOmega ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getSrcCons", Nvars, TC_Nvars_ ) ;

  ans       = 0 ;
  rhomix    =  scal[0] ;
  Yspec     = &scal[1] ;
  omegaspec = &omega[1] ;

  /* get temperature, cpmix, and mixture molecular weight */
  ans = TC_getTmixMs(scal,Nvars,&temperature) ;
  if ( ans != 0 ) return ( ans ) ;
  scal[0] = temperature ;

  /* get species molar reaction rates [kmol/(m3.s)]*/
  ans = TC_getTY2RRml(scal,Nvars,omegaspec) ;
  if ( ans != 0 ) return ( ans ) ;

  /* get cpmix [J/(kg.K)] and mixture molecular weight [kg/kmol] */
  ans = TC_getMs2CpMixMs(scal,Nvars,&cpmix) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix(Yspec,TC_Nspec_,&Wmix) ;
  if ( ans != 0 ) return ( ans ) ;

  /* get species enthalpies [J/kg] */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;
  if ( ans != 0 ) return ( ans ) ;

#ifdef NONNEG
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) 
    if ( ( Yspec[i]<0 ) && ( omegaspec[i] < 0 ) ) 
      omegaspec[i] = 0.0 ;
#endif

  /* sum molar reaction rates (for density source term) */
  sumOmega = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) sumOmega += omegaspec[i] ;

  /* transform molar reaction rate to mass reaction rate(*Wi) [kg/(m3.s)]*/
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omegaspec[i] *= TC_sMass_[i] ;

  /* compute source term for density */
  omega[0] = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[0] += omegaspec[i]*TC_hks[i] ;
  omega[0] /= (cpmix*temperature) ;
  omega[0] -= Wmix*sumOmega ;

  return ( ans ) ;

}
