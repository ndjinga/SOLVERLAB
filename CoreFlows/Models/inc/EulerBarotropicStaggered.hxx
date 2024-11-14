//============================================================================
/**
 * \file EulerBarotropicStaggered.hxx
 * \author Esteban COIFFIER
 * \version 1.0
 * \date 01 janv. 2018
 * \briefThe Wave system
 * */
//============================================================================

/*! \class EulerBarotropicStaggered EulerBarotropicStaggered.hxx "EulerBarotropicStaggered.hxx"
 *  \brief The Wave system
 *  \details Wave system
 */
#ifndef EULERBAROTROPICSTAGGERED_HXX_
#define EULERBAROTROPICSTAGGERED_HXX_

#include "WaveStaggered.hxx"
#include "StiffenedGas.hxx"
#include "BarotropicLaw.hxx"
#include "Node.hxx"
#include "utilitaire_algebre.h"
#include "Mesh.hxx"

//! enumeration phaseType
/*! The material phase can be Gas or liquid  */
enum phaseTypeStaggered
{
    LiquidStaggered,
    GasStaggered
};

class EulerBarotropicStaggered : public WaveStaggered{
public :
	/** \fn EulerBarotropicStaggered
	 * \param [in] phaseType : \ref Liquid or \ref Gas
	 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
	 * \param [in] int : mesh dimension
	 *  */
	EulerBarotropicStaggered(phaseTypeStaggered fluid, pressureEstimate pEstimate, double a, double gamma, int dim);


	//! system initialisation
	void initialize();

	/** \fn terminate
     * \brief empties the memory
     * @param void
     *  */
    void terminate();
    double computeTimeStep(bool & stop);

    BarotropicLaw getBarotropicEOS(int i){
        BarotropicLaw * result = dynamic_cast<BarotropicLaw*>(_fluides[i]); 
        if(result)
             return *result;
        else
            throw CdmathException("Fluid EOS is not barotropic");
    }

    void AssembleMetricsMatrices();
    void UpdateDualDensity();
    bool iterateTimeStep(bool &converged);
    void testConservation();
     std::vector<double>  getTimeEvol();
    
protected :
 /** Fluid equation of state **/
    vector<    Fluide* > _fluides;//
	BarotropicLaw   *_compressibleFluid;
	double _Tref, _Pref; //EOS reference temperature &pressure

    PetscReal _rhoMax, _uMax;
    Vec _DualDensity ;
	Mat _InvVolPrim, _InvVolDual,_Conv, _DivRhoU, _LaplacianVelocity, _InvDualDensity  ;
	double _c, _ConvectiveMax;
	std::vector<double> _Entropy, _Time;
				

};
#endif /* EULERBAROTROPICSTAGGERED_HXX_*/
