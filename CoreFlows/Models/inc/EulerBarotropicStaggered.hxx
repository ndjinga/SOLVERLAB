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
#include "Node.hxx"
#include "utilitaire_algebre.h"
#include "Mesh.hxx"

//! enumeration phaseType
/*! The material phase can be Gas or liquid  */
enum phaseType
{
    Liquid,/**< Material considered is Liquid */
    Gas/**< Material considered is Gas */
};

class EulerBarotropicStaggered : public WaveStaggered{
public :
	/** \fn EulerBarotropicStaggered
	 * \param [in] phaseType : \ref Liquid or \ref Gas
	 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
	 * \param [in] int : mesh dimension
	 *  */
	EulerBarotropicStaggered(phaseType fluid, pressureEstimate pEstimate, int dim);


	//! system initialisation
	void initialize();

	/** \fn terminate
     * \brief empties the memory
     * @param void
     *  */
    void terminate();


	/** \fn computeTimeStep
     * \brief Proposes a value for the next time step to be solved using mesh data and cfl coefficient
     *  \return  double dt the proposed time step
     *  \return  bool stop, true if the calculation should not be continued (stationary state, maximum time or time step numer reached)
     *  */
    double computeTimeStep(bool & stop);

	/** \fn computeNewtonVariation
	 * \brief Builds and solves the linear system to obtain the variation Vkp1-Vk in a Newton scheme using primitive variables
	 * @param
	 * */
	void computeNewtonVariation();


	/** \fn getStiffenedGasEOS
     * \brief return the stiffened gas law associated to fluid i
     * @param int i : the index of the fluid
     * @return throws an exception if the fluid with index i does not follow a stiffened gas law.
     * */
    StiffenedGas getStiffenedGasEOS(int i)
    {
        StiffenedGas * result = dynamic_cast<StiffenedGas*>(_fluides[i]); 
        if(result)
             return *result;
        else
            throw CdmathException("ProblemFluid::getStiffenedGasEOS() : fluid EOS is not a stiffened gas law");
    }


protected :
 /** Fluid equation of state **/
    vector<    Fluide* > _fluides;//
	CompressibleFluid *_compressibleFluid;
	double _Tref; //EOS reference temperature
    double _Pref; //EOS reference pressure

	Mat _Conv, _DivRhoU, _LaplacianVelocity  ;
	double _c;
	std::vector<double> _Entropy;
				

};
#endif /* EULERBAROTROPICSTAGGERED_HXX_*/
