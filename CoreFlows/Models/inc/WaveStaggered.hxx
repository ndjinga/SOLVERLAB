//============================================================================
/**
 * \file WaveStaggered.hxx
 * \author Esteban COIFFIER
 * \version 1.0
 * \date 01 janv. 2018
 * \briefThe Wave system
 * */
//============================================================================

/*! \class WaveStaggered WaveStaggered.hxx "WaveStaggered.hxx"
 *  \brief The Wave system
 *  \details Wave system
 */
#ifndef WAVESTAGGERED_HXX_
#define WAVESTAGGERED_HXX_

#include "ProblemFluid.hxx"
#include "Node.hxx"

class WaveStaggered : public ProblemFluid{
public :
	/** \fn WaveStaggered
	 * \brief Constructor for the Navier-Stokes system
	 * \param [in] phaseType : \ref Liquid or \ref Gas
	 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
	 * \param [in] int : mesh dimension
	 * \param [in] bool : There are two possible equations of state for the fluid
	 *  */
	WaveStaggered(phaseType fluid,int dim, double kappa, double rho);

	void setInitialField(const Field &field);

	//! system initialisation
	void initialize();

	/** \fn terminate
     * \brief empties the memory
     * @param void
     *  */
    void terminate();

	void save();

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

	/** \fn iterateTimeStep
	 * \brief calls computeNewtonVariation to perform one Newton iteration and tests the convergence
	 * @param
	 * @return boolean ok is true is the newton iteration gave a physically acceptable result
	 * */
	bool iterateTimeStep(bool &ok);

	 /** \fn validateTimeStep
     * \brief Validates the solution computed y solveTimeStep
     * \details updates the currens time t=t+dt, save unknown fields, resets the time step dt to 0, tests the stationnarity.
     * c It is a pure virtual function overloaded in each model
     * @param  void
     *  */
    void validateTimeStep();

	 /** \fn savePressure
     * \brief saves the Pressure field in a separate file 
     * @param bool
     * */
    void savePressure(bool save_p=true){
        _savePressure=save_p;
    }

protected :
	Field _Velocity, _Pressure ;
	Mat _Q; // matrice Q such that U^n+1 = (Id + dt V^-1 Q)U^n for explicit scheme
	double _kappa, _rho,  _c, _d;
	bool _savePressure;
	//Vec _boundaryPressure;
	std::map<int, double>  _boundaryPressure;
				

};
#endif /* WAVESTAGGERED_HXX_*/
