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

#include "ProblemCoreFlows.hxx"
#include "Node.hxx"
#include "utilitaire_algebre.h"

class WaveStaggered : public ProblemCoreFlows{
public :
	/** \fn WaveStaggered
	 * \brief Constructor for the Navier-Stokes system
	 * \param [in] phaseType : \ref Liquid or \ref Gas
	 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
	 * \param [in] int : mesh dimension
	 * \param [in] bool : There are two possible equations of state for the fluid
	 *  */
	WaveStaggered(int dim, double kappa, double rho, MPI_Comm comm = MPI_COMM_WORLD);

	void setInitialField(const Field &field);
	void setInitialFieldFunction(const Mesh& M, std::map<int, double> V, EntityType typeField, const string name); //TODO : à dégager

	//! system initialisation
	void initialize();

	/** \fn terminate
     * \brief empties the memory
     * @param void
     *  */
    void terminate();

	double getTimeStep();

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

    /** \fn saveVelocity
     * \brief saves the velocity field in a separate 3D file so that paraview can display the streamlines
     * @param bool
     * */
    void saveVelocity(bool save_v=true){
        _saveVelocity=save_v;
    }

	 /** \fn testConservation
     * \brief Teste et affiche la conservation de masse et de la quantité de mouvement
     * \Details la fonction est virtuelle pure, on la surcharge dans chacun des modèles
     * @param void
     * */
	void testConservation();

	std::map<int,double>  getboundaryPressure();
	void  setboundaryPressure(map< int, double> BoundaryPressure);
	void  setboundaryVelocity(map< int, double> BoundaryVelocity);

	void  abortTimeStep();
	bool  initTimeStep( double dt);
	vector<string> getInputFieldsNames();
	void setInputField(const string& nameField, Field& inputField );


protected :
	Field _Velocity, _Pressure ;
	Vec _newtonVariation, _primitiveVars;
	Mat _InvVol; // matrice Q such that U^n+1 = (Id + dt V^-1 _A)U^n for explicit scheme
	double _kappa, _rho,  _c, _d, _maxPerim, _minCell ;
	bool _savePressure, _saveVelocity;
	std::map<int, double>  _boundaryPressure;
	bool _facesBoundinit; // To ensure that the boundary velocity is initialized after the initial velocity 
				

};
#endif /* WAVESTAGGERED_HXX_*/
