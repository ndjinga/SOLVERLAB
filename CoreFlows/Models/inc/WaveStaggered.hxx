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
#include "Mesh.hxx"

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

	//! system initialisation
	virtual void initialize();

	/** \fn terminate
     * \brief empties the memory
     * @param void
     *  */
    virtual void terminate();



	void save();

	/** \fn computeTimeStep
     * \brief Proposes a value for the next time step to be solved using mesh data and cfl coefficient
     *  \return  double dt the proposed time step
     *  \return  bool stop, true if the calculation should not be continued (stationary state, maximum time or time step numer reached)
     *  */
    virtual double computeTimeStep(bool & stop);

	/** \fn computeNewtonVariation
	 * \brief Builds and solves the linear system to obtain the variation Vkp1-Vk in a Newton scheme using primitive variables
	 * @param
	 * */
	void computeNewtonVariation(); //TODO ok ?

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

	double getTimeStep();
	void  abortTimeStep();
	bool  initTimeStep( double dt);

		vector<string> getInputFieldsNames();
	void setInputField(const string& nameField, Field& inputField );

    void savePressure(bool save_p=true){
        _savePressure=save_p;
	}
    void saveVelocity(bool save_v=true){
        _saveVelocity=save_v;
    }


	void ComputeEnergyAtTimeT();
	/*******Periodicity related ********/
	void setPeriodicFaces(	Mesh &M, const char &Direction, int ncells);
	std::map<int,int>  getFacePeriodicMap() const;
	bool  IsFaceBoundaryNotComputedInPeriodic(int j );
	bool  IsFaceBoundaryComputedInPeriodic(int j );

	/*******Imposed boundary conditions related ********/
	void setWallBoundIndex(int j );
	void setSteggerBoundIndex(int j ); //Imposed pressure and velocity
	void setInteriorIndex(int j ); //To avoid complicated implementation in periodic

	std::map<int,double>  getboundaryPressure() const;
	void  setboundaryPressure(map< int, double> BoundaryPressure);
	void  setboundaryVelocity(map< int, double> BoundaryVelocity);
	
	/***********Orientation *************/
	double getOrientation(int l, Cell Cint);
	void  setOrientation(int j,std::vector<double> vec_normal_sigma);

	/***********Post Pro *************/
	std::vector<double> ErrorL2VelocityInfty(const Field &ExactVelocityInftyAtFaces, const Field &ExactVelocityInftyAtCells );
	void InterpolateFromFacesToCells(const Field &atFaces, Field &atCells);

	//TODO à supprimer ?
	void setExactVelocityFieldAtCells(const Field &atCells);
	void ComputeMinCellMaxPerim();
	void AssembleLocalMetricMatricsInterior(int j, Cell Ctemp1 , Cell Ctemp2);


protected :
	Field _Velocity, _Pressure, _Velocity_at_Cells, _DivVelocity, _ExactVelocityInftyAtCells, _ExactVelocityInftyInterpolate;
	Vec _newtonVariation, _primitiveVars,  _BoundaryTerms, _primitiveVars_seq ;
	Mat _InvVol,_InvSurface, _Div, _LaplacianPressure, _DivTranspose,  _GradDivTilde ; // matrice Q such that U^n+1 = (Id + dt V^-1 _A)U^n for explicit scheme
	double _kappa, _rho,  _c, _d, _maxPerim, _minCell ;
	double *_vec_normal;
	PetscScalar _pExt, _pInt;
	bool _savePressure, _saveVelocity;
	std::map<int, double>  _boundaryPressure;
	std::map<int, std::vector<double> > _vec_sigma; // arbitrary degree of liberty associated to a face
	std::map<int,int> _FacePeriodicMap;
	std::vector<int>_WallBoundFaceSet, _SteggerBoundFaceSet, _InteriorFaceSet; // map of perdiodic faces couples : only it->first is computed. it->second is avoided in the loop for matrices and is updated to it->first in save()
	bool _facesBoundinit,_indexFacePeriodicSet, _isWall; // To ensure that the boundary velocity is initialized after the initial velocity 
	std::vector<double> _Energy;
				

};
#endif /* WAVESTAGGERED_HXX_*/
