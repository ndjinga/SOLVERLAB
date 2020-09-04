//============================================================================
/**
 * \file StationaryDiffusionEquation.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date June 2019
 * \brief Stationary heat diffusion equation solved with either finite elements or finite volume method. 
 * -\lambda\Delta T=\Phi + \lambda_{sf} (T_{fluid}-T)
 * Dirichlet (imposed temperature) or Neumann (imposed normal flux) boundary conditions
 * */
//============================================================================

/*! \class StationaryDiffusionEquation StationaryDiffusionEquation.hxx "StationaryDiffusionEquation.hxx"
 *  \brief Scalar stationary heat equation solved with either finite elements or finite volume method. 
 *  \details see \ref StationaryDiffusionEqPage for more details
 * -\lambda\Delta T=\Phi(T) + \lambda_{sf} (T_{fluid}-T)
 */
#ifndef StationaryDiffusionEquation_HXX_
#define StationaryDiffusionEquation_HXX_

#include "ProblemCoreFlows.hxx"

/* for the laplacian spectrum */
#include <slepceps.h>
#include <slepcsvd.h>

using namespace std;

//! enumeration BoundaryType
/*! Boundary condition type  */
enum BoundaryTypeStationaryDiffusion	{ NeumannStationaryDiffusion, DirichletStationaryDiffusion, NoneBCStationaryDiffusion};

/** \struct LimitField
 * \brief value of some fields on the boundary  */
struct LimitFieldStationaryDiffusion{
	LimitFieldStationaryDiffusion(){bcType=NoneBCStationaryDiffusion; T=0; normalFlux=0;}
	LimitFieldStationaryDiffusion(BoundaryTypeStationaryDiffusion _bcType, double _T,	double _normalFlux){
		bcType=_bcType; T=_T; normalFlux=_normalFlux;
	}

	BoundaryTypeStationaryDiffusion bcType;
	double T; //for Dirichlet
	double normalFlux; //for Neumann
};

class StationaryDiffusionEquation
{

public :
	/** \fn StationaryDiffusionEquation
			 * \brief Constructor for the temperature diffusion in a solid
			 * \param [in] int : space dimension
			 * \param [in] double : solid conductivity
			 *  */

	StationaryDiffusionEquation( int dim,bool FECalculation=true,double lambda=1);

    void setMesh(const Mesh &M);
    void setFileName(string fileName){
	_fileName = fileName;
    }
    bool solveStationaryProblem();
    Field getOutputTemperatureField();
    
    //Linear system and spectrum
    void setLinearSolver(linearSolver kspType, preconditioner pcType);
    double getConditionNumber(bool isSingular=false, double tol=1e-6) const;
    std::vector< double > getEigenvalues (int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6) const;
    std::vector< Vector > getEigenvectors(int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6) const;
    Field getEigenvectorsField(int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6) const;

	//Gestion du calcul
	void initialize();
	void terminate();//vide la mémoire et enregistre le résultat final
	double computeDiffusionMatrix(bool & stop);
    double computeTimeStep(bool & stop);//For coupling calculations
	bool iterateNewtonStep(bool &ok);
	void save();

    /* Boundary conditions */
	void setBoundaryFields(map<string, LimitFieldStationaryDiffusion> boundaryFields){
		_limitField = boundaryFields;
    };
	/** \fn setDirichletBoundaryCondition
			 * \brief adds a new boundary condition of type Dirichlet
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the temperature at the boundary
			 * \param [out] void
			 *  */
	void setDirichletBoundaryCondition(string groupName,double Temperature){
		_limitField[groupName]=LimitFieldStationaryDiffusion(DirichletStationaryDiffusion,Temperature,-1);
	};

	/** \fn setNeumannBoundaryCondition
			 * \brief adds a new boundary condition of type Neumann
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [out] void
			 *  */
	void setNeumannBoundaryCondition(string groupName, double normalFlux=0){
		_limitField[groupName]=LimitFieldStationaryDiffusion(NeumannStationaryDiffusion,-1, normalFlux);
	};

	void setDirichletValues(map< int, double> dirichletBoundaryValues);
	void setNeumannValues  (map< int, double>   neumannBoundaryValues);
	
	void setConductivity(double conductivite){
		_conductivity=conductivite;
	};
	void setFluidTemperatureField(Field coupledTemperatureField){
		_fluidTemperatureField=coupledTemperatureField;
		_fluidTemperatureFieldSet=true;
	};
	void setFluidTemperature(double fluidTemperature){
	_fluidTemperature=fluidTemperature;
	}
	Field& getRodTemperatureField(){
		return _VV;
	}
	Field& getFluidTemperatureField(){
		return _fluidTemperatureField;
	}
	void setDiffusiontensor(Matrix DiffusionTensor){
		_DiffusionTensor=DiffusionTensor;
	};

	/** \fn setHeatPowerField
	 * \brief set the heat power field (variable in space)
	 * \details
	 * \param [in] Field
	 * \param [out] void
	 *  */
	void setHeatPowerField(Field heatPower){
		_heatPowerField=heatPower;
		_heatPowerFieldSet=true;
	}

	/** \fn setHeatPowerField
	 * \brief set the heat power field (variable in space)
	 * \details
	 * \param [in] string fileName (including file path)
	 * \param [in] string fieldName
	 * \param [out] void
	 *  */
	void setHeatPowerField(string fileName, string fieldName){
		_heatPowerField=Field(fileName, CELLS,fieldName);
		_heatPowerFieldSet=true;
	}

protected :
	//Main unknown field
	Field _VV;

	int _Ndim;//space dimension
	int _nVar;//Number of equations to solve=1

    //Mesh data
	Mesh _mesh;
    bool _meshSet;
	bool _initializedMemory;
	int _Nmailles;//number of cells for FV calculation
	int _neibMaxNbCells;//maximum number of cells around a cell
    
	double _precision;
	double _precision_Newton;
	double _erreur_rel;//norme(Uk+1-Uk)
    bool _computationCompletedSuccessfully;
    
	//Linear solver and petsc
	KSP _ksp;
	KSPType _ksptype;
	PC _pc;
	PCType _pctype;
	string _pc_hypre;
	int _maxPetscIts;//nombre maximum d'iteration gmres autorisé au cours d'une resolution de système lineaire
	int _maxNewtonIts;//nombre maximum d'iteration de Newton autorise au cours de la resolution d'un pas de temps
	int _PetscIts;//the number of iterations of the linear solver
	int _NEWTON_its;
	Mat  _A;//Linear system matrix
	Vec _b;//Linear system right hand side
	double _MaxIterLinearSolver;//nombre maximum d'iteration gmres obtenu au cours par les resolution de systemes lineaires au cours d'un pas de tmeps
	bool _conditionNumber;//computes an estimate of the condition number

	map<string, LimitFieldStationaryDiffusion> _limitField;
    bool _onlyNeumannBC;//if true then the linear system is singular and should be solved up to a constant vector
    
	bool _diffusionMatrixSet;
	Vector _normale;
	Matrix _DiffusionTensor;
	Vec _deltaT, _Tk, _Tkm1, _b0;
	double _dt_src;
    
	//Heat transfert variables
	double _conductivity, _fluidTemperature;
	Field _heatPowerField, _fluidTemperatureField;
	bool _heatPowerFieldSet, _fluidTemperatureFieldSet;
	double _heatTransfertCoeff, _heatSource;

	//Display variables
	bool _verbose, _system;
	ofstream * _runLogFile;//for creation of a log file to save the history of the simulation
    //saving parameters
	string _fileName;//name of the calculation
	string _path;//path to execution directory used for saving results
	saveFormat _saveFormat;//file saving format : MED, VTK or CSV

	double computeRHS(bool & stop);
	double computeDiffusionMatrixFV(bool & stop);

    /************ Data for FE calculation *************/
    bool _FECalculation;
	int _Nnodes;/* number of nodes for FE calculation */
	int _neibMaxNbNodes;/* maximum number of nodes around a node */
	int _NunknownNodes;/* number of unknown nodes for FE calculation */
	int _NboundaryNodes;/* total number of boundary nodes */
	int _NdirichletNodes;/* number of boundary nodes with Dirichlet BC for FE calculation */
    std::vector< int > _boundaryNodeIds;/* List of boundary nodes */
    std::vector< int > _dirichletNodeIds;/* List of boundary nodes with Dirichlet BC */

    /*********** Functions for finite element method ***********/
    Vector gradientNodal(Matrix M, vector< double > v);//gradient of nodal shape functions
	double computeDiffusionMatrixFE(bool & stop);
    int fact(int n);
    int unknownNodeIndex(int globalIndex, std::vector< int > dirichletNodes) const;
    int globalNodeIndex(int unknownIndex, std::vector< int > dirichletNodes) const;

    /********* Possibility to set a boundary field as DirichletNeumann boundary condition *********/
    bool _dirichletValuesSet;
    bool _neumannValuesSet;
    std::map< int, double> _dirichletBoundaryValues;
    std::map< int, double> _neumannBoundaryValues;
};

#endif /* StationaryDiffusionEquation_HXX_ */
