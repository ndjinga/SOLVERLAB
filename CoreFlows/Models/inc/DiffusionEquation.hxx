//============================================================================
/**
 * \file DiffusionEquation.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date 24 March 2015
 * \brief Heat diffusion equation
 * rho*cp*dT/dt - \lambda\Delta T = \Phi + \lambda_{sf} (T_{fluid}-T)
 * Dirichlet (imposed temperature) or Neumann (imposed normal flux) boundary conditions.
 * */
//============================================================================

/*! \class DiffusionEquation DiffusionEquation.hxx "DiffusionEquation.hxx"
 *  \brief Scalar heat equation for the Uranium rods temperature
 *  \details see \ref DiffusionEqPage for more details
 * rho*cp*dT/dt - \lambda\Delta T = \Phi + \lambda_{sf} (T_{fluid}-T)
 */
#ifndef DiffusionEquation_HXX_
#define DiffusionEquation_HXX_

#include "ProblemCoreFlows.hxx"
#include "Node.hxx"

using namespace std;

//! enumeration BoundaryType
/*! Boundary condition type  */
enum BoundaryTypeDiffusion	{ NeumannDiffusion, DirichletDiffusion, NoneBCDiffusion};

/** \struct LimitField
 * \brief value of some fields on the boundary  */
struct LimitFieldDiffusion{
	LimitFieldDiffusion(){bcType=NoneBCDiffusion; T=0; normalFlux=0;}
	LimitFieldDiffusion(BoundaryTypeDiffusion _bcType, double _T,	double _normalFlux){
		bcType=_bcType; T=_T; normalFlux=_normalFlux;
	}

	BoundaryTypeDiffusion bcType;
	double T; //for Dirichlet
	double normalFlux; //for Neumann
};

class DiffusionEquation: public ProblemCoreFlows
{

public :
	/** \fn DiffusionEquation
			 * \brief Constructor for the temperature diffusion in a solid
			 * \param [in] int : space dimension
			 * \param [in] double : solid density
			 * \param [in] double : solid specific heat at constant pressure
			 * \param [in] double : solid conductivity
			 *  */

	DiffusionEquation( int dim,bool FECalculation=true,double rho=10000,double cp=300,double lambda=5, MPI_Comm comm = MPI_COMM_WORLD);

	//Gestion du calcul (ICoCo)
	void initialize();
	void terminate();//vide la mémoire et enregistre le résultat final
	bool initTimeStep(double dt);
	double computeTimeStep(bool & stop);//propose un pas de temps pour le calcul. Celà nécessite de discrétiser les opérateur (convection, diffusion, sources) et pour chacun d'employer la condition cfl. En cas de problème durant ce calcul (exemple t=tmax), renvoie stop=true
	void abortTimeStep();//efface les inconnues calculées par solveTimeStep() et reinitialise dt à 0
	bool iterateTimeStep(bool &ok);
	void save();
	void validateTimeStep();

    /* Boundary conditions */
	void setBoundaryFields(map<string, LimitFieldDiffusion> boundaryFields){
		_limitField = boundaryFields;
    };
	/** \fn setDirichletBoundaryCondition
			 * \brief adds a new boundary condition of type Dirichlet on a boundary group
			 * \details Constant bounadry value on the boundary group
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the temperature at the boundary
			 * \param [out] void
			 *  */
	void setDirichletBoundaryCondition(string groupName,double Temperature=0){
		_limitField[groupName]=LimitFieldDiffusion(DirichletDiffusion,Temperature,-1);
	};
	/** \fn setDirichletBoundaryCondition
			 * \brief adds a new boundary condition of type Dirichlet on a boundary group
			 * \details Reads the boundary field in a med file
			 * \param [in] string : the name of the boundary
			 * \param [in] string : the file name
			 * \param [in] string : the field name
			 * \param [in] int : the time step number
			 * \param [in] int : int corresponding to the enum CELLS or NODES
			 * \param [out] void
			 *  */
	void setDirichletBoundaryCondition(string groupName, string fileName, string fieldName, int timeStepNumber, int order, int meshLevel, EntityType field_support_type);
	void setDirichletBoundaryCondition(string groupName, Field bc_field){
		_limitField[groupName]=LimitFieldDiffusion(DirichletDiffusion,0,-1);//This line will be deleted when variable BC are properly treated in solverlab 
	}
	/** \fn setNeumannBoundaryCondition
			 * \brief adds a new boundary condition of type Neumann on a boundary group
			 * \details Constant normal gradient at the boundary group
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the normal gradient
			 * \param [out] void
			 *  */
	void setNeumannBoundaryCondition(string groupName, double normalGradient=0){
		_limitField[groupName]=LimitFieldDiffusion(NeumannDiffusion,-1, normalGradient);
	};
	/** \fn setNeumannBoundaryCondition
			 * \brief adds a new boundary condition of type Neumann on a boundary group
			 * \details Reads the boundary field in a med file
			 * \param [in] string : the name of the boundary
			 * \param [in] string : the file name
			 * \param [in] string : the field name
			 * \param [in] int : the time step number
			 * \param [in] int : int corresponding to the enum CELLS or NODES 
			 * \param [out] void
			 *  */
	void setNeumannBoundaryCondition(string groupName, string fileName, string fieldName, int timeStepNumber, int order, int meshLevel, EntityType field_support_type);
	void setNeumannBoundaryCondition(string groupName, Field BC_Field){
		_limitField[groupName]=LimitFieldDiffusion(NeumannDiffusion,-1, 0);//This line will be deleted when variable BC are properly treated in solverlab 
	};
	/* Set all Dirichlet boundary values */
	void setDirichletValues(map< int, double> dirichletBoundaryValues);
	/* Set all Neumann boundary values */
	void setNeumannValues(map< int, double> neumannBoundaryValues);

	void setRodDensity(double rho){
		_rho=rho;
	};
	void setConductivity(double conductivite){
		_conductivity=conductivite;
	};
	void setDiffusiontensor(Matrix DiffusionTensor){
		_DiffusionTensor=DiffusionTensor;
	};


	/** Set input fields to prepare the simulation or coupling **/
	vector<string> getInputFieldsNames();
	void setInputField(const string& nameField, Field& inputField );//supply of a required input field

	void setFluidTemperatureField(Field coupledTemperatureField);
	void setFluidTemperature(double fluidTemperature){	_fluidTemperature=fluidTemperature;	}
	Field& getFluidTemperatureField(){  return _fluidTemperatureField;	}
	
	/*** get output fields names for postprocessing or coupling ***/
	vector<string> getOutputFieldsNames() ;//liste tous les champs que peut fournir le code pour le postraitement
	Field&         getOutputField(const string& nameField );//Renvoie un champs pour le postraitement

	Field& getOutputTemperatureField();//Return the main unknown if present (initialize() should be called first)
	Field& getRodTemperatureField();//Return the main unknown if present (initialize() should be called first)

    /*********** Generic functions for finite element method ***********/
    static Vector gradientNodal(Matrix M, vector< double > v);//gradient of nodal shape functions
    static int fact(int n);
    static int unknownNodeIndex(int globalIndex, std::vector< int > dirichletNodes);
    static int globalNodeIndex(int unknownIndex, std::vector< int > dirichletNodes);

protected :
	double computeDiffusionMatrix(bool & stop);
	double computeDiffusionMatrixFV(bool & stop);
	double computeDiffusionMatrixFE(bool & stop);
	double computeRHS(bool & stop);

	Field _fluidTemperatureField;
	bool _fluidTemperatureFieldSet, _diffusionMatrixSet;
	double _conductivity,_diffusivity, _fluidTemperature;
	double _rho;
	double _cp;
	Vector _normale;
	Matrix _DiffusionTensor;
	Vec _Tn, _deltaT, _Tk, _Tkm1, _b0;
	Vec _Tn_seq; // Local sequential copy of the parallel vector _Tn, used for saving result files

	double _dt_diffusion, _dt_src;
	map<string, LimitFieldDiffusion> _limitField;

    
    /************ Data for FE calculation *************/
	int _NunknownNodes;/* number of unknown nodes for FE calculation */
	int _NboundaryNodes;/* total number of boundary nodes */
	int _NdirichletNodes;/* number of boundary nodes with Dirichlet BC for FE calculation */
    std::vector< int > _boundaryNodeIds;/* List of boundary nodes */
    std::vector< int > _dirichletNodeIds;/* List of boundary nodes with Dirichlet BC */

    /********* Possibility to set a boundary field as DirichletNeumann boundary condition *********/
    std::map< int, double> _dirichletBoundaryValues;
    std::map< int, double> _neumannBoundaryValues;
};

#endif /* DiffusionEquation_HXX_ */
