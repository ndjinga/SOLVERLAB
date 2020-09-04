//============================================================================
/**
 * \file DiffusionEquation.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date 24 March 2015
 * \brief Heat diffusion equation
 * */
//============================================================================

/*! \class DiffusionEquation DiffusionEquation.hxx "DiffusionEquation.hxx"
 *  \brief Scalar heat equation for the Uranium rods temperature
 *  \details see \ref DiffusionEqPage for more details
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

	DiffusionEquation( int dim,bool FECalculation=true,double rho=10000,double cp=300,double lambda=5);

	//Gestion du calcul
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
			 * \brief adds a new boundary condition of type Dirichlet
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the temperature at the boundary
			 * \param [out] void
			 *  */
	void setDirichletBoundaryCondition(string groupName,double Temperature){
		_limitField[groupName]=LimitFieldDiffusion(DirichletDiffusion,Temperature,-1);
	};
	/** \fn setNeumannBoundaryCondition
			 * \brief adds a new boundary condition of type Neumann
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [out] void
			 *  */
	void setNeumannBoundaryCondition(string groupName, double normalFlux=0){
		_limitField[groupName]=LimitFieldDiffusion(NeumannDiffusion,-1, normalFlux);
	};

	void setRodDensity(double rho){
		_rho=rho;
	};
	void setConductivity(double conductivite){
		_conductivity=conductivite;
	};
	void setFluidTemperatureField(Field coupledTemperatureField){
		_fluidTemperatureField=coupledTemperatureField;
		_fluidTemperatureFieldSet=true;
	};

	void setDiffusiontensor(Matrix DiffusionTensor){
		_DiffusionTensor=DiffusionTensor;
	};

	void setFluidTemperature(double fluidTemperature){
	_fluidTemperature=fluidTemperature;
	}

	//get output fields for postprocessing or coupling
	vector<string> getOutputFieldsNames() ;//liste tous les champs que peut fournir le code pour le postraitement
	Field&         getOutputField(const string& nameField );//Renvoie un champs pour le postraitement

	Field& getRodTemperatureField(){
		return _VV;
	}
	Field& getFluidTemperatureField(){
		return _fluidTemperatureField;
	}

protected :
	double computeDiffusionMatrix(bool & stop);
	double computeDiffusionMatrixFV(bool & stop);
	double computeRHS(bool & stop);

	Field _fluidTemperatureField;
	bool _fluidTemperatureFieldSet, _diffusionMatrixSet;
	double _conductivity,_diffusivity, _fluidTemperature;
	double _rho;
	double _cp;
	Vector _normale;
	Matrix _DiffusionTensor;
	Vec _Tn, _deltaT, _Tk, _Tkm1, _b0;
	double _dt_diffusion, _dt_src;
    
    /************ Data for FE calculation *************/
    bool _FECalculation;
	int _NunknownNodes;/* number of unknown nodes for FE calculation */
	int _NboundaryNodes;/* total number of boundary nodes */
	int _NdirichletNodes;/* number of boundary nodes with Dirichlet BC for FE calculation */
    std::vector< int > _boundaryNodeIds;/* List of boundary nodes */
    std::vector< int > _dirichletNodeIds;/* List of boundary nodes with Dirichlet BC */

    /*********** Functions for finite element method ***********/
    Vector gradientNodal(Matrix M, vector< double > v);//gradient of nodal shape functions
	double computeDiffusionMatrixFE(bool & stop);
    int fact(int n);
    int unknownNodeIndex(int globalIndex, std::vector< int > dirichletNodes);
    int globalNodeIndex(int unknownIndex, std::vector< int > dirichletNodes);

	TimeScheme _timeScheme;
	map<string, LimitFieldDiffusion> _limitField;
};

#endif /* DiffusionEquation_HXX_ */
