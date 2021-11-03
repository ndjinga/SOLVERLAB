//============================================================================
/**
 * \file LinearElasticityModel.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date August 2020
 * \brief Stationary linear elasticity model  
 * -div \sigma = f 
 * with the stress \sigma given by the Hooke's law 
 * \sigma = 2 \mu e(u) + \lambda Tr(e(u)) I_d
 * solved with either finite element or finite volume method
 * Dirichlet (fixed boundary) or Neumann (free boundary) boundary conditions.
 * */
//============================================================================

/*! \class LinearElasticityModel LinearElasticityModel.hxx "LinearElasticityModel.hxx"
 *  \brief Linear Elasticity Model solved with either finite element or finite volume method. 
 * -div \sigma = f 
 * \sigma = 2 \mu e(u) + \lambda Tr(e(u)) I_d
 */
 
#ifndef LinearElasticityModel_HXX_
#define LinearElasticityModel_HXX_

#include "ProblemCoreFlows.hxx"

using namespace std;

/*! Boundary condition type  */
enum BoundaryTypeLinearElasticity	{ NeumannLinearElasticity, DirichletLinearElasticity, NoneBCLinearElasticity};

/** \struct LimitField
 * \brief value of some fields on the boundary  */
struct LimitFieldLinearElasticity{
	LimitFieldLinearElasticity(){bcType=NoneBCLinearElasticity; displacement=0; normalForce=0;}
	LimitFieldLinearElasticity(BoundaryTypeLinearElasticity _bcType, double _displacement,	double _normalForce){
		bcType=_bcType; displacement=_displacement; normalForce=_normalForce;
	}

	BoundaryTypeLinearElasticity bcType;
	double displacement; //for Dirichlet
	double normalForce; //for Neumann
};

class LinearElasticityModel
{

public :
	/** \fn LinearElasticityModel
			 * \brief Constructor for the linear elasticity in a solid
			 * \param [in] int : space dimension
			 * \param [in] bool : numerical method
			 * \param [in] double : solid density
			 * \param [in] double : first Lamé coefficient
			 * \param [in] double : second  Lamé coefficient
			 *  */

	LinearElasticityModel( int dim, bool FECalculation=true, double rho, double lambda, double mu,MPI_Comm comm = MPI_COMM_WORLD);

	void setConstantDensity(double rho) { _rho=rho; }
	void setDensityField(Field densityField) { _densityField=densityField; _densityFieldSet=true;}
	void setLameCoefficients(double lambda, double mu) { _lambda = lambda; _mu = mu;}
	void setYoungAndPoissonModuli(double E, double nu) { _lambda = E*nu/(1+nu)/(1-2*nu); _mu = E/2/(1+nu);}
	void setGravity(Vector gravite ) { _gravity=gravite; }

    void setMesh(const Mesh &M);
    void setFileName(string fileName){	_fileName = fileName;    }
    bool solveStationaryProblem();
    Field getOutputDisplacementField();
    
    //Linear system and spectrum
    void setLinearSolver(linearSolver kspType, preconditioner pcType);
    double getConditionNumber(bool isSingular=false, double tol=1e-6) const;

	//Gestion du calcul
	void initialize();
	void terminate();//vide la mémoire et enregistre le résultat final
	double computeStiffnessMatrix(bool & stop);
	bool solveLinearSystem();//return true if resolution successfull
	void save();

    /* Boundary conditions */
	void setBoundaryFields(map<string, LimitFieldLinearElasticity> boundaryFields){
		_limitField = boundaryFields;
    };
	/** \fn setDirichletBoundaryCondition
			 * \brief adds a new boundary condition of type Dirichlet
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the displacement at the boundary
			 * \param [out] void
			 *  */
	void setDirichletBoundaryCondition(string groupName,double displacement){
		_limitField[groupName]=LimitFieldLinearElasticity(DirichletLinearElasticity,displacement,-1);
	};

	/** \fn setNeumannBoundaryCondition
			 * \brief adds a new boundary condition of type Neumann
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : outward normal force
			 * \param [out] void
			 *  */
	void setNeumannBoundaryCondition(string groupName, double normalForce=0){
		_limitField[groupName]=LimitFieldLinearElasticity(DirichletLinearElasticity,-1, normalForce);
	};

	void setDirichletValues(map< int, double> dirichletBoundaryValues);
	void setNeumannValues  (map< int, double> neumannBoundaryValues);
	

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
	int _PetscIts;//the number of iterations of the linear solver
	Mat  _A;//Linear system matrix
	Vec _b;//Linear system right hand side
	double _MaxIterLinearSolver;//nombre maximum d'iteration gmres obtenu au cours par les resolution de systemes lineaires au cours d'un pas de tmeps
	bool _conditionNumber;//computes an estimate of the condition number

	map<string, LimitFieldLinearElasticity> _limitField;
    bool _onlyNeumannBC;//if true then the linear system is singular and should be solved up to a constant vector
    
	Vector _normale;
    Vec _displacements;//unknown of the linear system
    
	//Physical parameterss
	double _lambda, _mu;//Lamé coefficients
	double _rho;//constantDensity
	Field _densityField;//For non constant density field
	bool _densityFieldSet;
	Vector _gravity;

	//Display variables
	bool _verbose, _system;
	ofstream * _runLogFile;//for creation of a log file to save the history of the simulation
    //saving parameters
	string _fileName;//name of the calculation
	string _path;//path to execution directory used for saving results
	saveFormat _saveFormat;//file saving format : MED, VTK or CSV

	double computeRHS(bool & stop);
	double computeStiffnessMatrixFV(bool & stop);
	double computeStiffnessMatrixFE(bool & stop);

    /************ Data for FE calculation *************/
    bool _FECalculation;
	int _Nnodes;/* number of nodes for FE calculation */
	int _neibMaxNbNodes;/* maximum number of nodes around a node */
	int _NunknownNodes;/* number of unknown nodes for FE calculation */
	int _NboundaryNodes;/* total number of boundary nodes */
	int _NdirichletNodes;/* number of boundary nodes with Dirichlet BC for FE calculation */
    std::vector< int > _boundaryNodeIds;/* List of boundary nodes */
    std::vector< int > _dirichletNodeIds;/* List of boundary nodes with Dirichlet BC */

    /********* Possibility to set a boundary field as Dirichlet boundary condition *********/
    bool _dirichletValuesSet;
    bool _neumannValuesSet;
    std::map< int, double> _dirichletBoundaryValues;
    std::map< int, double> _neumannBoundaryValues;

	//MPI variables
	PetscMPIInt    _size;        /* size of communicator */
	PetscMPIInt    _rank;        /* processor rank */
};

#endif /* LinearElasticityModel_HXX_ */
