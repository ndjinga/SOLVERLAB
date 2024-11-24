//============================================================================
/**
 * \file StationaryDiffusionEquation.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date June 2019
 * \brief Stationary heat diffusion equation solved with either finite elements or finite volume method. 
 * -\lambda\Delta T=\Phi + \lambda_{sf} (T_{fluid}-T)
 * Dirichlet (imposed temperature) or Neumann (imposed normal flux) boundary conditions.
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

using namespace std;

/*! Boundary condition type  */
enum BoundaryTypeStationaryDiffusion    { NeumannStationaryDiffusion, DirichletStationaryDiffusion, NoneBCStationaryDiffusion};

/** \struct LimitField
 * \brief value of some fields on the boundary  */
struct LimitFieldStationaryDiffusion{
    LimitFieldStationaryDiffusion(){bcType=NoneBCStationaryDiffusion; T=0; normalFlux=0;}
    LimitFieldStationaryDiffusion(BoundaryTypeStationaryDiffusion _bcType, double _T,    double _normalFlux){
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
             * \param [in] bool : numerical method
             * \param [in] double : solid conductivity
             *  */

    StationaryDiffusionEquation( int dim,bool FECalculation=true,double lambda=1,MPI_Comm comm = MPI_COMM_WORLD);

    void setMesh(const Mesh &M);
    void setFileName(string fileName){
    _fileName = fileName;
    }
    bool solveStationaryProblem();
    
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

    /* option 1 : Boundary conditions via group name */
    void setBoundaryFields(map<string, LimitFieldStationaryDiffusion> boundaryFields){//Almost never used
        _limitField = boundaryFields;
    };
    /** \fn setDirichletBoundaryCondition
             * \brief adds a new boundary condition of type Dirichlet
             * \details specify the boundary type via the boundary name
             * \param [in] string : the name of the boundary
             * \param [in] double : the value of the temperature at the boundary
             * \param [out] void
             *  */
    void setDirichletBoundaryCondition(string groupName,double Temperature){
        _limitField[groupName]=LimitFieldStationaryDiffusion(DirichletStationaryDiffusion,Temperature,-1);
    };
    /** \fn setNeumannBoundaryCondition
             * \brief adds a new boundary condition of type Neumann
             * \details specify the boundary type via the boundary name
             * \param [in] string : the name of the boundary
             * \param [in] double : outward normal flux
             * \param [out] void
             *  */
    void setNeumannBoundaryCondition(string groupName, double normalFlux=0){
        _limitField[groupName]=LimitFieldStationaryDiffusion(NeumannStationaryDiffusion,-1, normalFlux);
    };

    /* option 2 : Boundary conditions via list of values (do not mix with option 1) */
    void setDirichletValues(map< int, double> dirichletBoundaryValues);
    void setNeumannValues  (map< int, double>   neumannBoundaryValues);
    
    /* option 3 : Boundary conditions via boundary field */
    /** \fn setDirichletBoundaryCondition
             * \brief adds a new boundary condition of type Dirichlet
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
        _limitField[groupName]=LimitFieldStationaryDiffusion(DirichletStationaryDiffusion, 0, -1);
    };

    /** \fn setNeumannBoundaryCondition
             * \brief adds a new boundary condition of type Neumann
             * \details Reads the boundary field in a med file
             * \param [in] string : the name of the boundary
             * \param [in] string : the file name
             * \param [in] string : the field name
             * \param [in] int : the time step number
             * \param [in] int : int corresponding to the enum CELLS or NODES 
             * \param [out] void
             *  */
    void setNeumannBoundaryCondition(string groupName, string fileName, string fieldName, int timeStepNumber, int order, int meshLevel, EntityType field_support_type);
    void setNeumannBoundaryCondition(string groupName, Field bc_field){
        _limitField[groupName]=LimitFieldStationaryDiffusion(NeumannStationaryDiffusion,-1, 0);
    };

    void setConductivity(double conductivite){
        _conductivity=conductivite;
    };
    void setDiffusiontensor(Matrix DiffusionTensor){
        _DiffusionTensor=DiffusionTensor;
    };


    //get input fields to prepare the simulation or coupling
    vector<string> getInputFieldsNames();
    void setInputField(const string& nameField, Field& inputField );//supply of a required input field
    
    void setFluidTemperatureField(Field coupledTemperatureField);
    void setFluidTemperature(double fluidTemperature){    _fluidTemperature=fluidTemperature;    }
    Field& getFluidTemperatureField(){    return _fluidTemperatureField;    }
    
    /** \fn setHeatPowerField
     * \brief set the heat power field (variable in space)
     * \details
     * \param [in] Field
     * \param [out] void
     *  */
    void setHeatPowerField(Field heatPower);

    /** \fn setHeatPowerField
     * \brief set the heat power field (variable in space)
     * \details
     * \param [in] string fileName (including file path)
     * \param [in] string fieldName
     * \param [out] void
     *  */
    void setHeatPowerField(string fileName, string fieldName, int iteration = 0, int order = 0, int meshLevel=0);

    /** \fn getHeatPowerField
     * \brief returns the heat power field
     * \details
     * \param [in] void
     * \param [out] Field
     *  */
    Field getHeatPowerField(){
        return _heatPowerField;
    }
    //get output fields names for postprocessing or coupling
    vector<string> getOutputFieldsNames() ;//liste tous les champs que peut fournir le code pour le postraitement
    Field&         getOutputField(const string& nameField );//Renvoie un champs pour le postraitement

    Field& getOutputTemperatureField();
    Field& getRodTemperatureField();

    /** \fn setVerbose
     * \brief Updates display options
     * \details
     * \param [in] bool
     * \param [in] bool
     * \param [out] void
     *  */
    void setVerbose(bool verbose,  bool system=false)
    {
        _verbose = verbose;
        _system = system;
    };

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
    Vec _Tk_seq; // Local sequential copy of the parallel vector _Tk, used for saving result files

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
    double computeDiffusionMatrixFE(bool & stop);

    /************ Data for FE calculation *************/
    bool _FECalculation;
    int _Nnodes;/* number of nodes for FE calculation */
    int _neibMaxNbNodes;/* maximum number of nodes around a node */
    int _NunknownNodes;/* number of unknown nodes for FE calculation */
    int _NboundaryNodes;/* total number of boundary nodes */
    int _NdirichletNodes;/* number of boundary nodes with Dirichlet BC for FE calculation */
    std::vector< int > _boundaryNodeIds;/* List of boundary nodes */
    std::vector< int > _dirichletNodeIds;/* List of boundary nodes with Dirichlet BC */

    /********* Possibility to set a boundary field as Dirichlet/Neumann boundary condition *********/
    bool _dirichletValuesSet;
    bool _neumannValuesSet;
    std::map< int, double> _dirichletBoundaryValues;
    std::map< int, double> _neumannBoundaryValues;

    /**** MPI related variables ***/
    PetscMPIInt    _mpi_size;        /* size of communicator */
    PetscMPIInt    _mpi_rank;        /* processor rank */
    VecScatter _scat;            /* For the distribution of a local vector */
    int _globalNbUnknowns, _localNbUnknowns;
    int _d_nnz, _o_nnz;            /* local and "non local" numbers of non zeros corfficients */
};

#endif /* StationaryDiffusionEquation_HXX_ */
