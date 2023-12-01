//============================================================================
// Name        : ProblemCoreFlows
// Author      : M. Ndjinga
// Version     :
// Copyright   : CEA Saclay 2014
// Description : Generic class for PDEs problems
//============================================================================
/* A ProblemCoreFlows class */

/*! \class ProblemCoreFlows ProblemCoreFlows.hxx "ProblemCoreFlows.hxx"
 *  \brief Generic class for thermal hydraulics problems
 *  \details  Common functions to CoreFlows models
 */

#ifndef PROBLEMCOREFLOWS_HXX_
#define PROBLEMCOREFLOWS_HXX_

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <string>
#include <map>

#include <petscksp.h>
#include <slepceps.h>
#include <slepcsvd.h>

#include "Field.hxx"
#include "Mesh.hxx"
#include "Cell.hxx"
#include "Face.hxx"
#include "CdmathException.hxx"

#include "MEDCouplingFieldDouble.hxx"

using namespace std;

//! enumeration linearSolver
/*! the linearSolver can be GMRES or BiCGStab (see Petsc documentation) */
enum linearSolver
{
	GMRES,/**< linearSolver is GMRES */
	BCGS,/**< linearSolver is BiCGSstab */
	CG,/**< linearSolver is Conjugate Gradient */
	CGNE,/**< linearSolver is Conjugate Gradient for Normal Equations */
	FGMRES/**< linearSolver is GMRES */
};

//! enumeration preconditioner
/*! the preconditioner can be ILU or LU  (see Petsc documentation) */
enum preconditioner
{
	ILU,/**< preconditioner is ILU(0) */
	LU,/**< preconditioner is actually a direct solver (LU factorisation)*/
	NOPC,/**< no preconditioner used */
	ICC,/**< preconditioner is ICC(0) */
	CHOLESKY,/**< preconditioner is actually a direct solver for symmetric matrices (CHOLESKY factorisation)*/
	GAMG,/**< multigrid preconditioner */
	QR,/**< preconditioner is actually a direct solver (QR factorisation)*/
	Svd/**< preconditioner is actually a direct solver (SVD decomposition)*/
};

//! enumeration nonLinearSolver
/*! the nonLinearSolver can be Newton_SOLVERLAB or using PETSc, Newton_PETSC_LINESEARCH, Newton_PETSC_TRUSTREGION (see Petsc documentation) */
enum nonLinearSolver
{
	Newton_SOLVERLAB,/**< nonLinearSolver is Newton_SOLVERLAB */
	Newton_PETSC_LINESEARCH,/**< nonLinearSolver is Newton_PETSC_LINESEARCH */
	Newton_PETSC_LINESEARCH_BASIC,/**< nonLinearSolver is Newton_PETSC_LINESEARCH_BASIC */
	Newton_PETSC_LINESEARCH_BT,/**< nonLinearSolver is Newton_PETSC_LINESEARCH_BT */
	Newton_PETSC_LINESEARCH_SECANT,/**< nonLinearSolver is Newton_PETSC_LINESEARCH_SECANT */
	Newton_PETSC_LINESEARCH_NLEQERR,/**< nonLinearSolver is Newton_PETSC_LINESEARCH_LEQERR */
	Newton_PETSC_TRUSTREGION,/**< nonLinearSolver is Newton_PETSC_TRUSTREGION */
	Newton_PETSC_NGMRES,/**< nonLinearSolver is Newton_PETSC_NGMRES */
	Newton_PETSC_ASPIN/**< nonLinearSolver is Newton_PETSC_ASPIN */
};

//! enumeration saveFormat
/*! the numerical results are saved using MED, VTK or CSV format */
enum saveFormat
{
	MED,/**< MED format is used  */
	VTK,/**< VTK format is used */
	CSV,/**< CSV format is used */
	NOSAVING/**< Do not save results to hard disk */
};

//! enumeration TimeScheme
/*! The numerical method can be Explicit or Implicit  */
enum TimeScheme
{
	Explicit,/**<  Explicit numerical scheme */
	Implicit/**< Implicit numerical scheme */
};

//! enumeration pressureEstimate
/*! the pressure estimate needed to fit physical parameters  */
enum pressureEstimate
{
	around1bar300K,/**< pressure is around 1 bar and temperature around 300K (for TransportEquation, SinglePhase and IsothermalTwoFluid) or 373 K (saturation for DriftModel and FiveEqsTwoFluid) */
	around155bars600K/**< pressure is around 155 bars  and temperature around 618 K (saturation) */
};

class ProblemCoreFlows
{
public :
	//! Constructeur par défaut
	ProblemCoreFlows(MPI_Comm comm = MPI_COMM_WORLD);
	virtual ~ProblemCoreFlows();
	
	// -*-*-*- Gestion du calcul (interface ICoCo) -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

	/** \fn initialize
	 * \brief Alloue la mémoire et vérifie que  le maillage et les conditions limites/initiales sont bien définis
	 * \Details c'est une fonction virtuelle pure, on la surcharge dans le problem Fluide .
	 * @param  void
	 *  */
	virtual void initialize()=0;

	/** \fn terminate
	 * \brief vide la mémoire et enregistre le résultat final
	 * \Details on la surcharge dans le problem fluid
	 *  @param void
	 *  */
	virtual void terminate()=0;

	/** \fn computeTimeStep
	 * \brief Propose un pas de temps pour le calcul
	 * \details Pour proposer un pas de temps la fonction nécessite de discrétiser les opérateurs
	 *  convection, diffusion et sources. Pour chacun on emploie la condition cfl.
	 *   En cas de problème durant ce calcul (exemple t=tmax), la fonction renvoie stop=true sinon stop = false
	 *   Cest une fonction virtuelle pure, on la surcharge dans Problème Fulide
	 *  @param Stop booléen correspond a True si ya un problème dans le code (exemple t=tmax) False sinon
	 *  \return  dt le pas de temps
	 *  */
	virtual double computeTimeStep(bool & stop)=0;

	/** \fn initTimeStep
	 * \brief Enregistre le nouveau pas de temps dt et l'intègre aux opérateurs
	 *  discrets (convection, diffusion, sources)
	 *  \details c'est une fonction virtuelle pure, on la surcharge dans le problemfluid
	 *  @param  dt est le nouvel pas de temps (double)
	 *  \return false si dt <0 et True sinon
	 *  			   */
	virtual bool initTimeStep(double dt)=0;

	/** \fn solveTimeStep
	 * \brief calcule les valeurs inconnues au pas de temps +1 .
	 *  \details c'est une fonction virtuelle
	 *  @param  void
	 *  \return Renvoie false en cas de problème durant le calcul (valeurs non physiques..)
	 *  */
	virtual bool solveTimeStep();//

	/** \fn validateTimeStep
	 * \brief Valide le calcule au temps courant
	 * \details met à jour le temps présent t=t+dt, sauvegarde les champs inconnus
	 * et reinitialise dt à 0, teste la stationnarité .
	 * c'est une fonction virtuel , on la surchage dans chacun des modèles
	 * @param  void
	 * \return  Renvoie false en cas de problème durant le calcul
	 *  */
	virtual void validateTimeStep()=0;//

	/** \fn abortTimeStep
	 * \brief efface les inconnues calculées par solveTimeStep() et reinitialise dt à 0
	 * \details c'est une fonction virtuelle pure, elle est surchargée dans ProblemFluid, TransportEquation et DiffusionEquation
	 *  */
	virtual void abortTimeStep()=0;

	/** \fn run
	 * \brief Vérifie tous les paramètres et lance le code en supposant que la cfl ne changera pas
	 * \details  Renvoie "false" en cas de problème durant le calcul
	 * C'est une fonction virtuelle pure, on la surcharge dans problème Fluide
	 * @param void
	 * \return Renvoie "false" en cas de problème durant le calcul
	 *  */
	virtual bool run();//

	/** \fn iterateTimeStep
	 * \brief Calcul d'une sous-itération du pas de temps en cours, typiquement une itération de Newton pour un schéma implicite
	 * \details Deux paramètres booléen (converged et ok) sont retournés
	 * converged, Vaut true si on peut passer au pas de temps suivant (shéma explicite ou convergence du schéma de Newton dans le schéma implicite)
	 * ok  vaut true si le calcul n'a pas rencontré d'erreur ou de problème particulier dans la mise à jour des variables physiques
	 * \param [in] bool, passage par reférence.
	 * \param [out] bool
	 *  */
	virtual bool iterateTimeStep(bool &converged) = 0; 

	/** \fn isStationary
	 * \brief vérifie la stationnairité du problème .
	 * \details Renvoie "true" si le problème atteint l'état stationnaire et "false" sinon
	 * \param [in] void
	 * \param [out] bool
	 *  */
	virtual bool isStationary() const;

	/** \fn presentTime
	 * \brief Calcule la valeur du temps courant .
	 * \details
	 * \param [in] void
	 * \param [out] double
	 *  */
	virtual double presentTime() const;

	/** \fn setStationaryMode
	 * \brief Perform the search of a stationary regime
	 * \details
	 * \param [in] bool
	 * \param [out] 
	 *  */
	virtual void setStationaryMode(bool stationaryMode){ _stationaryMode=stationaryMode;};

	/** \fn getStationaryMode
	 * \brief Tells if we are seeking a stationary regime
	 * \details
	 * \param [in] 
	 * \param [out] bool
	 *  */
	virtual bool getStationaryMode(){return _stationaryMode;};

	/** \fn resetTime
	 * \brief sets the current time (typically to start a new calculation)
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void resetTime (double time);

	/*
	//Coupling interface
	virtual void getInputMEDDoubleFieldTemplate(const std::string& name, MEDDoubleField& afield) const;//Renvoie le format de champs attendu (maillage, composantes etc)
	virtual vector<string> getOutputFieldsNames()=0 ;//liste tous les champs que peut fournir le code pour le postraitement
	virtual void getOutputMEDDoubleField(const std::string& name, MEDDoubleField& afield)//Renvoie un champs pour le postraitement ou le couplage
	 */

	/** Set input fields to prepare the simulation or coupling **/
	virtual vector<string> getInputFieldsNames()=0;
	virtual void setInputField(const string& nameField, Field& inputField )=0;//supply of a required input field

     /*! @brief (Optional) Provide the code with a scalar double data.
     *
     * See Problem documentation for more details on the time semantic of a scalar value.
     *
     * @param[in] name name of the scalar value that is given to the code.
     * @param[in] val value passed to the code.
     * @throws ICoCo::WrongArgument exception if the scalar name ('name' parameter) is invalid.
     */
    /* virtual void setInputDoubleValue(const std::string& name, const double& val); */

    /*! @brief (Optional) Retrieve a scalar double value from the code.
     *
     * See Problem documentation for more details on the time semantic of a scalar value.
     *
     * @param[in] name name of the scalar value to be read from the code.
     * @return the double value read from the code.
     * @throws ICoCo::WrongArgument exception if the scalar name ('name' parameter) is invalid.
     */
    /* virtual double getOutputDoubleValue(const std::string& name) const; */

    /*! @brief (Optional) Similar to setInputDoubleValue() but for an int value.
     * @sa setInputDoubleValue()
     */
    /* virtual void setInputIntValue(const std::string& name, const int& val); */

    /*! @brief (Optional) Similar to getOutputDoubleValue() but for an int value.
     * @sa getOutputDoubleValue()
     */
    /* virtual int getOutputIntValue(const std::string& name) const; */

    /*! @brief (Optional) Similar to setInputDoubleValue() but for an string value.
     * @sa setInputDoubleValue()
     */
    /* virtual void setInputStringValue(const std::string& name, const std::string& val); */

    /*! @brief (Optional) Similar to getOutputDoubleValue() but for an string value.
     * @sa getOutputDoubleValue()
     */
    /* virtual std::string getOutputStringValue(const std::string& name) const; */
    
   /*! @brief Return ICoCo interface major version number.
     * @return ICoCo interface major version number (2 at present)
     */
    static int GetICoCoMajorVersion() { return 2; }

    /*! @brief (Optional) Get MEDCoupling major version, if the code was built with MEDCoupling support.
     *
     * This can be used to assess compatibility between codes when coupling them.
     *
     * @return the MEDCoupling major version number (typically 7, 8, 9, ...)
     */
    virtual int getMEDCouplingMajorVersion() const{ return MEDCOUPLING_VERSION_MAJOR; };

    /*! @brief (Optional) Indicate whether the code was built with a 64-bits version of MEDCoupling.
     *
     * Implemented if the code was built with MEDCoupling support.
     * This can be used to assess compatibility between codes when coupling them.
     *
     * @return the MEDCoupling major version number
     */
    virtual bool isMEDCoupling64Bits() const;

    /*! @brief (Optional) Get the list of input scalars accepted by the code.
     *
     * @return the list of scalar names that represent inputs of the code
     * @throws ICoCo::WrongContext exception if called before initialize() or after terminate().
     */
    /* virtual std::vector<std::string> getInputValuesNames() const; */

    /*! @brief (Optional) Get the list of output scalars that can be provided by the code.
     *
     * @return the list of scalar names that can be returned by the code
     * @throws ICoCo::WrongContext exception if called before initialize() or after terminate().
     */
    /* virtual std::vector<std::string> getOutputValuesNames() const; */

	Field getUnknownField() const;
	
	//paramètres du calcul -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

	/** \fn setMaxNbOfTimeStep
	 * \brief met à jour _maxNbOfTimeStep ( le nombre maximum d'itération du calcul )
	 * \details
	 * \param [in] int
	 * \param [out] void
	 *  */
	void setMaxNbOfTimeStep(int maxNbOfTimeStep);

	/** \fn setTimeMax
	 * \brief met à jour _timeMax (Le temps maximum du calcul)
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void setTimeMax(double timeMax);

	/** \fn setCFL
	 * \brief met à jour la _CFL
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void setCFL(double cfl);

	/** \fn setPrecision
	 * \brief met à jour _precision (la précision du calcule)
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void setPrecision(double precision);

	/** \fn setInitialField
	 * \brief sets the initial field
	 * \details
	 * \param [in] Field
	 * \param [out] void
	 *  */
	void setInitialField(const Field &VV);

	/** \fn setInitialField
	 * \brief sets the initial field
	 * \details
	 * \param [in] MEDCouplingField*
	 * \param [out] void
	 *  */
	void setInitialField( MEDCoupling::MEDCouplingFieldDouble* myMEDCouplingield );

	/** \fn setInitialField
	 * \brief sets the initial field
	 * \details
	 * \param [in] MCAuto<MEDCoupling::MEDCouplingFieldDouble>
	 * \param [out] void
	 *  */
	void setInitialField( const MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> myMEDCouplingield );

	/** \fn setInitialField
	 * \brief sets the initial field from a field in a med file
	 * \details
	 * \param [in] string : the file name
	 * \param [in] string : the field name
	 * \param [in] int : the time step number
	 * \param [in] EntityType : CELLS, NODES or FACES
	 * \param [out] void
	 *  */
	void setInitialField(string fileName, string fieldName, int timeStepNumber, int order = 0, int meshLevel=0, EntityType typeField = CELLS);

	/** \fn setInitialFieldConstant
	 * \brief sets a constant initial field on a mesh stored in a med file
	 * \details
	 * \param [in] string : the file name
	 * \param [in] vector<double> : the value in each cell
	 * \param [in] EntityType : CELLS, NODES or FACES
	 * \param [out] void
	 *  */
	void setInitialFieldConstant(string fileName, const vector<double> Vconstant, EntityType typeField = CELLS);

	/** \fn setInitialFieldConstant
	 * \brief sets a constant initial field 
	 * \details
	 * \param [in] Mesh 
	 * \param [in] Vector
	 * \param [in] EntityType : CELLS, NODES or FACES
	 * \param [out] void
	 *  */
	void setInitialFieldConstant(const Mesh& M, const Vector Vconstant, EntityType typeField = CELLS);

	/** \fn setInitialFieldConstant
	 * \brief sets a constant initial field
	 * \details
	 * \param [in] Mesh
	 * \param [in] vector<double>
	 * \param [in] EntityType : CELLS, NODES or FACES
	 * \param [out] void
	 *  */
	void setInitialFieldConstant(const Mesh& M, const vector<double> Vconstant, EntityType typeField = CELLS);

	/** \fn setInitialFieldConstant
	 * \brief sets a constant initial field
	 * \details
	 * \param [in] int the space dimension
	 * \param [in] vector<double> the value in each cell
	 * \param [in] double the lowest value in the x direction
	 * \param [in] double the highest value in the x direction
	 * \param [in] string name of the left boundary
	 * \param [in] string name of the right boundary
	 * \param [in] double the lowest value in the y direction
	 * \param [in] double the highest value in the y direction
	 * \param [in] string name of the back boundary
	 * \param [in] string name of the front boundary
	 * \param [in] double the lowest value in the z direction
	 * \param [in] double the highest value in the z direction
	 * \param [in] string name of the bottom boundary
	 * \param [in] string name of the top boundary
	 * \param [in] EntityType : CELLS, NODES or FACES
	 * \param [out] void
	 *  */
	void setInitialFieldConstant( int nDim, const vector<double> Vconstant, double xmin, double xmax,int nx, string leftSide, string rightSide,
			double ymin=0, double ymax=0, int ny=0, string backSide="", string frontSide="",
			double zmin=0, double zmax=0, int nz=0, string bottomSide="", string topSide="", EntityType typeField = CELLS);

	/** \fn setInitialFieldStepFunction
	 * \brief sets a step function initial field (Riemann problem)
	 * \details
	 * \param [in] Mesh
	 * \param [in] Vector
	 * \param [in] Vector
	 * \param [in] double position of the discontinuity on one of the three axis
	 * \param [in] int direction (axis carrying the discontinuity) : 0 for x, 1 for y, 2 for z
	 * \param [in] EntityType : CELLS, NODES or FACES
	 * \param [out] void
	 *  */
	void setInitialFieldStepFunction(const Mesh M, const Vector Vleft, const Vector Vright, double disc_pos, int direction=0, EntityType typeField = CELLS);

	/** \fn setInitialFieldStepFunction
	 * \brief sets a constant initial field
	 * \details
	 * \param [in] int the space dimension
	 * \param [in] vector<double> the value left of the discontinuity
	 * \param [in] vector<double> the value right of the discontinuity
	 * \param [in] double the position of the discontinuity in the x direction
	 * \param [in] double the lowest value in the x direction
	 * \param [in] double the highest value in the x direction
	 * \param [in] string name of the left boundary
	 * \param [in] string name of the right boundary
	 * \param [in] double the lowest value in the y direction
	 * \param [in] double the highest value in the y direction
	 * \param [in] string name of the back boundary
	 * \param [in] string name of the front boundary
	 * \param [in] double the lowest value in the z direction
	 * \param [in] double the highest value in the z direction
	 * \param [in] string name of the bottom boundary
	 * \param [in] string name of the top boundary
	 * \param [in] EntityType : CELLS, NODES or FACES
	 * \param [out] void
	 *  */
	void setInitialFieldStepFunction( int nDim, const vector<double> VV_Left, vector<double> VV_Right, double xstep,
			double xmin, double xmax,int nx, string leftSide, string rightSide,
			double ymin=0, double ymax=0, int ny=0, string backSide="", string frontSide="",
			double zmin=0, double zmax=0, int nz=0, string bottomSide="", string topSide="", EntityType typeField = CELLS);

	/** \fn setInitialFieldSphericalStepFunction
	 * \brief sets a step function initial field with value Vin inside the ball with radius Radius and Vout outside
	 * \details
	 * \param [in] Mesh
	 * \param [in] Vector Vin, value inside the ball
	 * \param [in] Vector Vout, value outside the ball
	 * \param [in] double radius of the ball
	 * \param [in] Vector Center, coordinates of the ball center
	 * \param [in] EntityType : CELLS, NODES or FACES
	 * \param [out] void
	 *  */
	void setInitialFieldSphericalStepFunction(const Mesh M, const Vector Vin, const Vector Vout, double Radius, Vector Center, EntityType typeField = CELLS);

	/** \fn getTime
	 * \brief renvoie _time (le temps courant du calcul)
	 * \details
	 * \param [in] void
	 * \param [out] double
	 *  */
	double getTime();

	/** \fn getNbTimeStep
	 * \brief renvoie _nbTimeStep le Numéro d'itération courant
	 * \details
	 * \param [in] void
	 * \param [out] unsigned
	 *  */
	unsigned getNbTimeStep();

	/** \fn getCFL
	 * \brief renvoie la _CFL
	 * \details
	 * \param [in] void
	 * \param [out] double
	 *  */
	double getCFL();

	/** \fn getPrecision
	 * \brief renvoie _precision (la précision du calcul)
	 * \details
	 * \param [in] void
	 * \param [out] double
	 *  */
	double getPrecision();

	/** \fn getMesh
	 * \brief renvoie _Mesh (le maillage du problème)
	 * \details
	 * \param [in] void
	 * \param [out] Mesh
	 *  */
	Mesh getMesh();

	/** \fn setFileName
	 * \brief met à jour _fileName le nom du fichier
	 * \details
	 * \param [in]  string
	 * \param [out] void
	 *  */
	void setFileName(string fileName);

	/** \fn setFreqSave
	 * \brief met à jour _FreqSave (la fréquence du sauvgarde de la solution)
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void setFreqSave(int freqSave);

	/** \fn save
	 * \brief sauvgarde les données dans des fichiers MED ou VTK
	 * \details c'est une fonction virtuelle pure , on la surcharge
	 * dans chacun des modèles
	 * @param  void
	 */
	virtual void save() = 0;

	/** \fn getLinearSolver
	 * \brief renvoie _ksptype (le type du solveur linéaire utilisé)
	 * \details
	 * \param [in] void
	 * \param [out] string
	 *  */
	string getLinearSolver() {
		return _ksptype;
	};

	/** \fn getNonLinearSolver
	 * \brief renvoie _nonLinearSolver (le type du solveur de Newton utilisé)
	 * \details
	 * \param [in] void
	 * \param [out] string
	 *  */
	nonLinearSolver getNonLinearSolver() {
		return _nonLinearSolver;
	};

	/** \fn getNumberOfVariables
	 * \brief le nombre d'inconnues du problème
	 * \details
	 * @param void
	 * \return renvoie _nVar (le nombre d'inconnues du problème)
	 *  */
	int getNumberOfVariables(){
		return _nVar;
	};

	/** \fn setWellBalancedCorrection
	 * \brief include a well balanced correction to treat stiff source terms
	 * @param boolean that is true if a well balanced correction should be applied
	 * */
	void setWellBalancedCorrection(bool wellBalancedCorr){
		_wellBalancedCorrection=wellBalancedCorr;
	}

	/** \fn setLinearSolver
	 * \brief Legacy function that sets the linear solver and preconditioner
	 * @param Three choices of linear solvers (GMRES, BICGSTAB or CG)
	 * @param Five choices of preconditioner (ILU,LU, ICC, CHOLESKY or NOPC)
	 */
	void setLinearSolver(linearSolver solverName, preconditioner pcType, double maxIts=50);

	/** \fn setLinearSolver
	 * \brief sets the linear solver and preconditioner
	 * @param Any available solver in PETSc (See PETSc KSPType)
	 * @param Any available preconditioner in PETSc (See PETSc PCType)
	 */
	void setLinearSolver(std::string solverName, std::string pcName, double maxIts=50);

	/** \fn setNewtonSolver
	 * \brief sets the Newton type algorithm for solving the nonlinear algebraic system arising from the discretisation of the PDE
	 * \param [in] double : precision required for the convergence of the newton scheme
	 * \param [in] int : maximum number of newton iterations
	 * \param [in] nonLinearSolver : the algorithm to be used to solve the nonlinear system
	 * \param [out] void
	 *  */
	void setNewtonSolver(double precision, int iterations=20, nonLinearSolver solverName=Newton_SOLVERLAB);

	/** \fn displayConditionNumber
	 * \brief display the condition number of the preconditioned linear systems
	 */
	void displayConditionNumber(bool display=true){
		_conditionNumber=display;
	}

	/** \fn setSaveFileFormat
	 * \brief sets the numerical results file format (MED, VTK or CSV)
	 * \details
	 * \param [in] saveFormat
	 * \param [out] void
	 *  */
	void setSaveFileFormat(saveFormat saveFileFormat){
		_saveFormat=saveFileFormat;
	}

	/** \fn setResultDirectory
	 * \brief sets the directory where the results will be saved
	 * \details
	 * \param [in] resultsPath
	 * \param [out] void
	 *  */
	void setResultDirectory(string resultsPath){
		_path=resultsPath;
	}

	/** \fn getTimeScheme
	 * \brief returns the  time scheme name
	 * \param [in] void
	 * \param [out] enum TimeScheme (explicit or implicit)
	 *  */
	TimeScheme getTimeScheme();

	/** \fn setNumericalScheme
	 * \brief sets the numerical method ( explicit vs implicit )
	 * \details
	 * \param [in] TimeScheme
	 * \param [out] void
	 *  */
	void setTimeScheme( TimeScheme method);

	//Couplages Thermohydraulique-thermique-neutronique *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

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

	/** \fn setHeatSource
	 * \brief sets a constant heat power field
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void setHeatSource(double phi){
		_heatSource=phi;
		_isStationary=false;//Source term may be changed after previously reaching a stationary state
	}

	/** \fn getHeatPowerField
	 * \brief returns the heat power field
	 * \details
	 * \param [in] void
	 * \param [out] Field
	 *  */
	Field getHeatPowerField(){
		return _heatPowerField;
	}

	/** \fn setHeatTransfertCoeff
	 * \brief set the heat transfert coefficient for heat exchange between fluid and solid
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void setHeatTransfertCoeff(double heatTransfertCoeff){
		_heatTransfertCoeff=heatTransfertCoeff;
		_isStationary=false;//Source term may be changed after previously reaching a stationary state
	}

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

    //Spectral analysis
    double getConditionNumber(bool isSingular=false, double tol=1e-6) const;
    std::vector< double > getEigenvalues (int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6, EPSType type = EPSKRYLOVSCHUR, bool viewEigenvaluesInXWindows=false, double pause_lenght=0) const;
    std::vector< Vector > getEigenvectors(int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6) const;
    Field getEigenvectorsField(int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6) const;
    std::vector< double > getSingularValues( int nsv, SVDWhich which=SVD_SMALLEST, double tol=1e-6, SVDType type = SVDCYCLIC, bool viewSingularValuesInXWindows=false, double pause_lenght=0) const;
    std::vector< Vector > getSingularVectors(int nsv, SVDWhich which=SVD_SMALLEST, double tol=1e-6) const;

	//  some supplementary functions

	/** \fn displayMatrix
	 * \brief displays a matrix of size "size x size" for profiling
	 * @param  matrix is a pointer of size "size"
	 * @param size, size of the matrix
	 * @param name, string, name or description of the matrix
	 * @return displays the matrix on the terminal
	 *  */
	static void displayMatrix(double *matrix, int size, string name="Matrix coefficients :");

	/** \fn displayMatrix
	 * \brief displays a vector of size "size" for profiling
	 * @param  vector is a pointer of size "size"
	 * @param size, size of the vector
	 * @param name, string, name or description of the vector
	 * @return displays the vector on the terminal
	 *  */
	static void displayVector(double *vector, int size, string name="Vector coefficients :");

protected :

	// Mesh info
	int _Ndim;//space dimension
	int _nVar;//Number of equations to solve
	int _Nmailles;//number of cells
	int _Nnodes;//number of nodes
	int _Nfaces;//number of faces
	int _neibMaxNbCells;//maximum number of neighbours around a cell
	int _neibMaxNbNodes;/* maximum number of nodes around a node */
	Mesh _mesh;
	Field _perimeters;

	//Main unknown field
	Field _VV;

	//Numerical method
	double _dt;
	double _cfl;
	double _maxvp;//valeur propre max pour calcul cfl
	double _minl;//minimum cell diameter
    bool _FECalculation;
	/** boolean used to specify that a well balanced correction should be used */
	bool _wellBalancedCorrection;
	TimeScheme _timeScheme;

	//Linear solver and petsc
	KSP _ksp;
	KSPType _ksptype;
	PC _pc;
	PCType _pctype;
	string _pc_hypre;
	nonLinearSolver _nonLinearSolver;
	int _maxPetscIts;//nombre maximum d'iteration gmres autorise au cours d'une resolution de systeme lineaire
	int _PetscIts;//the number of iterations of the linear solver
	int _maxNewtonIts;//nombre maximum d'iteration de Newton autorise au cours de la resolution d'un pas de temps
	int _NEWTON_its;
	Mat _A;//Linear system matrix
	Vec _b;//Linear system right hand side
	int _MaxIterLinearSolver;//nombre maximum d'iteration gmres obtenu au cours par les resolution de systemes lineaires au cours d'un pas de tmeps
	bool _conditionNumber;//computes an estimate of the condition number
	/** \fn createKSP
	 * \brief Create PETSc solver and preconditioner structures
	 *  */
	void createKSP();

	// Simulation monitoring variables
	bool _isStationary;
	bool _initialDataSet;
	bool _initializedMemory;
	bool _stationaryMode;//ICoCo V2
	bool _restartWithNewTimeScheme;
	bool _restartWithNewFileName;
	double _timeMax,_time;
	int _maxNbOfTimeStep,_nbTimeStep;
	double _precision;
	double _precision_Newton;
	double _erreur_rel;//norme(Un+1-Un)
	string _fileName;//name of the calculation
	int _freqSave;
	ofstream * _runLogFile;//for creation of a log file to save the history of the simulation

	//Heat transfert variables
	Field _heatPowerField;
	bool _heatPowerFieldSet;
	double _heatTransfertCoeff;
	double _heatSource;
	double _hsatv, _hsatl;//all models appart from DiffusionEquation will need this

	//Display variables
	bool _verbose, _system;

	string _path;//path to execution directory used for saving results
	saveFormat _saveFormat;//file saving format : MED, VTK or CSV
	
	//MPI related variables
	PetscMPIInt    _mpi_size;        /* size of communicator */
	PetscMPIInt    _mpi_rank;        /* processor rank */
	VecScatter	   _scat;			/* For the distribution of a local vector */
	int _globalNbUnknowns, _localNbUnknowns;
	int _d_nnz, _o_nnz;			/* local and "non local" numbers of non zeros corfficients */
};

#endif /* PROBLEMCOREFLOWS_HXX_ */
