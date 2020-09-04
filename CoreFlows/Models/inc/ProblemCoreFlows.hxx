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

#include <petsc.h>
#include <slepceps.h>
#include <slepcsvd.h>

#include "Field.hxx"
#include "Mesh.hxx"
#include "Cell.hxx"
#include "Face.hxx"
#include "CdmathException.hxx"

using namespace std;

//! enumeration linearSolver
/*! the linearSolver can be GMRES or BiCGStab (see Petsc documentation) */
enum linearSolver
{
	GMRES,/**< linearSolver is GMRES */
	BCGS,/**< linearSolver is BiCGSstab */
	CG/**< linearSolver is CG */
};

//! enumeration preconditioner
/*! the preconditioner can be ILU or LU  (see Petsc documentation) */
enum preconditioner
{
	ILU,/**< preconditioner is ILU(0) */
	LU,/**< preconditioner is actually a direct solver (LU factorisation)*/
	NONE,/**< no preconditioner used */
	ICC,/**< preconditioner is ICC(0) */
	CHOLESKY/**< preconditioner is actually a direct solver for symmetric matrices (CHOLESKY factorisation)*/
};

//! enumeration saveFormat
/*! the numerical results are saved using MED, VTK or CSV format */
enum saveFormat
{
	MED,/**< MED format is used  */
	VTK,/**< VTK format is used */
	CSV/**< CSV format is used */
};

//! enumeration TimeScheme
/*! The numerical method can be Explicit or Implicit  */
enum TimeScheme
{
	Explicit,/**<  Explicit numerical scheme */
	Implicit/**< Implicit numerical scheme */
};

class ProblemCoreFlows
{
public :
	//! Constructeur par défaut
	ProblemCoreFlows();
	virtual ~ProblemCoreFlows();
	// -*-*-*- Gestion du calcul -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*


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
	 *  \details c'est une fonction virtuelle ,
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
	virtual bool iterateTimeStep(bool &converged) = 0; //??

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

	/*
	//Coupling interface
	virtual vector<string> getInputFieldsNames()=0 ;//Renvoie les noms des champs dont le problème a besoin (données initiales)
	virtual  Field& getInputFieldTemplate(const string& name)=0;//Renvoie le format de champs attendu (maillage, composantes etc)
	virtual void setInputField(const string& name, const Field& afield)=0;//enregistre les valeurs d'une donnée initiale
	virtual vector<string> getOutputFieldsNames()=0 ;//liste tous les champs que peut fournir le code pour le postraitement
	virtual Field& getOutputField(const string& nameField )=0;//Renvoie un champs pour le postraitement
	 */

	//paramètres du calcul -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

	/** \fn setPresentTime
	 * \brief met à jour _time (le temps courant du calcul)
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void setPresentTime (double time);

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
	 * \brief sets the initial field from a field in a med file
	 * \details
	 * \param [in] string : the file name
	 * \param [in] string : the field name
	 * \param [in] int : the time step number
	 * \param [out] void
	 *  */
	void setInitialField(string fileName, string fieldName, int timeStepNumber);

	/** \fn setInitialFieldConstant
	 * \brief sets a constant initial field on a mesh stored in a med file
	 * \details
	 * \param [in] string : the file name
	 * \param [in] vector<double> : the value in each cell
	 * \param [out] void
	 *  */
	void setInitialFieldConstant(string fileName, const vector<double> Vconstant);

	/** \fn setInitialFieldConstant
	 * \brief sets a constant initial field 
	 * \details
	 * \param [in] Mesh 
	 * \param [in] Vector
	 * \param [out] void
	 *  */
	void setInitialFieldConstant(const Mesh& M, const Vector Vconstant);

	/** \fn setInitialFieldConstant
	 * \brief sets a constant initial field
	 * \details
	 * \param [in] Mesh
	 * \param [in] vector<double>
	 * \param [out] void
	 *  */
	void setInitialFieldConstant(const Mesh& M, const vector<double> Vconstant);

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
	 * \param [out] void
	 *  */
	void setInitialFieldConstant( int nDim, const vector<double> Vconstant, double xmin, double xmax,int nx, string leftSide, string rightSide,
			double ymin=0, double ymax=0, int ny=0, string backSide="", string frontSide="",
			double zmin=0, double zmax=0, int nz=0, string bottomSide="", string topSide="");

	/** \fn setInitialFieldStepFunction
	 * \brief sets a step function initial field (Riemann problem)
	 * \details
	 * \param [in] Mesh
	 * \param [in] Vector
	 * \param [in] Vector
	 * \param [in] double position of the discontinuity on one of the three axis
	 * \param [in] int direction (axis carrying the discontinuity) : 0 for x, 1 for y, 2 for z
	 * \param [out] void
	 *  */
	void setInitialFieldStepFunction(const Mesh M, const Vector Vleft, const Vector Vright, double disc_pos, int direction=0);

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
	 * \param [out] void
	 *  */
	void setInitialFieldStepFunction( int nDim, const vector<double> VV_Left, vector<double> VV_Right, double xstep,
			double xmin, double xmax,int nx, string leftSide, string rightSide,
			double ymin=0, double ymax=0, int ny=0, string backSide="", string frontSide="",
			double zmin=0, double zmax=0, int nz=0, string bottomSide="", string topSide="");

	/** \fn setInitialFieldSphericalStepFunction
	 * \brief sets a step function initial field with value Vin inside the ball with radius Radius and Vout outside
	 * \details
	 * \param [in] Mesh
	 * \param [in] Vector Vin, value inside the ball
	 * \param [in] Vector Vout, value outside the ball
	 * \param [in] double radius of the ball
	 * \param [in] Vector Center, coordinates of the ball center
	 * \param [out] void
	 *  */
	void setInitialFieldSphericalStepFunction(const Mesh M, const Vector Vin, const Vector Vout, double Radius, Vector Center);

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
	 * \brief sets the linear solver and preconditioner
	 * \details virtual function overloaded by intanciable classes
	 * @param kspType linear solver type (GMRES or BICGSTAB)
	 * @param pcType preconditioner (ILU,LU or NONE)
	 */
	void setLinearSolver(linearSolver solverName, preconditioner pcType);

	/** \fn setNewtonSolver
	 * \brief set the Newton algorithm parameters
	 * \param [in] int maximum number of newton iterations
	 * \param [in] double precision required for the convergence of the newton scheme
	 * \param [out] void
	 *  */
	void setNewtonSolver(double precision,int iterations=20)
	{
		_maxNewtonIts=iterations;
		_precision_Newton=precision;
	};

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

	//Couplages Thermohydraulique-thermique-neutronique *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

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

	/** \fn setHeatSource
	 * \brief met à jour la puissance thermique ( _phi )
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void setHeatSource(double phi){
		_heatSource=phi;
	}

	/** \fn getHeatPowerField
	 * \brief renvoie le champs ?? ( _heatPowerField )
	 * \details
	 * \param [in] void
	 * \param [out] Field
	 *  */
	Field getHeatPowerField(){
		return _heatPowerField;
	}

	/** \fn setRodTemperatureField ??
	 * \brief
	 * \details
	 * \param [in] Field
	 * \param [out] void
	 *  */
	void setRodTemperatureField(Field rodTemperature){
		_rodTemperatureField=rodTemperature;
		_rodTemperatureFieldSet=true;
	}

	/** \fn setRodTemperature ??
	 * \brief
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void setRodTemperature(double rodTemp){
		_rodTemperature=rodTemp;
	}

	/** \fn getRodTemperatureField
	 * \brief
	 * \details
	 * \param [in] void
	 * \param [out] Field
	 *  */
	virtual Field& getRodTemperatureField(){ // ?? je ne retrouve pas cet attribut dans le file.cxx
		return _rodTemperatureField;
	}

	/** \fn setHeatTransfertCoeff
	 * \brief set the heat transfert coefficient for heat exchange between fluid and solid
	 * \details
	 * \param [in] double
	 * \param [out] void
	 *  */
	void setHeatTransfertCoeff(double heatTransfertCoeff){
		_heatTransfertCoeff=heatTransfertCoeff;
	}

	/** \fn setDISPLAY
	 * \brief met à jour les paramètres de l'affichage
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

    //Spectrum analysis
    double getConditionNumber(bool isSingular=false, double tol=1e-6) const;
    std::vector< double > getEigenvalues (int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6) const;
    std::vector< Vector > getEigenvectors(int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6) const;
    Field getEigenvectorsField(int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6) const;

	//  some supplementary functions

	/** \fn displayMatrix
	 * \brief displays a matrix of size "size x size" for profiling
	 * @param  matrix is a pointer of size "size"
	 * @param size, size of the matrix
	 * @param name, string, name or description of the matrix
	 * @return displays the matrix on the terminal
	 *  */
	void displayMatrix(double *matrix, int size, string name);

	/** \fn displayMatrix
	 * \brief displays a vector of size "size" for profiling
	 * @param  vector is a pointer of size "size"
	 * @param size, size of the vector
	 * @param name, string, name or description of the vector
	 * @return displays the vector on the terminal
	 *  */
	void displayVector(double *vector, int size, string name);

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


protected :

	int _Ndim;//space dimension
	int _nVar;//Number of equations to sole
	int _Nmailles;//number of cells
	int _Nnodes;//number of nodes
	int _Nfaces;//number of faces
	int _neibMaxNb;//maximum number of neighbours around a cell
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
	int _maxPetscIts;//nombre maximum d'iteration gmres autorise au cours d'une resolution de systeme lineaire
	int _PetscIts;//the number of iterations of the linear solver
	int _maxNewtonIts;//nombre maximum d'iteration de Newton autorise au cours de la resolution d'un pas de temps
	int _NEWTON_its;
	Mat  _A;//Linear system matrix
	Vec _b;//Linear system right hand side
	double _MaxIterLinearSolver;//nombre maximum d'iteration gmres obtenu au cours par les resolution de systemes lineaires au cours d'un pas de tmeps
	bool _conditionNumber;//computes an estimate of the condition number

	//simulation monitoring variables
	bool _isStationary;
	bool _initialDataSet;
	bool _initializedMemory;
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
	Field _heatPowerField, _rodTemperatureField;
	bool _heatPowerFieldSet, _rodTemperatureFieldSet;
	double _heatTransfertCoeff;
	double _heatSource, _rodTemperature;
	double _hsatv, _hsatl;//all models appart from DiffusionEquation will need this

	//Display variables
	bool _verbose, _system;

	string _path;//path to execution directory used for saving results
	saveFormat _saveFormat;//file saving format : MED, VTK or CSV
};

#endif /* PROBLEMCOREFLOWS_HXX_ */
