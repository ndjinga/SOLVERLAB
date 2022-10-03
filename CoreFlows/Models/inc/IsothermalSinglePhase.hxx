//============================================================================
/**
 * \file IsothermalSinglePhase.hxx
 * \author Michael NDJINGA, Esteban Coiffier
 * \date 07 July 2022
 * \brief The isothermal single phase model
 * */
//============================================================================

/*! \class IsothermalSinglePhase IsothermalSinglePhase.hxx "IsothermalSinglePhase.hxx"
 *  \brief Isothermal single phase model
 *  \details The model consists in two phasic mass equations, two phasic momentum equations, see \ref IsothermalPage for more details
 */
#ifndef IsothermalSinglePhase_HXX_
#define IsothermalSinglePhase_HXX_

#include "ProblemFluid.hxx"

class IsothermalSinglePhase : public ProblemFluid{
public :
	/** \fn IsothermalSinglePhase
			 * \brief Constructor for isothermal single phase model
			 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
			 * \param [in] int : mesh dimension
			 *  */
	IsothermalSinglePhase(phaseType fluid, pressureEstimate pEstimate, int dim, bool isCompressibleFluid=true);
	//!initialisation du systeme (allocations mémoire)
	void initialize();
	//!libération de la mémoire
	void terminate();

	void testConservation();

	void save();

	// Boundary conditions
	/** \fn setIntletBoundaryCondition
			 * \brief adds a new boundary condition of type Inlet
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] vector<double> : the values of the x component of the 2 velocities at the boundary
			 * \param [in] vector<double> : the values of the y component of the 2 velocities at the boundary
			 * \param [in] vector<double> : the values of the z component of the 2 velocities at the boundary
			 * \param [out] void
			 *  */
	void setInletBoundaryCondition(string groupName, vector<double> v_x=vector<double>(3,0), vector<double> v_y=vector<double>(3,0), vector<double> v_z=vector<double>(3,0)){
		_limitField[groupName]=LimitField(Inlet,-1,v_x,v_y,v_z,-1,-1,-1,-1);
	};
	/** \fn setIntletPressureBoundaryCondition
			 * \brief adds a new boundary condition of type InletPressure
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the pressure at the boundary
			 * \param [out] void
			 *  */
	void setInletPressureBoundaryCondition(string groupName, double Pressure){
		_limitField[groupName]=LimitField(InletPressure,Pressure,vector<double>(0,0),vector<double>(0,0),vector<double>(0,0),-1,-1,-1,-1);
	};
	/** \fn setWallBoundaryCondition
	 * \brief adds a new boundary condition of type Wall
	 * \details
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [in] double : the value of the x component of the velocity at the boundary
	 * \param [in] double : the value of the y component of the velocity at the boundary
	 * \param [in] double : the value of the z component of the velocity at the boundary
	 * \param [out] void
	 *  */
	void setWallBoundaryCondition(string groupName, double v_x=0, double v_y=0, double v_z=0){
		_limitField[groupName]=LimitField(Wall,-1,vector<double>(1,v_x),vector<double>(1,v_y),vector<double>(1,v_z),getReferenceTemperature(),-1,-1,-1);
	};
	/** \fn addConvectionToSecondMember
	 * \brief Adds the contribution of the convection to the system right hand side for a face (i,j) inside the domain
	 * @param i left cell number
	 * @param j right cell number
	 * @param isBord is a boolean that is true if the interface (i,j) is a boundary interface
	 * @param groupname : is a string that may be used when isBord is true to specify which boundary the face (i,j) belongs to
	 * */
	void addConvectionToSecondMember(const int &i,const int &j,bool isBord, string groupname="");

	/** \fn addSourceTermToSecondMember
	 * \brief Adds the contribution of source terms to the right hand side of the system: gravity,
	 *  phase change, heat power and friction
	 * @param i,j : left and right cell number
	 * @param nbNeighboursi, integer giving the number of neighbours of cell i
	 * @param nbNeighboursj, integer giving the number of neighbours of cell j
	 * @param boolean isBoundary is true for a boundary face (i,j) and false otherwise
	 * @param double mesureFace the lenght or area of the face
	 * */
	void addSourceTermToSecondMember(const int i, int nbNeighboursi,const int j, int nbNeighboursj,bool isBoundary, int ij, double mesureFace);

	//EOS functions
	double getReferenceTemperature() { return _Temperature; };
	
	/* Get output fields for postprocessing or coupling */
	vector<string> getOutputFieldsNames() ;//liste tous les champs que peut fournir le code pour le postraitement
	Field&         getOutputField(const string& nameField );//Renvoie un champ pour le postraitement
	Field& getPressureField();
	Field& getVelocityField();
	Field& getVelocityXField();
	Field& getVelocityYField();
	Field& getVelocityZField();
	Field& getDensityField();
	Field& getMomentumField();
	Field& getMachNumberField();

protected :
	//Thermodynamical quantities
	double _Temperature, _internalEnergy;
	double  _drho_sur_dp;
	//Saving results
	Field _Pressure, _Density, _Momentum, _Vitesse, _VitesseX, _VitesseY, _VitesseZ, _MachNumber;
	bool _saveAllFields;
	
	//Vecteurs nécessaires pour utilisation variables primitives dans schéma de Newton
	double * _Vextdiff, *_Vext, *_Vdiff;
	
	//!calcule l'etat de Roe de deux etats
	void convectionState( const long &i, const long &j, const bool &IsBord);
	//!calcule la matrice de convection de l'etat interfacial entre deux cellules voisinnes
	void convectionMatrices();
	//!Calcule le flux pour un état, une porosité et une normale donnés
	Vector convectionFlux(Vector U,Vector V, Vector normale, double porosity);
	//!Computation of the Roe matrix in primitive variables
	void convectionMatrixPrimitiveVariables(double u_n);
	//!Computation of the Roe matrix in conservative variables
	void convectionMatrixConservativeVariables(double u_n);
	//!Computation of the staggered Roe upwinding matrix in primitive variables
	void staggeredRoeUpwindingMatrixPrimitiveVariables( double u_n);
	//!calcule la matrice de diffusion de l'etat interface pour la diffusion
	void diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord);
	//!Ajoute au second membre la contribution de la gravite et du frottement
	void sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i);
	//!Computes the pressure loss associated to the face ij
	void pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj);
	//!Computes the contribution of the porosity gradient associated to the face ij to the source term
	void porosityGradientSourceVector();
	//matrice de gravite
	// void gravityMatrix();
	//!Calcule la jacobienne de la CL convection
	void jacobian(const int &j, string nameOfGroup,double * normale);
	//!Calcule la jacobienne de la CL de diffusion
	void jacobianDiff(const int &j, string nameOfGroup);
	//!Calcule l'etat fictif à la frontiere
	void setBoundaryState(string nameOfGroup, const int &j,double *normale);
	//!Ajoute au second membre la contribution de la diffusion
	void addDiffusionToSecondMember(const int &i,const int &j,bool isBoundary);
	//!Compute the corrected interfacial state for lowMach, pressureCorrection and staggered versions of the VFRoe formulation
	void applyVFRoeLowMachCorrections(bool isBord, string groupname="");
	//!remplit les vecteurs de scaling
	void computeScaling(double offset);

	// Functions of equations of states
	bool _isCompressibleFluid;
	CompressibleFluid *_compressibleFluid;//This class works with both compressible and incompressible fluids
	void consToPrim(const double *Ucons, double* Vprim,double porosity=1);
	void primToCons(const double *V, const int &i, double *U, const int &j);
	void primToConsJacobianMatrix(double *V);
	void primToConsJacobianMatrix(double rho, double* velocity, double invSoundSpeedsquared);
	/** \fn getDensityDerivatives
	 * \brief Computes the partial derivatives of rho, with regard to the primitive variables  p
	 * @param pressure
	*/
	void getDensityDerivatives( double pressure);
	
	//Functions not yet implemented (generate exceptions)
	Vector staggeredVFFCFlux();
	void entropicShift(double* n);

	bool _isSingularSystem;
	Vec _constantPressureVector;
};

#endif /* IsothermalSinglePhase_HXX_ */
