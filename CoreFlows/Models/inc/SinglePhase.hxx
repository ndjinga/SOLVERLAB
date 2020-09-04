//============================================================================
/**
 * \file DriftModel.hxx
 * \author Michael NDJINGA, Kieu Nguyen
 * \version 1.0
 * \date 01 janv. 2016
 * \brief The compressible Navier-Stokes equations
 * */
//============================================================================

/*! \class SinglePhase SinglePhase.hxx "SinglePhase.hxx"
 *  \brief The compressible Navier-Stokes equations
 *  \details The model consists in one mass, one momentum and one energy equation, see \ref NSModelsPage for more details
 */
#ifndef SINGLEPHASE_HXX_
#define SINGLEPHASE_HXX_

#include "ProblemFluid.hxx"

class SinglePhase : public ProblemFluid{
public :
	/** \fn SinglePhase
	 * \brief Constructor for the Navier-Stokes system
	 * \param [in] phaseType : \ref Liquid or \ref Gas
	 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
	 * \param [in] int : mesh dimension
	 * \param [in] bool : There are two possible equations of state for the fluid
	 *  */
	SinglePhase(phaseType fluid, pressureEstimate pEstimate,int dim,bool useDellacherieEOS=false);

	/** \fn setViscosity
	 * \brief sets the viscosity
	 * @param viscosite : value of the dynamic viscosity
	 * 	 * */
	void setViscosityConstant( double viscosite ){
		_fluides[0]->setViscosity(viscosite);
	};

	//! system initialisation
	void initialize();

	//fonctions d'echange de flux
	//	void getOutputField(const Vec &Flux, const string Champ, const int numBord)=0;//, PetscInt *indices_Flux, PetscInt *indices_Bord, const long range)=0;
	//	double trace(const int &numBord, Vec &out)=0;
	void testConservation();

	void save();

	/** \fn setIntletBoundaryCondition
	 * \brief adds a new boundary condition of type Inlet
	 * \details
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [in] double : the value of the x component of the velocity at the boundary
	 * \param [in] double : the value of the y component of the velocity at the boundary
	 * \param [in] double : the value of the z component of the velocity at the boundary
	 * \param [out] void
	 *  */
	void setInletBoundaryCondition(string groupName,double Temperature,double v_x=0, double v_y=0, double v_z=0){
		_limitField[groupName]=LimitField(Inlet,-1,vector<double>(1,v_x),vector<double>(1,v_y),vector<double>(1,v_z),Temperature,-1,-1,-1);
	};
	/** \fn setIntletRotationBoundaryCondition
	 * \brief adds a new boundary condition of type InletRotationVelocity
	 * \details
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [in] double : the components of the rotational of the inlet velocity
	 * \param [out] void
	 *  */
	void setInletRotationBoundaryCondition(string groupName,double Temperature,double omega_x=0, double omega_y=0, double omega_z=0){
		_limitField[groupName]=LimitField(InletRotationVelocity,0,vector<double>(1,omega_x),vector<double>(1,omega_y),vector<double>(1,omega_z),Temperature,-1,-1,-1);
	};
	/** \fn setIntletPressureBoundaryCondition
	 * \brief adds a new boundary condition of type InletPressure
	 * \details
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the pressure at the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [in] double : the components of the rotational of the inlet velocity
	 * \param [out] void
	 *  */
	void setInletPressureBoundaryCondition(string groupName, double pressure,double Temperature,double omega_x=0, double omega_y=0, double omega_z=0){
		_limitField[groupName]=LimitField(InletPressure,pressure,vector<double>(1,omega_x),vector<double>(1,omega_y),vector<double>(1,omega_z),Temperature,-1,-1,-1);
	};
	/** \fn setIntletPressureBoundaryCondition
	 * \brief adds a new boundary condition of type InletPressure taking into account the hydrostatic pressure variations
	 * \details The pressure is not constant on the boundary but varies linearly with a slope given by the gravity vector
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the pressure at the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [in] vector<double> : reference_point position on the boundary where the value Pressure will be imposed
	 * \param [out] void
	 *  */
	void setInletPressureBoundaryCondition(string groupName, double pressure,double Temperature, vector<double> reference_point){
		/* On the boundary we have P-Pref=rho g\cdot(x-xref) hence P=Pref-g\cdot xref + g\cdot x */
		pressure-=reference_point[0]*_GravityField3d[0];
		if(_Ndim>1){
			pressure-=reference_point[1]*_GravityField3d[1];
			if(_Ndim>2)
				pressure-=reference_point[2]*_GravityField3d[2];
		}

		_limitField[groupName]=LimitField(InletPressure,pressure,vector<double>(0,0),vector<double>(0,0),vector<double>(0,0),Temperature,-1,-1,-1);
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
	void setWallBoundaryCondition(string groupName,double Temperature=-1, double v_x=0, double v_y=0, double v_z=0){
		if(Temperature!=-1)
			_limitField[groupName]=LimitField(Wall,-1,vector<double>(1,v_x),vector<double>(1,v_y),vector<double>(1,v_z),Temperature,-1,-1,-1);
		else
			_limitField[groupName]=LimitField(Wall,-1,vector<double>(1,v_x),vector<double>(1,v_y),vector<double>(1,v_z),getReferenceTemperature(),-1,-1,-1);
	};

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

	double getReferencePressure()    { return _Pref; };
	double getReferenceTemperature() { return _Tref; };
	
	//get output fields for postprocessing or coupling
	vector<string> getOutputFieldsNames() ;//liste tous les champs que peut fournir le code pour le postraitement
	Field&         getOutputField(const string& nameField );//Renvoie un champs pour le postraitement
	Field& getPressureField();
	Field& getVelocityField();
	Field& getVelocityXField();
	Field& getTemperatureField();
	Field& getDensityField();
	Field& getMomentumField();
	Field& getTotalEnergyField();
	Field& getEnthalpyField();

protected :
	double  _drho_sur_dp,   _drho_sur_dT;//derivatives of the density rho wrt cv, p, T
	double  _drhoE_sur_dp,  _drhoE_sur_dT;//derivatives of the total energy rho E wrt cv, p, T
	bool _useDellacherieEOS;
	double _Tref; //EOS reference temperature
	double _Pref; //EOS reference pressure

	//!calcule l'etat de Roe de deux etats
	void convectionState( const long &i, const long &j, const bool &IsBord);
	//!calcule la matrice de convection de l'etat interfacial entre deux cellules voisinnes
	void convectionMatrices();
	//!Calcule le flux pour un état et une porosité et une normale donnés
	Vector convectionFlux(Vector U,Vector V, Vector normale, double porosity);
	//!calcule la matrice de diffusion de l'etat interface pour la diffusion
	void diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord);
	//!Computes the source vector associated to the cell i
	void sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i);
	//!Computes the pressure loss associated to the face ij
	void pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj);
	//!Computes the contribution of the porosity gradient associated to the face ij to the source term
	void porosityGradientSourceVector();
	//!Calcule la jacobienne de la CL convection
	void jacobian(const int &j, string nameOfGroup,double * normale);
	//!Calcule la jacobienne de la CL de diffusion
	void jacobianDiff(const int &j, string nameOfGroup);
	//!Calcule l'etat fictif a la frontiere
	void setBoundaryState(string nameOfGroup, const int &j,double *normale);// delete &nf Kieu
	//!Adds the contribution of diffusion to the RHS
	void addDiffusionToSecondMember(const int &i,const int &j,bool isBord);
	//!Computation of the Roe matrix
	void RoeMatrixConservativeVariables(double u_n, double total_enthalpy,Vector velocity, double k, double K);
	void convectionMatrixPrimitiveVariables( double rho, double u_n, double H,Vector velocity);
	//!Computation of the staggered Roe upwinding matrix in conservative variables
	void staggeredRoeUpwindingMatrixConservativeVariables(  double u_n, double total_enthalpy, Vector velocity, double k, double K);
	//!Computation of the staggered Roe upwinding matrix in primitive variables
	void staggeredRoeUpwindingMatrixPrimitiveVariables(double density, double u_n,double total_enthalpy, Vector velocity);
	/**** staggered VFFC scheme has some specific features (no need for Roe matrix), hence some specific functions ******/
	//!Computes the interfacial flux for the VFFC formulation of the staggered upwinding
	Vector staggeredVFFCFlux();
	//!Computes the matrices A^+ and A^- for the VFFC formulation of the staggered upwinding
	void staggeredVFFCMatricesConservativeVariables(double u_n);
	//!Computes the matrices A^+ and A^- for the VFFC formulation of the staggered upwinding using Primitive Variables
	void staggeredVFFCMatricesPrimitiveVariables(double u_n);
	//!Compute the corrected interfacial state for lowMach, pressureCorrection and staggered versions of the VFRoe formulation
	void applyVFRoeLowMachCorrections(bool isBord, string groupname="");
	//!Special preconditioner based on a matrix scaling strategy
	void computeScaling(double offset);
	//!Calcule les saut de valeurs propres pour la correction entropique
	void entropicShift(double* n);
	// Fonctions utilisant la loi d'etat
	void consToPrim(const double *Ucons, double* Vprim,double porosity=1);
	void primToCons(const double *V, const int &i, double *U, const int &j);
	void primToConsJacobianMatrix(double *V);
	/** \fn getDensityDerivatives
	 * \brief Computes the partial derivatives of rho, and rho E with regard to the primitive variables  p and  T
	 * @param pressure
	 * @param temperature
	 * @param square of the velocity vector
	*/
	void getDensityDerivatives( double pressure, double temperature, double v2);

	bool _saveAllFields;
	Field _Enthalpy, _Pressure, _Density, _Temperature, _Momentum, _TotalEnergy, _Vitesse, _VitesseX, _VitesseY, _VitesseZ;

	};
#endif /* SINGLEPHASE_HXX_*/
