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
	IsothermalSinglePhase(phaseType fluid, pressureEstimate pEstimate, int dim);
	//initialisation du systeme
	void initialize();

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
			 * \param [in] vector<double> : the values of the x component of the 2 velocities at the boundary
			 * \param [in] vector<double> : the values of the y component of the 2 velocities at the boundary
			 * \param [in] vector<double> : the values of the z component of the 2 velocities at the boundary
			 * \param [out] void
			 *  */
	void setWallBoundaryCondition(string groupName,vector<double> v_x=vector<double>(3,0), vector<double> v_y=vector<double>(3,0), vector<double> v_z=vector<double>(3,0)){
		_limitField[groupName]=LimitField(Wall,-1,v_x,v_y,v_z,-1,-1,-1,-1);
	};

protected :
	double _Temperature, _internalEnergy;
	double  _drho_sur_dp;
	Field _Pressure, _Density, _Momentum, _Vitesse, _VitesseX, _VitesseY, _VitesseZ, _MachNumber;
	bool _saveAllFields;
	
	//!calcule l'etat de Roe de deux etats
	void convectionState( const long &i, const long &j, const bool &IsBord);
	//!calcule la matrice de convection de l'etat interfacial entre deux cellules voisinnes
	void convectionMatrices();
	//!Calcule le flux pour un état, une porosité et une normale donnés
	Vector convectionFlux(Vector U,Vector V, Vector normale, double porosity);
	//!Computation of the Roe matrix
	void RoeMatrixConservativeVariables(double u_n, double total_enthalpy,Vector velocity, double k, double K);
	void convectionMatrixPrimitiveVariables( double rho, double u_n, double H,Vector velocity);
	//!Computation of the staggered Roe upwinding matrix in conservative variables
	void staggeredRoeUpwindingMatrixConservativeVariables(  double u_n, double total_enthalpy, Vector velocity, double k, double K);
	//!Computation of the staggered Roe upwinding matrix in primitive variables
	void staggeredRoeUpwindingMatrixPrimitiveVariables(double density, double u_n,double total_enthalpy, Vector velocity);
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
	void jacobianConvGhostState(const int &j, string nameOfGroup,double * normale);
	//!Calcule la jacobienne de la CL de diffusion
	void jacobianDiffGhostState(const int &j, string nameOfGroup);
	//!Calcule l'etat fictif à la frontiere
	void setBoundaryState(string nameOfGroup, const int &j,double *normale);
	//!Ajoute au second membre la contribution de la diffusion
	void addDiffusionToSecondMember(const int &i,const int &j,bool isBoundary);
	//!Compute the corrected interfacial state for lowMach, pressureCorrection and staggered versions of the VFRoe formulation
	void applyVFRoeLowMachCorrections(bool isBord, string groupname="");
	//!remplit les vecteurs de scaling
	void computeScaling(double offset);

	// Functions of equations of states
	void consToPrim(const double *Ucons, double* Vprim,double porosity=1);
	void primToCons(const double *V, const int &i, double *U, const int &j);
	void primToConsJacobianMatrix(double *V);
	/** \fn getDensityDerivatives
	 * \brief Computes the partial derivatives of rho, with regard to the primitive variables  p
	 * @param pressure
	*/
	void getDensityDerivatives( double pressure);

};

#endif /* IsothermalSinglePhase_HXX_ */
