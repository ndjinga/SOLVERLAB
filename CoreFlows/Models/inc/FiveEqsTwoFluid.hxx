//============================================================================
/**
 * \file DriftModel.hxx
 * \author Kieu Nuyen, Michael NDJINGA
 * \version 1.0
 * \date 01 Sept. 2014
 * \brief Five equation two phase flow model
 * */
//============================================================================

/*! \class FiveEqsTwoFluid FiveEqsTwoFluid.hxx "FiveEqsTwoFluid.hxx"
 *  \brief The model consists in the phasic mass and momentum balance equations and one mixture total energy balance equation.
 *  \details The model consists in two phasic mass equations, two phasic momentum equations, one mixture energy equation, see \ref FiveEqPage for more details
 */

#ifndef FiveEqsTwoFluid_HXX_
#define FiveEqsTwoFluid_HXX_

#include "ProblemFluid.hxx"

class FiveEqsTwoFluid : public ProblemFluid{
  public :
	/** \fn FiveEqsTwoFluid
			 * \brief Constructor for the five equation two-fluid model with two velocities and one temperature
			 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
			 * \param [in] int : mesh dimension
			 *  */
	FiveEqsTwoFluid(pressureEstimate pEstimate, int dim);
	//initialisation du systeme
	void initialize();

	void testConservation();

	void save();

	// Boundary conditions
	/** \fn setIntletBoundaryCondition
			 * \brief adds a new boundary condition of type Inlet
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the vapour volume fraction at the boundary
			 * \param [in] double : the value of the temperature at the boundary
			 * \param [in] vector<double> : the values of the x component of the 2 velocities at the boundary
			 * \param [in] vector<double> : the values of the y component of the 2 velocities at the boundary
			 * \param [in] vector<double> : the values of the z component of the 2 velocities at the boundary
			 * \param [out] void
			 *  */
	void setInletBoundaryCondition(string groupName,double alpha,double Temperature,vector<double> v_x, vector<double> v_y=vector<double>(3,0), vector<double> v_z=vector<double>(3,0)){
		_limitField[groupName]=LimitField(Inlet,-1,v_x,v_y,v_z,Temperature,-1,alpha,-1);
	};
	/** \fn setIntletPressureBoundaryCondition
			 * \brief adds a new boundary condition of type InletPressure
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the vapour volume fraction at the boundary
			 * \param [in] double : the value of the Pressure at the boundary
			 * \param [in] double : the value of the temperature at the boundary
			 * \param [out] void
			 *  */
	void setInletPressureBoundaryCondition(string groupName,double alpha,double Pressure,double Temperature){
		_limitField[groupName]=LimitField(InletPressure,Pressure,vector<double>(0,0),vector<double>(0,0),vector<double>(0,0),Temperature,-1,alpha,-1);
	};
	/** \fn setWallBoundaryCondition
			 * \brief adds a new boundary condition of type Wall
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the temperature at the boundary
			 * \param [in] vector<double> : the values of the x component of the 2 velocities at the boundary
			 * \param [in] vector<double> : the values of the y component of the 2 velocities at the boundary
			 * \param [in] vector<double> : the values of the z component of the 2 velocities at the boundary
			 * \param [out] void
			 *  */
	void setWallBoundaryCondition(string groupName,double Temperature,vector<double> v_x, vector<double> v_y=vector<double>(3,0), vector<double> v_z=vector<double>(3,0)){
		_limitField[groupName]=LimitField(Wall,-1,v_x,v_y,v_z,Temperature,-1,-1,-1);
	};

	/** \fn setIntPressCoeff
			 * \brief sets a value for the interfacial pressure default coefficient
			 * \details
			 * \param [in] double : the value for the interfacial pressure default coefficient
			 * \param [out] void
			 *  */
	void setIntPressCoeff(double delta){
		_intPressCoeff=delta;
	}

  protected :
	Field _Vitesse1,_Vitesse2;
	PetscScalar *_lCon, *_rCon;	// left and right conservative vectors
	PetscScalar * _JacoMat; //Jacobian matrix of the convection fluxes, used to compute the entropic corrections for the 5eqs two-fluid model
	PetscReal *_realPart, *_imagPart;
	double _intPressCoeff;
	//!calcule l'etat de Roe de deux etats
	void convectionState( const long &i, const long &j, const bool &IsBord);
	//!calcule la matrice de jacobienne de la convection de l'etat associé à une cellule
	void convectionJacobianMatrix(double *V, double *n);
	//!calcule la matrice de convection de l'etat interfacial entre deux cellules voisinnes
	void convectionMatrices();
	//!Calcule le flux pour un état et une porosité et une normale donnés
	Vector convectionFlux(Vector U,Vector V, Vector normale, double porosity);
	//!calcule la matrice de diffusion de l'etat interface pour la diffusion
	void diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord);
	//!Ajoute au second membre la contribution de la gravite, chgt phase, chauffage et frottement
	void sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i);
	//!Computes the pressure loss associated to the face ij
	void pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj);
	//!Computes the contribution of the porosity gradient associated to the face ij to the source term
	void porosityGradientSourceVector();
	//!Calcule la jacobienne de la CL convection
	void jacobian(const int &j, string nameOfGroup,double * normale);
	//!Calcule la jacobienne de la CL de diffusion
	void jacobianDiff(const int &j, string nameOfGroup);
	//!Calcule l'etat fictif à la frontiere
	void setBoundaryState(string nameOfGroup, const int &j,double *normale);
	//!Ajoute au second membre la contribution de la diffusion
	void addDiffusionToSecondMember(const int &i,const int &j,bool isBoundary);
	//!Computes the interfacial flux for the VFFC formulation of the staggered upwinding
	Vector staggeredVFFCFlux();
	//!Compute the corrected interfacial state for lowMach, pressureCorrection and staggered versions of the VFRoe formulation
	void applyVFRoeLowMachCorrections(bool isBord, string groupname="");

	//!Special preconditioner based on a matrix scaling strategy
	void computeScaling(double offset);
	//!Calcule les saut de valeurs propres pour la correction entropique
	void entropicShift(double* n);

	// Functions of equations of states
	void consToPrim(const double *Ucons, double* Vprim,double porosity=1);
	void primToCons(const double *V, const int &i, double *U, const int &j);
	void primToConsJacobianMatrix(double *V);


	double intPressDef(double alpha, double ur_n, double rho1, double rho2, double temperature);
	void entropicShift(double*n, double& vpcorr0, double& vpcorr1);
};

#endif /* FiveEqsTwoFluid_HXX_ */
