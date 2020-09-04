//============================================================================
/**
 * \file DriftModel.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date 01 janv. 2018
 * \brief The compressible Navier-Stokes equations with an ICE scheme on staggered meshes
 * */
//============================================================================

/*! \class SinglePhase SinglePhase.hxx "SinglePhase.hxx"
 *  \brief The compressible Navier-Stokes equations
 *  \details The model consists in one mass, one momentum and one energy equation, see \ref NSModelsPage for more details
 */
#ifndef SINGLEPHASESTAGGERED_HXX_
#define SINGLEPHASESTAGGERED_HXX_

#include "ProblemFluid.hxx"

class SinglePhaseStaggered : public ProblemFluid{
public :
	/** \fn SinglePhaseStaggered
	 * \brief Constructor for the Navier-Stokes system
	 * \param [in] phaseType : \ref Liquid or \ref Gas
	 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
	 * \param [in] int : mesh dimension
	 * \param [in] bool : There are two possible equations of state for the fluid
	 *  */
	SinglePhaseStaggered(phaseType fluid, pressureEstimate pEstimate,int dim,bool useDellacherieEOS=false);
	//! system initialisation
	void initialize();

	//fonctions d'echange de flux
	//	void getOutputField(const Vec &Flux, const string Champ, const int numBord)=0;//, PetscInt *indices_Flux, PetscInt *indices_Bord, const long range)=0;
	//	double trace(const int &numBord, Vec &out)=0;

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
	/** \fn setIntletPressureBoundaryCondition
	 * \brief adds a new boundary condition of type InletPressure
	 * \details
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the pressure at the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [out] void
	 *  */
	void setInletPressureBoundaryCondition(string groupName, double pressure,double Temperature){
		_limitField[groupName]=LimitField(InletPressure,pressure,vector<double>(0,0),vector<double>(0,0),vector<double>(0,0),Temperature,-1,-1,-1);
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
	void setWallBoundaryCondition(string groupName,double Temperature,double v_x, double v_y=0, double v_z=0){
		_limitField[groupName]=LimitField(Wall,-1,vector<double>(1,v_x),vector<double>(1,v_y),vector<double>(1,v_z),Temperature,-1,-1,-1);
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

	void computeVelocityMCells(const Field& velocity,
                      Field& velocityMCells)
protected :
	Field _Vitesse;
	double  _drho_sur_dp,   _drho_sur_dT;//derivatives of the density rho wrt cv, p, T
	double  _drhoE_sur_dp,  _drhoE_sur_dT;//derivatives of the total energy rho E wrt cv, p, T
	bool _useDellacherieEOS;

	//!calcule l'etat de Roe de deux etats
	void convectionState( const long &i, const long &j, const bool &IsBord);
	//!calcule la matrice de convection de l'etat interfacial entre deux cellules voisinnes
	void convectionMatrices();
	//!Calcule le flux pour un état et une porosité et une normale donnés
	Vector convectionFlux(Vector U,Vector V, Vector normale, double porosity);
	//!Computes the source vector associated to the cell i
	void sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i);
	//!Computes the pressure loss associated to the face ij
	void pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj);
	//!Computes the contribution of the porosity gradient associated to the face ij to the source term
	void porosityGradientSourceVector();
	//!Calcule la jacobienne de la CL convection
	void jacobian(const int &j, string nameOfGroup,double * normale);
	//!Calcule l'etat fictif a la frontiere
	void setBoundaryState(string nameOfGroup, const int &j,double *normale);// delete &nf Kieu
	//!Adds the contribution of diffusion to the RHS
	void convectionMatrixPrimitiveVariables( double rho, double u_n, double H,Vector velocity);
	/** \fn getDensityDerivatives
	 * \brief Computes the partial derivatives of rho, and rho E with regard to the primitive variables  p and  T
	 * @param pressure
	 * @param temperature
	 * @param square of the velocity vector
	*/
	void getDensityDerivatives( double pressure, double temperature, double v2);

};
#endif /* SINGLEPHASESTAGGERED_HXX_*/
