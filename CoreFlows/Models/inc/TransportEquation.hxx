//============================================================================
/**
 * \file TransportEquation.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date 24 March 2015
 * \brief Fluid enthalpy transport equation
 * */
//============================================================================

/*! \class TransportEquation TransportEquation.hxx "TransportEquation.hxx"
 *  \brief Scalar advection equation for a fluid enthalpy
 *  \details see \ref TransportEqPage for more details
 */
#ifndef TransportEquation_HXX_
#define TransportEquation_HXX_

#include "ProblemCoreFlows.hxx"

using namespace std;


//! enumeration phase
/*! The fluid type can be LiquidPhase or water  */
enum phase
{
	LiquidPhase,/**< Fluid considered is GasPhase */
	GasPhase/**< Fluid considered is Gas */
};

//! enumeration pressureEstimate
/*! the pressure estimate needed to fit physical parameters  */
enum pressureMagnitude
{
	around1bar300KTransport,/**< pressure is around 1 bar and temperature around 300K (for TransportEquation, SinglePhase and IsothermalTwoFluid) or 373 K (saturation for DriftModel and FiveEqsTwoFluid) */
	around155bars600KTransport/**< pressure is around 155 bars  and temperature around 618 K (saturation) */
};

//! enumeration BoundaryType
/*! Boundary condition type  */
enum BoundaryTypeTransport	{InletTransport,  OutletTransport, NeumannTransport, DirichletTransport, NoneBCTransport};//Actually Inlet=Dirichlet and Outlet=Neumann

/** \struct LimitField
 * \brief value of some fields on the boundary  */
struct LimitFieldTransport{
	LimitFieldTransport(){bcType=NoneBCTransport; T=0; h=0; flux=0; }
	LimitFieldTransport(BoundaryTypeTransport _bcType, double _T,	double _h,double _flux	){
		bcType=_bcType; T=_T; h=_h; flux=_flux;
	}

	BoundaryTypeTransport bcType;
	double T; //for inlet or Dirichlet
	double h; //for inlet or Dirichlet
	double flux; //for Neumann or outlet
};

class TransportEquation: public ProblemCoreFlows
{

public :
	/** \fn TransportEquation
			 * \brief Constructor for the enthalpy transport in a fluid
			 * \param [in] phase : \ref Liquid or \ref Gas
			 * \param [in] pressureMagnitude : \ref around1bar or \ref around155bars
			 * \param [in] vector<double> : fluid velocity (assumed constant)
			 *  */
	TransportEquation(phase fluid, pressureMagnitude pEstimate,vector<double> vitesseTransport);

	//Gestion du calcul
	virtual void initialize();
	virtual void terminate();//vide la mémoire et enregistre le résultat final
	bool initTimeStep(double dt);
	double computeTimeStep(bool & stop);//propose un pas de temps pour le calcul. Celà nécessite de discrétiser les opérateur (convection, diffusion, sources) et pour chacun d'employer la condition cfl. En cas de problème durant ce calcul (exemple t=tmax), renvoie stop=true
	void abortTimeStep();//efface les inconnues calculées par solveTimeStep() et reinitialise dt à 0
	bool iterateTimeStep(bool &ok);
	virtual void save();
	virtual void validateTimeStep();

	/* Boundary conditions */
	/** \fn setIntletBoundaryCondition
			 * \brief adds a new boundary condition of type Inlet
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the temperature at the boundary
			 * \param [out] void
			 *  */
	void setInletBoundaryCondition(string groupName,double enthalpy){
		_limitField[groupName]=LimitFieldTransport(InletTransport,-1,enthalpy,-1);
	};

	/** \fn setNeumannBoundaryCondition
	 * \brief adds a new boundary condition of type Neumann
	 * \details
	 * \param [in] string the name of the boundary
	 * \param [out] void
	 *  */
	void setNeumannBoundaryCondition(string groupName, double flux=0){
		_limitField[groupName]=LimitFieldTransport(NeumannTransport,-1,flux,-1);
	};

	/** \fn setBoundaryFields
	 * \brief met à jour  _limitField  ( le type de condition limite )
	 * \details
	 * \param [in] string
	 * \param [out] void
	 *  */
	void setBoundaryFields(map<string, LimitFieldTransport> boundaryFields){
		_limitField = boundaryFields;
	};


	/*Physical parameters*/
	void setLiqSatEnthalpy(double hsatl){
		_hsatl=hsatl;
	};
	void setVapSatEnthalpy(double hsatv){
		_hsatv=hsatv;
	};
	void setLiqSatDensity(double rhosatl){
		_rhosatl=rhosatl;
	};
	void setVapSatDensity(double rhosatv){
		_rhosatv=rhosatv;
	};
	void setTransportVelocity(Vector v){
		_vitesseTransport=v;
	};

	//get output fields for postprocessing or coupling
	vector<string> getOutputFieldsNames() ;//liste tous les champs que peut fournir le code pour le postraitement
	Field&         getOutputField(const string& nameField );//Renvoie un champs pour le postraitement

	Field& getFluidTemperatureField(){
		return _TT;
	}

	Field& getEnthalpyField(){
		return _VV;
	}

	/** \fn getTimeScheme
	 * \brief returns the  time scheme name
	 * \param [in] void
	 * \param [out] enum TimeScheme (explicit or implicit)
	 *  */
	TimeScheme getTimeScheme();

protected :
	double computeTransportMatrix();
	double computeRHS();
	void updatePrimitives();
	double temperature(double h){
		return _Tref+(h-_href)/_cpref;
	};
	double voidFraction(double h){
		double titre=(h-_href)/(_hsatv-_hsatl);
		if (titre<0)
			return 0;
		else if (titre>1)
			return 1;
		else return titre*_rhosatl/(titre*_rhosatl+(1-titre)*_rhosatv);
	};
	double density(double alpha){
		return alpha*_rhosatv+(1-alpha)*_rhosatl;
	};

	Field   _TT, _Alpha, _Rho;//Fields of temperature and coupled temperature
	double _rhosatv, _rhosatl;
	double _Tref, _href, _cpref;
	Vector _vitesseTransport, _normale;
	bool _transportMatrixSet;
	Vec _Hn, _deltaH, _Hk, _Hkm1, _b0;
	double _dt_transport, _dt_src;

	map<string, LimitFieldTransport> _limitField;
};

#endif /* TransportEquation_HXX_ */
