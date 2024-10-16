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


//! enumeration BoundaryType
/*! Boundary condition type  */
enum BoundaryTypeTransport    {InletTransport,  OutletTransport, NeumannTransport, DirichletTransport, NoneBCTransport};//Actually Inlet=Dirichlet and Outlet=Neumann

//! enumeration phaseType
/*! The material phase can be Solid, Gas or liquid  */
enum FluidMaterial
{
    Air,/**< Material considered is air */
    Water/**< Material considered is water */
};

/** \struct LimitField
 * \brief value of some fields on the boundary  */
struct LimitFieldTransport{
    LimitFieldTransport(){bcType=NoneBCTransport; T=0; h=0; flux=0; }
    LimitFieldTransport(BoundaryTypeTransport _bcType, double _T,    double _h,double _flux    ){
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
             * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
             * \param [in] vector<double> : fluid velocity (assumed constant)
             *  */
    TransportEquation(FluidMaterial fluid, pressureEstimate pEstimate,vector<double> vitesseTransport, MPI_Comm comm = MPI_COMM_WORLD);

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
             * \details same as setDirichletBoundaryCondition
             * \param [in] string : the name of the boundary
             * \param [in] double : the value of the enthalpy at the boundary
             * \param [out] void
             *  */
    void setInletBoundaryCondition(string groupName,double enthalpy){
        _limitField[groupName]=LimitFieldTransport(InletTransport,-1,enthalpy,-1);
    };
    /** \fn setDirichletBoundaryCondition 
             * \brief adds a new boundary condition of type Dirichlet
             * \details same as setInletBoundaryCondition
             * \param [in] string : the name of the boundary
             * \param [in] double : the value of the enthalpy at the boundary
             * \param [out] void
             *  */
    void setDirichletBoundaryCondition(string groupName,double enthalpy){
        _limitField[groupName]=LimitFieldTransport(DirichletTransport,-1,enthalpy,-1);
    };
    /** \fn setDirichletBoundaryCondition
             * \brief adds a new boundary condition of type Inlet
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
        _limitField[groupName]=LimitFieldTransport(DirichletTransport, -1, 0, -1);
    };

    /** \fn setNeumannBoundaryCondition 
     * \brief adds a new boundary condition of type Neumann
     * \details same as setOutletBoundaryCondition
     * \param [in] string the name of the boundary
     * \param [in] double : the value of the enthalpy flux at the boundary
     * \param [out] void
     *  */
    void setNeumannBoundaryCondition(string groupName, double flux=0){
        _limitField[groupName]=LimitFieldTransport(NeumannTransport,-1,-1,flux);
    };
    /** \fn setOutletBoundaryCondition 
     * \brief adds a new boundary condition of type Outlet
     * \details same as setNeumannBoundaryCondition
     * \param [in] string the name of the boundary
     * \param [in] double : the value of the enthalpy flux at the boundary
     * \param [out] void
     *  */
    void setOutletBoundaryCondition(string groupName, double flux=0){
        _limitField[groupName]=LimitFieldTransport(OutletTransport,-1,-1,flux);
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
        _limitField[groupName]=LimitFieldTransport(NeumannTransport,-1,-1, 0);
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

    /* set input fields to prepare the simulation */
    vector<string> getInputFieldsNames();
    void setInputField(const string& nameField, Field& inputField );//supply of a required input field
    
    /** \fn setRodTemperatureField
     * \brief Set the rod temperature field
     * \details
     * \param [in] Field
     * \param [out] void
     *  */
    void setRodTemperatureField(Field rodTemperature);

    /** \fn setRodTemperature 
     * \brief Set a constant rod temperature field
     * \details
     * \param [in] double
     * \param [out] void
     *  */
    void setRodTemperature(double rodTemp){
        _rodTemperature=rodTemp;
        _isStationary=false;//Source term may be changed after previously reaching a stationary state
    }

    /** \fn getRodTemperatureField
     * \brief
     * \details
     * \param [in] void
     * \param [out] Field
     *  */
    Field& getRodTemperatureField(){ // ?? je ne retrouve pas cet attribut dans le file.cxx
        return _rodTemperatureField;
    }

    /* get output fields for postprocessing or coupling */
    vector<string> getOutputFieldsNames() ;//liste tous les champs que peut fournir le code pour le postraitement
    Field&         getOutputField(const string& nameField );//Renvoie un champs pour le postraitement

    Field& getFluidTemperatureField(){
        return _TT;
    }

    Field& getEnthalpyField(){
        return _VV;
    }

    Field& getVoidFractionField(){
        return _Alpha;
    }

    Field& getDensityField(){
        return _Rho;
    }

protected :
    double computeTransportMatrix();
    double computeRHS();
    void updatePrimitives();

    /* Postprocessing fields */
    Field   _TT, _Alpha, _Rho;//Fields of temperature, void fraction, density. Unknown field is enthalpy (_VV)
    double _rhosatv, _rhosatl;
    double _Tref, _href, _cpref;

    double temperature(double h){
        return _Tref+(h-_href)/_cpref;
    };
    double voidFraction(double h){
        double titre=(h-_hsatl)/(_hsatv-_hsatl);
        if (titre<0)
            return 0;
        else if (titre>1)
            return 1;
        else return titre*_rhosatl/(titre*_rhosatl+(1-titre)*_rhosatv);
    };
    double density(double alpha){
        return alpha*_rhosatv+(1-alpha)*_rhosatl;
    };

    Vector _vitesseTransport, _normale;
    bool _transportMatrixSet;
    Vec _Hn, _deltaH, _Hk, _Hkm1, _b0;
    Vec _Hn_seq; // Local sequential copy of the parallel vector _Hn, used for saving result files
    double _dt_transport, _dt_src;

    map<string, LimitFieldTransport> _limitField;
    
    /* source terms */
    bool   _rodTemperatureFieldSet;
    Field  _rodTemperatureField;
    double _rodTemperature;
};

#endif /* TransportEquation_HXX_ */
