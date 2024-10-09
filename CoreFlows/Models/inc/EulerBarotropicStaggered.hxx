//============================================================================
/**
 * \file EulerBarotropicStaggered.hxx
 * \author Esteban COIFFIER
 * \version 1.0
 * \date 01 janv. 2018
 * \briefThe Wave system
 * */
//============================================================================

/*! \class EulerBarotropicStaggered EulerBarotropicStaggered.hxx "EulerBarotropicStaggered.hxx"
 *  \brief The Wave system
 *  \details Wave system
 */
#ifndef EULERBAROTROPICSTAGGERED_HXX_
#define EULERBAROTROPICSTAGGERED_HXX_

#include "ProblemCoreFlows.hxx"
#include "StiffenedGas.hxx"
#include "Node.hxx"
#include "utilitaire_algebre.h"
#include "Mesh.hxx"

//! enumeration phaseType
/*! The material phase can be Gas or liquid  */
enum phaseType
{
    Liquid,/**< Material considered is Liquid */
    Gas/**< Material considered is Gas */
};

class EulerBarotropicStaggered : public ProblemCoreFlows{
public :
	/** \fn EulerBarotropicStaggered
	 * \param [in] phaseType : \ref Liquid or \ref Gas
	 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
	 * \param [in] int : mesh dimension
	 *  */
	EulerBarotropicStaggered(phaseType fluid, pressureEstimate pEstimate, int dim);

	void setInitialField(const Field &field);

	//! system initialisation
	void initialize();

	/** \fn terminate
     * \brief empties the memory
     * @param void
     *  */
    void terminate();

	double getTimeStep();

	void save();

	/** \fn computeTimeStep
     * \brief Proposes a value for the next time step to be solved using mesh data and cfl coefficient
     *  \return  double dt the proposed time step
     *  \return  bool stop, true if the calculation should not be continued (stationary state, maximum time or time step numer reached)
     *  */
    double computeTimeStep(bool & stop);

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

	 /** \fn validateTimeStep
     * \brief Validates the solution computed y solveTimeStep
     * \details updates the currens time t=t+dt, save unknown fields, resets the time step dt to 0, tests the stationnarity.
     * c It is a pure virtual function overloaded in each model
     * @param  void
     *  */
    void validateTimeStep();

	 /** \fn saveDensity
     * \brief saves the Density field in a separate file 
     * @param bool
     * */
    void saveDensity(bool save_p=true){
        _saveDensity=save_p;
    }

    /** \fn saveVelocity
     * \brief saves the velocity field in a separate 3D file so that paraview can display the streamlines
     * @param bool
     * */
    void saveVelocity(bool save_v=true){
        _saveVelocity=save_v;
    }

	/** \fn getStiffenedGasEOS
     * \brief return the stiffened gas law associated to fluid i
     * @param int i : the index of the fluid
     * @return throws an exception if the fluid with index i does not follow a stiffened gas law.
     * */
    StiffenedGas getStiffenedGasEOS(int i)
    {
        StiffenedGas * result = dynamic_cast<StiffenedGas*>(_fluides[i]); 
        if(result)
             return *result;
        else
            throw CdmathException("ProblemFluid::getStiffenedGasEOS() : fluid EOS is not a stiffened gas law");
    }

	void setVerticalPeriodicFaces();
	void setHorizontalPeriodicFaces();

	double getOrientation(int j, Cell Cint);
	void  setOrientation(int j,std::vector<double> vec_normal_sigma);
	void setWallBoundIndex(int j );

	std::map<int,double>  getboundaryDensity() const;
	void  setboundaryDensity(map< int, double> BoundaryDensity);
	void  setboundaryVelocity(map< int, double> BoundaryVelocity);

	void  abortTimeStep();
	bool  initTimeStep( double dt);
	vector<string> getInputFieldsNames();
	void setInputField(const string& nameField, Field& inputField );

	void InterpolateFromFacesToCells(const Field &atFaces, Field &atCells);


protected :
 /** Fluid equation of state **/
    vector<    Fluide* > _fluides;//
	CompressibleFluid *_compressibleFluid;
	double _Tref; //EOS reference temperature
    double _Pref; //EOS reference pressure

	Field _Velocity, _Density, _Velocity_at_Cells;
	Vec _newtonVariation, _primitiveVars,  _BoundaryTerms, _primitiveVars_seq;
	Mat _InvVol,_InvSurface,  _Div,  _LaplacianDensity,_Conv, _MinusGrad,_LaplacianVelocity  ;
	double _c, _maxPerim, _minCell ;
	double *_vec_normal;
	bool _saveDensity, _saveVelocity;
	std::map<int, double>  _boundaryDensity;
	std::map<int, std::vector<double> > _vec_sigma; // arbitrary degree of liberty associated to a face

	std::map<int,int> _indexFacePeriodicMap;
	std::vector<int>_indexWallBoundFaceSet; // map of perdiodic faces couples : only it->first is computed. it->second is avoided in the loop for matrices and is updated to it->first in save()
	bool _facesBoundinit,_indexFacePeriodicSet, _isWall; // To ensure that the boundary velocity is initialized after the initial velocity 
	std::vector<double> _Entropy;
				

};
#endif /* EULERBAROTROPICSTAGGERED_HXX_*/
