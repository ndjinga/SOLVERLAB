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

#include "WaveStaggered.hxx"
#include "StiffenedGas.hxx"
#include "BarotropicLaw.hxx"
#include "Node.hxx"
#include "utilitaire_algebre.h"
#include "Mesh.hxx"

//! enumeration phaseType
/*! The material phase can be Gas or liquid  */
enum phaseTypeStaggered
{
    LiquidStaggered,
    GasStaggered
};

class EulerBarotropicStaggered : public WaveStaggered{
public :
	/** \fn EulerBarotropicStaggered
	 * \param [in] phaseType : \ref Liquid or \ref Gas
	 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
	 * \param [in] int : mesh dimension
	 *  */
	EulerBarotropicStaggered(phaseTypeStaggered fluid, pressureEstimate pEstimate, double a, double gamma, int dim);


	//! system initialisation
	void initialize();

	/** \fn terminate
     * \brief empties the memory
     * @param void
     *  */
    void terminate();
    double computeTimeStep(bool & stop);

    BarotropicLaw getBarotropicEOS(int i){
        BarotropicLaw * result = dynamic_cast<BarotropicLaw*>(_fluides[i]); 
        if(result)
             return *result;
        else
            throw CdmathException("Fluid EOS is not barotropic");
    }

    void AssembleMetricsMatrices();
    void UpdateDualDensity();
    bool iterateTimeStep(bool &converged);
    std::vector<double>  getTimeEvol();
    void Rhomax_Umax_Cmax();

    //********* Raviart-Thomas related functions ***********//
    std::vector<double> ReferenceBasisFunctionRaviartThomas(const int &i, const Point &Xhat, const std::vector<Node> &K_Nodes );
    std::vector<double>  Gradient_ReferenceBasisFunctionRaviartThomas(int i, const std::vector<Node> &K_Nodes );
    Point xToxhat(const Cell &K, const  Point &X, const std::vector<Node> & K_Nodes); 
    std::vector<double>  JacobianTransfor_K_X(const Point &X, const std::vector<Node> &K_Nodes);

    //In order to find the corresponding basis function on the reference element we test its image by Piola transform 
    //and select the only one that is non-zero when taken against n_sigma. 
    bool FindlocalBasis(const int &m,const Face &Facej, const int &j,const  Cell& K, const std::vector<Node> &K_Nodes );
    double MassLumping(const Cell &K, const int &idcell, const Face & Facej, const int &j);

    // K is the cell on which we evaluate the basis function
    // Support is the table containing the support (so only K when we compute in the cell)
    // Facej is the face of the basis function and j its number
    // X the point inwhich it is evaluated
    std::vector<double> PhysicalBasisFunctionRaviartThomas(Cell K, int idcell, std::vector<Cell> Support, Face Facej,int j, Point X);
    std::vector<double> Gradient_PhysicalBasisFunctionRaviartThomas(Cell K, int idcell,  std::vector<Cell> Support, Face Facej, int j, Point X);
    std::vector<double> VelocityRaviartThomas_at_point_X(Cell K, int idcell, Point X);

    // operation on matrices
    double Jacobian(const std::vector<double> & mat);
    std::vector<double> TensorProduct(std::vector<double> &u, std::vector<double> &v); //returns u tenso v
    double Contraction(std::vector<double> &u, std::vector<double> &v); // returns contraction of two order 2 tensors
    std::vector<double> InvTranspose(std::vector<double> &u); // returns (u^{-1})^t
     
protected :
 /** Fluid equation of state **/
    vector<    Fluide* > _fluides;//
	BarotropicLaw   *_compressibleFluid;
	double _Tref, _Pref; //EOS reference temperature &pressure

    PetscReal _rhoMax, _uMax;
    Vec _DualDensity,_Conv ;
	Mat _InvVolPrim, _InvVolDual, _DivRhoU, _LaplacianVelocity, _InvDualDensity  ;
	double _c;
	std::vector<double> _Entropy, _Time;
    std::map< int , std::map< int, std::vector< std::pair<std::vector<double>, std::vector<double> > >  >> _PhysicalPsif,_GradientPhysicalPsif;
    //  map< idcell, map< idface, pair< K, [(x_0, nabla Psi_sigma_|K (x_0) ), ..., (x_f, nabla Psi_sigma_|K (x_f) )] >>>
    

};
#endif /* EULERBAROTROPICSTAGGERED_HXX_*/
