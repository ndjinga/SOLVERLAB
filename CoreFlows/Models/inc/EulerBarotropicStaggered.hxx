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
    void save();
    void terminate();
    double computeTimeStep(bool & stop);
    void computeNewtonVariation();

    BarotropicLaw getBarotropicEOS(int i){
        BarotropicLaw * result = dynamic_cast<BarotropicLaw*>(_fluides[i]); 
        if(result)
             return *result;
        else
            throw CdmathException("Fluid EOS is not barotropic");
    }

    void UpdateDualDensity();
    //void computeNewtonVariation(); 
    bool iterateTimeStep(bool &converged);
    std::vector<double>  getTimeEvol();
    void Rhomax_Umax_Cmax();

    //********* Raviart-Thomas related functions ***********//
    std::vector<double> HessienneTransfo(const int component, const std::vector<Node> &K_Nodes);
    std::vector<double> Gradient_ReferenceBasisFunctionRaviartThomas(int i, const std::vector<Node> &K_Nodes, const Point & Xhat );
    std::vector<double> Gradient_PhysicalBasisFunctionRaviartThomas(Cell K, int idcell,  const std::vector<Cell>& Support, Face Facej, int j, Point X);
    std::vector<double> MomentumRaviartThomas_at_point_X(Cell K, int idcell, Point X);

    // operation on matrices
    std::vector<double> TensorProduct(std::vector<double> &u, std::vector<double> &v); //returns u tenso v
    double Contraction(std::vector<double> &u, std::vector<double> &v); // returns contraction of two order 2 tensors
    std::vector<double> Inverse(std::vector<double> &u); // returns (u^{-1})^t

    double getOrientationNode(int n, int j); //n is a node, j a face, gives back sign( (x_n - x_j). n_sigma^perp    )

    std::map<int, std::vector<double> >  getboundaryVelocityVector() const;
	void setboundaryVelocityVector(int j,  std::vector<double>  boundaryVelocityVector);
    std::vector<double> H_1DensitySemi_Norm__H_divVelocitySemi_Norm();
    void computeOrder2Density(const double& rho_b, const double& Mach );
     
protected :
    /** Fluid equation of state **/
    vector<    Fluide* > _fluides;//
	BarotropicLaw   *_compressibleFluid;
    Field _MachNumber;
	
    Vec _DualDensity, _Conv, _GradPressure ;
    Mat _JacobianMatrix;
	double _c, _rhoMax, _uMax;
	std::vector<double> _Entropy, _Time;
    std::map< int , std::map< int, std::vector< std::pair<std::vector<double>, std::vector<double> > >  >> _PhysicalPsif,_GradientPhysicalPsif;
    std::map< int ,std::vector<double> > _boundaryVelocityVector;
    //  map< idcell, map< idface, pair< K, [(x_0, nabla Psi_sigma_|K (x_0) ), ..., (x_f, nabla Psi_sigma_|K (x_f) )] >>>
    

};
#endif /* EULERBAROTROPICSTAGGERED_HXX_*/
