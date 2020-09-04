#ifndef utilitaire_algebre_h
#define utilitaire_algebre_h

/*      rpoly.cpp -- Jenkins-Traub real polynomial root finder.
 *
 *      (C) 2000, C. Bond.  All rights reserved.
 *
 *      Translation of TOMS493 from FORTRAN to C. This
 *      implementation of Jenkins-Traub partially adapts
 *      the original code to a C environment by restruction
 *      many of the 'goto' controls to better fit a block
 *      structured form. It also eliminates the global memory
 *      allocation in favor of local, dynamic memory management.
 *
 *      The calling conventions are slightly modified to return
 *      the number of roots found as the function value.
 *
 *      INPUT:
 *      op - double precision vector of coefficients in order of
 *              decreasing powers.
 *      degree - integer degree of polynomial
 *
 *      OUTPUT:
 *      zeror,zeroi - output double precision vectors of the
 *              real and imaginary parts of the zeros.
 *
 *      RETURN:
 *      returnval:   -1 if leading coefficient is zero, otherwise
 *                  number of roots found. 
 */

/* A roots_polynoms class */

/*! \class roots_polynoms utilitaire_algebre.hxx "utilitaire_algebre.hxx"
 *  \brief  Computes the roots of a given polynomial
 *  \details  Translation of TOMS493 from FORTRAN to C. This
 *      implementation of Jenkins-Traub partially adapts
 *      the original code to a C environment by restruction
 *      many of the 'goto' controls to better fit a block
 *      structured form. It also eliminates the global memory
 *      allocation in favor of local, dynamic memory management. 
 *      The calling conventions are slightly modified to return
 *      the number of roots found as the function value.
 */
#include <math.h>
#include <vector>
#include <complex>
#include <iostream>

using namespace std;

class roots_polynoms 
{
public:

	roots_polynoms();

	/** \fn quad
	 * \brief calcule les racines d'un polynome de degré 2
	 * \Details Calculate the zeros of the quadratic a*x² + b*x + c.
	 *  The quadratic formula, modified to avoid overflow, is used
	 *  to find the larger zero if the zeros are real and both
	 *  are complex. The smaller real zero is found directly from
	 *  the product of the zeros c/a
	 * @param a double correspond au coefficient de x²
	 * @param b double correspond au coefficient de x
	 * @param c double correspond au coefficient de 1
	 * @param *sr double correspond à la pratie réelle de la premiere racine
	 * @param *si double correspond à la pratie imaginaire de la premiere racine
	 * @param *lr double correspond à la pratie réelle de la deuxième racine
	 * @param *li double correspond à la pratie imaginaire de la deuxième racine
	 * @return calcule les racine du polynome dans sr+i*si et lr+i*li
	 * */
	void quad(double a,double b,double c,double *sr,double *si, double *lr,double *li);


	/** \fn fxshfr
	 * \Details  Computes up to L2 fixed shift k-polynomials,
	 *  testing for convergence in the linear or quadratic
	 *  case. Initiates one of the variable shift
	 *  iterations and returns with the number of zeros found
	 * @param l2 entier
	 * @param *nz entier
	 * @return */
	void fxshfr(int l2, int *nz);

	/** \fn quadit
	 * \brief
	 * \Details /*  Variable-shift k-polynomial iteration for a
	 *  quadratic factor converges only if the zeros are
	 *  equimodular or nearly so.
	 * @param *uu, *vv double coefficients of starting quadratic
	 * @param *nz entier number of zeros found
	 * @return */
	void quadit(double *uu,double *vv,int *nz);

	/** \fn realit
	 * \brief Variable-shift H polynomial iteration for a real zero
	 * @param  *sss starting iterate (double)
	 * @param *nz number of zeros found (int)
	 * @param *iflag flag to indicate a pair of zeros near real axis (int)
	 * @return */
	void realit(double *sss, int *nz, int *iflag);

	/** \fn calcsc
	 * \Details This routine calculates scalar quantities used to
	 *  compute the next k polynomial and new estimates of
	 *  the quadratic coefficients.
	 * @param *type integer variable set here
	 * indicating how the calculations are normalized to avoid overflow
	 * @return */
	void calcsc(int *type);

	/** \fn nextk
	 * \brief Computes the next k polynomials using scalars computed in calcsc().
	 * @param *type integer variable set here
	 * indicating how the calculations are normalized to avoid overflow
	 * @return */
	void nextk(int *type);

	/** \fn newest
	 * \brief Compute new estimates of the quadratic coefficients
	 *  using the scalars computed in calcsc ().
	 * @param type
	 * @param *uu, *vv double coefficients of starting quadratic
	 * @return */
	void newest(int type,double *uu,double *vv);

	/** \fn quadsd
	 * \brief Divides p by the quadratic 1,u,v placing the quotient
	 *  in q and the remainder in a,b.
	 * @return */
	void quadsd(int n,double *u,double *v,double *p,double *q,
			double *a,double *b);

	/** \fn rpoly
	 * \brief
	 * @param *op
	 * @param degree
	 * @param *zeror
	 * @param *zeroi
	 * @return */
	int rpoly(double *op, int degree, double *zeror, double *zeroi);

private:
	double infin;
	double smalno;
	double eta;
	double base;
	double *p,*qp,*k,*qk,*svk;
	double sr,si,u,v,a,b,c,d,a1,a2;
	double a3,a6,a7,e,f,g,h,szr,szi,lzr,lzi;
	double are,mre;
	int n,nmi;
};

class Polynoms
{
public:

	void abs_par_interp_directe(int n,  vector< complex< double > > vp,   double * A, int sizeA, double epsilon, double *result,vector< complex< double > > y);

	/** \fn polynome_caracteristique
	 * \brief Calcule le polynome caracteristique
	 * @return */
	vector< double > polynome_caracteristique(double alpha_v, double alpha_l, double Uv, double Ul, double rhov, double rhol, double   cv_2, double cl_2, double dPiv, double dPil);//c= inverse vitesse du son!
	/** \fn polynome_caracteristique
	 * \brief
	 * @param alpha1
	 * @param alpha2
	 * @param u1
	 * @param u2
	 * @param rho1
	 * @param rho2
	 * @param invc1sq
	 * @param invc2sq
	 * @param dpi1
	 * @param dpi2
	 * @param g2press
	 * @param g2alpha
	 * @param g2
	 * @param epsilon
	 * @return */
	vector< double > polynome_caracteristique(double alpha1, double alpha2, double u1, double u2, double rho1, double rho2, double invc1sq, double invc2sq, double dpi1, double dpi2, double g2press, double g2alpha, double g2, double epsilon);

	/** \fn module
	 * \brief calcule le module d'un nombre complexe
	 * @param z est un nombre complexe
	 * @return  calcule le module de z*/
	double module(complex< double > z);
	/** \fn modulec calcule le module² de z */
	double modulec(complex< double > z);
	/** \fn abs_generalise
	 * \brief calcule la valeur absolue d'un nombre complexe
	 * \Details calcule la valeur absolue d'un nombre complexe en prenant en compte
	 * que la partie réelle c-a-d si z= a+ib abs_generalize() renvoie |a|+ib
	 * @param z
	 * @return si z = a+ib elle renvoie |a|+ib */
	complex< double > abs_generalise(complex< double > z);

	int new_tri_selectif(vector< complex< double > > &L, int n, double epsilon);

	/** \fn tri_selectif
	 * \brief
	 * @param L
	 * @param n
	 * @param epsilon
	 * @return */
	template< typename T >
	int tri_selectif(T& L, int n, double epsilon);

	/** \fn matrixProduct
	 * \brief calcule le produit de 2 matrices .
	 * @param *Matrix1 la 1ere matrice
	 * @param R1,C1 la taille de la matrice1 Nbr de ligne et le nbr de colonne
	 * @param Matrix2 la 2ieme matrice
	 * @param R2,C2 la taille de la matrice2 Nbr de ligne et le nbr de colonne
	 * @param out le resultat de Matrix1*Matrix2
	 * @return */
	void matrixProduct
	(
			const double *Matrix1,
			const int &R1,
			const int &C1,
			const double *Matrix2,
			const int &R2,
			const int &C2,
			double *out
	);

	/** \fn matrixProdVec
	 * \brief Calcule le produit Matrice vecteur
	 * @param *Matrix correspond à la matrice du produit
	 * @param R1,C1 la taille de la matrice( Nbr de ligne et le nbr de colonne)
	 * @param Vector correspond au vecteur du produit
	 * @param out le resultat de Matrix*Vector
	 * @return   */
	void matrixProdVec
	(	const double *Matrix,
			const int &R1,
			const int &C1,
			const double *Vector,
			double *out
	);

	bool belongTo(complex< double > a , vector< complex <double > > v, double epsilon);

private:

/** \fn add
 * \brief Calcule la somme de 2vecteurs ( U=U+V)
 * @param *U,*V les 2 vecteurs
 * @param N taille des 2 vecteurs
 * @return la somme de U et V dans U*/
void add
(
		double *U,		//vecteur auquel on ajoute
		int N,			//taille du vecteur
		const double *V	//ajout
);
/** \fn tensor
 * \brief Calcule le tenseur de 2 vecteurs
 * @param *a le 1er vecteur (*double)
 * @param N la taille du premier vecteur (int)
 * @param *b le 2ieme vecteur (*double)
 * @param M la taille du 2ieme vecteur (int)
 * @param *ab le resultat
 * @return */

void tensor
(	const double *a,	//vecteur gauche
		int N,			//taille du vecteur gauche
		const double *b,	//vecteur droit
		int M,			//taille du vecteur droit
		double *ab		//conteneur de sortie
);

/** \fn shift_diagonal
 * \brief Calcule la trace d'une matrice carrée
 * @param *Matrix correspond  à la matrice
 * @param N taille de la matrice
 * @param shift résultat
 * @return renvoie la trace de la matrice Matrix dans shift */
void shift_diagonal( double *Matrix, const int N, double shift);
//qques outils d'algebre des polynomes

/** \fn remplir_trinome
 * \brief remplie un trinome
 */
void remplir_trinome(double coeff, double u, double alpha, vector< double >& trinome);

/** \fn additionne
 * \brief Additionne 2 vecteurs de tailles differentes
 * @param P,Q les 2 vecteurs à additionner
 * @param n,m les 2 tailles respectives des vecteurs P et Q
 * @return un vecteur qui correspond a P+Q*/
vector< double > additionne(const vector< double >& P, const vector< double >& Q, int n, int m);

/** \fn multiplie
 * \brief Calcule le produit scalaire de 2 vecteurs de tailles differentes
 * @param P,Q les 2 vecteurs
 * @param n,m les 2 tailles respectives des 2 vecteurs
 * @return un vecteur correspond à (P,Q) */
vector< double > multiplie(const vector< double >& P, const vector< double >& Q, int n, int m);

//Pour le tri des valeurs propres

void ordre_croissant_abs(vector< complex< double > > &L, int n);

//Calcul des coefficients du polynome d'interpolation x->y dans la base de Newton
template<class T>

/** \fn dif_div
 * \brief
 * @param n
 * @param x
 * @param y
 * @param epsilon
 * @return */
vector< complex< double > > dif_div(int n, const vector< complex< double > >& x, const vector< T >& y, double epsilon);//n=nb valeurs à interpoler

//attention, Les vp complexes conjugu�es doivent se suivre dans la liste x
void appliquer_dif_div(int n, const vector< complex< double > >& dif_div, const vector< complex< double > >& x, const double* A, const int sizeA, double epsilon, double *p);//p=result

double avr(double a, double b);

};
#endif 
