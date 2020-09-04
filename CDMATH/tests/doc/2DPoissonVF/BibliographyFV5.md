## Some bibliographical remarks about the two points finite volume scheme

- Order 1 convergence on orthogonal meshes : neighbouring cells  $C_i$ and $C_j$ must be separated by a face (or edge in 2D) $f_{ij}$ that is perpendicular to the straight line connecting the center of masses $x_i$ of $C_i$ and $x_j$ of $C_j$  
  *R. Eymard, T. Gallouët, R. Herbin, Finite Volume Methods, Handbook for Numerical Analysis, Ph. Ciarlet, J.L. Lions eds, North Holland, 2000, 715-1022. 

- Order 1 convergence on not too deformed triangular meshes : the triangles edges must be in $O(h)$ and the triangle areas must be in $O(h^2)$ (angles must not shrink to $0{}^\circ$ nor $180{}^\circ$)   
  *R. Herbin, An error estimate for a four point finite volume scheme for the convection-diffusion equation on a triangular mesh, Num. Meth. P.D.E., 165-173, 1995.*


- Order 2 convergence on triangular meshes, provided 
    - the center of the circumscribed circle is used instead of the center of mass in each cell for the evaluation of the source term and analytical solution
    - the Delaunay conditions are satisfied (no neighboring cell is included in the circumscribed circle of an arbitrary cell)


- Non convergence on flat degenerating triangular meshes  
  *K. Domelevo, P. Omnes, A finite volume method for the Laplace equation on almost arbitrary 2D grids, Mathematical Modelling and Numerical Analysis, 2005*


- Order 1 if the mesh is conforming except on a line  
  *J. Droniou, C. Le Potier, Construction and Convergence Study of Schemes Preserving the Elliptic Local Maximum Principle, SIAM Journal on Numerical Analysis, 2011*


- Order 2 on triangular meshes provided 1) Delaunay type conditions are satisfied and 2) $f\in H^1$ and meshes are generated from an initial mesh either by subdivisions,symmetry or translation  
  *J. Droniou, Improved L^2 estimate for gradient schemes and super-convergence of the TPFA finite volume scheme, IMA Journal of Numerical Analysis, 2018*


- It is possible to converge with order 1 on the gradient, but only order 1 on the function ie there is no equivalent of the Aubin-Nitsche lemma in the finite volume context  
  *P. Omnes, Error estimates for a finite volume method for the Laplace equation in dimension one through discrete Green functions. International Journal on Finite Volumes 6(1), 18p., electronic only, 2009*
