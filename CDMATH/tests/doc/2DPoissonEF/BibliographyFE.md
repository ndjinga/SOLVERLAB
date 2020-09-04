- In 2D the minimum angle condition (1968)  
  $\alpha_K\geq\alpha_0>0,\,\forall K\in \mathcal{T}$ where $\alpha_K$ is the minimal angle of $K$  
  is sufficient for the convergence of the 2D linear finite element method
  $$
  \lim_{h\to 0} || u - u_h||_{H^1} =0.
  $$ 
  Moreover, if the exact solution $u\in H^{k+1}$, then
  $$
  \exists C,\forall h,\quad || u - u_h||_{H^1}\geq C h^k || u ||_{H^{k+1}}
  $$
  See theorem 6.3.13 and remark 6.3.12 in  
  *Grégoire Allaire, Numerical analysis and optimization, Oxford University Press, 2007*

- In dimension greater than 2, the minimum angle condition can be replaced by the condition
  $$
  \exists C,\forall h,\forall K\in\mathcal{T}\quad diam(K)\leq C \rho(K)
  $$
  where $diam(K)$ is the diameter of the cell $K$ and $\rho(K)$ is the diameter of the largest inscribed circle.  
  See definition 6.3.11 in  
  *Grégoire Allaire, Numerical analysis and optimization, Oxford University Press, 2007*  
  *J. Brandts, S. Korotov, M. Křı́žek, On the equivalence of ball conditions for simplicial finite elements in $\mathbb{R}^d$ , Appl. Math. Lett. 22 (2009), 1210–1212*

- The maximum angle condition (1976)  
  $\gamma_K\leq\gamma_0<\pi,\,\forall K\in \mathcal{T}$ where $\gamma_K$ is the maximal angle of $K$  
  is sufficient for the convergence of the quadratic finite element method
  $$
  \exists C,\quad || u - u_h||_1 \leq C h ||u||_2
  $$ 
  as $h\to 0$.

- The maximum angle condition is not necessary for convergence  
  *A. Hannukainen, S. Korotov, M. Křı́žek, Maximum angle condition is not necessary for convergence of the finite element method,  M. Numer. Math. (2012) 120: 79*

- Acute triangles (angles less than $\pi/2$) yield a maximum principle thanks to a diagonally dominant stiffness matrices  
  See exercice 6.3.13 in
  *Grégoire Allaire, Numerical analysis and optimization, Oxford University Press, 2007*  
  *J. Karátson, S. Korotov, M. Křı́žek, On discrete maximum principles for nonlinear elliptic problems, Math. Comput. Simulation 76 (2007), 99–108*  
  *J. Brandts, S. Korotov, M. Křı́žek, The discrete maximum principle for linear simplicial finite element approximations of a reaction-diffusion problem, Linear Algebra Appl. 429 (2008), 2344–2357*
