The numerical methods 
=====================

CoreFlows proposes a variety of finite volume methods (see \ref leveque for an introduction). The method  can be explicit or implicit (enum \ref TimeScheme), and the convection fluxes can be approximated using the  \ref roe,  \ref vfroe or \ref vffc formulations (enum \ref NonLinearFormulation).

The basic method for fluid models is the \ref roe scheme with entropic correction (see \ref kieu) and source upwinding (see\ref wbsourceterms). In order to increase precision it is possible to use centered, staggered, pressureCorrection or lowMach schemes through the enum \ref SpaceScheme.

The finite volume discretisation allows an easy handling of general geometries and meshes generated by \ref salome .

Explicit schemes are used in general for fast dynamics solved with small time steps while implicit schemes allow the use of large time step to quickly reach the stationary regime. The implicit schemes result in nonlinear systems that are solved using a Newton type method.

The upwind scheme is the basic scheme but options are available to use a centered scheme (second order in space) or entropic corrections.

Our models can be written in generic form as a nonlinear system of balance laws:

$$
\frac{\partial U}{\partial t} + \nabla \cdot \left(\mathcal{F}^{conv}(U)\right)+\nabla \cdot \left(\mathcal{F}^{diff}(U)\right) = S(U,x), 
$$

where 
- $U$ is the vector of conservative unknowns, 
- $\mathcal{F}^{conv}$ is the convective flux 
- and $\mathcal{F}^{diff}$ the diffusive flux.

We decompose the computational domain into $N$ disjoint cells $C_i$ with volume $v_i$.
+ Two neighboring cells $C_i$ 
  and $C_j$ have a common boundary 
  $\partial C_{ij}$ 
  with area $s_{ij}$.
+ We denote $N(i)$ the set of neighbors 
  of a given cell $C_i$ 
  and $\vec n_{ij}$ the exterior unit normal vector 
  of $\partial C_{ij}$ . 

Integrating the system (\ref NStokesEq) over $C_{i}$ and setting 
$U_i(t)=  \frac{1}{v_i} \int_{C_i} U(x,t) dx$, the semi-discrete equations can be written:

$$
 \frac{\mathrm{d} U_i}{\mathrm{d} t} + \sum_{j \in N(i)} \frac{s_{ij}}{v_i}\left(\overrightarrow \Phi^{conv}_{ij} +  \overrightarrow{\Phi}^{diff}_{ij}\right)= S_i(U,x).
$$

with: 
- the numerical convection flux 

$$
\overrightarrow{\Phi}_{ij}^{conv}= \frac{1}{s_{ij}}\int_{\partial C_{ij}}\mathcal F^{conv}(U)\cdot\vec n_{ij}ds,
$$

- the numerical  diffusion flux 

$$
\overrightarrow{\Phi}_{ij}^{diff}= \frac{1}{s_{ij}}\int_{\partial C_{ij}}\mathcal {F}^{diff}(U)\cdot\vec n_{ij}ds.
$$

To approximate the convection numerical flux $\overrightarrow{\Phi}^{conv}_{ij}$ we solve an  approximate Riemann problem 
at the interface $\partial C_{ij}$. There are three possible formulations for the convection fluxes. 
- Using the \ref roe local linearisation of the fluxes, we obtain the following formula:

$$
\overrightarrow{\Phi}^{conv, Roe}_{ij}= \frac{\mathcal{F}^{conv}(U_i) + \mathcal{F}^{conv}(U_j)}{2} \vec{n}_{ij}- \mathcal{D}(U_i,U_j) \frac{U_j-U_i}{2}\\
=\mathcal{F}^{conv}(U_i) \vec{n}_{ij} + A^-(U_i,U_j) (U_j - U_i),
$$

- Using the \ref vfroe local linearisation of the fluxes, we obtain the following formula:

$$
\overrightarrow{\Phi}^{conv, VFRoe}_{ij}=\mathcal{F}^{conv}\left(\frac{U_i + U_j}{2} - \mathcal{D}(U_i,U_j) \frac{U_j-U_i}{2}\right)\vec{n}_{ij},
$$

- Using the \ref vffc local linearisation of the fluxes, we obtain the following formula:

$$
\overrightarrow{\Phi}^{conv, VFFC}_{ij}= \frac{\mathcal{F}^{conv}(U_i) + \mathcal{F}^{conv}(U_j)}{2} \vec{n}_{ij}- \mathcal{D}(U_i,U_j) \frac{\mathcal{F}^{conv}(U_j)-\mathcal{F}^{conv}(U_i)}{2} \vec{n}_{ij},
$$

where 
- $\mathcal{D}$ is an upwinding matrix,
- $A(U_i,U_j)$ the Roe matrix
- and $A^-= \frac{A - \mathcal{D}}{2}$.

 The choice $\mathcal{D}= 0$ gives the \ref centered upwinding, 
 $\mathcal{D}= |A|$ for the \ref roe formulation 
 and  $\mathcal{D}= sign(A)$ for the \ref vfroe and \ref vffc formulations give the full \ref upwind upwinding. The \ref lowMach, \ref pressureCorrection and \ref staggered upwindings allow more precision for  Mach number flows. 

The diffusion numerical flux $\overrightarrow\Phi_{ij}^{diff}$ is approximated on structured meshes using the formula:

$$
\overrightarrow \Phi_{ij}^{diff}= D (\frac{U_i+U_j}{2},\vec{n}_{ij})(U_j-U_i),
$$

where 
- the numerical diffusion tensor is 

  $$
  D(U,\vec{n}_{ij})=\nabla\mathcal{F}^{diff}(U) \cdot \vec{n}_{ij}.
  $$
  
  The expression of $\Phi_{ij}^{diff}$ above is not accurate for highly non structured or non conforming meshes. 
  However, since we are mainly interested in convection driven flows, we do not ask for a very precise scheme.

Finally, since $\sum_{j \in N(i)}\mathcal {F}^{conv}(U_i)\cdot \vec{n}_{ij}=0$, 
using *(\ref{eq:flux roe})* and *(\ref{eq:flux diff})* the equation *(\ref{eq:numer scheme})* of the semi-discrete scheme becomes:

$$
\frac{\mathrm{d} U_{i}}{\mathrm{d} t} + \sum_{j\in N(i)} {\frac{s_{ij}}{v_i}\{(A^-+ D)(U_i,U_j)\}(U_j-U_i)} = S_i(U,x), 
$$

The source term in *(\ref{eq:reduced scheme})* can be approximated using either a 

$$
 \textrm{ Centered source } S_i=S(U_i)\nonumber
$$

or an

$$
 \textrm{ Upwind source } S_i=\frac{1}{2}(Id-sign(A^{Roe}_{i,i+1}))\frac{S(U_i)+S(U_{i+1})}{2}
			      +\frac{1}{2}(Id+sign(A^{Roe}_{i-1,i}))\frac{S(U_{i-1})+S(U_i)}{2}.
$$



Explicit schemes
----------------

In explicit schemes, in order to compute the values $U_i^{n+1}$, 
the convection flux $\overrightarrow\Phi_{ij}^{conv}$, 
the diffusion flux $\overrightarrow\Phi_{ij}^{diff}$ 
and the source term $S(U,x)$ in *(\ref{eq:numer scheme})* 
are evaluated at time $n$ as :

$$
\frac{U_{i}^{n+1} - U_{i}^{n}}{\Delta t} + \sum_{j\in N(i)} \frac{s_{ij}}{v_i}\left(\frac{1}{2}(\mathcal{F}^{conv}(U_i^n) + \mathcal{F}^{conv}(U_j^n))\cdot \vec{n}_{ij}- \mathcal{D}(U_i^n,U_j^n,\vec{n}_{ij}) \frac{U_j^n-U_i^n}{2}\right)
+\frac{s_{ij}}{v_i} D (\frac{U_i+U_j}{2},\vec{n}_{ij})(U_j-U_i)= S(U^n,x_i),
$$

or equivalently using *(\ref{eq:flux roe})*$ 
and *(\ref{eq:flux diff})*

$$
\frac{U_{i}^{n+1} - U_{i}^{n}}{\Delta t} + \sum_{j\in N(i)} {\frac{s_{ij}}{v_i}\{(A^-+ D)(U_i^{n},U_j^{n},\vec{n}_{ij})\}(U^{n}_j-U^{n}_i)} =  S(U^n,x_i). 
$$

From the system *(\ref{explicitscheme})* 
we can obtain $U_i^{n+1}$ easily using matrix-vector products and vector sum. 
However the time step of explicit methods is constrained by the CFL condition for stability reasons.

Implicit schemes
----------------

In implicit schemes, in order to compute the values $U_i^{n+1}$, 
the fluxes $\Phi^{conv}_{ij}$, 
$\Phi^{diff}_{ij}$ 
and the source term $S(U,x)$ 
in *(\ref{eq:numer scheme})* 
are evaluated at time $n+1$ :

$$
\frac{U_{i}^{n+1} - U_{i}^{n}}{\Delta t} + \sum_{j\in N(i)} \frac{s_{ij}}{v_i}\left(\frac{1}{2}(\mathcal{F}^{conv}(U_i^{n+1}) + \mathcal{F}^{conv}(U_j^{n+1})). \vec{n}_{ij}- \mathcal{D}(U_i^{n+1},U_j^{n+1},\vec{n}_{ij}) \frac{U_j^{n+1}-U_i^{n+1}}{2}\right)
+\frac{s_{ij}}{v_i} D (\frac{U_i+U_j}{2},\vec{n}_{ij})(U_j-U_i) = S(U^{n+1},x_i),
$$

or equivalently using *(\ref{eq:flux roe})* and *(\ref{eq:flux diff})*

$$
\frac{U_{i}^{n+1} - U_{i}^{n}}{\Delta t} + \sum_{j\in N(i)} {\frac{s_{ij}}{v_i}\{(A^-+ D)(U_i^{n+1},U_j^{n+1},\vec{n}_{ij})\}(U^{n+1}_j-U^{n+1}_i)} =  S(U^{n+1},x_i). 
$$

The system *(\ref{implicitscheme})* is nonlinear. The computation of $U_i^n$ is more expensive but we can expect to use larger time steps than would be allowed with the explicit scheme.

We use the following Newton iterative method to obtain the required solutions:

$$
\frac{\delta U_i^{k+1}}{\Delta t} + \sum_{j \in N(i)} \frac{s_{ij}}{v_i} \left[( A^-+ D)(U_i^k,U_j^k) \right] \left(\delta U_j^{k+1}- \delta U_i^{k+1} \right)\\
 = - \frac{U^k_i-U^n_i}{\Delta t} - \sum_{j \in N(i)} \frac{s_{ij}}{v_i} \left[( A^-+ D)(U_i^k,U_j^k)\right] (U^k_j-U^k_i),
$$

where :
- $\delta U_i^{k+1} = U_i^{k+1} - U_i^{k}$ is the variation 
  of the $k$-th iterate that approximate the solution 
  at time $n+1$. 

Defining the unknown vector $\mathcal U = (U_1,\dots,U_N)^t$, 
each Newton iteration for the computation of $\mathcal U$ 
at time step $n+1$ requires the numerical solution of the following linear system:

$$
 \mathcal A(\mathcal U^k)\delta \mathcal U^{k+1} =  b(\mathcal U^n, \mathcal U^{k}).
$$



Numerical scheme for the Navier-Stokes equations
------------------------------------------------

For the Navier-Stokes equations $U=(\rho, \vec q, \rho E)^t$ and the fluxes 
in *(\ref{Navierstokes})* write

$$
\mathcal{F}^{conv}(U)= \left(\begin{array}{c}
\vec{q} \\
\vec{q} \otimes \frac{\vec{q}}{\rho} + p {I}_d \\
\left( \rho E + p \right) \frac{\vec{q}}{\rho} 
\end{array}
 \right) ,\quad
\mathcal{F}^{diff}(U)=\left(\begin{array}{c}
0\\
-\nu \vec \nabla (\frac{\vec{q}}{\rho}) \\
-\lambda \vec \nabla T 
\end{array}
 \right).
$$

For the Euler equations, we can build the \ref roe matrix $A(U_i,U_j)* explicitly using the Roe averaged state $U_{Roe}(U_i,U_j)=(\tilde\rho, \tilde \rho\tilde u, \tilde\rho \tilde E=\tilde \rho\tilde H-\tilde p)^t$ defined by

$$
\tilde{\rho}=\sqrt{\rho_{i}\rho_{j}}
$$

$$
\tilde{u}=\frac{\sqrt{\rho_{i}}u_{i}+\sqrt{\rho_{j}}u_{j}}{\sqrt{\rho_{i}}+\sqrt{\rho_{j}}}
$$

$$
\tilde{H}=\frac{\sqrt{\rho_{i}}H_{i}+\sqrt{\rho_{j}}H_{j}}{\sqrt{\rho_{i}}+\sqrt{\rho_{j}}}.
$$

The Roe matrix writes (see \ref leveque)

$$
A_{Roe}(U_i,U_j)=\nabla\mathcal{F}^{conv}(U_{Roe}(U_i,U_j))\vec{n}_{ij}=
\left(\begin{array}{ccc}
       0  & 1 & 0\\
       \tilde{\chi}+(\frac{1}{2}\tilde{\kappa} -1)\tilde{u}^2 & (2-\tilde{\kappa})\tilde{u} & \tilde{\kappa}\\
       (\tilde{\chi}+\frac{1}{2}\tilde{\kappa} \tilde{u}^2- \tilde{H})\tilde{u} & \tilde{H}-\tilde{\kappa} \tilde{u}^2 & (\tilde{\kappa}+1)\tilde{u}
      \end{array}\right)
$$

The diffusion numerical flux $\overrightarrow\Phi_{ij}^{diff}$ is approximated with the formula:

$$
\overrightarrow \Phi_{ij}^{diff}= D (\frac{U_i+U_j}{2})(U_j-U_i)
$$

with the matrix 

$$
D(U)= 
\left(
\begin{array}{ccc}
 0&\vec 0& 0\\
\frac{\nu \vec q}{\rho^2}& \frac{-\nu}{\rho} I_d&0\\
\frac{\lambda}{c_v}\left(\frac{c_vT}{\rho}-\frac{||\vec q||^2}{2\rho^3}\right)&\quad \frac{{\vec q}^{\;t} \lambda}{\rho^2 c_v}&\quad-\frac{\lambda}{c_v  \rho}
\end{array}
\right) , 
$$ 

where $c_v$ is the heat capacity at constant volume.

