The two-phase flow models
=========================

We present the homogeneised two phase flow models implemented in CoreFlows. 

This models are obtained by averaging the balance equations for each separated phase or for the mixture, using space, time or ensemble averaged quantities (\ref ishii and \ref Drew ). 

\ref The drift model is used in the thermal hydraulics softwares \ref flica4 and \ref flocal, whilst the two-fluid models are used in \ref cathare , \ref neptuneCFD, \ref CobraTF , \ref relap5 .


The Drift model
---------------

The drift model is a system of four nonlinear equations taking the following conservative form
$$
\left\{\begin{array}{lll}
         \partial_t(\alpha_g\rho_g+\alpha_l\rho_l)&+\nabla\cdot(\alpha_g\rho_g{}^t\vec{u}_g+\alpha_l\rho_l{}^t\vec{u}_l)&=0\\
         \partial_t(\alpha_g\rho_g)&+\nabla\cdot(\alpha_g\rho_g{}^t\vec{u}_g)&=\Gamma_g(h_m,\Phi)\\
         \partial_t(\alpha_g\rho_g\vec{u}_g+\alpha_l\rho_l\vec{u}_l)&+\nabla\cdot(\alpha_g\rho_g\vec{u}_g\otimes\vec{u}_g+\alpha_l\rho_l\vec{u}_l\otimes\vec{u}_l+p {I}_d)&=\rho_m\vec{g}-K_g\alpha_g\rho_g||\vec{u}_g||\vec{u}_g-K_l\alpha_l\rho_l||\vec{u}_l||\vec{u}_l\\
         \partial_t(\alpha_g\rho_g E_g+\alpha_l\rho_l E_l)&+\nabla\cdot(\alpha_g\rho_g H_g{}^t\vec{u}_g+\alpha_l\rho_l H_l{}^t\vec{u}_l)&=\Phi+\rho\vec{g}\cdot\vec{u}-K_g\alpha_g\rho_g||\vec{u}_g||^3-K_l\alpha_l\rho_l||\vec{u}_l||^3
        \end{array}\right.,
$$
where the total energy and total enthalpy are defined by
$$
E_k=e_k+\frac{1}{2}|\vec{u}_k|^2,\quad H_k=h_k+\frac{1}{2}|\vec{u}_k|^2,\qquad k=v,l,
$$
where $e_k$ is the internal energy, and $h_k=e_k+\frac{p}{\rho_k}$ the enthalpy associated to phase $k$ and
$$
\begin{array}{lll}
\rho_m&=&\alpha_g\rho_g+\alpha_l\rho_l\\
\vec{u}_m&=&\frac{\alpha_g\rho_g\vec{u}_g+\alpha_l\rho_l\vec{u}_l}{\alpha_g\rho_g+\alpha_l\rho_l}\\
h_m&=&\frac{\alpha_g\rho_g h_g+\alpha_l\rho_l h_l}{\alpha_g\rho_g+\alpha_l\rho_l}.
\end{array}
$$
We need a drift correlation for the relative velocity:
$$
\vec{u}_r=\vec{u}_g-\vec{u}_l=\vec{f}_r(c_g,\vec{u}_m,\rho_m).
$$
The phase change is modeled using the formula
$$
 \Gamma_g=\left\{\begin{array}{cc}
         \frac{\Phi}{\mathcal{L}}&\textrm{ if } h_l^{sat}\leq h< h_g^{sat} \textrm{ and } 0<\alpha_g<1\\[1.5ex]
         0& \textrm{ otherwise }
        \end{array}\right..
$$
The parameters $\lambda_k, \nu_k,\vec g, K_k $ and $\Phi$ can be set by the user.

[More details about the drift model are available here](TwoPhase/DriftModelPage.ipynb)

	
The isothermal two-fluid model
-----------------------------------------------

The model consists in the phasic mass and momentum balance equations.

The main unknowns are $\alpha$, $P$, $\vec{u}_g$, $\vec{u}_l$. The model uses stiffened gas laws $p_g(\rho_g)$ and  $p_l(\rho_l)$ for a contant temperature $T_0$ provided by the user.

The subscript $k$ stands for $l$ for the liquid phase and $g$ for the gas phase. The common
averaged pressure of the two phases is denoted by $p$. 

In our model, pressure equilibrium between the two phases is postulated, and the resulting system to solve is:
$$
\left\{
\begin{array}{ccll}
 \frac{\partial m_g}{\partial t}& +& \nabla \cdot \vec{q}_g &= 0,\\[1.5ex]
\frac{\partial m_l}{\partial t} &+ &\nabla \cdot \vec{q}_l &= 0,\\[1.5ex]
\frac{\partial \vec{q}_g}{\partial t}& +& \nabla \cdot (\vec{q}_g\otimes\frac{\vec{q}_g}{m_g})+ \alpha_g \vec\nabla p&\\[1.5ex] 
 &+&\Delta p \nabla \alpha_g -\nu_g\Delta \vec{u}_g &= m_g\vec{g}-K_gm_g||\vec{u}_g||\vec{u}_g\\[1.5ex]
\frac{\partial \vec{q}_l}{\partial t}& +& \nabla \cdot (\vec{q}_l\otimes\frac{\vec{q}_l}{m_l})+ \alpha_l \vec\nabla p&\\[1.5ex]
&+&\Delta p \nabla \alpha_l -\nu_l\Delta \vec{u}_l &= m_l\vec{g}-K_lm_l||\vec{u}_l||\vec{u}_l,\\
\end{array}
\right.
$$

Here :
- $\nu_k$ is the viscosity of phase $k$,
- $\Delta p$ denotes the pressure default $p-p_k$ between the bulk average pressure and the interfacial average pressure.

where 
$$ 
\left\{\begin{array}{clc}
	\alpha_g +\alpha_l &=& 1 \\[1.5ex]
	m_k &=& \alpha_k \rho_k \\[1.5ex]
	\vec{q}_k &=& \alpha_k \rho_k \vec{u}_k \\[1.5ex]
        	\end{array}\right..
$$

The parameters $\lambda_k, \nu_k,\vec g, K_k $ and $\Phi$ can be set by the user.

[More details about the isothermal two-fluid model are available here](IsothermalPage.ipynb)


The five equation two-fluid model
-----------------------------------------------


The model consists in the phasic mass and momentum balance equations and one mixture total energy balance equation. 

The main unknowns are $\alpha$,$P$,$\vec{u}_g$,$\vec{u}_l$ and $T=T_g=T_l$. 

The model uses stiffened gas laws $p_g(\rho_g,T)$ and  $p_l(\rho_l,T)$.

$$
\left\{
\begin{array}{ccll}
 \frac{\partial m_g}{\partial t}& +& \nabla \cdot \vec{q}_g &= \Gamma_g(h_g,\Phi),\\[1.5ex]
\frac{\partial m_l}{\partial t} &+ &\nabla \cdot \vec{q}_l &= \Gamma_l(h_l,\Phi),\\[1.5ex]
\frac{\partial \vec{q}_g}{\partial t}& +& \nabla \cdot (\vec{q}_g\otimes\frac{\vec{q}_g}{m_g})+ \alpha_g \nabla p&\\[1.5ex] 
 &+&\Delta p \nabla \alpha_g -\nu_g(\Delta \frac{\vec{q}_g}{m_g}) &= m_g\vec{g}-K_gm_g||\vec{u}_g||\vec{u}_g\\[1.5ex]
\frac{\partial \vec{q}_l}{\partial t}& +& \nabla \cdot (\vec{q}_l\otimes\frac{\vec{q}_l}{m_l})+ \alpha_l \nabla p&\\[1.5ex]
&+&\Delta p \nabla \alpha_l -\nu_l(\Delta \frac{\vec{q}_l}{m_l}) &= m_l\vec{g}-K_lm_l||\vec{u}_l||\vec{u}_l,\\[1.5ex]
\partial_t\rho_mE_m&+&\nabla\cdot(\alpha_g\rho_g H_g{}^t\vec{u}_g+\alpha_l\rho_l H_l{}^t\vec{u}_l)&=\Phi+\rho\vec{g}\cdot\vec{u}-K_gm_g||\vec{u}_g||^3-K_lm_l||\vec{u}_l||^3
\end{array}
\right. \nonumber
$$
where
$$
\begin{array}{ccl}
\rho_m&=&\alpha_g\rho_g+\alpha_l\rho_l\\
 E_m&=&\frac{\alpha_g\rho_g E_g+\alpha_l\rho_l E_l}{\alpha_g\rho_g+\alpha_l\rho_l}.
\end{array}
$$

The phase change is modeled using the formula

$$
\Gamma_g=\left\{\begin{array}{cc}
         \frac{\Phi}{\mathcal{L}}&\textrm{ if } h_l^{sat}\leq h< h_g^{sat} \textrm{ and } 0<\alpha_g<1\\[1.5ex]
         0& \textrm{ otherwise }
        \end{array}\right..
$$

The parameters $\lambda_k, \nu_k,\vec g, K_k $ and $\Phi$ can be set by the user.

[More details about the five equation two-fluid model are available here](TwoPhase/FiveEqPage.ipynb)

