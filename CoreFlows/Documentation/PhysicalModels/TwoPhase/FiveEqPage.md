The five equation two-fluid model
=================================

The model consists in five balance laws :
- the gas mass balance
- the liquid mass balance
- the gas momentum balance
- the liquid momentum balance
- the mixture energy balance. 

The main unknowns are $\alpha, P, \vec{u}_g, \vec{u}_l$ and $T=T_g=T_l$. 

The model uses stiffened gas laws $p_g(\rho_g,T)$ and  $p_l(\rho_l,T)$ and linearised internal energy laws $e_g(T)$ and  $e_l(T)$ fitted by either around 1 bar and 373K or around 155 bars and 618K (see pressureEstimate).

$$
 \frac{\partial m_g}{\partial t} + \nabla \cdot \vec{q}_g = \Gamma_g(h_g,\Phi),
$$
 
$$
\frac{\partial m_l}{\partial t} + \nabla \cdot \vec{q}_l = \Gamma_l(h_l,\Phi),
$$
 
$$
\frac{\partial \vec{q}_g}{\partial t} + \nabla \cdot (\vec{q}_g\otimes\frac{\vec{q}_g}{m_g})+ \alpha_g \nabla p 
 +\Delta p \nabla \alpha_g -\nu_g(\Delta \frac{\vec{q}_g}{m_g}) = m_g\vec{g}-K\rho_m||\vec{u}_g-\vec{u}_l||(\vec{u}_g-\vec{u}_l)-K_s\delta(x)m_g||\vec{u}_g||\vec{u}_g
$$
 
$$
\frac{\partial \vec{q}_l}{\partial t} + \nabla \cdot (\vec{q}_l\otimes\frac{\vec{q}_l}{m_l})+ \alpha_l \nabla p
+\Delta p \nabla \alpha_l -\nu_l(\Delta \frac{\vec{q}_l}{m_l}) = m_l\vec{g}-K\rho_m||\vec{u}_l-\vec{u}_g||(\vec{u}_l-\vec{u}_g)-K_s\delta(x)m_l||\vec{u}_l||\vec{u}_l
$$
 
$$
\partial_t\rho_m E_m+\nabla\cdot(\alpha_g\rho_g H_g{}^t\vec{u}_g+\alpha_l\rho_l H_l{}^t\vec{u}_l)=\Phi+\rho\vec{g}\cdot\vec{u}_m-K\rho_m||\vec{u}_g-\vec{u}_l||^3-K_s\delta(x)(m_g||\vec{u}_g||^3+m_l||\vec{u}||^3)
$$

where the mixture quantities are defined by

$$
\rho_m=\alpha_g\rho_g+\alpha_l\rho_l
$$
 
$$
\vec{u}_m=\frac{\alpha_g\rho_g\vec{u}_g+\alpha_l\rho_l\vec{u}_l}{\alpha_g\rho_g+\alpha_l\rho_l}
$$
 
$$
E_m=\alpha_g\rho_g E_g+\alpha_l\rho_l E_l,
$$

whereas the quantities associated to each to phase $k=g,l$ are defined as
- $\alpha_k$ is the phasic volumic presence rate,
- $\rho_k$ is the phasic density,
- $m_k=\alpha_k\rho_k$ is the phasic partial density,
- $\vec u_k$ the phasic velocity,
- $\vec q_k = \rho \vec u$ the phasic momentum,
- $p$ the common phasic pressure,
- $e_k$ the phasic internal energy,
- $E_k = e_k + \frac{||\vec u||^2}{2}$ the phasic total energy,
- $h_k=e_k+\frac{p}{\rho_k}$ the phasic enthalpy
- $H_k=h_k+\frac{1}{2}|\vec{u}_k|^2$ the phasic total enthalpy
- $T$ the common absolute temperature,
- $\nu_k$ the viscosity ([FiveEqsTwoFluid](../../../Models/inc/FiveEqsTwoFluid.hxx)::setViscosity),
- $\lambda_k$ the thermal conductivity ([FiveEqsTwoFluid](../../../Models/inc/FiveEqsTwoFluid.hxx)::setConductivity),
- $K$ the interphase friction coefficient ([FiveEqsTwoFluid](../../../Models/inc/FiveEqsTwoFluid.hxx)::setDragCoeffs),

Geometric and physical source terms are
- $\vec g$ the gravity vector ([FiveEqsTwoFluid](../../../Models/inc/FiveEqsTwoFluid.hxx)::setGravity)
- $\Phi(\vec x)$ the heat power received by the fluid ([FiveEqsTwoFluid](../../../Models/inc/FiveEqsTwoFluid.hxx)::setHeatPowerField),
- $K_s(\vec x)$ the singular friction function, $\delta_s(\vec x)$ the Dirac delta function with support on the set $s$ ([FiveEqsTwoFluid](../../../Models/inc/FiveEqsTwoFluid.hxx)::setPressureLossField).

We close the this system with two stiffened gas lasw $p = (\gamma_k -1) \rho_k e_k -\gamma_k p_{0k}$ for each phase 
and a linearised internal energy law $h_k(T)$ 
valid around the points $(P=1 bar, T=300K)$ 
or $(P=155 bars, T=618K)$ depending on the value of the enum pressureEstimate.

For the sake of simplicity, for the moment we consider constant viscosity and conductivity, and neglect the contribution of viscous forces in the energy equation. The constant parameters $\lambda_k, \nu_k,\vec g, K_k$ 
and the fields $\phi(\vec x), \Phi(\vec x), K_s(\vec x)$ can be set by the user. 
The default value for $\phi$ is $\phi=1$.

The phase change is modeled using the formula

$$
 \Gamma_g=\frac{\Phi}{\mathcal{L}}\textrm{ if } T^{sat}\leq T \textrm{ and } 0<\alpha_g<1
$$

$$
  \Gamma_g= 0 \textrm{ otherwise }
$$

The parameters $\lambda_k, \nu_k,\vec g, K$ 
and $\Phi$ can be set by the user.

The class : [FiveEqsTwoFluid](../../../Models/inc/FiveEqsTwoFluid.hxx) implements the equal temperature two fluid model  

\subpage Example5EqPage "Here are C and Python example scripts using the five equation two-fluid model "

