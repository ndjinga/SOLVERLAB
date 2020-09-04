The Drift model
===============	

The model consists in the steam mass balance equation together with the mixture mass conservation, the mixture momentum balance and mixture energy balance equations. The main unknowns are the steam mass concentration $c_v$, the pressure $P$, the mixture velocity $\vec{u}_m$, and the common temperature $T$. The model uses stiffened gas laws $p_g(\rho_g,T)$ and  $p_l(\rho_l,T)$ as well as  linearised internal energy law $e_k(T)$ valid around the saturation points $(P=1 bar, T=373K)$ or $(P=155 bars, T=618K)$ depending on the value of the enum \ref pressureEstimate.

The drift model is a system of four nonlinear equations taking the following conservative form
$$
\left\{\begin{array}{lll}
         \partial_t(\phi \rho_m) &+\nabla\cdot(\phi\rho_m\vec{u}_m)&=0\\
         \partial_t(\phi m_g)&+\nabla\cdot\phi(m_g\vec{u}_g)&=\phi\Gamma_g(h_m,\Phi)\\
         \partial_t(\phi\rho_m\vec{u}_m)&+\nabla\cdot\phi(m_g\vec{u}_g\otimes\vec{u}_g+ m_l\vec{u}_l\otimes\vec{u}_l)+\vec{\nabla}(p\phi)&=p\vec{\nabla}\phi+\phi\rho_m\vec{g}- K_g\phi m_g||\vec{u}_g||\vec{u}_g- K_l\phi m_l||\vec{u}_l||\vec{u}_l- K_s\delta_s(x)\phi\rho_m||\vec{u}_m||\vec{u}_m\\
         \partial_t\phi (\rho_m E_m)&+\nabla\cdot\phi(m_g H_g{}^t\vec{u}_g+m_l H_l{}^t\vec{u}_l)&=\Phi+\phi\rho_m\vec{g}\cdot\vec{u}_m- K_g\phi m_g||\vec{u}_g||^3- K_l\phi m_l||\vec{u}_l||^3- K_s\delta_s(x)\phi\rho_m||\vec{u}_m||^3
        \end{array}\right.,
$$
where the mixture quantities are defined by
$$
\begin{array}{lll}
\rho_m&=&\alpha_g\rho_g+\alpha_l\rho_l\\
\vec{u}_m&=&\frac{\alpha_g\rho_g\vec{u}_g+\alpha_l\rho_l\vec{u}_l}{\alpha_g\rho_g+\alpha_l\rho_l}\\
E_m&=&\alpha_g\rho_g E_g+\alpha_l\rho_l E_l,
\end{array}
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
- $\nu_k$ the viscosity ([DriftModel](../../../Models/inc/DriftModel.hxx)::setViscosity),
- $\lambda_k$ the thermal conductivity ([DriftModel](../../../Models/inc/DriftModel.hxx)::setConductivity),
- $K_k$ the phasic regular friction coefficient ([DriftModel](../../../Models/inc/DriftModel.hxx)::setDragCoeffs),

Geometric and physical source terms are
- $\vec g$ the gravity vector ([DriftModel](../../../Models/inc/DriftModel.hxx)::setGravity)
- $\Phi(\vec x)$ the heat power received by the fluid ([DriftModel](../../../Models/inc/DriftModel.hxx)::setHeatPowerField),
- $\phi(\vec x)$ the volumic porosity field ([DriftModel](../../../Models/inc/DriftModel.hxx)::setPorosityField),
- $K_s(\vec x)$ the singular friction function, $\delta_s(\vec x)$ the Dirac delta function with support on the set $s$ ([DriftModel](../../../Models/inc/DriftModel.hxx)::setPressureLossField).

We close the Drift-Model system with a stiffened gas law $p = (\gamma_k -1) \rho_k e_k -\gamma_k p_{0k}$ for each phase and a linearised enthalpy law $h_k(T)$ valid around the points $(P=1 bar, T=300K)$ or $(P=155 bars, T=618K)$ depending on the value of the enum \ref pressureEstimate.

For the sake of simplicity, for the moment we consider constant viscosity and conductivity, and neglect the contribution of viscous forces in the energy equation.

The constant parameters $\lambda_k, \nu_k,\vec g, K_k$ and the fields $\phi(\vec x),\: \Phi(\vec x),\: K_s(\vec x)$ can be set by the user. The default value for $\phi$ is $\phi=1$.


To close the system we need a drift correlation for the relative velocity:
$$
\vec{u}_r=\vec{u}_g-\vec{u}_l=\vec{f}_r(c_g,\vec{u}_m,\rho_m).
$$
For the moment the only drift correlation available is $\vec{u}_g=\vec{u}_l$.

The phase change is modeled using the formula
$$
 \Gamma_g=\left\{\begin{array}{cc}
         \frac{\Phi}{\mathcal{L}}&\textrm{ if } T^{sat}\leq T \textrm{ and } 0<\alpha_g<1\\[1.5ex]
         0& \textrm{ otherwise }
        \end{array}\right..
$$

For the moment the boiling temperature $T^{sat}$ is constant and can be changed using the function DriftModell::setSatTemp.

The class : [DriftModel](../../../Models/inc/DriftModel.hxx) implements the 4 equation drift model  

\subpage ExampleDriftModelPage "Here are C and Python example scripts using the Drift Model  "


