The Navier-Stokes equations
===========================
The model consists of the following three balance laws for the mass, the momentum and the energy:
 
$$
\left\{\begin{array}{cclclcc}
\frac{\partial \phi\rho}{\partial t}&+&\nabla\cdot\vec{\phi q} & & &= &0\\[1.5ex]
\frac{\partial \phi\vec{q}}{\partial t}&+&\nabla\cdot\left(\phi\vec{q} \otimes \frac{\vec{q}}{\rho})+\phi\vec{\nabla} p \right) &-& \nu \nabla\cdot(\phi\vec{\nabla}\vec {u})&=&p\vec{\nabla\phi}+\phi\rho\vec{g}- (K_r+K_s\delta_s(x))\phi \rho||\vec{u}||\vec{u}\\[1.5ex]
\frac{\partial(\phi\rho E)}{\partial t} &+&\nabla\cdot\left[\phi(\rho E + p) \frac{\vec{q}}{\rho}\right]&-&\lambda \nabla\cdot(\phi\vec{\nabla} T)&=&\Phi+\phi\rho\vec{g}\cdot\vec{u}-(K_r+K_s\delta_s(x))\phi \rho||\vec{u}||^3
\end{array}\right.,
$$

where 
- $\rho$ is the density,
- $\vec u$ the velocity,
- $\vec q = \rho \vec u$ the momentum,
- $p$ the pressure,
- $\rho e$ the volumic internal energy,
- $\rho E = \rho e + \frac{||\vec q||^2}{2 \rho}$ the volumic total energy,
- $T$ the absolute temperature,
- $\Phi(\vec x)$ the heat power received by the fluid  ([SinglePhase](../../Models/inc/SinglePhase.hxx)::setHeatPowerField),
- $\phi(\vec x)$ the volumic porosity field ([SinglePhase](../../Models/inc/SinglePhase.hxx)::setPorosityField),
- $\vec g$ the gravity vector ([SinglePhase](../../Models/inc/SinglePhase.hxx)::setGravity)
- $\nu$ the viscosity ([SinglePhase](../../Models/inc/SinglePhase.hxx)::setViscosity),
- $\lambda$ the thermal conductivity ([SinglePhase](../../Models/inc/SinglePhase.hxx)::setConductivity),
- $K_r$ the regular friction coefficient ([SinglePhase](../../Models/inc/SinglePhase.hxx)::setDragCoeffs),
- $K_s(\vec x)$ the singular friction function, $\delta_s(\vec x)$ the Dirac delta function with support on the set $s$ ([SinglePhase](../../Models/inc/SinglePhase.hxx)::setPressureLossField).

We close the Navier-Stokes system by the ideal gas law $p = (\gamma -1) \rho e$ for steam water and a stiffened gas law $p = (\gamma -1) \rho e -\gamma p_0$ for liquid water and a linearised internal energy law valid around the points $(P=1 bar, T=300K)$ or $(P=155 bars, T=618K)$ depending on the value of the enum \ref pressureEstimate.

For the sake of simplicity, for the moment we consider constant viscosity and conductivity, and neglect the contribution of viscous forces in the energy equation.

The constant parameters $\lambda, \nu,\vec g, K_r$ and the fields $\phi(\vec x),\: \Phi(\vec x),\: K_s(\vec x)$ can be set by the user.


The class : [SinglePhase](../../Models/inc/SinglePhase.hxx) implements the single phase model  

\subpage ExampleSinglePhase "Here are C and Python example scripts using the single phase model "

