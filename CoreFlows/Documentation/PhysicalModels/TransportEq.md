The transport equation
======================
 
$$
 \partial_t h + \vec{u}\cdot\vec{\nabla} h = \Phi+\lambda_{sf}(T_s-T)
$$

where

- $ h $ the main unknown is the fluid enthalpy field
- $ \vec{u} $ is the constant transport velocity set by the user
- $ \Phi $ is the heat source term if explicitely known
- $ T_s $ is the rod temperature field provided by the user
- $ T=T_0+\frac{h-h_0}{c_p}$ is the fluid temperature field
- $ \lambda_{sf}$ is the fluid-rod heat transfer coefficient provided by the user
- $ c_p $ is the fluid specific heat, provided by the user and assumed constant



The class [TransportEquation](../../Models/inc/TransportEquation.hxx) implements a scalar advection equation for the enthalpy of a fluid. The fluid can be either steam or liquid water around 1 bar or 155 bars.  

\subpage ExampleTransportEqPage "Here are C and Python example scripts"	


