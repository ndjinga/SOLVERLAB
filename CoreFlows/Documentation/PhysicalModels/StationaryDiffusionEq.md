The stationary diffusion equation
=================================

$$
 -\lambda\triangle T = \Phi+\lambda_{sf}(T_f-T)
$$
where
- $T$ the main unknown is the solid temperature field
- $\lambda$ is the solid thermal conductivity possibly set by the user (default value is 1)
- $\Phi$ is the heat source term possibly set by the user (default value is 0)
- $\lambda_{sf}$ is the fluid-solid heat transfer coefficient set by the user (default value is 0)
- $T_f$ is the fluid temperature field provided by the user

The class [StationaryDiffusionEquation](../../Models/inc/StationaryDiffusionEquation.hxx) implements a scalar stationary diffusion equation for the temperature in a solid.  


\subpage ExampleStationaryDiffusionEqPage "Here are C and Python example scripts"


