The diffusion equation
======================

$$
 \partial_t T =d\triangle T +\frac{ \Phi+\lambda_{sf}(T_f-T)}{\rho c_p}
$$
where
- $T$ the main unknown is the solid temperature field
- $\rho$ is the solid density assumed constant and possibly set by the user
- $c_p$ is the solid specific heat, possibly set by the user and assumed constant
- $\lambda$ is the solid thermal conductivity possibly set by the user
- $d=\frac{\lambda}{\rho c_p}$ is the solid diffusivity
- $\lambda_{sf}$ is the fluid-solid heat transfer coefficient provided by the user
- $\Phi$ is the heat source term set by the user
- $T_f$ is the fluid temperature field provided by the user

The class [DiffusionEquation](../../Models/inc/DiffusionEquation.hxx) implementing a scalar diffusion equation for the temperature in a solid. The default values for $\rho, c_p, \lambda$ are those of Uranium oxyde at $900 K$.  


\subpage ExampleDiffusionEqPage "Here are C and Python example scripts"


