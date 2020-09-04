The scalar models
=================

The stationary diffusion equation
---------------------------------

$$
 -\lambda\triangle T = \Phi+\lambda_{sf}(T_f-T)
$$
where
- $T$ the main unknown is the solid temperature field
- $\lambda$ is the solid thermal conductivity possibly set by the user (default value is 1)
- $\Phi$ is the heat source term possibly set by the user (default value is 0)
- $\lambda_{sf}$ is the fluid-solid heat transfer coefficient set by the user (default value is 0)
- $T_f$ is the fluid temperature field provided by the user

See the [Stationary diffusion equation page](StationaryDiffusionEq.ipynb)

The diffusion equation	
----------------------
$$
 \partial_t T =d\triangle T +\frac{ \Phi+\lambda_{sf}(T_f-T)}{\rho c_p}
$$
where
- $T$ the main unknown is the rod temperature field
- $\rho$ is the rod density assumed constant (default value 10000)
- $c_p$ is the rod specific heat, provided by the user and assumed constant (default value 300)
- $d=\frac{\lambda}{\rho c_p}$ is the rod diffusivity  (default value 5/(10000*300))
- $\lambda_{sf}$ is the fluid-rod heat transfer coefficient provided by the user (default value 0)
- $\Phi$ is the heat source term if explicitely known (default value 0)
- $T_f$ is the fluid temperature field provided by the user

See the [Diffusion equation page](DiffusionEq.ipynb)

The transport equation	
----------------------
 
$$
 \partial_t H + \vec{u}\cdot\vec{\nabla} H = \Phi+\lambda_{sf}(T_s-T)
$$

where

- $ H $ the main unknown is the fluid enthalpy field
- $ \vec{u} $ is the constant transport velocity
- $ \Phi $ is the heat source term if explicitely known
- $ T_s $ is the rod temperature field provided by the user
- $ T=T_0+\frac{H-H_0}{c_p}$ is the fluid temperature field
- $ \lambda_{sf}$ is the fluid-rod heat transfer coefficient provided by the user
- $ c_p $ is the fluid specific heat, provided by the user and assumed constant

See the [Transport equation page](TransportEq.ipynb)


