## The Diffusion problem on the square

We consider the following diffusion problem with Dirichlet boundary conditions

$$
\left\{\begin{array}{c}
-\left(\partial_{xx} u + K\partial_{yy} u\right) = f \textrm{ on } \Omega\\
u=0 \textrm{ on } \partial\Omega
\end{array}\right.
$$

on the square domain $\Omega= [0,1]\times [0,1]$ with  

$$f=(1+K)\pi^2 sin(\pi x) sin(\pi y).$$  
The unique solution of the problem is  

$$
u=sin(\pi x) sin(\pi y).
$$

The Diffusion equation can be written in a matrix form
$$
-\nabla\cdot(D\vec\nabla u)=f
$$
and the associated diffusion matrix is
$$
D=\left(\begin{array}{cc}
    1 & 0\\
    0 & K
  \end{array}\right)
$$

We are interested in case where $K\gg 1$. In the following numerical results we take the value $K=10^4$.
