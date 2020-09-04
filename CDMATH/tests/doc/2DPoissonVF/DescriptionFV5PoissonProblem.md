## The FV5 scheme for the Laplace equation

The domain $\Omega$ is decomposed into cells $C_i$.

$|C_i|$ is the measure of the cell $C_i$.

$f_{ij}$ is the interface between two cells $C_i$ and $C_j$. 

$s_{ij}$ is the measure of the interface $f_{ij}$.

$d_{ij}$ is the distance between the centers of mass of the two cells $C_i$ and $C_j$.

The discrete Poisson problem is
$$
-\frac{1}{|C_i|} \sum s_{ij}F_{ij}=f_i,
$$
where
$u_i$ is the approximation of $u$ in the cell $C_i$,

$f_i$ is the approximation of $f$ in the cell $C_i$,

$F_{ij}$ is a numerical approximation of the outward normal diffusion flux from cell $i$ to cell $j$.

In the case of the scheme FV5, the flux formula are
$$
F_{ij}=\frac{u_j-u_i}{d_{ij}},
$$
for two cells $i$ and $j$ inside the domain,

and
$$
F_{boundary}=\frac{u(x_f)-u_i}{d_{if}},
$$
for a boundary face with center $x_f$, inner cell $i$ and distance between face and cell centers $d_{if}$
