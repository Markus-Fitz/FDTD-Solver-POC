# FDTD-Solver-POC
Proof of concept of FDTD algorithm

Simple Jupyter-Notebook containing the FDTD-algorithm and surrounding structures.
This is meant as a verification of the algorithm and not optimized for usability or performance.
The algorithm only considers vacuum as a medium.

## Derivation of the discretized update equations

## Maxwell's equations

$$
\nabla \cdot E = \frac{\rho}{\epsilon_0}
$$

$$
\nabla \cdot B = 0
$$

$$
\nabla \times E = - \frac{\partial B}{\partial t}
$$

$$
\nabla \times B = \mu_0 (J + \epsilon_0 \frac{\partial E}{\partial t})
$$

## Yee grid

The Yee grid FDTD method is based on a staggered grid with the electric and magnetic fields stored at different places within each cell.

The electric fields are separated into the x-, y- and z-components and stored on the edges of each cubic cell. The directions of the respective field is aligned with the cell edge on which it is placed.

The magnetic fields are stored in the middle of each cell face. The magnetic fields are placed on the cell faces for which the field direction and the normal vector of the cell face align. The By field components are stored in the middle of the cell-faces in y+ and y- directions of the cubic cell.

## Discretized operators

The functions specifying the time-evolution for the electric and magnetic field components require the evaluation of the cross-product of both the electric and the magnetic fields.

To come up with the discretized partial differential equations, the curl of a vector field

$$
\nabla \times A
$$

can first be expressed in three separate equations

$$
(\nabla \times A)_x = \frac{\partial A_y}{\partial z} - \frac{\partial A_z}{\partial y}
$$

$$
(\nabla \times A)_y = \frac{\partial A_z}{\partial x} - \frac{\partial A_x}{\partial z}
$$

$$
(\nabla \times A)_z = \frac{\partial A_x}{\partial y} - \frac{\partial A_y}{\partial x}
$$

. Each partial differential acting on a field $F$ at a point specified by the index i, j, and k for the three spacial directions (from now on noted as a superscript on the corresponding field) can be expressed by the approximations

$$
\frac{\partial F}{\partial x} |_{i, j, k} \approx \frac{F^{i+\frac{1}{2}, j, k} - F^{i-\frac{1}{2}, j, k}}{\Delta x}
$$

$$
\frac{\partial F}{\partial y} |_{i, j, k} \approx \frac{F^{i, j+\frac{1}{2}, k} - F^{i, j-\frac{1}{2}, k}}{\Delta y}
$$

$$
\frac{\partial F}{\partial x} |_{i, j, k} \approx \frac{F^{i, j, k+\frac{1}{2}} - F^{i, j, k-\frac{1}{2}}}{\Delta z}
$$

, where the $\Delta x$, $\Delta y$ and $\Delta z$ describe the length of the cells in each direction (which is equal to the distance in between each field quantity).

With this approximation, the cross product of the curl operator acting on a vector field can be expressed as

$$
(\nabla \times A)_x^{i, j, k} \approx \frac{A_y^{i, j, k+\frac{1}{2}} - A_y^{i, j, k-\frac{1}{2}}}{\Delta z} - \frac{A_z^{i, j+\frac{1}{2}, k} - A_z^{i, j-\frac{1}{2}, k}}{\Delta y}
$$

$$
(\nabla \times A)_y^{i, j, k} \approx \frac{A_z^{i+\frac{1}{2}, j, k} - A_z^{i-\frac{1}{2}, j, k}}{\Delta x} - \frac{A_x^{i, j, k+\frac{1}{2}} - A_x^{i, j, k-\frac{1}{2}}}{\Delta z}
$$

$$
(\nabla \times A)_z^{i, j, k} \approx \frac{A_x^{i, j+\frac{1}{2}, k} - A_x^{i, j-\frac{1}{2}, k}}{\Delta y} - \frac{A_y^{i+\frac{1}{2}, j, k} - A_y^{i-\frac{1}{2}, j, k}}{\Delta x}
$$

. The nature of the staggered Yee grid is very convenient for evaluating these approximations to the partial differential equations. As is seen in the approximations above, the field values shifted by half a cell distance is needed to evaluate the central differences. This fits perfectly with the staggered grid proposed by Yee in his paper "Numerical Solution of Initial Boundary Value Problems Involving Maxwell's Equations in Isotropic Media". The values of $E_x$, $E_y$, $E_z$ as well as $B_x$, $B_y$ and $B_z$ are placed in such a way that the central difference scheme described above comes naturally.

## Time evolution of fields

Assuming that the initial fields adhere to all four equations, the induction law of Faraday (equation nr. 3) and Ampère's law (equation nr. 4) describe the time-evolution of both the electric and magnetic fields. The time-step becomes relevant now, and is noted by the time-index $n$, used to describe the $n$-th timestep $t_n$.

The time-derivative $\frac{1}{\partial t}$ can be approximated in a similar way to the spacial derivatives. For a field $F(t)$, the approximation

$$
\frac{\partial F}{\partial t} \approx \frac{F(t_{n+1}) - F(t_n)}{\Delta t}
$$

is used, with $\Delta t$ being one time-step in the discretized time domain.

With this approximation, the Induction law of Faraday,

$$
\nabla \times E = - \frac{\partial B}{\partial t}
$$

, can be rewritten as the three equations

$$
\frac{E_y^{i, j, k+\frac{1}{2}} - E_y^{i, j, k-\frac{1}{2}}}{\Delta z} - \frac{E_z^{i, j+\frac{1}{2}, k} - E_z^{i, j-\frac{1}{2}, k}}{\Delta y} = \frac{B_x^{i, j, k}(t_{n+1}) - B_x^{i, j, k}(t_n)}{\Delta t}
$$

$$
\frac{E_z^{i+\frac{1}{2}, j, k} - E_z^{i-\frac{1}{2}, j, k}}{\Delta x} - \frac{E_x^{i, j, k+\frac{1}{2}} - E_x^{i, j, k-\frac{1}{2}}}{\Delta z} = \frac{B_y^{i, j, k}(t_{n+1}) - B_y^{i, j, k}(t_n)}{\Delta t}
$$

$$
\frac{E_x^{i, j+\frac{1}{2}, k} - E_x^{i, j-\frac{1}{2}, k}}{\Delta y} - \frac{E_y^{i+\frac{1}{2}, j, k} - E_y^{i-\frac{1}{2}, j, k}}{\Delta x} = \frac{B_z^{i, j, k}(t_{n+1}) - B_z^{i, j, k}(t_n)}{\Delta t}
$$

. This directly leads to an expression for the three vector fields $B_x$, $B_y$ and $B_z$ at the time-step $t_{n+1}$, based only on information we have at time-step $t_n$. The following three equations therefore are our update-equations for the B-field at the next timestep

$$
B_x^{i, j, k}(t_{n+1}) = B_z^{i, j, k}(t_n) + \Delta t \left( \frac{E_y^{i, j, k+\frac{1}{2}} - E_y^{i, j, k-\frac{1}{2}}}{\Delta z} - \frac{E_z^{i, j+\frac{1}{2}, k} - E_z^{i, j-\frac{1}{2}, k}}{\Delta y} \right)
$$

$$
B_y^{i, j, k}(t_{n+1}) = B_y^{i, j, k}(t_n) + \Delta t \left( \frac{E_z^{i+\frac{1}{2}, j, k} - E_z^{i-\frac{1}{2}, j, k}}{\Delta x} - \frac{E_x^{i, j, k+\frac{1}{2}} - E_x^{i, j, k-\frac{1}{2}}}{\Delta z} \right)
$$

$$
B_z^{i, j, k}(t_{n+1}) = B_z^{i, j, k}(t_n) + \Delta t \left( \frac{E_x^{i, j+\frac{1}{2}, k} - E_x^{i, j-\frac{1}{2}, k}}{\Delta y} - \frac{E_y^{i+\frac{1}{2}, j, k} - E_y^{i-\frac{1}{2}, j, k}}{\Delta x} \right)
$$

. Similarly, Ampère's law can be approximated in the following three equations

$$
\frac{B_y^{i, j, k+\frac{1}{2}} - B_y^{i, j, k-\frac{1}{2}}}{\Delta z} - \frac{B_z^{i, j+\frac{1}{2}, k} - B_z^{i, j-\frac{1}{2}, k}}{\Delta y} = \mu_0 \left( J_x^{i, j, k} + \epsilon_0 \frac{E_x^{i, j, k}(t_{n+1}) - E_x^{i, j, k}(t_n)}{\Delta t} \right)
$$

$$
\frac{B_z^{i+\frac{1}{2}, j, k} - B_z^{i-\frac{1}{2}, j, k}}{\Delta x} - \frac{B_x^{i, j, k+\frac{1}{2}} - B_x^{i, j, k-\frac{1}{2}}}{\Delta z} = \mu_0 \left( J_y^{i, j, k} + \epsilon_0 \frac{E_y^{i, j, k}(t_{n+1}) - E_y^{i, j, k}(t_n)}{\Delta t} \right)
$$

$$
\frac{B_x^{i, j+\frac{1}{2}, k} - B_x^{i, j-\frac{1}{2}, k}}{\Delta y} - \frac{B_y^{i+\frac{1}{2}, j, k} - B_y^{i-\frac{1}{2}, j, k}}{\Delta x} = \mu_0 \left( J_z^{i, j, k} + \epsilon_0 \frac{E_z^{i, j, k}(t_{n+1}) - E_z^{i, j, k}(t_n)}{\Delta t} \right)
$$

, leading to the update equations for the E-fields

$$
E_x^{i, j, k}(t_{n+1}) = E_x^{i, j, k}(t_n) + \frac{\Delta t}{\epsilon_0} \left( \frac{1}{\mu_0} \left( \frac{B_y^{i, j, k+\frac{1}{2}} - B_y^{i, j, k-\frac{1}{2}}}{\Delta z} - \frac{B_z^{i, j+\frac{1}{2}, k} - B_z^{i, j-\frac{1}{2}, k}}{\Delta y} \right) - J_x^{i, j, k} \right)
$$

$$
E_y^{i, j, k}(t_{n+1}) = E_y^{i, j, k}(t_n) + \frac{\Delta t}{\epsilon_0} \left( \frac{1}{\mu_0} \left( \frac{B_z^{i+\frac{1}{2}, j, k} - B_z^{i-\frac{1}{2}, j, k}}{\Delta x} - \frac{B_x^{i, j, k+\frac{1}{2}} - B_x^{i, j, k-\frac{1}{2}}}{\Delta z} \right) - J_y^{i, j, k} \right)
$$

$$
E_z^{i, j, k}(t_{n+1}) = E_z^{i, j, k}(t_n) +  \frac{\Delta t}{\epsilon_0} \left( \frac{1}{\mu_0} \left( \frac{B_x^{i, j+\frac{1}{2}, k} - B_x^{i, j-\frac{1}{2}, k}}{\Delta y} - \frac{B_y^{i+\frac{1}{2}, j, k} - B_y^{i-\frac{1}{2}, j, k}}{\Delta x} \right) - J_z^{i, j, k} \right)
$$
. 
