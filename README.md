# fokker-plank-equation
This module contains a 1st order numerical discretisation of the Fokker-Plank equation with constant drift and diffusion. The current implementation uses a central difference diffusion term and a 1st order upwind scheme for the advective term of the PDE.

Temporal and spatial discretisation of the domain has been done using the principles of Finite Volume. In order to address the boundary node, a ghost-node approach has been carried out.

Grid independence has been shown and the results look similar to how they should analytically.

Future work:

- write up simple analytical solution and compare the L2-norm error of the two methods.
- Implement a simple TVD scheme to prevent overshoot, reduce numerical dispersion and increase order of accuracy. This could be using either Albada or UMIST.
- Extend for non-constant drift and diffusion terms.
