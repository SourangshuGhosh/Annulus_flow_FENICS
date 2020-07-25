# Annulus_flow_FENICS
 
 Annulus_flow, a FENICS code which models the flow of a fluid, governed by the time dependent Navier-Stokes equations, in a 2D eccentric annulus.
 
To run this code under docker(), I issue the commands:

 ```html
- cd $HOME/fenicsproject/annulus_flow
- fenicsproject run
- python3 annulus_flow.py
- exit
```
## Reference:
Hans Petter Langtangen, Anders Logg,
Solving PDEs in Python - The FEniCS Tutorial Volume 1.

## Source Code:

The code normally takes 1000 time steps. For documentation, we only take 10 steps.

- annulus_flow.py a code which sets up and solves the problem.
- annulus_flow.sh runs all the tests.
- annulus_flow.txt the output file.
- annulus_flow_mesh.png a plot of the mesh.

## The pressure scalar field is saved on each step.

The following files can be found inside the pressure solutions directory
- pressure_solution000000.vtu
- pressure_solution000001.vtu
- pressure_solution000002.vtu
- pressure_solution000003.vtu
- pressure_solution000004.vtu
- pressure_solution000005.vtu
- pressure_solution000006.vtu
- pressure_solution000007.vtu
- pressure_solution000008.vtu
- pressure_solution000009.vtu

## The velocity vector field is saved on each step.

The following files can be found inside the velocity solutions directory
- velocity_solution000000.vtu
- velocity_solution000001.vtu
- velocity_solution000002.vtu
- velocity_solution000003.vtu
- velocity_solution000004.vtu
- velocity_solution000005.vtu
- velocity_solution000006.vtu
- velocity_solution000007.vtu
- velocity_solution000008.vtu
- velocity_solution000009.vtu
