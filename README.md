# Shape optimisation of horn

This notebooks are experimental work to get optimal shape for a horn.

# Notebooks

1. horn1-mesh

Let's generate a mesh first with [gmsh](http://gmsh.info/).

2. horn1-fem

Let's model the Helmotz equation with some boundary conditions.
I use [firedrake](https://www.firedrakeproject.org) to solve the equation.

3. horn1-optim

Let's use some shape optimisation. I use [fireshape](https://link.springer.com/article/10.1007/s00158-020-02813-y) as the engine.

# References

## Horn Introduction

- Horn Theory: An Introduction, Part 1 & 2 By Bjørn Kolbrek. Audio Xpress. 2008.

## Shape Optimisation

- Shape optimization of an acoustic horn. E. Bängtsson et al. Computational Methods in Apply Mechanics and Engineering. 192 (2003) 1533–1571.
- Optimization of an acoustic horn with respect to efficiency and directivity. Udawalpola R., Berggren M. (2008) Int. J. Numer. Meth. Engng., 73(11):1571-1606, doi:10.1002/nme.2132.
- Shape and topology op- timization of an acoustic horn-lens combination. Wadbro E., Udawalpola R., Berggren M. (2009) Published online in J. Comput. Appl. Math., doi:10.1016/j.cam.2009.08.028.
- Optimization of vari- able mouth acoustic horn. Submitted to Int. J. Numer. Meth. Engng. Udawalpola R., Wadbro E., Berggren M. (2010) 

## FEM software

- Fireshape: a shape optimization toolbox for Firedrake. Paganini, A., Wechsung, F. Struct Multidisc Optim 63, 2553–2569 (2021).
- Automated shape differentiation in the Unified Form Language. David A. Ham1 · Lawrence Mitchell2 · Alberto Paganini3 · Florian Wechsung.

## Math

- On a class of preconditioners for solving the Helmholtz equation. Y.A. Erlangga C. Vuik C.W. Oosterlee. Department of Applied Mathematical Analysis, Delft University of Technology, Mekelweg 4, 2628 CD, Delft, The Netherlands. 2004.
- How large a shift is needed in the shifted Helmholtz preconditioner for its effective inversion by multigrid? Pierre-Henri Cocquet and Martin J. Gander. SIAM J. SCI. COMPUT. c 2017 Society for Industrial and Applied Mathematics. Vol. 39, No. 2, pp. A438–A478
- Efficient model reduction of parametrized systems by matrix discrete empirical interpolation. Federico Negri, Andrea Manzoni, David Amsallem.  MATHICSE Technical Report Nr. 02.2015 February 2015.

