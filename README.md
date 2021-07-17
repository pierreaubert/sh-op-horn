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

1. E. Bängtsson et al. Shape optimization of an acoustic horn. Computational Methods in Apply Mechanics and Engineering. 192 (2003) 1533–1571.
2. Paganini, A., Wechsung, F. Fireshape: a shape optimization toolbox for Firedrake. Struct Multidisc Optim 63, 2553–2569 (2021).
3. Udawalpola R., Berggren M. (2008) Optimization of an acoustic horn with respect to efficiency and directivity. Int. J. Numer. Meth. Engng., 73(11):1571-1606, doi:10.1002/nme.2132.
4. Wadbro E., Udawalpola R., Berggren M. (2009) Shape and topology op- timization of an acoustic horn-lens combination. Published online in J. Comput. Appl. Math., doi:10.1016/j.cam.2009.08.028.
5. Udawalpola R., Wadbro E., Berggren M. (2010) Optimization of vari- able mouth acoustic horn. Submitted to Int. J. Numer. Meth. Engng.
