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