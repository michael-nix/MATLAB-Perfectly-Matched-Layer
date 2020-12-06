[![View Perfectly Matched Layer for a Standard Wave Equation on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/77521-perfectly-matched-layer-for-a-standard-wave-equation)

## Perfectly Matched Layer for a Standard Wave Equation
How do you add decent absorbing boundary conditions so that you can pretend you're simulating real electromagnetic phenomenon except inside of a computer?  How do you do this when you're not solving Maxwell's equations, but wave equations for potentials and not fields?  Well look no further, 'cause some ongoing work is being done all up in this folder.  So far this includes a, "standard," analytic continuation of spatial coordinates into the complex domain, and then discretized and solved using a few different techniques:

- A fully explicit finite difference method using first order equations via an auxiliary differential equation,
- A fully explicit finite difference method using second order equations,
- A semi-implicit finite difference method using first order equations via an auxiliary differential equation.

The nice thing about these methods is that the exact same files should work exactly the same in 3D (tweaked to add in the third dimension, albeit quite slow and memory intensive) because MATLAB is rad like that.  Since my larger research requires absolute stability, hopefully a fully implicit method will appear soon.  This requires a completely different approach to thinking about integrating the wave equation so there will be a quick writeup for it as well.

With the exception of `setupPML.m` each of these files is a standalone file that you should be able to run to see how things play out for a standard oscillating source charge distribution.
