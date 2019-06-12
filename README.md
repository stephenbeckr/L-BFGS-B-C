# L-BFGS-B-C
L-BFGS-B, converted from Fortran to C with Matlab wrapper

This is a C version of the well-known [L-BFGS-B code](http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html), version 3.0.

It was created with f2c, then hand-coded to remove dependences on the f2c library

There is a Matlab mex wrapper (mex files and .m files, with example). This was the main
motivation for converting to C, since compiling C and Fortran from Matlab is a pain,
especially since many standard users don't have a Fortran compiler (especially for Windows).

This is an update of my previous wrapper that was on the [Mathworks file-exchange](http://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb--l-bfgs-b--mex-wrapper) from 2012--2015.  This code is completely re-done. I no longer post on the mathworks file-exchange since they have a restrictive license that, e.g., prevents one from using that code with a Matlab alternative such as [Octave](https://www.gnu.org/software/octave/). This current code has not been tested under Octave but there is no reason why it should not be able to work without major modification.


More info on the algorithm is available at the [L-BFGS-B wikipedia page](http://en.wikipedia.org/wiki/L-BFGS-B:_Optimization_subject_to_simple_bounds). References for the algorithm:

* R. H. Byrd, P. Lu and J. Nocedal. A Limited Memory Algorithm for Bound Constrained Optimization, (1995), SIAM Journal on Scientific and Statistical Computing , 16, 5, pp. 1190-1208.
* C. Zhu, R. H. Byrd and J. Nocedal. L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization (1997), ACM Transactions on Mathematical Software, Vol 23, Num. 4, pp. 550 - 560.
* J.L. Morales and J. Nocedal. L-BFGS-B: Remark on Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization (2011), ACM Transactions on Mathematical Software, Vol 38, Num. 1.

There are wrappers for L-BFGS-B in other languages (most link to one of the 2.x versions, but some, like this [R package](http://cran.r-project.org/web/packages/lbfgsb3/), link to the current version 3.0).

# Installation

To use in C, go to the `src/` subdirectory and type `make`. This is unnecessary if you just want to use the Matlab wrapper. For an example of how to call the library from C, see the `driver1.c` file. The included Makefile includes a test of the installation using the problem defined in `driver1.c`.

To use in Matlab, you need to compile the mex files. You can either go to `Matlab/` and type `make` from a shell (assumes you are using `gcc`, and using linux/Mac or Windows+cygwin), or from Matlab, go to the `Matlab/` subdirectory and run `compile_mex.m` which will install the mex using whatever compiler you set up during `mex -setup` (should work for all linux/Mac/Windows). The Matlab install file also runs some basic test to verify the installation.

# License

L-BFGS-B is released under the BSD 3-clause license, and I am releasing this software under the same license. See LICENSE for details

The L-BFGS-B website requests that you cite them. From their website:
"Condition for Use: This software is freely available, but we expect that all publications describing  work using this software , or all commercial products using it, quote at least one of the references given below. This software is released under the "New BSD License" (aka "Modified BSD License" or "3-clause license"). "

It would be nice to cite this website as well since it took a significant amount of work...

# Authors
This C version and Matlab wrapper are written by [Stephen Becker](http://amath.colorado.edu/faculty/becker/), stephen.becker@colorado.edu

The L-BFGS-B algorithm was written in the 1990s (mainly 1994, some revisions 1996) by Ciyou Zhu (in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal)

L-BFGS-B Version 3.0 is an algorithmic update from 2011, with coding changes by J. L. Morales.

