Rusty Dierckx
===

Rusty Dierckx provide the 1D routines for spline interpolation and fitting 
in Rust. The functions are translated from the original Fortran77 code [FITPACK](http://www.netlib.org/dierckx) by Paul Dierckx.
This packages provides almost the same interface as the [SciPy](http://www.scipy.org) wrapper for FITPACK. 
In concrete terms, the package implements the two functions, `splrep` and `splev`.
 
 
References
===

Paul Dierckx, *Curve and Surface Fitting with Splines*, Oxford University Press, 1993.