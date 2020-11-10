Rusty FITPACK
===

Rusty FITPACK provides the 1D routines for spline interpolation and fitting 
in Rust. The functions are translated from the original Fortran77 code [FITPACK](http://www.netlib.org/dierckx) by Paul Dierckx.
This packages provides almost the same interface as the [SciPy](http://www.scipy.org) wrapper for FITPACK. 
In concrete terms, the package implements three functions, `splrep`, `splev` and `splev_uniform`.
 

Please read the [Documentation](http://jhoche.de/Rusty-FITPACK/rusty_fitpack/)

References
----------
Based on algorithms described by Paul Dierckx in Ref [1-4]:<br>

[1] P. Dierckx, "An algorithm for smoothing, differentiation and integration of experimental data using spline functions", J.Comp.Appl.Maths 1 (1975) 165-184.

[2] P. Dierckx, "A fast algorithm for smoothing data on a rectangular grid while using spline functions", SIAM J.Numer.Anal. 19 (1982) 1286-1304.

[3] P. Dierckx, "An improved algorithm for curve fitting with spline functions", report tw54, Dept. Computer Science,K.U. Leuven, 1981.

[4] P. Dierckx, "Curve and surface fitting with splines", Monographs on Numerical Analysis, Oxford University Press, 1993.