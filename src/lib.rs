//! Rusty FITPACK provides the 1D routines for spline interpolation and fitting
//! in Rust. The functions are translated from the original Fortran77 code [FITPACK](http://www.netlib.org/dierckx) by Paul Dierckx.
//! This packages provides almost the same interface as the [SciPy](http://www.scipy.org) wrapper for FITPACK.
//! In concrete terms, the package implements the two functions, `splrep` and `splev` (not yet supported).
//!
//!
//!
//! References
//! ----------
//! Based on algorithms described by Paul Dierckx in Ref [1-4]:<br>
//!
//! [1] P. Dierckx, "An algorithm for smoothing, differentiation and integration of experimental data using spline functions", J.Comp.Appl.Maths 1 (1975) 165-184.
//!
//! [2] P. Dierckx, "A fast algorithm for smoothing data on a rectangular grid while using spline functions", SIAM J.Numer.Anal. 19 (1982) 1286-1304.
//!
//! [3] P. Dierckx, "An improved algorithm for curve fitting with spline functions", report tw54, Dept. Computer Science,K.U. Leuven, 1981.
//!
//! [4] P. Dierckx, "Curve and surface fitting with splines", Monographs on Numerical Analysis, Oxford University Press, 1993.

use std::cmp::max;
mod curfit;
mod fpchec;
mod fpcurf;
mod fpdisc;
mod fpgivs;
//mod fpknot;
mod fpback;
mod fpbspl;
mod fprati;
mod fprota;

/// Find the B-spline representation of a 1-D curve.
/// Given the set of data points ``(x[i], y[i])`` determine a smooth spline
/// approximation of degree k on the interval ``xb <= x <= xe``.
///
/// #### Example
/// Simple example of spline interpolation
/// ```
/// use rusty_fitpack::splrep;
/// let x = vec![0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
/// let y = vec![0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0];
///
/// let (t, c, k) = splrep(x, y, None, None, None, None, None, None, None, None, None, None);
/// ```
///
/// #### Parameters
/// ----------
/// `x, y` : The data points defining a curve `y = f(x)`. <br> <br>
/// `w` : Strictly positive `Vec<f64>` of weights the same length as `x` and `y`.
/// The weights are used in computing the weighted least-squares spline
/// fit. If the errors in the `y` values have standard-deviation given by the
/// vector `d`, then w should be `1/d`. Default is ones(len(x)). <br> <br>
/// `xb, xe` : The interval to fit.  If None, these default to `x[0]` and `x[-1]`
/// respectively. <br> <br>
/// `k` : The degree of the spline fit. It is recommended to use cubic splines.
/// Even values of `k` should be avoided especially with small `s` values.
/// `1 <= k <= 5` <br> <br>
/// `task` : `{0, -1}` If `task==0` find `t` and `c` for a given smoothing factor, `s`. <br>
/// If `task=-1` find the weighted least square spline for a given set of
/// knots, `t`. These should be interior knots as knots on the ends will be
/// added automatically.<br> <br>
/// `s` : A smoothing condition. The amount of smoothness is determined by
/// satisfying the conditions: `sum((w * (y - g)).powi(2),axis=0) <= s` where g(x)
/// is the smoothed interpolation of `(x,y)`. The user can use s to control
/// the tradeoff between closeness and smoothness of fit. Larger s means
/// more smoothing while smaller values of s indicate less smoothing.
/// Recommended values of s depend on the weights, w. If the weights
/// represent the inverse of the standard-deviation of y, then a good s
/// value should be found in the range `(m-(2*m).sqrt(),m+(2*m).sqrt())` where m is
/// the number of datapoints in `x, y,` and `w`. default : `s=m-(2*m).sqrt()` if
/// weights are supplied. `s = 0.0` (interpolating) if no weights are
/// supplied. <br> <br>
/// `t` : The knots needed for `task=-1`. If given then task is automatically set
/// to -1. <br> <br>
/// `full_output` Should be None. Feature is not implemented yet. <br> <br>
/// `per` : Should be None. Periodic spline approximations are not supported yet. <br> <br>
/// `quiet`: Should be None. Feature is not implemented yet. <br> <br>
///

pub fn splrep(
    x: Vec<f64>,
    y: Vec<f64>,
    w: Option<Vec<f64>>,
    xb: Option<f64>,
    xe: Option<f64>,
    k: Option<usize>,
    task: Option<i8>,
    s: Option<f64>,
    t: Option<Vec<f64>>,
    full_output: Option<bool>,
    per: Option<bool>,
    quiet: Option<bool>,
    //) -> (Vec<f64>, Vec<f64>, usize) {
) -> (Vec<f64>, Vec<f64>, usize) {
    let m: usize = x.len();
    let s: f64 = match s {
        None => {
            if w == None {
                0.0
            } else {
                (m as f64) - (2.0 * m as f64).sqrt()
            }
        }
        Some(value) => value,
    };
    let w: Vec<f64> = w.unwrap_or(vec![1.0; m]);
    let k: usize = k.unwrap_or(3);
    let mut task: i8 = task.unwrap_or(0);
    let _full_output: bool = full_output.unwrap_or(false);
    let _per: bool = per.unwrap_or(false);
    let _quiet: bool = quiet.unwrap_or(true);

    assert_eq!(w.len(), m, "length of w is not equal to length of x");
    assert_eq!(x.len(), y.len(), "length of x is not equal to length of y");
    assert!(
        1 <= k && k <= 5,
        "Given degree of the spline is not supported (1<=k<=5)."
    );
    assert!(m > k, "m > must hold");
    let xb: f64 = xb.unwrap_or(x[0]);
    let xe: f64 = xe.unwrap_or(x[x.len() - 1]);
    assert!(0 <= task && task <= 1, "task must be 0, or 1");
    if t.is_some() {
        task = -1;
    }
    let (t, nest): (Vec<f64>, usize) = if task == -1 {
        assert!(t.is_some(), "knots must be given for task = -1");
        let numknots: usize = t.clone().unwrap().len();
        let nest: usize = numknots + 2 * k + 2;
        let mut new_t: Vec<f64> = vec![0.0; nest];
        for (i, value) in t.unwrap().iter().enumerate() {
            new_t[k + 1 + i] = *value;
        }
        (new_t, nest)
    } else if task == 0 {
        let nest: usize = max(m + k + 1, 2 * k + 3);
        (vec![0.0; nest], nest)
    } else {
        (vec![0.0; 1], 0)
    };
    let wrk: Vec<f64> = vec![0.0; m * (k + 1) + nest * (7 + 3 * k)];

    let (t, n, c, _fp, _ier): (Vec<f64>, usize, Vec<f64>, f64, i8) =
        curfit::curfit(task, x, y, w, xb, xe, k, s, nest, t.clone(), wrk);

    let tck = (t[..n].to_vec(), c[..n].to_vec(), k);
    return tck;
}


///  The function splev evaluates a number of points $x(i)$ with $i=1,2,...,m$
///  a spline $s(x)$ of degree $k$, given in its b-spline representation.
///
///  calling sequence:
///     call splev(t,n,c,k,x,y,m,e,ier)
///
///  input parameters:
///    t    : array,length n, which contains the position of the knots.
///    n    : integer, giving the total number of knots of s(x).
///    c    : array,length n, which contains the b-spline coefficients.
///    k    : integer, giving the degree of s(x).
///    x    : array,length m, which contains the points where s(x) must
///           be evaluated.
///    m    : integer, giving the number of points where s(x) must be
///           evaluated.
///    e    : integer, if 0 the spline is extrapolated from the end
///           spans for points not in the support, if 1 the spline
///           evaluates to zero for those points, if 2 ier is set to
///           1 and the subroutine returns, and if 3 the spline evaluates
///           to the value of the nearest boundary point.
///
///  output parameter:
///    y    : array,length m, giving the value of s(x) at the different
///           points.
///    ier  : error flag
///      ier = 0 : normal return
///      ier = 1 : argument out of bounds and e == 2
///      ier =10 : invalid input data (see restrictions)
///
///  restrictions:
///    m >= 1
///--    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
///
///  other subroutines required: fpbspl.
///
///  references :
///    de boor c  : on calculating with b-splines, j. approximation theory
///                 6 (1972) 50-62.
///    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
///                 applics 10 (1972) 134-149.
///    dierckx p. : curve and surface fitting with splines, monographs on
///                 numerical analysis, oxford university press, 1993.
pub fn splev() -> f64 {

    return 1.0;
}


#[cfg(test)]
mod tests {
    use crate::splrep;
    /// The reference values were calculated using the SciPy interface to Fitpack
    /// Python: 3.7.6 (conda, GCC 7.3.0), NumPy: 1.18.1, SciPy: 1.4.1
    /// the following code was used
    /// ```python
    /// from scipy.interpolate import splrep
    /// import numpy as np
    /// x = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    /// y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0]
    ///
    /// spl = splrep(x, y)
    /// np.set_printoptions(16)
    /// print(spl[0])
    /// print(spl[1])
    /// print(spl[2])
    /// ```
    #[test]
    fn simple_spline_interpolation() {
        let x = vec![0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];
        let y = vec![
            0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0,
        ];
        let (t, c, k) = splrep(
            x, y, None, None, None, None, None, None, None, None, None, None,
        );
        let t_ref: Vec<f64> = vec![
            0.5, 0.5, 0.5, 0.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 11.0, 11.0, 11.0,
        ];
        let c_ref: Vec<f64> = vec![
            -2.6611993517399935e-17,
            9.0032063142935170e-01,
            2.7876063631714838e+00,
            8.6767816283663706e+00,
            1.5663956371192338e+01,
            2.4667392886864288e+01,
            3.5666472081350541e+01,
            4.8666718787733515e+01,
            6.3666652767715512e+01,
            8.6333342599300764e+01,
            1.0633332870034963e+02,
            1.2100000000000000e+02,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
        ];
        assert_eq!(t, t_ref);
        assert_eq!(c, c_ref);
        assert_eq!(k, 3);
    }
    /// The reference values were calculated using the SciPy interface to Fitpack
    /// Python: 3.7.6 (conda, GCC 7.3.0), NumPy: 1.18.1, SciPy: 1.4.1
    /// the following code was used
    /// ```python
    /// from scipy.interpolate import splrep
    /// import numpy as np
    /// x = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    /// y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0]
    /// spl = splrep(x, y, s=0.5)
    /// np.set_printoptions(16)
    /// print(spl[0])
    /// print(spl[1])
    /// print(spl[2])
    /// ```
    #[test]
    fn simple_spline_fit() {
        let x = vec![0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];
        let y = vec![
            0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0,
        ];
        let (t, c, k) = splrep(
            x,
            y,
            None,
            None,
            None,
            None,
            None,
            Some(0.5),
            None,
            None,
            None,
            None,
        );
        let t_ref: Vec<f64> = vec![0.5, 0.5, 0.5, 0.5, 11.0, 11.0, 11.0, 11.0];
        let c_ref: Vec<f64> = vec![
            9.6048864019998764e-02,
            3.9719628633913779e+00,
            4.3868995123976298e+01,
            1.2101960437915093e+02,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
        ];
        assert_eq!(t, t_ref);
        assert_eq!(c, c_ref);
        assert_eq!(k, 3);
    }

    /// The reference values were calculated using the SciPy interface to Fitpack
    /// Python: 3.7.6 (conda, GCC 7.3.0), NumPy: 1.18.1, SciPy: 1.4.1
    /// the following code was used
    /// ```python
    /// from scipy.interpolate import splrep
    /// import numpy as np
    /// x = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    /// y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0]
    /// w = [0.5, 1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 0.1, 8.0, 0.2, 3.0, 2.0]
    /// spl = splrep(x, y, w=wm s=0)
    /// np.set_printoptions(16)
    /// print(spl[0])
    /// print(spl[1])
    /// print(spl[2])
    /// ```
    #[test]
    fn spline_interpolation_with_weights() {
        let x = vec![0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];
        let y = vec![
            0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0,
        ];
        let w: Vec<f64> = vec![0.5, 1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 0.1, 8.0, 0.2, 3.0, 2.0];
        let (t, c, k) = splrep(
            x,
            y,
            Some(w),
            None,
            None,
            None,
            None,
            Some(0.0),
            None,
            None,
            None,
            None,
        );
        let t_ref: Vec<f64> = vec![
            0.5, 0.5, 0.5, 0.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 11.0, 11.0, 11.0,
        ];
        let c_ref: Vec<f64> = vec![
            3.5816830069903035e-17,
            9.0032063142935159e-01,
            2.7876063631714838e+00,
            8.6767816283663706e+00,
            1.5663956371192336e+01,
            2.4667392886864263e+01,
            3.5666472081350641e+01,
            4.8666718787733217e+01,
            6.3666652767715597e+01,
            8.6333342599300707e+01,
            1.0633332870034965e+02,
            1.2100000000000000e+02,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
        ];
        assert_eq!(t, t_ref);
        assert_eq!(c, c_ref);
        assert_eq!(k, 3);
    }

    /// The reference values were calculated using the SciPy interface to Fitpack
    /// Python: 3.7.6 (conda, GCC 7.3.0), NumPy: 1.18.1, SciPy: 1.4.1
    /// the following code was used
    /// ```python
    /// from scipy.interpolate import splrep
    /// import numpy as np
    /// x = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    /// y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0]
    /// w = [0.5, 1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 0.1, 8.0, 0.2, 3.0, 2.0]
    /// spl = splrep(x, y, w=w)
    /// np.set_printoptions(16)
    /// print(spl[0])
    /// print(spl[1])
    /// print(spl[2])
    /// ```
    #[test]
    fn spline_fit_with_weights() {
        let x = vec![0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];
        let y = vec![
            0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0,
        ];
        let w: Vec<f64> = vec![0.5, 1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 0.1, 8.0, 0.2, 3.0, 2.0];
        let (t, c, k) = splrep(
            x,
            y,
            Some(w),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        let t_ref: Vec<f64> = vec![0.5, 0.5, 0.5, 0.5, 11., 11., 11., 11.];
        let c_ref: Vec<f64> = vec![
            0.20514750870201512,
            3.79852892458793200,
            43.9786486869686500,
            121.003719055670530,
            0.00000000000000000,
            0.00000000000000000,
            0.00000000000000000,
            0.00000000000000000,
        ];
        assert_eq!(t, t_ref);
        assert_eq!(c, c_ref);
        assert_eq!(k, 3);
    }

    /// The reference values were calculated using the SciPy interface to Fitpack
    /// Python: 3.7.6 (conda, GCC 7.3.0), NumPy: 1.18.1, SciPy: 1.4.1
    /// the following code was used
    /// ```python
    /// from scipy.interpolate import splrep
    /// import numpy as np
    /// x = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    /// y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0]
    /// t = [2.5, 3.5, 4.5]
    /// spl = splrep(x, y, t=t)
    /// np.set_printoptions(16)
    /// print(spl[0])
    /// print(spl[1])
    /// print(spl[2])
    /// ```
    #[test]
    fn spline_interpolation_with_specified_knots() {
        let x = vec![0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];
        let y = vec![
            0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0,
        ];
        let t: Vec<f64> = vec![2.5, 3.5, 4.5];
        let (t, c, k) = splrep(
            x,
            y,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(t),
            None,
            None,
            None,
        );
        let t_ref: Vec<f64> = vec![0.5, 0.5, 0.5, 0.5, 2.5, 3.5, 4.5, 11., 11., 11., 11.];
        let c_ref: Vec<f64> = vec![
            4.3456188984254242e-03,
            1.1369046801251863e+00,
            3.8089766153250060e+00,
            1.1942469486307893e+01,
            3.4554109906393840e+01,
            7.3350195466184161e+01,
            1.2099781223876728e+02,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
        ];
        assert_eq!(t, t_ref);
        assert_eq!(c, c_ref);
        assert_eq!(k, 3);
    }

    /// The reference values were calculated using the SciPy interface to Fitpack
    /// Python: 3.7.6 (conda, GCC 7.3.0), NumPy: 1.18.1, SciPy: 1.4.1
    /// the following code was used
    /// ```python
    /// from scipy.interpolate import splrep, splev
    /// import numpy as np
    /// x = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    /// y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0]
    /// w = [0.5, 1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 0.1, 8.0, 0.2, 3.0, 2.0]
    /// t = [2.5, 3.5, 4.5]
    /// spl = splrep(x, y, w=w, xb=0.0, xe=16.0, t=t, s=0.8)
    /// np.set_printoptions(16)
    /// print(spl[0])
    /// print(spl[1])
    /// print(spl[2])
    /// ```
    #[test]
    fn spline_fit_limits_knots_weights() {
        let x = vec![0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];
        let y = vec![
            0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0,
        ];
        let w: Vec<f64> = vec![0.5, 1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 0.1, 8.0, 0.2, 3.0, 2.0];
        let t: Vec<f64> = vec![2.5, 3.5, 4.5];
        let (t, c, k) = splrep(
            x,
            y,
            Some(w),
            Some(0.0),
            Some(16.0),
            None,
            None,
            None,
            Some(t),
            None,
            None,
            None,
        );
        let t_ref: Vec<f64> = vec![0., 0., 0., 0., 2.5, 3.5, 4.5, 16., 16., 16., 16.];
        let c_ref: Vec<f64> = vec![
            -0.6664597168700137,
            0.3028630978060425,
            2.7945126095902090,
            11.928341312803433,
            47.885102645361140,
            133.39937455474050,
            255.87891800341376,
            0.0000000000000000,
            0.0000000000000000,
            0.0000000000000000,
            0.0000000000000000,
        ];
        assert_eq!(t, t_ref);
        assert_eq!(c, c_ref);
        assert_eq!(k, 3);
    }
    /// The reference values were calculated using the SciPy interface to Fitpack
    /// Python: 3.7.6 (conda, GCC 7.3.0), NumPy: 1.18.1, SciPy: 1.4.1
    /// the following code was used
    /// ```python
    /// from scipy.interpolate import splrep
    /// import numpy as np
    /// x = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    /// y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0]
    ///
    /// spl = splrep(x, y, k=5)
    /// np.set_printoptions(16)
    /// print(spl[0])
    /// print(spl[1])
    /// print(spl[2])
    /// ```
    #[test]
    fn spline_interpolation_fifth_order() {
        let x = vec![0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];
        let y = vec![
            0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0,
        ];
        let (t, c, k) = splrep(
            x,
            y,
            None,
            None,
            None,
            Some(5),
            None,
            None,
            None,
            None,
            None,
            None,
        );
        let t_ref: Vec<f64> = vec![
            0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 3., 4., 5., 6., 7., 8., 11., 11., 11., 11., 11., 11.,
        ];
        let c_ref: Vec<f64> = vec![
            -3.1525642173709244e-17,
            9.6187233587515120e-01,
            2.2151448816620882e+00,
            5.9694686107264907e+00,
            1.2784814311098723e+01,
            2.4504570963883427e+01,
            3.5498043000268858e+01,
            5.3701926216065736e+01,
            7.2898267193716578e+01,
            9.1401339513703903e+01,
            1.0779928004881573e+02,
            1.2100000000000000e+02,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
        ];
        assert_eq!(t, t_ref);
        assert_eq!(c, c_ref);
        assert_eq!(k, 5);
    }
}
