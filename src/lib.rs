use std::cmp::max;
mod splines;
pub mod curfit;

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

    let (t, n, c, fp, ier): (Vec<f64>, usize, Vec<f64>, f64, i8) =
        splines::curfit(task, x, y, w, xb, xe, k, s, nest, t.clone(), wrk);

    let tck = (t[..n].to_vec(), c[..n].to_vec(), k);
    return tck;
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
            0.0000000000000000e+00
        ];
        assert_eq!(t, t_ref);
        assert_eq!(c, c_ref);
        assert_eq!(k, 3);
    }
}
