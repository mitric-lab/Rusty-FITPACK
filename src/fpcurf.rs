use crate::{fpback, fpbspl, fpdisc, fpgivs, fprati, fprota};
use ndarray::prelude::*;
use ndarray::Array2;
use std::cmp::{max, min};

pub fn fpcurf(
    iopt: i8,
    x: Vec<f64>,
    y: Vec<f64>,
    w: Vec<f64>,
    m: usize,
    xb: f64,
    xe: f64,
    k: usize,
    s: f64,
    nest: usize,
    k1: usize,
    k2: usize,
    n: usize,
    t: Vec<f64>,
) -> (Vec<f64>, usize, Vec<f64>, f64, i8) {
    let nmin: usize = 2 * k as usize + 2;

    let tol: f64 = 1e-3;
    let acc: f64 = tol * s;
    let mut n: usize = n;
    let mut ier: i8 = 0;
    let mut t: Vec<f64> = t;
    let maxit: usize = 20;
    let mut nrdata: Vec<usize> = vec![0; nest];
    let mut fpint: Vec<f64> = vec![0.0; nest];
    let mut c: Vec<f64> = vec![0.0; nest];
    let mut z: Vec<f64> = vec![0.0; nest];
    let nmax: usize = m + k1;
    let mut h: Vec<f64> = vec![0.0; 20];
    let mut fpms: f64 = 0.0;
    let nk1: usize;
    let mut yi: f64;
    let mut i2: usize;
    let mut fp: f64 = 0.0;

    let mut finished: bool = false;

    let mut piv: f64;
    let mut p: f64 = 0.0;
    let mut p1: f64 = 0.0;
    let mut p2: f64;
    let mut p3: f64;
    let mut f2: f64;
    let mut ich3: usize;

    // create arrays with Fortran (column-major) memory layout filled with zeros
    let mut a: Array2<f64> = Array2::zeros((nest, k2).f());
    let mut b: Array2<f64> = Array2::zeros((nest, k2).f());
    let mut g: Array2<f64> = Array2::zeros((nest, k2).f());
    let mut q: Array2<f64> = Array2::zeros((m, k1).f());

    let mut fpold: f64 = 0.0;
    let mut fp0: f64 = 0.0;
    let mut nplus: usize = 0;

    if s > 0.0 && iopt >= 0 {
        if iopt != 0 || n != nmin {
            fp0 = fpint[n - 1];
            fpold = fpint[n - 2];
            nplus = nrdata[n - 1];
        }
        if fp0 <= s {
            n = nmin;
            fpold = 0.0;
            nplus = 0;
            nrdata[0] = m - 2;
        }
    } else if s == 0.0 && iopt >= 0 {
        n = nmax;
        assert!(
            nmax <= nest,
            "the storage space exceeds available space, try to increase nest"
        );
        // find the position of the interior knots in case of interpolation
        let mk1 = m - k1;
        if mk1 != 0 {
            let mut i: usize = k2;
            let mut j: usize = k / 2 + 2;
            if k % 2 == 0 {
                for _l in 1..(mk1 + 1) {
                    t[i - 1] = (x[j - 1] + x[j - 2]) * 0.5;
                    i = i + 1;
                    j = j + 1;
                }
            } else {
                for _l in 1..(mk1 + 1) {
                    t[i - 1] = x[j - 1];
                    i = i + 1;
                    j = j + 1;
                }
            }
        }
    }
    nk1 = n - k1;
    // main loop for the different sets of knots. m is a save upper bound
    // for the number of trials.
    // loop 200 in original code
    'main_loop: for _iter in 1..(m + 1) {
        if n == nmin {
            ier = -2;
        };
        // find nrint, tne number of knot intervals.
        let nrint: usize = n - nmin + 1;
        // find the position of the additional knots which are needed for
        // the b-spline representation of s[x].
        let mut i: usize = n;
        for j in 1..(k1 + 1) {
            t[j - 1] = xb;
            t[i - 1] = xe;
            i = i - 1;
        }
        //  compute the b-spline coefficients of the least-squares spline
        //  sinf(x). the observation matrix a is built up row by row and
        //  reduced to upper triangular form by Givens transformations.
        //  at the same time fp=f(p=inf) is computed.
        // initialize the observation matrix a
        for i in 1..(nk1 + 1) {
            z[i - 1] = 0.0;
            // TODO: we dont need this?
            for j in 1..(k1 + 1) {
                a[[i - 1, j - 1]] = 0.0;
            }
        }
        let mut l: usize = k1;
        // loop 130 in original code
        for it in 1..(m + 1) {
            // fetch the current data point x[it], y[it]
            let xi: f64 = x[it - 1];
            let wi: f64 = w[it - 1];

            yi = y[it - 1] * wi;

            // search for knot interval t[l] <= xi < t[l+1].
            while xi >= t[l] && l != nk1 {
                l = l + 1;
            }
            // evaluate the (k+1) non-zero b-splines at xi and store them in q
            h = fpbspl::fpbspl(xi, &t, k, l);
            for i in 1..(k1 + 1) {
                q[[it - 1, i - 1]] = h[i - 1];
                h[i - 1] = h[i - 1] * wi;
            }

            let mut j: i32 = l as i32 - k1 as i32;
            // loop 110 in original code
            'baby_loop: for i in 1..(k1 + 1) {
                j = j + 1;
                piv = h[i - 1];
                if piv != 0.0 {
                    // calculate the parameters of the Givens transformation
                    let (ww, cos, sin): (f64, f64, f64) =
                        fpgivs::fpgivs(piv, a[[j as usize - 1, 0]]);
                    a[[j as usize - 1, 0]] = ww;
                    // transformations to right hand side.
                    let (tmp_a, tmp_b): (f64, f64) =
                        fprota::fprota(cos, sin, yi, z[j as usize - 1]);
                    yi = tmp_a;
                    z[j as usize - 1] = tmp_b;
                    if i == k1 {
                        break 'baby_loop;
                    }
                    i2 = 1;
                    let i3: usize = i + 1;
                    for i1 in i3..(k1 + 1) {
                        i2 = i2 + 1;
                        // transformations to left hand side
                        let (tmp_a, tmp_b): (f64, f64) =
                            fprota::fprota(cos, sin, h[i1 - 1], a[[j as usize - 1, i2 - 1]]);
                        h[i1 - 1] = tmp_a;
                        a[[j as usize - 1, i2 - 1]] = tmp_b;
                    }
                }
            }
            //  add contribution of this row to the sum of squares of residual
            //  right hand sides.
            fp = fp + yi * yi;
        }

        if ier == -2 {
            fp0 = fp;
        }

        fpint[n - 1] = fp0;
        fpint[n - 2] = fpold;
        nrdata[n - 1] = nplus;
        //  backward substitution to obtain the b-spline coefficients.
        c = fpback::fpback(a.view(), z.clone(), nk1, k1, c.clone());

        //  test whether the approximation sinf(x) is an acceptable solution .
        fpms = fp - s;
        if iopt < 0 || fpms.abs() < acc || n == nmax || fpms < 0.0 {
            finished = true;
            break 'main_loop;
        }

        // if f(p=inf) < s accept the choice of knots
        if fpms >= 0.0 {
            //  if n = nmax, sinf(x) is an interpolating spline.
            // if(n.eq.nmax) go to 430
            // increase the number of knots.
            // if n=nest we cannot increase the number of knots because of
            //  the storage capacity limitation.
            if ier != 0 {
                nplus = 1;
                ier = 0;
            } else {
                let mut npl1: f64 = (nplus * 2) as f64;
                let rn: f64 = nplus as f64;
                if fpold - fp > acc {
                    npl1 = rn * fpms / (fpold - fp);
                    //;
                }
                nplus = min(nplus * 2, max(max(npl1 as usize, nplus / 2), 1 as usize));
            }
            fpold = fp;
            //  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
            //  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
            let mut fpart: f64 = 0.0;
            i = 1;
            l = k2;
            let mut new: usize = 0;
            for it in 1..(m + 1) {
                if x[it - 1] >= t[l - 1] && l <= nk1 {
                    new = 1;
                    l = l + 1;
                }
                let mut term: f64 = 0.0;
                let mut l0 = l - k2;
                for j in 1..(k1 + 1) {
                    l0 = l0 + 1;
                    term = term + c[l0 - 1] * q[[it - 1, j - 1]];
                }
                term = (w[it - 1] * (term - y[it - 1])).powi(2);
                fpart = fpart + term;
                if new != 0 {
                    let store = term * 0.50;
                    fpint[i - 1] = fpart - store;
                    i = i + 1;
                    fpart = store;
                    new = 0;
                }
                fpint[nrint - 1] = fpart;
                for _l in 1..(nplus + 1) {
                    println!("WE SHOULD ADD A KNOT");
                    // TODO: call of fpknot is missing
                    //  add a new knot.
                    // call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
                    //  if n=nmax we locate the knots as for interpolation.
                    // if(n.eq.nmax) go to 10
                    //  test whether we cannot further increase the number of knots.
                    // if(n.eq.nest) go to 200
                }
            }
        }
    }
    if ier == -2 {
        finished = true;
    }
    // PART 2
    // initial value for p
    b = fpdisc::fpdisc(t.clone(), k2, n, b);
    let mut f1: f64 = fp0 - s;
    let mut f3: f64 = fpms;
    p3 = -1.0;
    for i in 1..(nk1 + 1) {
        p = p + a[[i - 1, 0]]
    }
    let rn: usize = nk1;
    p = rn as f64 / p;
    let mut ich1: usize = 0;
    ich3 = 0;
    let n8: usize = n - nmin;
    if finished == false {
        // iteration process to find the root of f[p] = s
        'outer: for iter in 1..(maxit + 1) {
            //  the rows of matrix b with weight 1/p are rotated into the
            //  triangularised observation matrix a which is stored in g.
            let pinv: f64 = 1.0 / p;
            for i in 1..(nk1 + 1) {
                c[i - 1] = z[i - 1];
                g[[i - 1, k2 - 1]] = 0.0;
                for j in 1..(k1 + 1) {
                    g[[i - 1, j - 1]] = a[[i - 1, j - 1]];
                }
                g = a.clone();
                for it in 1..(n8 + 1) {
                    for i in 1..(k2 + 1) {
                        h[i - 1] = b[[it - 1, i - 1]] * pinv;
                    }
                    yi = 0.0;
                    for j in it..(nk1 + 1) {
                        piv = h[0];
                        //  calculate the parameters of the Givens transformation.
                        //  call fpgivs(piv,g(j,1),cos,sin)
                        let (ww, cos, sin): (f64, f64, f64) =
                            fpgivs::fpgivs(piv, g[[j as usize - 1, 0]]);
                        g[[j as usize - 1, 0]] = ww;
                        //  transformations to right hand side.
                        //  call fprota(cos,sin,yi,c(j))
                        let (tmp_a, tmp_b): (f64, f64) =
                            fprota::fprota(cos, sin, yi, c[j as usize - 1]);
                        yi = tmp_a;
                        c[j as usize - 1] = tmp_b;
                        if j != nk1 {
                            if j > n8 {
                                i2 = nk1 - j;
                            } else {
                                i2 = k1;
                            }
                            for i in 1..(i2 + 1) {
                                // transformations to the left hand side
                                let i1 = i + 1;
                                let (tmp_a, tmp_b): (f64, f64) = fprota::fprota(
                                    cos,
                                    sin,
                                    h[i1 - 1],
                                    g[[j as usize - 1, i1 - 1]],
                                );
                                h[i1 - 1] = tmp_a;
                                g[[j as usize - 1, i1 - 1]] = tmp_b;
                                h[i - 1] = h[i1 - 1];
                            }
                            h[i2] = 0.0;
                        }
                    }
                }
                // backward substitution to obtain the b-spline coefficients
                // computation of f[p]
                c = fpback::fpback(g.view(), c.clone(), nk1, k1, c.clone());
                let mut l = k2;
                for it in 1..(m + 1) {
                    if x[it - 1] >= t[l - 1] && l <= nk1 {
                        l = l + 1;
                    }
                    let mut l0 = l - k2;
                    let mut term = 0.0;
                    for j in 1..(k1 + 1) {
                        l0 = l0 + 1;
                        term = term + c[l0 - 1] * q[[it - 1, j - 1]];
                    }
                    fp = fp + (w[it - 1] * (term - y[it - 1])).powi(2);
                }
                // test whether the approximation sp(x) is an acceptable solution
                // or the maximal number of iterations is reached
                fpms = fp - s;
                if fpms < acc {
                    break 'outer;
                }
                assert_ne!(
                iter, maxit,
                "Error. The maximal number of iterations maxit (set to 20 by the program) allowed
                 for finding a smoothing spline with fp=s has been reached.
                 Probably cause: s is too small."
            );
                // carry out one more step of the iteration process
                p2 = p;
                f2 = fpms;
                if ich3 == 0 {
                    // our initial choice of p is too large
                    if (f2 - f3) <= acc {
                        p3 = p2;
                        f3 = f2;
                        p = p * 0.4;
                        if p <= p1 {
                            p = p1 * 0.9 + p2 * 0.1;
                        }
                        continue;
                    } else {
                        if f2 < 0.0 {
                            ich3 = 1;
                        }
                        // our initial choice of p is too small
                        if ich1 == 0 && (f1 - f2) <= acc {
                            p1 = p2;
                            f1 = f2;
                            p = p * 0.4;
                            if p >= p3 {
                                //p = p2 * 0.1 + p3 * 0.9
                            }
                            if f2 > 0.0 {
                                ich1 = 1
                            }
                        }
                    }
                }
                //test whether the iteration process proceeds as theoretically
                //expected
                if f2 >= f1 || f2 <= f3 {
                    panic!(
                        "Error. A theoretically impossible result was found during the iteration process
                         for finding a smoothing spline with fp = s.
                         Probably cause: s is too small."
                    )
                }
                // find the new value for p
                let (tmp_p, tmp_p1, tmp_f1, tmp_p3, tmp_f3): (f64, f64, f64, f64, f64) =
                    fprati::fprati(p1, f1, p2, f2, p3, f3);
                p = tmp_p;
                p1 = tmp_p1;
                f1 = tmp_f1;
                p3 = tmp_p3;
                f3 = tmp_f3;
                // error codes and messages
            }
        }
    }
    return (t, n, c, fp, ier);
}
