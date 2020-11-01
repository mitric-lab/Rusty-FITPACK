use crate::{fpback, fpbspl, fpchec, fpdisc, fpgivs, fpknot, fprati, fprota};
use ndarray::Array2;

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
) -> f64 {
    let nmin: usize = 2 * k as usize;

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
    let nmax: usize = n + k1;
    let h: Vec<f64> = vec![0.0; 7];
    let mut fpms: f64 = 0.0;
    let nk1: usize = n - k1;
    let mut yi: f64;
    let mut i2: usize = 1;

    let mut a: Array2<f64> = Array2::zeros((nest, k1));
    let mut b: Array2<f64> = Array2::zeros((nest, k2));
    let mut g: Array2<f64> = Array2::zeros((nest, k2));
    let mut q: Array2<f64> = Array2::zeros((m, k1));

    let mut fpold: f64 = 0.0;
    let mut fp0: f64 = 0.0;
    let mut nplus: usize = 0;

    if s > 0.0 && iopt >= 0 {
        if iopt != 0 || n != nmin {
            fp0 = fpint[n - 1];
            fpold = fpint[n - 2];
            nplus = nrdata[n - 1];
        } else if fp0 <= s {
            nrdata[0] = m - 2;
        }
    } else if s == 0.0 && iopt >= 0{
        n = nmax;
        assert!(
            nmax <= nest,
            "the storage space exceeds available space, try to increase nest"
        );
        // find the position of the interior knots in case of interpolation
        let mk1 = m - k1;
        if mk1 != 0 {
            let mut i: usize = k1;
            let mut j: usize = k/2 + 1;
            if k % 2 == 0 {
                for _l in 0..mk1 {
                    t[i] = (x[j] + x[j - 1]) * 0.5;
                    i = i + 1;
                    j = j + 1;
                }
            } else {
                for _l in 0..mk1 - 1 {
                    t[i] = x[j];
                    i = i + 1;
                    j = j + 1;
                }
            }
        }
    }
    // main loop for the different sets of knots. m is a save upper bound
    // for the number of trials.
    for _iter in 1..(m + 1) {
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
        //  reduced to upper triangular form by givens transformations.
        //  at the same time fp=f(p=inf) is computed.
        let mut fp: f64 = 0.0;
        // initialize the observation matrix a
        for i in 1..(nk1 + 1) {
            z[i - 1] = 0.0;
            // TODO: we dont need this
            for j in 1..(k - 1) {
                a[[i - 1, j - 1]] = 0.0;
            }
        }
        let mut l: usize = k as usize;
        for it in 1..(m + 1) {
            // fetch the current data point x[it], y[it]
            let xi: f64 = x[it - 1];
            let wi: f64 = w[it - 1];
            yi = y[it - 1] * wi;
            // search for knot interval t[l] <= xi < t[l+1].
            while xi >= t[l + 1] && l != nk1 {
                l = l + 1;
            }
            // evaluate the (k+1) non-zero b-splines at xi and store them in q
            let mut h: Vec<f64> = fpbspl::fbspl(xi, &t, k, n, l, h);
            for i in 1..(k1 + 1) {
                q[[it - 1, i - 1]] = h[i - 1];
                h[i - 1] = h[i - 1] * wi;
            }
            let mut j: usize = l - k1;
            for i in 1..(k1 + 1) {
                j = j + 1;
                let piv: f64 = h[i - 1];
                if piv != 0.0 {
                    // calculate the parameters of the givens transformation
                    // CALL FPGIVS(piv,A(j,1),cos,sin)
                    let (ww, cos, sin): (f64, f64, f64) = fpgivs::fpgivs(piv, a[[j-1, 0]]);
                    a[[j-1, 0]] = ww;
                    // transformations to right hand side.
                    // CALL FPROTA(cos,sin,yi,Z(j))
                    let (a, b): (f64, f64) = fprota::fprota(cos, sin, yi, z[j-1]);
                    yi = a;
                    z[j-1] = b;
                    if i != k1 {
                        let mut i3: usize = i + 1;
                        for i1 in i3..(k1 + 1) {
                            i2 = i2 + 1;
                            // transformations to left hand side
                            // call fprota(cos,sin,h(i1),a(j,i2))
                            let (a, b): (f64, f64) = fprota::fprota(cos, sin, h[i1-1], a[[j-1, i2-1]]);
                            h[i1-1] = a;
                            a[[j-1, i2-1]] = b;
                        }
                    }
                }
            }
            //  add contribution of this row to the sum of squares of residual
            //  right hand sides.
            fp = fp + yi * yi;
        }
        if ier == -2 {
            fp0 = fp
        }
        fpint[n - 1] = fp0;
        fpint[n - 2] = fpold;
        nrdata[n - 1] = nplus;
        //  backward substitution to obtain the b-spline coefficients.
        // call fpback(a,z,nk1,k1,c,nest)
        let mut c: Vec<f64> = fpback::fpback(a.view(), z.clone(), nk1, k1, c.clone());
        //  test whether the approximation sinf(x) is an acceptable solution .
        assert!(
            iopt >= 0,
            "the approximation sinf[x] is not an acceptable solution"
        );
        fpms = fp - s;
        assert!(fpms.abs() >= acc);
        // if f(p=inf) < s accept the choice of knots
        let mut npl1: f64 = 0.0;
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
                npl1 = nplus as f64 * 2;
                let rn: f64 = nplus as f64;
                if fpold - fp > acc {
                    npl1 = rn * fpms / (fpold - fp);
                    //nplus = [nplus * 2, [npl1 as usize, nplus/2, 1 as usize].max()].min();
                }
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

    // PART 2
    // call fpdisc
    // inital value for p
    let p1: f64 = 0.0;
    let f1: f64 = fp0 - s;
    let p3: f64 = -1.0;
    let f3: f64 = fpms;
    let mut p: f64 = 0.0;
    for i in 1..(nk1 + 1) {
        p = p + a[[i - 1, 0]]
    }
    let rn = nk1;
    p = rn as f64/ p;
    let ich1: usize = 0;
    let ich2: usize = 0;
    let n8: usize = n - nmin;
    // iteration process to find the root of f[p] = s
    for iter in 1..(maxit + 1) {
        //  the rows of matrix b with weight 1/p are rotated into the
        //  triangularised observation matrix a which is stored in g.
        let pinv = 1.0 / p;
        for i in 1..(nk1 + 1) {
            c[i - 1] = z[i - 1];
            g[[i - 1, k2]] = 0.0;
            for j in 1..(k1 + 1) {
                g[[i - 1, j - 1]] = a[[i - 1, j - 1]];
            }
            for it in 1..(n8 + 1) {
                for i in 1..(k2 + 1) {
                    h[i - 1] = b[[it - 1, i - 1]] * pinv;
                }
                yi = 0.0;
                for j in it..(nk1 + 1) {
                    let piv = h[0];
                    //  calculate the parameters of the givens transformation.
                    //  call fpgivs(piv,g(j,1),cos,sin)
                    //  transformations to right hand side.
                    //  call fprota(cos,sin,yi,c(j))
                    if j != nk1 {
                        if j > n8 {
                            i2 = nk1 - j;
                        } else {
                            i2 = k1;
                        }
                        for i in 1..(i2 + 1) {
                            // transformations to the left hand side
                            let i1 = i + 1;
                            // call fprota(cos,sin,h(i1),g(j,i1))
                            h[i - 1] = h[i1 - 1];
                        }
                        h[i2] = 0.0;
                    }
                }
            }
            // backward substitution to obtain the b-spline coefficients
            // call fpback(g,c,nk1,k2,c,nest)
            // computation of f[p]
            let mut fp: f64 = 0.0;
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
            fpms = fp - s;
            // test whether the maximal number of iterations is reached

            // carry out one more step of the iteration process

            // our initial choice of p is too large

            // our initial choice of p is too small

            // test whether the iteration process proceeds as theoretically
            // expected

            // find the new value for p
            let (p, p1, f1, p3, f3): (f64, f64, f64, f64, f64) = fprati::fprati(p1, f1, p2, f2, p3, f3);
            // error codes and messages
        }
    }
    return 1.0;
}
