use crate::fpchec;
use crate::fpcurf;

pub fn curfit(
    iopt: i8,
    m: usize,
    x: Vec<f64>,
    y: Vec<f64>,
    w: Vec<f64>,
    xb: f64,
    xe: f64,
    k: u8,
    s: f64,
    nest: usize,
    n: usize,
    t: Vec<f64>,
    c: Vec<f64>,
    fp: f64,
    wrk: Vec<f64>,
    lwrk: usize,
    iwrk: Vec<usize>,
    ier: i8,
) -> Option<f64> {
    let maxit: usize = 20;
    assert!(k >= 0 && k <= 5, "the degree k must be within 1 and 5");
    let k1: u8 = k + 1;
    let k2: u8 = k1 + 1;
    assert!(iopt >= -1 && iopt <= 1, "iopt must be within -1 and 1");
    let nmin: usize = 2 * k1 as usize;
    assert!(
        m >= k1 as usize && nest >= nmin,
        "number of data points (m) must be greater than k"
    );
    let lwest: usize = m * k1 as usize + nest * (7 + 3 * k as usize);
    assert!(lwrk >= lwest, "lwrk is to small");
    assert!(xb <= x[0] && xe >= x[m - 1]);
    //assert!(x.is_sorted(), "x has to be sorted in ascending order");
    if iopt >= 0 {
        assert!(s >= 0.0, "smoothing factor cannot be negative");
        if s == 0.0 && nest <= m + k1 {
            panic!()
        }
    }
    assert!(
        n >= nmin && ng <= nest,
        "total number of knots must be within nmin and nest"
    );

    let mut j: usize = n;
    let mut t: Vec<f64> = t;
    for i in 1..(k1 + 1) {
        t[(i - 1) as usize] = xb;
        t[j - 1] = xe;
        j = j - 1;
    }

    // verify the number and position of knots
    let mut ier: u8 = fpchec::check_knots(&x, &t, k, m, n);
    assert_eq!(ier, 0, "The knots dont fullfill all five conditions");

    // we partition the working space and determine the spline approximation
    let ifp: usize = 0;
    let iz: usize = ifp + nest;
    let ia: usize = iz + nest;
    let ib: usize = ia + nest * k1;
    let ig: usize = ib + nest * k2;
    let iq: usize = ig + nest * k2;

    fpcurf::fpcurf(
        iopt, x, y, w, m, xb, xe, k, s, nest, tol, maxit, k1, k2, n, t, c, fp, wrk[ifp-1], wrk[iz-1],
        wrk[ia-1], wrk[ib-1], wrk[ig-1], wrk[iq-1], iwrk, ier,
    );
    return Some(1.0);
}
