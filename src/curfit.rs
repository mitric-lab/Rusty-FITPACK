use crate::fpchec::fpchec;
use crate::fpcurf::fpcurf;

pub fn curfit(
    iopt: i8,
    x: Vec<f64>,
    y: Vec<f64>,
    w: Vec<f64>,
    xb: f64,
    xe: f64,
    k: usize,
    s: f64,
    nest: usize,
    t: Vec<f64>,
    wrk: Vec<f64>,
) -> (Vec<f64>, usize, Vec<f64>, f64, i8) {
    let m: usize = x.len();
    let n: usize = nest.clone();
    let lwrk: usize = wrk.len();
    assert!(k <= 5, "the degree k must be within 1 and 5");
    let k1: usize = k + 1;
    let k2: usize = k1 + 1;
    assert!(iopt >= -1 && iopt <= 1, "iopt must be within -1 and 1");
    let nmin: usize = 2 * k1 as usize;
    assert!(
        m >= k1 as usize && nest >= nmin,
        "number of data points (m) must be greater than k"
    );
    let lwest: usize = m * k1 as usize + nest * (7 + 3 * k as usize);
    assert!(lwrk >= lwest, "lwrk is to small");
    assert!(xb <= x[0] && xe >= x[m - 1]);
    if iopt >= 0 {
        assert!(s >= 0.0, "smoothing factor cannot be negative");
        if s == 0.0 && nest <= m + k {
            panic!(
                "s: {}, nest: {}, m+k+1: {}, nest cannot be smaller than m+k+1",
                s,
                nest,
                m + k1
            )
        }
    }
    assert!(
        n >= nmin && n <= nest,
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
    let ier: u8 = 0; //
    //fpchec::check_knots(&x, &t, k, m, n);
    assert_eq!(ier, 0, "The knots dont fullfill all five conditions");

    let (t, n, c, fp, ier) = fpcurf(iopt, x, y, w, m, xb, xe, k, s, nest, k1, k2, n, t);
    return (t, n, c, fp, ier);
}
