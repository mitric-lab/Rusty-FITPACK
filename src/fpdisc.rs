///  The function fpdisc calculates the discontinuity jumps of the kth
///  derivative of the b-splines of degree k at the knots t[k+2]..t[n-k-1]
///  The original subroutine in FITPACK by Paul Dierckx is named fpdisc
fn fpdisc(t: Vec<f64>, k2: usize, n: usize, b: Array2<f64>) -> f64 {
    let mut b: Array2<f64> = b;

    let k1: usize = k2 - 1;
    let k: usize = k1 - 1;
    let nk1: usize = n - k1;
    let nrint: usize = nk1 - k;
    let an: f64 = nrint as f64;
    let fac: f64 = an / (t[nk1] - t[k1 - 1]);
    for l in k2..(nk1 + 1) {
        let lmk: usize = l - k1;
        for j in 1..(k1 + 1) {
            let ik: usize = j + k1;
            let lj: usize = l + j;
            let lk: usize = lj - k2;
            h[j - 1] = t[l - 1] - t[lk - 1];
            h[ik - 1] = t[l - 1] - t[lj - 1];
        }
        let mut lp: usize = lmk;
        for j in 1..(k2 + 1) {
            let mut jk: usize = j;
            let mut prod = h[j - 1];
            for i in 1..(k + 1) {
                jk = jk + 1;
                prod = prod * h[jk - 1] * fac;
            }
            let lk: usize = lp + k1;
            b[[lmk - 1, j - 1]] = (t[lk - 1] - t[lp - 1]) / prod;
            lp = lp + 1;
        }
    }
    return 1.0;
}
