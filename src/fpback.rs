///  Function fpback calculates the solution of the system of
///  equations $a * c = z$ with a a n x n upper triangular matrix
///  of bandwidth k.
fn fpback(a: ArrayView2<f64>, z: Vec<f64>, n: usize, k: u8, c: Vec<f64>) -> Vec<f64> {
    let mut c: Vec<f64> = c;

    let k1: usize = k as usize - 1;
    c[n - 1] = z[n - 1] / a[[n - 1, 0]];
    let mut i: usize = n - 1;
    if i != 0 {
        for j in 2..(n + 1) {
            let mut store: f64 = z[i - 1];
            let mut i1: usize = k1;
            if j <= k1 {
                i1 = j - 1;
            }
            let mut m: usize = i;
            for l in 1..(i1 + 1) {
                m = m + 1;
                store = store - c[m - 1] * a[[i - 1, l]];
            }
            c[i - 1] = store / a[[i - 1, 0]];
            i = i - 1;
        }
    }
    return c;
}