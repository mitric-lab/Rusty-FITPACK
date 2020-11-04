///  The function check_knots verifies the number and the position of the knots
///  t[j],j=0,1,2,...,n-1 of a spline of degree k, in relation to the number
///  and the position of the data points x[i],i=0,1,2,...,m-1. if all of the
///  following conditions are fulfilled, the function will return zero.
///  if one of the conditions is violated the function returns 1.
///      1) k+1 <= n-k-1 <= m
///      2) t[0] <= t[1] <= ... <= t[k]
///         t[n-k] <= t[n-k+1] <= ... <= t[n-1]
///      3) t[k] < t[k+1] < ... < t[n-k-1]
///      4) t[k] <= x[i] <= t[n-k-1]
///      5) the conditions specified by schoenberg and whitney must hold
///         for at least one subset of data points, i.e. there must be a
///         subset of data points y(j) such that
///             t[j] < y[j] < t[j+k+1], j=0,1,2,...,n-k-1
///
///  The original subroutine in FITPACK by Paul Dierckx is named fpchec
pub(crate) fn check_knots(x: &Vec<f64>, t: &Vec<f64>, k: usize, m: usize, n: usize) -> u8 {
    let k1: usize = k + 1;
    let k2: usize = k1 + 1;
    let nk1: usize = n - k1;
    let nk2: usize = nk1 + 1;
    // check condition no 1
    if nk1 < k1 || nk1 > m {
        panic!(
            "condition 1 is not satisfied. nk1: {}, k1: {}, m: {}",
            nk1, k1, m
        )
    }
    // check condition no 2
    let mut j: usize = n;
    for i in 1..(k + 1) {
        if t[(i - 1) as usize] > t[i as usize] {
            panic!(
                "condition 2 is not satisfied. t[i-1]: {}, t[i]: {}",
                t[i - 1],
                t[i]
            );
        }
        if t[j - 1] < t[j - 2] {
            panic!(
                "condition 2 is not satisfied. t[j-1]: {}, t[j-2]: {}",
                t[j - 1],
                t[j - 2]
            );
        }
        j = j - 1;
    }
    // check condition no 3
    for i in k1..(nk1 + 1) {
        if t[i - 1] <= t[i - 2] {
            panic!(
                "condition 3 is not satisfied. t[i-1]: {}, t[i-2]: {}",
                t[i - 1],
                t[i - 2]
            );
        }
    }
    // check condition no 4
    if x[0] < t[k1 - 1] || x[m - 1] > t[nk2 - 1] {
        panic!(
            "condition 4 is not satisfied. x[0]: {}, t[k1-1]: {}, x[m-1]: {}, t[nk2-1]: {}",
            x[0],
            t[k1 - 1],
            x[m - 1],
            t[nk2 - 1]
        );
    }
    //check condition no 5
    if x[0] > t[k2 - 1] || x[m - 1] <= t[nk1 - 1] {
        panic!(
            "condition 5 is not satisfied. x[0]: {}, t[k2-1]: {}, x[m-1]: {}, t[nk1-1]: {}",
            x[0],
            t[k2 - 1],
            x[m - 1],
            t[nk1 - 1]
        );
    }
    let mut i: usize = 1;
    let mut l: usize = k2;
    let nk3: usize = nk1 - 1;
    if nk3 < 2 {
        return 0;
    }
    for j in 2..(nk3 + 1) {
        l = l + 1;
        i = i + 1;
        if i > m {
            panic!("check failed after five conditions. i: {}, m: {}", i, m);
        }
        while x[i - 1] <= t[j - 1] {
            i = i + 1;
            if i > m {
                panic!("check failed after five conditions. i: {}, m: {}", i, m);
            }
        }
        if x[i - 1] >= t[l - 1] {
            panic!(
                "check failed after five conditions. x[i - 1]: {}, t[l - 1]: {}",
                x[i - 1],
                t[l - 1]
            );
        }
    }
    return 0;
}
