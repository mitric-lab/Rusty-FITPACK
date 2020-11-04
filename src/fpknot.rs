///  The function fpknot locates an additional knot for a spline of degree
///  k and adjusts the corresponding parameters,i.e.
///    t     : the position of the knots.
///    n     : the number of knots.
///    nrint : the number of knot intervals.
///    fpint : the sum of squares of residual right hand sides
///            for each knot interval.
///    nrdata: the number of data points inside each knot interval.
///  istart indicates that the smallest data point at which the new knot
///  may be added is x(istart+1)
/// The original subroutine in FITPACK by Paul Dierckx is named fpknot
fn fpknot(
    x: Vec<f64>,
    t: Vec<f64>,
    n: usize,
    fpint: Vec<f64>,
    nrdata: Vec<usize>,
    nrint: usize,
    istart: usize,
) -> f64 {
    let k: usize = (n - nrint - 1) / 2;
    let mut t: Vec<f64> = t;
    let mut number: usize = 0;
    let mut fpint: Vec<f64> = fpint;
    let mut nrdata: Vec<usize> = nrdata;
    let mut nrint: usize = nrint;
    let mut n: usize = n;
    //  search for knot interval t(number+k) <= x <= t(number+k+1) where
    //  fpint(number) is maximal on the condition that nrdata(number)
    //  not equals zero.
    let mut fpmax: f64 = 0.0;
    let mut jbegin: usize = istart;
    let mut maxpt: usize = 0;
    let mut maxbeg: usize = 0;
    for j in 1..(nrint + 1) {
        let jpoint: usize = nrdata[j - 1];
        if fpmax < fpint[j - 1] && jpoint != 0 {
            fpmax = fpint[j - 1];
            number = j;
            maxpt = jpoint;
            maxbeg = jbegin;
        }
        jbegin = jbegin + jpoint + 1;
    }
    //  let coincide the new knot t(number+k+1) with a data point x(nrx)
    //  inside the old knot interval t(number+k) <= x <= t(number+k+1).
    let ihalf: usize = maxpt / 2 + 1;
    let nrx: usize = maxbeg + ihalf;
    let next: usize = number + 1;
    if next <= nrint {
        // adjust different parameters
        for j in next..nrint {
            let jj: usize = next + nrint - j;
            fpint[jj] = fpint[jj - 1];
            nrdata[jj] = nrdata[jj - 1];
            let jk: usize = jj + k;
            t[jk] = t[jk - 1];
        }
    }
    nrdata[number - 1] = ihalf - 1;
    nrdata[next - 1] = maxpt - ihalf;
    let am = maxpt;
    let an = nrdata[number - 1];
    fpint[number] = fpmax * an as f64 / am as f64;
    let an = nrdata[next - 1];
    fpint[next - 1] = fpmax * an as f64 / am as f64;
    let jk = next + k;
    t[jk - 1] = x[nrx - 1];
    n = n + 1;
    nrint = nrint + 1;
    return 1.0;
}
