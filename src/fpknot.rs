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
    mut n: usize,
    fpint: Vec<f64>,
    nrdata: Vec<f64>,
    mut nrint: usize,
    istart: usize,
) -> f64 {
    let k: usize = (n - nrint - 1) / 2;
    //  search for knot interval t(number+k) <= x <= t(number+k+1) where
    //  fpint(number) is maximal on the condition that nrdata(number)
    //  not equals zero.
    let mut fpmax: f64 = 0.0;
    let mut jbegin: usize = istart;
    let mut maxpt: f64 = 0.0;
    let mut maxbeg: usize = 0;
    for j in 1..(nrint + 1) {
        let jpoint = nrdata[j - 1];
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
    let ihalf = maxpt / 2 + 1;
    let nrx = maxbeg + ihalf;
    let next = number + 1;
    if next <= nrint {
        // adjust different parameters
        for j in next..nrint {
            jj = next + nrint - j;
            fpint[jj] = fpint[jj - 1];
            nrdata[jj] = nrdata[jj - 1];
            jk = jj + k;
            t[jk] = t[jk - 1];
        }
    }
    nrdata[number - 1] = ihalf - 1;
    nrdata[next - 1] = maxpt - ihalf;
    am = maxpt;
    an = nrdata[number - 1];
    fpint[number] = fpmax * an / am;
    an = nrdata[next - 1];
    fpint[next - 1] = fpmax * an / am;
    jk = next + k;
    t[jk - 1] = x[nrx - 1];
    n = n + 1;
    nrint = nrint + 1;
    return 1.0;
}
