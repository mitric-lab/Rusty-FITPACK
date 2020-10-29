/// given three points (p1, f1), (p2, f2) and (p3, f3), function fprati
/// gives the value of p such that the rational interpolating function
/// of the form r(p) = (u*p+v)/(p+w) equals zero at p.
/// The original subroutine in FITPACK by Paul Dierckx is named fprati
fn fprati(
    mut p1: f64,
    mut f1: f64,
    p2: f64,
    f2: f64,
    mut p3: f64,
    mut f3: f64,
) -> (f64, f64, f64, f64, f64) {
    let p: f64 = if p3 <= 0.0 {
        (p1 * (f1 - f3) * f2 - p2 * (f2 - f3) * f1) / ((f1 - f2) * f3)
    } else {
        let h1: f64 = f1 * (f2 - f3);
        let h2: f64 = f2 * (f3 - f1);
        let h3: f64 = f3 * (f1 - f2);
        -(p1 * p2 * h3 + p2 * p3 * h1 + p3 * p1 * h2) / (p1 * h1 + p2 * h2 + p3 * h3)
    };
    // adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0
    if f2 >= 0.0 {
        p1 = p2;
        f1 = f2;
    } else {
        p3 = p2;
        f3 = f2;
    }
    return (p, p1, f1, p3, f3);
}
