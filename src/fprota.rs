/// Function fprota applies a Givens rotation to a and b.
/// For details on the Givens rotation see for example (Wikipedia)[https://en.wikipedia.org/wiki/Givens_rotation]
/// The original subroutine in FITPACK by Paul Dierckx is named fprota
fn fprota(cos: f64, sin: f64, a: f64, b: f64) -> (f64, f64) {
    let stor1: f64 = a;
    let stor2: f64 = b;
    let b: f64 = cos * stor2 + sin * stor1;
    let a: f64 = cos * stor1 - sin * stor2;

    return (a, b);
}
