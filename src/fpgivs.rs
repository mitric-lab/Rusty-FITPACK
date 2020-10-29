/// Function fpgivs calculates the parameters of a givens transformation
/// The original subroutine in FITPACK by Paul Dierckx is named fpgivs
fn fpgivs(piv: f64, ww: f64) -> (f64, f64, f64) {
    let store: f64 = piv.abs();
    let dd: f64 = if store >= ww {
        store * (1.0 + (ww / piv).powi(2)).sqrt()
    } else {
        ww * (1 + (piv / ww).powi(2)).sqrt()
    };
    let cos: f64 = ww / dd;
    let sin: f64 = piv / dd;
    let ww: f64 = dd;

    return (ww, cos, sin);
}