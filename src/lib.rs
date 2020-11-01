use std::cmp::max;

mod curfit;
mod fpback;
mod fpbspl;
mod fpchec;
mod fpcurf;
mod fpdisc;
mod fpgivs;
mod fpknot;
mod fprati;
mod fprota;

pub fn splrep(
    x: Vec<f64>,
    y: Vec<f64>,
    w: Option<Vec<f64>>,
    xb: Option<f64>,
    xe: Option<f64>,
    k: Option<usize>,
    task: Option<i8>,
    s: Option<f64>,
    t: Option<Vec<f64>>,
    full_output: Option<bool>,
    per: Option<bool>,
    quiet: Option<bool>,
    //) -> (Vec<f64>, Vec<f64>, usize) {
) -> () {
    let m: usize = x.len();
    let s: f64 = match s {
        None => {
            if w == None {
                0.0
            } else {
                (m as f64) - (2.0 * m as f64).sqrt()
            }
        }
        Some(value) => value,
    };
    let w: Vec<f64> = w.unwrap_or(vec![1.0; m]);
    let k: usize = k.unwrap_or(3);
    let mut task: i8 = task.unwrap_or(0);
    let full_output: bool = full_output.unwrap_or(false);
    let per: bool = per.unwrap_or(false);
    let quiet: bool = quiet.unwrap_or(true);

    assert_eq!(w.len(), m, "length of w is not equal to length of x");
    assert_eq!(x.len(), y.len(), "length of x is not equal to length of y");
    assert!(
        1 <= k && k <= 5,
        "Given degree of the spline is not supported (1<=k<=5)."
    );
    assert!(m > k, "m > must hold");
    let xb: f64 = xb.unwrap_or(x[0]);
    let xe: f64 = xe.unwrap_or(x[-1]);
    assert!(0 <= task && task <= 1, "task must be 0, or 1");
    if t == Some(T) {
        task = -1;
    }
    let (t, nest): (Vec<f64>, usize) = if task == -1 {
        assert_eq!(t, Some(T), "knots must be given for task = -1");
        let numknots: usize = t.unwrap().len();
        let nest: usize = numknots + 2 * k + 2;
        let mut new_t: Vec<f64> = vec![0.0; numknots];
        for (i, value) in t.unwrap().iter().enumerate() {
            if i < (numknots - k - 2) {
                new_t[k + 1 + i] = *value;
            }
        }
        (new_t, nest)
    } else if task == 0 {
        (vec![0.0; nest], max(m + k + 1, 2 * k + 3))
    };
    let wrk: Vec<f64> = vec![0.0; m * (k + 1) + nest * (7 + 3 * k)];
    let iwrk: Vec<i32> = vec![0; nest];

    let (n, c, fp, ier) = curfit::curfit(task, x, y, w, xb, xe, k, s, nest, t, wrk, iwrk);
    // curfit call here

    //let tck = (t[..n], c[..n], k);
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
