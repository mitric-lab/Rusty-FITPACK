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
) -> (Vec<f64>, Vec<f64>, usize) {
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
    let xe: f64 = xe.unwrap_or(x[x.len() - 1]);
    assert!(0 <= task && task <= 1, "task must be 0, or 1");
    if t.is_some() {
        task = -1;
    }
    let (t, nest): (Vec<f64>, usize) = if task == -1 {
        assert!(t.is_some(), "knots must be given for task = -1");
        let numknots: usize = t.clone().unwrap().len();
        let nest: usize = numknots + 2 * k + 2;
        let mut new_t: Vec<f64> = vec![0.0; numknots];
        for (i, value) in t.unwrap().iter().enumerate() {
            if i < (numknots - k - 2) {
                new_t[k + 1 + i] = *value;
            }
        }
        (new_t, nest)
    } else if task == 0 {
        let nest: usize = max(m + k + 1, 2 * k + 3);
        (vec![0.0; nest], nest)
    } else {
        (vec![0.0; 1], 0)
    };
    let wrk: Vec<f64> = vec![0.0; m * (k + 1) + nest * (7 + 3 * k)];
    let iwrk: Vec<i32> = vec![0; nest];
    println!("len wrk {}", wrk.len());
    println!("nest {} iwrk {:?}", nest, iwrk);
    //assert_ne!(1, 1);
    let (t, n, c, fp, ier): (Vec<f64>, usize, Vec<f64>, f64, i8) =
        curfit::curfit(task, x, y, w, xb, xe, k, s, nest, t.clone(), wrk, iwrk);
    // curfit call here

    let tck = (t[..n].to_vec(), c[..n].to_vec(), k);
    return tck;
}

#[cfg(test)]
mod tests {
    use crate::splrep;

    #[test]
    fn it_works() {
        let x = vec![0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];
        let y = vec![0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 72.0, 81.0, 100.0, 121.0];
        let (a, b, c) = splrep(
            x, y, None, None, None, None, None, None, None, None, None, None,
        );
        println!("{:?}", a);
        println!("{:?}", b);
        println!("{}", c);
        let t_ref: Vec<f64> = vec![0.5, 0.5, 0.5, 0.5, 2.,  3.,  4.,  6.,  6.,  6.,  6.];
        assert_eq!(a, t_ref);
        panic!();
    }
}
