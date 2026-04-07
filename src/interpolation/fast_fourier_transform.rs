use num_complex::Complex;
use std::f64::consts::{PI, E};

pub fn fast_fourier_transform(
    xs: Vec<Complex<f64>>,
) -> Vec<Complex<f64>> {
    let n = xs.len();
    if n <= 1 {
        return xs
    };
    let Xs_even = fast_fourier_transform(
        xs.iter().step_by(2).cloned().collect()
    );
    let Xs_odd = fast_fourier_transform(
        xs.iter().skip(1) .step_by(2).cloned().collect()
    );

    let mut Xs: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); n];
    let base = Complex::new(E, 0.0);
    for i in 0..n/2 {
        let w = base.powc(Complex::new(0.0, -2.0) * PI * (i as f64)/(n as f64));
        Xs[i] = Xs_even[i] + w * Xs_odd[i];
        Xs[i + n/2] = Xs_even[i] - w * Xs_odd[i];
    }
    Xs
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fast_fourier_transform() {
        let test_xs = vec![
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0)
        ];

        let result = fast_fourier_transform(test_xs);
        println!("{:?}", result);
        assert_eq!(result, vec![
            Complex::new(10.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(-2.0, 0.0),
            Complex::new(0.0, 0.0),
        ])
    }
}