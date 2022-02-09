use ark_ff::{PrimeField, Zero};
use crate::Poly;
use ark_poly::UVPolynomial;

pub(crate) fn constant_poly<F: PrimeField>(c: F) -> Poly<F> {
    Poly::from_coefficients_vec(vec![c])
}

fn interpolate<F: PrimeField>(xs: &[F], ys: &[F]) -> Poly<F> {
    let x1 = xs[0];
    let mut l = crate::utils::z_of_point(&x1);
    for &xj in xs.iter().skip(1) {
        let q = crate::utils::z_of_point(&xj);
        l = &l * &q;
    }

    let mut ws = vec![];
    for xj in xs {
        let mut wj = F::one();
        for xk in xs {
            if xk != xj {
                let d = *xj - xk;
                wj *= d;
            }
        }
        ws.push(wj);
    }
    ark_ff::batch_inversion(&mut ws);

    let mut res = Poly::zero();
    for ((&wi, &xi), &yi) in ws.iter().zip(xs).zip(ys) {
        let d = crate::utils::z_of_point(&xi);
        let mut z = &l / &d;
        z = &z * wi;
        z = &z * yi;
        res = res + z;
    }
    res
}


/// returns `r(zeta)` and `z(zeta)`, where `r` is the interpolting poly, and `z` is the vanishing poly. 
pub(crate) fn interpolate_evaluate<F: PrimeField>(xs: &[F], ys: &[F], zeta: &F) -> (F, F) {
    assert_eq!(xs.len(), ys.len());

    let zeta_minus_xs = ark_std::iter::repeat(zeta).zip(xs.iter())
        .map(|(&zeta, xi)| zeta - xi)
        .collect::<Vec<_>>();

    let l_at_zeta = zeta_minus_xs.iter().cloned()
        .reduce(|acc, item| item * acc)
        .expect("TODO");

    let mut ws = vec![];
    for xj in xs {
        let mut wj = F::one();
        for xk in xs {
            if xk != xj {
                let d = *xj - xk;
                wj *= d;
            }
        }
        ws.push(wj);
    }

    let mut denominator = ws.into_iter().zip(zeta_minus_xs.iter())
        .map(|(a, b)| a * b)
        .collect::<Vec<_>>();

    ark_ff::batch_inversion(&mut denominator);

    let sum = denominator.into_iter().zip(ys.iter()).map(|(a, b)| a * b).sum::<F>();
    (sum * l_at_zeta, l_at_zeta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;
    use crate::tests::BenchField;
    use ark_ff::UniformRand;
    use ark_poly::Polynomial;
    use crate::utils::z_of_set;

    #[test]
    fn test_interpolation() {
        let rng = &mut test_rng();

        let d = 15;
        let (xs, ys): (Vec<_>, Vec<_>) = (0..d + 1)
            .map(|_| (BenchField::rand(rng), BenchField::rand(rng)))
            .unzip();

        let poly = interpolate(&xs, &ys);

        assert_eq!(poly.degree(), d);
        assert!(xs.iter().zip(ys.iter()).all(|(x, &y)| poly.evaluate(x) == y));

        for _ in 0..10 {
            let zeta = BenchField::rand(rng);
            let (r_at_zeta, z_at_zeta) = interpolate_evaluate(&xs, &ys, &zeta);
            assert_eq!(r_at_zeta, poly.evaluate(&zeta));
            assert_eq!(z_at_zeta, z_of_set(&xs).evaluate(&zeta));
        }
    }
}