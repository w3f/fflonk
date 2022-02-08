use ark_ff::{PrimeField, Zero};
use crate::Poly;

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

pub(crate) fn interpolate_evaluate<F: PrimeField>(xs: &[F], ys: &[F], zeta: &F) -> F {
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

    denominator.into_iter().zip(ys.iter()).map(|(a, b)| a * b).sum::<F>() * l_at_zeta
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;
    use crate::tests::TestField;
    use ark_ff::UniformRand;
    use ark_poly::Polynomial;

    #[test]
    fn test_interpolation() {
        let rng = &mut test_rng();

        let d = 15;
        let (xs, ys): (Vec<_>, Vec<_>) = (0..d + 1)
            .map(|_| (TestField::rand(rng), TestField::rand(rng)))
            .unzip();

        let poly = interpolate(&xs, &ys);

        assert_eq!(poly.degree(), d);
        assert!(xs.iter().zip(ys.iter()).all(|(x, &y)| poly.evaluate(x) == y));

        for _ in 0..10 {
            let zeta = TestField::rand(rng);
            let poly_at_zeta = interpolate_evaluate(&xs, &ys, &zeta);
            assert_eq!(poly_at_zeta, poly.evaluate(&zeta));
        }
    }
}