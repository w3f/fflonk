use ark_ff::{FftField, Field, PrimeField, Zero, BigInteger};
use ark_poly::{Polynomial, UVPolynomial};
use ark_poly::univariate::DensePolynomial;
use ark_ec::{AffineCurve, ProjectiveCurve};

/// (max_exp+1)-sized vec: 1, base, base^2,... ,base^{max_exp}
pub fn powers<F: Field>(base: F, max_exp: usize) -> Vec<F> {
    let mut result = Vec::with_capacity(max_exp + 1);
    result.push(F::one());
    if max_exp > 0 {
        result.push(base);
    }
    let mut curr = base;
    for _ in 1..max_exp {
        curr *= base;
        result.push(curr);
    };
    result
}

pub fn randomize<P, F>(
    r: F,
    polys: &[P],
) -> P
    where
        F: Field,
        P: Polynomial<F> {
    let mut res = P::zero();
    if polys.is_empty() {
        return res;
    }
    let powers = powers(r, polys.len() - 1);

    powers.into_iter().zip(polys).for_each(|(r, p)| {
        res += (r, p);
    });
    res
}

// Z(X) = X - x
pub fn z_of_point<F: Field>(x: &F) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(vec![x.neg(), F::one()])
}


pub fn z_of_set<'a, F: FftField>(xs: impl IntoIterator<Item=&'a F>) -> DensePolynomial<F> {
    xs.into_iter()
        .map(|xi| z_of_point(xi))
        .reduce(|p, pi| &p * &pi)
        .unwrap()
}


/// Performs a small multi-exponentiation operation.
/// Uses the double-and-add algorithm with doublings shared across points.
// adopted from https://github.com/zcash/halo2/pull/20
pub fn small_multiexp<G: AffineCurve>(coeffs: &[G::ScalarField], bases: &[G]) -> G::Projective {
    let bytes_in_repr = <G::ScalarField as PrimeField>::BigInt::NUM_LIMBS * 8;
    let coeffs: Vec<_> = coeffs.iter().map(|c| c.into_repr().to_bytes_le()).collect();

    let mut acc = G::Projective::zero();

    // for byte idx
    for byte_idx in (0..bytes_in_repr).rev() {
        // for bit idx
        for bit_idx in (0..8).rev() {
            acc.double_in_place();
            // for each coeff
            for coeff_idx in 0..coeffs.len() {
                let byte = coeffs[coeff_idx][byte_idx];
                if ((byte >> bit_idx) & 1) != 0 {
                    acc.add_assign_mixed(&bases[coeff_idx]);
                }
            }
        }
    }

    acc
}


pub fn naive_multiexp<G: AffineCurve>(coeffs: &[G::ScalarField], bases: &[G]) -> G::Projective {
    bases.iter().zip(coeffs.iter()).map(|(b, &c)| b.mul(c)).sum()
}


#[cfg(test)]
pub mod tests {
    use super::*;
    use ark_std::test_rng;
    use ark_ff::UniformRand;
    use ark_ec::PairingEngine;
    use ark_bw6_761::{Fr, G1Affine};


    pub fn curve_name<E: PairingEngine>() -> &'static str {
        // ark_ec::models::bw6::BW6<ark_bw6_761::curves::Parameters>
        let full_name = std::any::type_name::<E>();
        full_name.split_once("<").unwrap().1
            .split_once(":").unwrap().0
    }


    #[test]
    fn test_multiexp() {
        let rng = &mut test_rng();

        let n = 10;
        let exps = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let bases = (0..n).map(|_| G1Affine::rand(rng)).collect::<Vec<_>>();

        let naive = naive_multiexp(&exps, &bases);
        let small = small_multiexp(&exps, &bases);
        assert_eq!(naive, small)
    }
}