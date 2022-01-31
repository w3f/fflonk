use ark_ff::{FftField, Field};
use ark_poly::{Polynomial, UVPolynomial};
use ark_poly::univariate::DensePolynomial;
use ark_ec::PairingEngine;

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

#[cfg(test)]
pub fn curve_name<E: PairingEngine>() -> &'static str {
    // ark_ec::models::bw6::BW6<ark_bw6_761::curves::Parameters>
    let full_name = std::any::type_name::<E>();
    full_name.split_once("<").unwrap().1
        .split_once(":").unwrap().0
}