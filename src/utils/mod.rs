pub mod ec;
pub(crate) mod poly;

use ark_ff::{FftField, Field};
use ark_poly::{Polynomial, UVPolynomial};
use ark_poly::univariate::DensePolynomial;


pub fn powers<F: Field>(base: F) -> impl Iterator<Item=F> {
    ark_std::iter::successors(Some(F::one()), move |power| Some(base * power))
}


pub fn randomize<P, F>(
    r: F,
    polys: &[P],
) -> P
    where
        F: Field,
        P: Polynomial<F> {

    let mut res = P::zero();
    polys.iter().zip(powers(r)).for_each(|(p, power_of_r)| {
        res += (power_of_r, p);
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
pub fn curve_name<E: ark_ec::PairingEngine>() -> &'static str {
    // ark_ec::models::bw6::BW6<ark_bw6_761::curves::Parameters>
    let full_name = std::any::type_name::<E>();
    full_name.split_once("<").unwrap().1
        .split_once(":").unwrap().0
}