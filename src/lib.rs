extern crate ark_ff;
extern crate ark_poly;
extern crate ark_std;

use ark_ff::{Field, FftField, PrimeField};
use ark_poly::{Polynomial, UVPolynomial};
use ark_poly::univariate::DensePolynomial;
use ark_std::ops::{Sub, Add};

pub mod shplonk;
pub mod fflonk;

pub trait AdditiveCommitment<F: PrimeField>:
    Sized
    + Add<Self, Output=Self>
    + Sub<Self, Output=Self>
    + core::iter::Sum<Self>
{
    fn mul(&self, by: F) -> Self; //TODO: Self::Projective;
}

pub trait CommitmentScheme<F: PrimeField, P> {
    type G: AdditiveCommitment<F>;
    fn commit(&self, poly: &P) -> Self::G;
    fn commit_const(&self, c: &F) -> Self::G;
}

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


fn z_of_set<'a, F: FftField>(xs: impl IntoIterator<Item=&'a F>) -> DensePolynomial<F> {
    xs.into_iter()
        .map(|xi| z_of_point(xi))
        .reduce(|p, pi| &p * &pi)
        .unwrap()
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::marker::PhantomData;

    pub struct WrappedPolynomial<F: Field, P: UVPolynomial<F>>(pub P, PhantomData<F>);

    impl<F: PrimeField, P: UVPolynomial<F>> Add<Self> for WrappedPolynomial<F, P> {
        type Output = WrappedPolynomial<F, P>;

        fn add(self, other: WrappedPolynomial<F, P>) -> Self::Output {
            WrappedPolynomial(self.0 + other.0, PhantomData)
        }
    }

    impl<F: PrimeField, P: UVPolynomial<F>> Sub<Self> for WrappedPolynomial<F, P> {
        type Output = WrappedPolynomial<F, P>;

        fn sub(self, other: WrappedPolynomial<F, P>) -> Self::Output {
            let mut temp = self.0;
            temp -= &other.0; //TODO
            WrappedPolynomial(temp, PhantomData)
        }
    }

    impl<F: PrimeField, P: UVPolynomial<F>> core::iter::Sum<Self> for WrappedPolynomial<F, P> {
        fn sum<I: Iterator<Item=Self>>(iter: I) -> Self {
            iter.reduce(|a, b| a + b).unwrap()
        }
    }

    impl<F: PrimeField, P: UVPolynomial<F>> AdditiveCommitment<F> for WrappedPolynomial<F, P> {
        fn mul(&self, by: F) -> Self {
            let mut temp = P::zero(); //TODO
            temp += (by, &self.0);
            WrappedPolynomial(temp, PhantomData)
        }
    }

    pub struct IdentityCommitment {}

    impl<F: PrimeField, P: UVPolynomial<F>> CommitmentScheme<F, P> for IdentityCommitment {
        type G = WrappedPolynomial<F, P>;

        fn commit(&self, poly: &P) -> Self::G {
            WrappedPolynomial(poly.clone(), PhantomData)
        }

        fn commit_const(&self, c: &F) -> Self::G {
            WrappedPolynomial(P::from_coefficients_vec(vec![*c]), PhantomData)
        }
    }
}