use ark_ff::PrimeField;
use ark_std::ops::{Add, Sub};
use ark_std::iter::Sum;

pub trait AdditiveCommitment<F: PrimeField>:
Sized
+ Clone
+ Add<Self, Output=Self>
+ Sub<Self, Output=Self>
+ Sum<Self>
{
    fn mul(&self, by: F) -> Self;
}

/// Polynomial commitment scheme.
pub trait CommitmentScheme<F: PrimeField, P> {
    type G: AdditiveCommitment<F>;
    fn commit(&self, poly: &P) -> Self::G;
    fn commit_const(&self, c: &F) -> Self::G;
    /// Verifies that `p(x) = y` given the commitment to `p` and a proof.
    fn verify(&self, commitment: &Self::G, x: &F, y: F, proof: &Self::G) -> bool; // TODO: generalize proof type
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use ark_ff::Field;
    use ark_poly::UVPolynomial;
    use ark_std::marker::PhantomData;

    #[derive(Clone)]
    pub struct WrappedPolynomial<F: Field, P: UVPolynomial<F>>(pub P, PhantomData<F>);

    impl<F: Field, P: UVPolynomial<F>> WrappedPolynomial<F, P> {
        fn evaluate(&self, x: &F) -> F {
            self.0.evaluate(x)
        }
    }

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

        fn verify(&self, commitment: &Self::G, x: &F, y: F, _proof: &Self::G) -> bool {
            commitment.evaluate(x) == y
        }
    }
}