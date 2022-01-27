pub mod kzg;

use ark_ff::PrimeField;
use ark_std::iter::Sum;
use ark_std::ops::{Add, Sub};

use crate::Poly;
use ark_serialize::*;
use ark_std::io::{Read, Write};
use ark_std::fmt::Debug;
use ark_std::rand::Rng;


pub trait CommitmentSpace<F: PrimeField>:
Sized
+ Clone
+ Add<Self, Output=Self>
+ Sub<Self, Output=Self>
+ Sum<Self>
{
    fn mul(&self, by: F) -> Self;
}


/// Can be used to commit and open commitments to DensePolynomial<F> of degree up to max_degree.
pub trait CommitterKey: Clone + Debug + CanonicalSerialize + CanonicalDeserialize {
    /// Maximal degree of a polynomial supported.
    fn max_degree(&self) -> usize;

    /// Maximal number of evaluations supported when committing in the Lagrangian base.
    fn max_evals(&self) -> usize {
        self.max_degree() + 1
    }
}


/// Can be used to verify openings to commitments.
pub trait VerifierKey: Clone + Debug {
    /// Maximal number of openings of the same commitment that can be verified.
    fn max_points(&self) -> usize {
        1
    }
}


pub trait PcsParams {
    type CommitterKey: CommitterKey;
    type VerifierKey: VerifierKey;
    fn ck(&self) -> Self::CommitterKey; //TODO: trim
    fn rk(&self) -> Self::VerifierKey;
}


/// Polynomial commitment scheme.
pub trait PCS<F: PrimeField> {
    type G: CommitmentSpace<F>;
    type Params: PcsParams;
    type Proof;

    fn setup<R: Rng>(max_degree: usize, rng: &mut R) -> Self::Params;

    fn commit(ck: &<Self::Params as PcsParams>::CommitterKey, p: &Poly<F>) -> Self::G;

    fn open(ck: &<Self::Params as PcsParams>::CommitterKey, p: &Poly<F>, x: F) -> Self::Proof; //TODO: eval?

    fn verify(pvk: &<Self::Params as PcsParams>::VerifierKey, c: &Self::G, x: F, z: F, proof: Self::Proof) -> bool;

    /// Commit to degree-0 polynomials (see shplonk scheme #2 ot Halo infinite 4.2).
    fn commit_to_one(pvk: &<Self::Params as PcsParams>::VerifierKey) -> Self::G;
}


#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    use ark_ff::Zero;
    use ark_poly::UVPolynomial;
    use ark_poly::Polynomial;

    use crate::Poly;

    #[derive(Clone)]
    pub struct WrappedPolynomial<F: PrimeField>(pub Poly<F>);

    impl<F: PrimeField> WrappedPolynomial<F> {
        fn evaluate(&self, x: &F) -> F {
            self.0.evaluate(x)
        }
    }

    impl<F: PrimeField> Add<Self> for WrappedPolynomial<F> {
        type Output = WrappedPolynomial<F>;

        fn add(self, other: WrappedPolynomial<F>) -> Self::Output {
            WrappedPolynomial(self.0 + other.0)
        }
    }

    impl<F: PrimeField> Sub<Self> for WrappedPolynomial<F> {
        type Output = WrappedPolynomial<F>;

        fn sub(self, other: WrappedPolynomial<F>) -> Self::Output {
            let mut temp = self.0;
            temp -= &other.0; //TODO
            WrappedPolynomial(temp)
        }
    }

    impl<F: PrimeField> core::iter::Sum<Self> for WrappedPolynomial<F> {
        fn sum<I: Iterator<Item=Self>>(iter: I) -> Self {
            iter.reduce(|a, b| a + b).unwrap()
        }
    }

    impl<F: PrimeField> CommitmentSpace<F> for WrappedPolynomial<F> {
        fn mul(&self, by: F) -> Self {
            let mut temp = Poly::zero(); //TODO
            temp += (by, &self.0);
            WrappedPolynomial(temp)
        }
    }


    impl CommitterKey for () {
        fn max_degree(&self) -> usize {
            usize::MAX >> 1
        }
    }

    impl VerifierKey for () {
        fn max_points(&self) -> usize {
            1
        }
    }


    impl PcsParams for () {
        type CommitterKey = ();
        type VerifierKey = ();

        fn ck(&self) -> Self::CommitterKey {
            ()
        }

        fn rk(&self) -> Self::VerifierKey {
            ()
        }
    }


    pub struct IdentityCommitment {}

    impl<F: PrimeField> PCS<F> for IdentityCommitment {
        type G = WrappedPolynomial<F>;
        type Params = ();
        type Proof = ();

        fn setup<R: Rng>(max_degree: usize, rng: &mut R) -> Self::Params {
            ()
        }

        fn commit(ck: &(), p: &Poly<F>) -> Self::G {
            WrappedPolynomial(p.clone())
        }

        fn open(ck: &(), p: &Poly<F>, x: F) -> Self::Proof {
            ()
        }

        fn verify(pvk: &(), c: &Self::G, x: F, z: F, proof: Self::Proof) -> bool {
            c.evaluate(&x) == z
        }

        fn commit_to_one(pvk: &()) -> Self::G {
            WrappedPolynomial(Poly::from_coefficients_slice(&[F::one()]))
        }
    }
}