pub mod kzg;
pub mod aggregation;

use ark_ff::PrimeField;

use ark_std::iter::Sum;
use ark_std::ops::{Add, Sub};
use ark_std::fmt::Debug;
use ark_std::rand::Rng;
use ark_serialize::*;

use crate::Poly;


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
    /// Maximal number of openings that can be verified.
    fn max_points(&self) -> usize {
        1
    }
}


pub trait PcsParams<CK, VK> {
    fn ck(&self) -> CK; //TODO: trim
    fn vk(&self) -> VK;
}


/// Polynomial commitment scheme.
pub trait PCS<F: PrimeField> {
    type G: CommitmentSpace<F>;

    type Proof;

    type CK: CommitterKey;
    type VK: VerifierKey + Into<Self::CK>;
    type Params: PcsParams<Self::CK, Self::VK>;

    fn setup<R: Rng>(max_degree: usize, rng: &mut R) -> Self::Params;

    fn commit(ck: &Self::CK, p: &Poly<F>) -> Self::G;

    fn open(ck: &Self::CK, p: &Poly<F>, x: F) -> Self::Proof; //TODO: eval?

    fn verify(vk: &Self::VK, c: Self::G, x: F, z: F, proof: Self::Proof) -> bool;
}


#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    use ark_ff::Zero;
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


    impl PcsParams<(), ()> for () {
        fn ck(&self) -> () {
            ()
        }

        fn vk(&self) -> () {
            ()
        }
    }


    pub struct IdentityCommitment {}

    impl<F: PrimeField> PCS<F> for IdentityCommitment {
        type G = WrappedPolynomial<F>;
        type Params = ();
        type Proof = ();
        type CK = ();
        type VK = ();

        fn setup<R: Rng>(_max_degree: usize, _rng: &mut R) -> Self::Params {
            ()
        }

        fn commit(_ck: &(), p: &Poly<F>) -> Self::G {
            WrappedPolynomial(p.clone())
        }

        fn open(_ck: &(), _p: &Poly<F>, _x: F) -> Self::Proof {
            ()
        }

        fn verify(_vk: &(), c: Self::G, x: F, z: F, _proof: Self::Proof) -> bool {
            c.evaluate(&x) == z
        }
    }
}