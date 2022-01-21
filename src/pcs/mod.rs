use ark_ff::{PrimeField, Field};
use ark_std::ops::{Add, Sub};
use ark_std::iter::Sum;
use ark_poly::{UVPolynomial, Polynomial};
use ark_poly::univariate::{DensePolynomial, DenseOrSparsePolynomial};
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

// pub trait EuclideanPolynomial<F: PrimeField>: UVPolynomial<F> {
//     fn divide_with_q_and_r(&self, divisor: &DensePolynomial<F>) -> (Self, Self);
// }
//
// impl<F: PrimeField> EuclideanPolynomial<F> for DensePolynomial<F> {
//     fn divide_with_q_and_r(&self, divisor: &DensePolynomial<F>) -> (Self, Self) {
//         let a: DenseOrSparsePolynomial<F> = self.into();
//         let b: DenseOrSparsePolynomial<F> = divisor.into();
//         a.divide_with_q_and_r(&b).unwrap()
//     }
// }

pub trait VerifierKey<F, G> {
    fn commit_to_one(&self) -> G;
}

/// Polynomial commitment scheme.
pub trait PCS<F: PrimeField> {
    type G: CommitmentSpace<F>;
    type CK;
    type OK;
    type VK: VerifierKey<F, Self::G>;
    type Proof;

    fn setup() -> (Self::CK, Self::OK, Self::VK);
    fn commit(ck: &Self::CK, p: &Poly<F>) -> Self::G;
    fn open(ok: &Self::OK, p: &Poly<F>, x: &F) -> Self::Proof; //TODO: eval?
    fn verify(vk: &Self::VK, c: &Self::G, x: &F, z: &F, proof: &Self::Proof) -> bool;
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use ark_ff::Field;
    use ark_poly::UVPolynomial;
    use ark_std::marker::PhantomData;
    use ark_ff::Zero;
    use ark_poly::univariate::DensePolynomial;
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


    impl<F: PrimeField> VerifierKey<F, WrappedPolynomial<F>> for () {
        fn commit_to_one(&self) -> WrappedPolynomial<F> {
            WrappedPolynomial(Poly::from_coefficients_slice(&[F::one()]))
        }
    }


    pub struct IdentityCommitment {}

    impl<F: PrimeField> PCS<F> for IdentityCommitment {
        type G = WrappedPolynomial<F>;
        type CK = ();
        type OK = ();
        type VK = ();
        type Proof = WrappedPolynomial<F>;

        fn setup() -> (Self::CK, Self::OK, Self::VK) {
            ((), (), ())
        }

        fn commit(ck: &Self::CK, p: &Poly<F>) -> Self::G {
            WrappedPolynomial(p.clone())
        }

        fn open(ok: &Self::OK, p: &Poly<F>, x: &F) -> Self::Proof {
            WrappedPolynomial(p.clone())
        }

        fn verify(vk: &Self::VK, c: &Self::G, x: &F, z: &F, proof: &Self::Proof) -> bool {
            c.evaluate(x) == *z
        }
    }
}