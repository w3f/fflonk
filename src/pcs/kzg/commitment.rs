use ark_ec::pairing::Pairing;
use crate::pcs::Commitment;
use ark_std::ops::{Add, Mul, Sub};
use ark_std::iter::Sum;

use ark_serialize::*;
use ark_std::io::{Read, Write};
use crate::utils::ec::{small_multiexp_affine};


/// KZG commitment to G1 represented in affine coordinates.
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct KzgCommitment<E: Pairing>(pub E::G1Affine);


impl <E: Pairing> Commitment<E::ScalarField> for KzgCommitment<E> {
    fn mul(&self, by: E::ScalarField) -> KzgCommitment<E> {
        KzgCommitment(self.0.mul(by).into())
    }

    fn combine(coeffs: &[<E as Pairing>::Fr], commitments: &[Self]) -> Self {
        let bases = commitments.iter().map(|c| c.0).collect::<Vec<_>>();
        let prod = small_multiexp_affine(coeffs, &bases);
        KzgCommitment(prod.into())
    }
}

impl<E: Pairing> Add<Self> for KzgCommitment<E> {
    type Output = KzgCommitment<E>;

    fn add(self, other: KzgCommitment<E>) -> KzgCommitment<E> {
        KzgCommitment(self.0 + other.0)
    }
}

impl<E: Pairing> Sub<Self> for KzgCommitment<E> {
    type Output = KzgCommitment<E>;

    fn sub(self, other: KzgCommitment<E>) -> KzgCommitment<E> {
        KzgCommitment(self.0 + -other.0)
    }
}

impl<E: Pairing> Sum<Self> for KzgCommitment<E> {
    fn sum<I: Iterator<Item=Self>>(iter: I) -> KzgCommitment<E> {
        KzgCommitment(iter.map(|c| c.0).sum())
    }
}