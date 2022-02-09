use ark_ec::PairingEngine;
use crate::pcs::Commitment;
use ark_ff::PrimeField;
use ark_std::ops::Add;
use std::ops::Sub;
use ark_std::iter::Sum;
use ark_ec::ProjectiveCurve;

use ark_serialize::*;
use ark_std::io::{Read, Write};
use crate::utils::ec::small_multiexp_proj;


/// KZG commitment to G1 represented in projective coordinates.
// Is serialized in affine, so muight make sense to batch-convert if there are many points.
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct KzgCommitment<E: PairingEngine>(pub E::G1Projective);

impl <E: PairingEngine> Commitment<E::Fr> for KzgCommitment<E> {
    fn mul(&self, by: E::Fr) -> KzgCommitment<E> {
        let x: E::G1Projective = self.0.mul(by.into_repr());
        KzgCommitment(x)
    }

    fn combine(coeffs: &[<E as PairingEngine>::Fr], commitments: &[Self]) -> Self {
        let bases = commitments.iter().map(|c| c.0).collect::<Vec<_>>();
        let prod = small_multiexp_proj(coeffs, bases.as_slice());
        KzgCommitment(prod)
    }
}

impl<E: PairingEngine> Add<Self> for KzgCommitment<E> {
    type Output = KzgCommitment<E>;

    fn add(self, other: KzgCommitment<E>) -> KzgCommitment<E> {
        KzgCommitment(self.0 + other.0)
    }
}

impl<E: PairingEngine> Sub<Self> for KzgCommitment<E> {
    type Output = KzgCommitment<E>;

    fn sub(self, other: KzgCommitment<E>) -> KzgCommitment<E> {
        KzgCommitment(self.0 - other.0)
    }
}

impl<E: PairingEngine> Sum<Self> for KzgCommitment<E> {
    fn sum<I: Iterator<Item=Self>>(iter: I) -> KzgCommitment<E> {
        let s: E::G1Projective = iter.map(|c| c.0).sum();
        KzgCommitment(s)
    }
}