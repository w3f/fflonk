use ark_ec::{PairingEngine, AffineCurve};
use crate::pcs::{PcsParams, CommitterKey, VerifierKey};
use ark_std::rand::RngCore;
use ark_ff::{UniformRand, PrimeField};
use crate::utils;
use crate::pcs::kzg::urs::URS;

use ark_serialize::*;
use ark_std::io::{Read, Write};



impl<E: PairingEngine> URS<E> {
    /// Non-prepared verifier key. Can be used for serialization.
    pub fn raw_rk(&self) -> RawKzgVerifierKey<E> {
        assert!(self.powers_in_g1.len() > 0, "no G1 generator");
        assert!(self.powers_in_g2.len() > 1, "{} powers in G2", self.powers_in_g2.len());

        RawKzgVerifierKey {
            g1: self.powers_in_g1[0],
            g2: self.powers_in_g2[0],
            tau_in_g2: self.powers_in_g2[1],
        }
    }
}

impl<E: PairingEngine> PcsParams for URS<E> {
    type CommitterKey = KzgCommitterKey<E::G1Affine>;
    type VerifierKey = KzgVerifierKey<E>;

    fn ck(&self) -> Self::CommitterKey {
        KzgCommitterKey {
            powers_in_g1: self.powers_in_g1.clone() //TODO: Cow?
        }
    }

    fn rk(&self) -> Self::VerifierKey {
        self.raw_rk().prepare()
    }
}




/// Used to commit to and to open univariate polynomials of degree up to self.max_degree().
#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct KzgCommitterKey<G: AffineCurve> {
    // G1, tau.G1, tau^2.G1, ..., tau^n1.G1
    pub(crate) powers_in_g1: Vec<G>,
}

impl<G: AffineCurve> CommitterKey for KzgCommitterKey<G> {
    fn max_degree(&self) -> usize {
        self.powers_in_g1.len() - 1
    }
}



/// Verifier key with G2 elements not "prepared". Exists only to be serializable.
/// KzgVerifierKey is used for verification.
#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct RawKzgVerifierKey<E: PairingEngine> {
    g1: E::G1Affine, // generator of G1
    g2: E::G2Affine, // generator of G2
    tau_in_g2: E::G2Affine, // tau.g2
}

impl<E: PairingEngine> RawKzgVerifierKey<E> {
    /// Returns the key that is used to verify openings in a single point. It has points in G2 "prepared".
    /// "Preparation" is a pre-computation that makes pairing computation with these points more efficient.
    /// At the same time usual arithmetic operations are not implemented for "prepared" points.
    pub fn prepare(&self) -> KzgVerifierKey<E> {
        KzgVerifierKey {
            g1: self.g1,
            g2: self.g2.into(),
            tau_in_g2: self.tau_in_g2.into(),
        }
    }
}



/// "Prepared" verifier key capable of verifying opening in a single point, given the commitment is in G1.
/// Use RawKzgVerifierKey for serialization.
#[derive(Clone, Debug)]
pub struct KzgVerifierKey<E: PairingEngine> {
    // generator of G1
    pub(crate) g1: E::G1Affine, // G1Prepared is just a wrapper around G1Affine // TODO: fixed-base precomputations
    // generator of G2, prepared
    pub(crate) g2: E::G2Prepared, // G2Prepared can be used as a pairing RHS only
    // tau.g2, prepared
    pub(crate) tau_in_g2: E::G2Prepared, // G2Prepared can be used as a pairing RHS only
}

impl<E: PairingEngine> VerifierKey for KzgVerifierKey<E> {
    fn max_points(&self) -> usize {
        1
    }
}
