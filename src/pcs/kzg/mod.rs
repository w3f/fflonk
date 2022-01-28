mod urs;
mod params;
mod commitment;

use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_std::marker::PhantomData;
use crate::pcs::{PCS, CommitterKey};
use crate::pcs::kzg::params::{KzgCommitterKey, KzgVerifierKey};
use crate::Poly;
use crate::pcs::kzg::commitment::KzgCommitment;
use crate::pcs::kzg::urs::URS;
use ark_poly::{Polynomial, UVPolynomial};
use ark_ec::msm::VariableBase;
use ark_ff::{PrimeField, One};
use ark_ec::AffineCurve;

use ark_std::{end_timer, start_timer};
use ark_std::rand::Rng;

pub struct KZG<E: PairingEngine> {
    _engine: PhantomData<E>,
}


impl<E: PairingEngine> KZG<E> {
    fn z(x: E::Fr) -> Poly<E::Fr> {
        Poly::from_coefficients_slice(&[-x, E::Fr::one()])
    }

    fn q(p: &Poly<E::Fr>, d: &Poly<E::Fr>) -> Poly<E::Fr> {
        p / d
    }

    fn compute_quotient(p: &Poly<E::Fr>, x: E::Fr) -> Poly<E::Fr> {
        Self::q(p, &Self::z(x))
    }

    fn opening(g1: &E::G1Affine, c: &E::G1Projective, x: E::Fr, z: E::Fr, proof: E::G1Affine) -> (E::G1Projective, E::G1Affine) {
        (g1.mul(z) - c - proof.mul(x), proof)
    }
}

impl<E: PairingEngine> PCS<E::Fr> for KZG<E> {
    type G = KzgCommitment<E>;
    type Params = URS<E>;
    type Proof = E::G1Affine;

    fn setup<R: Rng>(max_degree: usize, rng: &mut R) -> Self::Params {
        URS::<E>::generate(max_degree + 1, 2, rng)
    }


    fn commit(ck: &KzgCommitterKey<E::G1Affine>, p: &Poly<E::Fr>) -> Self::G {
        assert!(p.degree() <= ck.max_degree(), "Can't commit to degree {} polynomial using {} bases", p.degree(), ck.powers_in_g1.len());

        let coeffs = p.coeffs.iter().map(|c| c.into_repr()).collect::<Vec<_>>();
        let commitment = VariableBase::msm(
            &ck.powers_in_g1,
            &coeffs,
        );

        KzgCommitment(commitment)
    }

    fn open(ck: &KzgCommitterKey<E::G1Affine>, p: &Poly<E::Fr>, x: E::Fr) -> Self::Proof {
        let q = Self::compute_quotient(p, x);
        Self::commit(ck, &q).0.into_affine()
    }

    fn verify(pvk: &KzgVerifierKey<E>, c: &Self::G, x: E::Fr, z: E::Fr, proof: Self::Proof) -> bool {
        let (agg, proof) = Self::opening(&pvk.g1, &c.0, x, z, proof);
        E::product_of_pairings(&[
            (agg.into_affine().into(), pvk.g2.clone()),
            (proof.into(), pvk.tau_in_g2.clone()),
        ]).is_one()
    }

    fn commit_to_one(pvk: &KzgVerifierKey<E>) -> Self::G {
        KzgCommitment(pvk.g1.into_projective())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bw6_761::{BW6_761, G1Projective, Fr, G1Affine};
    use ark_std::test_rng;
    use crate::pcs::PcsParams;
    use ark_poly::UVPolynomial;
    use ark_ff::UniformRand;
    use ark_ec::ProjectiveCurve;
    use ark_std::rand::prelude::SliceRandom;
    use ark_ff::BigInteger;
    use ark_std::cmp::Ordering::Greater;

    fn _test_minimal_kzg<E: PairingEngine>(log_n: usize) {
        let rng = &mut test_rng();

        let t_setup = start_timer!(|| format!("KZG setup of size 2^{} on {}", log_n, std::any::type_name::<E>()));
        let max_degree = 1 << log_n - 1;
        let urs = KZG::<E>::setup(max_degree, rng);
        end_timer!(t_setup);

        let ck = urs.ck();
        let vk = urs.rk();

        let p = Poly::<E::Fr>::rand(ck.max_degree(), rng);
        let x = E::Fr::rand(rng);
        let z = p.evaluate(&x);

        let t_commit = start_timer!(|| format!("Committing to a dense degree-{} polynomial", ck.max_degree()));
        let c = KZG::<E>::commit(&ck, &p);
        end_timer!(t_commit);

        let t_prove = start_timer!(|| "Generating an opening proof for a single point");
        let proof = KZG::<E>::open(&ck, &p, x);
        end_timer!(t_prove);

        let t_verify = start_timer!(|| "Verification of a single-point opening");
        assert!(KZG::<E>::verify(&vk, &c, x, z, proof));
        end_timer!(t_verify);
    }

    #[test]
    fn test_minimal_kzg() {
        _test_minimal_kzg::<BW6_761>(10);
    }

    #[test]
    fn bench_minimal_kzg() {
        _test_minimal_kzg::<BW6_761>(16);
    }
}