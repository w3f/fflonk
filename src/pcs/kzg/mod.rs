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

use ark_std::rand::Rng;

pub struct KZG<E: PairingEngine> {
    _engine: PhantomData<E>,
}


/// Represents a claim that f(z) = v, for some f such that CS::commit(f) = c.
/// Notion of "some f" depends on the soundness properties of the commitment scheme.
struct Claim<F: PrimeField, CS: PCS<F>> {
    c: CS::G,
    z: F,
    v: F,
}

/// Represents an opening claim together with an alleged proof
struct Opening<F: PrimeField, CS: PCS<F>> {
    claim: Claim<F, CS>,
    proof: CS::Proof,
}

struct PreparedOpening<F: PrimeField, CS: PCS<F>, E: PairingEngine> {
    acc: E::G1Projective,
    acc_proof: CS::Proof,
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

    fn prepare(opening: Opening<E::Fr, Self>, pvk: &<KZG<E> as PCS<E::Fr>>::VK) -> PreparedOpening<E::Fr, Self, E> {
        let Opening { claim, proof } = opening;
        let Claim { c, z, v } = claim;
        let acc = pvk.g1.mul(v) - c.0 - proof.mul(z);
        PreparedOpening {
            acc,
            acc_proof: proof
        }
    }

    fn verify_prepared(opening: PreparedOpening<E::Fr, Self, E>, pvk: &<KZG<E> as PCS<E::Fr>>::VK) -> bool {
        E::product_of_pairings(&[
            (opening.acc.into_affine().into(), pvk.g2.clone()),
            (opening.acc_proof.into(), pvk.tau_in_g2.clone()),
        ]).is_one()
    }
}

impl<E: PairingEngine> PCS<E::Fr> for KZG<E> {
    type G = KzgCommitment<E>;
    type Proof = E::G1Affine;

    type VK = KzgVerifierKey<E>;
    type CK = KzgCommitterKey<E::G1Affine>;
    type Params = URS<E>;


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

    fn verify(vk: &KzgVerifierKey<E>, c: Self::G, z: E::Fr, v: E::Fr, proof: Self::Proof) -> bool {
        let claim = Claim{ c, z, v };
        let opening = Opening {claim, proof};
        let prepared = Self::prepare(opening, vk);
        Self::verify_prepared(prepared, vk)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bw6_761::BW6_761;
    use ark_std::test_rng;
    use crate::pcs::PcsParams;
    use crate::utils::tests::curve_name;
    use ark_poly::UVPolynomial;
    use ark_ff::UniformRand;

    use ark_std::{end_timer, start_timer};

    fn _test_minimal_kzg<E: PairingEngine>(log_n: usize) {
        let rng = &mut test_rng();

        let t_setup = start_timer!(|| format!("KZG setup of size 2^{} on {}", log_n, curve_name::<E>()));
        let max_degree = (1 << log_n) - 1;
        let urs = KZG::<E>::setup(max_degree, rng);
        end_timer!(t_setup);

        let ck = urs.ck();
        let vk = urs.vk();

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
        assert!(KZG::<E>::verify(&vk, c, x, z, proof));
        end_timer!(t_verify);
    }

    #[test]
    fn test_minimal_kzg() {
        _test_minimal_kzg::<BW6_761>(8);
    }

    #[test]
    #[ignore]
    fn bench_minimal_kzg() {
        _test_minimal_kzg::<BW6_761>(16);
    }
}