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
        let t_verify = start_timer!(|| "1-point KZG verification");

        let (agg, proof) = Self::opening(&pvk.g1, &c.0, x, z, proof);
        let pp = E::product_of_pairings(&[
            (agg.into_affine().into(), pvk.g2.clone()),
            (proof.into(), pvk.tau_in_g2.clone()),
        ]);
        let valid = pp.is_one();

        end_timer!(t_verify, || format!("valid = {}", valid));
        return valid;
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

    fn _test_commit<E: PairingEngine>() {
        let rng = &mut test_rng();

        let num_bases = 1 << 10;
        let urs = KZG::<E>::setup(num_bases - 1, rng);
        let ck = urs.ck();
        let vk = urs.rk();

        let p = Poly::<E::Fr>::rand(ck.max_degree(), rng);
        let x = E::Fr::rand(rng);
        let z = p.evaluate(&x);

        let c = KZG::<E>::commit(&ck, &p);
        let proof = KZG::<E>::open(&ck, &p, x);
        assert!(KZG::<E>::verify(&vk, &c, x, z, proof));
    }

    #[test]
    fn test_commit() {
        _test_commit::<BW6_761>();
    }

    // #[test]
    // fn bench_msm() {
    //     let n = 1 << 16;
    //     let rng = &mut test_rng();
    //
    //     type S = <Fr as PrimeField>::BigInt;
    //
    //     let bases = (0..2*n).map(|_| G1Projective::rand(rng).into_affine()).collect::<Vec<_>>();
    //     let scalars = (0..2*n).map(|_| S::rand(rng)).collect::<Vec<_>>();
    //
    //     let t_msm = start_timer!(|| "msm");
    //     VariableBase::msm_checked_len(&bases, &scalars);
    //     end_timer!(t_msm);
    //
    //     let zeroes = vec![S::from(0); 2*n];
    //     let t_msm = start_timer!(|| "msm by 0");
    //     VariableBase::msm_checked_len(&bases, &zeroes);
    //     end_timer!(t_msm);
    //
    //     let zeroes = vec![S::from(0); n];
    //     let scalars = (0..n).map(|_| S::rand(rng)).collect::<Vec<_>>();
    //     let mut half_zeroes: Vec<S> = [zeroes, scalars].concat();
    //
    //     let t_msm = start_timer!(|| "1/2 0");
    //     VariableBase::msm_checked_len(&bases, &half_zeroes);
    //     end_timer!(t_msm);
    //
    //     half_zeroes.shuffle(rng);
    //     let t_msm = start_timer!(|| "shuffled");
    //     VariableBase::msm_checked_len(&bases, &half_zeroes);
    //     end_timer!(t_msm);
    //
    //
    //     let t_filter = start_timer!(|| "filter");
    //     let mut pairs: Vec<(S, G1Affine)> = half_zeroes.into_iter().zip(bases.into_iter()).collect();
    //     pairs.retain(|(s, b)| !s.is_zero());
    //     assert_eq!(pairs.len(), n);
    //     let (scalars, bases): (Vec<S>, Vec<G1Affine>)  = pairs.into_iter().unzip();
    //     end_timer!(t_filter);
    //
    //     let t_msm = start_timer!(|| "filtered");
    //     VariableBase::msm_checked_len(&bases, &scalars);
    //     end_timer!(t_msm);
    // }
}