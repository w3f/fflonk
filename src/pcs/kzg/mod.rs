pub mod urs;
pub mod params;
mod commitment;

use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_std::ops::Mul;
use ark_std::marker::PhantomData;
use crate::pcs::{PCS, CommitterKey};
use crate::pcs::kzg::params::{KzgCommitterKey, KzgVerifierKey};
use crate::Poly;
use crate::pcs::kzg::commitment::KzgCommitment;
use crate::pcs::kzg::urs::URS;
use ark_poly::{Polynomial, DenseUVPolynomial};
use ark_ec::msm::VariableBaseMSM;
use ark_ff::{One, UniformRand};

use ark_std::rand::Rng;
use crate::utils::ec::{small_multiexp_proj, small_multiexp_affine};

pub struct KZG<E: PairingEngine> {
    _engine: PhantomData<E>,
}

/// e(acc, g2) = e(proof, tau.g2)
#[derive(Clone, Debug)]
pub struct AccumulatedOpening<E: PairingEngine> {
    pub acc: E::G1Affine,
    pub proof: E::G1Affine,
}

#[derive(Clone, Debug)]
pub struct KzgOpening<E: PairingEngine> {
    pub c: E::G1Affine,
    pub x: E::Fr,
    pub y: E::Fr,
    pub proof: E::G1Affine,
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

    fn parse(openings: Vec<KzgOpening<E>>) -> Vec<((E::G1Projective, E::G1Affine), E::Fr)> {
        openings.into_iter().map(|KzgOpening { c, x, y, proof }|
            ((proof.mul(x).add_mixed(&c), proof), y)
        ).collect()
    }

    pub fn accumulate(openings: Vec<KzgOpening<E>>, rs: &[E::Fr], vk: &KzgVerifierKey<E>) -> AccumulatedOpening<E> {
        assert_eq!(openings.len(), rs.len());
        let openings = Self::parse(openings);
        let ((accs, proofs), ys): ((Vec<E::G1Projective>, Vec<E::G1Affine>), Vec<E::Fr>) = openings.into_iter().unzip();
        let sum_ry = rs.iter().zip(ys.into_iter()).map(|(r, y)| y * r).sum::<E::Fr>();
        let acc = vk.g1.mul(sum_ry) - small_multiexp_proj(rs, &accs);
        let proof = small_multiexp_affine(rs, &proofs);
        E::G1Projective::batch_normalization(&mut [acc, proof]);
        let acc = acc.into_affine();
        let proof = proof.into_affine();
        AccumulatedOpening { acc, proof }
    }

    fn accumulate_single(opening: KzgOpening<E>, g1: &E::G1Affine) -> AccumulatedOpening<E> {
        let KzgOpening { c, x, y, proof } = opening;
        let acc = (g1.mul(y) - proof.mul(x).add_mixed(&c)).into_affine();
        AccumulatedOpening { acc, proof }
    }

    pub fn verify_accumulated(opening: AccumulatedOpening<E>, vk: &KzgVerifierKey<E>) -> bool {
        E::product_of_pairings(&[
            (opening.acc.into(), vk.g2.clone()),
            (opening.proof.into(), vk.tau_in_g2.clone()),
        ]).is_one()
    }

    pub fn verify_single(opening: KzgOpening<E>, vk: &KzgVerifierKey<E>) -> bool {
        let acc_opening = Self::accumulate_single(opening, &vk.g1);
        Self::verify_accumulated(acc_opening, vk)
    }

    pub fn verify_batch<R: Rng>(openings: Vec<KzgOpening<E>>, vk: &KzgVerifierKey<E>, rng: &mut R) -> bool {
        let one = std::iter::once(E::Fr::one());
        let coeffs: Vec<E::Fr> = one.chain((1..openings.len()).map(|_| u128::rand(rng).into())).collect();
        let acc_opening = Self::accumulate(openings, &coeffs, vk);
        Self::verify_accumulated(acc_opening, vk)
    }
}

impl<E: PairingEngine> PCS<E::Fr> for KZG<E> {
    type C = KzgCommitment<E>;
    type Proof = E::G1Affine;

    type CK = KzgCommitterKey<E::G1Affine>;
    type VK = KzgVerifierKey<E>;
    type Params = URS<E>;

    fn setup<R: Rng>(max_degree: usize, rng: &mut R) -> Self::Params {
        URS::<E>::generate(max_degree + 1, 2, rng)
    }

    fn commit(ck: &KzgCommitterKey<E::G1Affine>, p: &Poly<E::Fr>) -> Self::C {
        assert!(p.degree() <= ck.max_degree(), "Can't commit to degree {} polynomial using {} bases", p.degree(), ck.powers_in_g1.len());

        let commitment: E::G1Projective = VariableBaseMSM::msm(
            &ck.powers_in_g1,
            &p.coeffs,
        );

        KzgCommitment(commitment.into())
    }

    fn open(ck: &KzgCommitterKey<E::G1Affine>, p: &Poly<E::Fr>, x: E::Fr) -> Self::Proof {
        let q = Self::compute_quotient(p, x);
        Self::commit(ck, &q).0
    }

    fn verify(vk: &KzgVerifierKey<E>, c: Self::C, x: E::Fr, y: E::Fr, proof: Self::Proof) -> bool {
        let opening = KzgOpening { c: c.0, x, y, proof };
        Self::verify_single(opening, vk)
    }

    fn batch_verify<R: Rng>(vk: &KzgVerifierKey<E>, c: Vec<Self::C>, x: Vec<E::Fr>, y: Vec<E::Fr>, proof: Vec<Self::Proof>, rng: &mut R) -> bool {
        assert_eq!(c.len(), x.len());
        assert_eq!(c.len(), y.len());
        let openings = c.into_iter().zip(x.into_iter()).zip(y.into_iter()).zip(proof.into_iter())
            .map(|(((c, x), y), proof)| KzgOpening {c: c.0, x, y, proof})
            .collect();
        Self::verify_batch(openings, vk, rng)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;
    use crate::pcs::PcsParams;
    use ark_poly::DenseUVPolynomial;
    use ark_ff::UniformRand;

    use ark_std::{end_timer, start_timer};
    use crate::tests::{BenchCurve, TestCurve};

    fn _test_minimal_kzg<E: PairingEngine>(log_n: usize) {
        let rng = &mut test_rng();

        let max_degree = (1 << log_n) - 1;

        let t_setup = start_timer!(|| format!("KZG setup of size 2^{} on {}", log_n, crate::utils::curve_name::<E>()));
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

    fn random_openings<R: Rng, E: PairingEngine>(
        k: usize,
        ck: &KzgCommitterKey<E::G1Affine>,
        xs: Vec<E::Fr>,
        rng: &mut R,
    ) -> Vec<KzgOpening<E>>
    {
        assert_eq!(xs.len(), k);
        let d = ck.max_degree();

        (0..k).map(|i| {
            let f = Poly::<E::Fr>::rand(d, rng);
            let x = xs[i];
            let y = f.evaluate(&x);
            let c = KZG::<E>::commit(ck, &f).0;
            let proof = KZG::<E>::open(ck, &f, x);
            KzgOpening { c, x, y, proof }
        }).collect()
    }

    fn _test_batch_verification<E: PairingEngine>(log_n: usize, k: usize) {
        let rng = &mut test_rng();

        let max_degree = (1 << log_n) - 1;

        let urs = KZG::<E>::setup(max_degree, rng);
        let (ck, vk) = (urs.ck(), urs.vk());

        let xs = (0..k).map(|_| E::Fr::rand(rng)).collect();
        let openings = random_openings(k, &ck, xs, rng);
        let t_verify_batch = start_timer!(|| format!("Batch verification of {} openings of degree ~2^{} on {} with {}-bit xs", k, log_n, crate::utils::curve_name::<E>(), E::Fr::MODULUS_BIT_SIZE));
        assert!(KZG::<E>::verify_batch(openings, &vk, rng));
        end_timer!(t_verify_batch);

        let xs = (0..k).map(|_| E::Fr::from(u128::rand(rng))).collect();
        let openings = random_openings(k, &ck, xs, rng);
        let t_verify_batch = start_timer!(|| format!("Batch verification of {} openings of degree ~2^{} on {} with {}-bit xs", k, log_n, crate::utils::curve_name::<E>(), 128));
        assert!(KZG::<E>::verify_batch(openings, &vk, rng));
        end_timer!(t_verify_batch);
    }

    #[test]
    fn test_minimal_kzg() {
        _test_minimal_kzg::<TestCurve>(8);
    }

    #[test]
    #[ignore]
    fn bench_minimal_kzg() {
        _test_minimal_kzg::<BenchCurve>(16);
    }

    #[test]
    fn test_batch_verification() {
        _test_batch_verification::<TestCurve>(8, 4);
    }

    #[test]
    #[ignore]
    fn bench_batch_verification() {
        _test_batch_verification::<BenchCurve>(12, 5);
    }
}