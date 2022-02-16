use std::marker::PhantomData;
use ark_ff::{PrimeField, UniformRand, One};
use fflonk::pcs::{PCS, PcsParams, Commitment};
use crate::{DecoyPlonk, VanillaPlonkAssignments, random_polynomials};
use fflonk::Poly;
use ark_poly::{Radix2EvaluationDomain, Polynomial, EvaluationDomain, UVPolynomial};
use fflonk::shplonk::AggregateProof;
use ark_std::rand::Rng;
use ark_std::{end_timer, start_timer, test_rng};
use fflonk::utils::poly;
use fflonk::aggregation::single::aggregate_claims_multiexp;
use ark_bls12_381::{Fr, Bls12_381};
use fflonk::pcs::kzg::{KzgOpening, KZG};
use ark_ec::PairingEngine;

impl<F: PrimeField> VanillaPlonkAssignments<F> {
    fn constraints(&self) -> Vec<Poly<F>> {
        vec![
            self.arithmetic_constraint.clone(),
            self.permutation_constraint_1.clone(),
            self.permutation_constraint_2.clone(),
        ]
    }

    fn aggregated_constraints_quotient(&self, alpha: F) -> Poly<F> {
        let aggregate_constraint = poly::sum_with_powers(alpha, &self.constraints());
        self.quotient(&aggregate_constraint)
    }

    fn permutation_polynomial_at_zeta_omega(&self, zeta: F) -> F {
        let zeta_omega = zeta * self.omega;
        self.permutation_polynomial.evaluate(&zeta_omega)
    }
}

struct PlonkBatchKzgTest<F: PrimeField, CS: PCS<F>> {
    polys: VanillaPlonkAssignments<F>,

    linearization_polynomial: Poly<F>,
    // verifier challenges in order:
    // (beta, gamma): (F, F) // permutation argument challenges (aka "permutation challenges")
    alpha: F,
    // constraint aggregation challenge (aka "quotient challenge")
    zeta: F,
    // e evaluation challenge
    nus: Vec<F>,
    // polynomial aggregation challenge (aka "opening challenge")
    cs: PhantomData<CS>,
}

impl<F: PrimeField, CS: PCS<F>> PlonkBatchKzgTest<F, CS> {
    fn new(log_n: usize) -> Self {
        let rng = &mut test_rng();

        let n = 1 << log_n;

        let polys = VanillaPlonkAssignments::<F>::new(n, rng);
        let linearization_polynomial = Poly::rand(polys.max_degree, rng);

        let alpha: F = Self::get_128_bit_challenge(rng);
        let zeta: F = Self::get_128_bit_challenge(rng);
        let nu = std::iter::once(F::one());
        let nus = nu.chain((1..6).map(|_| Self::get_128_bit_challenge(rng))).collect();

        Self {
            polys,
            linearization_polynomial,
            alpha,
            zeta,
            nus,
            cs: PhantomData,
        }
    }

    fn get_128_bit_challenge<R: Rng>(rng: &mut R) -> F {
        u128::rand(rng).into()
    }

    fn commit_polynomials(&self, ck: &CS::CK, polys: &Vec<Poly<F>>) -> Vec<CS::C> {
        let t_commitment = start_timer!(|| format!("Committing to {} polynomials", polys.len()));
        let commitments = polys.iter().map(|p| CS::commit(ck, &p)).collect();
        end_timer!(t_commitment);

        commitments
    }

    //TODO
    fn polynomials_to_open_at_zeta(&self) -> Vec<Poly<F>> {
        let mut res = vec![self.linearization_polynomial.clone()];
        res.extend_from_slice(&self.polys.wire_polynomials);
        let preprocessed_permutation_polynomials_1_and_2 = &self.polys.preprocessed_polynomials[5..7]; // S_{sigma_1}, S_{sigma_2}
        res.extend_from_slice(preprocessed_permutation_polynomials_1_and_2);
        assert_eq!(res.len(), 6);
        assert_eq!(res.iter().map(|p| p.degree()).max().unwrap(), self.polys.max_degree);
        res
    }

    fn aggregated_polynomial_to_open_in_zeta(&self, nus: Vec<F>) -> Poly<F> {
        poly::sum_with_coeffs(nus, &self.polynomials_to_open_at_zeta())
    }

    fn evaluations_at_zeta(&self, zeta: F) -> Vec<F> {
        let polys = self.polynomials_to_open_at_zeta();
        assert_eq!(polys.len(), 6);
        polys.iter().map(|p| p.evaluate(&zeta)).collect()
    }
}

struct BatchyPlonkProof<F, C> {
    aggregate_kzg_proof_at_zeta: C,
    // [W_{\zeta}]_1
    permutation_polynomial_kzg_proof_at_zeta_omega: C,
    // [W_{\zeta\omega}]_1
    evals_at_zeta: Vec<F>,
    permutation_polynomial_at_zeta_omega: F,
}

impl<F: PrimeField, CS: PCS<F>> DecoyPlonk<F, CS> for PlonkBatchKzgTest<F, CS> {
    type Proof = BatchyPlonkProof<F, CS::Proof>;

    fn setup<R: Rng>(&mut self, rng: &mut R) -> (<CS as PCS<F>>::CK, <CS as PCS<F>>::VK) {
        let urs_degree = 123; //TODO

        let t_setup = start_timer!(|| format!("KZG setup of degree {} on {}",
                urs_degree, fflonk::utils::curve_name::<TestCurve>()));
        let params = CS::setup(urs_degree, rng);
        end_timer!(t_setup);
        (params.ck(), params.vk())
    }

    fn preprocess(&mut self, ck: &CS::CK) -> Vec<CS::C> {
        let t_preprocessing = start_timer!(|| "Preprocessing");
        let commitments = self.commit_polynomials(ck, &self.polys.preprocessed_polynomials);
        end_timer!(t_preprocessing);
        commitments
    }

    fn prove(&mut self, ck: &CS::CK) -> (Self::Proof, Vec<CS::C>) {
        let t_proving = start_timer!(|| "Proving");
        let mut proof_commitments = vec![];
        proof_commitments.extend(self.commit_polynomials(ck, &self.polys.wire_polynomials));
        proof_commitments.extend(self.commit_polynomials(ck, &vec![self.polys.permutation_polynomial.clone()]));
        let q = self.polys.aggregated_constraints_quotient(self.alpha);
        assert_eq!(q.degree() as u64, 3 * self.polys.domain.size - 4);
        proof_commitments.extend(self.commit_polynomials(ck, &vec![q]));

        let aggregate_kzg_proof_at_zeta = CS::open(ck, &self.aggregated_polynomial_to_open_in_zeta(self.nus.clone()), self.zeta);
        let permutation_polynomial_kzg_proof_at_zeta_omega = CS::open(ck, &self.polys.permutation_polynomial, self.zeta * self.polys.omega);
        end_timer!(t_proving);

        let evals_at_zeta = self.evaluations_at_zeta(self.zeta);
        let permutation_polynomial_at_zeta_omega = self.polys.permutation_polynomial_at_zeta_omega(self.zeta);

        let proof = BatchyPlonkProof {
            aggregate_kzg_proof_at_zeta,
            permutation_polynomial_kzg_proof_at_zeta_omega,
            evals_at_zeta,
            permutation_polynomial_at_zeta_omega,
        };

        proof_commitments.extend(self.commit_polynomials(ck, &vec![self.linearization_polynomial.clone()]));

        (proof, proof_commitments)
    }

    fn verify(&self, vk: &CS::VK, preprocessed_commitments: Vec<CS::C>, commitments_from_proof: Vec<CS::C>, proof: Self::Proof) -> bool {
        let t_kzg = start_timer!(|| "KZG batch verification");
        let (agg_comm, agg_eval) = {
            let t_aggregate_claims = start_timer!(|| "aggregate evaluation claims at zeta");

            let mut comms = vec![commitments_from_proof[5].clone()];
            comms.extend_from_slice(&commitments_from_proof[0..3]); // wire polynomials commitments: [a]_1, [b]_1, [b]_1
            comms.extend_from_slice(&preprocessed_commitments[5..7]); // [s_{sigma_1}]_1, [s_{sigma_2}]_1
            assert_eq!(comms.len(), self.nus.len());
            let agg_comms = CS::C::combine( &self.nus, &comms);

            let evals = proof.evals_at_zeta;
            assert_eq!(evals.len(), self.nus.len());
            let agg_evals = evals.into_iter().zip(self.nus.iter()).map(|(y, r)| y * r).sum();

            end_timer!(t_aggregate_claims);
            (agg_comms, agg_evals)
        };
        let permutation_polynomial_commitment = commitments_from_proof[3].clone();
        let t_kzg_batch_opening = start_timer!(|| "batched KZG openning");
        // let opening_at_zeta = KzgOpening {
        //     c: agg_comm,
        //     x: self.zeta,
        //     y: agg_eval,
        //     proof: proof.cs_proof.0,
        // };
        // let opening_at_zeta_omega = KzgOpening {
        //     c: permutation_polynomial_commitment,
        //     x: self.zeta * self.polys.omega,
        //     y: proof.r_zeta_omega,
        //     proof: proof.r_at_zeta_omega_proof,
        // };
        // let openings = vec![opening_at_zeta, opening_at_zeta_omega];
        let valid = CS::batch_verify(vk,
                                     vec![agg_comm, permutation_polynomial_commitment],
                                     vec![self.zeta, self.polys.omega],
                                     vec![agg_eval, proof.permutation_polynomial_at_zeta_omega],
                                     vec![proof.aggregate_kzg_proof_at_zeta, proof.permutation_polynomial_kzg_proof_at_zeta_omega],
                                     &mut test_rng());

        end_timer!(t_kzg_batch_opening);
        end_timer!(t_kzg);
        valid
    }
}