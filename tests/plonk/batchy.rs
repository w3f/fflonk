use std::marker::PhantomData;
use ark_ff::{PrimeField, UniformRand};
use fflonk::pcs::{PCS, PcsParams};
use crate::{DecoyPlonk, VanillaPlonkAssignments, random_polynomials};
use fflonk::Poly;
use ark_poly::{Radix2EvaluationDomain, Polynomial, EvaluationDomain, UVPolynomial};
use fflonk::shplonk::AggregateProof;
use ark_std::rand::Rng;
use ark_std::{end_timer, start_timer, test_rng};
use fflonk::utils::poly;

impl<F: PrimeField> VanillaPlonkAssignments<F, Radix2EvaluationDomain<F>> {
    fn constraints(&self) -> Vec<Poly<F>> {
        vec![
            self.arithmetic_constraint.clone(),
            self.permutation_constraint_1.clone(),
            self.permutation_constraint_2.clone(),
        ]
    }

    fn aggregate(constraints: &Vec<Poly<F>>) -> Poly<F> {
        let gamma: F = u128::rand(&mut test_rng()).into();
        poly::sum_with_powers(gamma, constraints)
    }

    fn aggregated_constraints_quotient(&self) -> Poly<F> {
        let aggregate_constraint = Self::aggregate(&self.constraints());
        self.quotient(&aggregate_constraint)
    }

    //TODO
    fn evaluations(&self) -> Vec<F> {
        vec![]
    }

    //TODO
    fn linearization_polynomial(&self) -> Poly<F> {
        let d = self.domain.size() - 1;
        Poly::rand(d, &mut test_rng())
    }

    fn zeta_omega(&self) -> (F, F) {
        let zeta = u128::rand(&mut test_rng()).into();
        let omega = self.domain.group_gen;
        (zeta, omega)
    }

    //TODO
    fn aggregated_polynomial_to_open_in_zeta(&self) -> Poly<F> {
        let d = self.domain.size() - 1;
        Poly::rand(d, &mut test_rng())
    }
}

struct PlonkBatchKzgTest<F: PrimeField, CS: PCS<F>> {
    polys: VanillaPlonkAssignments<F, Radix2EvaluationDomain<F>>,
    cs: PhantomData<CS>,
}

impl<F: PrimeField, CS: PCS<F>> PlonkBatchKzgTest<F, CS> {
    fn commit_polynomials(&self, ck: &CS::CK, polys: &Vec<Poly<F>>) -> Vec<CS::C> {
        let t_commitment = start_timer!(|| format!("Committing to {} polynomials", polys.len()));
        let commitments = polys.iter().map(|p| CS::commit(ck, &p)).collect();
        end_timer!(t_commitment);

        commitments
    }
}

impl<F: PrimeField, CS: PCS<F>> DecoyPlonk<F, CS> for PlonkBatchKzgTest<F, CS> {
    type Proof = (CS::Proof, CS::Proof);

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

    fn prove(&mut self, ck: &CS::CK) -> (Self::Proof, Vec<CS::C>, Vec<Vec<Vec<F>>>) {
        let empty_transcript = &mut merlin::Transcript::new(b"plonk-batch-kzg");
        let t_proving = start_timer!(|| "Proving");
        let mut proof_commitments = vec![];
        proof_commitments.extend(self.commit_polynomials(ck, &self.polys.wire_polynomials));
        proof_commitments.extend(self.commit_polynomials(ck, &vec![self.polys.permutation_polynomial.clone()]));
        let q = self.polys.aggregated_constraints_quotient();
        assert_eq!(q.degree() as u64, 3 * self.polys.domain.size - 4);
        proof_commitments.extend(self.commit_polynomials(ck, &vec![q]));

        let (zeta, omega) = self.polys.zeta_omega();
        let opening_in_zeta = CS::open(ck, &self.polys.aggregated_polynomial_to_open_in_zeta(), zeta);
        let opening_in_zeta_omega = CS::open(ck, &self.polys.permutation_polynomial, zeta * omega);
        end_timer!(t_proving);
        ((opening_in_zeta, opening_in_zeta_omega), proof_commitments, vec![])
    }

    fn verify(&self, vk: &<CS as PCS<F>>::VK, preprocessed_commitments: Vec<<CS as PCS<F>>::C>, commitments_from_proof: Vec<<CS as PCS<F>>::C>, evals_from_proof: Vec<Vec<Vec<F>>>, cs_proof: Self::Proof) -> bool {
        todo!()
    }
}