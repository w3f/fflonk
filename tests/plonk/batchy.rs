use std::marker::PhantomData;
use ark_ff::{PrimeField, UniformRand};
use fflonk::pcs::{PCS, PcsParams};
use crate::{DecoyPlonk, VanillaPlonkAssignments};
use fflonk::Poly;
use ark_poly::Radix2EvaluationDomain;
use fflonk::shplonk::AggregateProof;
use ark_std::rand::Rng;
use ark_std::{end_timer, start_timer};

impl<F: PrimeField> VanillaPlonkAssignments<F, Radix2EvaluationDomain<F>> {
    // fn combinations(&self) -> Vec<Combination<F>> {
    //     let zeta: F = u128::rand(&mut test_rng()).into();
    //     let omega = self.domain.group_gen;
    //     let t0 = self.quotient(&self.arithmetic_constraint);
    //     let t1 = self.quotient(&self.permutation_constraint_1);
    //     let t2 = self.quotient(&self.permutation_constraint_2);
    //     let z = self.permutation_polynomial.clone();
    //
    //     let fs0 = self.preprocessed_polynomials.clone();
    //     let mut fs1 = self.wire_polynomials.clone();
    //     fs1.push(t0);
    //     let fs2 = vec![z, t1, t2, Poly::zero()]; //TODO: zero is not strictly necessary
    //     vec![
    //         Combination { fs: fs0, roots_of_xs: vec![zeta] },
    //         Combination { fs: fs1, roots_of_xs: vec![zeta] },
    //         Combination { fs: fs2, roots_of_xs: vec![zeta, zeta * omega] },
    //     ]
    // }
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

    fn prove(&mut self, ck: &<CS as PCS<F>>::CK) -> (AggregateProof<F, CS>, Vec<<CS as PCS<F>>::C>, Vec<Vec<Vec<F>>>) {
        todo!()
    }

    fn verify(&self, vk: &<CS as PCS<F>>::VK, preprocessed_commitments: Vec<<CS as PCS<F>>::C>, commitments_from_proof: Vec<<CS as PCS<F>>::C>, evals_from_proof: Vec<Vec<Vec<F>>>, cs_proof: AggregateProof<F, CS>) -> bool {
        todo!()
    }
}