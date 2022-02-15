mod fflonky;
mod batchy;


use fflonk::{Poly, FflonkyKzg};
use ark_std::test_rng;
use ark_poly::{UVPolynomial, Polynomial};
use ark_std::rand::Rng;
use ark_ff::{PrimeField, UniformRand, Zero};
use ark_poly::EvaluationDomain;
use fflonk::pcs::{PCS, PcsParams};
use ark_std::{end_timer, start_timer};
use fflonk::fflonk::Fflonk;
use std::marker::PhantomData;
use ark_poly::Radix2EvaluationDomain;
use ark_ec::PairingEngine;
use fflonk::shplonk::AggregateProof;

struct VanillaPlonkAssignments<F: PrimeField, D: EvaluationDomain<F>> {
    domain: D,
    // [Poly<F>; 8], max_deg = d
    preprocessed_polynomials: Vec<Poly<F>>,
    // [Poly<F>; 3], max_deg = d
    wire_polynomials: Vec<Poly<F>>,
    // max_deg = d
    permutation_polynomial: Poly<F>,
    // max_deg = 3 * d
    arithmetic_constraint: Poly<F>,
    // max_deg = 2 * d
    permutation_constraint_1: Poly<F>,
    // // max_deg = 4 * d
    permutation_constraint_2: Poly<F>,
}

fn random_polynomials<F: PrimeField, R: Rng>(k: usize, degree: usize, rng: &mut R) -> Vec<Poly<F>> {
    (0..k).map(|_| Poly::rand(degree, rng)).collect()
}

impl<F: PrimeField> VanillaPlonkAssignments<F, Radix2EvaluationDomain<F>> {
    fn new<R: Rng>(n: usize, rng: &mut R) -> Self {
        let d = n - 1;
        let domain = Radix2EvaluationDomain::<F>::new(n).expect("TODO"); //TODO
        Self {
            domain,
            preprocessed_polynomials: random_polynomials(8, d, rng),
            wire_polynomials: random_polynomials(3, d, rng),
            permutation_polynomial: Poly::rand(d, rng),
            arithmetic_constraint: Poly::rand(3 * d, rng),
            permutation_constraint_1: Poly::rand(2 * d, rng),
            permutation_constraint_2: Poly::rand(4 * d, rng),
        }
    }

    fn quotient(&self, constraint: &Poly<F>) -> Poly<F> {
        constraint.divide_by_vanishing_poly(self.domain).unwrap().0
    }
}

trait DecoyPlonk<F: PrimeField, CS: PCS<F>> {
    fn setup<R: Rng>(&mut self, rng: &mut R) -> (CS::CK, CS::VK);
    fn preprocess(&mut self, ck: &CS::CK) -> Vec<CS::C>;
    fn prove(&mut self, ck: &CS::CK) -> (AggregateProof<F, CS>, Vec<CS::C>, Vec<Vec<Vec<F>>>);
    fn verify(&self, vk: &CS::VK, preprocessed_commitments: Vec<CS::C>, commitments_from_proof: Vec<CS::C>, evals_from_proof: Vec<Vec<Vec<F>>>, cs_proof: AggregateProof<F, CS>) -> bool;
}


