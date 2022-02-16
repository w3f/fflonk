use ark_ff::{PrimeField, UniformRand, Zero};
use ark_poly::{Radix2EvaluationDomain, Polynomial};
use ark_std::test_rng;
use fflonk::{Poly, FflonkyKzg};
use crate::{VanillaPlonkAssignments, DecoyPlonk};
use fflonk::pcs::PCS;
use ark_std::rand::Rng;
use fflonk::fflonk::Fflonk;
use ark_std::{end_timer, start_timer};
use fflonk::pcs::PcsParams;
use fflonk::shplonk::AggregateProof;
use std::marker::PhantomData;


impl<F: PrimeField> VanillaPlonkAssignments<F> {
    fn combinations(&self) -> Vec<Combination<F>> {
        let zeta: F = u128::rand(&mut test_rng()).into();
        let omega = self.domain.group_gen;
        let t0 = self.quotient(&self.arithmetic_constraint);
        let t1 = self.quotient(&self.permutation_constraint_1);
        let t2 = self.quotient(&self.permutation_constraint_2);
        let z = self.permutation_polynomial.clone();

        let fs0 = self.preprocessed_polynomials.clone();
        let mut fs1 = self.wire_polynomials.clone();
        fs1.push(t0);
        let fs2 = vec![z, t1, t2, Poly::zero()]; //TODO: zero is not strictly necessary
        vec![
            Combination { fs: fs0, roots_of_xs: vec![zeta] },
            Combination { fs: fs1, roots_of_xs: vec![zeta] },
            Combination { fs: fs2, roots_of_xs: vec![zeta, zeta * omega] },
        ]
    }
}

struct Combination<F: PrimeField> {
    fs: Vec<Poly<F>>,
    roots_of_xs: Vec<F>,
}

impl<F: PrimeField> Combination<F> {
    fn max_degree(&self) -> usize {
        self.fs.iter().map(|f| f.degree()).max().unwrap()
    }

    fn t(&self) -> usize {
        self.fs.len().next_power_of_two() //TODO: should work fine for other roots
    }

    fn max_combined_degree(&self) -> usize {
        self.t() * (self.max_degree() + 1) - 1
    }

    fn xs(&self) -> Vec<F> {
        self.roots_of_xs.iter() // opening points
            .map(|root| root.pow([self.t() as u64]))
            .collect()
    }

    fn yss(&self) -> Vec<Vec<F>> {
        self.xs().iter().map(|x|
            self.fs.iter().map(|f| f.evaluate(x)).collect()
        ).collect()
    }
}

struct PlonkWithFflonkTest<F: PrimeField, CS: PCS<F>> {
    combinations: Vec<Combination<F>>,
    cs: PhantomData<CS>,
}

impl<F: PrimeField, CS: PCS<F>> PlonkWithFflonkTest<F, CS> {
    fn new(combinations: Vec<Combination<F>>) -> Self {
        Self {
            combinations,
            cs: PhantomData,
        }
    }

    fn _commit_proof_polynomials(&self, ck: &CS::CK) -> Vec<CS::C> {
        let t_commitment = start_timer!(|| format!("Committing to {} proof polynomials", self.combinations.len() - 1));
        let commitments = self.combinations.iter().enumerate()
            .skip(1) // preprocessing
            .map(|(i, _)| self._commit_single(i, ck))
            .collect();
        end_timer!(t_commitment);
        commitments
    }

    fn _commit_single(&self, i: usize, ck: &CS::CK) -> CS::C {
        let combination = &self.combinations[i];
        let t_commit = start_timer!(|| format!("Committing to combination #{}", i));

        let t_combine = start_timer!(|| format!("combining {} polynomials: t = {}, max_degree = {}", combination.fs.len(), combination.t(), combination.max_degree()));
        let poly = Fflonk::combine(combination.t(), &combination.fs);
        end_timer!(t_combine);

        let t_commit_combined = start_timer!(|| format!("committing to the combined polynomial: degree = {}", poly.degree()));
        let commitment = CS::commit(ck, &poly);
        end_timer!(t_commit_combined);

        end_timer!(t_commit);
        commitment
    }

    fn _open(&self, transcript: &mut merlin::Transcript, ck: &CS::CK) -> AggregateProof<F, CS> {
        let (ts, (fss, xss)): (Vec<_>, (Vec<_>, Vec<_>)) =
            self.combinations.iter()
                .map(|c| (c.t(), (c.fs.clone(), c.roots_of_xs.clone())))
                .unzip();

        let t_open = start_timer!(|| "Opening");
        let proof = FflonkyKzg::<F, CS>::open(ck, &fss, &ts, &xss, transcript);
        end_timer!(t_open);
        proof
    }

    fn _evaluate(&self) -> Vec<Vec<Vec<F>>> {
        self.combinations.iter()
            .map(|c| c.yss()).collect()
    }
}

struct FflonkyPlonkProof<F: PrimeField, CS: PCS<F>> {
    cs_proof: AggregateProof<F, CS>,
    evals: Vec<Vec<Vec<F>>>,
}


impl<F: PrimeField, CS: PCS<F>> DecoyPlonk<F, CS> for PlonkWithFflonkTest<F, CS> {
    type Proof = FflonkyPlonkProof<F, CS>;

    fn setup<R: Rng>(&mut self, rng: &mut R) -> (CS::CK, CS::VK) {
        let urs_degree = self.combinations.iter()
            .map(|c| c.max_combined_degree())
            .max().unwrap();

        let t_setup = start_timer!(|| format!("KZG setup of degree {} on {}",
                urs_degree, fflonk::utils::curve_name::<TestCurve>()));
        let params = CS::setup(urs_degree, rng);
        end_timer!(t_setup);
        (params.ck(), params.vk())
    }

    fn preprocess(&mut self, ck: &CS::CK) -> Vec<CS::C> {
        let t_preprocessing = start_timer!(|| "Preprocessing");
        let commitment = self._commit_single(0, ck);
        end_timer!(t_preprocessing);
        vec![commitment]
    }

    fn prove(&mut self, ck: &CS::CK) -> (FflonkyPlonkProof<F, CS>, Vec<CS::C>) {
        let empty_transcript = &mut merlin::Transcript::new(b"plonk-fflonk-shplonk-kzg");
        let t_proving = start_timer!(|| "Proving");
        let commitments = self._commit_proof_polynomials(ck);
        let cs_proof = self._open(empty_transcript, ck);
        end_timer!(t_proving);
        let evals = self._evaluate();
        (FflonkyPlonkProof { cs_proof, evals }, commitments)
    }

    fn verify(&self, vk: &CS::VK, preprocessed_commitments: Vec<CS::C>, commitments_from_proof: Vec<CS::C>, proof: Self::Proof) -> bool {
        let empty_transcript = &mut merlin::Transcript::new(b"plonk-fflonk-shplonk-kzg");
        let (ts, xss): (Vec<_>, Vec<_>) =
            self.combinations.iter()
                .map(|c| (c.t(), c.roots_of_xs.clone()))
                .unzip();

        let commitments = [preprocessed_commitments, commitments_from_proof].concat();
        FflonkyKzg::<F, CS>::verify(vk, &commitments, &ts, proof.cs_proof, &xss, &proof.evals, empty_transcript)
    }
}

fn _test_vanilla_plonk_opening<F: PrimeField, CS: PCS<F>>(log_n: usize) {
    let rng = &mut test_rng();

    let n = 1 << log_n;

    let piop = VanillaPlonkAssignments::<F>::new(n, rng);
    let combinations = piop.combinations();

    let mut test = PlonkWithFflonkTest::<F, CS>::new(combinations);

    let (ck, vk) = test.setup(rng);
    let commitments_to_preprocessed_polynomials = test.preprocess(&ck);
    let (cs_proof, commitments) = test.prove(&ck);
    // assert_eq!(self.evals.len(), self.combinations.len());
    // assert_eq!(self.evals[0].len(), self.combinations[0].roots_of_xs.len());
    // assert_eq!(self.evals[0][0].len(), self.combinations[0].fs.len());
    assert!(test.verify(&vk, commitments_to_preprocessed_polynomials, commitments, cs_proof));
}

#[test]
fn test_vanilla_plonk_opening2() {
    assert!(cfg!(test));
    _test_vanilla_plonk_opening::<_, fflonk::pcs::kzg::KZG<ark_bls12_381::Bls12_381>>(8);
}
