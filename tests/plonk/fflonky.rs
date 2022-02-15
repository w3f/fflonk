use ark_ff::{PrimeField, UniformRand, Zero};
use ark_poly::{Radix2EvaluationDomain, Polynomial};
use ark_std::test_rng;
use fflonk::{Poly, FflonkyKzg};
use crate::{VanillaPlonkAssignments, PlonkTest};
use fflonk::pcs::PCS;
use ark_std::rand::Rng;
use fflonk::fflonk::Fflonk;
use ark_std::{end_timer, start_timer};
use fflonk::pcs::PcsParams;


impl<F: PrimeField> VanillaPlonkAssignments<F, Radix2EvaluationDomain<F>> {
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
    ck: Option<CS::CK>,
    vk: Option<CS::VK>,
    commitments: Vec<CS::C>,
    proof: Option<(CS::C, CS::Proof)>,
    evals: Vec<Vec<Vec<F>>>,
}

impl<F: PrimeField, CS: PCS<F>> PlonkWithFflonkTest<F, CS> {
    fn new(combinations: Vec<Combination<F>>) -> Self {
        Self {
            combinations,
            ck: None,
            vk: None,
            commitments: vec![],
            proof: None,
            evals: vec![],
        }
    }

    fn _commit_proof_polynomials(&self) -> Vec<CS::C> {
        let t_commitment = start_timer!(|| format!("Committing to {} proof polynomials", self.combinations.len() - 1));
        let commitments = self.combinations.iter().enumerate()
            .skip(1) // preprocessing
            .map(|(i, _)| self._commit_single(i))
            .collect();
        end_timer!(t_commitment);
        commitments
    }

    fn _commit_single(&self, i: usize) -> CS::C {
        let combination = &self.combinations[i];
        let t_commit = start_timer!(|| format!("Committing to combination #{}", i));

        let t_combine = start_timer!(|| format!("combining {} polynomials: t = {}, max_degree = {}", combination.fs.len(), combination.t(), combination.max_degree()));
        let poly = Fflonk::combine(combination.t(), &combination.fs);
        end_timer!(t_combine);

        let t_commit_combined = start_timer!(|| format!("committing to the combined polynomial: degree = {}", poly.degree()));
        let commitment = CS::commit(self.ck.as_ref().unwrap(), &poly);
        end_timer!(t_commit_combined);

        end_timer!(t_commit);
        commitment
    }

    fn _open(&self, transcript: &mut merlin::Transcript) -> (CS::C, CS::Proof) {
        let (ts, (fss, xss)): (Vec<_>, (Vec<_>, Vec<_>)) =
            self.combinations.iter()
                .map(|c| (c.t(), (c.fs.clone(), c.roots_of_xs.clone())))
                .unzip();

        let t_open = start_timer!(|| "Opening");
        let proof = FflonkyKzg::<F, CS>::open(&self.ck.as_ref().unwrap(), &fss, &ts, &xss, transcript);
        end_timer!(t_open);
        proof
    }

    fn _evaluate(&self) -> Vec<Vec<Vec<F>>> {
        self.combinations.iter()
            .map(|c| c.yss()).collect()
    }
}

impl<F: PrimeField, CS: PCS<F>> PlonkTest for PlonkWithFflonkTest<F, CS> {
    fn setup<R: Rng>(&mut self, rng: &mut R) {
        let urs_degree = self.combinations.iter()
            .map(|c| c.max_combined_degree())
            .max().unwrap();

        let t_setup = start_timer!(|| format!("KZG setup of degree {} on {}",
                urs_degree, fflonk::utils::curve_name::<TestCurve>()));
        let params = CS::setup(urs_degree, rng);
        end_timer!(t_setup);

        self.ck = Some(params.ck());
        self.vk = Some(params.vk());
    }

    fn preprocess(&mut self) {//-> CS::C {
        let t_preprocessing = start_timer!(|| "Preprocessing");
        let commitment = self._commit_single(0);
        end_timer!(t_preprocessing);

        self.commitments.push(commitment.clone());
        // commitment
    }

    fn prove(&mut self) {//-> (&Vec<CS::C>, &(CS::C, CS::Proof), &Vec<Vec<Vec<F>>>) {
        let empty_transcript = &mut merlin::Transcript::new(b"plonk-fflonk-shplonk-kzg");
        let t_proving = start_timer!(|| "Proving");
        let commitments = self._commit_proof_polynomials();
        let proof = self._open(empty_transcript);
        end_timer!(t_proving);
        self.evals = self._evaluate();

        assert_eq!(self.evals.len(), self.combinations.len());
        assert_eq!(self.evals[0].len(), self.combinations[0].roots_of_xs.len());
        assert_eq!(self.evals[0][0].len(), self.combinations[0].fs.len());

        self.commitments.extend_from_slice(&commitments);
        self.proof = Some(proof);
        // (&self.commitments, self.proof.as_ref().unwrap(), &self.evals)
    }

    fn verify(&self) {
        let empty_transcript = &mut merlin::Transcript::new(b"plonk-fflonk-shplonk-kzg");
        let (ts, xss): (Vec<_>, Vec<_>) =
            self.combinations.iter()
                .map(|c| (c.t(), c.roots_of_xs.clone()))
                .unzip();

        let result = FflonkyKzg::<F, CS>::verify(&self.vk.as_ref().unwrap(), &self.commitments, &ts, self.proof.as_ref().unwrap().clone(), &xss, &self.evals, empty_transcript);
        assert!(result);
    }
}

fn _test_vanilla_plonk_opening<F: PrimeField, CS: PCS<F>>(log_n: usize) {
    let rng = &mut test_rng();

    let n = 1 << log_n;

    let piop = VanillaPlonkAssignments::<F, _>::new(n, rng);
    let combinations = piop.combinations();

    let mut test = PlonkWithFflonkTest::<F, CS>::new(combinations);

    test.setup(rng);
    test.preprocess();
    test.prove();
    test.verify();
}

#[test]
fn test_vanilla_plonk_opening2() {
    assert!(cfg!(test));
    _test_vanilla_plonk_opening::<_, fflonk::pcs::kzg::KZG<ark_bls12_381::Bls12_381>>(8);
}
