#[cfg(test)]
mod tests {
    use crate::{Poly, FflonkyKzg};
    use ark_std::test_rng;
    use ark_poly::{UVPolynomial, Polynomial};
    use crate::utils::poly;
    use ark_std::rand::Rng;
    use ark_ff::{PrimeField, UniformRand};
    use ark_std::convert::TryInto;
    use ark_poly::Radix2EvaluationDomain;
    use ark_poly::EvaluationDomain;
    use crate::tests::{TestKzg, TestCurve, TestField};
    use crate::pcs::{PCS, PcsParams};
    use ark_std::{end_timer, start_timer};
    use crate::fflonk::Fflonk;

    struct VanillaPlonk<F: PrimeField, D: EvaluationDomain<F>> {
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

    impl<F: PrimeField> VanillaPlonk<F, Radix2EvaluationDomain<F>> {
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
            let fs2 = vec![z, t1, t2];
            vec![
                Combination { fs: fs0, xs_roots: vec![zeta] },
                Combination { fs: fs1, xs_roots: vec![zeta] },
                Combination { fs: fs2, xs_roots: vec![zeta, zeta * omega] },
            ]
        }
    }

    struct Combination<F: PrimeField> {
        fs: Vec<Poly<F>>,
        xs_roots: Vec<F>,
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
    }

    struct Prover<F: PrimeField, CS: PCS<F>> {
        ck: CS::CK,
    }

    impl<F: PrimeField, CS: PCS<F>> Prover<F, CS> {
        fn new(ck: CS::CK) -> Self {
            Self { ck }
        }

        fn prove(&self, combinations: Vec<Combination<F>>) {
            let comms = self.commit(&combinations);
            let proof = self.open(combinations);
        }

        fn commit(&self, combinations: &[Combination<F>]) -> Vec<CS::C> {
            let t_preprocessing = start_timer!(|| "Preprocessing");
            let preprocessed_c = self.commit_combination(&combinations[0]);
            end_timer!(t_preprocessing);

            let t_commitment = start_timer!(|| "Commitment");
            let fcs = combinations.iter().map(|c| self.commit_combination(c)).collect::<Vec<_>>();
            end_timer!(t_commitment);

            let mut res = vec![preprocessed_c];
            res.extend(fcs);
            res
        }

        fn commit_combination(&self, combination: &Combination<F>) -> CS::C {
            let t_commit = start_timer!(|| format!("{} polys of degree up to {}", combination.fs.len(), combination.max_degree()));
            let combined_p = Fflonk::combine(combination.t(), &combination.fs);
            let combined_c = CS::commit(&self.ck, &combined_p);
            end_timer!(t_commit, || format!("combined poly degree = {}", combined_p.degree()));
            combined_c
        }

        fn open(&self, combinations: Vec<Combination<F>>)  {
            let rng = &mut test_rng();
            let transcript = &mut (F::rand(rng), F::rand(rng));

            let (ts, (fss, xss)): (Vec<_>, (Vec<_>, Vec<_>)) =
                combinations.into_iter().map(|c| (c.t(), (c.fs, c.xs_roots))).unzip();

            let t_open = start_timer!(|| "Opening");
            let proof = FflonkyKzg::<F, CS>::open(&self.ck, &fss, &ts, &xss, transcript);
            end_timer!(t_open);
        }
    }

    fn _test_vanilla_plonk_opening<F: PrimeField, CS: PCS<F>>(log_n: usize) {
        let rng = &mut test_rng();

        let n = 1 << log_n;

        let piop = VanillaPlonk::<F, _>::new(n, rng);
        let combinations = piop.combinations();

        let urs_degree = combinations.iter().map(|c| c.max_combined_degree()).max().unwrap();

        let t_setup = start_timer!(|| format!("KZG setup of size {} on {}", urs_degree, crate::utils::curve_name::<TestCurve>()));
        let params = CS::setup(urs_degree, rng);
        end_timer!(t_setup);

        let (ck, vk) = (params.ck(), params.vk());

        let prover = Prover::<F, CS>::new(ck);
        prover.prove(combinations);
    }

    #[test]
    fn test_vanilla_plonk_opening() {
        _test_vanilla_plonk_opening::<TestField, TestKzg>(8);
    }
}