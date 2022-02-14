#[cfg(test)]
mod tests {
    use crate::Poly;
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
                Combination { fs: fs0, xs: vec![zeta] },
                Combination { fs: fs1, xs: vec![zeta] },
                Combination { fs: fs2, xs: vec![zeta, zeta * omega] },
            ]
        }
    }

    struct Combination<F: PrimeField> {
        fs: Vec<Poly<F>>,
        xs: Vec<F>,
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


    #[test]
    fn test_vanilla_plonk_opening() {
        let rng = &mut test_rng();

        let log_n = 8;
        let n = 1 << log_n;

        let piop = VanillaPlonk::new(n, rng);
        let combinations = piop.combinations();

        let urs_degree = combinations.iter().map(|c| c.max_combined_degree()).max().unwrap();

        let t_setup = start_timer!(|| format!("KZG setup of size {} on {}", urs_degree, crate::utils::curve_name::<TestCurve>()));
        let urs = TestKzg::setup(urs_degree, rng);
        end_timer!(t_setup);

        let (ck, vk) = (urs.ck(), urs.vk());

        let preprocessed = &combinations[0];
        let t_commit_preprocessed = start_timer!(|| format!("Preprocessing: committing to the combination of {} polynomials of degree up to {}", preprocessed.fs.len(), preprocessed.max_degree()));
        let preprocessed_combined = Fflonk::combine(preprocessed.t(), &preprocessed.fs);
        let preprocessed_combined_c = TestKzg::commit(&ck, &preprocessed_combined);
        end_timer!(t_commit_preprocessed, || format!("combined polynomial degree = {}!", preprocessed_combined.degree()));

        let t_commit = start_timer!(|| format!("Generating commitments to {} combinations", combinations[1..].len()));
        for combination in combinations {
            let t_commit_combination = start_timer!(|| format!("committing to to the combination of {} polynomials of degree up to {}", combination.fs.len(), combination.max_degree()));
            let combined = Fflonk::combine(combination.t(), &combination.fs);
            let combined_c = TestKzg::commit(&ck, &combined);
            end_timer!(t_commit_combination, || format!("combined polynomial degree = {}!", combined.degree()));
        }
        end_timer!(t_commit);
    }
}