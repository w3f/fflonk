#[cfg(test)]
mod tests {
    use crate::Poly;
    use ark_std::test_rng;
    use ark_poly::UVPolynomial;
    use crate::utils::poly;
    use ark_std::rand::Rng;
    use ark_ff::{PrimeField, UniformRand};
    use ark_std::convert::TryInto;
    use ark_poly::Radix2EvaluationDomain;
    use ark_poly::EvaluationDomain;
    use crate::tests::{TestKzg, TestCurve, TestField};
    use crate::pcs::{PCS, PcsParams};
    use ark_std::{end_timer, start_timer};

    struct VanillaPlonk<F: PrimeField> {
        d: usize,
        preprocessed_polynomials: [Poly<F>; 8],
        // max_deg = d
        wire_polynomials: [Poly<F>; 3],
        // max_deg = d
        permutation_polynomial: Poly<F>,
        // max_deg = d
        arithmetic_constraint: Poly<F>,
        // max_deg = 3 * d
        permutation_constraint_1: Poly<F>,
        // max_deg = 2 * d
        permutation_constraint_2: Poly<F>, // max_deg = 4 * d
    }

    fn random_polynomials<F: PrimeField, R: Rng, const N: usize>(degree: usize, rng: &mut R) -> [Poly<F>; N] {
        (0..N).map(|_| Poly::rand(degree, rng)).collect::<Vec<_>>().try_into().unwrap()
    }

    impl<F: PrimeField> VanillaPlonk<F> {
        fn new<R: Rng>(d: usize, rng: &mut R) -> Self {
            Self {
                d,
                preprocessed_polynomials: random_polynomials::<_, _, 8>(d, rng),
                wire_polynomials: random_polynomials::<_, _, 3>(d, rng),
                permutation_polynomial: Poly::rand(d, rng),
                arithmetic_constraint: Poly::rand(3 * d, rng),
                permutation_constraint_1: Poly::rand(2 * d, rng),
                permutation_constraint_2: Poly::rand(4 * d, rng),
            }
        }

        fn get_constraint_polynomials(&self) -> [Poly<F>; 3] {
            vec![
                &self.arithmetic_constraint,
                &self.permutation_constraint_1,
                &self.permutation_constraint_2,
            ].into_iter().cloned().collect::<Vec<_>>().try_into().unwrap()
        }

        fn get_quotient_polynomial(&self, constraint_polynomial: &Poly<F>) -> Poly<F> {
            let domain = Radix2EvaluationDomain::new(8).unwrap(); //TODO
            constraint_polynomial.divide_by_vanishing_poly(domain).unwrap().0
        }

        fn get_quotient_polynomials<const N: usize>(&self, constraint_polynomials: &[Poly<F>; N]) -> [Poly<F>; N] {
            constraint_polynomials.iter().map(|p| self.get_quotient_polynomial(p))
                .collect::<Vec<_>>().try_into().unwrap()
        }

        fn agg_quotient_polynomial(&self, constraint_polynomials: &[Poly<F>]) -> Poly<F> {
            let alpha = u128::rand(&mut test_rng()).into(); //TODO
            let agg_constraint_polynomial = poly::sum_with_powers(alpha, constraint_polynomials);
            self.get_quotient_polynomial(&agg_constraint_polynomial)
        }
    }


    trait CommitmentScheme {}

    trait Piop<F: PrimeField, CS: CommitmentScheme, const R: usize> {
        fn get_polynomials_to_commit_per_round(&self) -> [Vec<Poly<F>>; R];
    }

    struct BatchedKzg {}

    impl CommitmentScheme for BatchedKzg {}

    impl<F: PrimeField> Piop<F, BatchedKzg, 4> for VanillaPlonk<F> {
        fn get_polynomials_to_commit_per_round(&self) -> [Vec<Poly<F>>; 4] {
            let constraints = self.get_constraint_polynomials();
            let t = self.agg_quotient_polynomial(&constraints);
            let z = self.permutation_polynomial.clone();
            [
                self.preprocessed_polynomials.to_vec(),
                self.wire_polynomials.to_vec(),
                vec![z],
                vec![t],
            ]
        }
    }

    struct Fflonk {}

    impl CommitmentScheme for Fflonk {}

    impl<F: PrimeField> Piop<F, Fflonk, 3> for VanillaPlonk<F> {
        fn get_polynomials_to_commit_per_round(&self) -> [Vec<Poly<F>>; 3] {
            let constraints = self.get_constraint_polynomials();
            let quotients = self.get_quotient_polynomials(&constraints);
            let [t0, t1, t2] = self.get_quotient_polynomials(&constraints);
            let z = self.permutation_polynomial.clone();
            [
                self.preprocessed_polynomials.to_vec(),
                [self.wire_polynomials.as_slice(), &[t0]].concat(),
                [z, t1, t2].to_vec(),
            ]
        }
    }

    #[test]
    fn test_vanilla_plonk_opening() {
        let rng = &mut test_rng();
        let log_n = 8;
        let max_degree = (1 << log_n) - 1;

        let t_setup = start_timer!(|| format!("KZG setup of size 2^{} on {}", log_n, crate::utils::curve_name::<E>()));
        let urs = TestKzg::setup(max_degree, rng);
        end_timer!(t_setup);

        let (ck, vk) = (urs.ck(), urs.vk());

        let polys = VanillaPlonk::<TestField>::new(max_degree, rng);

        for round in Piop::<_, BatchedKzg, 4>::get_polynomials_to_commit_per_round(&polys) {
            for poly in round {
                let c = TestKzg::commit(&ck, &poly);
            }
        }
    }
}