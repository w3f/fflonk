use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_std::rand::RngCore;
use ark_ff::{UniformRand, PrimeField, BitIteratorBE, Zero, BitIteratorLE, FftField, FftParameters};
use crate::utils;
use ark_ec::msm::FixedBase;

use ark_serialize::*;
use ark_std::io::{Read, Write};

use ark_std::{end_timer, start_timer};


/// Updatable Universal References String
#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
// The bases are presented in affine as ark_ec::msm::VariableBaseMSM, though does double-and-add in projective, still enjoys mixed addition,
// and Celo's https://github.com/celo-org/zexe/blob/master/algebra-core/src/curves/batch_arith.rs uses affine ops followed with with batch inversions.
// See https://github.com/arkworks-rs/algebra/issues/60
pub struct URS<E: PairingEngine> {
    // g1, tau.g1, tau^2.g1, ..., tau^n1.g1, where g1 is a generator of G1
    pub powers_in_g1: Vec<E::G1Affine>,
    // g2, tau.g2, tau^2.g2, ..., tau^n2.g2, where g2 is a generator of G2
    pub powers_in_g2: Vec<E::G2Affine>,
}

impl<E: PairingEngine> URS<E> {
    // Multiply the same base by each scalar.
    fn single_base_msm<G: ProjectiveCurve>(scalars: &[G::ScalarField], g: G) -> Vec<G::Affine> {
        let num_scalars = scalars.len();
        let window_size = FixedBase::get_mul_window_size(num_scalars);
        let bits_in_scalar = G::ScalarField::size_in_bits();
        let table = FixedBase::get_window_table(bits_in_scalar, window_size, g);
        let scalars_in_g = FixedBase::msm(bits_in_scalar, window_size, &table, scalars);
        assert_eq!(scalars_in_g.len(), num_scalars);

        G::batch_normalization_into_affine(&scalars_in_g)
    }

    /// Generates URS of the form:
    /// g1, tau.g1, ..., tau^(n1-1).g1, g2, tau.g2, ..., tau^(n2-1).g2
    pub fn generate<R: RngCore>(n1: usize, n2: usize, rng: &mut R) -> Self {
        let tau = E::Fr::rand(rng);
        let n = n1.max(n2);
        assert!(n > 0, "nothing to generate");

        // Until ECFFT for more curves is implemented, see https://github.com/wborgeaud/ecfft-bn254
        assert!(n <= 1 << <E::Fr as FftField>::FftParams::TWO_ADICITY, "number of bases exceeds curve 2-adicity");

        let t_powers = start_timer!(|| "Compute scalars powers");
        let powers_of_tau = utils::powers(tau, n - 1);
        // tau^0, ..., tau^(n-1))
        end_timer!(t_powers);
        assert_eq!(powers_of_tau.len(), n);

        let g1 = E::G1Projective::rand(rng);
        let g2 = E::G2Projective::rand(rng);

        let t_msm_g1 = start_timer!(|| "Scalar muls in G1");
        let powers_in_g1 = Self::single_base_msm(&powers_of_tau[..n1], g1);
        end_timer!(t_msm_g1);

        let t_msm_g2 = start_timer!(|| "Scalar muls in G2");
        let powers_in_g2 = Self::single_base_msm(&powers_of_tau[..n2], g2);
        end_timer!(t_msm_g2);

        URS {
            powers_in_g1,
            powers_in_g2,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_bw6_761::{BW6_761, G1Projective, Fr};
    use ark_std::test_rng;
    use ark_ec::AffineCurve;

    fn _bench_urs_generation<E: PairingEngine>() {
        let log_num_bases = 8;
        let num_bases = 1 << log_num_bases;
        let t_generate = start_timer!(|| format!("Generate 2^{} G1 and 2^{} G2 bases for {}", log_num_bases, log_num_bases, std::any::type_name::<E>()));
        URS::<E>::generate(num_bases, num_bases, &mut test_rng());
        end_timer!(t_generate);
    }

    #[test]
    fn bench_urs_generation() {
        _bench_urs_generation::<BW6_761>();
    }

    #[test]
    fn test_urs_generation() {
        let n1 = 1 << 8;
        let n2 = 2;
        let urs = URS::<BW6_761>::generate(n1, n2, &mut test_rng());
        assert_eq!(urs.powers_in_g1.len(), n1);
        assert_eq!(urs.powers_in_g2.len(), n2);
    }

    #[test]
    #[should_panic]
    fn test_max_bases() {
        let max_bases = 1 << <Fr as FftField>::FftParams::TWO_ADICITY;
        URS::<BW6_761>::generate(max_bases + 1, 0, &mut test_rng());
    }
}