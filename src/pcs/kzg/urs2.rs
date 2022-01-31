use ark_ec::PairingEngine;

#[derive(Clone, Debug)]
pub struct ConstUrs<E: PairingEngine, const N1: usize, const N2: usize> {
    // g1, tau.g1, tau^2.g1, ..., tau^n1.g1, where g1 is a generator of G1
    pub powers_in_g1: [E::G1Affine; N1],
    // g2, tau.g2, tau^2.g2, ..., tau^n2.g2, where g2 is a generator of G2
    pub powers_in_g2: [E::G2Affine; N2],
}

impl<E: PairingEngine, const N1: usize, const N2: usize> ConstUrs<E, N1, N2> {
    fn trim_to<const N: usize>(&self) -> ConstUrs<E, N, N> {
        let mut powers_in_g1 = [E::G1Affine::default(); N];
        let mut powers_in_g2 = [E::G2Affine::default(); N];
        powers_in_g1.copy_from_slice(&self.powers_in_g1[..N]);
        powers_in_g2.copy_from_slice(&self.powers_in_g2[..N]);
        ConstUrs {
            powers_in_g1,
            powers_in_g2,
        }
    }
}

