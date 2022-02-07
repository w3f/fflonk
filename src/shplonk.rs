use ark_ff::PrimeField;
use ark_poly::{UVPolynomial, Polynomial};

use std::collections::HashSet;
use std::marker::PhantomData;

use crate::Poly;
use crate::pcs::PCS;
use crate::pcs::aggregation::{aggregate_polys, aggregate_claims, Claim};


pub trait ShplonkTranscript<F, G> {
    fn get_gamma(&mut self) -> F;
    fn commit_to_q(&mut self, q_comm: &G);
    fn get_zeta(&mut self) -> F;
}


pub struct Shplonk<F: PrimeField, CS: PCS<F>> {
    _field: PhantomData<F>,
    _pcs: PhantomData<CS>,
}

impl<F: PrimeField, CS: PCS<F>> Shplonk<F, CS> {
    pub fn open_many<T: ShplonkTranscript<F, CS::G>>(
        ck: &CS::CK,
        fs: &[Poly<F>],
        xss: &[HashSet<F>],
        transcript: &mut T,
    ) -> (CS::G, CS::Proof)
    {
        let (agg_poly, zeta, agg_proof) = aggregate_polys::<F, CS, T>(ck, fs, xss, transcript);
        assert!(agg_poly.evaluate(&zeta).is_zero());
        let opening_proof = CS::open(ck, &agg_poly, zeta);
        (agg_proof, opening_proof)
    }

    pub fn verify_many<T: ShplonkTranscript<F, CS::G>>(
        vk: &CS::VK,
        fcs: &[CS::G],
        proof: (CS::G, CS::Proof),
        xss: &Vec<Vec<F>>,
        yss: &Vec<Vec<F>>,
        transcript: &mut T,
    ) -> bool
    {
        let (agg_proof, opening_proof) = proof;
        let onec = CS::commit(&vk.clone().into(), &Poly::from_coefficients_slice(&[F::one()]));
        let claims = Self::group_by_commitment(fcs, xss, yss);
        let agg_claim = aggregate_claims::<F, CS, T>(claims, &agg_proof, &onec, transcript);
        CS::verify(vk, agg_claim.c, agg_claim.xs[0], agg_claim.ys[0], opening_proof)
    }

    fn group_by_commitment(
        fcs: &[CS::G],
        xss: &Vec<Vec<F>>,
        yss: &Vec<Vec<F>>,
    ) -> Vec<Claim<F, CS::G>> {
        fcs.iter().cloned()
            .zip(xss.iter().cloned())
            .zip(yss.iter().cloned())
            .map(|((c, xs), ys)| Claim { c, xs, ys })
            .collect()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::pcs::tests::IdentityCommitment;

    use ark_std::test_rng;
    use ark_std::iter::FromIterator;
    use ark_std::rand::Rng;
    use crate::Poly;
    use ark_bw6_761::{BW6_761, Fr};
    use crate::pcs::kzg::KZG;
    use crate::pcs::PcsParams;

    impl<F: PrimeField, G> ShplonkTranscript<F, G> for (F, F) {
        fn get_gamma(&mut self) -> F { self.0 }

        fn commit_to_q(&mut self, _q_comm: &G) {}

        fn get_zeta(&mut self) -> F { self.1 }
    }

    pub fn generate_test_data<R, F, CS>(
        rng: &mut R,
        ck: &CS::CK,
        d: usize, // degree of polynomials
        t: usize, // number of polynomials
        xss: &Vec<Vec<F>>, // vecs of opening points per polynomial
    ) -> (
        Vec<Poly<F>>, // polynomials
        Vec<CS::G>, // commitments
        Vec<Vec<F>>, // evaluations per polynomial
    ) where
        R: Rng,
        F: PrimeField,
        CS: PCS<F>,
    {
        // polynomials
        let fs: Vec<_> = (0..t)
            .map(|_| Poly::<F>::rand(d, rng))
            .collect();
        // commitments
        let fcs: Vec<_> = fs.iter()
            .map(|fi| CS::commit(&ck, fi))
            .collect();

        // evaluations per polynomial
        let yss: Vec<_> = fs.iter()
            .zip(xss.iter())
            .map(|(f, xs)|
                xs.iter().map(
                    |x| f.evaluate(x))
                    .collect::<Vec<_>>()
            ).collect();

        (fs, fcs, yss)
    }

    fn _test_shplonk<F: PrimeField, CS: PCS<F>>() {
        let rng = &mut test_rng();

        let params = CS::setup(123, rng);

        let t = 4; // number of polynomials
        let max_m = 3; // maximal number of opening points per polynomial

        let xss: Vec<_> = (0..t)
            .map(|_| (0..rng.gen_range(1..max_m))
                .map(|_| F::rand(rng)).collect::<Vec<_>>())
            .collect();

        let (fs, fcs, yss) =
            generate_test_data::<_, _, CS>(rng, &params.ck(), 15, 4, &xss);

        let sets_of_xss: Vec<HashSet<F>> = xss.iter()
            .map(|xs| HashSet::from_iter(xs.iter().cloned()))
            .collect();

        let transcript = &mut (F::rand(rng), F::rand(rng));

        let (qc, qlc) = Shplonk::<F, CS>::open_many(&params.ck(), &fs, sets_of_xss.as_slice(), transcript);

        assert!(Shplonk::<F, CS>::verify_many(&params.vk(), &fcs, (qc, qlc), &xss, &yss, transcript))
    }

    #[test]
    fn test_shplonk_id() {
        _test_shplonk::<Fr, IdentityCommitment>();
    }

    #[test]
    fn test_shplonk_kzg() {
        _test_shplonk::<Fr, KZG<BW6_761>>();
    }
}