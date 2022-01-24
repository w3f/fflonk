use ark_std::Zero;
use ark_ff::PrimeField;
use ark_poly::{UVPolynomial, Polynomial};
use ark_poly::univariate::DensePolynomial;

use std::collections::HashSet;

use crate::EuclideanPolynomial;
use crate::pcs::{PCS, CommitmentSpace, VerifierKey};


use std::marker::PhantomData;
use crate::Poly;


pub trait ShplonkTranscript<F, G> {
    fn get_gamma(&mut self) -> F;
    fn commit_to_q(&mut self, q_comm: &G);
    fn get_zeta(&mut self) -> F;
}


pub struct Shplonk<F: PrimeField, CS: PCS<F>> {
    _field: PhantomData<F>,
    _pcs: PhantomData<CS>,
}

impl<F: PrimeField, CS: PCS<F>> PCS<F> for Shplonk<F, CS> {
    type G = CS::G;
    type CK = CS::CK;
    type VK = CS::VK;
    type Proof = CS::Proof;

    fn setup() -> (Self::CK, Self::VK) {
        CS::setup()
    }

    fn commit(ck: &Self::CK, p: &Poly<F>) -> Self::G {
        CS::commit(ck, p)
    }

    fn open(ck: &Self::CK, p: &Poly<F>, x: &F) -> Self::Proof {
        CS::open(ck, p, x)
    }

    fn verify(vk: &Self::VK, c: &Self::G, x: &F, z: &F, proof: &Self::Proof) -> bool {
        CS::verify(vk, c, x, z, proof)
    }
}

impl<F: PrimeField, CS: PCS<F>> Shplonk<F, CS> {
    pub fn open_many<T: ShplonkTranscript<F, CS::G>>(
        ck: &<Shplonk<F, CS> as PCS<F>>::CK,
        fs: &[Poly<F>],
        xss: &[HashSet<F>],
        transcript: &mut T,
    ) -> (CS::G, CS::Proof)
    {
        assert_eq!(xss.len(), fs.len(), "{} opening sets specified for {} polynomials", xss.len(), fs.len());
        let mut opening_set = HashSet::new();
        for xs in xss {
            opening_set.extend(xs);
        }
        let z = crate::utils::z_of_set(&opening_set);

        let zs: Vec<_> = xss.iter()
            .map(|xs| crate::utils::z_of_set(xs))
            .collect();

        let (qs, rs): (Vec<_>, Vec<_>) = fs.iter().zip(&zs)
            .map(|(fi, zi)| fi.divide_with_q_and_r(zi))
            .unzip();

        let gamma = transcript.get_gamma();
        let q = crate::utils::randomize(gamma, &qs);
        let q_comm = CS::commit(ck, &q);//scheme.commit(&q);
        transcript.commit_to_q(&q_comm);
        let zeta = transcript.get_zeta();

        let z_zeta = z.evaluate(&zeta);
        let mut zs_zeta: Vec<_> = zs.iter().map(|zi| zi.evaluate(&zeta)).collect();
        let rs_zeta: Vec<_> = rs.iter().map(|ri| ri.evaluate(&zeta)).collect();
        ark_ff::batch_inversion(&mut zs_zeta);

        let gs = crate::utils::powers(gamma, fs.len() - 1);

        let mut l = DensePolynomial::zero();
        for (((fi, ri), zi_inv), gi) in fs.iter()
            .zip(rs_zeta)
            .zip(zs_zeta)
            .zip(gs) {
            l += (gi * zi_inv, &(fi - &DensePolynomial::from_coefficients_vec(vec![ri])));
        }
        let l = &(&l - &q) * z_zeta;

        // let z_of_zeta = crate::utils::z_of_point(&zeta);
        // let (q_of_l, r_of_l) = l.divide_with_q_and_r(&z_of_zeta);
        // assert!(r_of_l.is_zero());
        // let q_of_l1 = CS::commit(ck, &q_of_l);//scheme.commit(&q_of_l);
        let q_of_l1 = CS::open(ck, &l, &zeta);
        (q_comm, q_of_l1)
    }

    fn get_linearization_commitment<T: ShplonkTranscript<F, CS::G>>(
        fcs: &[<Shplonk<F, CS> as PCS<F>>::G],
        qc: &<Shplonk<F, CS> as PCS<F>>::G,
        onec: &<Shplonk<F, CS> as PCS<F>>::G,
        xss: &Vec<Vec<F>>,
        yss: &Vec<Vec<F>>,
        transcript: &mut T,
    ) -> <Shplonk<F, CS> as PCS<F>>::G
    {
        let gamma = transcript.get_gamma();
        transcript.commit_to_q(&qc);
        let zeta = transcript.get_zeta();

        let mut opening_set = HashSet::new();
        for xs in xss {
            opening_set.extend(xs);
        }

        let z_at_zeta = opening_set.iter()
            .map(|xi| zeta - xi)
            .reduce(|z, zi| z * zi)
            .unwrap();

        let rs = xss.iter().zip(yss).map(|(xs, ys)| Self::interpolate(&xs, &ys));
        let rs_at_zeta = rs.map(|ri| ri.evaluate(&zeta));

        let mut zs_at_zeta: Vec<_> = xss.iter().map(|xs|
            xs.iter()
                .map(|xi| zeta - xi)
                .reduce(|z, zi| z * zi)
                .unwrap()
        ).collect();

        ark_ff::batch_inversion(&mut zs_at_zeta);

        let gs = crate::utils::powers(gamma, fcs.len() - 1);

        let gzs: Vec<_> = gs.iter().zip(zs_at_zeta).map(|(&gi, zi_inv)| gi * zi_inv).collect();

        let fc: <Shplonk<F, CS> as PCS<F>>::G = fcs.iter().zip(&gzs)
            .map(|(f1, gzi)| f1.mul(z_at_zeta * gzi))
            .sum();

        let r = rs_at_zeta.zip(gzs).map(|(ri_at_zeta, gzi)| ri_at_zeta * &gzi).sum::<F>() * z_at_zeta;

        fc - onec.mul(r) - qc.mul(z_at_zeta)
    }

    pub fn verify_many<T: ShplonkTranscript<F, CS::G>>(
        vk: &<Shplonk<F, CS> as PCS<F>>::VK,
        fcs: &[<Shplonk<F, CS> as PCS<F>>::G],
        proof: &(CS::G, CS::Proof),
        xss: &Vec<Vec<F>>,
        yss: &Vec<Vec<F>>,
        transcript: &mut T,
    ) -> bool
    {
        let qc = &proof.0; // commitment to the original quotient
        let lqc = &proof.1; // commitment to the quotient of the linearization polynomial
        let onec = vk.commit_to_one();
        let lc = Self::get_linearization_commitment(fcs, qc, &onec, xss, yss, transcript);
        CS::verify(vk, &lc, &transcript.get_zeta(), &F::zero(), lqc)//verify(&lc, &transcript.get_zeta(), F::zero(), &lqc)
    }

    fn interpolate(xs: &[F], ys: &[F]) -> DensePolynomial<F> {
        let x1 = xs[0];
        let mut l = crate::utils::z_of_point(&x1);
        for &xj in xs.iter().skip(1) {
            let q = crate::utils::z_of_point(&xj);
            l = &l * &q;
        }

        let mut ws = vec![];
        for xj in xs {
            let mut wj = F::one();
            for xk in xs {
                if xk != xj {
                    let d = *xj - xk;
                    wj *= d;
                }
            }
            ws.push(wj);
        }
        ark_ff::batch_inversion(&mut ws);

        let mut res = DensePolynomial::zero();
        for ((&wi, &xi), &yi) in ws.iter().zip(xs).zip(ys) {
            let d = crate::utils::z_of_point(&xi);
            let mut z = &l / &d;
            z = &z * wi;
            z = &z * yi;
            res = res + z;
        }
        res
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::F;
    use crate::pcs::tests::IdentityCommitment;

    use ark_std::{test_rng, UniformRand};
    use ark_std::iter::FromIterator;
    use ark_std::rand::Rng;
    use crate::Poly;

    type TestShplonk = Shplonk<F, IdentityCommitment>;

    impl<G> ShplonkTranscript<F, G> for (F, F) {
        fn get_gamma(&mut self) -> F { self.0 }

        fn commit_to_q(&mut self, _q_comm: &G) {}

        fn get_zeta(&mut self) -> F { self.1 }
    }

    pub fn generate_test_data<R, F, CS>(
        rng: &mut R,
        ck: &<Shplonk<F, CS> as PCS<F>>::CK,
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

    #[test]
    fn test_shplonk() {
        let rng = &mut test_rng();

        let (ck, vk) = TestShplonk::setup();

        let t = 4; // number of polynomials
        let max_m = 3; // maximal number of opening points per polynomial

        let xss: Vec<_> = (0..t)
            .map(|_| (0..rng.gen_range(1..max_m))
                .map(|_| F::rand(rng)).collect::<Vec<_>>())
            .collect();

        let (fs, fcs, yss) =
            generate_test_data::<_, _, IdentityCommitment>(rng, &ck,15, 4, &xss);

        let sets_of_xss: Vec<_> = xss.iter()
            .map(|xs| HashSet::from_iter(xs.iter().cloned()))
            .collect();

        let transcript = &mut (F::rand(rng), F::rand(rng));

        let (qc, qlc) = TestShplonk::open_many(&ck, &fs, sets_of_xss.as_slice(), transcript);

        let onec = ().commit_to_one();
        let lc = TestShplonk::get_linearization_commitment(&fcs, &qc, &onec, &xss, &yss, transcript);
        assert!(lc.0.evaluate(&transcript.1).is_zero());

        assert!(TestShplonk::verify_many(&vk, &fcs, &(qc, qlc), &xss, &yss, transcript))
    }
}