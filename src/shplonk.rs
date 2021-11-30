use ark_std::Zero;
use ark_ff::{Field, FftField, PrimeField};
use ark_poly::{UVPolynomial, Polynomial};
use ark_poly::univariate::{DensePolynomial, DenseOrSparsePolynomial};

use std::collections::HashSet;

use crate::{AdditiveCommitment, CommitmentScheme};


pub trait ShplonkTranscript<F, G> {
    fn get_gamma(&mut self) -> F;
    fn commit_to_q(&mut self, q_comm: &G);
    fn get_zeta(&mut self) -> F;
}


trait EuclideanPolynomial: Sized {
    fn divide_with_q_and_r(&self, divisor: &Self) -> (Self, Self);
}

impl<F: Field> EuclideanPolynomial for DensePolynomial<F> {
    fn divide_with_q_and_r(&self, divisor: &Self) -> (Self, Self) {
        let a: DenseOrSparsePolynomial<F> = self.into();
        let b: DenseOrSparsePolynomial<F> = divisor.into();
        a.divide_with_q_and_r(&b).unwrap()
    }
}


pub fn open<F, C, T>(
    fs: &[DensePolynomial<F>],
    xss: &[HashSet<F>],
    scheme: &C,
    transcript: &mut T,
) -> (C::G, C::G)
    where
        F: PrimeField, //TODO
        C: CommitmentScheme<F, DensePolynomial<F>>,
        T: ShplonkTranscript<F, C::G>,
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
    let q_comm = scheme.commit(&q);
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

    let z_of_zeta = crate::utils::z_of_point(&zeta);
    let (q_of_l, r_of_l) = l.divide_with_q_and_r(&z_of_zeta);
    assert!(r_of_l.is_zero());
    let q_of_l1 = scheme.commit(&q_of_l);
    (q_comm, q_of_l1)
}

fn get_linearization_commitment<F, C, T>(
    fcs: &[C::G],
    qc: &C::G,
    xss: &Vec<Vec<F>>,
    yss: &Vec<Vec<F>>,
    scheme: &C,
    transcript: &mut T,
) -> C::G
    where
        F: PrimeField,
        C: CommitmentScheme<F, DensePolynomial<F>>,
        T: ShplonkTranscript<F, C::G>,
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

    let rs = xss.iter().zip(yss).map(|(xs, ys)| interpolate(&xs, &ys));
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

    let fc: C::G = fcs.iter().zip(&gzs)
        .map(|(f1, gzi)| f1.mul(z_at_zeta * gzi))
        .sum();

    let r = rs_at_zeta.zip(gzs).map(|(ri_at_zeta, gzi)| ri_at_zeta * &gzi).sum::<F>() * z_at_zeta;

    fc - scheme.commit_const(&r) - qc.mul(z_at_zeta)
}

pub fn verify<F, C, T>(
    fcs: &[C::G],
    proof: (C::G, C::G),
    xss: &Vec<Vec<F>>,
    yss: &Vec<Vec<F>>,
    scheme: &C,
    transcript: &mut T,
) -> bool
    where
        F: PrimeField,
        C: CommitmentScheme<F, DensePolynomial<F>>,
        T: ShplonkTranscript<F, C::G>,
{
    let qc = proof.0; // commitment to the original quotient
    let lqc = proof.1; // commitment to the quotient of the linearization polynomial
    let lc = get_linearization_commitment(fcs, &qc, xss, yss, scheme, transcript);
    scheme.verify(&lc, &transcript.get_zeta(), F::zero(), &lqc)
}

fn interpolate<F: FftField>(xs: &[F], ys: &[F]) -> DensePolynomial<F> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::{IdentityCommitment, F};

    use ark_std::{test_rng, UniformRand};
    use ark_std::iter::FromIterator;
    use ark_std::rand::Rng;
    use crate::Poly;

    impl<G> ShplonkTranscript<F, G> for (F, F) {
        fn get_gamma(&mut self) -> F { self.0 }

        fn commit_to_q(&mut self, _q_comm: &G) {}

        fn get_zeta(&mut self) -> F { self.1 }
    }

    pub fn generate_test_data<R, F, C>(
        rng: &mut R,
        d: usize, // degree of polynomials
        t: usize, // number of polynomials
        xss: &Vec<Vec<F>>, // vecs of opening points per polynomial
        scheme: &C, // commitment scheme
    ) -> (
        Vec<Poly<F>>, // polynomials
        Vec<C::G>, // commitments
        Vec<Vec<F>>, // evaluations per polynomial
    ) where
        R: Rng,
        F: PrimeField,
        C: CommitmentScheme<F, Poly<F>>
    {
        // polynomials
        let fs: Vec<Poly<F>> = (0..t)
            .map(|_| Poly::rand(d, rng))
            .collect();
        // commitments
        let fcs: Vec<_> = fs.iter()
            .map(|fi| scheme.commit(fi))
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
        let scheme = IdentityCommitment {};

        let t = 4; // number of polynomials
        let max_m = 3; // maximal number of opening points per polynomial

        let xss: Vec<_> = (0..t)
            .map(|_| (0..rng.gen_range(1..max_m))
                .map(|_| F::rand(rng)).collect::<Vec<_>>())
            .collect();

        let (fs, fcs, yss) =
            generate_test_data(rng, 15, 4, &xss, &scheme);

        let sets_of_xss: Vec<_> = xss.iter()
            .map(|xs| HashSet::from_iter(xs.iter().cloned()))
            .collect();

        let transcript = &mut (F::rand(rng), F::rand(rng));

        let (qc, qlc) = open(&fs, sets_of_xss.as_slice(), &scheme, transcript);

        let lc = get_linearization_commitment(&fcs, &qc, &xss, &yss, &scheme, transcript);
        assert!(lc.0.evaluate(&transcript.1).is_zero());

        assert!(verify(&fcs, (qc, qlc), &xss, &yss, &scheme, transcript))
    }
}