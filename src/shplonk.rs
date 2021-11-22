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
    let z = crate::z_of_set(&opening_set);

    let zs: Vec<_> = xss.iter()
        .map(|xs| crate::z_of_set(xs))
        .collect();

    let (qs, rs): (Vec<_>, Vec<_>) = fs.iter().zip(&zs)
        .map(|(fi, zi)| fi.divide_with_q_and_r(zi))
        .unzip();

    let gamma = transcript.get_gamma();
    let q = crate::randomize(gamma, &qs);
    let q_comm = scheme.commit(&q);
    transcript.commit_to_q(&q_comm);
    let zeta = transcript.get_zeta();

    let z_zeta = z.evaluate(&zeta);
    let mut zs_zeta: Vec<_> = zs.iter().map(|zi| zi.evaluate(&zeta)).collect();
    let rs_zeta: Vec<_> = rs.iter().map(|ri| ri.evaluate(&zeta)).collect();
    ark_ff::batch_inversion(&mut zs_zeta);

    let gs = crate::powers(gamma, fs.len() - 1);

    let mut l = DensePolynomial::zero();
    for (((fi, ri), zi_inv), gi) in fs.iter()
        .zip(rs_zeta)
        .zip(zs_zeta)
        .zip(gs) {
        l += (gi * zi_inv, &(fi - &DensePolynomial::from_coefficients_vec(vec![ri])));
    }
    let l = &(&l - &q) * z_zeta;

    let z_of_zeta = crate::z_of_point(&zeta);
    let (q_of_l, r_of_l) = l.divide_with_q_and_r(&z_of_zeta);
    assert!(r_of_l.is_zero());
    let q_of_l1 = scheme.commit(&q_of_l);
    (q_comm, q_of_l1)
}

pub fn verify<F, C, T>(
    f1s: &[C::G],
    q1: C::G,
    q_of_l1: C::G,
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
    transcript.commit_to_q(&q1);
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

    let gs = crate::powers(gamma, f1s.len() - 1);

    let gzs: Vec<_> = gs.iter().zip(zs_at_zeta).map(|(&gi, zi_inv)| gi * zi_inv).collect();

    let f_1: C::G = f1s.iter().zip(&gzs)
        .map(|(f1, gzi)| f1.mul(z_at_zeta * gzi))
        .sum();
        // .reduce(|a, b| a + b) //TODO: Sum
        // .unwrap();

    let r = rs_at_zeta.zip(gzs).map(|(ri_at_zeta, gzi)| ri_at_zeta * &gzi).sum::<F>() * z_at_zeta;

    let l_1 = f_1 - scheme.commit_const(&r) - q1.mul(z_at_zeta);

    // E::product_of_pairings(&[
    //     (l_1.into_affine().into(), pvk.prepared_h.clone()),
    //     ((-q_of_l1).into(), pvk.prepared_beta_h.clone()),
    // ]).is_one()
    l_1
}

fn interpolate<F: FftField>(xs: &[F], ys: &[F]) -> DensePolynomial<F> {
    let x1 = xs[0];
    let mut l = crate::z_of_point(&x1);
    for &xj in xs.iter().skip(1) {
        let q = crate::z_of_point(&xj);
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
        let d = crate::z_of_point(&xi);
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
    use crate::tests::IdentityCommitment;

    use ark_std::{test_rng, UniformRand};
    use ark_std::rand::{Rng, rngs::StdRng};
    use ark_std::iter::FromIterator;

    type F = ark_bw6_761::Fr;
    type P = DensePolynomial<F>; //TODO: to lib

    impl<G> ShplonkTranscript<F, G> for (F, F) {
        fn get_gamma(&mut self) -> F { self.0 }

        fn commit_to_q(&mut self, q_comm: &G) {}

        fn get_zeta(&mut self) -> F { self.1 }
    }

    #[test]
    fn test_shplonk() {
        let rng = &mut test_rng();
        let scheme = IdentityCommitment {};

        let d = 15; // degree of polynomials
        let t = 4; // number of polynomials
        let max_m = 3; // maximal number of opening points per polynomial
        let ms: Vec<usize> =  (0..t) // number of opening points per polynomial
            .map(|_| rng.gen_range(1..max_m))
            .collect();

        // polynomials
        let fs: Vec<P> = (0..t)
            .map(|_| P::rand(d, rng))
            .collect();
        // commitments
        let fcs: Vec<_> = fs.iter().map(|fi| scheme.commit(fi)).collect();
        // opening points per polynomial
        let xss: Vec<_> = ms.into_iter().map(|m| (0..m).map(|_| F::rand(rng)).collect::<Vec<_>>()).collect();
        let sets_of_xss: Vec<_> = xss.iter().map(|xs| HashSet::from_iter(xs.iter().cloned())).collect();
        // evaluations
        let yss: Vec<_> = fs.iter().zip(xss.iter()).map(|(f, xs)| {
            xs.iter().map(|x| f.evaluate(x)).collect::<Vec<_>>()
        }).collect();

        let transcript = &mut (F::rand(rng), F::rand(rng));

        let (qc, qlc) = open(&fs, sets_of_xss.as_slice(), &scheme, transcript);

        let lc = verify(&fcs, qc, qlc, &xss, &yss, &scheme, transcript);

        assert!(lc.0.evaluate(&transcript.1).is_zero());
    }
}