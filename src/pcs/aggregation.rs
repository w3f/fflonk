use crate::shplonk::ShplonkTranscript;
use crate::{Poly, EuclideanPolynomial};
use std::collections::HashSet;
use ark_ff::{PrimeField, Zero};
use crate::pcs::{PCS, CommitmentSpace};
use ark_poly::{Polynomial, UVPolynomial};

pub fn aggregate_polys<F: PrimeField, CS: PCS<F>, T: ShplonkTranscript<F, CS::G>>(
    ck: &CS::CK,
    fs: &[Poly<F>],
    xss: &[HashSet<F>],
    transcript: &mut T,
) -> (Poly<F>, F, CS::G) {
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

    let mut l = Poly::zero();
    for (((fi, ri), zi_inv), gi) in fs.iter()
        .zip(rs_zeta)
        .zip(zs_zeta)
        .zip(gs) {
        l += (gi * zi_inv, &(fi - &Poly::from_coefficients_vec(vec![ri])));
    }
    let l = &(&l - &q) * z_zeta;
    (l, zeta, q_comm)
}

pub struct Claim<F: PrimeField, C: CommitmentSpace<F>> {
    pub c: C,
    pub xs: Vec<F>,
    pub ys: Vec<F>,
}


pub fn aggregate_claims<F: PrimeField, CS: PCS<F>, T: ShplonkTranscript<F, CS::G>>(
    claims: Vec<Claim<F, CS::G>>,
    qc: &CS::G,
    onec: &CS::G,
    transcript: &mut T,
) -> Claim<F, CS::G>
{
    let gamma = transcript.get_gamma();
    transcript.commit_to_q(&qc);
    let zeta = transcript.get_zeta();

    let mut opening_set = HashSet::new();
    for xs in claims.iter().map(|claim| &claim.xs) {
        opening_set.extend(xs);
    }

    let z_at_zeta = opening_set.iter()
        .map(|xi| zeta - xi)
        .reduce(|z, zi| z * zi)
        .unwrap();

    let rs = claims.iter()
        .map(|Claim { c: _, xs, ys }| interpolate(xs, ys));
    let rs_at_zeta = rs.map(|ri| ri.evaluate(&zeta));

    let mut zs_at_zeta: Vec<F> = claims.iter().map(|claim|
        claim.xs.iter()
            .map(|xi| zeta - xi)
            .reduce(|z, zi| z * zi)
            .unwrap()
    ).collect();

    ark_ff::batch_inversion(&mut zs_at_zeta);

    let gs = crate::utils::powers(gamma, claims.len() - 1);
    assert_eq!(gs.len(), zs_at_zeta.len());
    let gzs: Vec<F> = gs.iter().zip(zs_at_zeta.iter()).map(|(&gi, zi_inv)| gi * zi_inv).collect();
    assert_eq!(claims.len(), gzs.len());

    let fc: CS::G = claims.iter().zip(gzs.iter())
        .map(|(claim, &gzi)| claim.c.mul(z_at_zeta * gzi))
        .sum();

    let r = rs_at_zeta.zip(gzs).map(|(ri_at_zeta, gzi)| ri_at_zeta * &gzi).sum::<F>() * z_at_zeta;

    let c = fc - onec.mul(r) - qc.mul(z_at_zeta);
    Claim { c, xs: vec![zeta], ys: vec![F::zero()] }
}


fn interpolate<F: PrimeField>(xs: &[F], ys: &[F]) -> Poly<F> {
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

    let mut res = Poly::zero();
    for ((&wi, &xi), &yi) in ws.iter().zip(xs).zip(ys) {
        let d = crate::utils::z_of_point(&xi);
        let mut z = &l / &d;
        z = &z * wi;
        z = &z * yi;
        res = res + z;
    }
    res
}