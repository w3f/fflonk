use std::collections::HashSet;

use ark_ff::{PrimeField, Zero};
use ark_poly::{Polynomial, UVPolynomial};

use crate::{EuclideanPolynomial, Poly};
use crate::pcs::{Commitment, PCS};

pub struct MultipointClaim<F: PrimeField, C: Commitment<F>> {
    pub c: C,
    pub xs: Vec<F>,
    pub ys: Vec<F>,
}

pub trait Transcript<F, G> {
    fn get_gamma(&mut self) -> F;
    fn commit_to_q(&mut self, q_comm: &G);
    fn get_zeta(&mut self) -> F;
}


pub fn aggregate_polys<F: PrimeField, CS: PCS<F>, T: Transcript<F, CS::C>>(
    ck: &CS::CK,
    fs: &[Poly<F>],
    xss: &[HashSet<F>],
    transcript: &mut T,
) -> (Poly<F>, F, CS::C) {
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

pub fn group_by_commitment<F: PrimeField, C: Commitment<F>>(
    fcs: &[C],
    xss: &Vec<Vec<F>>,
    yss: &Vec<Vec<F>>,
) -> Vec<MultipointClaim<F, C>> {
    fcs.iter().cloned()
        .zip(xss.iter().cloned())
        .zip(yss.iter().cloned())
        .map(|((c, xs), ys)| MultipointClaim { c, xs, ys })
        .collect()
}

pub fn aggregate_claims<F: PrimeField, CS: PCS<F>, T: Transcript<F, CS::C>>(
    claims: Vec<MultipointClaim<F, CS::C>>,
    qc: &CS::C,
    onec: &CS::C,
    transcript: &mut T,
) -> MultipointClaim<F, CS::C>
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
        .map(|MultipointClaim { c: _, xs, ys }| interpolate(xs, ys));
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

    let fc: CS::C = claims.iter().zip(gzs.iter())
        .map(|(claim, &gzi)| claim.c.mul(z_at_zeta * gzi))
        .sum();

    let r = rs_at_zeta.zip(gzs).map(|(ri_at_zeta, gzi)| ri_at_zeta * &gzi).sum::<F>() * z_at_zeta;

    let c = fc - onec.mul(r) - qc.mul(z_at_zeta);
    MultipointClaim { c, xs: vec![zeta], ys: vec![F::zero()] }
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


#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;
    use crate::pcs::tests::IdentityCommitment;
    use crate::shplonk::tests::{random_xss, random_opening};
    use ark_bw6_761::{Fr, BW6_761};
    use crate::pcs::kzg::KZG;
    use crate::pcs::PcsParams;
    use ark_std::iter::FromIterator;


    impl<F: PrimeField, G> Transcript<F, G> for (F, F) {
        fn get_gamma(&mut self) -> F { self.0 }

        fn commit_to_q(&mut self, _q_comm: &G) {}

        fn get_zeta(&mut self) -> F { self.1 }
    }


    fn _test_aggregation<F: PrimeField, CS: PCS<F>>() {
        let rng = &mut test_rng();

        let d = 15; // degree of polynomials
        let t = 4; // number of polynomials
        let max_m = 3; // maximal number of opening points per polynomial

        let params = CS::setup(d, rng);
        let (ck, vk) = (params.ck(), params.vk());

        let xss = random_xss(rng, t, max_m);
        let opening = random_opening::<_, _, CS>(rng, &ck, d, t, xss);

        let sets_of_xss: Vec<HashSet<F>> = opening.xss.iter()
            .map(|xs| HashSet::from_iter(xs.iter().cloned()))
            .collect();

        let transcript = &mut (F::rand(rng), F::rand(rng));

        let (agg_poly, zeta, agg_proof) = aggregate_polys::<_, CS, _>(&ck, &opening.fs, &sets_of_xss, transcript);
        let claims = group_by_commitment(&opening.fcs, &opening.xss, &opening.yss);
        let onec = CS::commit(&vk.clone().into(), &Poly::from_coefficients_slice(&[F::one()]));
        let agg_claim = aggregate_claims::<_, CS, _>(claims, &agg_proof, &onec, transcript);

        assert_eq!(CS::commit(&ck, &agg_poly), agg_claim.c);
        assert_eq!(zeta, agg_claim.xs[0]);
        assert_eq!(agg_poly.evaluate(&zeta), agg_claim.ys[0]);
        assert!(agg_claim.ys[0].is_zero());
    }

    #[test]
    fn test_aggregation_id() {
        _test_aggregation::<Fr, IdentityCommitment>();
    }

    #[test]
    fn test_aggregation_kzg() {
        _test_aggregation::<Fr, KZG<BW6_761>>();
    }
}