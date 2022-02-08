use std::collections::HashSet;

use ark_ff::{PrimeField, Zero};
use ark_poly::{Polynomial, UVPolynomial};

use crate::{EuclideanPolynomial, Poly};
use crate::pcs::{Commitment, PCS};
use ark_std::iterable::Iterable;
use crate::utils::poly::interpolate_evaluate;

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

    // Vanishing polynomial of the whole opening set evaluated at zeta
    let agg_z_at_zeta = opening_set.iter()
        .map(|x| zeta - x)
        .reduce(|a, b| a * b)
        .expect("TODO");

    // For each polynomial the opening claim {(xi, yi)} can be presented in polynomial form
    // as a pair of polynomials (r, z), where z is the vanishing polynomial of the set {xi},
    // and r is the interpolation polynomial of the set {(xi, yi)}.
    // rj(zeta), zj(zeta)
    let (rs_at_zeta, zs_at_zeta): (Vec<_>, Vec<_>) = claims.iter()
        .map(|MultipointClaim { c: _, xs, ys }| interpolate_evaluate(xs, ys, &zeta))
        .unzip();

    let mut zs_at_zeta_inv = zs_at_zeta;
    ark_ff::batch_inversion(&mut zs_at_zeta_inv);

    // 1, gamma, ..., gamma^{k-1}
    let powers = crate::utils::powers(gamma, claims.len() - 1);
    let coeffs: Vec<F> = powers.iter().zip(zs_at_zeta_inv.iter())
        .map(|(&gj, zj_inv)| gj * zj_inv * agg_z_at_zeta) // (g^j / z_j) * agg_z
        .collect();

    //TODO: multiexp
    let agg_c: CS::C = claims.iter().zip(coeffs.iter())
        .map(|(claim, &coeff)| claim.c.mul(coeff))
        .sum();

    let agg_r_at_zeta: F = rs_at_zeta.into_iter().zip(coeffs.iter())
        .map(|(ri_at_zeta, coeff)| ri_at_zeta * coeff)
        .sum();

    let c = agg_c - onec.mul(agg_r_at_zeta) - qc.mul(agg_z_at_zeta);
    MultipointClaim { c, xs: vec![zeta], ys: vec![F::zero()] }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;
    use crate::pcs::tests::IdentityCommitment;
    use crate::shplonk::tests::{random_xss, random_opening};
    use crate::pcs::PcsParams;
    use ark_std::iter::FromIterator;
    use crate::tests::{TestKzg, TestField};


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
        _test_aggregation::<TestField, IdentityCommitment>();
    }

    #[test]
    fn test_aggregation_kzg() {
        _test_aggregation::<TestField, TestKzg>();
    }
}