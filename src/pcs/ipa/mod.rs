use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::{batch_inversion, Field, One};
use ark_std::test_rng;
use ark_std::vec::Vec;
use ark_std::rand::Rng;
use ark_std::vec;

fn scalar_prod<F: Field>(a: &[F], b: &[F]) -> F {
    ark_std::cfg_iter!(a)
        .zip(b)
        .map(|(a, b)| *a * b)
        .sum()
}

// Computes `l + xr` pointwise.
fn fold_points<C: AffineRepr>(l: &[C], r: &[C], x: &C::ScalarField) -> Vec<C> {
    assert_eq!(l.len(), r.len());
    let proj: Vec<C::Group> = ark_std::cfg_iter!(l)
        .zip(r)
        .map(|(&l, &r)| (l + r * x))
        .collect();
    C::Group::normalize_batch(&proj)
}

// Computes `l + xr` pointwise.
fn fold_scalars<F: Field>(l: &[F], r: &[F], x: &F) -> Vec<F> {
    assert_eq!(l.len(), r.len());
    ark_std::cfg_iter!(l)
        .zip(r)
        .map(|(&l, &r)| (l + r * x))
        .collect()
}


// n = 2^m
// Folding elements V = [A1, ..., An] with scalars [x1, ..., xm] recursively m times using formula
// V = VL || VR, Vi = VL + xi * VR pointwise, V := Vi, i = 1,...,m
// results in V = Vm = [c1A1 + ... + cnAn], where ci = prod({xj | if j-th bit of i-1 is set}).
// This function computes these ci-s.
fn final_folding_exponents<F: Field>(xs: &[F]) -> Vec<F> {
    let m = xs.len();
    let mut n = 2usize.pow(m as u32);
    let mut res = vec![F::one(); n];
    for x in xs {
        n = n / 2;
        for chunk in res.rchunks_mut(n).step_by(2) {
            for elem in chunk.iter_mut() {
                *elem *= x;
            }
        }
    }
    res
}

pub struct Proof<C: AffineRepr> {
    ls: Vec<C>,
    rs: Vec<C>,
    final_p: C,
    final_a: C::ScalarField,
    final_b: C::ScalarField,
    xs: Vec<C::ScalarField>,
}

pub fn bullet_prove<C: AffineRepr>(log_n: usize, g: &[C], h: &[C], a: &[C::ScalarField], b: &[C::ScalarField], p: C, u: C) -> Proof<C> {
    let n = 2usize.pow(log_n as u32);
    assert_eq!(g.len(), n);
    assert_eq!(h.len(), n);
    assert_eq!(a.len(), n);
    assert_eq!(b.len(), n);


    let mut g_folded = g.to_vec();
    let mut h_folded = h.to_vec();
    let mut a_folded = a.to_vec();
    let mut b_folded = b.to_vec();

    let mut ls = Vec::<C>::with_capacity(log_n);
    let mut rs = Vec::<C>::with_capacity(log_n);

    let xs: Vec<_> = (0..log_n).map(|_| C::ScalarField::from(test_rng().gen::<u128>())).collect();

    let mut p1 = p;

    let mut n1 = n;
    for x in xs.iter() {
        n1 /= 2;

        let gl = &g_folded[..n1];
        let gr = &g_folded[n1..];
        let hl = &h_folded[..n1];
        let hr = &h_folded[n1..];
        let al = &a_folded[..n1];
        let ar = &a_folded[n1..];
        let bl = &b_folded[..n1];
        let br = &b_folded[n1..];

        let cl = scalar_prod(al, br);
        let cr = scalar_prod(ar, bl);

        let points = [gr, hl, &[u]].concat();
        let scalars = [al, br, &[cl]].concat();
        let l = C::Group::msm(&points, &scalars).unwrap();
        ls.push(l.into());

        let points = [gl, hr, &[u]].concat();
        let scalars = [ar, bl, &[cr]].concat();
        let r = C::Group::msm(&points, &scalars).unwrap();
        rs.push(r.into());

        let x_inv = x.inverse().unwrap();

        g_folded = fold_points(gl, gr, &x_inv);
        h_folded = fold_points(hl, hr, &x);
        a_folded = fold_scalars(al, ar, &x);
        b_folded = fold_scalars(bl, br, &x_inv);

        p1 = ((l * x_inv) + (r * x) + p1).into_affine();
    }

    let final_a = a_folded[0];
    let final_b = b_folded[0];

    Proof {
        ls,
        rs,
        final_p: p1,
        final_a,
        final_b,
        xs,
    }
}

pub fn verify<C: AffineRepr>(g: &[C], h: &[C], p :C, u: C, proof: Proof<C>) {
    let xs = proof.xs;
    let mut xs_inv = xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let g_exps = final_folding_exponents(&xs_inv);
    let h_exps = final_folding_exponents(&xs);

    let g_exps = g_exps.iter().map(|e| proof.final_a * e).collect();
    let h_exps = h_exps.iter().map(|e| proof.final_b * e).collect();


    let final_c = proof.final_a * proof.final_b;

    let points = [g, h, &[u]].concat();
    let scalars = [g_exps, h_exps, vec![final_c]].concat();
    let res1 = C::Group::msm(&points, &scalars).unwrap();

    let points = [&[p], proof.ls.as_slice(), proof.rs.as_slice()].concat();
    let scalars = [vec![C::ScalarField::one()], xs_inv, xs].concat();
    let res2 = C::Group::msm(&points, &scalars).unwrap();
    assert_eq!(res1, res2);
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::{Fr, G1Affine, G1Projective};
    use ark_ff::vec;
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn bullet() {
        let rng = &mut test_rng();

        let log_n = 2;
        let n = 2usize.pow(log_n as u32);

        let g = (0..n).map(|_| G1Affine::rand(rng)).collect::<Vec<_>>();
        let h = (0..n).map(|_| G1Affine::rand(rng)).collect::<Vec<_>>();
        let a = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let b = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let c = scalar_prod(&a, &b);
        let u = G1Affine::rand(rng);

        let points = [g.clone(), h.clone(), vec![u]].concat();
        let scalars = [a.clone(), b.clone(), vec![c]].concat();
        let p: G1Affine = G1Projective::msm(&points, &scalars).unwrap().into_affine();

        let proof = bullet_prove(log_n, &g, &h, &a, &b, p, u);
        verify(&g, &h, p, u, proof);
    }
}