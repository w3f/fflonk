use crate::pcs::ipa::{evaluate_final_poly, final_folding_exponents, fold_points, fold_scalars, scalar_prod};
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::{batch_inversion, Field, One};
use ark_std::rand::Rng;
use ark_std::test_rng;
use ark_std::vec::Vec;

pub struct Proof<C: AffineRepr> {
    l: Vec<C>,
    r: Vec<C>,
    g_final: C, // aka U
    f_final: C::ScalarField,
    xs: Vec<C::ScalarField>, //TODO:
}

pub fn open<C: AffineRepr>(log_n: usize, g: Vec<C>, h: C, f: Vec<C::ScalarField>, z: Vec<C::ScalarField>) -> Proof<C> {
    let zz = z[1];

    let n = 2usize.pow(log_n as u32);
    assert_eq!(g.len(), n);
    assert_eq!(f.len(), n);
    assert_eq!(z.len(), n);

    let mut g_folded = g;
    let mut f_folded = f;
    let mut z_folded = z;

    let mut l = Vec::<C>::with_capacity(log_n);
    let mut r = Vec::<C>::with_capacity(log_n);

    let xs: Vec<_> = (0..log_n).map(|_| C::ScalarField::from(test_rng().gen::<u128>())).collect();

    let mut n1 = n;
    for x in xs.iter() {
        n1 /= 2;

        let gl = &g_folded[..n1];
        let gr = &g_folded[n1..];
        let fl = &f_folded[..n1];
        let fr = &f_folded[n1..];
        let zl = &z_folded[..n1];
        let zr = &z_folded[n1..];

        let cl = scalar_prod(fl, zr);
        let cr = scalar_prod(fr, zl);

        let points = [gr, &[h]].concat();
        let scalars = [fl, &[cl]].concat();
        l.push(C::Group::msm(&points, &scalars).unwrap().into());

        let points = [gl, &[h]].concat();
        let scalars = [fr, &[cr]].concat();
        r.push(C::Group::msm(&points, &scalars).unwrap().into());

        let x_inv = x.inverse().unwrap();

        g_folded = fold_points(gl, gr, &x_inv);
        f_folded = fold_scalars(fl, fr, &x);
        z_folded = fold_scalars(zl, zr, &x_inv);
    }

    let g_final = g_folded[0];
    let f_final = f_folded[0];

    Proof {
        l,
        r,
        g_final,
        f_final,
        xs,
    }
}

// f(z) = v
pub fn check<C: AffineRepr>(g: Vec<C>, h: C, c: C, z: C::ScalarField, v: C::ScalarField, proof: Proof<C>) {
    let xs = proof.xs;
    let mut xs_inv = xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let final_exps = final_folding_exponents(&xs_inv);
    let z_final = evaluate_final_poly(&xs_inv, &z);

    let bases = [proof.l, proof.r].concat();
    let exps = [xs_inv, xs].concat();
    let res1 = C::Group::msm(&bases, &exps).unwrap();
    let res1  = res1 + c + (h * v);

    let res2 = proof.g_final * proof.f_final + h * (proof.f_final * z_final);
    assert_eq!(res1, res2);

    assert_eq!(proof.g_final, C::Group::msm(&g, &final_exps).unwrap().into_affine());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils;
    use ark_bls12_381::{Fr, G1Affine, G1Projective};
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn ipa_pc() {
        let rng = &mut test_rng();

        let log_n = 2;
        let n = 2usize.pow(log_n as u32);

        let g = (0..n).map(|_| G1Affine::rand(rng)).collect::<Vec<_>>();
        let h = G1Affine::rand(rng);

        let f = DensePolynomial::<Fr>::rand(n - 1, rng);
        let f_comm: G1Affine = G1Projective::msm(&g, &f.coeffs).unwrap().into_affine();

        let z = Fr::rand(rng);
        let z_powers: Vec<Fr> = utils::powers(z).take(n).collect();
        let v = f.evaluate(&z);

        let proof = open(log_n, g.clone(), h, f.coeffs, z_powers);
        check(g, h, f_comm, z, v, proof);
    }
}