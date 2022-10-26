use std::ops::Mul;
use criterion::{criterion_group, criterion_main, Criterion};

use ark_ff::UniformRand;
use ark_std::test_rng;

use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ec::pairing::Pairing;
use ark_bw6_761::{BW6_761};
use fflonk::utils::curve_name;


fn scalar_mul<E: Pairing>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("{}/scalar-mul", curve_name::<E>()));

    let rng = &mut test_rng();
    let n = 100;

    let mut exps = vec![];
    exps.resize_with(n, || E::ScalarField::rand(rng));
    // the timing depends on the exponent
    let bases_projective = vec![E::G1::rand(rng); n];
    let bases_affine = vec![E::G1Affine::rand(rng); n];

    let _res: E::G1 = bases_affine[0].mul(exps[0]); // result of affine mul is projective

    let mut i = 0;
    group.bench_function("proj", |b|
        b.iter_with_setup(
            || {
                let pair = (bases_projective[i], exps[i]);
                i = (i + 1) % n;
                pair
            },
            |(base, exp)| base.mul(exp),
        ));

    let mut i = 0;
    group.bench_function("aff", |b|
        b.iter_with_setup(
            || {
                let pair = (bases_affine[i], exps[i]);
                i = (i + 1) % n;
                pair
            },
            |(base, exp)| base.mul(exp),
        ));

    group.finish();
}

fn coordinates_conversion<E: Pairing>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("{}/into", curve_name::<E>()));
    let rng = &mut test_rng();
    let projective = E::G1::rand(rng);
    let affine = E::G1Affine::rand(rng);
    group.bench_function("affine", |b| b.iter(|| projective.into_affine()));
    group.bench_function("projective", |b| b.iter(|| affine.into_group()));
    group.finish();
}

fn additions<E: Pairing>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("{}/addition", curve_name::<E>()));
    let rng = &mut test_rng();
    let a_projective = E::G1::rand(rng);
    let b_projective = E::G1::rand(rng);
    let a_affine = E::G1Affine::rand(rng);
    let b_affine = E::G1Affine::rand(rng);
    group.bench_function("projective", |b| b.iter(|| a_projective + b_projective));
    group.bench_function("affine", |b| b.iter(|| a_affine + b_affine));
    group.bench_function("mixed", |b| b.iter(|| a_projective + &b_affine));
    group.bench_function("doubling", |b| b.iter(|| a_projective.double()));
    // group.bench_function("doubling", |b| b.iter(|| CurveGroup::double_in_place(&mut b_projective)));
    group.finish();
}


criterion_group!(benches, scalar_mul::<BW6_761>, coordinates_conversion::<BW6_761>, additions::<BW6_761>);
criterion_main!(benches);