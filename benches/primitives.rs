use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

use ark_ff::UniformRand;
use ark_std::test_rng;

use ark_ff::PrimeField;
use ark_ec::{PairingEngine, ProjectiveCurve, AffineCurve};
use ark_bw6_761::BW6_761;
use fflonk::utils;


fn scalar_mul<E: PairingEngine>(c: &mut Criterion) {
    let n = 100;

    let rng = &mut test_rng();

    // mul is not constant time, but timing depends on the exponent, not the base, i guess
    let mut exps = vec![];
    exps.resize_with(n, || E::Fr::rand(rng));

    let bases_projective = vec![E::G1Projective::rand(rng); n];
    let bases_affine = vec![E::G1Affine::rand(rng); n];

    let _res: E::G1Projective = bases_affine[0].mul(exps[0]); // result of affine mul is projective

    let mut i = 0;
    c.bench_function("mul projective", |b|
        b.iter_with_setup(
            || {
                let pair = (bases_projective[i], exps[i]);
                i = (i + 1) % n;
                pair
            },
            |(base, exp)| base.mul(exp.into_repr()),
        ));

    let mut i = 0;
    c.bench_function("mul affine", |b|
        b.iter_with_setup(
            || {
                let pair = (bases_affine[i], exps[i]);
                i = (i + 1) % n;
                pair
            },
            |(base, exp)| base.mul(exp),
        ));
}

fn coordinates_conversion<E: PairingEngine>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let projective = E::G1Projective::rand(rng);
    let affine = E::G1Affine::rand(rng);
    c.bench_function("into_affine", |b| b.iter(|| projective.into_affine()));
    c.bench_function("into_projective", |b| b.iter(|| affine.into_projective()));
}

fn additions<E: PairingEngine>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let a_projective = E::G1Projective::rand(rng);
    let b_projective = E::G1Projective::rand(rng);
    let a_affine = E::G1Affine::rand(rng);
    let b_affine = E::G1Affine::rand(rng);
    c.bench_function("add projective", |b| b.iter(|| a_projective + b_projective));
    c.bench_function("add affine", |b| b.iter(|| a_affine + b_affine));
    c.bench_function("add mixed", |b| b.iter(|| a_projective.add_mixed(&b_affine)));
}

fn short_mul<E: PairingEngine>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let projective = E::G1Projective::rand(rng);
    let affine = E::G1Affine::rand(rng);
    let exp: E::Fr = u128::rand(rng).into();
    c.bench_function("128-bit projective", |b| b.iter(|| projective.mul(exp.into_repr())));
    c.bench_function("128-bit affine", |b| b.iter(|| affine.mul(exp)));
}

fn small_multiexp<E: PairingEngine>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let n = 10;

    let bases = (0..n).map(|_| E::G1Affine::rand(rng)).collect::<Vec<_>>();
    let exps_full = (0..n).map(|_| E::Fr::rand(rng)).collect::<Vec<_>>();
    let exps_128 = (0..n).map(|_| E::Fr::from(u128::rand(rng))).collect::<Vec<_>>();

    let mut group = c.benchmark_group("msm");
    group.bench_with_input(BenchmarkId::new("small-multiexp-full", n), &n, |b, _n| {
        b.iter(|| utils::small_multiexp(&exps_full, &bases))
    });
    group.bench_with_input(BenchmarkId::new("naive-multiexp-full", n), &n, |b, _n| {
        b.iter(|| utils::naive_multiexp(&exps_full, &bases))
    });
    group.bench_with_input(BenchmarkId::new("small-multiexp-128", n), &n, |b, _n| {
        b.iter(|| utils::small_multiexp(&exps_128, &bases))
    });
    group.bench_with_input(BenchmarkId::new("naive-multiexp-128", n), &n, |b, _n| {
        b.iter(|| utils::naive_multiexp(&exps_128, &bases))
    });
    group.finish();
}

fn small_multiexp_proj<E: PairingEngine>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let n = 10;

    let bases = (0..n).map(|_| E::G1Projective::rand(rng)).collect::<Vec<_>>();
    let exps_full = (0..n).map(|_| E::Fr::rand(rng)).collect::<Vec<_>>();
    let exps_128 = (0..n).map(|_| E::Fr::from(u128::rand(rng))).collect::<Vec<_>>();

    let mut group = c.benchmark_group("msm");
    group.bench_with_input(BenchmarkId::new("small-multiexp-proj-full", n), &n, |b, _n| {
        b.iter(|| utils::small_multiexp_proj(&exps_full, &bases))
    });
    group.bench_with_input(BenchmarkId::new("small-multiexp-proj2-full", n), &n, |b, _n| {
        b.iter(|| utils::_small_multiexp_proj2(&exps_full, &bases))
    });
    group.bench_with_input(BenchmarkId::new("small-multiexp-proj-128", n), &n, |b, _n| {
        b.iter(|| utils::small_multiexp_proj(&exps_128, &bases))
    });
    group.bench_with_input(BenchmarkId::new("small-multiexp-proj2-128", n), &n, |b, _n| {
        b.iter(|| utils::_small_multiexp_proj2(&exps_128, &bases))
    });
    group.finish();
}

criterion_group!(
    benches,
    scalar_mul::<BW6_761>,
    short_mul::<BW6_761>,
    coordinates_conversion::<BW6_761>,
    additions::<BW6_761>,
    small_multiexp::<BW6_761>,
    small_multiexp_proj::<BW6_761>,
);
criterion_main!(benches);