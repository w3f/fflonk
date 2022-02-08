use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

use ark_ff::UniformRand;
use ark_ec::PairingEngine;
use ark_std::test_rng;

use fflonk::utils::ec;


fn small_multiexp_affine<E: PairingEngine>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let n = 10;

    let bases = (0..n).map(|_| E::G1Affine::rand(rng)).collect::<Vec<_>>();
    let exps_full = (0..n).map(|_| E::Fr::rand(rng)).collect::<Vec<_>>();
    let exps_128 = (0..n).map(|_| E::Fr::from(u128::rand(rng))).collect::<Vec<_>>();

    let mut group = c.benchmark_group("msm");
    group.bench_with_input(BenchmarkId::new("small-multiexp-full", n), &n, |b, _n| {
        b.iter(|| ec::small_multiexp_affine(&exps_full, &bases))
    });
    group.bench_with_input(BenchmarkId::new("naive-multiexp-full", n), &n, |b, _n| {
        b.iter(|| ec::naive_multiexp_affine(&exps_full, &bases))
    });
    group.bench_with_input(BenchmarkId::new("small-multiexp-128", n), &n, |b, _n| {
        b.iter(|| ec::small_multiexp_affine(&exps_128, &bases))
    });
    group.bench_with_input(BenchmarkId::new("naive-multiexp-128", n), &n, |b, _n| {
        b.iter(|| ec::naive_multiexp_affine(&exps_128, &bases))
    });
    group.finish();
}



criterion_group!(benches,
    small_multiexp_affine::<ark_bw6_761::BW6_761>,
);
criterion_main!(benches);