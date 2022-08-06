use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

use ark_ff::{UniformRand, PrimeField};
use ark_ec::PairingEngine;
use ark_std::test_rng;

use fflonk::utils::ec;
use ark_ec::msm::VariableBaseMSM;


fn small_multiexp_affine<E: PairingEngine>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let n = 10;

    let bases = (0..n).map(|_| E::G1Affine::rand(rng)).collect::<Vec<_>>();
    let exps_full = (0..n).map(|_| E::Fr::rand(rng)).collect::<Vec<_>>();
    let exps_128 = (0..n).map(|_| E::Fr::from(u128::rand(rng))).collect::<Vec<_>>();

    let mut group = c.benchmark_group("small-multiexp-affine");
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

fn small_multiexp_proj<E: PairingEngine>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let n = 10;

    let bases = (0..n).map(|_| E::G1Projective::rand(rng)).collect::<Vec<_>>();
    let exps_128 = (0..n).map(|_| E::Fr::from(u128::rand(rng))).collect::<Vec<_>>();

    let mut group = c.benchmark_group("small-multiexp-proj");
    group.bench_with_input(BenchmarkId::new("in_affine", n), &n, |b, _n| {
        b.iter(|| ec::small_multiexp_proj(&exps_128, &bases))
    });
    group.bench_with_input(BenchmarkId::new("in-proj", n), &n, |b, _n| {
        b.iter(|| ec::_small_multiexp_proj_2(&exps_128, &bases))
    });
    group.finish();
}

fn small_multiexp_vs_msm<E: PairingEngine>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let mut group = c.benchmark_group("small-multiexp-vs-msm");

    for n in [10, 20] {
        let bases = (0..n).map(|_| E::G1Affine::rand(rng)).collect::<Vec<_>>();

        let exps_full = (0..n).map(|_| E::Fr::rand(rng)).collect::<Vec<_>>();
        let exps_128 = (0..n).map(|_| E::Fr::from(u128::rand(rng))).collect::<Vec<_>>();

        let exps_full_repr = exps_full.iter().map(|exp| exp.into_bigint()).collect::<Vec<_>>();
        let exps_128_repr = exps_128.iter().map(|exp| exp.into_bigint()).collect::<Vec<_>>();


        group.bench_with_input(BenchmarkId::new("small-multiexp-full", n), &n, |b, _n| {
            b.iter(|| ec::small_multiexp_affine(&exps_full, &bases))
        });
        group.bench_with_input(BenchmarkId::new("var-base-msm-full", n), &n, |b, _n| {
            b.iter(|| <E::G1Projective as VariableBaseMSM>::msm_bigint(&bases, &exps_full_repr))
        });
        group.bench_with_input(BenchmarkId::new("small-multiexp-128", n), &n, |b, _n| {
            b.iter(|| ec::small_multiexp_affine(&exps_128, &bases))
        });
        group.bench_with_input(BenchmarkId::new("var-base-msm-128", n), &n, |b, _n| {
            b.iter(|| <E::G1Projective as VariableBaseMSM>::msm_bigint(&bases, &exps_128_repr))
        });
    }

    group.finish();
}


criterion_group!(benches,
    small_multiexp_affine::<ark_bw6_761::BW6_761>,
    small_multiexp_proj::<ark_bw6_761::BW6_761>,
    small_multiexp_vs_msm::<ark_bw6_761::BW6_761>,
);
criterion_main!(benches);