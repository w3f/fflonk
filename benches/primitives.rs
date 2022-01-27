use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ark_bw6_761::{G1Projective, BW6_761};
use ark_ff::UniformRand;
use ark_std::test_rng;
use ark_ec::PairingEngine;
use ark_ec::group::Group;
use ark_ec::ProjectiveCurve;
use ark_ec::AffineCurve;

use ark_ff::PrimeField;



fn scalar_mul<E: PairingEngine>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let exp = E::Fr::rand(rng);
    let base_projective = E::G1Projective::rand(rng);
    let base_affine = base_projective.into_affine();
    let res: E::G1Projective = base_affine.mul(exp);
    c.bench_function("mul projective", |b| b.iter(|| base_projective.mul(exp.into_repr())));
    c.bench_function("mul affine", |b| b.iter(|| base_affine.mul(exp)));
}

criterion_group!(benches, scalar_mul::<BW6_761>);
criterion_main!(benches);