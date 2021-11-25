use ark_ff::PrimeField;
use ark_std::ops::{Add, Sub};
use ark_poly::univariate::DensePolynomial;
use crate::fflonk::Fflonk;
use crate::shplonk::ShplonkTranscript;

pub mod shplonk;
pub mod fflonk;
mod utils;

pub trait AdditiveCommitment<F: PrimeField>:
Sized
+ Add<Self, Output=Self>
+ Sub<Self, Output=Self>
+ core::iter::Sum<Self>
{
    fn mul(&self, by: F) -> Self;
}

/// Polynomial commitment scheme.
pub trait CommitmentScheme<F: PrimeField, P> {
    type G: AdditiveCommitment<F>;
    fn commit(&self, poly: &P) -> Self::G;
    fn commit_const(&self, c: &F) -> Self::G;
    /// Verifies that `p(x) = y` given the commitment to `p` and a proof.
    fn verify(&self, commitment: &Self::G, x: &F, y: F, proof: &Self::G) -> bool; // TODO: generalize proof type
}

type Poly<F> = DensePolynomial<F>; // currently SparsePolynomial doesn't implement UVPolynomial anyway

pub fn open<F, C, T>(
    fss: &[Vec<Poly<F>>], // vecs of polynomials to combine
    ts: &[usize], // lengths of each combination
    // TODO: ts can be inferred from li := len(fss[i]) as ti = min(x : x >= li and x | p-1)
    rootss: &[Vec<F>], // sets of opening points per a combined polynomial presented as t-th roots
    scheme: &C,
    transcript: &mut T,
) -> (C::G, C::G)
    where
        F: PrimeField,
        C: CommitmentScheme<F, DensePolynomial<F>>,
        T: ShplonkTranscript<F, C::G>,
{
    let k = fss.len();
    assert_eq!(k, ts.len());
    assert_eq!(k, rootss.len());
    let gs: Vec<Poly<F>> = fss.iter()
        .zip(ts.iter())
        .map(|(fs, t)| Fflonk::combine(*t, fs))
        .collect();
    let xss: Vec<_> = rootss.iter()
        .zip(ts.iter())
        .map(|(roots, t)|
            roots.iter()
                .flat_map(|root| Fflonk::<F, Poly<F>>::roots(*t, *root))
                .collect()
        ).collect();

    shplonk::open(&gs, &xss, scheme, transcript)
}

pub fn verify<F, C, T>(
    gcs: &[C::G],
    ts: &[usize],
    proof: (C::G, C::G),
    rootss: &[Vec<F>],
    vss: &[Vec<Vec<F>>],
    scheme: &C,
    transcript: &mut T,
) -> bool
    where
        F: PrimeField,
        C: CommitmentScheme<F, DensePolynomial<F>>,
        T: ShplonkTranscript<F, C::G>,
{
    let (xss, yss) = rootss.iter()
        .zip(vss.iter())
        .zip(ts.iter())
        .map(|((roots, vs), t)|
            Fflonk::<F, Poly<F>>::multiopening(*t, roots, vs)
        ).unzip();

    shplonk::verify(&gcs, proof, &xss, &yss, scheme, transcript)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use ark_std::rand::Rng;
    use std::marker::PhantomData;

    use ark_ff::Field;
    use ark_poly::{UVPolynomial, Polynomial};

    pub type F = ark_bw6_761::Fr;

    pub struct WrappedPolynomial<F: Field, P: UVPolynomial<F>>(pub P, PhantomData<F>);

    impl<F: Field, P: UVPolynomial<F>> WrappedPolynomial<F, P> {
        fn evaluate(&self, x: &F) -> F {
            self.0.evaluate(x)
        }
    }

    impl<F: PrimeField, P: UVPolynomial<F>> Add<Self> for WrappedPolynomial<F, P> {
        type Output = WrappedPolynomial<F, P>;

        fn add(self, other: WrappedPolynomial<F, P>) -> Self::Output {
            WrappedPolynomial(self.0 + other.0, PhantomData)
        }
    }

    impl<F: PrimeField, P: UVPolynomial<F>> Sub<Self> for WrappedPolynomial<F, P> {
        type Output = WrappedPolynomial<F, P>;

        fn sub(self, other: WrappedPolynomial<F, P>) -> Self::Output {
            let mut temp = self.0;
            temp -= &other.0; //TODO
            WrappedPolynomial(temp, PhantomData)
        }
    }

    impl<F: PrimeField, P: UVPolynomial<F>> core::iter::Sum<Self> for WrappedPolynomial<F, P> {
        fn sum<I: Iterator<Item=Self>>(iter: I) -> Self {
            iter.reduce(|a, b| a + b).unwrap()
        }
    }

    impl<F: PrimeField, P: UVPolynomial<F>> AdditiveCommitment<F> for WrappedPolynomial<F, P> {
        fn mul(&self, by: F) -> Self {
            let mut temp = P::zero(); //TODO
            temp += (by, &self.0);
            WrappedPolynomial(temp, PhantomData)
        }
    }

    pub struct IdentityCommitment {}

    impl<F: PrimeField, P: UVPolynomial<F>> CommitmentScheme<F, P> for IdentityCommitment {
        type G = WrappedPolynomial<F, P>;

        fn commit(&self, poly: &P) -> Self::G {
            WrappedPolynomial(poly.clone(), PhantomData)
        }

        fn commit_const(&self, c: &F) -> Self::G {
            WrappedPolynomial(P::from_coefficients_vec(vec![*c]), PhantomData)
        }

        fn verify(&self, commitment: &Self::G, x: &F, y: F, _proof: &Self::G) -> bool {
            commitment.evaluate(x) == y
        }
    }

    pub fn generate_test_data<R, F, C>(
        rng: &mut R,
        d: usize, // degree of polynomials
        t: usize, // number of polynomials
        max_m: usize, // maximal number of opening points per polynomial
        scheme: &C, // commitment scheme
    ) -> (
        Vec<Poly<F>>, // polynomials
        Vec<C::G>, // commitments
        Vec<Vec<F>>, // opening points per polynomial
        Vec<Vec<F>>, // evaluations per polynomial
    ) where
        R: Rng,
        F: PrimeField,
        C: CommitmentScheme<F, Poly<F>>
    {
        let ms: Vec<usize> = (0..t) // number of opening points per polynomial
            .map(|_| rng.gen_range(1..max_m))
            .collect();

        // polynomials
        let fs: Vec<Poly<F>> = (0..t)
            .map(|_| Poly::rand(d, rng))
            .collect();
        // commitments
        let fcs: Vec<_> = fs.iter()
            .map(|fi| scheme.commit(fi))
            .collect();
        // opening points per polynomial
        let xss: Vec<_> = ms.into_iter()
            .map(|m|
                (0..m).map(|_| F::rand(rng)
                ).collect::<Vec<_>>())
            .collect();
        // evaluations per polynomial
        let yss: Vec<_> = fs.iter()
            .zip(xss.iter())
            .map(|(f, xs)|
                xs.iter().map(
                    |x| f.evaluate(x))
                    .collect::<Vec<_>>()
            ).collect();

        (fs, fcs, xss, yss)
    }
}