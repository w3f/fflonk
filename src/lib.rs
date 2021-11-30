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
+ Clone
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

pub fn open_single<F, C, T>(
    fs: &[Poly<F>], // polynomials to combine
    t: usize, // lengths of the combination
    roots: &[F], // set of opening points presented as t-th roots
    scheme: &C,
    transcript: &mut T,
) -> (C::G, C::G)
    where
        F: PrimeField,
        C: CommitmentScheme<F, DensePolynomial<F>>,
        T: ShplonkTranscript<F, C::G>,
{
    open(&[fs.to_vec()], &[t], &[roots.to_vec()], scheme, transcript)
}

pub fn verify_single<F, C, T>(
    gc: &C::G,
    t: usize,
    proof: (C::G, C::G),
    roots: &[F],
    vss: &[Vec<F>], // evaluations per point // TODO: shplonk provides API with evals per polynomial
    scheme: &C,
    transcript: &mut T,
) -> bool
    where
        F: PrimeField,
        C: CommitmentScheme<F, DensePolynomial<F>>,
        T: ShplonkTranscript<F, C::G>,
{
    verify(&[(*gc).clone()], &[t], proof, &[roots.to_vec()], &[vss.to_vec()], scheme, transcript)
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use ark_std::rand::Rng;
    use ark_std::marker::PhantomData;

    use ark_ff::Field;
    use ark_poly::{UVPolynomial, Polynomial};

    pub type F = ark_bw6_761::Fr;

    #[derive(Clone)]
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

    fn generate_test_data<R, F>(
        rng: &mut R,
        d: usize, // degree of polynomials
        t: usize, // number of polynomials
        m: usize, // number of opening points
    ) -> (
        Vec<Poly<F>>, // polynomials
        Vec<F>, // roots of evaluation points
        Vec<Vec<F>>, // evaluations per point
    ) where
        R: Rng,
        F: PrimeField,
    {
        // polynomials
        let fs: Vec<Poly<F>> = (0..t)
            .map(|_| Poly::rand(d, rng))
            .collect();

        let roots: Vec<_> = (0..m)
            .map(|_| F::rand(rng))
            .collect();

        let xs: Vec<F> = roots.iter() // opening points
            .map(|root| root.pow([t as u64]))
            .collect();

        // evaluations per point
        let vss: Vec<_> = xs.iter()
            .map(|x|
                fs.iter()
                    .map(|f| f.evaluate(x))
                    .collect::<Vec<_>>()
            ).collect();

        (fs, roots, vss)
    }

    #[test]
    fn test_fflonk_single() {
        let rng = &mut test_rng();
        let scheme = IdentityCommitment {};
        let transcript = &mut (F::rand(rng), F::rand(rng));

        let t = 4; // number of polynomials in a combination
        let m = 3; // number of opening points per a combination
        let d = 15;

        let (fs, roots, vss) = generate_test_data(rng, d, t, m);

        let g = Fflonk::combine(t, &fs);
        let gc = scheme.commit(&g);

        let proof = open_single(&fs, t, &roots, &scheme, transcript);
        assert!(verify_single(&gc, t, proof, &roots, &vss, &scheme, transcript));
    }

    #[test]
    fn test_fflonk() {
        let rng = &mut test_rng();
        let scheme = IdentityCommitment {};
        let transcript = &mut (F::rand(rng), F::rand(rng));

        let ds = [31, 15];
        let ts = [2, 4]; // number of polynomials in a combination
        let ms = [2, 2]; // number of opening points per a combination

        let mut fss = vec![];
        let mut rootss = vec![];
        let mut vsss = vec![];
        for ((d, t), m) in ds.into_iter()
            .zip(ts)
            .zip(ms)
        {
            let (fs, roots, vss) = generate_test_data(rng, d, t, m);
            fss.push(fs);
            rootss.push(roots);
            vsss.push(vss);
        }

        let gcs: Vec<_> = fss.iter()
            .zip(ts)
            .map(|(fs, t)| scheme.commit(&Fflonk::combine(t, &fs)))
            .collect();

        let proof = open(&fss, &ts, &rootss, &scheme, transcript);
        assert!(verify(&gcs, &ts, proof, &rootss, &vsss, &scheme, transcript));
    }
}