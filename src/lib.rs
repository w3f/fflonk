use ark_ff::PrimeField;
use ark_poly::univariate::{DensePolynomial, DenseOrSparsePolynomial};
use crate::fflonk::Fflonk;
use crate::shplonk::{ShplonkTranscript, Shplonk};
use crate::pcs::{PCS};
use std::marker::PhantomData;

pub mod shplonk;
pub mod fflonk;
pub mod pcs;
mod utils;


type Poly<F> = DensePolynomial<F>; // currently SparsePolynomial doesn't implement UVPolynomial anyway

pub trait EuclideanPolynomial<F: PrimeField> {
    fn divide_with_q_and_r(&self, divisor: &Poly<F>) -> (Poly<F>, Poly<F>);
}

impl<F: PrimeField> EuclideanPolynomial<F> for Poly<F> {
    fn divide_with_q_and_r(&self, divisor: &Poly<F>) -> (Poly<F>, Poly<F>) {
        let a: DenseOrSparsePolynomial<F> = self.into();
        let b: DenseOrSparsePolynomial<F> = divisor.into();
        a.divide_with_q_and_r(&b).unwrap()
    }
}


pub struct FflonkyKzg<F: PrimeField, CS: PCS<F>> {
    _field: PhantomData<F>,
    _pcs: PhantomData<CS>,
}

impl<F: PrimeField, CS: PCS<F>> FflonkyKzg<F, CS> {
    pub fn setup() -> (<Shplonk<F, CS> as PCS<F>>::CK, <Shplonk<F, CS> as PCS<F>>::VK) {
        Shplonk::<F, CS>::setup()
    }

    pub fn open<T: ShplonkTranscript<F, CS::G>>(
        ck: &<Shplonk<F, CS> as PCS<F>>::CK,
        fss: &[Vec<Poly<F>>], // vecs of polynomials to combine
        ts: &[usize], // lengths of each combination
        // TODO: ts can be inferred from li := len(fss[i]) as ti = min(x : x >= li and x | p-1)
        rootss: &[Vec<F>], // sets of opening points per a combined polynomial presented as t-th roots
        transcript: &mut T,
    ) -> (CS::G, CS::Proof)
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

        Shplonk::<F, CS>::open_many(ck, &gs, &xss, transcript)
    }

    pub fn verify<T: ShplonkTranscript<F, CS::G>>(
        vk: &<Shplonk<F, CS> as PCS<F>>::VK,
        gcs: &[CS::G],
        ts: &[usize],
        proof: (CS::G, CS::Proof),
        rootss: &[Vec<F>],
        vss: &[Vec<Vec<F>>],
        transcript: &mut T,
    ) -> bool
    {
        let (xss, yss) = rootss.iter()
            .zip(vss.iter())
            .zip(ts.iter())
            .map(|((roots, vs), t)|
                Fflonk::<F, Poly<F>>::multiopening(*t, roots, vs)
            ).unzip();

        Shplonk::<F, CS>::verify_many(vk, &gcs, &proof, &xss, &yss, transcript)
    }

    pub fn open_single<T: ShplonkTranscript<F, CS::G>>(
        ck: &<Shplonk<F, CS> as PCS<F>>::CK,
        fs: &[Poly<F>], // polynomials to combine
        t: usize, // lengths of the combination
        roots: &[F], // set of opening points presented as t-th roots
        transcript: &mut T,
    ) -> (CS::G, CS::Proof)
    {
        Self::open(ck, &[fs.to_vec()], &[t], &[roots.to_vec()], transcript)
    }

    pub fn verify_single<T: ShplonkTranscript<F, CS::G>>(
        vk: &<Shplonk<F, CS> as PCS<F>>::VK,
        gc: &CS::G,
        t: usize,
        proof: (CS::G, CS::Proof),
        roots: &[F],
        vss: &[Vec<F>], // evaluations per point // TODO: shplonk provides API with evals per polynomial
        transcript: &mut T,
    ) -> bool
    {
        Self::verify(vk, &[(*gc).clone()], &[t], proof, &[roots.to_vec()], &[vss.to_vec()], transcript)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use ark_std::rand::Rng;

    use ark_poly::{UVPolynomial, Polynomial};

    use crate::pcs::tests::IdentityCommitment;

    pub type F = ark_bw6_761::Fr;
    type TestFlonk = FflonkyKzg<F, IdentityCommitment>;

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
        let transcript = &mut (F::rand(rng), F::rand(rng));

        let (ck, vk) = TestFlonk::setup();

        let t = 4; // number of polynomials in a combination
        let m = 3; // number of opening points per a combination
        let d = 15;

        let (fs, roots, vss) = generate_test_data(rng, d, t, m);

        let g = Fflonk::combine(t, &fs);
        let gc = IdentityCommitment::commit(&(), &g);

        let proof = TestFlonk::open_single(&ck, &fs, t, &roots, transcript);
        assert!(TestFlonk::verify_single(&vk, &gc, t, proof, &roots, &vss, transcript));
    }

    #[test]
    fn test_fflonk() {
        let rng = &mut test_rng();
        let transcript = &mut (F::rand(rng), F::rand(rng));

        let (ck, vk) = TestFlonk::setup();

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
            .map(|(fs, t)| IdentityCommitment::commit(&(), &Fflonk::combine(t, &fs)))
            .collect();

        let proof = TestFlonk::open(&ck, &fss, &ts, &rootss, transcript);
        assert!(TestFlonk::verify(&vk, &gcs, &ts, proof, &rootss, &vsss, transcript));
    }
}