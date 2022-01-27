use ark_ff::PrimeField;
use ark_poly::univariate::{DensePolynomial, DenseOrSparsePolynomial};
use crate::fflonk::Fflonk;
use crate::shplonk::{ShplonkTranscript, Shplonk};
use crate::pcs::{PCS, PcsParams};
use std::marker::PhantomData;
use ark_std::rand::Rng;

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

    pub fn setup<R: Rng>(max_degree: usize, rng: &mut R) -> <Shplonk<F, CS> as PCS<F>>::Params {
        <Shplonk<F, CS> as PCS<F>>::setup(max_degree, rng)
    }

    pub fn open<T: ShplonkTranscript<F, CS::G>>(
        ck: &<<Shplonk<F, CS> as PCS<F>>::Params as PcsParams>::CommitterKey,
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
        vk: &<<Shplonk<F, CS> as PCS<F>>::Params as PcsParams>::VerifierKey,
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

        Shplonk::<F, CS>::verify_many(vk, &gcs, proof, &xss, &yss, transcript)
    }

    pub fn open_single<T: ShplonkTranscript<F, CS::G>>(
        ck: &<<Shplonk<F, CS> as PCS<F>>::Params as PcsParams>::CommitterKey,
        fs: &[Poly<F>], // polynomials to combine
        t: usize, // lengths of the combination
        roots: &[F], // set of opening points presented as t-th roots
        transcript: &mut T,
    ) -> (CS::G, CS::Proof)
    {
        Self::open(ck, &[fs.to_vec()], &[t], &[roots.to_vec()], transcript)
    }

    pub fn verify_single<T: ShplonkTranscript<F, CS::G>>(
        vk: &<<Shplonk<F, CS> as PCS<F>>::Params as PcsParams>::VerifierKey,
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
    use ark_bw6_761::{BW6_761, Fr};
    use crate::pcs::kzg::KZG;

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

    fn _test_fflonk_single<F: PrimeField, CS: PCS<F>>() {
        let rng = &mut test_rng();
        let transcript = &mut (F::rand(rng), F::rand(rng));

        let params = FflonkyKzg::<F, CS>::setup(123, rng);

        let t = 4; // number of polynomials in a combination
        let m = 3; // number of opening points per a combination
        let d = 15;

        let (fs, roots, vss) = generate_test_data(rng, d, t, m);

        let g = Fflonk::combine(t, &fs);
        let gc = CS::commit(&params.ck(), &g);

        let proof = FflonkyKzg::<F, CS>::open_single(&params.ck(), &fs, t, &roots, transcript);
        assert!(FflonkyKzg::<F, CS>::verify_single(&params.rk(), &gc, t, proof, &roots, &vss, transcript));
    }

    fn _test_fflonk<F: PrimeField, CS: PCS<F>>() {
        let rng = &mut test_rng();
        let transcript = &mut (F::rand(rng), F::rand(rng));

        let params = FflonkyKzg::<F, CS>::setup(123, rng);

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
            .map(|(fs, t)| CS::commit(&params.ck(), &Fflonk::combine(t, &fs)))
            .collect();

        let proof = FflonkyKzg::<F, CS>::open(&params.ck(), &fss, &ts, &rootss, transcript);
        assert!(FflonkyKzg::<F, CS>::verify(&params.rk(), &gcs, &ts, proof, &rootss, &vsss, transcript));
    }

    #[test]
    fn test_fflonk_single_id() {
        _test_fflonk_single::<Fr, IdentityCommitment>();
    }

    #[test]
    fn test_fflonk_single_kzg() {
        _test_fflonk_single::<Fr, KZG<BW6_761>>();
    }

    #[test]
    fn test_fflonk_id() {
        _test_fflonk::<Fr, IdentityCommitment>();
    }

    #[test]
    fn test_fflonk_kzg() {
        _test_fflonk::<Fr, KZG<BW6_761>>();
    }
}