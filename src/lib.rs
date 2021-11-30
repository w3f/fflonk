use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use crate::fflonk::Fflonk;
use crate::shplonk::ShplonkTranscript;
use crate::pcs::CommitmentScheme;

pub mod shplonk;
pub mod fflonk;
mod utils;
mod pcs;

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

    use ark_poly::{UVPolynomial, Polynomial};

    use crate::pcs::tests::IdentityCommitment;

    pub type F = ark_bw6_761::Fr;

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