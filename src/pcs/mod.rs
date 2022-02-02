pub mod kzg;

use ark_ff::PrimeField;

use ark_std::iter::Sum;
use ark_std::ops::{Add, Sub};
use ark_std::fmt::Debug;
use ark_std::rand::Rng;
use ark_serialize::*;

use crate::Poly;
use ark_std::Zero;
use ark_poly::Polynomial;


pub trait CommitmentSpace<F: PrimeField>:
Sized
+ Clone
+ Add<Self, Output=Self>
+ Sub<Self, Output=Self>
+ Sum<Self>
{
    fn mul(&self, by: F) -> Self;
}


/// Can be used to commit and open commitments to DensePolynomial<F> of degree up to max_degree.
pub trait CommitterKey: Clone + Debug + CanonicalSerialize + CanonicalDeserialize {
    /// Maximal degree of a polynomial supported.
    fn max_degree(&self) -> usize;

    /// Maximal number of evaluations supported when committing in the Lagrangian base.
    fn max_evals(&self) -> usize {
        self.max_degree() + 1
    }
}


/// Can be used to verify openings to commitments.
pub trait VerifierKey: Clone + Debug {
    /// Maximal number of openings that can be verified.
    fn max_points(&self) -> usize {
        1
    }
}


pub trait PcsParams<CK, VK> {
    fn ck(&self) -> CK;
    //TODO: trim
    fn vk(&self) -> VK;
}


/// Polynomial commitment scheme.
pub trait PCS<F: PrimeField> {
    type G: CommitmentSpace<F>;

    type Proof;

    type CK: CommitterKey;
    type VK: VerifierKey + Into<Self::CK>;
    type Params: PcsParams<Self::CK, Self::VK>;

    fn setup<R: Rng>(max_degree: usize, rng: &mut R) -> Self::Params;

    fn commit(ck: &Self::CK, p: &Poly<F>) -> Self::G;

    fn open(ck: &Self::CK, p: &Poly<F>, x: F) -> Self::Proof; //TODO: eval?

    fn verify(vk: &Self::VK, c: Self::G, x: F, z: F, proof: Self::Proof) -> bool;
}

/// Represents a claim that f(x) = y, for some f such that CS::commit(f) = c.
/// Notion of "some f" depends on the soundness properties of the commitment scheme.
#[derive(Clone)]
pub struct Claim<F: PrimeField, CS: PCS<F>> {
    c: CS::G,
    x: F,
    y: F,
}

/// Represents an opening claim together with an alleged proof.
#[derive(Clone)]
pub struct Opening<F: PrimeField, CS: PCS<F>> {
    claim: Claim<F, CS>,
    proof: CS::Proof,
}

impl<F: PrimeField, CS: PCS<F>> Claim<F, CS> {
    fn new(poly: &Poly<F>, at: F, ck: &CS::CK) -> Claim<F, CS> {
        Claim {
            c: CS::commit(ck, poly),
            x: at,
            y: poly.evaluate(&at),
        }
    }

    fn aggregate_at_same_x(claims: &[Claim<F, CS>], rs: &[F]) -> Claim<F, CS> {
        assert_eq!(claims.len(), rs.len());

        let mut iter_over_xs = claims.iter().map(|cl| cl.x);
        let same_x = iter_over_xs.next().expect("claims is empty");
        assert!(iter_over_xs.all(|x| x == same_x));

        let (rcs, rys): (Vec<CS::G>, Vec<F>) = claims.iter().zip(rs.iter())
            .map(|(cl, &r)| (cl.c.mul(r), r * cl.y)).unzip();

        Claim {
            c: rcs.into_iter().sum(),
            x: same_x,
            y: rys.iter().sum(),
        }
    }
}

fn aggregate_polys<F: PrimeField>(polys: &[Poly<F>], rs: &[F]) -> Poly<F> {
    assert_eq!(polys.len(), rs.len());
    polys.iter().zip(rs.iter())
        .map(|(p, &r)| p * r)
        .fold(Poly::zero(), |acc, p| acc + p)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    use ark_ff::Zero;
    use ark_poly::{Polynomial, UVPolynomial};

    use crate::Poly;
    use ark_std::test_rng;
    use crate::pcs::kzg::KZG;

    pub(crate) type TestCurve = ark_bls12_381::Bls12_381;
    pub(crate) type TestField = ark_bls12_381::Fr;
    pub(crate) type TestKzg = KZG::<TestCurve>;

    #[derive(Clone)]
    pub struct WrappedPolynomial<F: PrimeField>(pub Poly<F>);

    impl<F: PrimeField> WrappedPolynomial<F> {
        fn evaluate(&self, x: &F) -> F {
            self.0.evaluate(x)
        }
    }

    impl<F: PrimeField> Add<Self> for WrappedPolynomial<F> {
        type Output = WrappedPolynomial<F>;

        fn add(self, other: WrappedPolynomial<F>) -> Self::Output {
            WrappedPolynomial(self.0 + other.0)
        }
    }

    impl<F: PrimeField> Sub<Self> for WrappedPolynomial<F> {
        type Output = WrappedPolynomial<F>;

        fn sub(self, other: WrappedPolynomial<F>) -> Self::Output {
            let mut temp = self.0;
            temp -= &other.0; //TODO
            WrappedPolynomial(temp)
        }
    }

    impl<F: PrimeField> core::iter::Sum<Self> for WrappedPolynomial<F> {
        fn sum<I: Iterator<Item=Self>>(iter: I) -> Self {
            iter.reduce(|a, b| a + b).unwrap()
        }
    }

    impl<F: PrimeField> CommitmentSpace<F> for WrappedPolynomial<F> {
        fn mul(&self, by: F) -> Self {
            let mut temp = Poly::zero(); //TODO
            temp += (by, &self.0);
            WrappedPolynomial(temp)
        }
    }


    impl CommitterKey for () {
        fn max_degree(&self) -> usize {
            usize::MAX >> 1
        }
    }

    impl VerifierKey for () {
        fn max_points(&self) -> usize {
            1
        }
    }


    impl PcsParams<(), ()> for () {
        fn ck(&self) -> () {
            ()
        }

        fn vk(&self) -> () {
            ()
        }
    }


    pub struct IdentityCommitment {}

    impl<F: PrimeField> PCS<F> for IdentityCommitment {
        type G = WrappedPolynomial<F>;
        type Params = ();
        type Proof = ();
        type CK = ();
        type VK = ();

        fn setup<R: Rng>(_max_degree: usize, _rng: &mut R) -> Self::Params {
            ()
        }

        fn commit(_ck: &(), p: &Poly<F>) -> Self::G {
            WrappedPolynomial(p.clone())
        }

        fn open(_ck: &(), _p: &Poly<F>, _x: F) -> Self::Proof {
            ()
        }

        fn verify(_vk: &(), c: Self::G, x: F, z: F, _proof: Self::Proof) -> bool {
            c.evaluate(&x) == z
        }
    }

    fn _test_aggregate_polys<F: PrimeField, CS: PCS<F>>() {
        let rng = &mut test_rng();
        let d = 10;
        let params = CS::setup(d, rng);
        let (ck, vk) = (params.ck(), params.vk());

        assert!(aggregate_polys::<F>(&[], &[]).is_zero());

        let (polys, rs): (Vec<_>, Vec<_>) = (0..4)
            .map(|_| (Poly::<F>::rand(d, rng), F::rand(rng)))
            .unzip();
        let agg_poly = aggregate_polys(&polys, &rs);

        let same_x = F::rand(rng);
        let claims_at_same_x = polys.iter().map(|p|
            Claim::<F, CS>::new(p, same_x, &ck)
        ).collect::<Vec<_>>();
        let agg_claim = Claim::aggregate_at_same_x(&claims_at_same_x, &rs);

        let agg_proof = CS::open(&ck, &agg_poly, same_x);
        assert!(CS::verify(&vk, agg_claim.c, agg_claim.x, agg_claim.y, agg_proof));
    }

    #[test]
    fn test_aggregate_polys() {
        _test_aggregate_polys::<TestField, IdentityCommitment>();
        _test_aggregate_polys::<TestField, TestKzg>()
    }
}