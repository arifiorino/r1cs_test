use ark_marlin::Marlin;
use ark_std::One;
use ark_marlin::rng::FiatShamirRng;
use ark_bn254::{Bn254, Fr};
use ark_ff::UniformRand;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::marlin_pc::MarlinKZG10;
use ark_std::ops::MulAssign;
use blake2::Blake2s;
use rand_chacha::ChaChaRng;
use ark_ff::Field;
use ark_relations::{
  lc,
  r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError},
};
use ark_std::marker::PhantomData;
use rand::Rng;
use std::ops::AddAssign;
use ark_ff::PrimeField;
use rayon::prelude::*;

#[derive(Clone)]
struct Circuit<F: Field> {
  N: usize,
  A: Vec<Option<F>>,
  B: Vec<Option<F>>,
  C: Vec<Option<F>>,
}

impl<ConstraintF: Field> ConstraintSynthesizer<ConstraintF> for Circuit<ConstraintF> {
  fn generate_constraints(
    self,
    cs: ConstraintSystemRef<ConstraintF>,
  ) -> Result<(), SynthesisError> {
    assert!(self.A.len() == self.N * self.N);
    assert!(self.B.len() == self.N * self.N);
    let one = cs.new_input_variable(||Ok(ConstraintF::one())).unwrap();
    let A:Vec<_> = self.A.iter().map(|x| cs.new_witness_variable(|| x.ok_or(SynthesisError::AssignmentMissing)).unwrap()).collect();
    let B:Vec<_> = self.B.iter().map(|x| cs.new_witness_variable(|| x.ok_or(SynthesisError::AssignmentMissing)).unwrap()).collect();
    let C:Vec<_> = self.C.iter().map(|x| cs.new_witness_variable(|| x.ok_or(SynthesisError::AssignmentMissing)).unwrap()).collect();
    let mults:Vec<_> = (0..self.N * self.N * self.N).map(|idx|{
      let i = idx / self.N / self.N;
      let j = idx / self.N % self.N;
      let k = idx % self.N;
      cs.new_witness_variable(|| {
        let mut x = self.A[i * self.N + k].unwrap();
        let y = self.B[j * self.N + k].unwrap();
        x.mul_assign(y);
        Ok(x)
      }).unwrap()
    }).collect();
    for idx in 0..self.N * self.N * self.N{
      let i = idx / self.N / self.N;
      let j = idx / self.N % self.N;
      let k = idx % self.N;
      cs.enforce_constraint(lc!() + A[i * self.N + k], lc!() + B[j * self.N + k], lc!() + mults[idx]);
    }
    for idx in 0..self.N * self.N{
      let i = idx / self.N;
      let j = idx % self.N;
      let mut lc = lc!();
      for k in 0..self.N{
        lc = lc + mults[i * self.N * self.N + j * self.N + k];
      }
      cs.enforce_constraint(lc, lc!() + one, lc!() + C[idx]);
    }

    Ok(())
  }
}

type MultiPC = MarlinKZG10<Bn254, DensePolynomial<Fr>>;
type MarlinInst = Marlin<Fr, MultiPC, Blake2s>;

fn main(){
  let rng = &mut ark_std::test_rng();

  let N = 16;

  let universal_srs = MarlinInst::universal_setup(N * N * N + N * N, N * N * 3 + N * N * N, 4 * N * N * N, rng).unwrap();


  let A:Vec<_> = (0..N * N).into_par_iter().map_init(rand::thread_rng, |rng, _| Some(Fr::from(rng.gen_range(-100..100)))).collect();
  let B:Vec<_> = (0..N * N).into_par_iter().map_init(rand::thread_rng, |rng, _| Some(Fr::from(rng.gen_range(-100..100)))).collect();
  let C:Vec<_> = (0..N * N).into_par_iter().map(|idx|{
    let i = idx / N;
    let j = idx % N;
    let mut x = A[i * N].unwrap() * B[j * N].unwrap();
    for k in 1..N{
      let y = A[i * N + k].unwrap() * B[j * N + k].unwrap();
      x.add_assign(y);
    }
    Some(x)
  }).collect();
  //println!("A:{:?}",A.iter().map(|x|x.unwrap().into_repr()).collect::<Vec<_>>());
  //println!("B:{:?}",B.iter().map(|x|x.unwrap().into_repr()).collect::<Vec<_>>());
  //println!("C:{:?}",C.iter().map(|x|x.into_repr()).collect::<Vec<_>>());

  let circ = Circuit {
    N: N,
    A: A,
    B: B,
    C: C,
  };

  let circ2 = Circuit {
    N: N,
    A: vec![None;N * N],
    B: vec![None;N * N],
    C: vec![None;N * N],
  };

  println!("Setting up");
  let (index_pk, index_vk) = MarlinInst::index(&universal_srs, circ2).unwrap();
  println!("Called index");

  let time = std::time::Instant::now();
  let proof = MarlinInst::prove(&index_pk, circ, rng).unwrap();
  println!("Proved {:?}",time.elapsed());

  assert!(MarlinInst::verify(&index_vk, &[Fr::one()], &proof, rng).unwrap());
  println!("Called verifier");
}
