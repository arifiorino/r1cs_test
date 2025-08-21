#![allow(non_snake_case)]
use ark_bn254::{Bn254, Fr};
use ark_ff::Field;
use ark_marlin::Marlin;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::marlin_pc::MarlinKZG10;
use ark_relations::{
  lc,
  r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError},
};
use ark_std::One;
use blake2::Blake2s;
use rand::Rng;

#[derive(Clone)]
struct Circuit<F: Field> {
  N: usize,
  A: Vec<Option<F>>,
  B: Vec<Option<F>>,
  C: Vec<Option<F>>,
}

impl<ConstraintF: Field> ConstraintSynthesizer<ConstraintF> for Circuit<ConstraintF> {
  fn generate_constraints(self, cs: ConstraintSystemRef<ConstraintF>) -> Result<(), SynthesisError> {
    assert!(self.A.len() == self.N * self.N);
    assert!(self.B.len() == self.N * self.N);
    let one = cs.new_input_variable(|| Ok(ConstraintF::one())).unwrap();
    let A: Vec<_> = self.A.iter().map(|x| cs.new_witness_variable(|| x.ok_or(SynthesisError::AssignmentMissing)).unwrap()).collect();
    let B: Vec<_> = self.B.iter().map(|x| cs.new_witness_variable(|| x.ok_or(SynthesisError::AssignmentMissing)).unwrap()).collect();
    let C: Vec<_> = self.C.iter().map(|x| cs.new_witness_variable(|| x.ok_or(SynthesisError::AssignmentMissing)).unwrap()).collect();
    let mults: Vec<_> = (0..self.N * self.N * self.N)
      .map(|idx| {
        let i = idx / self.N / self.N;
        let j = idx / self.N % self.N;
        let k = idx % self.N;
        cs.new_witness_variable(|| Ok(self.A[i * self.N + k].unwrap() * self.B[j * self.N + k].unwrap())).unwrap()
      })
      .collect();
    for idx in 0..self.N * self.N * self.N {
      let i = idx / self.N / self.N;
      let j = idx / self.N % self.N;
      let k = idx % self.N;
      let _ = cs.enforce_constraint(lc!() + A[i * self.N + k], lc!() + B[j * self.N + k], lc!() + mults[idx]);
    }
    for idx in 0..self.N * self.N {
      let i = idx / self.N;
      let j = idx % self.N;
      let mut lc = lc!();
      for k in 0..self.N {
        lc = lc + mults[i * self.N * self.N + j * self.N + k];
      }
      let _ = cs.enforce_constraint(lc, lc!() + one, lc!() + C[idx]);
    }

    Ok(())
  }
}

type MultiPC = MarlinKZG10<Bn254, DensePolynomial<Fr>>;
type MarlinInst = Marlin<Fr, MultiPC, Blake2s>;

fn main() {
  let args: Vec<String> = std::env::args().collect();
  let N: usize = args.get(1).expect("Please provide a value for N").parse().expect("N must be a positive integer");

  let rng = &mut ark_std::test_rng();

  println!("Creating universal SRS...");
  let universal_srs = MarlinInst::universal_setup(N * N * N + N * N, 3 * N * N + N * N * N, 4 * N * N * N + 2 * N * N, rng).unwrap();

  println!("Creating A,B,C...");
  let A: Vec<_> = (0..N * N).map(|_| Some(Fr::from(rng.gen_range(-100..100)))).collect();
  let B: Vec<_> = (0..N * N).map(|_| Some(Fr::from(rng.gen_range(-100..100)))).collect();
  let C: Vec<_> = (0..N * N)
    .map(|idx| {
      let i = idx / N;
      let j = idx % N;
      let x = (0..N).map(|k| A[i * N + k].unwrap() * B[j * N + k].unwrap()).sum();
      Some(x)
    })
    .collect();

  let circ = Circuit {
    N: N,
    A: vec![None; N * N],
    B: vec![None; N * N],
    C: vec![None; N * N],
  };

  println!("Setting up circuit...");
  let (index_pk, index_vk) = MarlinInst::index(&universal_srs, circ).unwrap();

  let circ_with_witness = Circuit { N: N, A: A, B: B, C: C };

  println!("Proving...");
  let time = std::time::Instant::now();
  let proof = MarlinInst::prove(&index_pk, circ_with_witness, rng).unwrap();
  println!("Proof time: {:?}", time.elapsed());

  println!("Verifying...");
  assert!(MarlinInst::verify(&index_vk, &[Fr::one()], &proof, rng).unwrap());
  println!("Done");
}
