use std::process::exit;
use fftw::types::*;
use ::fhew::{
    *,
    BinGate::*,
    lwe::*,
    fhew::*
};
use rand::Rng;

fn help(cmd: &String) {
    eprintln!("\nusage: {} n\n", cmd);
    eprintln!("  Generate a secret key sk and evaluation key ek, and repeat the following test n times:");
    eprintln!("   - generate random bits b1,b2,b3,b4");
    eprintln!("   - compute ciphertexts c1, c2, c3 and c4 encrypting b1, b2, b3 and b4  under sk");
    eprintln!("   - homomorphically compute the encrypted (c1 NAND c2) NAND (c3 NAND c4)");
    eprintln!("   - decrypt all the intermediate results and check correctness");
    eprintln!("\n If any of the tests fails, print ERROR and stop immediately.\n");
    exit(0);
}
fn cleartext_gate(v1: bool, v2: bool, gate: BinGate) -> bool {
    match gate {
        OR => v1 || v2,
        AND => v1 && v2,
        NOR => !(v1 || v2),
        NAND => !(v1 && v2)
    }
}
fn eprint_gate(gate: BinGate) {
    match gate {
        OR => eprint!(" OR\t"),
        AND => eprint!(" AND\t"),
        NOR => eprint!(" NOT\t"),
        NAND => eprint!(" NAND\t")
    }
}
fn main() {
    // assert_eq!(q, 512);
    let mut rng = rand::thread_rng();
    let mut ffto: FFT = Default::default();

    // fftSetup(&mut ffto);
    let mut tTestMSB: RingFFT = Default::default();

    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        help(&args[0]);
    }
    let count: i32 = args[1].parse().unwrap();
    eprintln!("Setting up FHEW");
    fhew::setup(&mut ffto, &mut tTestMSB);
    eprint!("Generating secret key ... ");
    let mut lwe_sk: lwe::SecretKey = Default::default();
    lwe::keyGen(&mut lwe_sk, &mut rng);
    eprintln!("Done.\n");
    eprintln!("Generating evaluation key ...  this may take a while ... ");
    let mut ek: EvalKey = Default::default();
    fhew::keyGen(&mut ek, &lwe_sk, &mut rng, &mut ffto);
    eprintln!("Done.\n");
    eprintln!("Testing depth-2 homomorphic circuits {} times.", count);
    eprintln!("Circuit shape : (a GATE NOT(b)) GATE (c GATE d)\n");

    let (mut v1, mut v2): (i32, i32);
    let (mut sv1, mut sv2): (i32, i32) = (2, 2);
    let (mut se1, mut se2, mut e1, mut e2, mut e12): (CipherText, CipherText, CipherText, CipherText, CipherText)
        = Default::default();

    for i in 1..(3*count) {
        if i % 3 != 0 {
            v1 = rng.gen::<i32>() % 2;
            v2 = rng.gen::<i32>() % 2;
            lwe::encrypt(&mut e1, &lwe_sk, v1, &mut rng);
            lwe::encrypt(&mut e2, &lwe_sk, v2, &mut rng);
            if i % 3 == 1 {
                eprint!(" NOT\tEnc({}) = ", v2);
                let e2_temp = e2.clone();
                fhew::hom_not(&mut e2, &e2_temp);
                let notv2 = lwe::decrypt(&lwe_sk, &e2);
                eprintln!("Enc({})", v2);
                if !(notv2 == !v2) {
                    eprintln!("ERROR: incorrect NOT Homomorphic computation at iteration {}", i+1);
                    exit(1);
                }
                v2 = !v2;
            }
        } else {
            v1 = sv1;
            v2 = sv2;
            e1 = se1.clone();
            e2 = se2.clone();
        }
        let gate: BinGate = match rng.gen::<usize>() % 4 {
            0 => BinGate::OR,
            1 => BinGate::AND,
            2 => BinGate::NOR,
            3 => BinGate::NAND,
            _ => BinGate::OR
        };
        eprint!("Enc({})", v1);
        eprint_gate(gate);
        eprint!("Enc({}) = ", v2);

        fhew::hom_gate(&mut e12, gate, &ek, &e1, &e2, &mut ffto, &tTestMSB, &mut rng);
        let v12: i32 = lwe::decrypt(&lwe_sk, &e12);

        eprintln!("Enc({})", v12);

        match i % 3 {
            0 => eprintln!(""),
            1 => {
                sv1 = v12;
                se1 = e12.clone();
            },
            2 => {
                sv2 = v12;
                se2 = e12.clone();
            }
            _ => ()
        }

        if cleartext_gate(v1 != 0, v2 != 0, gate) != (v12 != 0) {
            eprintln!("\n ERROR: incorrect Homomorphic Gate computation at iteration {}", i+1);
            exit(1);
        }
    }
    eprintln!("\nPassed all tests!\n");
}