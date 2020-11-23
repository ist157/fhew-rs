use std::fs::File;
use std::process::exit;
// use fftw::types::*;
use ::fhew::{
    *,
    BinGate::*,
    lwe::*,
    fhew::*
};
use rand::Rng;

fn help(cmd: &String) {
    eprintln!("\nusage: {} <count>\n", cmd);
    eprintln!("  Generate a secret key sk and evaluation key ek, and repeat the following test <count> times:");
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
        OR => eprint!("OR"),
        AND => eprint!("AND"),
        NOR => eprint!("NOR"),
        NAND => eprint!("NAND")
    }
}
fn main() {
    // assert_eq!(q, 512);
    // let mut rng = rand::thread_rng();
    let mut ffto: FFT = Default::default();

    // fftSetup(&mut ffto);
    let mut t_test_msb: RingFFT = RingFFT();

    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        help(&args[0]);
    }
    let count: i32 = args[1].parse().unwrap();
    eprintln!("Setting up FHEW");
    fhew::setup(&mut ffto, &mut t_test_msb);
    eprint!("Generating secret key ... ");
    let mut lwe_sk: lwe::SecretKey = SecretKey();
    lwe::key_gen(&mut lwe_sk);
    // dbg!(&lwe_sk);
    eprintln!("Done.\n");
    eprintln!("Generating evaluation key ...  this may take a while ... ");
    let mut ek: EvalKey = Default::default();
    fhew::key_gen(&mut ek, &lwe_sk, &mut ffto);

    // let mut f = File::create("key.txt").unwrap();
    // fwrite_ek(&ek, &mut f);

    eprintln!("Done.\n");
    eprintln!("Testing depth-2 homomorphic circuits {} times.", count);
    eprintln!("Circuit shape : (a GATE NOT(b)) GATE (c GATE d)\n");

    let (mut v1, mut v2): (i32, i32);
    let (mut sv1, mut sv2): (i32, i32) = (2, 2);
    let (mut se1, mut se2, mut e1, mut e2, mut e12): (CipherText, CipherText, CipherText, CipherText, CipherText)
        = Default::default();

    for i in 1..(3*count+1) {
        // if i != 1 {break;}
        if i % 3 != 0 { // 1,2
            v1 = (rand::thread_rng().gen::<u32>() % 2) as i32;
            v2 = (rand::thread_rng().gen::<u32>() % 2) as i32;
            // v1 = 0;
            // v2 = 0;
            lwe::encrypt(&mut e1, &lwe_sk, v1);
            lwe::encrypt(&mut e2, &lwe_sk, v2);
            // for t in 0..n {
            //     println!("t = {}, lwe_sk[t] = {}, e1.a[t] = {}, e2.a[t] = {}", t, lwe_sk[t], e1.a[t], e2.a[t]);
            // }
            if i % 3 == 1 { // 1
                eprint!("\tNOT\tEnc({}) = ", v2);
                let e2_temp = e2.clone();
                fhew::hom_not(&mut e2, &e2_temp);
                let notv2 = lwe::decrypt(&lwe_sk, &e2);
                eprintln!("Enc({})", notv2);
                // dbg!(v2,notv2,!v2,!notv2);
                if !(notv2 != v2 && notv2 * v2 == 0) {
                    eprintln!("ERROR: incorrect NOT Homomorphic computation at iteration {}", i+1);
                    exit(1);
                }
                v2 = if v2 == 0 {1} else {0};
            }
        } else { // 3
            v1 = sv1;
            v2 = sv2;
            e1 = se1.clone();
            e2 = se2.clone();
        }
        // let gate: BinGate = BinGate::NAND;
        let gate: BinGate = match rand::thread_rng().gen::<usize>() % 4 {
            0 => BinGate::OR,
            1 => BinGate::AND,
            2 => BinGate::NOR,
            3 => BinGate::NAND,
            _ => BinGate::OR
        };
        
        lwe::encrypt(&mut e1, &lwe_sk, v1);
        lwe::encrypt(&mut e2, &lwe_sk, v2);
        
        fhew::hom_gate(&mut e12, gate, &ek, &e1, &e2, &mut ffto, &t_test_msb);       
        let v12: i32 = lwe::decrypt(&lwe_sk, &e12);
        
        eprint!("Enc({})\t", v1);
        eprint_gate(gate);
        eprint!("\tEnc({}) = ", v2);

        eprint!("Enc({})", v12);
        eprintln!("");

        // for j in 0..n {
        //     println!("i = {}, j = {}, e1.a[j] = {}, e2.a[j] = {}, e12.a[j] = {}", i, j, e1.a[j], e2.a[j], e12.a[j]);
        // }

        // println!("i = {}\ne1.a = {:?}\ne2.a = {:?}\ne12.a = {:?}", i, e1.a, e2.a, e12.a);
        // println!("e1.b = {}, e2.b = {}, e12.b = {}", e1.b, e2.b, e12.b);

        match i % 3 {
            1 => {
                sv1 = v12;
                se1 = e12.clone();
            },
            2 => {
                sv2 = v12;
                se2 = e12.clone();
            },
            _ => eprintln!("")
        }
        // println!("i = {}, v1 = {}, v2 = {}, v12 = {}", i, v1, v2, v12);
        if cleartext_gate(v1 != 0, v2 != 0, gate) != (v12 != 0) {
            eprintln!("\n ERROR: incorrect Homomorphic Gate computation at iteration {}", i+1);
            exit(1);
        }
    }
    eprintln!("\nPassed all tests!\n");
}