use rand::Rng;
use fftw::types::*;
use crate::{
    *,
    BinGate::*,
    lwe::*
};
use std::fs::File;
use std::io::{
    BufRead,
    Write
};
#[derive(Clone)]
struct CtModQ(Array1<CtModQ1>); // [CtModQ1; K2];
impl Default for CtModQ {
    fn default() -> Self {
        CtModQ(Array::default(K2))
    }
}
#[derive(Clone)]
struct CtModQ1(Array1<RingModQ>); // [RingModQ; 2];
impl Default for CtModQ1 {
    fn default() -> Self {
        CtModQ1(Array::default(2))
    }
}
#[derive(Clone)]
struct DctModQ(Array1<DctModQ1>); // [DctModQ1; K2];
impl Default for DctModQ {
    fn default() -> Self {
        DctModQ(Array::default(K2))
    }
}
#[derive(Clone)]
struct DctModQ1(Array1<RingModQ>); // [RingModQ; K2];
impl Default for DctModQ1 {
    fn default() -> Self {
        DctModQ1(Array::default(K2))
    }
}
#[derive(Clone)]
struct DctFFT(Array1<DctFFT1>); // [DctFFT1; K2];
impl Default for DctFFT {
    fn default() -> Self {
        DctFFT(Array::default(K2))
    }
}
#[derive(Clone)]
struct DctFFT1(Array1<RingFFT>); // [RingFFT; K2];
impl Default for DctFFT1 {
    fn default() -> Self {
        DctFFT1(Array::default(K2))
    }
}
#[derive(Clone)]
pub struct CtFFT(Array1<CtFFT1>); // [CtFFT1; K2];
impl Default for CtFFT {
    fn default() -> Self {
        CtFFT(Array::default(K2))
    }
}
#[derive(Clone)]
struct CtFFT1(Array1<RingFFT>); // [RingFFT; 2];
impl Default for CtFFT1 {
    fn default() -> Self {
        CtFFT1(Array::default(2))
    }
}
#[derive(Clone)]
pub struct BootstrappingKey(Array3<CtFFT>); // [[[CtFFT; BS_EXP]; BS_BASE]; n];
impl Default for BootstrappingKey {
    fn default() -> Self {
        BootstrappingKey(Array::default((n,BS_BASE,BS_EXP)))
    }
}
#[derive(Default,Clone)]
pub struct EvalKey {
    pub BSkey: BootstrappingKey,
    pub KSkey: lwe::SwitchingKey
}
pub fn setup(ffto: &mut FFT, tTestMSB: &mut RingFFT) {
    fft_setup(ffto);
    let mut tmsb: RingModQ = Default::default();
    tmsb.0[0] = Wrapping(-1);
    for i in 1..N {
        tmsb.0[i] = Wrapping(1);
    }
    fft_forward(ffto, &tmsb, tTestMSB);
}
fn fhewEncrypt(ct: &mut CtFFT, skFFT: &RingFFT, m: i32, rng: &mut rand::rngs::ThreadRng, ffto: &mut FFT) {
    let mut ai: RingFFT = Default::default();
    let mut res: CtModQ = Default::default();
    let qi = q as i32;
    let Ni = N as i32;
    let mut mm: i32 = (((m & qi) + qi) % qi) * (2*Ni/qi);
    let mut sign: i32 = 1;
    if mm >= Ni {
        mm -= Ni;
        sign = -1;
    }
    // println!("stage 1");
    for i in 0..K2 {
        for k in 0..(Ni as usize) {
            res.0[i].0[0].0[k] = rng.gen();
        }
        fft_forward(ffto, &((res.0[i]).0[0]), &mut ai); 
        for k in 0..(N2 as usize) {
            ai.0[k] = ai.0[k] * skFFT.0[k];
        }
        fft_backward(ffto, &ai, &mut res.0[i].0[1]);
        for k in 0..(Ni as usize) {
            // println!("i = {}, k = {}, res = {}", i, k, res.0[i].0[1].0[k]);
            res.0[i].0[1].0[k] += Wrapping(sample(&CHI1, rng) as i32);
        }
    }
    // println!("stage 2");
    for i in 0..K {
        res.0[2*i].0[0].0[mm as usize] += Wrapping(sign) * vgprime[i];
        res.0[2*i+1].0[1].0[mm as usize] += Wrapping(sign) * vgprime[i];
    }
    // println!("stage 3");
    for i in 0..K2 {
        for j in 0..2 {
            fft_forward(ffto, &res.0[i].0[j], &mut ct.0[i].0[j]);
        }
    }
    // println!("stage 4");
}
pub fn keyGen(ek: &mut EvalKey, lweSk: &lwe::SecretKey, rng: &mut rand::rngs::ThreadRng, ffto: &mut FFT) {
    let mut fhewSK: SecretKeyN = Default::default();
    keyGenN(&mut fhewSK, rng);
    switchingKeyGen(&mut ek.KSkey, &lweSk, &fhewSK, rng);

    let mut fhewSkFFT: RingFFT = Default::default();
    let mut fhew_rmq: RingModQ = Default::default();
    for i in 0..N {
        fhew_rmq.0[i] = Wrapping(fhewSK.0[i]);
    }
    fft_forward(ffto, &fhew_rmq, &mut fhewSkFFT);
    for i in 0..n {
        for j in 1..BS_BASE {
            for k in 0..BS_EXP {
                ek.BSkey.0[[i,j,k]] = Default::default();
                fhewEncrypt(&mut ek.BSkey.0[[i,j,k]], &fhewSkFFT, lweSk.0[i]*(j as i32)*(BS_TABLE[k] as i32), rng, ffto);
            }
        }
    }
}
fn addToAcc(acc: &mut CtFFT1, c: &CtFFT, ffto: &mut FFT) {
    let mut ct: CtModQ1 = Default::default();
    let mut dct: DctModQ1 = Default::default();
    let mut dctFFT: DctFFT1 = Default::default();
    for j in 0..2 {
        fft_backward(ffto, &acc.0[j], &mut ct.0[j]);
    }
    for j in 0..2 {
        for k in 0..N {
            let mut t: ZmodQ = ct.0[j].0[k] * v_inverse;
            for l in 0..K {
                let r: ZmodQ = Wrapping((t.0 << g_bits_32[l]) >> g_bits_32[l]);
                t = Wrapping((t - r).0 >> g_bits[l]);
                dct.0[j+2*1].0[k] = r;
            }
        }
    }
    for j in 0..K2 {
        fft_forward(ffto, &dct.0[j], &mut dctFFT.0[j]);
    }
    for j in 0..2 {
        for k in 0..N2 {
            acc.0[j].0[k] = c64::new(0.0,0.0);
            for l in 0..K2 {
                acc.0[j].0[k] += dctFFT.0[l].0[k] * c.0[l].0[j].0[k];
            }
        }
    }
}
fn initializeAcc(acc: &mut CtFFT1, m: i32, ffto: &mut FFT) {
    let mut res: CtModQ1 = Default::default();
    let qi = q as i32;
    let Ni = N as i32;
    let mut mm = (((m % qi) + qi) % qi) * (2*Ni/qi);
    let mut sign = 1;
    if mm >= Ni {
        mm -= Ni;
        sign = -1;
    }
    for j in 0..2 {
        for k in 0..N {
            res.0[j].0[k] = Wrapping(0);
        }
    }
    res.0[1].0[mm as usize] += Wrapping(sign) * vgprime[0];
    for j in 0..2 {
        fft_forward(ffto, &res.0[j], &mut acc.0[j]);
    }
}
fn memberTest(t: &RingFFT, c: &CtFFT1, ffto: &mut FFT) -> CipherTextQN {
    let mut temp: RingFFT = Default::default();
    let mut tempModQ: RingModQ = Default::default();
    let mut ct: CipherTextQN = Default::default();
    for i in 0..N2 {
        temp.0[i] = (c.0[0].0[i] * t.0[i]).conj();
    }
    fft_backward(ffto, &temp, &mut tempModQ);
    for i in 0..N {
        ct.a[i] = tempModQ.0[i];
    }
    for i in 0..N2 {
        temp.0[i] = c.0[1].0[i] * t.0[i];
    }
    fft_backward(ffto, &temp, &mut tempModQ);
    ct.b = v + tempModQ.0[0];
    ct
}
pub fn hom_gate(
    res: &mut CipherText, 
    gate: BinGate, 
    ek: &EvalKey, 
    ct1: &CipherText, 
    ct2: &CipherText, 
    ffto: &mut FFT, 
    tTestMSB: &RingFFT, 
    rng: &mut rand::rngs::ThreadRng
) {
    let mut e12: CipherText = Default::default();
    let qi = q as i32;
    for i in 0..n {
        if ((ct1.a[i] - ct2.a[i]) % qi) != 0 && ((ct1.a[i] + ct2.a[i]) % qi) != 0 {
            break;
        }
        if i == n - 1 {
            panic!("ERROR: Please only use independant ciphertexts as inputs.");
        }
    }
    for i in 0..n {
        e12.a[i] = (2*qi - (ct1.a[i] + ct2.a[i])) % qi;
    }
    e12.b = (GATE_CONST[gate as usize] as i32) - (ct1.b + ct2.b) % qi;
    let mut acc: CtFFT1 = Default::default();
    initializeAcc(&mut acc, (e12.b + qi/4) % qi, ffto);
    for i in 0..n {
        let mut a = (qi - e12.a[i] % qi) % qi;
        for k in 0..BS_EXP {
            let a0 = a % (BS_BASE as i32);
            if a0 != 0 {
                addToAcc(&mut acc, &ek.BSkey.0[[i,a0 as usize,k]], ffto);
            }
            a /= BS_BASE as i32;
        }
    }
    let eQN: CipherTextQN = memberTest(tTestMSB, &acc, ffto);
    let mut eQ: CipherTextQ = Default::default();
    keySwitch(&mut eQ, &ek.KSkey, &eQN);
    modSwitch(res, &eQ, rng);
}
pub fn hom_nand(
    res: &mut CipherText, 
    ek: &EvalKey, 
    ct1: &CipherText, 
    ct2: &CipherText, 
    ffto: &mut FFT, 
    tTestMSB: &RingFFT, 
    rng: &mut rand::rngs::ThreadRng) {
    hom_gate(res, NAND, ek, ct1, ct2, ffto, tTestMSB, rng)
}
pub fn hom_not(res: &mut CipherText, ct: &CipherText) {
    let qi = q as i32;
    for i in 0..n {
        res.a[i] = (qi - ct.a[i]) % qi;
    }
    res.b = (5 * qi / 4 - ct.b) % qi;
}

pub fn fwriteEK(ek: &EvalKey, f: &mut File) {
    for i in 0..n {
        for j in 1..BS_BASE {
            for k in 0..BS_EXP {
                for l in 0..K2 {
                    for m in 0..2 {
                        for n_ in 0..N2 {
                            writeln!(*f, "{}", ek.BSkey.0[[i,j,k]].0[l].0[m].0[n_]).unwrap();
                        }
                    }
                }
            }
        }
    }
    for i in 0..N {
        for j in 0..KS_BASE {
            for k in 0..KS_EXP {
                for l in 0..n {
                    writeln!(*f, "{}", ek.KSkey.0[[i,j,k]].a[l]).unwrap();
                }
                writeln!(*f, "{}", ek.KSkey.0[[i,j,k]].b).unwrap();
            }
        }
    }
    // f.sync_all().unwrap();
}
pub fn freadEK(f: &File) -> EvalKey {
    use std::str::FromStr;
    let mut ek: EvalKey = Default::default();
    let mut b = std::io::BufReader::new(f);
    let mut s = String::new();
    for i in 0..n {
        for j in 1..BS_BASE {
            for k in 0..BS_EXP {
                for l in 0..K2 {
                    for m in 0..2 {
                        for n_ in 0..N2 {
                            b.read_line(&mut s).unwrap();
                            ek.BSkey.0[[i,j,k]].0[l].0[m].0[n_] = c64::from_str(&s).unwrap();
                        }
                    }
                }
            }
        }
    }
    for i in 0..N {
        for j in 0..KS_BASE {
            for k in 0..KS_EXP {
                for l in 0..n {
                    b.read_line(&mut s).unwrap();
                    ek.KSkey.0[[i,j,k]].a[l] = Wrapping(i32::from_str(&s).unwrap());
                }
                b.read_line(&mut s).unwrap();
                ek.KSkey.0[[i,j,k]].b = Wrapping(i32::from_str(&s).unwrap());
            }
        }
    }
    // f.sync_all().unwrap();
    ek
}