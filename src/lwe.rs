use rand::Rng;
use crate::*;

#[derive(Clone)]
pub struct CipherText {
    pub a: Array1<i32>, // [isize; n],
    pub b: i32
}
impl Default for CipherText {
    fn default() -> Self {
        CipherText {
            a: Array::default(n),
            b: 0
        }
    }
}
#[derive(Clone)]
pub struct CipherTextQ {
    pub a: Array1<ZmodQ>, // [ZmodQ; n],
    pub b: ZmodQ
}
impl Default for CipherTextQ {
    fn default() -> Self {
        CipherTextQ {
            a: Array::default(n),
            b: Default::default()
        }
    }
}
#[derive(Clone)]
pub struct CipherTextQN {
    pub a: Array1<ZmodQ>, // [ZmodQ; n],
    pub b: ZmodQ
}
impl Default for CipherTextQN {
    fn default() -> Self {
        CipherTextQN {
            a: Array::default(n),
            b: Default::default()
        }
    }
}

#[derive(Clone)]
pub struct SecretKey(pub Array1<i32>); // [isize; n];
impl Default for SecretKey {
    fn default() -> Self {
        SecretKey(Array::default(n))
    }
}
#[derive(Clone)]
pub struct SecretKeyN(pub Array1<i32>); // [isize; N];
impl Default for SecretKeyN {
    fn default() -> Self {
        SecretKeyN(Array::default(N))
    }
}

const qi: i32 = q as i32;

pub fn keyGen(sk: &mut SecretKey, rng: &mut rand::rngs::ThreadRng) {
    loop {
        let mut s = 0;
        let mut ss = 0;
        for i in 0..n {
            sk.0[i] = sample(&CHI_BINARY, rng);
            s += sk.0[i];
            ss += (sk.0[i]).abs();
        }
        if s.abs() > 5 || (ss - (n as i32) / 2).abs() > 5 {
            continue;
        } else {
            break;
        }
    }
}
pub fn keyGenN(sk: &mut SecretKeyN, rng: &mut rand::rngs::ThreadRng) {
    for i in 0..N {
        sk.0[i] = sample(&CHI1, rng);
    }
}
pub fn encrypt(ct: &mut CipherText, sk: &SecretKey, m: i32, rng: &mut rand::rngs::ThreadRng) {
    ct.b = (m % 4) * qi / 4 + sample(&CHI3, rng);
    for i in 0..n {
        ct.a[i] = rng.gen::<i32>() % qi;
        ct.b = (ct.b + ct.a[i] * sk.0[i]) % qi;
    }
}
pub fn decrypt(sk: &SecretKey, ct: &CipherText) -> i32 {
    let mut r = ct.b;
    for i in 0..n {
        r -= ct.a[i] * sk.0[i];
    }
    r = ((r % qi) + qi + qi/8) % qi;
    4 * r / qi
}

#[derive(Clone)]
pub struct SwitchingKey(pub Array3<CipherTextQ>); // [[[CipherTextQ; KS_EXP]; KS_BASE]; N];
impl Default for SwitchingKey {
    fn default() -> Self {
        SwitchingKey(Array::default((N,KS_BASE,KS_EXP)))
    }
}

pub fn switchingKeyGen(res: &mut SwitchingKey, new_sk: &SecretKey, old_sk: &SecretKeyN, rng: &mut rand::rngs::ThreadRng) {
    for i in 0..N {
        for j in 0..KS_BASE {
            for k in 0..KS_EXP {
                let mut ct: CipherTextQ = Default::default();
                ct.b = - Wrapping(old_sk.0[i]) * Wrapping(j as i32) * Wrapping(KS_TABLE[k]) + Wrapping(sample(&CHI2, rng));
                for l in 0..n {
                    ct.a[l] = rng.gen();
                    let addend = ct.a[l] * Wrapping(new_sk.0[l]);
                    // dbg!(ct.b, addend, isize::MAX, isize::MIN);
                    // let upper = addend <= isize::MAX - ct.b;
                    // let lower = addend >= isize::MIN - ct.b;
                    // ct.b = if upper && lower {
                    //     ct.b + addend
                    // } else if !upper && lower {
                    //     isize::MIN + (addend - (isize::MAX - ct.b)) - 1
                    // } else {
                    //     isize::MAX - (- addend - (ct.b - isize::MIN)) + 1
                    // };
                    ct.b = ct.b + addend;
                    // dbg!(ct.b);
                }
                res.0[[i,j,k]] = ct;
            }
        }
    }
}
pub fn keySwitch(res: &mut CipherTextQ, ksk: &SwitchingKey, ct: &CipherTextQN) {
    for k in 0..n {
        res.a[k] = Wrapping(0);
    }
    res.b = ct.b;
    for i in 0..N {
        let a: UZmodQ = Wrapping(- ct.a[i].0 as u32);
        for j in 0..KS_BASE {
            let a0: UZmodQ = a % Wrapping(KS_BASE as u32);
            for k in 0..n {
                res.a[k] -= ksk.0[[i,a0.0 as usize,j]].a[k];
                res.b -= ksk.0[[i,a0.0 as usize,j]].b;
            }
        }
    }
}
pub fn modSwitch(ct: &mut CipherText, c: &CipherTextQ, rng: &mut rand::rngs::ThreadRng) {
    for i in 0..n {
        ct.a[i] = round_qQ(c.a[i]);
    }
    ct.b = round_qQ(c.b);
}
pub fn round_qQ(v_: ZmodQ) -> i32 {
    (0.5 + (v_.0 as f64) * (q as f64) / (Q as f64)).floor() as i32
}
