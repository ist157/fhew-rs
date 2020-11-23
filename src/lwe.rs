use rand::Rng;
use crate::*;
use ndarray::parallel::prelude::*;
use ndarray::linalg::*;
use rayon::prelude::*;

#[derive(Clone,Debug)]
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
#[derive(Clone,Debug)]
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
#[derive(Clone,Debug)]
pub struct CipherTextQN {
    pub a: Array1<ZmodQ>, // [ZmodQ; n],
    pub b: ZmodQ
}
impl Default for CipherTextQN {
    fn default() -> Self {
        CipherTextQN {
            a: Array::default(N),
            b: Default::default()
        }
    }
}

// #[derive(Clone,Debug)]
// pub struct SecretKey(pub Array1<i32>); // [isize; n];
// impl Default for SecretKey {
//     fn default() -> Self {
//         SecretKey(Array::default(n))
//     }
// }
pub type SecretKey = Array<i32,Ix1>;
pub fn SecretKey() -> SecretKey {
    Array::default(n)
}
// #[derive(Clone,Debug)]
// pub struct SecretKeyN(pub Array1<i32>); // [isize; N];
// impl Default for SecretKeyN {
//     fn default() -> Self {
//         SecretKeyN(Array::default(N))
//     }
// }
pub type SecretKeyN = Array<i32,Ix1>;
pub fn SecretKeyN() -> SecretKeyN {
    Array::default(N)
}

// #[derive(Clone,Debug)]
// pub struct SwitchingKey(pub Array3<CipherTextQ>); // [[[CipherTextQ; KS_EXP]; KS_BASE]; N];
// impl Default for SwitchingKey {
//     fn default() -> Self {
//         SwitchingKey(Array::default((N,KS_BASE,KS_EXP)))
//     }
// }
pub type SwitchingKey = Array<CipherTextQ,Ix3>;
pub fn SwitchingKey() -> SwitchingKey {
    let k = Array::default((N,KS_BASE,KS_EXP));
    eprintln!("evaluation key 2 -> switching key zeroing complete");
    k
}

const qi: i32 = q as i32;

pub fn key_gen(sk: &mut SecretKey) {
    // loop {
        let mut s = 0;
        let mut ss = 0;
        for i in 0..n {
            // sk[i] = sample(&CHI_BINARY);
            sk[i] = BINARY_TABLE[i%3].floor() as i32;
            s += sk[i];
            ss += (sk[i]).abs();
        }
    //     if s.abs() > 5 || (ss - (n as i32) / 2).abs() > 5 {
    //         continue;
    //     } else {
    //         break;
    //     }
    // }
}
pub fn key_gen_N(sk: &mut SecretKeyN) {
    for i in 0..N {
        // sk[i] = sample(&CHI1);
        sk[i] = (i as i32) % 2;
    }
}
pub fn encrypt(ct: &mut CipherText, sk: &SecretKey, m: i32) {
    ct.b = (m % 4) * qi / 4; // + sample(&CHI3);
    for i in 0..n {
        // ct.a[i] = rand::thread_rng().gen::<i32>() % qi;
        ct.a[i] = (i as i32) % qi;
        ct.b = (ct.b + ct.a[i] * sk[i]) % qi;
        // println!("i = {}, ct.b = {}", i, ct.b);
    }
    // println!("m = {}, ct.b = {}, ct->a = {:?}", m, ct.b, ct.a);
}
pub fn decrypt(sk: &SecretKey, ct: &CipherText) -> i32 {
    let mut r = ct.b;
    // dbg!(r);
    for i in 0..n {
        r -= ct.a[i] * sk[i];
        // dbg!(r);
    }
    r = ((r % qi) + qi + qi/8) % qi;
    // dbg!(r,qi);
    // dbg!(r, qi, 4*r/qi);
    4 * r / qi
}

pub fn switching_key_gen<'b,'c>(res: &mut SwitchingKey, new_sk: ArrayView<'b,i32,Ix1>, old_sk: ArrayView<'c,i32,Ix1>) {
    // dbg!(&res[[0,0,0]]);
    for i in 0..N {
        for j in 0..KS_BASE {
            for k in 0..KS_EXP {
                // dbg!("switching key",i,j,k);
                let mut ct: CipherTextQ = Default::default();
                // ct.a.par_map_inplace(|x| *x = rand::thread_rng().gen());
                ct.a.par_map_inplace(|x| *x = Wrapping(i as i32) * Wrapping(j as i32) * Wrapping(k as i32));
                ct.b = // Wrapping(sample(&CHI2))
                    - Wrapping(old_sk[i]) * Wrapping(j as i32) * Wrapping(KS_TABLE[k])
                    + ct.a.dot::<Array<Wrapping<i32>,Ix1>>(&new_sk.mapv(Wrapping));
                    // + Zip::from(&ct.a).and(new_sk).fold(Wrapping(0), |acc, &x, &y| acc + x * Wrapping(y));
                    // + ct.a.par_iter().zip(new_sk.par_iter()).map(|(x,y)| x*Wrapping(y)).sum();
                res[[i,j,k]] = ct;
            }
        }
    }
    // dbg!(&res[[0,0,0]]);
    eprintln!("evaluation key 3 -> switching key generation complete");
}
pub fn key_switch(res: &mut CipherTextQ, ksk: &SwitchingKey, ct: &CipherTextQN) {
    for k in 0..n {
        res.a[k] = Wrapping(0);
        // println!("res.a[{}] = {}", k, res.a[k]);
    }
    res.b = ct.b;
    // dbg!(res.b);
    for i in 0..N {
        let mut a: UZmodQ = Wrapping(0) - Wrapping(ct.a[i].0 as u32);
        // println!("i = {}, a = {}", i, a);
        for j in 0..KS_EXP {
            let a0: UZmodQ = a % Wrapping(KS_BASE as u32);
            // print!("i = {}, j = {}, ksk[i][{}][j].a = [", i, j, a0.0);
            for k in 0..n {
                res.a[k] -= ksk[[i,a0.0 as usize,j]].a[k];
                // print!("{},", ksk[[i,a0.0 as usize,j]].a[k]);
            }
            // println!("]");
            res.b -= ksk[[i,a0.0 as usize,j]].b;
            a /= Wrapping(KS_BASE as u32);
        }
    }
}
pub fn mod_switch(ct: &mut CipherText, c: &CipherTextQ) {
    ct.a = c.a.mapv(round_q_Q);
    ct.b = round_q_Q(c.b);
}
pub fn round_q_Q(v: ZmodQ) -> i32 {
    (0.5 + (v.0 as f64) * (q as f64) / (Q as f64)).floor() as i32
}
