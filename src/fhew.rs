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

// #[derive(Clone,Debug)]
// struct CtModQ(Array1<CtModQ1>); // [CtModQ1; K2];
// impl Default for CtModQ {
//     fn default() -> Self {
//         CtModQ(Array::default(K2))
//     }
// }
pub type CtModQ = Array<ZmodQ,Ix3>;
pub fn CtModQ() -> CtModQ {
    Array::default((K2,2,N))
}

// #[derive(Clone,Debug)]
// struct CtModQ1(Array1<RingModQ>); // [RingModQ; 2];
// impl Default for CtModQ1 {
//     fn default() -> Self {
//         CtModQ1(Array::default(2))
//     }
// }
pub type CtModQ1 = Array<ZmodQ,Ix2>;
pub fn CtModQ1() -> CtModQ1 {
    Array::default((2,N))
}

// #[derive(Clone,Debug)]
// struct DctModQ(Array1<DctModQ1>); // [DctModQ1; K2];
// impl Default for DctModQ {
//     fn default() -> Self {
//         DctModQ(Array::default(K2))
//     }
// }
pub type DctModQ = Array<ZmodQ,Ix3>;
pub fn DctModQ() -> DctModQ {
    Array::default((K2,K2,N))
}

// #[derive(Clone,Debug)]
// struct DctModQ1(Array1<RingModQ>); // [RingModQ; K2];
// impl Default for DctModQ1 {
//     fn default() -> Self {
//         DctModQ1(Array::default(K2))
//     }
// }
pub type DctModQ1 = Array<ZmodQ,Ix2>;
pub fn DctModQ1() -> DctModQ1 {
    Array::default((K2,N))
}

// #[derive(Clone,Debug)]
// struct DctFFT(Array1<DctFFT1>); // [DctFFT1; K2];
// impl Default for DctFFT {
//     fn default() -> Self {
//         DctFFT(Array::default(K2))
//     }
// }
pub type DctFFT = Array<c64,Ix3>;
pub fn DctFFT() -> DctFFT {
    Array::default((2,K2,N2))
}

// #[derive(Clone,Debug)]
// struct DctFFT1(Array1<RingFFT>); // [RingFFT; K2];
// impl Default for DctFFT1 {
//     fn default() -> Self {
//         DctFFT1(Array::default(K2))
//     }
// }
pub type DctFFT1 = Array<c64,Ix2>;
pub fn DctFFT1() -> DctFFT1 {
    Array::default((K2,N2))
}

// #[derive(Clone,Debug)]
// pub struct CtFFT(Array1<CtFFT1>); // [CtFFT1; K2];
// impl Default for CtFFT {
//     fn default() -> Self {
//         CtFFT(Array::default(K2))
//     }
// }
pub type CtFFT = Array<c64,Ix3>;
pub fn CtFFT() -> CtFFT {
    Array::default((K2,2,N2))
}

// #[derive(Clone,Debug)]
// struct CtFFT1(Array1<RingFFT>); // [RingFFT; 2];
// impl Default for CtFFT1 {
//     fn default() -> Self {
//         CtFFT1(Array::default(2))
//     }
// }
pub type CtFFT1 = Array<c64,Ix2>;
pub fn CtFFT1() -> CtFFT1 {
    Array::default((2,N2))
}

// #[derive(Clone,Debug)]
// pub struct BootstrappingKey(Array3<CtFFT>); // [[[CtFFT; BS_EXP]; BS_BASE]; n];
// impl Default for BootstrappingKey {
//     fn default() -> Self {
//         BootstrappingKey(Array::default((n,BS_BASE,BS_EXP)))
//     }
// }
pub type BootstrappingKey = Array<c64,Ix6>;
pub fn BootstrappingKey() -> BootstrappingKey {
    let k = Array::default((n,BS_BASE,BS_EXP,K2,2,N2));
    eprintln!("evaluation key 1 -> bootstrapping key zeroing complete");
    k
}

#[derive(Clone,Debug)]
pub struct EvalKey {
    pub bskey: BootstrappingKey,
    pub kskey: lwe::SwitchingKey
}
impl Default for EvalKey {
    fn default() -> Self {
        EvalKey {
            bskey: BootstrappingKey(),
            kskey: SwitchingKey()
        }
    }
}
pub fn setup(ffto: &mut FFT, t_test_msb: &mut RingFFT) {
    fft_setup(ffto);
    let mut tmsb: RingModQ = RingModQ();
    tmsb[0] = Wrapping(-1);
    for i in 1..N {
        tmsb[i] = Wrapping(1);
    }
    fft_forward(ffto, tmsb.view(), t_test_msb.view_mut());
    // for i in 0..N2 { // no issue here
    //     println!("t_test_msb[{}] = {}, {}", i, t_test_msb[i].re, t_test_msb[i].im);
    // }
}
fn fhew_encrypt<'a>(mut ct: ArrayViewMut<'a,c64,Ix3>, sk_fft: &RingFFT, m: i32, ffto: &mut FFT) {
    let mut ai: RingFFT = RingFFT();
    let mut res: CtModQ = CtModQ();
    let qi = q as i32;
    let Ni = N as i32;
    let mut mm: i32 = (((m % qi) + qi) % qi) * (2*Ni/qi);
    let old_mm = mm;
    let mut sign: i32 = 1;
    let old_sign = sign;
    if mm >= Ni {
        mm -= Ni;
        sign = -1;
    }
    // println!("m = {}, old_mm = {}, new_mm = {}, old_sign = {}, new_sign = {}", m, old_mm, mm, old_sign, sign);

    // println!("stage 1");
    for i in 0..K2 {
        for k in 0..(Ni as usize) {
            // res[[i,0,k]] = rand::thread_rng().gen();
            res[[i,0,k]] = Wrapping(k as i32);
            // println!("res[[{},0,{}]] = {}", i, k, res[[i,0,k]]);
        }
        fft_forward(ffto, res.slice::<Ix1>(s![i,0,..]), ai.view_mut()); 
        for k in 0..(N2 as usize) { // minute differences, mostly because of precision
            // print!("i = {}, k = {}, ai_old[k] = ({}, {}), , sk_fft[k] = ({}, {}), ", i, k, ai[k].re, ai[k].im, sk_fft[k].re, sk_fft[k].im);
            ai[k] = ai[k] * sk_fft[k];
            // println!("ai_new[k] = ({}, {})", ai[k].re, ai[k].im);
        }
        fft_backward(ffto, ai.view(), res.slice_mut::<Ix1>(s![i,1,..]));
        // for k in 0..(Ni as usize) {
        //     res[[i,1,k]] += Wrapping(sample(&CHI1) as i32);
        // }
        
        // for k in 0..(Ni as usize) { // res is clean
        //     println!("i = {}, k = {}, res[i][0][k] = {}, res[i][1][k] = {}", i, k, 
        //         res[[i,0,k]], res[[i,1,k]]);
        // }
    }
    // println!("stage 2");
    for i in 0..K {
        res[[2*i,0,mm as usize]] += Wrapping(sign) * VGPRIME[i];
        res[[2*i+1,1,mm as usize]] += Wrapping(sign) * VGPRIME[i];
        // println!("i = {}, mm = {}, res[2*i+0][0][{}] = {}, res[2*i+1][1][{}] = {}", i, mm, mm as usize, res[[2*i,0,mm as usize]], mm as usize, res[[2*i+1,1,mm as usize]]);
    }
    // println!("stage 3");
    for i in 0..K2 {
        for j in 0..2 {
            fft_forward(ffto, res.slice::<Ix1>(s![i,j,..]), ct.slice_mut::<Ix1>(s![i,j,..]));
            // for k in 0..N2 {
            //     ct[[i,j,k]] = c64::new(ct[[i,j,k]].re.round(), ct[[i,j,k]].im.round());
            // }
        }
    }
    // println!("stage 4");
}
pub fn key_gen(ek: &mut EvalKey, lwe_sk: &lwe::SecretKey, ffto: &mut FFT) {
    let mut fhew_sk: SecretKeyN = SecretKeyN();
    key_gen_N(&mut fhew_sk);
    switching_key_gen(&mut ek.kskey, lwe_sk.view(), fhew_sk.view());
    
    let mut fhew_sk_fft: RingFFT = RingFFT();
    let fhew_rmq: RingModQ = fhew_sk.mapv(Wrapping);
    fft_forward(ffto, fhew_rmq.view(), fhew_sk_fft.view_mut());
    // for i in 0..N2 {
    //     println!("fhew_sk_fft[{}] = ({}, {})", i, fhew_sk_fft[i].re, fhew_sk_fft[i].im);
    // }
    for i in 0..n {
        for j in 1..BS_BASE {
            for k in 0..BS_EXP {

                // dbg!("bootstrapping key",i,j,k);

                // println!("i0 = {}, i1 = {}, i2 = {}", i, j, k);
                fhew_encrypt(ek.bskey.slice_mut::<Ix3>(s![i,j,k,..,..,..]), &fhew_sk_fft, lwe_sk[i]*(j as i32)*(BS_TABLE[k] as i32), ffto); // wrong at 2,4,0,0,0,0
                
                // if i == 2 && j == 4 && k == 0 {
                //     println!("lwe_sk[i] = {}, j = {}, BS_TABLE[k] = {}, product = {}", 
                //         lwe_sk[i], j as i32, BS_TABLE[k] as i32, lwe_sk[i]*(j as i32)*(BS_TABLE[k] as i32));
                // }
            }
        }
    }

    // for i in 0..n { // bootstrapping key clear, except for precision issues
    //     for j in 1..BS_BASE {
    //         for k in 0..BS_EXP {
    //             for x in 0..K2 {
    //                 for y in 0..2 {
    //                     for z in 0..N2 {
    //                         let bskey = &ek.bskey[[i,j,k,x,y,z]];
    //                         println!("bskey[{}][{}][{}][{}][{}][{}] = ({}, {})", i, j, k, x, y, z, bskey.re, bskey.im);
    // }}}}}}
    // for i in 0..N { // switching key is clear
    //     for j in 0..KS_BASE {
    //         for k in 0..KS_EXP {
    //             let kskey = &ek.kskey[[i,j,k]];
    //             print!("i = {}, j = {}, k = {}, kskey.b = {}, kskey.a = [", i, j, k, kskey.b);
    //             for t in 0..n {print!("{}, ", kskey.a[t]);}
    //             print!("]\n");
    // }}}

    eprintln!("evaluation key 4 -> bootstrapping key generation complete");
}
fn add_to_acc<'a,'b>(mut acc: ArrayViewMut<'a,c64,Ix2>, c: ArrayView<'b,c64,Ix3>, ffto: &mut FFT) {
    let mut ct: CtModQ1 = CtModQ1();
    let mut dct: DctModQ1 = DctModQ1();
    let mut dct_fft: DctFFT1 = DctFFT1();
    for j in 0..2 {
        fft_backward(ffto, acc.slice(s![j,..]), ct.slice_mut(s![j,..])); 
        // for k in 0..N {
        //     print!("j = {}, k = {}, ct[j][k] = {}\t", j, k, ct[[j,k]]);
        //     if k < N2 {
        //         print!("acc[j][k] = ({}, {})", acc[[j,k]].re, acc[[j,k]].im);
        //     }
        //     println!("");
        // }
    }
    for j in 0..2 {
        for k in 0..N {
            let mut t: ZmodQ = ct[[j,k]] * V_INVERSE;
            // println!("j = {}, k = {}, t = {}, ct[[j,k]] = {}, V_INVERSE = {}", j, k, t, ct[[j,k]], V_INVERSE);
            for l in 0..K {
                let t_temp = t.clone();
                let r: ZmodQ = Wrapping((t.0 << G_BITS_32[l]) >> G_BITS_32[l]);
                t = Wrapping((t - r).0 >> G_BITS[l]);
                dct[[j+2*l,k]] = r;
                // println!("j = {}, k = {}, l = {}, ct[j][k] = {}, dct[j+2*l][k] = {}, t_old = {}, r = {}, t_new = {}", 
                //     j, k, l, ct[[j,k]], dct[[j+2*l,k]], t_temp, r, t);
            }
        }
    }
    
    for j in 0..K2 {
        fft_forward(ffto, dct.slice(s![j,..]), dct_fft.slice_mut(s![j,..]));
    }
    for j in 0..2 {
        for k in 0..N2 {
            acc[[j,k]] = c64::default();
            for l in 0..K2 {
                acc[[j,k]] += dct_fft[[l,k]] * c[[l,j,k]]; // here
                // println!("l = {}, j = {}, k = {}, dct_fft[l][k] = ({:.0}, {:.0}), c[l][j][k] = ({:.0}, {:.0})", l, j, k, 
                //     dct_fft[[l,k]].re, dct_fft[[l,k]].im, c[[l,j,k]].re, c[[l,j,k]].im);
            }
        }
        // if j == 0 {
        //     for k in 0..N {
        //         print!("j = {}, k = {}, ct[j][k] = {}\t", j, k, ct[[j,k]]);
        //         if k < N2 {
        //             print!("acc[j][k] = {}", acc[[j,k]]);
        //         }
        //         println!("");
        //     }
        // }
    }
}
fn initialize_acc(acc: &mut CtFFT1, m: i32, ffto: &mut FFT) {
    let mut res: CtModQ1 = CtModQ1();
    let qi = q as i32;
    let Ni = N as i32;
    let mut mm = (((m % qi) + qi) % qi) * (2*Ni/qi);
    let old_mm = mm;
    let mut sign: i32 = 1;
    let old_sign = sign;
    if mm >= Ni {
        mm -= Ni;
        sign = -1;
    }
    // println!("m = {}, old_mm = {}, new_mm = {}, old_sign = {}, new_sign = {}", m, old_mm, mm, old_sign, sign);

    for j in 0..2 {
        for k in 0..N {
            res[[j,k]] = Wrapping(0);
        }
    }
    res[[1,mm as usize]] += Wrapping(sign) * VGPRIME[0];

    // for j in 0..2 {
    //     for k in 0..N {
    //         println!("j = {}, k = {}, res[j][k] = {}", j, k, res[[j,k]]);
    //     }
    // }
    for j in 0..2 { // this is the culprit
        fft_forward(ffto, res.slice(s![j,..]), acc.slice_mut(s![j,..]));
    }
    // for j in 0..2 {
    //     for k in 0..N2 {
    //         println!("j = {}, k = {}, acc[j][k] = ({}, {})", j, k, acc[[j,k]].re, acc[[j,k]].im);
    //     }
    // }
}
fn member_test(t: &RingFFT, c: &CtFFT1, ffto: &mut FFT) -> CipherTextQN {
    let mut temp: RingFFT = RingFFT();
    let mut temp_mod_q: RingModQ = RingModQ();
    let mut ct: CipherTextQN = Default::default();
    // for i in 0..N2 {
    //         dbg!(i,t[[i]]);
    // }

    // for i in 0..2 {
    //     for j in 0..N2 {
    //         dbg!(i,j,c[[i,j]]);
    //     }
    // }    
    for i in 0..N2 {
        temp[i] = (c[[0,i]] * t[i]).conj();
        // print!("i = {}, temp[i] = {}, c[1][j] = {}, t[i] = {}\n", i, temp[i], c[[0,i]], t[i]);
    }
    // for i in 0..N2 {
    //         dbg!(i,temp[[i]]);
    // }
    // dbg!(&temp);
    fft_backward(ffto, temp.view(), temp_mod_q.view_mut());
    // for i in 0..N {
    //         dbg!(i,temp_mod_q[[i]]);
    // }
    for i in 0..N {
        ct.a[i] = temp_mod_q[i];
    }
    for i in 0..N2 {
        temp[i] = c[[1,i]] * t[i];
        // dbg!(i,temp[i],c[[1,i]],t[i]);
    }
    // for i in 0..N2 {
    //         dbg!(i,temp[[i]]);
    // }
    fft_backward(ffto, temp.view(), temp_mod_q.view_mut());
    // for i in 0..N {
    //     dbg!(i,temp_mod_q[[i]]);
    // }
    ct.b = V + temp_mod_q[0];
    // dbg!(ct.b,V,temp_mod_q[0]);
    ct
}
pub fn hom_gate(
    res: &mut CipherText, 
    gate: BinGate, 
    ek: &EvalKey, 
    ct1: &CipherText, 
    ct2: &CipherText, 
    ffto: &mut FFT, 
    t_test_msb: &RingFFT
) {
    let mut e12: CipherText = Default::default();
    let qi = q as i32;
    // for i in 0..n {
    //     if ((ct1.a[i] - ct2.a[i]) % qi) != 0 && ((ct1.a[i] + ct2.a[i]) % qi) != 0 {
    //         break;
    //     }
    //     if i == n - 1 {
    //         panic!("ERROR: Please only use independant ciphertexts as inputs.");
    //     }
    // }
    for i in 0..n {
        e12.a[i] = (2*qi - (ct1.a[i] + ct2.a[i])) % qi;
    }
    e12.b = (GATE_CONST[gate as usize] as i32) - ((ct1.b + ct2.b) % qi);

    // println!("e12.b = {}, e12.a = {}", e12.b, e12.a);
    // println!("ct1.b = {}, ct1.a = {}", ct1.b, ct1.a);
    // println!("ct2.b = {}, ct2.a = {}", ct2.b, ct2.a);

    let mut acc: CtFFT1 = CtFFT1();
    initialize_acc(&mut acc, (e12.b + qi/4) % qi, ffto); // all clear
    // for i in 0..2 {
    //     for j in 0..N2 {
    //         println!("acc[{}][{}] = {}, {}", i, j, acc[[i,j]].re, acc[[i,j]].im);
    //     }
    // }
    for i in 0..n {
        let mut a = (qi - e12.a[i] % qi) % qi;
        for k in 0..BS_EXP {
            let a0 = a % (BS_BASE as i32);
            // println!("i = {}, a0 = {}, k = {}", i, a0, k); // all clear
            if a0 != 0 {
                add_to_acc(acc.view_mut(), ek.bskey.slice(s![i,a0 as usize,k,..,..,..]), ffto); // trouble here
            }
            a /= BS_BASE as i32;
        }
    }
    // for i in 0..2 {
    //     for j in 0..N2 {
    //         println!("acc[{}][{}] = ({}, {})", i, j, acc[[i,j]].re, acc[[i,j]].im);
    //     }
    // }
    let e_qn: CipherTextQN = member_test(t_test_msb, &acc, ffto);
    // println!("e_qn.b = {}, e_qn.a = {}", e_qn.b, e_qn.a);
    let mut e_q: CipherTextQ = Default::default();
    key_switch(&mut e_q, &ek.kskey, &e_qn);
    // print!("e_q.b = {}, e_q.a = [", e_q.b);
    // for t in 0..n {println!("{}, ", e_q.a[t]);}
    // print!("]\n");
    mod_switch(res, &e_q);
    // dbg!(res);
}
pub fn hom_nand(
    res: &mut CipherText, 
    ek: &EvalKey, 
    ct1: &CipherText, 
    ct2: &CipherText, 
    ffto: &mut FFT, 
    t_test_msb: &RingFFT) {
    hom_gate(res, NAND, ek, ct1, ct2, ffto, t_test_msb)
}
pub fn hom_not(res: &mut CipherText, ct: &CipherText) {
    let qi = q as i32;
    for i in 0..n {
        res.a[i] = (qi - ct.a[i]) % qi;
    }
    res.b = (5 * qi / 4 - ct.b) % qi;
}

pub fn fwrite_ek(ek: &EvalKey, f: &mut File) {
    for i in 0..n {
        for j in 1..BS_BASE {
            for k in 0..BS_EXP {
                for l in 0..K2 {
                    for m in 0..2 {
                        for n_ in 0..N2 {
                            writeln!(*f, "{}", ek.bskey[[i,j,k,l,m,n_]]).unwrap();
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
                    writeln!(*f, "{}", ek.kskey[[i,j,k]].a[l]).unwrap();
                }
                writeln!(*f, "{}", ek.kskey[[i,j,k]].b).unwrap();
            }
        }
    }
    // f.sync_all().unwrap();
}
pub fn fread_ek(f: &File) -> EvalKey {
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
                            ek.bskey[[i,j,k,l,m,n_]] = c64::from_str(&s).unwrap();
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
                    ek.kskey[[i,j,k]].a[l] = Wrapping(i32::from_str(&s).unwrap());
                }
                b.read_line(&mut s).unwrap();
                ek.kskey[[i,j,k]].b = Wrapping(i32::from_str(&s).unwrap());
            }
        }
    }
    // f.sync_all().unwrap();
    ek
}