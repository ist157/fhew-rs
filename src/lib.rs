#![allow(dead_code)]

use fftw::{
    array::AlignedVec,
    plan::*,
    types::*
};
use ndarray::*;
use rand::Rng;
use std::num::Wrapping;

pub const n: usize = 10;

pub const N: usize = 1024;
pub const N2: usize = N/2+1;

pub const K: usize = 3; // K
pub const K2: usize = 6;

pub const Q: usize = 1 << 32; // Q
pub const q: usize = 512;
pub const q2: usize = q/2;

pub type ZmodQ = Wrapping<i32>;
type UZmodQ = Wrapping<u32>;
const V: ZmodQ = Wrapping((1 << 29) + 1);
const V_INVERSE: ZmodQ = Wrapping(-536870911); // 3758096385; 1/V mod Q

const VGPRIME: [ZmodQ; 3] = [Wrapping(V.0), Wrapping(V.0 << 11), Wrapping(V.0 << 22)]; // [V, V<<11, V<<22];
const G_BITS: [isize; 3] = [11, 11, 10];
const G_BITS_32: [isize; 3] = [21, 21, 22];

pub const KS_BASE: usize = 25;
pub const KS_EXP: usize = 7;
pub const KS_TABLE: [i32; 7] = [
    1,
    25,
    25*25,
    25*25*25,
    25*25*25*25,
    25*25*25*25*25,
    25*25*25*25*25*25
];

pub const BS_BASE: usize = 23;
pub const BS_EXP: usize = 2;
pub const BS_TABLE: [usize; 2] = [1, 23];

// #[derive(Clone,Debug)]
// pub struct RingModQ(pub Array1<ZmodQ>); // [ZmodQ; N];
// impl Default for RingModQ {
//     fn default() -> RingModQ {
//         RingModQ(Array::default(N))
//     }
// }
pub type RingModQ = Array<ZmodQ,Ix1>;
pub fn RingModQ() -> RingModQ {
    Array::default(N)
}
// #[derive(Clone,Debug)]
// pub struct RingFFT(pub Array1<c64>); // [c64; N2];
// impl Default for RingFFT {
//     fn default() -> RingFFT {
//         RingFFT(Array::default(N2))
//     }
// }
pub type RingFFT = Array<c64,Ix1>;
pub fn RingFFT() -> RingFFT {
    Array::default(N2)
}

#[derive(Copy,Clone,Debug)]
pub enum BinGate { OR, AND, NOR, NAND }
const GATE_CONST: [usize; 4] = [15*q/8, 9*q/8, 11*q/8, 13*q/8];

pub mod lwe;
pub mod fhew;

struct Distrib {
    std_dev: f64,
    max: i32,
    offset: i32,
    table: &'static [f64]
}

fn sample(chi: &Distrib) -> i32 {
    // dbg!(chi.std_dev);
    if chi.max != 0 { // CHI1, CHI_BINARY
        // println!("path 1");
        let r: f64 = rand::thread_rng().gen();
        for i in 0..chi.max {
            if r <= chi.table[i as usize] {
                return i - chi.offset;
            }
        }
        panic!("Sampling Error: distribution table ending before (double) 1.0");
    }
    let mut r: f64;
    let s = chi.std_dev;
    if s < 500.0 { // CHI3
        // println!("path 2");
        let mut x: i32;
        let maxx= (s*8.0).ceil() as i32;
        loop {
            x = rand::thread_rng().gen::<i32>() % (2*maxx + 1) - maxx;
            r = rand::thread_rng().gen();
            // println!("x = {}, y = {}, z = {}", r, x, s);
            if r < (- (x*x) as f64 / (2.0*s*s)).exp() {
                return x;
            }
        }
    } else { // CHI2
        // println!("path 3");
        let mut x: f64;
        loop {
            x = 16.0 * rand::thread_rng().gen::<f64>() - 8.0;
            r = rand::thread_rng().gen();
            // println!("r = {}\tx = {}\ts = {}", r, x, s);
            if r < (- x*x / 2.0).exp() {
                return (0.5 + x*s).floor() as i32;
            }
        }
    }
}

const CHI1_TABLE: [f64; 23] = [
    1.12011750313263e-14, 2.38717233762211e-12, 3.04966971020178e-10,
    2.34394541319773e-8, 1.08538196465647e-6, 0.0000303513786306856,
    0.000514575939439740, 0.00532464699317562, 0.0340111330223921,
    0.136723892128727, 0.357520614142345, 0.642479385857655,
    0.863276107871273, 0.965988866977608, 0.994675353006824,
    0.999485424060560, 0.999969648621369, 0.999998914618035,
    0.999999976560546, 0.999999999695033, 0.999999999997613,
    0.999999999999989, 1.00000000000000
];

const CHI1: Distrib = Distrib {
    std_dev: 1.4,
    max: 23,
    offset: 11,
    table: &CHI1_TABLE
};

const NO_TABLE: [f64; 1] = [1.0];

const CHI2: Distrib = Distrib {
    std_dev: (1 << 17) as f64,
    max: 0,
    offset: 0,
    table: &NO_TABLE
};

const CHI3: Distrib = Distrib {
    std_dev: 6.0,
    max: 0,
    offset: 0,
    table: &NO_TABLE
};

const BINARY_TABLE: [f64; 3] = [
    0.25,
    0.75,
    1.0
];

const CHI_BINARY: Distrib = Distrib {
    std_dev: 0.0,
    max: 3,
    offset: 1,
    table: &BINARY_TABLE
};

pub struct FFT {
    pub in_: AlignedVec<f64>,
    pub out: AlignedVec<c64>,
    pub plan_fft_forw: fftw::plan::R2CPlan64,
    pub plan_fft_back: fftw::plan::C2RPlan64
}
impl Default for FFT {
    fn default() -> Self {
        FFT {
            in_: AlignedVec::new(N*2),
            out: AlignedVec::new(N+1),
            plan_fft_forw: R2CPlan64::aligned(&[N*2], Flag::PATIENT).unwrap(),
            plan_fft_back: C2RPlan64::aligned(&[N*2], Flag::PATIENT | Flag::PRESERVEINPUT).unwrap()
        }
    }
}
pub fn fft_setup(ffto: &mut FFT) {
    *ffto = Default::default();
}
pub fn fft_forward<'a,'b>(ffto: &mut FFT, val: ArrayView<'a,ZmodQ,Ix1>, mut res: ArrayViewMut<'b,c64,Ix1>) {
    for k in 0..N {
        ffto.in_[k] = val[k].0 as f64;
        ffto.in_[k+N] = 0.0;
    }
    ffto.plan_fft_forw.r2c(&mut ffto.in_, &mut ffto.out).unwrap();
    for k in 0..(N2-1) {
        res[k] = ffto.out[2*k+1];
        // res[k] = ffto.out[k];
    }
}
pub fn fft_backward<'a,'b>(ffto: &mut FFT, val: ArrayView<'a,c64,Ix1>, mut res: ArrayViewMut<'b,ZmodQ,Ix1>) {
    for k in 0..N2 {
        ffto.out[2*k] = c64::new(0.0,0.0);
        if k < N2 - 1 {
            ffto.out[2*k+1] = val[k]/c64::new(N as f64,0.0);
        }
        // ffto.out[k] = val[k]; // /c64::new(N as f64,0.0);
    }
    ffto.plan_fft_back.c2r(&mut ffto.out, &mut ffto.in_).unwrap();
    for k in 0..N {
        // let max: i64 = i32::MAX as i64;
        // let min: i64 = i32::MIN as i64;
        // let div: i64 = max - min + 1;
        // let mut t = ffto.in_[k].round() as i64 % div;
        // t = if t > max { t - div } else if t < min { t + div } else { t };
        // res[k] = Wrapping(t as i32);

        res[k] = Wrapping(ffto.in_[k].round().rem_euclid(2f64.powi(32)) as u32 as i32);
    }
}
