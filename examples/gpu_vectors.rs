//! Test vectors for the CUDA kernel emulation harness (tools/gpu_emulate).
//! Expected scores come from the naive i64 reference (`sw_reference`) and are
//! cross-checked against the production CPU kernel (`score_many`) here.
use saber_rs::algorithm::reference::sw_reference;
use saber_rs::algorithm::sw_simd::{score_many, KernelCounters, SimdScratch};
use saber_rs::alphabet::ProteinCodes;
use saber_rs::scoring::stats::GapCosts;
use saber_rs::scoring::{ScoreMatrix, ScoringSystem};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn emit(q: &[u8], bias: Option<&[i8]>, wins: &[Vec<u8>], m: &[[i32; 24]; 24], go: i32, ge: i32) {
    let mut out = String::new();
    out += &format!("CASE {} {} {} {}\n", q.len(), wins.len(), go, ge);
    out += &q.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(" ");
    out += "\n";
    let zeros = vec![0i8; q.len()];
    let b = bias.unwrap_or(&zeros);
    out += &b.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(" ");
    out += "\n";
    out += &m.iter().flat_map(|r| r.iter()).map(|c| c.to_string()).collect::<Vec<_>>().join(" ");
    out += "\n";
    for w in wins {
        out += &format!("{} {}\n", w.len(), w.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(" "));
    }
    let views: Vec<ProteinCodes<'_>> = wins.iter().map(|w| ProteinCodes::new_unchecked(w)).collect();
    let mut sc = SimdScratch::new();
    let mut kc = KernelCounters::default();
    let fast = score_many(ProteinCodes::new_unchecked(q), bias, &views, m, GapCosts { open: go, extend: ge }, &mut sc, &mut kc);
    let mut exp = Vec::new();
    for (i, w) in wins.iter().enumerate() {
        let (r, _, _) = sw_reference(q, w, m, bias, go as i64, ge as i64);
        assert_eq!(r, fast[i] as i64, "CPU kernel disagrees with reference");
        exp.push(r.to_string());
    }
    out += &exp.join(" ");
    out += "\n";
    print!("{out}");
}

fn main() {
    let mut r = Rng(0xC0DA);
    let blosum = ScoreMatrix::new(ScoringSystem::Blosum62).matrix;
    let pam = ScoreMatrix::new(ScoringSystem::Pam30).matrix;
    let nt = ScoreMatrix::nucleotide(2, -3).matrix;
    for case in 0..24 {
        let (m, alpha, go, ge) = match case % 4 {
            0 => (&blosum, 24u64, 11, 1),
            1 => (&pam, 20u64, 10, 1),
            2 => (&nt, 5u64, 7, 2),
            _ => (&blosum, 24u64, 13, 2),
        };
        let ql = 1 + r.below(300) as usize;
        let q: Vec<u8> = (0..ql).map(|_| r.below(alpha) as u8).collect();
        let bias: Option<Vec<i8>> = if case % 3 == 0 { Some((0..ql).map(|_| r.below(5) as i8 - 2).collect()) } else { None };
        let nw = 1 + r.below(40) as usize;
        let wins: Vec<Vec<u8>> = (0..nw)
            .map(|k| {
                if k % 5 == 0 {
                    // related window: mutated copy of a query segment
                    let lo = r.below(ql as u64) as usize;
                    let hi = lo + r.below((ql - lo) as u64 + 1) as usize;
                    q[lo..hi].iter().map(|&c| if r.below(4) == 0 { r.below(alpha) as u8 } else { c }).collect()
                } else {
                    let l = r.below(400) as usize;
                    (0..l).map(|_| r.below(alpha) as u8).collect()
                }
            })
            .collect();
        emit(&q, bias.as_deref(), &wins, m, go, ge);
    }
    // Long identical high-scoring sequences: scores far above i16::MAX
    // (v1.1's int16 kernel would have clamped these to 32767).
    let w = 17u8; // W
    let long: Vec<u8> = vec![w; 4000];
    emit(&long, None, &[long.clone(), long[..3500].to_vec()], &blosum, 11, 1);
    // Out-of-range codes are scored as X.
    emit(&[0, 1, 2, 200, 4], None, &[vec![0, 1, 2, 22, 4], vec![255, 255]], &blosum, 11, 1);
}
