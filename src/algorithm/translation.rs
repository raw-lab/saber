//! Genetic code translation for nucleotide-to-protein conversion.
//!
//! v2.0 had two bugs:
//!   1. Stop codons mapped to `'X'` (which is the BLAST ambiguity character)
//!      and then `translate_frame` broke out of the loop at the first stop.
//!      A 6-frame search must keep going past internal stops.
//!   2. Reverse-complement frames called `translate_frame(reverse_comp, frame)`,
//!      which is correct only if you also report subject coordinates back in
//!      the original strand. The reporting is now handled by the helpers below.
//!
//! Stops are now `*`; downstream scoring uses BLOSUM `*` row/col like blastp.

use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TranslationTable {
    Standard,
    Vertebrate,
    Yeast,
    Mold,
    Invertebrate,
    Ciliate,
    Echinoderm,
    /// Bacterial, archaeal and plant plastid code (NCBI 11).
    Bacterial,
    AltYeast,
    /// Archaea use NCBI table 11, so this translates exactly like
    /// [`TranslationTable::Bacterial`]. (SABER 1.0/1.1 mapped `archaeal` to
    /// NCBI table 15, the Blepharisma nuclear code, which reads TAG as Q.)
    Archaeal,
    /// Blepharisma nuclear code (NCBI 15): TAG = Q.
    Blepharisma,
}

impl std::str::FromStr for TranslationTable {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "standard" | "1" => Ok(Self::Standard),
            "vertebrate" | "2" => Ok(Self::Vertebrate),
            "yeast" | "3" => Ok(Self::Yeast),
            "mold" | "4" => Ok(Self::Mold),
            "invertebrate" | "5" => Ok(Self::Invertebrate),
            "ciliate" | "6" => Ok(Self::Ciliate),
            "echinoderm" | "9" => Ok(Self::Echinoderm),
            "bacterial" | "11" => Ok(Self::Bacterial),
            "alt-yeast" | "altyeast" | "12" => Ok(Self::AltYeast),
            "archaeal" => Ok(Self::Archaeal),
            "blepharisma" | "15" => Ok(Self::Blepharisma),
            _ => Err(format!("Unknown genetic code: {}", s)),
        }
    }
}

/// Information about one translated reading frame.
#[derive(Debug, Clone)]
pub struct TranslatedFrame {
    /// Frame number (0..=5). 0..=2 forward, 3..=5 reverse-complement.
    pub frame: u8,
    /// Translated protein, with `*` for stop codons.
    pub protein: String,
    /// Length of the source nucleotide region that produced this frame.
    pub source_length: usize,
}

impl TranslatedFrame {
    /// Whether this frame is on the reverse-complement strand.
    pub fn is_reverse(&self) -> bool { self.frame >= 3 }
}

/// Translate one frame of a nucleotide sequence to protein. Stops emit `*`.
/// We do NOT break at stops — callers either split on `*` to get ORFs or pass
/// the whole sequence through Smith-Waterman with BLOSUM's `*` row/col.
pub fn translate_frame(seq: &[u8], table: TranslationTable, frame: usize) -> String {
    let codes = genetic_code(table);
    let mut out = String::with_capacity(seq.len() / 3);
    let mut i = frame;
    while i + 3 <= seq.len() {
        let codon = [
            seq[i].to_ascii_uppercase(),
            seq[i + 1].to_ascii_uppercase(),
            seq[i + 2].to_ascii_uppercase(),
        ];
        let aa = lookup_codon(&codon, codes);
        out.push(aa);
        i += 3;
    }
    out
}

/// Translate all six reading frames (3 forward + 3 reverse-complement).
pub fn translate_six_frames(seq: &[u8], table: TranslationTable) -> Vec<TranslatedFrame> {
    let mut out = Vec::with_capacity(6);
    for frame in 0..3usize {
        let p = translate_frame(seq, table, frame);
        if !p.is_empty() {
            out.push(TranslatedFrame { frame: frame as u8, protein: p, source_length: seq.len() });
        }
    }
    let rc = reverse_complement(seq);
    for frame in 0..3usize {
        let p = translate_frame(&rc, table, frame);
        if !p.is_empty() {
            out.push(TranslatedFrame { frame: (frame + 3) as u8, protein: p, source_length: seq.len() });
        }
    }
    out
}

/// Single-frame helper kept for backwards compatibility.
pub fn translate_dna_to_protein(dna: &[u8], table: TranslationTable) -> String {
    translate_frame(dna, table, 0)
}

/// Watson-Crick reverse complement. IUPAC ambiguity bases pass through as N.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' | b'U' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'N' => b'N',
            other => other,
        });
    }
    out
}

fn lookup_codon(codon: &[u8; 3], codes: &HashMap<&'static str, char>) -> char {
    if let Some(s) = std::str::from_utf8(codon).ok() {
        if let Some(&aa) = codes.get(s) {
            return aa;
        }
        // Codon contains ambiguity (e.g. N) — return X (unknown amino acid).
        // BLAST does this too.
    }
    'X'
}

fn genetic_code(table: TranslationTable) -> &'static HashMap<&'static str, char> {
    match table {
        TranslationTable::Standard => &STANDARD_GENETIC_CODE,
        TranslationTable::Vertebrate => &VERTEBRATE_MITOCHONDRIAL,
        TranslationTable::Yeast => &YEAST_MITOCHONDRIAL,
        TranslationTable::Mold => &MOLD_MITOCHONDRIAL,
        TranslationTable::Invertebrate => &INVERTEBRATE_MITOCHONDRIAL,
        TranslationTable::Ciliate => &CILIATE_NUCLEAR,
        TranslationTable::Echinoderm => &ECHINODERM_MITOCHONDRIAL,
        TranslationTable::Bacterial | TranslationTable::Archaeal => &BACTERIAL_CODE,
        TranslationTable::AltYeast => &ALTERNATIVE_YEAST,
        TranslationTable::Blepharisma => &BLEPHARISMA_NUCLEAR,
    }
}

lazy_static::lazy_static! {
    static ref STANDARD_GENETIC_CODE: HashMap<&'static str, char> = {
        let mut m = HashMap::new();
        m.insert("TTT", 'F'); m.insert("TTC", 'F'); m.insert("TTA", 'L'); m.insert("TTG", 'L');
        m.insert("TCT", 'S'); m.insert("TCC", 'S'); m.insert("TCA", 'S'); m.insert("TCG", 'S');
        m.insert("TAT", 'Y'); m.insert("TAC", 'Y'); m.insert("TAA", '*'); m.insert("TAG", '*');
        m.insert("TGT", 'C'); m.insert("TGC", 'C'); m.insert("TGA", '*'); m.insert("TGG", 'W');
        m.insert("CTT", 'L'); m.insert("CTC", 'L'); m.insert("CTA", 'L'); m.insert("CTG", 'L');
        m.insert("CCT", 'P'); m.insert("CCC", 'P'); m.insert("CCA", 'P'); m.insert("CCG", 'P');
        m.insert("CAT", 'H'); m.insert("CAC", 'H'); m.insert("CAA", 'Q'); m.insert("CAG", 'Q');
        m.insert("CGT", 'R'); m.insert("CGC", 'R'); m.insert("CGA", 'R'); m.insert("CGG", 'R');
        m.insert("ATT", 'I'); m.insert("ATC", 'I'); m.insert("ATA", 'I'); m.insert("ATG", 'M');
        m.insert("ACT", 'T'); m.insert("ACC", 'T'); m.insert("ACA", 'T'); m.insert("ACG", 'T');
        m.insert("AAT", 'N'); m.insert("AAC", 'N'); m.insert("AAA", 'K'); m.insert("AAG", 'K');
        m.insert("AGT", 'S'); m.insert("AGC", 'S'); m.insert("AGA", 'R'); m.insert("AGG", 'R');
        m.insert("GTT", 'V'); m.insert("GTC", 'V'); m.insert("GTA", 'V'); m.insert("GTG", 'V');
        m.insert("GCT", 'A'); m.insert("GCC", 'A'); m.insert("GCA", 'A'); m.insert("GCG", 'A');
        m.insert("GAT", 'D'); m.insert("GAC", 'D'); m.insert("GAA", 'E'); m.insert("GAG", 'E');
        m.insert("GGT", 'G'); m.insert("GGC", 'G'); m.insert("GGA", 'G'); m.insert("GGG", 'G');
        m
    };

    /// Bacterial (NCBI 11): TGA still stop, but alternative start codons (handled at start-call time, not here).
    static ref BACTERIAL_CODE: HashMap<&'static str, char> = STANDARD_GENETIC_CODE.clone();

    /// Vertebrate mitochondrial (NCBI 2): AGA, AGG = stop; ATA = M; TGA = W.
    static ref VERTEBRATE_MITOCHONDRIAL: HashMap<&'static str, char> = {
        let mut m = STANDARD_GENETIC_CODE.clone();
        m.insert("AGA", '*'); m.insert("AGG", '*');
        m.insert("ATA", 'M'); m.insert("TGA", 'W');
        m
    };

    /// Yeast mitochondrial (NCBI 3): CTN = T, ATA = M, TGA = W, CGA/CGC absent (rare).
    static ref YEAST_MITOCHONDRIAL: HashMap<&'static str, char> = {
        let mut m = STANDARD_GENETIC_CODE.clone();
        m.insert("CTT", 'T'); m.insert("CTC", 'T');
        m.insert("CTA", 'T'); m.insert("CTG", 'T');
        m.insert("ATA", 'M'); m.insert("TGA", 'W');
        m
    };

    /// Mold/Protozoan/Coelenterate mitochondrial (NCBI 4): TGA = W.
    static ref MOLD_MITOCHONDRIAL: HashMap<&'static str, char> = {
        let mut m = STANDARD_GENETIC_CODE.clone();
        m.insert("TGA", 'W');
        m
    };

    /// Invertebrate mitochondrial (NCBI 5): AGA/AGG = S, ATA = M, TGA = W.
    static ref INVERTEBRATE_MITOCHONDRIAL: HashMap<&'static str, char> = {
        let mut m = STANDARD_GENETIC_CODE.clone();
        m.insert("AGA", 'S'); m.insert("AGG", 'S');
        m.insert("ATA", 'M'); m.insert("TGA", 'W');
        m
    };

    /// Ciliate (NCBI 6): TAA and TAG = Q.
    static ref CILIATE_NUCLEAR: HashMap<&'static str, char> = {
        let mut m = STANDARD_GENETIC_CODE.clone();
        m.insert("TAA", 'Q'); m.insert("TAG", 'Q');
        m
    };

    /// Echinoderm/Flatworm mitochondrial (NCBI 9): AAA = N, AGA/AGG = S, TGA = W.
    static ref ECHINODERM_MITOCHONDRIAL: HashMap<&'static str, char> = {
        let mut m = STANDARD_GENETIC_CODE.clone();
        m.insert("AAA", 'N'); m.insert("AGA", 'S'); m.insert("AGG", 'S');
        m.insert("TGA", 'W');
        m
    };

    /// Alternative yeast (NCBI 12): CTG = S.
    static ref ALTERNATIVE_YEAST: HashMap<&'static str, char> = {
        let mut m = STANDARD_GENETIC_CODE.clone();
        m.insert("CTG", 'S');
        m
    };

    /// Blepharisma nuclear (NCBI 15): TAG = Q.
    static ref BLEPHARISMA_NUCLEAR: HashMap<&'static str, char> = {
        let mut m = STANDARD_GENETIC_CODE.clone();
        m.insert("TAG", 'Q');
        m
    };
}

/// Table-driven translator (64-entry codon table per genetic code) producing
/// protein *codes* directly. Codons containing anything but A/C/G/T/U give X.
#[derive(Clone)]
pub struct CodonTable {
    aa: [u8; 64],
    table: TranslationTable,
}

impl CodonTable {
    pub fn new(table: TranslationTable) -> Self {
        let codes = genetic_code(table);
        let mut aa = [crate::alphabet::AA_X; 64];
        let nts = [b'A', b'C', b'G', b'T'];
        for a in 0..4 {
            for b in 0..4 {
                for c in 0..4 {
                    let codon = [nts[a], nts[b], nts[c]];
                    let ch = lookup_codon(&codon, codes) as u8;
                    aa[a * 16 + b * 4 + c] = crate::alphabet::encode_protein_byte(ch);
                }
            }
        }
        Self { aa, table }
    }

    pub fn genetic_code(&self) -> TranslationTable {
        self.table
    }

    /// Translate frame `frame` (0..=5; 3..=5 on the reverse complement) of a
    /// nucleotide sequence given as *nucleotide codes* (0..=3, 4 = N) into
    /// protein codes appended to `out` (cleared first).
    pub fn translate_codes(&self, nt: &[u8], frame: usize, out: &mut Vec<u8>) {
        out.clear();
        let n = nt.len();
        if frame < 3 {
            let mut i = frame;
            while i + 3 <= n {
                out.push(self.codon(nt[i], nt[i + 1], nt[i + 2]));
                i += 3;
            }
        } else {
            // reverse complement read right-to-left: complement of code c is 3-c
            let off = frame - 3;
            let mut k = off;
            while k + 3 <= n {
                let a = nt[n - 1 - k];
                let b = nt[n - 2 - k];
                let c = nt[n - 3 - k];
                out.push(self.codon(comp(a), comp(b), comp(c)));
                k += 3;
            }
        }
    }

    #[inline]
    fn codon(&self, a: u8, b: u8, c: u8) -> u8 {
        if a > 3 || b > 3 || c > 3 {
            crate::alphabet::AA_X
        } else {
            self.aa[(a as usize) * 16 + (b as usize) * 4 + c as usize]
        }
    }
}

#[inline]
fn comp(c: u8) -> u8 {
    if c > 3 { c } else { 3 - c }
}

/// Signed frame label (+1..+3 forward, -1..-3 reverse) for frame 0..=5.
pub fn frame_label(frame: usize) -> i8 {
    if frame < 3 { frame as i8 + 1 } else { -((frame - 3) as i8 + 1) }
}

/// Map 1-based inclusive protein coordinates of `frame` back to 1-based
/// nucleotide coordinates on the original strand. Minus-strand hits are
/// reported BLAST-style with start > end.
pub fn protein_to_nt(frame: usize, p_start: usize, p_end: usize, nt_len: usize) -> (usize, usize) {
    if frame < 3 {
        let s = frame + 3 * (p_start - 1) + 1;
        let e = frame + 3 * p_end;
        (s, e)
    } else {
        let off = frame - 3;
        let rs = off + 3 * (p_start - 1) + 1; // on the reverse complement
        let re = off + 3 * p_end;
        (nt_len - rs + 1, nt_len - re + 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codon_table_matches_hashmap_translation_on_all_frames() {
        let seq = b"ATGGCGTTTAAACCCGGGTTTAAATAGCGNATGCATGCTAGCTAGGATCCAAGCTTGAATTC";
        let codes = crate::alphabet::encode_nucleotide(seq);
        let ct = CodonTable::new(TranslationTable::Standard);
        let mut out = Vec::new();
        for f in translate_six_frames(seq, TranslationTable::Standard) {
            ct.translate_codes(&codes, f.frame as usize, &mut out);
            let want = crate::alphabet::EncodedProtein::from_ascii(f.protein.as_bytes());
            assert_eq!(out, want.as_slice(), "frame {}", f.frame);
        }
    }

    #[test]
    fn protein_to_nt_coordinates() {
        // forward frame +2: protein residue 1 = nt 2..4
        assert_eq!(protein_to_nt(1, 1, 1, 30), (2, 4));
        // reverse frame -1 (index 3): protein residue 1 = last codon, reported 30..28
        assert_eq!(protein_to_nt(3, 1, 1, 30), (30, 28));
        assert_eq!(protein_to_nt(4, 2, 3, 30), (26, 21));
        assert_eq!(frame_label(0), 1);
        assert_eq!(frame_label(5), -3);
    }

    #[test]
    fn reverse_complement_atcg() {
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT".to_vec());
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT".to_vec()); // palindrome
    }

    #[test]
    fn standard_code_translates_methionine() {
        let p = translate_frame(b"ATG", TranslationTable::Standard, 0);
        assert_eq!(p, "M");
    }

    #[test]
    fn stop_codons_become_asterisks_and_we_keep_reading() {
        // ATG (M) TAA (*) GCT (A) → "M*A" not "M"
        let p = translate_frame(b"ATGTAAGCT", TranslationTable::Standard, 0);
        assert_eq!(p, "M*A");
    }

    #[test]
    fn six_frames_yield_six_strings() {
        let frames = translate_six_frames(b"ATGCCCGGGTTTAAA", TranslationTable::Standard);
        assert_eq!(frames.len(), 6);
        assert!(frames.iter().any(|f| f.is_reverse()));
    }

    #[test]
    fn ambiguous_codon_yields_x() {
        let p = translate_frame(b"ATN", TranslationTable::Standard, 0);
        assert_eq!(p, "X");
    }

    #[test]
    fn mitochondrial_tga_is_tryptophan_not_stop() {
        let standard = translate_frame(b"TGA", TranslationTable::Standard, 0);
        let mito = translate_frame(b"TGA", TranslationTable::Vertebrate, 0);
        assert_eq!(standard, "*");
        assert_eq!(mito, "W");
    }

    #[test]
    fn genetic_code_parsing_accepts_numbers() {
        assert_eq!("11".parse::<TranslationTable>().unwrap(), TranslationTable::Bacterial);
        assert_eq!("standard".parse::<TranslationTable>().unwrap(), TranslationTable::Standard);
    }

    #[test]
    fn archaeal_uses_table_11_and_table_15_is_blepharisma() {
        // 1.0/1.1 mapped "archaeal" to NCBI 15 (TAG = Q); archaea use table 11.
        assert_eq!("archaeal".parse::<TranslationTable>().unwrap(), TranslationTable::Archaeal);
        assert_eq!("15".parse::<TranslationTable>().unwrap(), TranslationTable::Blepharisma);
        assert_eq!("blepharisma".parse::<TranslationTable>().unwrap(), TranslationTable::Blepharisma);
        let tag = b"ATGTAGTGG";
        assert_eq!(translate_frame(tag, TranslationTable::Archaeal, 0), translate_frame(tag, TranslationTable::Bacterial, 0));
        assert_eq!(translate_frame(tag, TranslationTable::Archaeal, 0), "M*W");
        assert_eq!(translate_frame(tag, TranslationTable::Blepharisma, 0), "MQW");
        let codes: Vec<u8> = tag.iter().map(|&b| match b { b'A' => 0, b'C' => 1, b'G' => 2, _ => 3 }).collect();
        let (mut a, mut b) = (Vec::new(), Vec::new());
        CodonTable::new(TranslationTable::Archaeal).translate_codes(&codes, 0, &mut a);
        CodonTable::new(TranslationTable::Bacterial).translate_codes(&codes, 0, &mut b);
        assert_eq!(a, b);
    }
}
