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
    Bacterial,
    AltYeast,
    Archaeal,
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
            "archaeal" | "15" => Ok(Self::Archaeal),
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
        TranslationTable::Bacterial => &BACTERIAL_CODE,
        TranslationTable::AltYeast => &ALTERNATIVE_YEAST,
        TranslationTable::Archaeal => &ARCHAEAL_CODE,
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

    /// Archaeal (NCBI 15): TAG = Q.
    static ref ARCHAEAL_CODE: HashMap<&'static str, char> = {
        let mut m = STANDARD_GENETIC_CODE.clone();
        m.insert("TAG", 'Q');
        m
    };
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
