//! Canonical residue encodings shared by indexing, scalar scoring, SIMD and CUDA.
//!
//! Every stage of SABER that touches residues uses exactly one of the encodings
//! defined here, and the types make the distinction between *ASCII* sequences
//! (what FASTA files contain) and *encoded* sequences (what the kernels consume)
//! explicit:
//!
//! * [`ProteinCodes`] / [`EncodedProtein`] — protein residues as indices into
//!   the 24-letter NCBI matrix alphabet `ARNDCQEGHILKMFPSTWYVBZX*`
//!   (codes `0..=23`). Unknown letters (`J`, `U`, `O`, digits, …) become `X`
//!   (22), never the stop symbol. `*` is the only input that becomes 23.
//! * [`NucleotideCodes`] — `A C G T/U → 0..3`, anything else `N` (4).
//!
//! Seed eligibility is a separate concept from the alignment encoding: only the
//! 20 standard amino acids (codes `0..20`) may participate in a seed, and the
//! seeding layer may additionally exclude masked (low-complexity) positions.
//! Before v1.2 the index stored ambiguous residues as `255` and the SIMD kernel
//! clamped that to 23, silently scoring every `X`/`B`/`Z` in the database as a
//! stop codon; that cannot happen with these types.

/// The 24-letter protein alignment alphabet, in NCBI matrix order.
pub const PROTEIN_ALPHABET: &[u8; 24] = b"ARNDCQEGHILKMFPSTWYVBZX*";
/// Number of letters in the protein alignment alphabet.
pub const PROTEIN_ALPHABET_SIZE: usize = 24;
/// Number of standard amino acids (codes `0..20`).
pub const N_STANDARD_AA: u8 = 20;
/// Code of `B` (Asx).
pub const AA_B: u8 = 20;
/// Code of `Z` (Glx).
pub const AA_Z: u8 = 21;
/// Code of `X` (unknown residue).
pub const AA_X: u8 = 22;
/// Code of `*` (translation stop).
pub const AA_STOP: u8 = 23;

/// Matrix index for an alignment code: codes outside the 24-letter alphabet
/// (which valid encodings never produce) are scored as X, everywhere.
#[inline(always)]
pub fn code_index(c: u8) -> usize {
    if c < 24 { c as usize } else { AA_X as usize }
}

/// Nucleotide alphabet `ACGTN`.
pub const NUCLEOTIDE_ALPHABET: &[u8; 5] = b"ACGTN";
/// Code of an ambiguous nucleotide.
pub const NT_N: u8 = 4;

const fn build_protein_table() -> [u8; 256] {
    let mut t = [AA_X; 256];
    let mut i = 0;
    while i < 24 {
        let c = PROTEIN_ALPHABET[i];
        t[c as usize] = i as u8;
        if c >= b'A' && c <= b'Z' {
            t[(c + 32) as usize] = i as u8;
        }
        i += 1;
    }
    t
}

const fn build_nucleotide_table() -> [u8; 256] {
    let mut t = [NT_N; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
}

/// ASCII → protein code lookup table (unknown → X).
pub static PROTEIN_ENCODE: [u8; 256] = build_protein_table();
/// ASCII → nucleotide code lookup table (unknown → N).
pub static NUCLEOTIDE_ENCODE: [u8; 256] = build_nucleotide_table();

/// Encode one ASCII protein letter. Unknown letters map to `X` (22).
#[inline]
pub fn encode_protein_byte(c: u8) -> u8 {
    PROTEIN_ENCODE[c as usize]
}

/// Decode a protein code back to its upper-case ASCII letter.
/// Out-of-range codes decode to `X`.
#[inline]
pub fn decode_protein(code: u8) -> u8 {
    if (code as usize) < PROTEIN_ALPHABET_SIZE {
        PROTEIN_ALPHABET[code as usize]
    } else {
        b'X'
    }
}

/// Encode one ASCII nucleotide. Anything but `ACGTU` (any case) maps to `N`.
#[inline]
pub fn encode_nucleotide_byte(c: u8) -> u8 {
    NUCLEOTIDE_ENCODE[c as usize]
}

/// True when `code` is one of the 20 standard amino acids, i.e. may be part of
/// a seed. Alignment still scores `B`, `Z`, `X` and `*` with their matrix rows.
#[inline]
pub fn is_seedable_aa(code: u8) -> bool {
    code < N_STANDARD_AA
}

/// Count of ASCII bytes that did not map to a canonical letter. Used to warn
/// users whose FASTA contains unexpected characters (gaps, digits, `J/U/O`).
pub fn count_noncanonical_protein(ascii: &[u8]) -> usize {
    ascii
        .iter()
        .filter(|&&c| {
            let code = PROTEIN_ENCODE[c as usize];
            code == AA_X && c != b'X' && c != b'x'
        })
        .count()
}

/// Borrowed view of a protein sequence in alignment encoding (codes 0..=23).
///
/// The only ways to build one are [`ProteinCodes::new`] (validated) and the
/// accessors on [`EncodedProtein`]; ASCII can never be passed by accident.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ProteinCodes<'a>(&'a [u8]);

impl<'a> ProteinCodes<'a> {
    /// Wrap already-encoded codes, verifying every code is `< 24`.
    pub fn new(codes: &'a [u8]) -> Result<Self, InvalidCode> {
        if let Some(pos) = codes.iter().position(|&c| c as usize >= PROTEIN_ALPHABET_SIZE) {
            return Err(InvalidCode { position: pos, code: codes[pos] });
        }
        Ok(Self(codes))
    }

    /// Wrap codes without validation. Callers must guarantee codes < 24;
    /// used on data that was validated when the index was opened.
    #[inline]
    pub fn new_unchecked(codes: &'a [u8]) -> Self {
        debug_assert!(codes.iter().all(|&c| (c as usize) < PROTEIN_ALPHABET_SIZE));
        Self(codes)
    }

    #[inline]
    pub fn as_slice(&self) -> &'a [u8] {
        self.0
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Sub-range view.
    #[inline]
    pub fn slice(&self, lo: usize, hi: usize) -> ProteinCodes<'a> {
        ProteinCodes(&self.0[lo..hi])
    }

    /// Decode to ASCII (upper case).
    pub fn to_ascii(&self) -> Vec<u8> {
        self.0.iter().map(|&c| decode_protein(c)).collect()
    }
}

/// Owned protein sequence in alignment encoding.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EncodedProtein(Vec<u8>);

impl EncodedProtein {
    /// Encode an ASCII sequence.
    pub fn from_ascii(ascii: &[u8]) -> Self {
        Self(ascii.iter().map(|&c| PROTEIN_ENCODE[c as usize]).collect())
    }

    /// Encode into an existing buffer (re-uses its allocation).
    pub fn encode_into(ascii: &[u8], out: &mut Vec<u8>) {
        out.clear();
        out.extend(ascii.iter().map(|&c| PROTEIN_ENCODE[c as usize]));
    }

    /// Take ownership of codes, validating them.
    pub fn from_codes(codes: Vec<u8>) -> Result<Self, InvalidCode> {
        ProteinCodes::new(&codes)?;
        Ok(Self(codes))
    }

    #[inline]
    pub fn codes(&self) -> ProteinCodes<'_> {
        ProteinCodes(&self.0)
    }

    #[inline]
    pub fn as_slice(&self) -> &[u8] {
        &self.0
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// Borrowed view of a nucleotide sequence in `ACGTN → 0..4` encoding.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NucleotideCodes<'a>(&'a [u8]);

impl<'a> NucleotideCodes<'a> {
    pub fn new(codes: &'a [u8]) -> Result<Self, InvalidCode> {
        if let Some(pos) = codes.iter().position(|&c| c > NT_N) {
            return Err(InvalidCode { position: pos, code: codes[pos] });
        }
        Ok(Self(codes))
    }
    #[inline]
    pub fn new_unchecked(codes: &'a [u8]) -> Self {
        debug_assert!(codes.iter().all(|&c| c <= NT_N));
        Self(codes)
    }
    #[inline]
    pub fn as_slice(&self) -> &'a [u8] {
        self.0
    }
    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// Encode an ASCII nucleotide sequence.
pub fn encode_nucleotide(ascii: &[u8]) -> Vec<u8> {
    ascii.iter().map(|&c| NUCLEOTIDE_ENCODE[c as usize]).collect()
}

/// Error returned when a byte is not a valid code for the requested encoding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvalidCode {
    pub position: usize,
    pub code: u8,
}

impl std::fmt::Display for InvalidCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "invalid residue code {} at position {} (expected an encoded residue, \
             not ASCII — did an ASCII sequence reach an encoded-only interface?)",
            self.code, self.position
        )
    }
}

impl std::error::Error for InvalidCode {}

// ---------------------------------------------------------------------------
// Reduced alphabets for seeding
// ---------------------------------------------------------------------------

/// A reduced amino-acid alphabet used only for *seed* matching. Alignment
/// always uses the full 24-letter encoding.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedAlphabet {
    /// Name stored in the index header.
    pub name: String,
    /// Number of reduced letters.
    pub size: u8,
    /// `map[code]` = reduced letter for codes 0..20, 255 for non-seedable.
    pub map: [u8; 24],
}

impl SeedAlphabet {
    /// Identity alphabet over the 20 standard amino acids.
    pub fn full() -> Self {
        let mut map = [255u8; 24];
        for (i, m) in map.iter_mut().enumerate().take(20) {
            *m = i as u8;
        }
        Self { name: "full".into(), size: 20, map }
    }

    /// Build from a space-separated group definition, e.g.
    /// `"A KR EDNQ C G H ILVM FYW P ST"` (Murphy et al. 2000, 10 groups — the
    /// reduction DIAMOND uses for most sensitivity modes).
    pub fn from_groups(name: &str, groups: &str) -> Result<Self, String> {
        let mut map = [255u8; 24];
        let mut n = 0u8;
        for g in groups.split_whitespace() {
            for c in g.bytes() {
                let code = encode_protein_byte(c);
                if code >= N_STANDARD_AA {
                    return Err(format!("reduced alphabet group '{g}' contains non-standard letter '{}'", c as char));
                }
                if map[code as usize] != 255 {
                    return Err(format!("letter '{}' appears in two groups", c as char));
                }
                map[code as usize] = n;
            }
            n += 1;
        }
        if map.iter().take(20).any(|&m| m == 255) {
            return Err("reduced alphabet must assign all 20 standard amino acids".into());
        }
        Ok(Self { name: name.into(), size: n, map })
    }

    /// Murphy, Wallqvist & Levy (2000) 10-letter alphabet, as used by DIAMOND.
    pub fn murphy10() -> Self {
        Self::from_groups("murphy10", "A KR EDNQ C G H ILVM FYW P ST").expect("static definition")
    }

    /// Parse a CLI/name specification: `full`, `murphy10`, or a custom group
    /// string prefixed by `custom:`.
    pub fn parse(spec: &str) -> Result<Self, String> {
        match spec.to_ascii_lowercase().as_str() {
            "full" | "full20" | "20" => Ok(Self::full()),
            "murphy10" | "murphy" | "10" => Ok(Self::murphy10()),
            s if s.starts_with("custom:") => Self::from_groups(spec, &spec[7..]),
            _ => Err(format!("unknown seed alphabet '{spec}' (expected full, murphy10, or custom:<groups>)")),
        }
    }

    /// Reduced letter for an alignment code, or `None` when not seedable.
    #[inline]
    pub fn reduce(&self, code: u8) -> Option<u8> {
        let m = self.map[(code as usize).min(23)];
        if m == 255 { None } else { Some(m) }
    }
}

/// A seed shape: positions of the "care" letters within a window. `"1111"` is
/// a contiguous 4-mer; `"11011"` is a spaced seed of weight 4 and span 5.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedShape {
    pub pattern: String,
    /// Offsets (from the window start) of the care positions.
    pub offsets: Vec<u16>,
    /// Window length (span).
    pub span: usize,
}

impl SeedShape {
    pub fn contiguous(k: usize) -> Self {
        Self::parse(&"1".repeat(k)).expect("contiguous shape")
    }

    pub fn parse(pattern: &str) -> Result<Self, String> {
        let p = pattern.trim();
        if p.is_empty() || !p.bytes().all(|b| b == b'0' || b == b'1') {
            return Err(format!("seed shape '{p}' must consist of 0/1"));
        }
        if !p.starts_with('1') || !p.ends_with('1') {
            return Err(format!("seed shape '{p}' must start and end with 1"));
        }
        let offsets: Vec<u16> = p.bytes().enumerate().filter(|(_, b)| *b == b'1').map(|(i, _)| i as u16).collect();
        if offsets.len() < 2 || offsets.len() > 13 {
            return Err(format!("seed shape '{p}' weight must be in [2, 13]"));
        }
        if p.len() > 64 {
            return Err(format!("seed shape '{p}' span must be <= 64"));
        }
        Ok(Self { pattern: p.to_string(), offsets, span: p.len() })
    }

    #[inline]
    pub fn weight(&self) -> usize {
        self.offsets.len()
    }

    #[inline]
    pub fn is_contiguous(&self) -> bool {
        self.span == self.offsets.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_round_trip_is_lossless_for_the_24_letters() {
        let enc = EncodedProtein::from_ascii(PROTEIN_ALPHABET);
        for (i, &c) in enc.as_slice().iter().enumerate() {
            assert_eq!(c as usize, i);
        }
        assert_eq!(enc.codes().to_ascii(), PROTEIN_ALPHABET.to_vec());
    }

    #[test]
    fn unknown_letters_become_x_not_stop() {
        for c in [b'J', b'U', b'O', b'1', b'-', b'.', b'j', 0u8, 255u8] {
            assert_eq!(encode_protein_byte(c), AA_X, "byte {c}");
        }
        assert_eq!(encode_protein_byte(b'*'), AA_STOP);
        assert_eq!(encode_protein_byte(b'b'), AA_B);
        assert_eq!(encode_protein_byte(b'z'), AA_Z);
    }

    #[test]
    fn ascii_is_rejected_by_encoded_interface() {
        assert!(ProteinCodes::new(b"MKV").is_err());
        assert!(ProteinCodes::new(&[0, 5, 23]).is_ok());
        assert!(ProteinCodes::new(&[24]).is_err());
    }

    #[test]
    fn seedability_excludes_ambiguous_codes() {
        assert!(is_seedable_aa(0));
        assert!(is_seedable_aa(19));
        assert!(!is_seedable_aa(AA_B));
        assert!(!is_seedable_aa(AA_Z));
        assert!(!is_seedable_aa(AA_X));
        assert!(!is_seedable_aa(AA_STOP));
    }

    #[test]
    fn murphy10_has_ten_groups_and_covers_all_standard_residues() {
        let m = SeedAlphabet::murphy10();
        assert_eq!(m.size, 10);
        assert_eq!(m.reduce(encode_protein_byte(b'K')), m.reduce(encode_protein_byte(b'R')));
        assert_eq!(m.reduce(encode_protein_byte(b'I')), m.reduce(encode_protein_byte(b'L')));
        assert_ne!(m.reduce(encode_protein_byte(b'A')), m.reduce(encode_protein_byte(b'S')));
        assert_eq!(m.reduce(AA_X), None);
    }

    #[test]
    fn seed_shape_parsing() {
        let s = SeedShape::parse("11011").unwrap();
        assert_eq!(s.offsets, vec![0, 1, 3, 4]);
        assert_eq!(s.span, 5);
        assert_eq!(s.weight(), 4);
        assert!(!s.is_contiguous());
        assert!(SeedShape::contiguous(4).is_contiguous());
        assert!(SeedShape::parse("0110").is_err());
        assert!(SeedShape::parse("1").is_err());
    }

    #[test]
    fn nucleotide_encoding() {
        assert_eq!(encode_nucleotide(b"ACGTUNacgtx"), vec![0, 1, 2, 3, 3, 4, 0, 1, 2, 3, 4]);
        assert!(NucleotideCodes::new(&[0, 1, 2, 3, 4]).is_ok());
        assert!(NucleotideCodes::new(b"A").is_err());
    }
}
