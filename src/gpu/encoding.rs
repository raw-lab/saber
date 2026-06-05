//! Host-side sequence encoding for GPU kernels.
//!
//! The CUDA kernels in `kernels/sw.cu` expect sequences as packed bytes
//! where each byte is an index into the substitution-matrix alphabet,
//! not an ASCII character. This module converts ASCII FASTA bytes into
//! the encoded form.
//!
//! The encoding MUST match the alphabet ordering used in
//! `scoring::matrices::ALPHABET`, which is "ARNDCQEGHILKMFPSTWYVBZX*".

use crate::scoring::matrices::ScoreMatrix;

/// Protein alphabet size (BLOSUM/PAM 24-letter alphabet).
pub const PROTEIN_ALPHABET_SIZE: usize = 24;

/// Nucleotide alphabet size (ACGT + N).
pub const NUCLEOTIDE_ALPHABET_SIZE: usize = 5;

/// Encode a protein sequence to GPU-ready bytes (0..23).
///
/// Anything not in the 24-letter alphabet maps to X (22). The result is
/// the same length as the input.
pub fn encode_protein(seq: &[u8]) -> Vec<u8> {
    let table = build_protein_table();
    seq.iter().map(|&c| table[c as usize]).collect()
}

/// Encode a nucleotide sequence to GPU-ready bytes (0..4).
///
/// A=0, C=1, G=2, T=3, U=3 (RNA), anything else=4 (N).
pub fn encode_nucleotide(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&c| nuc_encode_byte(c)).collect()
}

#[inline]
fn nuc_encode_byte(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4,
    }
}

fn build_protein_table() -> [u8; 256] {
    // Mirror of scoring::matrices::ALPHABET ordering.
    let alpha = b"ARNDCQEGHILKMFPSTWYVBZX*";
    let mut tbl = [22u8; 256]; // default: X
    for (idx, &c) in alpha.iter().enumerate() {
        tbl[c as usize] = idx as u8;
        tbl[c.to_ascii_lowercase() as usize] = idx as u8;
    }
    tbl
}

/// Pack a list of sequences into one contiguous buffer + offsets array.
///
/// Returns `(packed, offsets)` where:
/// - `packed[offsets[i]..offsets[i+1]]` is the i-th sequence
/// - `offsets.len() == sequences.len() + 1`
pub fn pack_sequences<T, F>(sequences: &[T], encode: F) -> (Vec<u8>, Vec<u32>)
where
    F: Fn(&T) -> Vec<u8>,
{
    let mut packed = Vec::with_capacity(sequences.iter().map(|s| encode(s).len()).sum());
    let mut offsets = Vec::with_capacity(sequences.len() + 1);
    offsets.push(0u32);
    for s in sequences {
        let enc = encode(s);
        packed.extend_from_slice(&enc);
        offsets.push(packed.len() as u32);
    }
    (packed, offsets)
}

/// Flatten a 24x24 substitution matrix into row-major i8 layout for the GPU.
pub fn flatten_matrix(matrix: &ScoreMatrix) -> Vec<i8> {
    let mut out = Vec::with_capacity(PROTEIN_ALPHABET_SIZE * PROTEIN_ALPHABET_SIZE);
    for i in 0..PROTEIN_ALPHABET_SIZE {
        for j in 0..PROTEIN_ALPHABET_SIZE {
            // Clamp to i8 range; standard BLOSUM/PAM scores are in [-17, +17].
            let v = matrix.score_by_index(i, j).clamp(-128, 127);
            out.push(v as i8);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_alanine_encodes_to_zero() {
        assert_eq!(encode_protein(b"A"), vec![0]);
        assert_eq!(encode_protein(b"a"), vec![0]);
    }

    #[test]
    fn protein_unknown_maps_to_x() {
        // 'J' is not in the 24-letter alphabet → X (22).
        assert_eq!(encode_protein(b"J")[0], 22);
    }

    #[test]
    fn protein_full_alphabet_encodes_in_order() {
        let enc = encode_protein(b"ARNDCQEGHILKMFPSTWYVBZX*");
        for (i, &b) in enc.iter().enumerate() {
            assert_eq!(b as usize, i);
        }
    }

    #[test]
    fn nucleotide_encodes_acgtn() {
        assert_eq!(encode_nucleotide(b"ACGTN"), vec![0, 1, 2, 3, 4]);
        assert_eq!(encode_nucleotide(b"acgtn"), vec![0, 1, 2, 3, 4]);
        assert_eq!(encode_nucleotide(b"U"), vec![3]); // RNA
    }

    #[test]
    fn pack_sequences_concatenates_with_offsets() {
        let seqs: Vec<&[u8]> = vec![b"ACG", b"TT", b"GATTACA"];
        let (packed, offsets) = pack_sequences(&seqs, |s| encode_nucleotide(s));
        assert_eq!(offsets, vec![0, 3, 5, 12]);
        assert_eq!(packed.len(), 12);
        // First sequence is ACG = [0,1,2]
        assert_eq!(&packed[offsets[0] as usize..offsets[1] as usize], &[0, 1, 2]);
        // Third is GATTACA = [2,0,3,3,0,1,0]
        assert_eq!(
            &packed[offsets[2] as usize..offsets[3] as usize],
            &[2, 0, 3, 3, 0, 1, 0]
        );
    }
}
