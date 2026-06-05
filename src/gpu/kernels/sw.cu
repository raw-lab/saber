// SABER GPU Smith-Waterman kernels
// =================================
// Inter-task parallelism: each CUDA thread computes the SW score of one
// (query, subject) pair. Subjects are processed in parallel across the grid.
//
// References:
//   - Liu, Y. et al. "CUDASW++: optimizing Smith-Waterman sequence database
//     searches for CUDA-enabled graphics processing units." BMC Res. Notes (2009)
//   - Korpar, M. & Sikic, M. "SW#–GPU-enabled exact alignments on genome scale."
//     Bioinformatics (2013)
//
// Memory model:
//   * Query: regular global memory; broadcast read pattern is L1/L2 friendly.
//     (Constant memory would be ~5% faster but introduces 64KB length limit.)
//   * Substitution matrix: shared memory (24x24 = 576 bytes per block).
//   * Per-thread DP rows H[query_len], E[query_len]: global memory.
//     A 5000-aa query × 100k subjects = ~2GB of DP state. Use batching
//     on the host side to stay under device memory.
//
// Algorithm: Gotoh affine-gap Smith-Waterman with local zero-floor.
//   Recurrence (subject as outer dim i, query as inner dim j):
//     H[i][j] = max(H[i-1][j-1] + sub(s_i, q_j),    // diagonal: match/mismatch
//                   E[i][j],                          // gap in query
//                   F[i][j],                          // gap in subject
//                   0)                                // local: never negative
//     E[i][j] = max(E[i-1][j] - gap_extend,
//                   H[i-1][j] - gap_open)
//     F[i][j] = max(F[i][j-1] - gap_extend,
//                   H[i][j-1] - gap_open)
//
// We track H and E across subject rows (i direction). F is computed inline
// during the inner query loop because it propagates only within a single row.
//
// Alphabet encoding (must match host-side encoding in src/gpu/encoding.rs):
//   Protein 24-letter: ARNDCQEGHILKMFPSTWYVBZX* → 0..23
//   Nucleotide 5-letter: ACGTN → 0..4 (we still use the 24x24 matrix shape;
//   only the top-left 5x5 region matters for nucleotide mode).

#include <stdint.h>

#define SABER_ALPHABET_SIZE 24
#define SABER_NEG_INF      ((int16_t)-30000)

// ---------------------------------------------------------------------------
// Protein Smith-Waterman, score-only.
//
// Grid: launch with num_subjects threads total.
// Block size: tunable (128 is a reasonable default; see CudaExecutor::launch).
//
// Input layout:
//   - query[0..query_len]: encoded query (0..23 per residue)
//   - subjects[0..total_residues]: concatenated encoded subjects
//   - subj_offsets[0..num_subjects+1]: prefix-sum, subject i is
//     subjects[subj_offsets[i]..subj_offsets[i+1]]
//   - sub_matrix[24*24]: row-major substitution matrix (int8)
//   - gap_open, gap_extend: positive penalty values
//   - work_H, work_E: scratch buffers, sized [num_subjects * query_len] each
//   - scores: output, one int32 per subject
//
// Returned score is the maximum SW score across the entire DP table for
// that (query, subject) pair. Traceback is NOT performed on the GPU —
// the host runs CPU traceback on the top-K hits.
// ---------------------------------------------------------------------------
extern "C" __global__ void saber_sw_protein_score(
    const uint8_t* __restrict__ query,
    int32_t                     query_len,
    const uint8_t* __restrict__ subjects,
    const uint32_t* __restrict__ subj_offsets,
    int32_t                     num_subjects,
    const int8_t* __restrict__  sub_matrix,
    int32_t                     gap_open,
    int32_t                     gap_extend,
    int16_t* __restrict__       work_H,
    int16_t* __restrict__       work_E,
    int32_t* __restrict__       scores
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Cooperative load of the substitution matrix into shared memory.
    __shared__ int8_t s_sub[SABER_ALPHABET_SIZE * SABER_ALPHABET_SIZE];
    for (int i = threadIdx.x; i < SABER_ALPHABET_SIZE * SABER_ALPHABET_SIZE; i += blockDim.x) {
        s_sub[i] = sub_matrix[i];
    }
    __syncthreads();

    if (tid >= num_subjects) return;

    const uint32_t subj_start = subj_offsets[tid];
    const uint32_t subj_end   = subj_offsets[tid + 1];
    const int      subj_len   = (int)(subj_end - subj_start);

    // Per-thread DP rows. work_H and work_E are sized [num_subjects][query_len].
    int16_t* H = work_H + (size_t)tid * query_len;
    int16_t* E = work_E + (size_t)tid * query_len;

    // Initialize H[0..query_len] = 0, E[0..query_len] = NEG_INF.
    // E = NEG_INF prevents spurious gap-opens at i = 0.
    for (int j = 0; j < query_len; j++) {
        H[j] = 0;
        E[j] = SABER_NEG_INF;
    }

    int32_t max_score = 0;

    const int16_t go = (int16_t)gap_open;
    const int16_t ge = (int16_t)gap_extend;

    for (int i = 0; i < subj_len; i++) {
        // Load and bounds-check the encoded subject residue.
        uint8_t s_raw = subjects[subj_start + i];
        // Out-of-range encoding maps to X (22) — defensive against malformed input.
        const int s_enc = (s_raw < SABER_ALPHABET_SIZE) ? (int)s_raw : 22;

        int16_t H_diag = 0;          // H[i-1][j-1] going into iteration j
        int16_t H_left = 0;          // H[i][j-1]
        int16_t F      = SABER_NEG_INF;  // gap-in-subject state, resets each row

        // Pointer to this subject's row in the substitution matrix.
        const int8_t* sub_row = &s_sub[s_enc * SABER_ALPHABET_SIZE];

        for (int j = 0; j < query_len; j++) {
            const uint8_t q_raw = query[j];
            const int q_enc = (q_raw < SABER_ALPHABET_SIZE) ? (int)q_raw : 22;
            const int16_t match_score = (int16_t)sub_row[q_enc];

            // Save H[i-1][j] before we overwrite — it becomes H_diag for j+1.
            const int16_t next_diag = H[j];

            // E[j]: gap-in-query state, propagates down (across subject rows).
            //   open from H[i-1][j], extend from E[i-1][j].
            int16_t e_ext  = E[j] - ge;
            int16_t e_open = H[j] - go;
            int16_t e      = (e_ext > e_open) ? e_ext : e_open;

            // F: gap-in-subject state, propagates right (within current row).
            int16_t f_ext  = F - ge;
            int16_t f_open = H_left - go;
            int16_t f      = (f_ext > f_open) ? f_ext : f_open;

            // H[i][j] = max(diag + sub, e, f, 0).
            int16_t h = H_diag + match_score;
            if (e > h) h = e;
            if (f > h) h = f;
            if (h < 0) h = 0;

            // Commit.
            H[j]    = h;
            E[j]    = e;
            F       = f;
            H_diag  = next_diag;
            H_left  = h;

            if ((int32_t)h > max_score) max_score = (int32_t)h;
        }
    }

    scores[tid] = max_score;
}

// ---------------------------------------------------------------------------
// Nucleotide Smith-Waterman, score-only.
//
// Identical structure to the protein kernel, but uses simple match/mismatch
// scoring inline rather than a substitution matrix lookup. This is ~10%
// faster than the protein path due to one fewer memory load per cell.
// ---------------------------------------------------------------------------
extern "C" __global__ void saber_sw_nucleotide_score(
    const uint8_t* __restrict__ query,
    int32_t                     query_len,
    const uint8_t* __restrict__ subjects,
    const uint32_t* __restrict__ subj_offsets,
    int32_t                     num_subjects,
    int32_t                     match_score,
    int32_t                     mismatch_penalty,  // negative
    int32_t                     gap_open,
    int32_t                     gap_extend,
    int16_t* __restrict__       work_H,
    int16_t* __restrict__       work_E,
    int32_t* __restrict__       scores
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_subjects) return;

    const uint32_t subj_start = subj_offsets[tid];
    const uint32_t subj_end   = subj_offsets[tid + 1];
    const int      subj_len   = (int)(subj_end - subj_start);

    int16_t* H = work_H + (size_t)tid * query_len;
    int16_t* E = work_E + (size_t)tid * query_len;
    for (int j = 0; j < query_len; j++) {
        H[j] = 0;
        E[j] = SABER_NEG_INF;
    }

    int32_t max_score = 0;

    const int16_t go     = (int16_t)gap_open;
    const int16_t ge     = (int16_t)gap_extend;
    const int16_t ms     = (int16_t)match_score;
    const int16_t mm     = (int16_t)mismatch_penalty;

    for (int i = 0; i < subj_len; i++) {
        const uint8_t s = subjects[subj_start + i];
        int16_t H_diag = 0;
        int16_t H_left = 0;
        int16_t F      = SABER_NEG_INF;

        for (int j = 0; j < query_len; j++) {
            const uint8_t q = query[j];
            // N (encoded as 4) treats as mismatch against everything.
            const int16_t score_here = (s == q && s < 4) ? ms : mm;

            const int16_t next_diag = H[j];

            int16_t e_ext  = E[j] - ge;
            int16_t e_open = H[j] - go;
            int16_t e      = (e_ext > e_open) ? e_ext : e_open;

            int16_t f_ext  = F - ge;
            int16_t f_open = H_left - go;
            int16_t f      = (f_ext > f_open) ? f_ext : f_open;

            int16_t h = H_diag + score_here;
            if (e > h) h = e;
            if (f > h) h = f;
            if (h < 0) h = 0;

            H[j]    = h;
            E[j]    = e;
            F       = f;
            H_diag  = next_diag;
            H_left  = h;

            if ((int32_t)h > max_score) max_score = (int32_t)h;
        }
    }

    scores[tid] = max_score;
}

// ---------------------------------------------------------------------------
// K-mer prefilter kernel: count how many query k-mer hashes appear in
// each subject. Used by the host to rank candidates before running full SW.
//
// Inputs:
//   - query_kmer_hashes: sorted array of unique k-mer hashes from the query
//   - subjects, subj_offsets, num_subjects: as above (but raw, not encoded)
//   - k: k-mer length (typically 3-5 for protein, 11-15 for nucleotide)
//   - alphabet_size: 24 for protein, 5 for nucleotide
//
// Output: hit_counts[num_subjects] = # of query k-mers that occur in subject.
//
// Implementation: each thread handles one subject. For each k-mer in the
// subject, it hashes and binary-searches the query k-mer table. O(L log N)
// per subject where L = subject length, N = # query k-mers. For typical
// protein queries (N ~= 100-1000), this is fast on GPU.
// ---------------------------------------------------------------------------
extern "C" __global__ void saber_prefilter_count(
    const uint64_t* __restrict__ query_kmer_hashes,
    int32_t                      num_query_kmers,
    const uint8_t* __restrict__  subjects,
    const uint32_t* __restrict__ subj_offsets,
    int32_t                      num_subjects,
    int32_t                      k,
    int32_t                      alphabet_size,
    uint32_t* __restrict__       hit_counts
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_subjects) return;

    const uint32_t subj_start = subj_offsets[tid];
    const uint32_t subj_end   = subj_offsets[tid + 1];
    const int      subj_len   = (int)(subj_end - subj_start);
    if (subj_len < k) {
        hit_counts[tid] = 0;
        return;
    }

    uint32_t hits = 0;

    for (int pos = 0; pos <= subj_len - k; pos++) {
        // Compute polynomial rolling hash matching host-side encoding.
        uint64_t h = 0;
        for (int kk = 0; kk < k; kk++) {
            uint8_t b = subjects[subj_start + pos + kk];
            h = h * (uint64_t)alphabet_size + (uint64_t)b;
        }

        // Binary search in sorted query_kmer_hashes.
        int lo = 0, hi = num_query_kmers - 1;
        while (lo <= hi) {
            int mid = (lo + hi) >> 1;
            uint64_t v = query_kmer_hashes[mid];
            if (v == h) { hits++; break; }
            if (v < h) lo = mid + 1;
            else hi = mid - 1;
        }
    }

    hit_counts[tid] = hits;
}
