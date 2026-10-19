// SABER GPU Smith-Waterman score kernel (v1.2)
// ============================================
// Inter-task parallelism: each CUDA thread computes the local-alignment score
// of the query against one subject window.
//
// Contract (must match src/gpu/mod.rs and the CPU kernel exactly):
//   * query / subjects are canonical alignment codes (protein ARNDCQEGHILKMFPSTWYVBZX*
//     -> 0..23, nucleotide ACGTN -> 0..4); codes >= 24 are scored as X (22).
//     The host never re-encodes and neither does this kernel.
//   * substitution matrix: 24x24 int32, row = query code, column = subject code.
//   * gaps (SABER convention): a gap of length k costs gap_open + (k-1)*gap_extend.
//   * bias[j] (int8, one per query position, zeros when unused) is added to every
//     aligned column involving query position j.
//   * all DP values are int32: no saturation (v1.1 used int16 and silently clamped).
//
// Work arrays are interleaved [query_len][num_subjects] so that, at a given
// query position, consecutive threads touch consecutive words (coalesced).
//
// Recurrence (subject outer index i, query inner index j):
//   E[i][j] = max(E[i-1][j] - ge, H[i-1][j] - go)        gap in query (vertical)
//   F[i][j] = max(F[i][j-1] - ge, H[i][j-1] - go)        gap in subject (horizontal)
//   H[i][j] = max(0, H[i-1][j-1] + S(q_j, s_i) + bias_j, E[i][j], F[i][j])
//
// The file is also compiled as plain C++ by gpu/emulate/ (with a small shim for
// the CUDA built-ins) to test the kernel logic against the CPU reference on
// machines without a GPU.

typedef unsigned char      uint8_t;
typedef unsigned int       uint32_t;
typedef signed char        int8_t;
typedef int                int32_t;

#define SABER_ALPHABET_SIZE 24
#define SABER_NEG_INF      (-(1 << 29))

extern "C" __global__ void saber_sw_score_i32(
    const uint8_t*  __restrict__ query,
    int32_t                      query_len,
    const int8_t*   __restrict__ bias,
    const uint8_t*  __restrict__ subjects,
    const uint32_t* __restrict__ subj_offsets,
    int32_t                      num_subjects,
    const int32_t*  __restrict__ sub_matrix,
    int32_t                      gap_open,
    int32_t                      gap_extend,
    int32_t*        __restrict__ work_H,
    int32_t*        __restrict__ work_E,
    int32_t*        __restrict__ scores
) {
    __shared__ int32_t s_sub[SABER_ALPHABET_SIZE * SABER_ALPHABET_SIZE];
    for (int i = threadIdx.x; i < SABER_ALPHABET_SIZE * SABER_ALPHABET_SIZE; i += blockDim.x) {
        s_sub[i] = sub_matrix[i];
    }
    __syncthreads();

    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_subjects) return;

    const uint32_t s0 = subj_offsets[tid];
    const int      sl = (int)(subj_offsets[tid + 1] - s0);
    const int      n  = num_subjects;

    for (int j = 0; j < query_len; j++) {
        work_H[(long long)j * n + tid] = 0;
        work_E[(long long)j * n + tid] = SABER_NEG_INF;
    }

    int32_t best = 0;
    for (int i = 0; i < sl; i++) {
        const uint8_t sr = subjects[s0 + i];
        const int sc = (sr < SABER_ALPHABET_SIZE) ? (int)sr : 22;
        int32_t h_diag = 0;              // H[i-1][j-1]
        int32_t h_left = 0;              // H[i][j-1]
        int32_t f      = SABER_NEG_INF;  // F[i][j-1]
        for (int j = 0; j < query_len; j++) {
            const long long at = (long long)j * n + tid;
            const uint8_t qr = query[j];
            const int qc = (qr < SABER_ALPHABET_SIZE) ? (int)qr : 22;
            const int32_t h_up = work_H[at];  // H[i-1][j]
            int32_t e = work_E[at] - gap_extend;
            const int32_t eo = h_up - gap_open;
            if (eo > e) e = eo;
            int32_t ff = f - gap_extend;
            const int32_t fo = h_left - gap_open;
            if (fo > ff) ff = fo;
            int32_t h = h_diag + s_sub[qc * SABER_ALPHABET_SIZE + sc] + (int32_t)bias[j];
            if (e > h) h = e;
            if (ff > h) h = ff;
            if (h < 0) h = 0;
            work_H[at] = h;
            work_E[at] = e;
            f = ff;
            h_diag = h_up;
            h_left = h;
            if (h > best) best = h;
        }
    }
    scores[tid] = best;
}
