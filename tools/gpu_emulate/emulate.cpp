// CPU emulation harness for src/gpu/kernels/sw.cu.
//
// Compiles the CUDA kernel source as plain C++ (shimming the CUDA built-ins),
// runs it one thread per block (so the cooperative shared-memory load is
// complete before use, as __syncthreads() guarantees on hardware) and compares
// every score with the expected values produced by the Rust CPU reference
// (examples/gpu_vectors.rs). This verifies the kernel *logic* (recurrence,
// gap convention, bias, encoding clamp, int32 range, interleaved work layout);
// it does not exercise the CUDA driver, NVRTC or the host code.
//
//   cargo run --release --example gpu_vectors > vectors.txt
//   g++ -O2 -std=c++17 -o emulate tools/gpu_emulate/emulate.cpp
//   ./emulate vectors.txt
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>

struct Dim3 { unsigned x, y, z; };
static Dim3 threadIdx, blockIdx, blockDim;
#define __global__
#define __shared__ static
#define __restrict__ __restrict
static inline void __syncthreads() {}

#include "../../src/gpu/kernels/sw.cu"

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: emulate vectors.txt\n"); return 2; }
    FILE* f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 2; }
    int ncases = 0, nwin_total = 0, bad = 0;
    char tag[16];
    while (fscanf(f, "%15s", tag) == 1) {
        if (std::string(tag) != "CASE") { fprintf(stderr, "bad tag %s\n", tag); return 2; }
        int qlen, nwin, go, ge;
        if (fscanf(f, "%d %d %d %d", &qlen, &nwin, &go, &ge) != 4) return 2;
        std::vector<uint8_t> q(qlen); std::vector<int8_t> bias(qlen);
        for (int i = 0; i < qlen; i++) { int v; if (fscanf(f, "%d", &v) != 1) return 2; q[i] = (uint8_t)v; }
        for (int i = 0; i < qlen; i++) { int v; if (fscanf(f, "%d", &v) != 1) return 2; bias[i] = (int8_t)v; }
        std::vector<int32_t> m(576);
        for (int i = 0; i < 576; i++) if (fscanf(f, "%d", &m[i]) != 1) return 2;
        std::vector<uint8_t> packed; std::vector<uint32_t> off(1, 0);
        for (int w = 0; w < nwin; w++) {
            int len; if (fscanf(f, "%d", &len) != 1) return 2;
            for (int i = 0; i < len; i++) { int v; if (fscanf(f, "%d", &v) != 1) return 2; packed.push_back((uint8_t)v); }
            off.push_back((uint32_t)packed.size());
        }
        std::vector<long long> expect(nwin);
        for (int w = 0; w < nwin; w++) if (fscanf(f, "%lld", &expect[w]) != 1) return 2;
        if (packed.empty()) packed.push_back(22);
        std::vector<int32_t> H((size_t)qlen * nwin), E((size_t)qlen * nwin), sc(nwin, -1);
        blockDim = {1, 1, 1};
        for (int b = 0; b < nwin; b++) {
            blockIdx = {(unsigned)b, 0, 0}; threadIdx = {0, 0, 0};
            saber_sw_score_i32(q.data(), qlen, bias.data(), packed.data(), off.data(), nwin, m.data(), go, ge, H.data(), E.data(), sc.data());
        }
        for (int w = 0; w < nwin; w++) {
            if ((long long)sc[w] != expect[w]) {
                if (bad < 10) fprintf(stderr, "case %d window %d: kernel %d expected %lld\n", ncases, w, sc[w], expect[w]);
                bad++;
            }
        }
        ncases++; nwin_total += nwin;
    }
    printf("{\"cases\": %d, \"windows\": %d, \"mismatches\": %d, \"result\": \"%s\"}\n", ncases, nwin_total, bad, bad ? "FAIL" : "PASS");
    return bad ? 1 : 0;
}
