# examples

Shell walkthroughs that run on the small FASTA files in this directory. Build
first (`cargo build --release`), then run any script from anywhere; outputs go to
a new temporary directory unless `OUT=/some/dir` is set, and `SABER=/path/to/saber`
selects another binary.

| Script | Shows |
|---|---|
| `example_protein_search.sh` | blastp-like search: default, `--preset sensitive`, BLOSUM45, m0/JSON output, `--makedb` + `--db` |
| `example_nucleotide_search.sh` | blastn-like search (`--mode nn`): defaults, identity filter, alternative scoring, SAM |
| `example_translated_search.sh` | blastx-like (`--mode nx`) and tblastn-like (`--mode pn`) searches, genetic codes, a translated index |
| `example_sam_output.sh` | SAM output for nucleotide and protein alignments; BAM conversion with samtools when installed |
| `example_gpu_search.sh` | `--doctor`, `--gpu` (falls back to the CPU kernels in builds without `--features gpu-cuda`) |

| Data | Content |
|---|---|
| `protein_query.fa`, `protein_db.fa` | human haemoglobin α and a kinase motif; database with globin homologues at decreasing identity plus unrelated and low-complexity decoys |
| `nucleotide_query.fa`, `nucleotide_db.fa` | haemoglobin α coding sequence and an AT repeat; matching, mutated and unrelated database entries |

The scripts written for SABER ≤ 1.1.0 used options that neither 1.0.2 nor
1.1.0 accepted (`-i`, `-v`, `-a`, `-c`, `--gpu-info`, `--streaming`, BAM output)
and input files that were never shipped; they were rewritten for 1.2.

Rust examples (`cargo run --release --example NAME`): `bench_banded`,
`bench_kernels`, `bench_simd` (kernel micro-benchmarks), `engine_probe` (stage
counters for one query set), `mask_probs` (per-residue tantan probabilities, used
by `bench/scripts/validate_masking.sh`), `mask_bench` (single-thread masker
speed and a checksum of the masks, `mask_bench FASTA [max_residues]`),
`gpu_vectors` (GPU test vectors).
