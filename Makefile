.PHONY: help build release test clean docs install bench profile gpu all

RUST_VERSION := 1.75.0
SHELL := /bin/bash

help:
	@echo "SABER Build System"
	@echo "====================="
	@echo ""
	@echo "Available targets:"
	@echo "  make build          - Build debug binary"
	@echo "  make release        - Build optimized release binary"
	@echo "  make test           - Run all tests"
	@echo "  make demo           - Run a smoke test against bundled example FASTA"
	@echo "  make benchmark      - Reproducible SABER 1.0.2/1.1.0/1.2/DIAMOND suite (bench/README.md)"
	@echo "  make benchmark-blast - Quick SABER vs NCBI blastp run (blastp on PATH)"
	@echo "  make bench          - Alignment-kernel micro-benchmarks"
	@echo "  make examples       - Run the example scripts in examples/"
	@echo "  make docs           - Generate documentation"
	@echo "  make install        - Install saber into ~/.cargo/bin (or system bin)"
	@echo "  make uninstall      - Remove saber binary from PATH"
	@echo "  make clean          - Clean build artifacts"
	@echo "  make gpu            - Build with GPU support"
	@echo "  make profile        - Build with profiling symbols"
	@echo "  make all            - Build all variants"
	@echo ""

check-deps:
	@echo "Checking dependencies..."
	@command -v cargo >/dev/null 2>&1 || (echo "Rust/Cargo not found. Install from https://rustup.rs/" && exit 1)
	@echo "✓ Rust toolchain found"

build: check-deps
	@echo "Building debug binary..."
	cargo build --verbose
	@echo "✓ Debug binary: target/debug/saber"

release: check-deps
	@echo "Building optimized release binary..."
	cargo build --release --verbose
	@echo "✓ Release binary: target/release/saber"
	@ls -lh target/release/saber

test: check-deps
	@echo "Running all tests..."
	cargo test --all --verbose
	@echo "✓ All tests passed"

test-release: check-deps
	@echo "Running tests in release mode..."
	cargo test --release --all --verbose
	@echo "✓ Release tests passed"

bench: check-deps
	@echo "Running alignment-kernel micro-benchmarks..."
	cargo run --release --example bench_banded
	cargo run --release --example bench_kernels
	@echo "✓ Benchmarks complete"

docs: check-deps
	@echo "Generating documentation..."
	cargo doc --no-deps --open --all-features
	@echo "✓ Documentation generated"

install: release
	@echo "Installing SABER..."
	cargo install --path . --force
	@echo "✓ Installation complete"
	@echo "  Binary location: $$(which saber)"

uninstall:
	@echo "Uninstalling SABER..."
	@if command -v cargo >/dev/null 2>&1; then \
		cargo uninstall saber 2>/dev/null || true ; \
	fi
	@for d in $$HOME/.cargo/bin /usr/local/bin /usr/bin ; do \
		if [ -f $$d/saber ]; then \
			echo "  removing $$d/saber"; \
			rm -f $$d/saber || sudo rm -f $$d/saber ; \
		fi ; \
	done
	@echo "✓ Uninstall complete"

# Quick functional smoke-test against the bundled FASTA files
demo: release
	@echo "=== blastp (pp) ==="
	./target/release/saber -q examples/protein_query.fa    -d examples/protein_db.fa     --mode pp --outfmt m8 --quiet
	@echo ""
	@echo "=== blastn (nn) ==="
	./target/release/saber -q examples/nucleotide_query.fa -d examples/nucleotide_db.fa  --mode nn --outfmt m8 --quiet
	@echo ""
	@echo "=== blastx (nx, translated query) ==="
	./target/release/saber -q examples/nucleotide_query.fa -d examples/protein_db.fa     --mode nx --outfmt m8 --quiet | head -6
	@echo ""
	@echo "=== tblastn (pn, translated db) ==="
	./target/release/saber -q examples/protein_query.fa    -d examples/nucleotide_db.fa  --mode pn --outfmt m8 --quiet | head -6

# Reproducible comparison with SABER 1.0.2 / 1.1.0 and DIAMOND (bench/README.md).
# Needs the public datasets and the other binaries; set V102 V110 DMD DATA DER
# (see the header of bench/scripts/run_all.sh), optionally SECTIONS="scop decoy"
# and BENCH_OUT (results directory).
BENCH_OUT ?= $(CURDIR)/bench/out
benchmark: release
	OUT=$(BENCH_OUT) V12=$${V12:-$(CURDIR)/target/release/saber} bash bench/scripts/run_all.sh $(SECTIONS)
	python3 -I bench/scripts/summarize.py $(BENCH_OUT)

# Quick side-by-side run against NCBI BLAST+ (blastp on PATH) on the 1.x
# benchmark files in benchmarks/results/; outputs go to benchmarks/out/.
benchmark-blast: release
	@command -v blastp >/dev/null 2>&1 || (echo "blastp not found — apt-get install ncbi-blast+" && exit 1)
	@mkdir -p benchmarks/out
	@echo "=== SABER blastp ==="
	@time ./target/release/saber -q benchmarks/results/bench_query.fa -d benchmarks/results/bench_db_5003seq.fa --mode pp -E 1e-5 --outfmt m8 --quiet > benchmarks/out/saber_bench5k.tsv
	@wc -l benchmarks/out/saber_bench5k.tsv
	@echo "=== NCBI blastp ==="
	@time blastp -query benchmarks/results/bench_query.fa -subject benchmarks/results/bench_db_5003seq.fa -outfmt 6 -evalue 1e-5 > benchmarks/out/blast_bench5k.tsv
	@wc -l benchmarks/out/blast_bench5k.tsv

gpu: check-deps
	@echo "Building with CUDA GPU support..."
	@echo "(Requires CUDA toolkit 11.0+ with libnvrtc available)"
	cargo build --release --features gpu-cuda --verbose
	@echo "✓ GPU-enabled binary: target/release/saber"
	@echo ""
	@echo "Smoke test: make gpu-test (or saber --doctor --gpu for a GPU/CPU score comparison)"

gpu-test: gpu
	@echo "Running GPU smoke test..."
	./target/release/saber -q examples/protein_query.fa -d examples/protein_db.fa \
		--mode pp --gpu --outfmt m8
	./target/release/saber --doctor --gpu
	@echo ""
	@echo "'region scoring on GPU backend ...' on stderr means the GPU path ran;"
	@echo "'--gpu requested but no GPU backend is available' means it fell back to"
	@echo "the CPU — check nvidia-smi and that libnvrtc is on LD_LIBRARY_PATH."

profile: check-deps
	@echo "Building with profiling symbols..."
	CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release
	@echo "✓ Profiling-enabled binary: target/release/saber"

# Development targets
dev-setup: check-deps
	@echo "Setting up development environment..."
	rustup update stable
	rustup component add clippy rustfmt
	@echo "✓ Development environment ready"

fmt:
	@echo "Formatting code..."
	cargo fmt --all

lint:
	@echo "Running clippy lint..."
	cargo clippy --all --all-targets --all-features -- -D warnings

check: fmt lint test
	@echo "✓ All checks passed"

# Build all variants
all: build release gpu profile
	@echo "✓ All variants built successfully"
	@echo ""
	@echo "Build artifacts:"
	@ls -lh target/*/saber 2>/dev/null || true
	@ls -lh target/release/saber 2>/dev/null || true

clean:
	@echo "Cleaning build artifacts..."
	cargo clean
	@echo "✓ Clean complete"

# Performance analysis (Q= and DB= select the inputs; the bundled examples are tiny)
Q ?= examples/protein_query.fa
DB ?= examples/protein_db.fa
perf-record: release
	@echo "Recording performance with Linux perf..."
	perf record -g ./target/release/saber -q $(Q) -d $(DB) -o /dev/null
	@echo "✓ Performance data recorded: perf.data"

perf-report:
	@echo "Analyzing performance..."
	perf report

flamegraph: release
	@echo "Building flamegraph..."
	cargo install flamegraph
	cargo flamegraph --release --bin saber -- -q $(Q) -d $(DB) -o /dev/null
	@echo "✓ Flamegraph generated: flamegraph.svg"

# Coverage
coverage: test
	@echo "Running coverage analysis..."
	cargo tarpaulin --out Html --output-dir coverage
	@echo "✓ Coverage report: coverage/tarpaulin-report.html"

# Example runs
examples: release
	@echo "Running example scripts..."
	@chmod +x examples/*.sh
	@for script in examples/example_*.sh; do \
	    echo ""; \
	    echo "Running $$script..."; \
	    bash "$$script" || echo "⚠ Script exited with error"; \
	done

# Version management
version:
	@cargo --version
	@rustc --version

update-deps:
	@echo "Updating dependencies..."
	cargo update
	@echo "✓ Dependencies updated"

# CI/CD helpers
ci-test: lint test test-release bench
	@echo "✓ CI pipeline passed"

# Publish to crates.io
publish: check
	@echo "Publishing to crates.io..."
	cargo publish

publish-dry: check
	@echo "Publishing (dry run)..."
	cargo publish --dry-run

# Help targets
show-features:
	@echo "Available Cargo features:"
	@grep -A 5 '\[features\]' Cargo.toml

show-deps:
	@echo "Project dependencies:"
	@cargo tree --depth 1

# Default target
.DEFAULT_GOAL := help

# Silent targets
.SILENT: help version

# Phony declarations
.PHONY: help check-deps build release test bench docs install uninstall demo benchmark benchmark-blast clean gpu gpu-test profile \
        dev-setup fmt lint check all perf-record \
        perf-report flamegraph coverage examples version update-deps ci-test \
        publish publish-dry show-features show-deps test-release
