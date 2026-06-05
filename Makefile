.PHONY: help build release test clean docs install bench profile gpu all

RUST_VERSION := 1.70.0
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
	@echo "  make benchmark      - Run SABER vs NCBI BLAST side-by-side"
	@echo "  make bench          - Run cargo bench"
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

bench: check-deps release
	@echo "Running benchmarks..."
	cargo bench --release --all
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

# Side-by-side comparison against BLAST+ (must be on PATH)
benchmark: release
	@command -v blastp >/dev/null 2>&1 || (echo "blastp not found — apt-get install ncbi-blast+" && exit 1)
	@mkdir -p benchmarks
	@echo "=== SABER blastp ==="
	@time ./target/release/saber -q benchmarks/bench_query.fa -d benchmarks/bench_db_5003seq.fa --mode pp -E 1e-5 --outfmt m8 --quiet > benchmarks/saber_bench5k.tsv
	@wc -l benchmarks/saber_bench5k.tsv
	@echo "=== NCBI blastp ==="
	@time blastp -query benchmarks/bench_query.fa -subject benchmarks/bench_db_5003seq.fa -outfmt 6 -evalue 1e-5 > benchmarks/blast_bench5k.tsv
	@wc -l benchmarks/blast_bench5k.tsv
	@echo "✓ Results in benchmarks/. See benchmarks/REPORT.md for analysis."

gpu: check-deps
	@echo "Building with CUDA GPU support..."
	@echo "(Requires CUDA toolkit 11.0+ with libnvrtc available)"
	cargo build --release --features gpu-cuda --verbose
	@echo "✓ GPU-enabled binary: target/release/saber"
	@echo ""
	@echo "Smoke test (--gpu flag should detect CUDA device and use it):"
	@echo "  ./target/release/saber -q examples/protein_query.fa -d examples/protein_db.fa \\"
	@echo "        --mode pp --gpu --prefilter --outfmt m8"

gpu-test: gpu
	@echo "Running GPU smoke test..."
	./target/release/saber -q examples/protein_query.fa -d examples/protein_db.fa \
		--mode pp --gpu --prefilter --outfmt m8 --verbose
	@echo ""
	@echo "If the output above shows 'Alignment backend: cuda', the GPU path ran."
	@echo "If it shows 'Alignment backend: cpu-fallback', the CUDA device was not"
	@echo "detected — check nvidia-smi and that libnvrtc is on LD_LIBRARY_PATH."

profile: check-deps
	@echo "Building with profiling symbols..."
	CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release
	@echo "✓ Profiling-enabled binary: target/release/saber"

# Development targets
dev-setup: check-deps
	@echo "Setting up development environment..."
	rustup update stable
	rustup install clippy
	rustup install rustfmt
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

# Docker build
docker-build:
	@echo "Building Docker image..."
	docker build -t saber:latest .
	@echo "✓ Docker image built"

docker-run:
	docker run --rm -it saber:latest saber -h

# Performance analysis
perf-record: release
	@echo "Recording performance with Linux perf..."
	perf record -g ./target/release/saber -i examples/queries.fasta -d examples/database.fasta
	@echo "✓ Performance data recorded: perf.data"

perf-report:
	@echo "Analyzing performance..."
	perf report

flamegraph: release
	@echo "Building flamegraph..."
	cargo install flamegraph
	cargo flamegraph --release -- -i examples/queries.fasta -d examples/database.fasta
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
.PHONY: help check-deps build release test bench docs install uninstall demo benchmark clean gpu gpu-test profile \
        dev-setup fmt lint check all docker-build docker-run perf-record \
        perf-report flamegraph coverage examples version update-deps ci-test \
        publish publish-dry show-features show-deps test-release
