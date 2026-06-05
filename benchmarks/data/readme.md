# Provenance

`bio_queries100.fa` is 100 sequences sampled (seed=42) from MMseqs2's bundled
benchmark queries (`MMseqs2/examples/QUERY.fasta`, 500 total sequences).

The reference DB used in the benchmark is `MMseqs2/examples/DB.fasta`
(20,000 sequences, ~9M residues). Not included here due to size; obtain via:

    git clone --depth=1 https://github.com/soedinglab/MMseqs2.git
    cp MMseqs2/examples/DB.fasta benchmarks/bio/data/bio_db.fa

These are real protein sequences spanning a broad taxonomic range
(bacterial, archaeal, viral, plant, mammalian) curated by the MMseqs2
authors for benchmarking purposes.
