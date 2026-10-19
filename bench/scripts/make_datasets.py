#!/usr/bin/env python3
"""Derive every benchmark dataset from the two public sources (deterministic).

Sources (see bench/README.md for URLs, commits and SHA-256):
  scop40.fa / scop40.labels.tsv  SCOPe 2.08 40% domains (reseek test_data/scop40.cal.gz)
  mmseqs_QUERY.fasta / mmseqs_DB.fasta  MMseqs2 examples/

Outputs (in OUTDIR):
  scop40.shuf.fa      each domain globally shuffled (composition-preserving decoys)
  scop40.wshuf10.fa   shuffled inside 10-residue windows (local composition kept)
  scop40.bt.fna       back-translated domains (random synonymous codons, random
                      0-60 nt flanks, random strand) for blastx/tblastn tests
  scop40.frag.fa      one random 40-residue window per domain (>= 60 aa): fragments
  scop40.multi.fa     400 two-domain queries "A+B" (different classes, GSGSGS linker)
  lowc.query.fa       100 low-complexity / tandem-repeat queries (adversarial)
  lowc.db.fa          mmseqs DB + 2,000 low-complexity / repeat sequences
  mm.db25.fa mm.db50.fa  every 4th / every 2nd mmseqs DB sequence
  mm.db400.fa         DB x4: original + 3 copies at 15% substitution (renamed)
  mm.q50.fa mm.q100.fa mm.q250.fa  first N queries
  mm.long.fa          mmseqs DB sequences longer than 1,500 aa (long-protein queries)
"""
import os, random, sys

AA = "ARNDCQEGHILKMFPSTWYV"
CODONS = {
    'A': ["GCT", "GCC", "GCA", "GCG"], 'R': ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], 'N': ["AAT", "AAC"],
    'D': ["GAT", "GAC"], 'C': ["TGT", "TGC"], 'Q': ["CAA", "CAG"], 'E': ["GAA", "GAG"],
    'G': ["GGT", "GGC", "GGA", "GGG"], 'H': ["CAT", "CAC"], 'I': ["ATT", "ATC", "ATA"],
    'L': ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], 'K': ["AAA", "AAG"], 'M': ["ATG"], 'F': ["TTT", "TTC"],
    'P': ["CCT", "CCC", "CCA", "CCG"], 'S': ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    'T': ["ACT", "ACC", "ACA", "ACG"], 'W': ["TGG"], 'Y': ["TAT", "TAC"], 'V': ["GTT", "GTC", "GTA", "GTG"],
}

def read_fasta(path):
    recs = []
    hdr = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(seq)))
                hdr = line[1:]
                seq = []
            elif line:
                seq.append(line.strip())
    if hdr is not None:
        recs.append((hdr, "".join(seq)))
    return recs

def write_fasta(path, recs):
    with open(path, "w") as f:
        for h, s in recs:
            f.write(f">{h}\n")
            for i in range(0, len(s), 80):
                f.write(s[i:i + 80] + "\n")

def rid(h):
    return h.split()[0]

def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

def main():
    data, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    scop = read_fasta(os.path.join(data, "scop40.fa"))
    labels = {}
    with open(os.path.join(data, "scop40.labels.tsv")) as f:
        next(f)
        for line in f:
            d, cls, fold, sf, fam, ln = line.rstrip("\n").split("\t")
            labels[d] = (cls, fold, sf, fam)
    mq = read_fasta(os.path.join(data, "mmseqs_QUERY.fasta"))
    mdb = read_fasta(os.path.join(data, "mmseqs_DB.fasta"))

    r = random.Random(20261018)
    shuf = []
    for h, s in scop:
        l = list(s)
        r.shuffle(l)
        shuf.append((rid(h) + "_shuf", "".join(l)))
    write_fasta(os.path.join(out, "scop40.shuf.fa"), shuf)

    w = []
    for h, s in scop:
        l = list(s)
        for i in range(0, len(l), 10):
            seg = l[i:i + 10]
            r.shuffle(seg)
            l[i:i + 10] = seg
        w.append((rid(h) + "_wshuf", "".join(l)))
    write_fasta(os.path.join(out, "scop40.wshuf10.fa"), w)

    bt = []
    for h, s in scop:
        cds = "".join(r.choice(CODONS.get(c, ["NNN"])) for c in s)
        dna = "".join(r.choice("ACGT") for _ in range(r.randint(0, 60))) + cds + \
              "".join(r.choice("ACGT") for _ in range(r.randint(0, 60)))
        if r.random() < 0.5:
            dna = revcomp(dna)
        bt.append((rid(h), dna))
    write_fasta(os.path.join(out, "scop40.bt.fna"), bt)

    frag = []
    for h, s in scop:
        if len(s) >= 60:
            st = r.randint(0, len(s) - 40)
            frag.append((rid(h), s[st:st + 40]))
    write_fasta(os.path.join(out, "scop40.frag.fa"), frag)

    ids = [rid(h) for h, _ in scop]
    seqs = {rid(h): s for h, s in scop}
    multi = []
    seen = set()
    while len(multi) < 400:
        a, b = r.choice(ids), r.choice(ids)
        if labels[a][0] == labels[b][0] or (a, b) in seen:
            continue
        seen.add((a, b))
        multi.append((f"{a}+{b}", seqs[a] + "GSGSGS" + seqs[b]))
    write_fasta(os.path.join(out, "scop40.multi.fa"), multi)

    def lowc(n, rr):
        kind = rr.randrange(6)
        if kind == 0:
            return rr.choice("QNAESGP") * n
        if kind == 1:
            unit = "".join(rr.choice("PESTQ") for _ in range(rr.randint(2, 6)))
        elif kind == 2:
            unit = "GPP"
        elif kind == 3:
            unit = "".join(rr.choice(AA) for _ in range(7))  # heptad-like repeat
        elif kind == 4:
            unit = "".join(rr.choice("KEDR") for _ in range(rr.randint(3, 8)))
        else:
            return "".join(rr.choice("SSSGGGNNQQ") for _ in range(n))
        return (unit * (n // len(unit) + 1))[:n]
    lq = []
    for i in range(100):
        n = r.randint(300, 3000)
        s = lowc(n, r)
        if i % 4 == 0:  # embed a real domain in the low-complexity sequence
            dseq = seqs[r.choice(ids)]
            p = r.randint(0, len(s))
            s = s[:p] + dseq + s[p:]
        lq.append((f"lowc{i}", s))
    write_fasta(os.path.join(out, "lowc.query.fa"), lq)
    ldb = [(rid(h), s) for h, s in mdb] + [(f"lowcdb{i}", lowc(r.randint(100, 2500), r)) for i in range(2000)]
    write_fasta(os.path.join(out, "lowc.db.fa"), ldb)

    write_fasta(os.path.join(out, "mm.db25.fa"), [(rid(h), s) for k, (h, s) in enumerate(mdb) if k % 4 == 0])
    write_fasta(os.path.join(out, "mm.db50.fa"), [(rid(h), s) for k, (h, s) in enumerate(mdb) if k % 2 == 0])
    x4 = [(rid(h), s) for h, s in mdb]
    for c in range(1, 4):
        for h, s in mdb:
            x4.append((f"{rid(h)}_v{c}", "".join(ch if r.random() >= 0.15 else r.choice(AA) for ch in s)))
    write_fasta(os.path.join(out, "mm.db400.fa"), x4)
    for n in (50, 100, 250):
        write_fasta(os.path.join(out, f"mm.q{n}.fa"), [(rid(h), s) for h, s in mq[:n]])
    write_fasta(os.path.join(out, "mm.long.fa"), [(rid(h), s) for h, s in mdb if len(s) > 1500])
    for fn in sorted(os.listdir(out)):
        p = os.path.join(out, fn)
        n = sum(1 for line in open(p) if line.startswith(">"))
        print(f"{fn}\t{n} sequences\t{os.path.getsize(p)} bytes")

if __name__ == "__main__":
    main()
