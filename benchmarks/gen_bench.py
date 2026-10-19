#!/usr/bin/env python3
"""Generate a realistic protein benchmark: UniProt-like length distribution,
real background AA composition, and planted homolog families at controlled
identity so recall is measurable (ground truth known)."""
import random, math, sys

random.seed(20260805)

# UniProt/Swiss-Prot background amino-acid frequencies
AA = "ARNDCQEGHILKMFPSTWYV"
FREQ = [.0826,.0553,.0406,.0546,.0137,.0393,.0674,.0708,.0227,.0593,
        .0966,.0582,.0241,.0386,.0472,.0661,.0535,.0109,.0292,.0687]
CUM = []
s = 0.0
for f in FREQ:
    s += f
    CUM.append(s)
for i in range(len(CUM)):
    CUM[i] /= s

# BLOSUM62-ish conservative substitution groups (for realistic homolog drift)
GROUPS = ["AGST","ILMV","FWY","DEQN","KRH","C","P"]
SUBMAP = {}
for g in GROUPS:
    for c in g:
        SUBMAP[c] = g

def rand_aa():
    r = random.random()
    lo, hi = 0, len(CUM)-1
    while lo < hi:
        mid = (lo+hi)//2
        if r > CUM[mid]: lo = mid+1
        else: hi = mid
    return AA[lo]

def rand_seq(n):
    return "".join(rand_aa() for _ in range(n))

def rand_len():
    # log-normal, mean ~350 aa, matches UniProt shape
    L = int(random.lognormvariate(5.6, 0.55))
    return max(50, min(L, 3000))

def mutate(seq, identity):
    """Produce a homolog at approximately `identity` fraction identical.
    Uses conservative substitutions (biologically realistic) plus indels."""
    out = []
    for c in seq:
        r = random.random()
        if r < identity:
            out.append(c)                       # conserved
        elif r < identity + (1-identity)*0.80:
            g = SUBMAP.get(c, AA)               # conservative substitution
            out.append(random.choice(g))
        elif r < identity + (1-identity)*0.93:
            out.append(rand_aa())               # random substitution
        elif r < identity + (1-identity)*0.97:
            pass                                # deletion
        else:
            out.append(c); out.append(rand_aa())  # insertion
    return "".join(out)

def write_fa(path, recs):
    with open(path, "w") as f:
        for name, seq in recs:
            f.write(">%s\n" % name)
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

n_db     = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
n_query  = int(sys.argv[2]) if len(sys.argv) > 2 else 200
prefix   = sys.argv[3] if len(sys.argv) > 3 else "bench"

# --- Build database: 85% unrelated background, 15% in homolog families ---
db = []
families = []          # family seed sequences we will query with
n_fam = max(1, n_query)
fam_sizes = [random.randint(2, 6) for _ in range(n_fam)]

for fi in range(n_fam):
    seed = rand_seq(rand_len())
    families.append(seed)
    for mi in range(fam_sizes[fi]):
        # family members span 25%-90% identity: the hard remote-homology range
        ident = random.uniform(0.25, 0.90)
        db.append(("fam%d_m%d_id%d" % (fi, mi, int(ident*100)), mutate(seed, ident)))

while len(db) < n_db:
    db.append(("bg%d" % len(db), rand_seq(rand_len())))

random.shuffle(db)
write_fa(prefix + "_db.fa", db)

# --- Queries: the family seeds (each has known homologs in the DB) ---
queries = [("query%d" % i, families[i]) for i in range(n_query)]
write_fa(prefix + "_query.fa", queries)

tot = sum(len(s) for _, s in db)
print("db: %d seqs, %d residues (mean %.0f aa)" % (len(db), tot, tot/len(db)))
print("query: %d seqs, %d residues" % (len(queries), sum(len(s) for _, s in queries)))
print("planted homologs: %d across %d families" % (sum(fam_sizes), n_fam))
