from pathlib import Path

script = r'''set -euo pipefail

THREADS=20
HAP1_FASTA="hap1_hic_output.fasta"
HAP2_FASTA="hap2_hic_output.fasta"
HAP1_BUSCO="hap1_full_table.tsv"
HAP2_BUSCO="hap2_full_table.tsv"
SCAFFOLD_MAP="hap2_to_hap1_scaffold.tsv"
FULL_PAF="hap1_to_hap2.full.paf"

samtools faidx "$HAP1_FASTA"
samtools faidx "$HAP2_FASTA"

minimap2 -t "$THREADS" -c -x asm5 --cs "$HAP2_FASTA" "$HAP1_FASTA" > "$FULL_PAF"

python3 - "$HAP1_BUSCO" "$HAP2_BUSCO" "${HAP2_FASTA}.fai" <<'PY'
import sys
from collections import Counter, defaultdict

h1_file, h2_file, fai_file = sys.argv[1:4]
ORDER = ("Complete", "Duplicated", "Fragmented", "Missing")

def parse_busco(path):
    d = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 2 or x[1] not in ORDER:
                continue
            bid, status = x[0], x[1]
            loc = None
            if status != "Missing" and len(x) >= 5:
                try:
                    s = min(int(x[3]), int(x[4])) - 1
                    e = max(int(x[3]), int(x[4]))
                    loc = (x[2], s, e)
                except ValueError:
                    pass
            d[bid].append((status, loc))

    out = {}
    for bid, rows in d.items():
        statuses = [s for s, _ in rows]
        if "Duplicated" in statuses or statuses.count("Complete") > 1:
            state = "Duplicated"
        elif "Complete" in statuses:
            state = "Complete"
        elif "Fragmented" in statuses:
            state = "Fragmented"
        else:
            state = "Missing"

        coords = [loc for st, loc in rows if st == "Complete" and loc is not None]
        representative = max(coords, key=lambda z: z[2] - z[1]) if coords else None
        out[bid] = (state, representative)
    return out

h1 = parse_busco(h1_file)
h2 = parse_busco(h2_file)
ids = sorted(set(h1) | set(h2))

matrix = Counter()
recovered = []
for bid in ids:
    s1, _ = h1.get(bid, ("Missing", None))
    s2, loc2 = h2.get(bid, ("Missing", None))
    matrix[(s1, s2)] += 1
    if s1 == "Missing" and s2 == "Complete" and loc2 is not None:
        recovered.append((bid, *loc2))

with open("matrix_hap1_hap2.tsv", "w") as out:
    out.write("hap1_status\thap2_status\torthologs\n")
    for a in ORDER:
        for b in ORDER:
            out.write(f"{a}\t{b}\t{matrix[(a,b)]}\n")

with open("missing_to_complete_hap2_coords.tsv", "w") as out:
    out.write("busco_id\thap2_scaffold\tstart_0based\tend_exclusive\n")
    for row in recovered:
        out.write("\t".join(map(str, row)) + "\n")

lengths = {}
with open(fai_file) as f:
    for line in f:
        x = line.rstrip("\n").split("\t")
        lengths[x[0]] = int(x[1])

counts = Counter(scaf for _, scaf, _, _ in recovered)

with open("recovered_buscos_by_contig.tsv", "w") as out:
    out.write("contig\trecovered_BUSCOs\n")
    for scaf in sorted(counts):
        out.write(f"{scaf}\t{counts[scaf]}\n")

with open("recovered_buscos_density_perMb.tsv", "w") as out:
    out.write("Contig\tLength_bp\tBUSCOs\tBUSCOs_per_Mb\n")
    for scaf in sorted(counts):
        length = lengths[scaf]
        n = counts[scaf]
        density = n / (length / 1_000_000)
        out.write(f"{scaf}\t{length}\t{n}\t{density:.6f}\n")
PY

python3 - "$FULL_PAF" "$SCAFFOLD_MAP" "missing_to_complete_hap2_coords.tsv" <<'PY'
import sys

paf, map_file, coords_file = sys.argv[1:4]

targets = set()
with open(coords_file) as f:
    next(f)
    for line in f:
        if line.strip():
            targets.add(line.split("\t")[1])

chrom_map = {}
with open(map_file) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        x = line.split()
        if len(x) < 2:
            continue
        if x[0].lower() == "hap2_scaffold":
            continue
        chrom_map[x[0]] = x[1]

with open(paf) as inp, open("busco_relevant_alignments.paf", "w") as out:
    for line in inp:
        if not line.strip():
            continue
        x = line.rstrip("\n").split("\t")
        if len(x) < 12:
            continue
        query = x[0]
        target = x[5]
        tags = x[12:]
        if target not in targets:
            continue
        if not any(t == "tp:A:P" for t in tags):
            continue
        expected = chrom_map.get(target)
        if expected is not None and query != expected:
            continue
        out.write(line)
PY

gzip -f busco_relevant_alignments.paf

python3 - "missing_to_complete_hap2_coords.tsv" "${HAP2_FASTA}.fai" "$SCAFFOLD_MAP" "busco_relevant_alignments.paf.gz" <<'PY'
import sys, gzip
from collections import Counter, defaultdict

coords_file, fai_file, map_file, paf_gz = sys.argv[1:5]

lengths = {}
with open(fai_file) as f:
    for line in f:
        x = line.rstrip("\n").split("\t")
        lengths[x[0]] = int(x[1])

chrom_map = {}
with open(map_file) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        x = line.split()
        if len(x) < 2 or x[0].lower() == "hap2_scaffold":
            continue
        chrom_map[x[0]] = x[1]

alignments = defaultdict(list)
with gzip.open(paf_gz, "rt") as f:
    for line in f:
        if not line.strip():
            continue
        x = line.rstrip("\n").split("\t")
        if len(x) < 12:
            continue
        query = x[0]
        target = x[5]
        tstart = int(x[7])
        tend = int(x[8])
        alignments[target].append((query, tstart, tend))

def union_coverage(intervals, start, end):
    clipped = []
    for a, b in intervals:
        a = max(a, start)
        b = min(b, end)
        if a < b:
            clipped.append((a, b))
    if not clipped:
        return 0
    clipped.sort()
    total = 0
    cs, ce = clipped[0]
    for a, b in clipped[1:]:
        if a <= ce:
            ce = max(ce, b)
        else:
            total += ce - cs
            cs, ce = a, b
    total += ce - cs
    return total

rows = []
category_counts = defaultdict(Counter)

with open(coords_file) as f:
    next(f)
    for line in f:
        if not line.strip():
            continue
        bid, scaf, s, e = line.rstrip("\n").split("\t")
        s, e = int(s), int(e)
        expected = chrom_map.get(scaf, "UNRESTRICTED")

        ints = []
        for query, a, b in alignments.get(scaf, []):
            if expected == "UNRESTRICTED" or query == expected:
                ints.append((a, b))

        flank_s = max(0, s - 50000)
        flank_e = min(lengths[scaf], e + 50000)

        locus_bp = union_coverage(ints, s, e)
        flank_bp = union_coverage(ints, flank_s, flank_e)

        locus_frac = locus_bp / max(1, e - s)
        flank_frac = flank_bp / max(1, flank_e - flank_s)

        if locus_bp > 0:
            category = "LOCUS_ALIGNED"
        elif flank_bp > 0:
            category = "FLANK_ONLY"
        else:
            category = "NO_HOMOLOGY_50KB"

        category_counts[scaf][category] += 1
        rows.append((bid, scaf, s, e, expected, locus_frac, flank_frac, category))

with open("BUSCO_locus_flank_classification.tsv", "w") as out:
    out.write("busco_id\thap2_scaffold\tstart_0based\tend_exclusive\texpected_hap1_scaffold\tlocus_coverage_fraction\tflank50kb_coverage_fraction\tcategory\n")
    for bid, scaf, s, e, expected, lc, fc, cat in rows:
        out.write(f"{bid}\t{scaf}\t{s}\t{e}\t{expected}\t{lc:.6f}\t{fc:.6f}\t{cat}\n")

with open("BUSCO_categories_by_scaffold.tsv", "w") as out:
    out.write("hap2_scaffold\tLOCUS_ALIGNED\tFLANK_ONLY\tNO_HOMOLOGY_50KB\tTOTAL\n")
    for scaf in sorted(category_counts):
        c = category_counts[scaf]
        total = c["LOCUS_ALIGNED"] + c["FLANK_ONLY"] + c["NO_HOMOLOGY_50KB"]
        out.write(
            f"{scaf}\t{c['LOCUS_ALIGNED']}\t{c['FLANK_ONLY']}\t"
            f"{c['NO_HOMOLOGY_50KB']}\t{total}\n"
        )
PY

rm -f "$FULL_PAF"
'''

out = Path("/mnt/data/BUSCO_haplotype_comparison_commands.sh")
out.write_text(script, encoding="utf-8")
print(out)
