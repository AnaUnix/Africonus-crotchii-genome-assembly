#!/usr/bin/env python3
import re, collections, zipfile
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

HAP_FILES = {
    "Hap1": "braker.gtf",
    "Hap2": "h2_braker.gtf",
}

def parse_attr(s):
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', s):
        d[m.group(1)] = m.group(2)
    for part in s.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            d.setdefault(k.strip(), v.strip())
    return d

def base_gene(x):
    return re.sub(r'\.t\d+$', '', x) if x else None

def transcript_from_line(fields, attrs):
    feat = fields[2].lower()
    if feat in ("transcript", "mrna"):
        raw = fields[8].strip().split(";")[0].strip()
        if raw and "=" not in raw and " " not in raw and '"' not in raw:
            return raw
        for k in ("ID", "transcript_id", "transcript"):
            if k in attrs:
                return attrs[k]
    for k in ("transcript_id", "transcript", "ID"):
        if k in attrs:
            return attrs[k]
    if "Parent" in attrs and re.search(r'\.t\d+$', attrs["Parent"]):
        return attrs["Parent"]
    return None

def gene_from_line(fields, attrs, tx=None):
    raw = fields[8].strip().split(";")[0].strip()
    if fields[2].lower() == "gene" and raw and "=" not in raw and " " not in raw and '"' not in raw:
        return base_gene(raw)
    if tx:
        return base_gene(tx)
    for k in ("gene_id", "gene", "Parent"):
        if k in attrs:
            return base_gene(attrs[k])
    return base_gene(raw) if raw and "=" not in raw and " " not in raw and '"' not in raw else None

def process_file(filename, hap):
    selected_tx = {}
    gene_coords = {}
    tx_exons = collections.defaultdict(list)
    counts = collections.Counter()

    with open(filename, errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2].lower()
            counts[feature] += 1
            chrom = fields[0]
            start, end = int(fields[3]), int(fields[4])
            attrs = parse_attr(fields[8])
            tx = transcript_from_line(fields, attrs)
            gid = gene_from_line(fields, attrs, tx)

            if feature == "gene" and gid:
                gene_coords.setdefault(gid, (chrom, start, end))

            if feature in ("transcript", "mrna") and tx:
                gid = base_gene(tx)
                # First transcript/model encountered for each gene wins.
                if gid not in selected_tx:
                    selected_tx[gid] = tx

            if feature == "exon" and tx:
                gid = base_gene(tx)
                if selected_tx.get(gid) == tx:
                    tx_exons[tx].append((chrom, start, end))

    gene_rows, exon_rows, intron_rows = [], [], []

    for gid, tx in selected_tx.items():
        exons = sorted(tx_exons[tx], key=lambda x: (x[0], x[1], x[2]))
        if not exons:
            continue

        if gid in gene_coords:
            chrom, gstart, gend = gene_coords[gid]
        else:
            chrom = exons[0][0]
            gstart = min(x[1] for x in exons)
            gend = max(x[2] for x in exons)

        gene_rows.append({
            "haplotype": hap, "gene_id": gid, "selected_transcript": tx,
            "chromosome": chrom, "gene_start": gstart, "gene_end": gend,
            "gene_length_bp": gend - gstart + 1,
            "n_exons": len(exons), "n_introns": max(0, len(exons) - 1)
        })

        for i, (c, s, e) in enumerate(exons, 1):
            exon_rows.append({
                "haplotype": hap, "gene_id": gid, "selected_transcript": tx,
                "exon_number": i, "chromosome": c, "start": s, "end": e,
                "length_bp": e - s + 1
            })

        for i in range(len(exons) - 1):
            intron_len = exons[i+1][1] - exons[i][2] - 1
            if intron_len > 0:
                intron_rows.append({
                    "haplotype": hap, "gene_id": gid, "selected_transcript": tx,
                    "intron_number": i + 1, "chromosome": chrom,
                    "start": exons[i][2] + 1, "end": exons[i+1][1] - 1,
                    "length_bp": intron_len
                })

    return (pd.DataFrame(gene_rows), pd.DataFrame(exon_rows),
            pd.DataFrame(intron_rows), selected_tx, counts)

def make_boxplot(ax, data, title):
    bp = ax.boxplot(data, labels=["Hap1", "Hap2"], patch_artist=True,
                    showfliers=False, widths=0.6)
    bp["boxes"][0].set_facecolor("green")
    bp["boxes"][1].set_facecolor("purple")
    for box in bp["boxes"]:
        box.set_alpha(0.65)
    ax.set_title(title)
    ax.set_ylabel("Length (bp)")
    ax.grid(axis="y", alpha=0.25)

def main():
    outdir = Path("boxplots_hap1_hap2")
    outdir.mkdir(exist_ok=True)

    res = {}
    for hap, fn in HAP_FILES.items():
        res[hap] = process_file(fn, hap)

    gene_df = pd.concat([res["Hap1"][0], res["Hap2"][0]], ignore_index=True)
    exon_df = pd.concat([res["Hap1"][1], res["Hap2"][1]], ignore_index=True)
    intron_df = pd.concat([res["Hap1"][2], res["Hap2"][2]], ignore_index=True)

    gene_df.to_csv(outdir/"gene_lengths.tsv", sep="\t", index=False)
    exon_df.to_csv(outdir/"exon_lengths.tsv", sep="\t", index=False)
    intron_df.to_csv(outdir/"intron_lengths.tsv", sep="\t", index=False)

    metrics = []
    for hap in ("Hap1", "Hap2"):
        for feature, df, col in (
            ("gene", res[hap][0], "gene_length_bp"),
            ("exon", res[hap][1], "length_bp"),
            ("intron", res[hap][2], "length_bp"),
        ):
            metrics.extend((hap, feature, int(x)) for x in df[col].dropna())
    pd.DataFrame(metrics, columns=["haplotype", "feature", "length_bp"]).to_csv(
        outdir/"all_feature_lengths.tsv", sep="\t", index=False)

    summary = []
    for hap in ("Hap1", "Hap2"):
        for feature, df, col in (
            ("gene", res[hap][0], "gene_length_bp"),
            ("exon", res[hap][1], "length_bp"),
            ("intron", res[hap][2], "length_bp"),
        ):
            s = df[col].dropna()
            summary.append({
                "haplotype": hap, "feature": feature, "n": len(s),
                "mean_bp": s.mean(), "median_bp": s.median(),
                "q1_bp": s.quantile(.25), "q3_bp": s.quantile(.75),
                "min_bp": s.min(), "max_bp": s.max()
            })
    pd.DataFrame(summary).to_csv(outdir/"length_summary.tsv", sep="\t", index=False)

    datasets = {
        "Gene length": [gene_df.loc[gene_df.haplotype == h, "gene_length_bp"] for h in ("Hap1","Hap2")],
        "Exon length": [exon_df.loc[exon_df.haplotype == h, "length_bp"] for h in ("Hap1","Hap2")],
        "Intron length": [intron_df.loc[intron_df.haplotype == h, "length_bp"] for h in ("Hap1","Hap2")],
    }

    for title, data in datasets.items():
        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        make_boxplot(ax, data, title)
        fig.tight_layout()
        stem = title.lower().replace(" ", "_")
        fig.savefig(outdir/f"{stem}_Hap1_vs_Hap2.png", dpi=300, bbox_inches="tight")
        fig.savefig(outdir/f"{stem}_Hap1_vs_Hap2.pdf", bbox_inches="tight")
        plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (title, data) in zip(axes, datasets.items()):
        make_boxplot(ax, data, title)
    fig.tight_layout()
    fig.savefig(outdir/"boxplots_gene_exon_intron_lengths.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir/"boxplots_gene_exon_intron_lengths.pdf", bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    main()
