# BUSCO haplotype comparison

This directory contains the final outputs of the comparative BUSCO analysis between the two chromosome-scale *Africonus crotchii* haplotype assemblies.

The analysis was designed to compare BUSCO status between Hap1 and Hap2 and to investigate BUSCOs classified as **Missing in Hap1 but Complete in Hap2**. For these loci, whole-genome Hap1-to-Hap2 alignments were used to determine whether the corresponding Hap2 BUSCO locus, and its surrounding ±50 kb region, showed homologous sequence coverage in Hap1.

The BUSCO inputs were generated in genome mode using the same `mollusca_odb12` lineage dataset for both haplotypes.

## Files

```text
BUSCO_comparison/
├── matrix_hap1_hap2.tsv
├── BUSCO_locus_flank_classification.tsv
├── missing_to_complete_hap2_coords.tsv
├── BUSCO_categories_by_scaffold.tsv
├── recovered_buscos_by_contig.tsv
├── recovered_buscos_density_perMb.tsv
├── hap2_to_hap1_scaffold.tsv
└── busco_relevant_alignments.paf.gz
```

### `matrix_hap1_hap2.tsv`

Transition matrix comparing BUSCO categories between Hap1 and Hap2.

BUSCOs are classified as Complete, Duplicated, Fragmented, or Missing in each haplotype, allowing direct identification of category changes such as Complete → Complete, Fragmented → Complete, Missing → Complete, and Complete → Missing.

### `BUSCO_locus_flank_classification.tsv`

Final classification of BUSCOs that are **Missing in Hap1 and Complete in Hap2**.

For each BUSCO, the analysis evaluates Hap1-to-Hap2 whole-genome alignments over the BUSCO locus in Hap2 and the BUSCO locus plus 50 kb of flanking sequence on each side.

This is the main output used to distinguish BUSCO differences associated with missing homologous assembled sequence from cases where homologous sequence is present, but BUSCO was not detected in Hap1.

### `missing_to_complete_hap2_coords.tsv`

Coordinates in the Hap2 assembly of BUSCOs classified as Missing in Hap1 and Complete in Hap2.

### `BUSCO_categories_by_scaffold.tsv`

Distribution of BUSCO categories across Hap2 scaffolds/pseudochromosomes.

### `recovered_buscos_by_contig.tsv`

Counts of BUSCOs recovered in Hap2, relative to Hap1, summarized by scaffold/contig.

### `recovered_buscos_density_perMb.tsv`

Density of recovered BUSCOs per megabase for each scaffold/contig, allowing comparison among scaffolds of different sizes.

### `hap2_to_hap1_scaffold.tsv`

Correspondence table between homologous Hap2 and Hap1 scaffolds/pseudochromosomes. We used this mapping to restrict BUSCO-locus comparisons to the expected homologous chromosome pair.

## Analysis logic

```text
BUSCO Hap1 + BUSCO Hap2
          |
          v
  BUSCO status comparison
          |
          v
 Identify Hap1 Missing /
      Hap2 Complete
          |
          v
 Obtain Hap2 BUSCO coordinates
          |
          v
 Compare with Hap1 → Hap2
 whole-genome alignments
          |
          +----------------------+
          |                      |
          v                      v
    BUSCO locus             ±50 kb flanks
          |                      |
          +----------+-----------+
                     |
                     v
       Locus/flank classification
```
