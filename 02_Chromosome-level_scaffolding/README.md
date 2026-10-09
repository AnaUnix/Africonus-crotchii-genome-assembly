# Chromosome-level scaffolding

This directory documents the Omni-C / Dovetail HiRise scaffolding of the two *Africonus crotchii* haplotype assemblies.

## Summary

HiRise increased scaffold N50 from 168,081 bp to 60,234,464 bp in Hap1 and from 241,767 bp to 78,584,641 bp in Hap2. The final project chromosome sets comprise 35 pseudochromosomes per haplotype.

## Omni-C / HiRise scaffolding

Chromosome-scale scaffolding was performed using Omni-C proximity-ligation data and the HiRise workflow.

The two haplotype-specific pre-HiRise assemblies were scaffolded separately.

Omni-C library reads were aligned to the draft assemblies using BWA. Read-pair separations within draft scaffolds were used by HiRise to build a likelihood model of genomic distance. HiRise then used this model to identify and break putative misjoins, score prospective joins, and make joins above a threshold.

Reads with mapping quality >50 were used for scaffolding.

### Hap1

- Input assembly length: 2,274,659,305 bp
- HiRise assembly length: 2,276,477,405 bp
- Input N50: 168,081 bp
- HiRise N50: 60,234,464 bp
- Input scaffolds: 19,653
- HiRise scaffolds: 3,491
- Largest HiRise scaffold: 114,775,838 bp
- Joins made: 18,181
- Breaks made: 2,019
- Read-pairs: 38,322,733

### Hap2

- Input assembly length: 3,173,374,471 bp
- HiRise assembly length: 3,175,209,571 bp
- Input N50: 241,767 bp
- HiRise N50: 78,584,641 bp
- Input scaffolds: 19,824
- HiRise scaffolds: 3,873
- Largest HiRise scaffold: 147,317,800 bp
- Joins made: 18,351
- Breaks made: 2,400
- Read-pairs: 50,388,525

The increase in assembly length exactly matches 100 bp per reported HiRise join in both haplotypes, consistent with the insertion of fixed-length gap sequence during scaffolding.

<img width="2641" height="88" alt="image" src="https://github.com/user-attachments/assets/0a41ea95-2ccb-4042-982e-0c12e620a099" />

## Pseudochromosomes

The final project chromosome set contains 35 pseudochromosomes for each haplotype. Scaffold_1–Scaffold_35 identifiers were used for the chromosome-scale assemblies.

