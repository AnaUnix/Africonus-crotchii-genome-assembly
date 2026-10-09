# Genome assembly workflow

## Input

- Sample: MUS4382
- HiFi reads: 2,931,519
- Yield: 55.4 Gb
- Genome-size estimate: 3.5 Gb
- Approximate reported coverage: 16x

## Workflow

```text
PacBio CCS / HiFi reads
        |
        v
Hifiasm v0.15.4-r347
(default parameters)
        |
        v
Hifiasm output assembly
        |
        +--> BLAST against nt
        |        |
        |        v
        |   BlobTools v1.1.1
        |        |
        |        v
        |   removal of possible contamination
        |
        v
haplotype resolved contig aasembly
```

## Haplotype-specific pre-HiRise assemblies

Separate Hap1 and Hap2 assemblies were subsequently supplied to Dovetail HiRise:

- Hap1 input: 2,274,659,305 bp; N50 168,081 bp; 19,653 sequences.
- Hap2 input: 3,173,374,471 bp; N50 241,767 bp; 19,824 sequences.
