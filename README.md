# Minipig ASE and ASM analysis

Scripts for sequencing data preprocessing and allele-specific expression (ASE) and methylation (ASM) analyses in Korean minipigs.

## Analyses

- WGS, RNA-seq and WGBS preprocessing
- ASE: haplotype counting, recurrence, tissue enrichment and gene-level PigGTEx comparisons
- ASM: single-CpG testing, common-anchor alignment, directional recurrence and tissue comparisons
- Gene annotation, GO enrichment and ASM–ASE integration

## Reference

- FASTA: `Sus_scrofa.Sscrofa11.1.dna.toplevel.fa`
- GTF: `Sus_scrofa.Sscrofa11.1.115.gtf`

## Usage

Configure input paths using the examples in `config/`. See the [environment requirements](docs/ENVIRONMENTS.md) for installation and external tools.

```sh
python -m pip install -r environments/requirements-portable-tested.txt
python scripts/run_analysis.py list
python scripts/run_analysis.py ase-count --help
```

Instructions: [Preprocessing](docs/PREPROCESSING.md) · [ASE](docs/ASE.md) · [ASM](docs/ASM.md) · [Integration](docs/INTEGRATION.md) · [Validation status](docs/VALIDATION.md)
