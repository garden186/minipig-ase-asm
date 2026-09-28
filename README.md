# Minipig ASE and ASM analysis scripts

Analysis code for the Korean minipig study, reconstructed from the local analysis sources and adopted final-result implementations. This repository contains preprocessing, statistical analysis and table-join code.

**Release: 0.2.0-rc2, a public code curation candidate.** Existing scientific kernels are retained wherever possible. New execution interfaces have been added. This is a collection of analysis scripts with documented input contracts, **not yet a verified single-command FASTQ-to-manuscript workflow**. See [validation and remaining boundaries](docs/VALIDATION.md).

## Reference inputs

- Assembly FASTA: `Sus_scrofa.Sscrofa11.1.dna.toplevel.fa`
- Gene annotation: `Sus_scrofa.Sscrofa11.1.115.gtf` (Ensembl release 115)
- Autosomes: `1` through `18`; use consistent chromosome names throughout.

FASTA, GTF and downloaded PigGTEx/GO resources are external inputs. Filenames identify the requested references but do not prove file identity. The annotation engine also records and checks a GTF SHA256. No reference sequence or annotation data are redistributed.

## Analysis map

| Topic | Included code |
|---|---|
| WGS | fastp, BWA-MEM, duplicate marking, DeepVariant, heterozygous-SNV filtering, WhatsHap, QC |
| RNA-seq | Trim Galore, two-pass STAR/WASP, duplicate marking and phase-ready BAMs |
| phASER | Joint WGS/RNA phasing and WGS-only allele counting required for WGS balance QC |
| WGBS | N-masked references, trimming, Bismark, deduplication and coordinate sorting |
| ASE | Unique-fragment counting, local block QC, final classification and cohort aggregation |
| ASE follow-up | Recurrence, tissue enrichment, SNP-level analyses, independent GATK validation |
| PigGTEx | Matched-tissue enrichment, exact ASE-site support, gene-eQTL and fine-mapping joins |
| ASM | Original single-CpG fragment extraction, GQ20 reassignment, individual exact tests and complete BH families |
| ASM follow-up | Original tested-count cache generation, common-anchor preparation, directional recurrence, paired tissue comparisons |
| Annotation/GO | Ensembl115 transcript/promoter annotation and native clusterProfiler BP enrichment |
| ASM–ASE integration | Same-animal/tissue/phase-set and same-transcript joins, signed effects |

```sh
python scripts/run_analysis.py list
python scripts/run_analysis.py ase-count --help
python scripts/run_analysis.py asm-count --help
```

## Start here

1. Use Linux for NGS processing. Install the external tools described in [environment notes](docs/ENVIRONMENTS.md).
2. Copy the example configuration and sample sheet outside the checkout. Set your own input and work locations.
3. Review [preprocessing](docs/PREPROCESSING.md), then [ASE analysis](docs/ASE.md) and [ASM analysis](docs/ASM.md).
4. Keep biological inputs and outputs outside this repository. Repository-relative imports are intentional; execution does not depend on where the repository is installed.

```sh
python -m pip install -r environments/requirements-portable-tested.txt
python scripts/verify_release.py
python scripts/run_tests.py --output /your/new/test-directory
```

The default test run does not require study datasets. Some tests require installed NGS tools and report a skip when unavailable. Additional frozen-data ASM regression tests use private external assets; their requirements are documented separately.

## Final-analysis contracts

The adopted ASM analysis uses minimum MAPQ **0**, while retaining the other original read/fragment filters and assignment GQ ≥ 20. Only single-CpG final analysis entry points are documented. Historical helper filenames containing `cpgk` are retained because the single-CpG engine imports their parsing/fragment functions; they are not an additional CpG4 analysis workflow.

Final tissue-enriched ASE selection uses ≥5 informative animals and pair-level global BH q ≤ 0.05 together with primary recurrence. The separate selector does not impose a gene-omnibus gate. Use the final commands in the documentation, rather than assuming every historical module default is the final primary setting.

H1/H2 labels refer to within-animal local phase sets. They do not establish parental origin or globally corresponding haplotypes. Common-anchor alignment is required before cross-animal signed ASM comparisons. Absence from a significant-only PigGTEx resource is not a tested negative.

Source lineage is recorded in `provenance/source_curation.json`; current file identities are fixed by `release_manifest.json`. Software licence terms have not been selected by the authors; no third-party tools or data are redistributed under an invented licence.
