# Stage 1: WGS preprocessing and local phasing

## Purpose

Stage 1 converts paired-end minipig WGS reads into a quality-controlled autosomal heterozygous-SNV set and a locally phased VCF. The phased VCF is the genotype backbone used by downstream ASE and ASM workflows.

## Workflow and retained analysis choices

| Step | Script | Main choices |
|---|---|---|
| Read QC and trimming | `run_fastp_wgs.sh` | explicit TruSeq adapters; right-end quality trimming; minimum length 50 bp; poly-G trimming |
| Alignment and QC | `wgs_align.sh` | BWA-MEM `-K 100000000 -Y`; coordinate sorting; GATK MarkDuplicates; mosdepth; flagstat |
| Phase-ready BAM | `wgs_align.sh` | proper pairs, MAPQ >= 30; excludes unmapped, secondary, QC-failed, duplicate, and supplementary records |
| Variant calling | `run_deepvariant.sh` | DeepVariant 1.10.0 WGS model; autosomes 1-18 |
| Het-SNV preparation | `prepare_het_vcf.sh` | autosomal, PASS, biallelic SNV, heterozygous, DP >= 10; multiallelic sites excluded |
| Phasing | `run_whatshap_wgs_honest.sh` | WhatsHap 2.8; SNVs only; local read-backed phasing; PS tags; MAPQ >= 30 |

The suffix `honest` is retained for compatibility with existing downstream files. It is an internal label for the stringent read-filtered baseline and is not a statistical guarantee.

## Required inputs

1. Paired-end WGS FASTQ files listed in `config/samples.tsv`.
2. Sscrofa11.1 reference FASTA with:
   - `samtools faidx` index (`.fai`)
   - BWA index (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`)
   - sequence dictionary (`.dict`)
3. A working Docker installation for the DeepVariant image.

The sample sheet is deliberately local and git-ignored because real file paths can reveal server layout or controlled-data locations.

## Main outputs

| Output | Role | Retention |
|---|---|---|
| `results/wgs/trim/{sample}/{sample}_R{1,2}.fastq.gz` | trimmed WGS reads | reproducibility intermediate |
| `results/wgs/align/{sample}.dedup.bam` | coordinate-sorted duplicate-marked WGS alignment | important intermediate |
| `results/wgs/align/{sample}.dedup.phaseready.bam` | high-quality input used by WhatsHap | archive; required only to rerun phasing |
| `results/wgs/variants/{sample}.deepvariant.vcf.gz` | raw autosomal DeepVariant calls | important intermediate |
| `results/wgs/variants/{sample}.het.snp.vcf.gz` | primary Stage 1 heterozygous-SNV set | downstream input |
| `results/phasing/{sample}.phased.wgs.honest.vcf.gz` | locally phased autosomal SNVs | primary deliverable |
| `*.stats.tsv`, `*.blocks.tsv`, `*.summary.tsv` | phasing QC and block summaries | manuscript/QC evidence |
| `*.used_reads.txt` | WhatsHap-selected phase-informative read entries | provenance/QC |

The phase-ready BAM does not need to remain on the compute server after successful phasing, provided the BAM and index are archived with checksums. It must be restored only when phasing is rerun or read-level phasing provenance is re-examined.

## Restart and repair behavior

Every computational step validates its own outputs and skips only complete output sets. The orchestrator does not use the phased VCF alone as a completion sentinel.

If a phased VCF already exists but `stats.tsv`, `blocks.tsv`, or `summary.tsv` is missing, rerun:

```bash
bash scripts/stage1/run_whatshap_wgs_honest.sh \
  --sample SAMPLE_001 --config config/stage1.env --summary-only
```

This reconstructs the QC files from the phased VCF and does not require the phase-ready BAM. It is also the intended repair path for a historical header-only summary.

## QC definitions

- `mean_autosome_depth` is length-weighted across chromosomes 1-18. It is not the `total` row from mosdepth, which can include unplaced scaffolds.
- Coverage breadth is calculated across autosomal chromosome lengths from the mosdepth cumulative distribution.
- `phased_fraction` is the fraction of heterozygous SNVs represented with phased genotypes; it does not measure chromosome-scale continuity.
- `block_span_n50_bp` is computed over spans in WhatsHap `blocks.tsv`. WhatsHap's reported `block_n50` can instead be genome-length-based NG50 and may be zero when blocks cover less than half of the supplied genome length.
- WhatsHap `blocks` excludes singleton phase sets, while the block-list file can include singleton rows. Both are reported explicitly.

The audited pilot run yielded approximately 40.7x autosomal mean depth, 8.19 million primary heterozygous SNVs, 98.04% phased SNVs, and a block-span N50 of 5,779 bp. These values are a reference for detecting gross drift, not hard-coded acceptance criteria for every sample.

## Interpretation limits

### Cross-species variant calling

DeepVariant's WGS model is trained on human truth data. The pipeline records its use exactly because it generated the analyzed cohort, but pig genotypes should not be described as species-benchmarked calls. Sensitivity analyses using genotype quality and allele balance, comparison with another caller, or orthogonal validation strengthen claims based on rare or sample-specific sites.

### Short-read phasing

The phased VCF is locally read-backed. `0|1` and `1|0` are meaningful only within the same phase set (`PS`). Orientation must not be propagated across different phase sets, and the result must not be described as parental or chromosome-resolved haplotypes.

### Primary versus sensitivity filters

To reproduce existing results, the primary het-SNV output applies DP >= 10 without a hidden GQ or allele-balance cutoff. `prepare_het_vcf.sh` writes GQ and allele-balance QC summaries. Any GQ >= 20/30 or allele-balance sensitivity set should be generated under a distinct filename and reported as a sensitivity analysis, not substituted silently for the primary set.

## Version provenance

The orchestrator writes `logs/{sample}/stage1_software_versions.txt`. The audited materials directly establish fastp 1.3.3, WhatsHap 2.8, and DeepVariant 1.10.0. Other exact package builds should be taken from the run-specific version file or environment export rather than inferred retrospectively.
