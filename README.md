# Minipig ASE/ASM

Reproducible workflows for integrative allele-specific expression (ASE) and allele-specific methylation (ASM) analyses in minipigs.

The repository is being released stage by stage. The current version contains the finalized Stage 1 workflow for autosomal WGS preprocessing, germline variant calling, and local read-backed phasing.

## Stage 1 at a glance

```text
paired-end WGS FASTQ
  -> fastp
  -> BWA-MEM
  -> GATK MarkDuplicates
  -> mosdepth / samtools QC
  -> DeepVariant
  -> PASS biallelic heterozygous SNVs (DP >= 10)
  -> WhatsHap local read-backed phasing
```

Stage 1 intentionally retains the filtering and phasing settings used for the analyzed cohort. Paths and sample-specific file names are supplied through local configuration files and are not embedded in the code.

## Quick start

1. Create the software environment.

   ```bash
   mamba env create -f envs/stage1.yml
   mamba activate minipig-stage1
   ```

2. Copy and edit the local configuration files.

   ```bash
   cp config/stage1.env.example config/stage1.env
   cp config/samples.tsv.example config/samples.tsv
   ```

3. Run one sample.

   ```bash
   bash scripts/stage1/stage1_wgs_preprocessing.sh \
     --sample 0235 \
     --config config/stage1.env \
     --threads 16
   ```

4. Validate an existing or completed result set.

   ```bash
   bash scripts/stage1/validate_stage1_outputs.sh \
     --sample 0235 \
     --config config/stage1.env
   ```

Detailed inputs, outputs, restart behavior, QC definitions, and interpretation limits are in [docs/stage1.md](docs/stage1.md).

## Cohort QC and manuscript outputs

After every sample has passed Stage 1 validation, build the concise cohort table:

```bash
python3 scripts/stage1/collect_stage1_cohort_qc.py \
  --project-dir /path/to/minipig_ase_asm \
  --sample-sheet /path/to/minipig_ase_asm/config/samples.tsv \
  --build-missing-dp-qc \
  --output /path/to/minipig_ase_asm/results/qc/stage1_cohort_qc.tsv
```

Generate the two supplementary figures, reporting table, data dictionary, descriptive statistics, and legends:

```bash
mamba env create -f envs/stage1-reporting.yml
mamba activate minipig-stage1-reporting
python3 scripts/stage1/make_stage1_manuscript_outputs.py \
  --cohort-qc /path/to/minipig_ase_asm/results/qc/stage1_cohort_qc.tsv \
  --output-dir /path/to/minipig_ase_asm/results/stage1_reporting
```

## Repository layout

```text
config/          local configuration templates
docs/            workflow and interpretation documentation
envs/            software environment definitions
scripts/stage1/  Stage 1 executable scripts
tests/           lightweight code and summary tests
```

## Scope and limitations

- Only Sus scrofa autosomes 1-18 on Sscrofa11.1 are processed in Stage 1.
- The DeepVariant WGS model is human-trained. Its use in minipig is a cross-species application and is not equivalent to a species-specific benchmark.
- Short-read WhatsHap results represent local phase sets. Phase orientation is not transferable between phase sets, and the output is not chromosome-scale or parental haplotype resolution.
- The primary heterozygous-SNV set reproduces the cohort analysis (`PASS`, biallelic SNV, heterozygous genotype, `DP >= 10`). GQ and allele-balance distributions are reported for sensitivity assessment rather than silently changing the primary call set.

## Status

Stage 1 is operationally finalized and reproduces the validated 0235 workflow. Later ASE, ASM, integration, and annotation stages will be added after their corresponding output audits are complete.

## License

No reuse license has been assigned yet. A license should be selected after institutional and data-governance review.
