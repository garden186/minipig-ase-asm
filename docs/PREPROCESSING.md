# Preprocessing

Use `preprocessing/run.py` on Linux. All stages print their exact commands by default; `--execute` runs them and records the plan and logs. There are no raw data downloads or hidden server connections.

Copy `config/preprocessing.example.json` and `config/samples.example.tsv` to your own configuration directory. Paths in the JSON resolve relative to that JSON; FASTQ paths resolve relative to the sample sheet. `tools` may map executable names to installed binaries. Tool arguments remain separate from executable paths.

The sample sheet has one row per library with `sample`, `assay` (`WGS`, `RNA`, `WGBS`), `tissue`, `read1`, `read2`, and `optical_duplicate_distance`. Use one WGS library per animal. Combine lanes before this interface, as appropriate for the original library; do not create duplicate sample/assay/tissue keys. Use normalized tissue names consistently.

The optical duplicate distance is explicit for WGS: the historical script inferred 100 or 2500 from the instrument prefix. The example's 2500 is illustrative, not authenticated per-library metadata. RNA duplicate marking retains the original 2500 setting.

```sh
python preprocessing/run.py reference --config /my/config/preprocessing.json --samples /my/config/samples.tsv
python preprocessing/run.py wgs --config /my/config/preprocessing.json --samples /my/config/samples.tsv --sample S01
python preprocessing/run.py rna --config /my/config/preprocessing.json --samples /my/config/samples.tsv --sample S01
python preprocessing/run.py phaser --config /my/config/preprocessing.json --samples /my/config/samples.tsv --sample S01
python preprocessing/run.py wgbs --config /my/config/preprocessing.json --samples /my/config/samples.tsv --sample S01
```

After reviewing commands, append `--execute`. Build the shared reference once, then run animals separately in this order. Start with a new work directory. A stage log directory cannot be reused automatically; inspect failed logs before choosing a new destination. These interfaces have not been validated by a complete Linux NGS rerun.

## Outputs and connections

- WGS writes `work/wgs/SAMPLE/{marked.bam,phase_ready.bam,variants.vcf.gz,heterozygous.vcf.gz,phased.vcf.gz}` plus QC.
- RNA writes `work/rna/SAMPLE/TISSUE/phase_ready.bam`. First-pass junctions are combined across the animal's RNA tissues for the second pass. STAR/WASP `vW=2..7` alignments and secondary/supplementary alignments are removed; downstream counting applies additional original filters.
- phASER writes joint WGS/RNA results under `work/phaser/SAMPLE/rna.*` and WGS-only results under `wgs.*`. WGS contributes to phasing but is excluded from RNA haplotypic counts. WGS-only allele counts are required by custom ASE WGS-balance QC. The final pipeline does not use `phaser_gene_ae` as its gene-count endpoint.
- WGBS writes the per-animal allele list and N-masked genome, then `work/wgbs/SAMPLE/TISSUE/asm_input.bam` and its index. This is the unsplit deduplicated BAM for custom CpG counting. SNPsplit BAM separation and per-haplotype methylation extraction are not part of this final path.

These are **new interface output names**, not claims about historical server filenames. The portable WGBS interface explicitly coordinate-sorts the deduplicated BAM; recovering the old `.srt.bam` filename or original sorting shell command is unnecessary.

## Independent GATK counts

Create a TSV with `sample`, `tissue`, `bam`, `vcf`. Supply phase-ready RNA BAMs and the matching heterozygous WGS SNV VCFs. Paths in this interface are explicit paths (relative to the working directory or absolute).

```sh
python analysis/ase_gatk_calls/count_reads.py --manifest gatk_inputs.tsv --reference reference.fa --output new_gatk_run
```

Append `--execute` after command review. This reproduces the original BQ10/MAPQ255/proper-pair/fragment-length and overlapping-mate counting options. Its `input_manifest.tsv` and `counts/` feed `ase_gatk_calls/call_gatk_snp_level_ase.py`. Statistical calling remains a separate step.
