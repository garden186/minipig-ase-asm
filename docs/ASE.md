# ASE analysis

Run entry points as separate processes. Their native `common`/validator module names are local to each package. `scripts/run_analysis.py list` lists the principal entry points; each accepts `--help`.

## Dependencies and order

1. `ase_fragment_counts/call_phaser_ase_candidates.py`: joint phASER variant connections, haplotypes and haplotypic counts; matching phased WGS VCF; WGS-only phASER allelic counts; annotation GTF; RNA BAM manifest. The caller recounts distinct QNAME fragments instead of summing per-SNP counts.
2. `ase_reclassification/reclassify_ase_v026.py`: saved block results and original QC/validation metadata. Preserve technical conflict and local-phase QC fields.
3. `ase_cohort/build_ase_v026_cohort.py`: all reclassified animal outputs. Exclude Cranial as in the final cohort; expected study scope is 103 animal–tissue units and 11 tissues.
4. `ase_recurrence/core_eligible_recurrence.py`: complete eligible cohort and canonical cohort core/recurrent outputs, using the primary fragment threshold 15. The additional threshold options in the original source are not required for the primary release.
5. `ase_tissue_enrichment/ase_tissue_specificity.py`, followed by the primary pair selector below.

Core fragment settings retained from the original launcher include minimum total RNA fragments 15, each haplotype 3, WGS GQ20, WGS total 15, each WGS allele 5, and WGS ALT fraction 0.30–0.70. RNA MAPQ255, BQ10, maximum insert size 1,000,000, proper-pair requirements, duplicate/QC filtering and reverse-stranded counting are retained. Do not replace these with the defaults of an unrelated external tool.

## Final tissue-enriched pairs

```sh
python scripts/run_analysis.py ase-tissues \
  --recurrence-dir RECURRENCE_DIRECTORY --out-dir NEW_CONTRAST_DIRECTORY

python scripts/run_analysis.py ase-select-tissues \
  --contrasts NEW_CONTRAST_DIRECTORY/threshold_15_tissue_contrasts_all.tsv.gz \
  --metadata NEW_CONTRAST_DIRECTORY/analysis_metadata.json \
  --output primary_tissue_enriched_pairs.tsv
```

The wrapper explicitly applies threshold 15, ≥5 informative animals, ≥2 informative target-core animals and FDR 0.05. The selector uses pair-level global BH and primary recurrence, with no gene-omnibus gate. The adopted original result has 39 pairs and 37 genes. The recurrence result has 646 gene–tissue pairs and 439 genes. These numbers are regression expectations for the original inputs, not enforced outcomes for other datasets.

## SNP and PigGTEx analyses

`ase_snp_level` → `ase_snp_recurrence` → `ase_snp_loci`/`ase_snp_annotation` provides block-linked SNP support. `ase_gatk_calls/count_reads.py` → `ase_gatk_calls/call_gatk_snp_level_ase.py` → `ase_gatk_recurrence` provides the independent GATK route; `ase_gatk_blocks` checks correspondence to primary blocks.

`piggtex` handles gene/tissue enrichment; `piggtex_sites` handles exact position/allele matching. `piggtex_pair_support` summarizes independent SNP support for tissue-enriched pairs. `piggtex_multilayer` joins that support with primary pairs, the tissue map and PigGTEx gene profiles. Its `--inputs` JSON maps the four keys `primary_39`, `tissue_map`, `piggtex_gene_profiles`, `independent_snp_support` to input TSVs; `--output` supplies a new directory. Its original study-specific 39-pair checks remain intentional.

`piggtex_eqtl/prepare_full_universe.py` and `audit_full_gene_eqtl_finemap.py` handle the full candidate universe and external archives. `piggtex_evidence/build_evidence_audit.py` joins canonical recurrent SNP, GATK and external evidence. The argument name `server-*` in this preserved module identifies an input table role; it does not initiate a server connection.

External inputs must include the same PigGTEx release, tissue mapping, gene mapping, significant ASE sites, gene-eQTL and fine-mapping archives used in the study. These datasets are not included. Missing entries in significant-only resources mean no observed positive evidence, not a tested negative. Retain the complete tested gene background for enrichment.
