# ASE analysis

Run entry points as separate processes. Their native `common`/validator module names are local to each package. `scripts/run_analysis.py list` lists the principal entry points; each accepts `--help`.

## Dependencies and order

1. `ase_fragment_counts/call_phaser_ase_candidates.py`: joint phASER variant connections, haplotypes and haplotypic counts; matching phased WGS VCF; WGS-only phASER allelic counts; annotation GTF; RNA BAM manifest. The caller recounts distinct QNAME fragments instead of summing per-SNP counts.
2. `ase_reclassification/reclassify_ase.py`: saved block results and original QC/validation metadata. Preserve technical conflict and local-phase QC fields.
3. `ase_cohort/build_ase_cohort.py`: all reclassified animal outputs. Exclude Cranial as in the final cohort; expected study scope is 103 animal–tissue units and 11 tissues.
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

## PigGTEx gene-level comparisons

`piggtex/piggtex_validation.py` handles matched-tissue gene support and enrichment against PigGTEx ASE and cis-eQTL resources, using the complete internal recurrence hypothesis universe as the background.

External inputs must include the study's PigGTEx release, tissue and gene mappings, significant ASE catalogue and permutation cis-eQTL tables. These datasets are not included. Missing entries in significant-only resources mean no observed positive evidence. Retain the complete tested gene background for enrichment.

## Analysis scope

The public workflow uses haplotype-block ASE and gene-level PigGTEx comparisons. Standalone SNP-level ASE analyses are outside its scope. WGS variant preparation, SNP balance QC, fragment assignment and phase/transcript evidence remain necessary inputs to haplotype ASE and ASM.
