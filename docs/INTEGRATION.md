# Promoter ASM–ASE integration

`analysis/integration/phase_helpers.py` and `link_engine.py` contain the original computational functions extracted from the adopted source chain. `link.py` supplies explicit paths instead of looking in dated author directories.

```sh
python scripts/run_analysis.py integration --cpgs promoter_profiles.tsv.gz \
  --blocks reclassified_blocks.tsv.gz --sites RNA_sites.tsv.gz \
  --genes gene_status.tsv --annotation-index annotation_index --output NEW_LINK_DIRECTORY
```

The input tables must use the original field schemas. They are biological/intermediate data and are not embedded in the repository:

- `cpgs`: eligible promoter profiles containing `sample`, `tissue`, `chrom`, `cpg_pos1`, `ps`; H1/H2 methylated/unmethylated counts and depths; `delta_M`; `fisher_two_sided_p`; `source_BH_q_within_animal`, `source_q_eligibility`; annotation gene/name/profile IDs; `nominal_ge_3`, `nominal_ge_5`, `nominal_ge_8`, `nominal_ge_10` depth flags. The misleading historical `nominal_ge_*` names mean depth only, not P-value thresholds.
- `blocks`: original block results with v0.2.6 QC/classification columns, WGS phase status/orientation, balanced SNP membership, haplotype alleles, fragment counts and ASE P/q.
- `sites`: original site-level fragment results including `all_parent_block_ids`, `matched_promoter_gene_ids`, REF/ALT counts, exact variant IDs and site-level ASE fields.
- `genes`: v0.2.6 gene status with sample/tissue/gene identifiers.

The join requires matching animal, tissue, chromosome and nonmissing exact WGS PS; technical QC; coherent DIRECT/FLIPPED SNP orientation; and compatible promoter/exonic transcript evidence. It emits block links, SNP links, CpG audit and rejected/accepted block attempts. Review `local_phase_status` and same-transcript support when selecting the final eligible comparisons. It does not silently pool disconnected or failed blocks.

RNA effect is `(H1-H2)/(H1+H2)` after WGS orientation; methylation effect is the H1 minus H2 methylation fraction. To express both against an ALT/REF anchor, flip **both** signs together when ALT is H2. The integration helper does not infer a shared allele label across animals without an anchor.

COX7C's six Heart observations came from the original full promoter linkage; THYN1's four Heart observations came from the later concordant-promoter selection. The THYN1 0377 second block failed technical/local-phase criteria and must not be pooled with the accepted block. These distinct source selections remain relevant even though plotting code is excluded.
