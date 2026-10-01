# Promoter ASM–ASE integration

`analysis/integration/phase_helpers.py` and `link_engine.py` retain the adopted haplotype-block phase, effect and transcript-support rules. `link.py` supplies explicit input paths. Integration reports block-level ASE effects.

```sh
python scripts/run_analysis.py integration --cpgs promoter_profiles.tsv.gz \
  --blocks reclassified_blocks.tsv.gz --sites RNA_sites.tsv.gz \
  --genes gene_status.tsv --annotation-index annotation_index --output NEW_LINK_DIRECTORY
```

The input tables must use the original field schemas. They are biological/intermediate data and are not embedded in the repository:

- `cpgs`: eligible promoter profiles containing `sample`, `tissue`, `chrom`, `cpg_pos1`, `ps`; H1/H2 methylated/unmethylated counts and depths; `delta_M`; `fisher_two_sided_p`; `source_BH_q_within_animal`, `source_q_eligibility`; annotation gene/name/profile IDs; `nominal_ge_3`, `nominal_ge_5`, `nominal_ge_8`, `nominal_ge_10` depth flags. The misleading historical `nominal_ge_*` names mean depth only, not P-value thresholds.
- `blocks`: original block results with v0.2.6 QC/classification columns, WGS phase status/orientation, balanced SNP membership, haplotype alleles, fragment counts and ASE P/q.
- `sites`: RNA SNP coverage evidence with `sample`, `tissue`, `chrom`, `position` (1-based), `variant_id` (`chrom_position_REF_ALT`), `all_parent_block_ids`, `matched_promoter_gene_ids` and `total_informative_fragments`. Membership fields contain comma-separated IDs. Keep observed fragment coverage and the original parent-block/gene mapping; a variant's presence in the WGS VCF alone does not establish RNA support. Extra fields in historical evidence tables are ignored; SNP-level ASE P/q values, calls and REF/ALT effect estimates are not consumed.
- `genes`: v0.2.6 gene status with sample/tissue/gene identifiers.

The join requires matching animal, tissue, chromosome and nonmissing exact WGS PS; technical QC; coherent DIRECT/FLIPPED SNP orientation; and compatible promoter/exonic transcript evidence. It emits `block_links.tsv`, `cpg_audit.tsv` and `block_attempts.tsv`. SNP coverage is used only to record transcript support for a block. Review `local_phase_status` and same-transcript support when selecting the final eligible comparisons. It does not silently pool disconnected or failed blocks.

RNA effect is `(H1-H2)/(H1+H2)` after WGS orientation; methylation effect is the H1 minus H2 methylation fraction. To express both against an ALT/REF anchor, flip **both** signs together when ALT is H2. The integration helper does not infer a shared allele label across animals without an anchor.

COX7C's six Heart observations came from the original full promoter linkage; THYN1's four Heart observations came from the later concordant-promoter selection. The THYN1 0377 second block failed technical/local-phase criteria and must not be pooled with the accepted block. These distinct source selections remain relevant even though plotting code is excluded.
