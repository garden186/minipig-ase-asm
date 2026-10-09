# Adopted ASM source lineage

A source audit on 2026-09-28 traced the adopted MAPQ0 results using locally received server execution and validation records. This was a provenance audit, not a new whole-genome analysis or live server inspection.

| Stage | Original implementation | Evidence |
|---|---|---|
| Tested-count cache | September 17 `asm_inference_prep_v0.1.0/cache.py` | Local package manifest matches the executed server release; all 30 cache gate hashes and 600 registered payload identities match the subsequent anchor inputs |
| Common-anchor preparation | September 24 `asm_allele_prep_v0.1.0` | Completed source gate, run identity and source configuration match the directional inputs |
| Directional recurrence | September 24 `asm_directional_v0.1.0` | Completed 54-family inference; historical r2 totals Heart 977, Kidney 1,084, Liver 1,451, union 2,197 |
| Paired tissue comparison | September 25 `asm_tissue_aligned_v0.1.0` | Completed adopted analysis; upstream identity, six root metadata files and all 236 tile registrations match |

The cache source covers 59,071,985 tested phase observations. Large chromosome arrays are not included in the compact received receipts; the audit compared their recorded hashes across stages, not all array bytes. Existing independent final-result reviews recorded passing completion gates for these stages.

The initial public curation omitted the original cache generator. The `asm_cache` addition corrects that omission. Historical execution success does not imply that every newly written public interface has completed a full study rerun. Remaining interface requirements are documented in [ASM analysis](ASM.md).
