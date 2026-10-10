# Tested-CpG cache preparation

`analysis/asm_cache/cache.py` and `common.py` are the original cache implementation from `asm_inference_prep_v0.1.0`. The final common-anchor run consumed these caches. Their scientific code is unchanged; the new `run.py` supplies a relocatable configuration interface. The historical inventory and exploratory statistical pilot are not included or executed. Historical policy descriptions in the original shared module describe its source package, not new inference by this cache-only command.

```sh
python scripts/run_analysis.py asm-cache --config cache_inputs.json --output new_cache
```

The JSON configuration has a `units` array. Each item has `sample`, `tissue` and `annotated`. `annotated` contains `sample`, `tissue`, `path`, the actual compressed-file `sha256`, and integer `rows`. If an original annotation validation record exists, also supply `gate` and `gate_sha256`; these are verified. Paths resolve beside the configuration. Output locations are supplied at execution, with no author-specific server paths.

Each source is a gzip-compressed TSV of all eligible tested CpG/PS rows for that unit, **including non-significant rows**. Required columns are:

- `sample`, `tissue`, `chrom`, `cpg_pos1`, `ps`, `genotype_disrupted`, `test_status`;
- `H1_methylated`, `H1_unmethylated`, `H2_methylated`, `H2_unmethylated`, `H1_depth`, `H2_depth`;
- `H1_methylation`, `H2_methylation`, `delta_M`, `fisher_two_sided_p`;
- `BH_q_within_animal`, `BH_q_tissue10`, `BH_q_all30`;
- at least one original `annotation_` column. All columns are retained in the source sidecars.

The original checks require `genotype_disrupted=0`, `test_status=TESTED`, autosomes 1–18, each haplotype depth ≥3, combined depth ≥10, consistent counts/effects, valid P/q, and unique exact CpG/PS keys. P/q are copied, not recalculated. Exact PS strings remain distinct, including `001` versus `1`.

Outputs include chromosome arrays, the phase dictionary, original-column sidecars, validation hashes and `cache_registry.json`. Its unit-level `annotated` and `cache` objects supply those roles in the native anchor configuration. Full raw opportunity tables, WGS phase inputs, reference inputs and the remaining anchor configuration fields must still be supplied. The final anchor controller expects all 30 units and all 18 autosomes; a partial cache run does not constitute the final analysis.

The new individual interface also emits untested rows; do not pass that output directly as a tested-only source or invent annotation columns. Complete the tested-row annotation step without significance/effect selection first. The automatic adapter from that new interface is not supplied here.

The five synthetic tests include preservation of counts/P/q/exact PS, duplicate and arithmetic rejection, changed-source rejection on restart, path relocation, and actual consumption of the exported cache by the original anchor availability indexer. They do not substitute for a full study rerun.
