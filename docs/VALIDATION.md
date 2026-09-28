# Validation scope

This release separates unchanged scientific kernels from new execution interfaces. File identities and source lineage are recorded in the release/provenance manifests.

Local checks cover Python parsing/importable command help, synthetic preprocessing plans, shell syntax, original ASE tests, single-CpG fragment/evidence tests, exact statistical kernels, and frozen-data ASM regression suites. The machine-readable public summary is `provenance/validation_summary.json`. Private biological fixtures, absolute paths and detailed study data are not included in that summary.

## What has not been established

- A complete Linux FASTQ-to-final-result rerun of the new interfaces.
- Historical execution identity for every preprocessing software version or server-modified wrapper.
- Automatic tested-row annotation of the new individual ASM output and complete construction of downstream native stage configurations. The original annotated-table-to-cache generator is now included and its output was consumed by the original anchor indexer in a synthetic integration test.
- Public availability of all external data, frozen regression inputs and PigGTEx/GO annotations.
- Reproduction of every final manuscript table by this new checkout, or a stable archival DOI.

The native tests retain source-hash, exact-phase, count arithmetic, full-family multiple-testing and checkpoint checks. Their success supports the tested scientific functions and fixtures; it does not certify unexecuted full-genome interfaces. A test skipped for missing SAMtools/BCFtools is recorded as skipped, not passed.

`0.2.0-rc2` remains a code release candidate. The rc1 validation summary is retained as the original validation record; `provenance/asm_cache_validation.json` records the five additional cache tests. No unchanged scientific package identity or historical test record was rewritten. Do not cite this as a complete end-to-end reproduction until the remaining interface/data/environment boundaries have been verified.
