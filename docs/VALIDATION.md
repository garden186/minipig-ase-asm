# Validation scope

This release separates unchanged scientific kernels from new execution interfaces. File identities and source lineage are recorded in the release/provenance manifests.

Local checks cover Python parsing/importable command help, synthetic preprocessing plans, shell syntax, original ASE tests, single-CpG fragment/evidence tests, exact statistical kernels, and frozen-data ASM regression suites. The machine-readable public summary is `provenance/validation_summary.json`. Private biological fixtures, absolute paths and detailed study data are not included in that summary.

## What has not been established

- A complete Linux FASTQ-to-final-result rerun of the new interfaces.
- Historical execution identity for every preprocessing software version or server-modified wrapper.
- Automatic conversion of the new individual ASM tables to the native anchor/recurrence controllers' cache and receipt formats.
- Public availability of all external data, frozen regression inputs and PigGTEx/GO annotations.
- Reproduction of every final manuscript table by this new checkout, or a stable archival DOI.

The native tests retain source-hash, exact-phase, count arithmetic, full-family multiple-testing and checkpoint checks. Their success supports the tested scientific functions and fixtures; it does not certify unexecuted full-genome interfaces. A test skipped for missing SAMtools/BCFtools is recorded as skipped, not passed.

`0.2.0-rc1` is therefore a code release candidate. Do not cite it as a complete end-to-end reproduction until the remaining interface/data/environment boundaries have been verified.
