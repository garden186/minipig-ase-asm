# Native ASM configuration contracts

The native controllers were designed around authenticated intermediate files. Their path configuration is now external, but their data formats and validation records remain strict. This file describes the required roles; it is not a substitute for complete source manifests.

| Stage | Main JSON roles |
|---|---|
| `asm_allele_prep` | `ids`, `tissues`, `chromosomes`; `units` with `sample`, `tissue`, `raw`, `annotated`, `cache`; `wgs`; `reference`; `frozen_anchors`, `frozen_profiles`; `expected`; `pilot_per_stratum`; `provenance` |
| `asm_directional` | `stage1_root`, `stage1_identity`, stage-1 validation identities; `locations`, `registry`, `source_inputs`; `tile_bp`; expected full/pilot counts and regression identities |
| `asm_tissue_comparison` | `stage2_root`, `stage2_identity`; `root_checks`, `tiles`, `fixture`, `expected`, `policy` |

Exact key usage is visible in each package's `run.py`, `pipeline.py`/`pilot.py`, and `indexing.py`/`data.py`. Paths, checksums, row counts and validation identities must all describe the supplied inputs. Most source-file descriptors use `path` and `sha256`. Cache directories also contain `phase_dictionary.json`, per-chromosome NumPy arrays and a validation receipt.

The original compact tested-count cache has chromosome, CpG position, exact-PS dictionary index, four H1/H2 methylated/unmethylated counts, original P/global q and source-row identity. Availability additionally includes raw CpG/PS opportunities. Both are required for P-independent anchor ranking; supplying only significant or testable sites loses the original selection background.

The `selftest.py` files generate synthetic small execution fixtures in temporary directories, and additionally use the original checksum-registered regression assets. Those biological fixtures are not included in this code-only repository. Default `scripts/run_tests.py` runs the data-free suites; native regression tests can be added with `--asm-config PRIVATE_TEST_CONFIG.json`, where each of `asm_allele_prep`, `asm_directional`, `asm_tissue_comparison` maps to `assets` and `inputs` paths.

A selftest record must match the exact native package `release.json`. This is enforced by the controller. Relocating an old checkpoint can change path-bearing identities: create new stage outputs from correctly described inputs; do not rewrite a receipt to claim prior validation.
