"""Discover and run a final-analysis entry point in an isolated Python process."""

import argparse, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMMANDS = {
    "ase-validate-count": "analysis/ase_fragment_counts/validate_fragment_outputs.py",
    "ase-validate-classification": "analysis/ase_reclassification/validate_ase_outputs.py",
    "ase-validate-cohort": "analysis/ase_cohort/validate_ase_cohort.py",
    "ase-validate-recurrence": "analysis/ase_recurrence/validate_outputs.py",
    "ase-validate-tissues": "analysis/ase_tissue_enrichment/validate_outputs.py",
    "piggtex-validate": "analysis/piggtex/validate_outputs.py",
    "preprocessing-qc": "preprocessing/helpers/stage2_qc_collect.py",
    "go": "analysis/go/enrich_final.R",
    "ase-count": "analysis/ase_fragment_counts/call_phaser_ase_candidates.py",
    "ase-classify": "analysis/ase_reclassification/reclassify_ase.py",
    "ase-cohort": "analysis/ase_cohort/build_ase_cohort.py",
    "ase-recurrence": "analysis/ase_recurrence/core_eligible_recurrence.py",
    "ase-tissues": "analysis/ase_tissue_enrichment/ase_tissue_specificity.py",
    "ase-select-tissues": "scripts/select_primary_tissue_pairs.py",
    "piggtex": "analysis/piggtex/piggtex_validation.py",
    "asm-count": "analysis/asm_counting/count.py",
    "asm-reassign": "analysis/asm_counting/engine/reassign.py",
    "asm-gq20": "analysis/asm_counting/filter_gq20.py",
    "asm-individual": "analysis/asm_individual/infer.py",
    "asm-cache": "analysis/asm_cache/run.py",
    "asm-anchors": "analysis/asm_allele_prep/run.py",
    "asm-recurrence": "analysis/asm_directional/run.py",
    "asm-tissues": "analysis/asm_tissue_comparison/run.py",
    "annotation": "analysis/annotation/build_index.py",
    "annotate-cpgs": "analysis/annotation/annotate_table.py",
    "asm-composition": "analysis/annotation/individual_composition.py",
    "external-overlap": "analysis/external_overlap/point_overlap.py",
    "integration": "analysis/integration/link.py",
}


def main():
    if sys.version_info < (3, 12):
        raise SystemExit("Python >= 3.12 is required")
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("analysis", choices=["list"] + list(COMMANDS))
    p.add_argument("arguments", nargs=argparse.REMAINDER)
    a = p.parse_args()
    if a.analysis == "list":
        for key, path in COMMANDS.items():
            print(key + "\t" + path)
        return 0
    argv = a.arguments
    if argv[:1] == ["--"]:
        argv = argv[1:]
    # Tissue-enrichment thresholds.
    if a.analysis == "ase-tissues" and not any(x in argv for x in ["-h", "--help"]):
        for option, value in [
            ("--thresholds", "15"),
            ("--min-informative-animals", "5"),
            ("--min-informative-target-core", "2"),
            ("--fdr", "0.05"),
        ]:
            if option not in argv and not any(x.startswith(option + "=") for x in argv):
                argv += [option, value]
    target = ROOT / COMMANDS[a.analysis]
    if target.suffix == ".R":
        if any(x in argv for x in ["-h", "--help"]):
            print(
                "Usage: run_analysis.py go specification.json (requires Rscript, jsonlite, clusterProfiler)"
            )
            return 0
        import shutil

        executable = shutil.which("Rscript")
        if executable is None:
            raise SystemExit("Rscript is required for GO analysis")
        return subprocess.call([executable, str(target), *argv])
    return subprocess.call([sys.executable, "-X", "utf8", str(target), *argv])


if __name__ == "__main__":
    raise SystemExit(main())
