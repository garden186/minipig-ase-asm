"""Discover and run a final-analysis entry point in an isolated Python process."""
import argparse,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
COMMANDS={
 'ase-count':'analysis/ase_fragment_counts/call_phaser_ase_candidates.py',
 'ase-classify':'analysis/ase_reclassification/reclassify_ase.py',
 'ase-cohort':'analysis/ase_cohort/build_ase_cohort.py',
 'ase-recurrence':'analysis/ase_recurrence/core_eligible_recurrence.py',
 'ase-tissues':'analysis/ase_tissue_enrichment/ase_tissue_specificity.py',
 'ase-select-tissues':'scripts/select_primary_tissue_pairs.py',
 'piggtex':'analysis/piggtex/piggtex_validation.py',
 'asm-count':'analysis/asm_counting/count.py',
 'asm-gq20':'analysis/asm_counting/filter_gq20.py',
 'asm-individual':'analysis/asm_individual/infer.py',
 'asm-cache':'analysis/asm_cache/run.py',
 'asm-anchors':'analysis/asm_allele_prep/run.py',
 'asm-recurrence':'analysis/asm_directional/run.py',
 'asm-tissues':'analysis/asm_tissue_comparison/run.py',
 'annotation':'analysis/annotation/build_index.py',
 'annotate-cpgs':'analysis/annotation/annotate_table.py',
 'asm-composition':'analysis/annotation/individual_composition.py',
 'external-overlap':'analysis/external_overlap/point_overlap.py',
 'integration':'analysis/integration/link.py',
}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('analysis',choices=['list']+list(COMMANDS));p.add_argument('arguments',nargs=argparse.REMAINDER);a=p.parse_args()
    if a.analysis=='list':
        for key,path in COMMANDS.items():print(key+'\t'+path)
        return 0
    argv=a.arguments
    if argv[:1]==['--']:argv=argv[1:]
    # Adopted final tissue-enrichment settings, overriding historical defaults.
    if a.analysis=='ase-tissues' and not any(x in argv for x in ['-h','--help']):
        for option,value in [('--thresholds','15'),('--min-informative-animals','5'),('--min-informative-target-core','2'),('--fdr','0.05')]:
            if option not in argv and not any(x.startswith(option+'=') for x in argv):argv += [option,value]
    return subprocess.call([sys.executable,'-X','utf8',str(ROOT/COMMANDS[a.analysis]),*argv])
if __name__=='__main__':raise SystemExit(main())
