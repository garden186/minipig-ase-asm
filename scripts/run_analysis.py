"""Discover and run a final-analysis entry point in an isolated Python process."""
import argparse,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
COMMANDS={
 'ase-count':'analysis/ase_fragment_counts/call_phaser_ase_candidates.py',
 'ase-classify':'analysis/ase_reclassification/reclassify_ase_v026.py',
 'ase-cohort':'analysis/ase_cohort/build_ase_v026_cohort.py',
 'ase-recurrence':'analysis/ase_recurrence/core_eligible_recurrence.py',
 'ase-tissues':'analysis/ase_tissue_enrichment/ase_tissue_specificity.py',
 'ase-select-tissues':'scripts/select_primary_tissue_pairs.py',
 'ase-sites':'analysis/ase_snp_level/call_snp_level_ase.py',
 'ase-site-recurrence':'analysis/ase_snp_recurrence/call_snp_recurrence.py',
 'ase-site-loci':'analysis/ase_snp_loci/collapse_recurrent_snp_loci.py',
 'ase-site-annotation':'analysis/ase_snp_annotation/annotate_snp_features.py',
 'gatk-ase':'analysis/ase_gatk_calls/call_gatk_snp_level_ase.py',
 'gatk-count':'analysis/ase_gatk_calls/count_reads.py',
 'gatk-recurrence':'analysis/ase_gatk_recurrence/call_gatk_snp_recurrence.py',
 'gatk-blocks':'analysis/ase_gatk_blocks/validate_gatk_against_primary_blocks.py',
 'piggtex':'analysis/piggtex/piggtex_validation.py',
 'piggtex-sites':'analysis/piggtex_sites/integrate_gatk_with_piggtex_sites.py',
 'piggtex-eqtl-universe':'analysis/piggtex_eqtl/prepare_full_universe.py',
 'piggtex-eqtl':'analysis/piggtex_eqtl/audit_full_gene_eqtl_finemap.py',
 'piggtex-pair-support':'analysis/piggtex_pair_support/summarize_pair_support.py',
 'piggtex-multilayer':'analysis/piggtex_multilayer/build_multilayer.py',
 'piggtex-evidence':'analysis/piggtex_evidence/build_evidence_audit.py',
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
