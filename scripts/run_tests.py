"""Run native test suites as separate processes; preserve exact outcomes."""
from pathlib import Path
import argparse,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--asm-config',type=Path);a=p.parse_args()
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    suites=[('portability',ROOT/'tests','test_*.py')]
    for name,pattern in [('ase_reclassification','test_reclassify_ase_v026.py'),('ase_recurrence','test_core_eligible_recurrence.py'),('ase_tissue_enrichment','test_ase_tissue_specificity.py'),('piggtex','test_piggtex_validation.py'),('ase_snp_recurrence','test_snp_recurrence.py')]:suites.append((name,ROOT/'analysis'/name,pattern))
    suites.extend([('ase_fragment_logic',ROOT/'analysis/ase_fragment_counts','test_v025_logic.py'),('ase_fragment_candidates',ROOT/'analysis/ase_fragment_counts','test_call_phaser_ase_candidates.py'),('ase_cohort',ROOT/'analysis/ase_cohort','test_ase_v026_cohort.py')])
    suites.append(('asm_fragment_evidence',ROOT/'analysis/asm_counting/engine/tests','test_*.py'))
    for name in ['ase_snp_level','ase_snp_loci','ase_snp_annotation','ase_gatk_calls','ase_gatk_recurrence','ase_gatk_blocks','piggtex_sites','piggtex_eqtl','piggtex_pair_support','piggtex_evidence']:
        suites.append((name,ROOT/'analysis'/name,'test_*.py'))
    jobs=[(name,[sys.executable,'-X','utf8','-m','unittest','discover','-s',str(folder),'-p',pattern],os.environ.copy()) for name,folder,pattern in suites]
    if a.asm_config:
        config=json.loads(a.asm_config.read_text(encoding='utf-8'))
        for name in ['asm_allele_prep','asm_directional','asm_tissue_comparison']:
            env=os.environ.copy();env['MINIPIG_ASSETS']=str(Path(config[name]['assets']).resolve());env['MINIPIG_ANALYSIS_CONFIG']=str(Path(config[name]['inputs']).resolve())
            jobs.append((name,[sys.executable,'-X','utf8',str(ROOT/'analysis'/name/'selftest.py'),'--record',str(a.output/(name+'.json'))],env))
    results=[]
    for name,command,env in jobs:
        start=time.monotonic()
        with (a.output/(name+'.log')).open('w',encoding='utf-8') as f:r=subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        results.append(dict(suite=name,exit_code=r.returncode,seconds=time.monotonic()-start));print(name+': '+('PASS' if r.returncode==0 else 'FAILED'),flush=True)
    result=dict(status='PASS' if all(r['exit_code']==0 for r in results) else 'FAILED',ASM_included=bool(a.asm_config),suites=results)
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return 0 if result['status']=='PASS' else 1

if __name__=='__main__':raise SystemExit(main())
