"""Exact-support calibration, frozen real counts, and real multiprocessing/resume."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import argparse,collections,copy,hashlib,itertools,json,math,sys,tempfile,time,unittest
from pathlib import Path
import numpy as np
from scipy.stats import chi2,false_discovery_control
from common import *
from schema import PROFILE,SITE,STATE,GT
from data import *
from exact import pair_test,reference,sign,paired_difference
from fractions import Fraction
from pipeline import analyze,family_job
from stats import partial_conjunction,bh_log
from run import execute

EVIDENCE={}
WORK=None

def real_arrays():
    rr=list(rows(asset_path('assets/real_profiles.tsv.gz')));bychrom=collections.defaultdict(list)
    for r in rr:bychrom[int(r['chrom'])].append(r)
    result={}
    for chrom,records in bychrom.items():
        positions=sorted({int(r['cpg_pos1']) for r in records});lookup={p:i for i,p in enumerate(positions)}
        sites=np.zeros(len(positions),dtype=SITE);sites['pos']=positions
        pr=np.zeros(len(records),dtype=PROFILE);ps={s:{v:i+1 for i,v in enumerate(sorted({r['anchor_PS'] for r in records if r['sample']==s and r['anchor_PS']!='.'}))} for s in IDS}
        original=set()
        for i,r in enumerate(records):
            pos=int(r['cpg_pos1']);j=lookup[pos];sample=IDS.index(r['sample']);tissue=TISSUES.index(r['tissue']);z=pr[i]
            z['pos']=pos;z['sample']=sample;z['tissue']=tissue;z['source_row']=int(r['source_row']);z['ps']=i+1
            for dest,src in zip('abcd',['H1_methylated','H1_unmethylated','H2_methylated','H2_unmethylated']):z[dest]=int(r[src])
            z['original_p']=float(r['original_fisher_two_sided_p']);z['original_q']=float(r['original_global30_q'])
            z['state']=STATE.index(r['orientation_status']);z['GT']=GT.index(r['anchor_GT'])
            z['anchor_ps']=ps[r['sample']].get(r['anchor_PS'],0);z['alt_hap']=2 if r['anchor_GT']=='0|1' else 1 if r['anchor_GT']=='1|0' else 0
            if r['anchor_pos1']!='.':sites['anchor'][j]=int(r['anchor_pos1'])*16+'ACGT'.index(r['anchor_ref'])*4+'ACGT'.index(r['anchor_alt'])
            original.add((j,sample,tissue))
            if z['state']==0:sites['N'][j,tissue]+=1
        for j,sample,tissue in original:sites['N_before'][j,tissue]+=1
        result[chrom]=(sites,pr,records)
    return result

def fixture_config(directory):
    config=copy.deepcopy(read(input_config()));config['fixture']=True;config['root_checks']={};config['tiles']=[]
    counts=dict(tiles=0,sites=0,profiles=0,aligned=0)
    for chrom,(sites,profiles,_) in sorted(real_arrays().items()):
        if chrom==1:
            # Explicit synthetic alternative exercises BH10 discovery export at
            # every r. It is never part of the frozen biological regression.
            selfpos=42
            assert selfpos not in sites['pos']
            added=np.zeros(1,dtype=SITE);added['pos']=selfpos;added['anchor']=41*16+3
            added['N']=10;added['N_before']=10;sites=np.concatenate([added,sites]);sites.sort(order='pos')
            extra=np.zeros(30,dtype=PROFILE)
            for animal in range(10):
                for tissue,synthetic_counts in enumerate([(60,0,0,60),(45,15,15,45),(0,60,60,0)]):
                    z=extra[animal*3+tissue];z['pos']=selfpos;z['sample']=animal;z['tissue']=tissue;z['source_row']=100000000+animal*3+tissue
                    z['state']=0;z['GT']=2;z['alt_hap']=1;z['anchor_ps']=999999;z['ps']=tissue+1;z['original_p']=1.;z['original_q']=1.
                    for k,v in zip('abcd',synthetic_counts):z[k]=v
            profiles=np.concatenate([extra,profiles])
        p=directory/f'chr{chrom}';p.mkdir(parents=True)
        np.save(p/'sites.npy',sites,allow_pickle=False);np.save(p/'profiles.npy',profiles,allow_pickle=False)
        spec=dict(chrom=chrom,start=int(sites['pos'][0]),end=int(sites['pos'][-1])+1,path=p.as_posix(),sites=len(sites),profiles=len(profiles),aligned=int(np.count_nonzero(profiles['state']==0)))
        spec['checks']={(p/n).as_posix():sha(p/n) for n in ['sites.npy','profiles.npy']};config['tiles'].append(spec)
        counts['tiles']+=1;counts['sites']+=len(sites);counts['profiles']+=len(profiles);counts['aligned']+=spec['aligned']
    config['expected']=counts;return config

class Tests(unittest.TestCase):
    def test_01_integer_oracle_and_conditional_size(self):
        checked=groups_checked=0;worst=0.
        for na,nr,nb,ns in [(3,3,3,3),(3,7,5,5),(4,6,8,3),(3,15,3,17)]:
            groups=collections.defaultdict(list)
            for a,c,e,g in itertools.product(range(na+1),range(nr+1),range(nb+1),range(ns+1)):
                one=(a,na-a,c,nr-c);two=(e,nb-e,g,ns-g)
                actual=pair_test(one,two);ref=reference(one,two)
                self.assertAlmostEqual(actual[0],math.log(ref['numerator'])-math.log(ref['denominator']),places=11)
                self.assertEqual(actual[1],ref['full_states'])
                weight=math.comb(na,a)*math.comb(nr,c)*math.comb(nb,e)*math.comb(ns,g)
                groups[(a+c,e+g,a+e)].append((math.exp(actual[0]),weight))
                checked+=1
            for values in groups.values():
                denominator=sum(w for _,w in values);cdf=0
                for p in sorted({p for p,_ in values}):
                    cdf=sum(w for q,w in values if q<=p)/denominator
                    worst=max(worst,cdf-p);self.assertLessEqual(cdf,p+1e-12)
                groups_checked+=1
        EVIDENCE['exact_enumeration']=dict(observed_table_pairs=checked,conditional_null_support_groups=groups_checked,max_CDF_minus_P=worst)

    def test_02_symmetry_extreme_and_underflow(self):
        for first,second in [((30,1,1,30),(1,30,30,1)),((10,5,3,9),(4,8,6,7)),((3,3,3,3),(9,0,0,9))]:
            a=pair_test(first,second)
            b=pair_test(first[2:]+first[:2],second[2:]+second[:2]);c=pair_test(second,first)
            self.assertAlmostEqual(a[0],b[0],places=11);self.assertAlmostEqual(a[0],c[0],places=11)
            self.assertEqual(a[1],b[1]);self.assertEqual(a[1],c[1])
        huge=pair_test((2000,0,0,2000),(0,2000,2000,0))
        self.assertTrue(math.isfinite(huge[0]) and huge[0]<-745)
        with self.assertRaises(ValueError):pair_test((3,2,1,-1),(3,2,1,3))
        with self.assertRaises(ValueError):pair_test((3,3,3,3),(3,3,3,3),1)
        EVIDENCE['underflow_example_logP']=huge[0]

    def test_03_frozen_profiles(self):
        checked=0;pairs=0;max_error=0.;pc_checks=0;started=time.monotonic()
        for chrom,(sites,pr,rr) in real_arrays().items():
            counts,phase,links=orient(pr,sites);delta=deltas(counts)
            for r in rr:
                if r['orientation_status']!='ALIGNED':continue
                j=np.searchsorted(sites['pos'],int(r['cpg_pos1']));a=IDS.index(r['sample']);t=TISSUES.index(r['tissue'])
                expected=[int(r[k]) for k in ['ALT_methylated','ALT_unmethylated','REF_methylated','REF_unmethylated']]
                self.assertEqual(counts[j,a,t].tolist(),expected);self.assertAlmostEqual(delta[j,a,t],float(r['delta_ALT_REF']),places=13);checked+=1
            f,z,individual=analyze(counts,sites,np.zeros(len(sites),dtype=np.uint32))
            self.assertEqual(len(z),int(np.count_nonzero(f['paired_N']>=2)))
            for row,result in enumerate(z):
                j=np.searchsorted(sites['pos'],result['pos']);u,v=PAIRS[int(result['pair'])]
                for animal in np.flatnonzero(np.isfinite(individual[row]['logp'])):
                    ref=reference(counts[j,animal,u],counts[j,animal,v]);lp=math.log(ref['numerator'])-math.log(ref['denominator'])
                    max_error=max(max_error,abs(lp-individual[row]['logp'][animal]));pairs+=1
                valid=np.isfinite(individual[row]['logp']);values=np.sort(individual[row]['logp'][valid]);N=len(values)
                self.assertEqual(N,int(result['N']))
                for r in range(2,N+1):
                    expected=chi2.logsf(-2*values[r-1:].sum(),2*(N-r+1))
                    self.assertAlmostEqual(float(result['log_pc'][r-2]),float(expected),places=10);pc_checks+=1
                measured=individual[row]['D'][valid]
                self.assertAlmostEqual(float(result['mean_D']),float(np.mean(measured)),places=13)
                self.assertAlmostEqual(float(result['median_D']),float(np.median(measured)),places=13)
        self.assertGreater(checked,1000);self.assertGreater(pairs,10);self.assertLess(max_error,1e-10)
        EVIDENCE['frozen_real']=dict(aligned_profiles=checked,paired_animal_tests=pairs,PC_oracles=pc_checks,OR_logP_max_error=max_error,seconds_including_independent_oracle=time.monotonic()-started)

    def test_04_keep_opposing_ties_and_ignore_third_tissue(self):
        x=np.full((1,10,3,4),-1,dtype=np.int32)
        x[0,:2,0]=[8,2,2,8];x[0,:2,2]=[2,8,8,2]
        # Third tissue mixes; Heart-Liver still passes.
        x[0,0,1]=[8,2,2,8];x[0,1,1]=[2,8,8,2]
        f,_,_=filters(x,[100]);self.assertEqual(int(f['reason'][0,1]),0)
        # Opposing and tied Heart-only animals do not remove paired animals.
        x[0,2,0]=[2,8,8,2];f,_,_=filters(x,[100]);self.assertEqual(int(f['reason'][0,1]),0)
        x[0,2,0]=[5,5,5,5];f,_,_=filters(x,[100]);self.assertEqual(int(f['reason'][0,1]),0)
        self.assertEqual(int(f['paired_N'][0,1]),2)
        x[0,1,0]=[5,5,5,5];x[0,1,2]=[5,5,5,5]
        sites=np.zeros(1,dtype=SITE);sites['pos']=100;sites['anchor']=101*16+1
        _,z,individual=analyze(x,sites,np.zeros(1,dtype=np.uint32))
        row=int(np.flatnonzero(z['pair']==1)[0])
        self.assertEqual(int(z['N'][row]),2);self.assertEqual(int(z['D_sign_counts'][row,2]),1)
        self.assertEqual(float(individual['D'][row,1]),0.)
        self.assertTrue(np.isfinite(individual['logp'][row,1]))
        self.assertFalse(np.isfinite(individual['logp'][row,2]))
        x[0,2,0]=-1;x[0,1,2]=-1;f,_,_=filters(x,[100]);self.assertEqual(int(f['reason'][0,1]),1)

    def test_05_phase_duplicate_and_depth_refusal(self):
        sites,pr,_=next(iter(real_arrays().values()));align=np.flatnonzero(pr['state']==0)[0]
        with self.assertRaises(ValueError):orient(np.concatenate([pr,pr[[align]]]),sites)
        bad=pr.copy();bad['alt_hap'][align]=3
        with self.assertRaises(ValueError):orient(bad,sites)
        bad=pr.copy();bad['a'][align]=0;bad['b'][align]=0
        with self.assertRaises(ValueError):orient(bad,sites)
        bad=pr.copy();bad['pos'][align]=sites['pos'][-1]+1
        with self.assertRaises(ValueError):orient(bad,sites)

    def test_06_PC_and_complete_family_correction(self):
        rng=np.random.default_rng(1729)
        for N in range(2,11):
            p=rng.uniform(.00001,1,(64,N));lp=np.log(p)
            for r in range(2,N+1):
                expected=chi2.logsf(-2*np.log(np.sort(p,axis=1)[:,r-1:]).sum(axis=1),2*(N-r+1))
                np.testing.assert_allclose(partial_conjunction(lp,r),expected,atol=1e-12,rtol=1e-12)
        p=np.array([.0001,.001,.03,.08,.5,1.]);q=np.exp(bh_log(np.log(p)))
        np.testing.assert_allclose(q,false_discovery_control(p,method='bh'),atol=1e-15)
        for fullM in [6,37,2000]:
            harmonic=sum(1/i for i in range(1,fullM+1))
            got=np.minimum(1.,q*(fullM/len(p))*harmonic)
            oracle=false_discovery_control(np.r_[p,np.ones(fullM-len(p))],method='by')[:len(p)]
            np.testing.assert_allclose(got,oracle,atol=1e-14)
        self.assertTrue(np.all(np.isfinite(partial_conjunction(np.array([[-3000.,-2000.,-1000.]]),2))))

    def test_07_end_to_end_parallel_resume_tamper(self):
        source=WORK/'input';config=fixture_config(source);out=WORK/'output';rh=sha(ROOT/'release.json')
        first=execute(config,out,rh,workers=2,fixture=True)
        locations=read(out/'result_locations.json');final=Path(locations['final'])
        hashes={p.name:sha(p) for p in final.iterdir() if p.is_file()}
        second=execute(config,out,rh,workers=2,fixture=True)
        self.assertEqual(first,second);self.assertEqual(hashes,{p.name:sha(p) for p in final.iterdir() if p.is_file()})
        self.assertEqual(len(locations['families']),27)
        self.assertGreater(first['BH05_hypotheses'],0)
        self.assertTrue(any(r['cpg_pos1']=='42' for r in rows(final/'tissue_contrasts_BH10.tsv.gz')))
        # Direct all-family N and SciPy complete-family BH/BY replay.
        direct_rows=0
        for path in locations['families']:
            a=np.load(Path(path)/'family.npy',allow_pickle=False);s=read(Path(path)/'summary.json');direct_rows+=len(a)
            if not len(a):continue
            p=np.exp(a['logp']);np.testing.assert_allclose(np.exp(a['logq_BH']),false_discovery_control(p,method='bh'),atol=3e-14)
            self.assertEqual(s['M'],len(p))
            expected=sum(int(np.count_nonzero(np.load(Path(t)/'availability.npy',allow_pickle=False)['paired_N'][:,s['pair_index']]>=s['r'])) for t in locations['tiles'])
            self.assertEqual(len(a),expected)
            np.testing.assert_allclose(np.exp(a['logq_BY']),false_discovery_control(p,method='by'),atol=3e-14)
        profiles=list(rows(final/'all_eligible_profiles.tsv.gz'));paired=list(rows(final/'paired_individual_results.tsv.gz'))
        lookup={(r['chrom'],r['cpg_pos1'],r['sample'],r['tissue']):r for r in profiles if r['orientation_status']=='ALIGNED'}
        for r in paired:
            if r['paired']!='1':continue
            u,v=PAIR_NAMES.index(r['pair']),None;one,two=[TISSUES[t] for t in PAIRS[u]]
            one=lookup[(r['chrom'],r['cpg_pos1'],r['sample'],one)];two=lookup[(r['chrom'],r['cpg_pos1'],r['sample'],two)]
            ref=reference(tuple(int(one[k]) for k in ['ALT_M','ALT_U','REF_M','REF_U']),tuple(int(two[k]) for k in ['ALT_M','ALT_U','REF_M','REF_U']))
            self.assertAlmostEqual(float(r['log_OR_P']),math.log(ref['numerator'])-math.log(ref['denominator']),places=11)
            d,sgn=paired_difference(tuple(int(one[k]) for k in ['ALT_M','ALT_U','REF_M','REF_U']),tuple(int(two[k]) for k in ['ALT_M','ALT_U','REF_M','REF_U']))
            self.assertEqual(int(r['D_sign']),sgn);self.assertAlmostEqual(float(r['D_pp']),100*d,places=12)
        plot=list(rows(final/'pairwise_plot_source_BH05.tsv'))
        hitrows=list(rows(final/'tissue_contrasts_BH10.tsv.gz'))
        for plotted in plot:
            key=(plotted['pair'],plotted['chrom'],plotted['cpg_pos1'])
            tested=[r for r in hitrows if (r['pair'],r['chrom'],r['cpg_pos1'])==key and r['BH05']=='1']
            self.assertEqual(int(plotted['max_significant_r']),max(int(r['r']) for r in tested))
            observed=[r for r in paired if (r['pair'],r['chrom'],r['cpg_pos1'])==key and r['paired']=='1']
            self.assertEqual(int(plotted['N_paired']),len(observed))
            self.assertAlmostEqual(float(plotted['median_D_pp']),np.median([float(r['D_pp']) for r in observed]),places=11)
        # Tampered source is refused even when checkpoints already exist.
        p=Path(config['tiles'][0]['path'])/'sites.npy';original=p.read_bytes();p.write_bytes(original+b'x')
        try:
            with self.assertRaises(ValueError):execute(config,out,rh,workers=2,fixture=True)
            self.assertFalse((out/'validation.json').exists())
        finally:p.write_bytes(original)
        restored=execute(config,out,rh,workers=2,fixture=True)
        self.assertEqual(first,restored)
        self.assertEqual(hashes,{p.name:sha(p) for p in final.iterdir() if p.is_file()})
        EVIDENCE['end_to_end']=dict(input=config['expected'],family_rows=direct_rows,returned_profiles=len(profiles),paired_rows=len(paired),resume_identical=True,tamper_refused=True,summary=first)

    def test_08_resource_lock_and_stop(self):
        p=WORK/'locks';p.mkdir()
        with lock(p):
            with self.assertRaises(AlreadyRunning):
                with lock(p):pass
        write(p/'STOP_REQUESTED',dict(reason='test'))
        with self.assertRaises(StopRequested):parallel(abs,[],p,'test',workers=1,fixture=True)
        r=resources();self.assertGreaterEqual(r['cpus'],1);self.assertGreaterEqual(r['reserve'],2*GIB)
        self.assertEqual(os.environ['OPENBLAS_NUM_THREADS'],'1')

    def test_09_empty_family_and_checkpoint_tamper(self):
        p=WORK/'empty_source';p.mkdir();f=np.zeros(1,dtype=FILTER);f['pos']=1;f['paired_N']=1;f['reason']=1
        np.save(p/'availability.npy',f,allow_pickle=False);np.save(p/'results.npy',np.zeros(0,dtype=RESULT),allow_pickle=False)
        write(p/'summary.json',dict(chrom=1))
        job=dict(parent=str(WORK/'empty_family'),identity='fixture_empty',pair=0,r=2,tiles=[str(p)],fixture=True)
        out=Path(family_job(job));s=read(out/'summary.json')
        self.assertEqual(s['M'],0);self.assertEqual(s['BH05'],0)
        with (out/'family.npy').open('ab') as f:f.write(b'x')
        with self.assertRaises(ValueError):family_job(job)

    def test_10_haplotype_orientation_invariance_and_unknown(self):
        sites,pr,_=next(iter(real_arrays().values()));base,_,_=orient(pr,sites)
        changed=pr.copy();good=changed['state']==0
        original=np.column_stack([pr[k][good] for k in 'abcd'])
        for k,value in zip('abcd',original[:,[2,3,0,1]].T):changed[k][good]=value
        changed['GT'][good]=3-changed['GT'][good];changed['alt_hap'][good]=3-changed['alt_hap'][good]
        counts,_,_=orient(changed,sites);np.testing.assert_array_equal(counts,base)
        aligned=np.flatnonzero(good);one=int(aligned[0]);p=pr['pos'][one];a=pr['sample'][one]
        same=np.flatnonzero(good&(pr['pos']==p)&(pr['sample']==a))
        if len(same)>=2:
            bad=pr.copy();bad['anchor_ps'][same[0]]+=100000
            with self.assertRaises(ValueError):orient(bad,sites)
        j=int(np.searchsorted(sites['pos'],p));t=int(pr['tissue'][one])
        changed=pr.copy();changed['state'][one]=2;ss=sites.copy();ss['N'][j,t]-=1
        missing,_,_=orient(changed,ss)
        self.assertTrue(np.all(missing[j,int(a),t]==-1))
        self.assertTrue(np.isnan(deltas(missing)[j,int(a),t]))

    def test_11_exact_fraction_sign_and_constant_not_excluded(self):
        rng=np.random.default_rng(92125)
        for _ in range(300):
            one=tuple(map(int,rng.integers(1,100,4)));two=tuple(map(int,rng.integers(1,100,4)))
            a,b,c,d=one;e,f,g,h=two
            exact=Fraction(a,a+b)-Fraction(c,c+d)-Fraction(e,e+f)+Fraction(g,g+h)
            value,sgn=paired_difference(one,two)
            self.assertEqual(value,float(exact));self.assertEqual(sgn,(exact>0)-(exact<0))
        # Equal nonzero differences across animals remain eligible count tables.
        x=np.full((1,10,3,4),-1,dtype=np.int32);x[0,:3,0]=[8,2,2,8];x[0,:3,2]=[5,5,5,5]
        sites=np.zeros(1,dtype=SITE);sites['pos']=100;sites['anchor']=101*16+1
        _,z,individual=analyze(x,sites,np.zeros(1,dtype=np.uint32))
        self.assertEqual(len(z),1);self.assertEqual(int(z['N'][0]),3)
        self.assertTrue(np.all(np.isfinite(z['log_pc'][0,:2])))
        self.assertEqual(z['min_D'][0],z['max_D'][0])

    def test_12_missing_eligible_family_is_refused(self):
        p=WORK/'lost_source';p.mkdir();f=np.zeros(1,dtype=FILTER);f['pos']=100;f['paired_N']=4
        np.save(p/'availability.npy',f,allow_pickle=False);np.save(p/'results.npy',np.zeros(0,dtype=RESULT),allow_pickle=False)
        write(p/'summary.json',dict(chrom=1))
        with self.assertRaises(ValueError):family_job(dict(parent=str(WORK/'lost_family'),identity='lost',pair=0,r=2,tiles=[str(p)],fixture=True))

def main():
    global WORK
    p=argparse.ArgumentParser();p.add_argument('--record',required=True);p.add_argument('--work');a=p.parse_args()
    release_hash=release();WORK=Path(a.work).resolve() if a.work else Path(tempfile.mkdtemp(prefix='asm_aligned_selftest_')).resolve();WORK.mkdir(parents=True,exist_ok=True)
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
    evidence=dict(status='PASS' if result.wasSuccessful() else 'FAIL',release_sha256=release_hash,tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),evidence=EVIDENCE,work=str(WORK),python=sys.version)
    write(a.record,evidence);sys.exit(0 if result.wasSuccessful() else 1)
if __name__=='__main__':main()
