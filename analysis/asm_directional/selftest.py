"""Exact integer oracles, direction/PC/BH semantics, optimized anchors and resume."""
import os
for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[n]='1'
import argparse,collections,csv,gzip,itertools,math,platform,sys,tempfile,time,unittest
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch
import numpy as np
from scipy.stats import fisher_exact,chi2,hypergeom
import common,anchors,stats,pipeline,reporting,run
from common import *

WORK=None
VAR=np.dtype([('pos','<u4'),('ref','u1'),('alt','u1'),('gt','u1'),('gq','<f8'),('ps','<u4'),('dup','u1'),('status','u1'),('alt_hap','u1'),('gq_source','u1'),('conversion','u1')])
CACHE=np.dtype([('chrom','u1'),('pos','<u4'),('ps','<u4'),('a','<u4'),('b','<u4'),('c','<u4'),('d','<u4'),('p','<f8'),('q_unit','<f8'),('q_tissue','<f8'),('q_all','<f8'),('source_row','<u4')])
AVAIL=np.dtype([('pos','<u4'),('ps','<u4')])

def tmp(name):
    p=WORK/(name+'_'+str(time.time_ns()));p.mkdir();return p

def gate(path,summary=None):
    if summary is not None:write(path/'summary.json',summary)
    files={p.name:sha(p) for p in path.iterdir() if p.is_file() and p.name!='validation.json'}
    write(path/'validation.json',dict(status='PASS_FIXTURE',errors=[],identity='fixture',sha256=files))
    return dict(gate_sha256=sha(path/'validation.json'),sha256=files)

def fixture(root):
    stage1=root/'stage1';stage1.mkdir();registry={};loc=dict(units={},phase_index={},pilot={},variant_sources={});specs=[]
    expected_profiles=[];allprofiles=[]
    for i,s in enumerate(IDS):
        phased=stage1/('phased_'+s);phased.mkdir();write(phased/'phase_dictionary.json',['.','10','20']);registry[str(phased)]=gate(phased,{})
        loc['variant_sources'][s+'-phased']=str(phased)
        pp=stage1/('phase_'+s);pp.mkdir()
        variants=[(90,0,2,1 if i==0 else 2,30.,1,0,0,2 if i==0 else 1,1,1),(290,0,2,1 if i==0 else 2,30.,2,0,0,2 if i==0 else 1,1,1)] if i<2 else []
        if i==0:variants.insert(1,(190,0,3,1,30.,2,0,0,2,1,0))
        a=np.array(variants,dtype=VAR);np.save(pp/'variants.npy',a);np.save(pp/'eligible_by_ps.npy',a[np.argsort(a['ps'],kind='stable')]);registry[str(pp)]=gate(pp,dict(phase_dictionary=str(phased/'phase_dictionary.json')));loc['phase_index'][s+'-chr1']=str(pp)
        for t in TISSUES:
            unit=s+'-'+t;cp=root/('cache_'+unit);cp.mkdir();up=stage1/('available_'+unit);up.mkdir()
            # Different cache and raw PS dictionaries deliberately exercise mapping.
            cps=['30','20','10'];rawps=['10','20','30'];write(cp/'phase_dictionary.json',cps);write(up/'phase_dictionary.json',rawps)
            entries=[]
            if i<2:entries.append((100,'10',*( (9,1,1,9) if i==0 else (1,9,9,1))))
            if i==0 and t=='Heart':entries.extend([(100,'20',7,3,1,9),(200,'20',6,4,1,9)])
            if i==1 and t=='Heart':entries.append((300,'30',9,1,1,9))
            entries.sort(key=lambda v:(v[0],cps.index(v[1])));ar=[];source=[]
            for j,e in enumerate(entries):
                pos,ps,A,B,C,D=e;p=float(fisher_exact([[A,B],[C,D]]).pvalue);ar.append((1,pos,cps.index(ps),A,B,C,D,p,1.,1.,1.,j+1))
                src=dict(sample=s,tissue=t,chrom=1,cpg_pos1=pos,ps=ps,source_row=j+1,H1_methylated=A,H1_unmethylated=B,H2_methylated=C,H2_unmethylated=D,annotation_gene_name='SYNTHETIC_FIXTURE',annotation_primary_feature='fixture')
                source.append(src)
                state='NO_ELIGIBLE_SITE_ANCHOR' if pos==300 else 'DIFFERENT_EXACT_PS' if pos==100 and ps=='20' else 'ALIGNED'
                O=(A,B,C,D) if i==1 else (C,D,A,B)
                expected_profiles.append(dict(sample=s,tissue=t,chrom=1,cpg_pos1=pos,ps=ps,groups='fixture',H1_methylated=A,H1_unmethylated=B,H2_methylated=C,H2_unmethylated=D,
                    original_fisher_two_sided_p=p,original_q_all30=1.,orientation_status=state,ALT_methylation=O[0]/(O[0]+O[1]) if state=='ALIGNED' else '.',REF_methylation=O[2]/(O[2]+O[3]) if state=='ALIGNED' else '.'))
            np.save(cp/'chr1.npy',np.array(ar,dtype=CACHE));fields=list(source[0]) if source else ['sample','tissue','chrom','cpg_pos1','ps','source_row','annotation_gene_name','annotation_primary_feature']
            pipeline.table(cp/'chr1.source.tsv.gz',source,fields);cpgate=gate(cp,dict(rows=len(ar)))
            raw=[(e[0],rawps.index(e[1])) for e in entries]
            if i==1 and t=='Heart':raw.append((200,rawps.index('20')))
            av=np.array(sorted(raw),dtype=AVAIL);np.save(up/'chr1.available.npy',av);registry[str(up)]=gate(up,{})
            loc['units'][unit]=str(up)
            specs.append(dict(sample=s,tissue=t,cache=dict(path=str(cp),gate=str(cp/'validation.json'),gate_sha256=cpgate['gate_sha256'],sha256={n:h for n,h in cpgate['sha256'].items() if not n.endswith('.tsv.gz')})))
    pilot=stage1/'pilot';pilot.mkdir();np.save(pilot/'eligible_positions.npy',np.array([100,200,300],dtype='<u4'))
    es=[]
    for pos,ap,nt,na,cand in [(100,90,2,2,3),(200,290,1,2,2),(300,'.',0,0,0)]:
        row=dict(chrom=1,cpg_pos1=pos,anchor_pos1=ap,anchor_ref='A' if ap!='.' else '.',anchor_alt='G' if ap!='.' else '.',anchor_testable_animals=nt,anchor_available_animals=na,candidates=cand)
        for t in TISSUES:row['N_before_'+t]=2 if pos==100 else 1 if t=='Heart' else 0;row['N_aligned_'+t]=2 if pos==100 else 1 if pos==200 and t=='Heart' else 0
        es.append(row)
    pipeline.table(pilot/'pilot_sites.tsv.gz',es,list(es[0]));pipeline.table(pilot/'pilot_profiles.tsv.gz',expected_profiles,list(expected_profiles[0]));registry[str(pilot)]=gate(pilot,{})
    loc['pilot']['chr1']=str(pilot);write(stage1/'run_identity.json',dict(identity='fixture-stage1'));write(stage1/'summary.json',dict(fixture=True));rootgate=gate(stage1)
    return dict(stage1_root=str(stage1),stage1_identity='fixture-stage1',stage1_gate_sha256=rootgate['gate_sha256'],stage1_root_payloads=rootgate['sha256'],locations=loc,registry=registry,
        source_inputs=dict(units=specs,expected=dict(full_eligible_coordinates=3,tested_rows=9)),tile_bp=150,expected_pilot_sites=3,expected_pilot_profiles=9)

def brute(pos,available,tested,phase):
    candidates={}
    for i,(av,te,a) in enumerate(zip(available,tested,phase)):
        for v in a:
            if int(v['ps']) in av:
                k=int(v['pos'])*16+int(v['ref'])*4+int(v['alt']);x=candidates.setdefault(k,[set(),set()]);x[0].add(i)
                if int(v['ps']) in te:x[1].add(i)
    if not candidates:return 0,0,0,0
    key=min(candidates,key=lambda k:(-len(candidates[k][1]),-len(candidates[k][0]),abs(k//16-pos),k))
    return key,len(candidates[key][1]),len(candidates[key][0]),len(candidates)

class Tests(unittest.TestCase):
    def test_real_count_tables_against_integer_oracle(self):
        data=list(pipeline.rows(asset_path('real_count_oracle.tsv.gz')))
        self.assertEqual(len(data),1500)
        counts=np.array([[int(r[k]) for k in ['a','b','c','d']] for r in data],dtype=np.int64)
        expected=[];two=[]
        for a,b,c,d in counts.tolist():
            N=a+b+c+d;K=a+c;n=a+b;lo=max(0,n-(N-K));hi=min(K,n)
            weights=[math.comb(K,x)*math.comb(N-K,n-x) for x in range(lo,hi+1)]
            den=math.comb(N,n);observed=weights[a-lo]
            expected.append([sum(weights[a-lo:])/den,sum(weights[:a-lo+1])/den])
            two.append(sum(w for w in weights if w<=observed)/den)
        lp,fallback=stats.directional_logp(counts)
        np.testing.assert_allclose(np.exp(lp),expected,rtol=2e-12,atol=0)
        np.testing.assert_allclose(two,[float(r['original_two_sided_p']) for r in data],rtol=2e-8,atol=0)
        swapped,_=stats.directional_logp(counts[:,[2,3,0,1]])
        np.testing.assert_allclose(lp[:,::-1],swapped,rtol=2e-12,atol=2e-12)

    def test_exact_integer_tails_and_swap(self):
        tables=[];oracle=[]
        for a,b,c,d in itertools.product(range(6),repeat=4):
            if not (a+b and c+d):continue
            N=a+b+c+d;K=a+c;n=a+b;lo=max(0,n-(N-K));hi=min(K,n);den=math.comb(N,n)
            up=sum(math.comb(K,x)*math.comb(N-K,n-x) for x in range(a,hi+1))/den
            down=sum(math.comb(K,x)*math.comb(N-K,n-x) for x in range(lo,a+1))/den
            tables.append([a,b,c,d]);oracle.append([up,down])
        a=np.array(tables);lp,n=stats.directional_logp(a);np.testing.assert_allclose(np.exp(lp),oracle,rtol=2e-13,atol=1e-15)
        swapped,_=stats.directional_logp(a[:,[2,3,0,1]]);np.testing.assert_allclose(lp[:,::-1],swapped,rtol=2e-13,atol=1e-13)

    def test_component_null_calibration_all_margins(self):
        tables=[];groups=[];weights=[]
        for n1,n2 in [(3,3),(3,7),(7,11),(20,20)]:
            N=n1+n2
            for K in range(N+1):
                lo=max(0,n1-(N-K));hi=min(K,n1);first=len(tables)
                for x in range(lo,hi+1):tables.append([x,n1-x,K-x,n2-K+x]);weights.append(math.comb(K,x)*math.comb(N-K,n1-x)/math.comb(N,n1))
                groups.append((first,len(tables)))
        lp,_=stats.directional_logp(np.array(tables));w=np.array(weights)
        for lo,hi in groups:
            for d in range(2):
                for alpha in [.001,.01,.025,.05,.1]:self.assertLessEqual(w[lo:hi][lp[lo:hi,d]<=math.log(alpha)].sum(),alpha+1e-12)

    def test_underflow_is_finite_log_and_no_fabricated_zero(self):
        a=np.array([[2000,0,0,2000],[0,2000,2000,0],[0,10,0,10],[10,0,10,0]],dtype=np.int64);lp,n=stats.directional_logp(a)
        exact=-math.log(math.comb(4000,2000));self.assertGreaterEqual(n,2);self.assertAlmostEqual(lp[0,0],exact,places=8);self.assertAlmostEqual(lp[1,1],exact,places=8)
        self.assertTrue(np.all(np.isfinite(lp)));self.assertTrue(np.all(stats.display(lp)>0));np.testing.assert_equal(lp[2:],0.)

    def test_partial_conjunction_largest_p_and_rN(self):
        logs=np.log([[1e-9,.02,.3,.8],[.8,.2,.5,1.]])
        for r in range(1,5):
            expected=chi2.logsf(-2*np.sort(logs,axis=1)[:,r-1:].sum(axis=1),2*(5-r))
            np.testing.assert_allclose(stats.partial_conjunction(logs,r),expected,atol=1e-13)
        np.testing.assert_allclose(stats.partial_conjunction(logs,4),np.max(logs,axis=1),atol=1e-14)
        with self.assertRaises(ValueError):stats.partial_conjunction(np.array([[-.5,np.nan]]),2)

    def test_PC_recurrence_null_with_r_minus_one_true(self):
        rng=np.random.default_rng(72401);u=rng.uniform(size=(60000,10));logs=np.log(u)
        for r in [2,5,10]:
            z=logs.copy();z[:,:r-1]=-1000
            p=stats.partial_conjunction(z,r);rate=float((p<=math.log(.05)).mean());self.assertLess(rate,.054)
        # Opposite animal subsets may each establish same-direction r=2.
        plus=np.log([[1e-12,1e-12,1e-12,1,1,1]]);minus=plus[:,::-1]
        self.assertLess(stats.partial_conjunction(plus,2)[0],math.log(.05));self.assertLess(stats.partial_conjunction(minus,2)[0],math.log(.05))

    def test_BH_full_family_ties_and_direction_budget(self):
        p=np.array([.04,.001,.2,.001,1.]);q=stats.bh_log(np.log(p));expected=np.array([.04*5/3,.0025,.25,.0025,1.])
        np.testing.assert_allclose(np.exp(q),expected,rtol=1e-13);np.testing.assert_allclose(np.exp(stats.adjusted_direction_logq(q)),np.minimum(1,2*expected))
        self.assertGreater(np.exp(q[1]),.001) # Full family includes nonselected rows.
        self.assertEqual(len(stats.bh_log(np.empty(0))),0)

    def test_optimized_anchor_against_independent_sets(self):
        rng=np.random.default_rng(9);n=250;positions=np.arange(100,100+n,dtype='<u4');phase=[];av=[];te=[];raw=[];test=[]
        for i in range(10):
            records=[]
            for ps in range(1,5):
                for pos in range(90+ps*10,390,37):records.append((pos,0,2 if (pos+i)%2 else 3,1,30.,ps,0,0,2,1,1))
            # A real SNP has one PS in an animal; avoid duplicate variant keys.
            seen=set();clean=[]
            for r in records:
                k=r[:3]
                if k not in seen:seen.add(k);clean.append(r)
            a=np.array(clean,dtype=VAR);a=a[np.argsort(a['ps'],kind='stable')];phase.append(a);pa=[];pt=[];sets=[];ts=[]
            for j in range(n):
                values=set(map(int,rng.choice(np.arange(1,5),size=int(rng.integers(0,3)),replace=False)));tested={x for x in values if rng.random()<.7}
                sets.append(values);ts.append(tested);pa.extend((j<<32)|x for x in values);pt.extend((j<<32)|x for x in tested)
            av.append(np.array(sorted(pa),dtype='<u8'));te.append(np.array(sorted(pt),dtype='<u8'));raw.append(sets);test.append(ts)
        k,nt,na,nc,perf=anchors.assign(positions,av,te,phase)
        for j,p in enumerate(positions):self.assertEqual((int(k[j]),int(nt[j]),int(na[j]),int(nc[j])),brute(int(p),[x[j] for x in raw],[x[j] for x in test],phase))

    def test_distance_tie_uses_variant_key(self):
        w=np.array([90*16+3,90*16+2,110*16+1],dtype='<u8');result=anchors.nearest(w,np.array([100,109,90]))
        self.assertEqual(result.tolist(),[90*16+2,110*16+1,90*16+2])

    def test_full_fixture_all54families_resume_annotations(self):
        d=tmp('whole');cfg=fixture(d);out=d/'output';s=run.execute(cfg,out,'fixture-release',workers=2,fixture=True)
        self.assertEqual(s['families'],54);self.assertEqual(s['original_eligible_profiles'],9);self.assertEqual(s['pilot_profiles_replayed'],9)
        self.assertGreater(s['BH05_rows'],0);loc=read(out/'result_locations.json');self.assertEqual(len(loc['families']),54);self.assertGreater(len(loc['annotations']),0)
        p=Path(loc['tiles'][0]);h=sha(p/'profiles.npy');s2=run.execute(cfg,out,'fixture-release',workers=2,fixture=True);self.assertEqual(sha(p/'profiles.npy'),h);self.assertEqual(s['BH05_rows'],s2['BH05_rows'])
        run.export(out)
        cp=Path(cfg['source_inputs']['units'][0]['cache']['path'])/'chr1.npy';cp.write_bytes(cp.read_bytes()+b'corrupt')
        with self.assertRaises(ValueError):run.execute(cfg,out,'fixture-release',workers=2,fixture=True)

    def test_stop_request_no_jobs_started(self):
        d=tmp('stop');write(d/'STOP_REQUESTED',{})
        with self.assertRaises(StopRequested):parallel(str,[dict(name='never',memory_estimate=1024)],d,'test')

    def test_gate_corruption_and_memory_rejection(self):
        d=tmp('gate');out,_=begin(d,'x');(out/'value').write_text('a');seal(d,out,'x',{})
        (out/'value').write_text('b')
        with self.assertRaises(ValueError):current(d,'x')
        small=dict(cpus=32,total=GIB,available=GIB,reserve=2*GIB,usable=0)
        with patch.object(common,'resources',return_value=small):
            with self.assertRaises(ValueError):guard(d)

def main():
    global WORK
    p=argparse.ArgumentParser();p.add_argument('--record',required=True);a=p.parse_args();WORK=Path(tempfile.mkdtemp(prefix='asmD_'));start=time.monotonic()
    tests=unittest.defaultTestLoader.loadTestsFromTestCase(Tests);r=unittest.TextTestRunner(verbosity=2).run(tests)
    write(a.record,dict(status='PASS' if r.wasSuccessful() else 'FAIL',tests=r.testsRun,failures=len(r.failures),errors=len(r.errors),skipped=len(r.skipped),
        release_sha256=release() if (ROOT/'release.json').exists() else 'UNPACKAGED_DEVELOPMENT',platform=platform.platform(),python=sys.version,seconds=time.monotonic()-start,work=str(WORK)))
    raise SystemExit(0 if r.wasSuccessful() else 1)

if __name__=='__main__':main()
