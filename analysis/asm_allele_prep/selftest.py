"""Portable scientific and execution tests; server launcher runs these again."""
from __future__ import annotations
import os
for variable in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[variable]='1'
import argparse, csv, gzip, io, json, platform, sys, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import common, indexing, pilot, run
from common import *
from indexing import *

CACHE_DTYPE=np.dtype([('chrom','u1'),('pos','<u4'),('ps','<u4'),('a','<u4'),('b','<u4'),('c','<u4'),('d','<u4'),
                     ('p','<f8'),('q_unit','<f8'),('q_tissue','<f8'),('q_all','<f8'),('source_row','<u4')])
TEST_ROOT=None

def case_dir(label):
    p=TEST_ROOT/(label+'_'+str(time.time_ns()));p.mkdir(parents=True);return p

def file_spec(p):return dict(path=str(p.resolve()),sha256=sha(p))

def make_vcf(p,sample,records):
    with gzip.open(p,'wt',encoding='utf-8',newline='') as f:
        f.write('##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')
        for pos,ref,alt,gt,gq,ps in records:f.write(f'1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:GQ:PS\t{gt}:{gq}:{ps}\n')

def fixture(root):
    ids=['0235','0242'];tissues=['Heart','Kidney','Liver'];units=[];wgs=[];expected=[]
    fa=root/'reference.fa';fa.write_bytes(b'>1\n'+b'A'*400+b'\n');fai=root/'reference.fa.fai';fai.write_text('1\t400\t3\t400\t401\n')
    reference=dict(**file_spec(fa),index_path=str(fai.resolve()),index_sha256=sha(fai))
    for s in ids:
        rec=[(90,'A','G','0|1' if s==ids[0] else '1|0',35,10),(290,'A','G','0|1' if s==ids[0] else '1|0','.' if s==ids[0] else 35,20)]
        if s==ids[0]:rec.insert(1,(190,'A','T','0|1',35,20))
        p=root/(s+'.phased.vcf.gz');make_vcf(p,s,rec)
        g=root/(s+'.genotype.vcf.gz');make_vcf(g,s,[(v[0],v[1],v[2],'0/1',35,v[5]) for v in rec])
        sp=root/(s+'.snpsplit.txt');sp.write_text(''.join(f'x\t1\t{v[0]}\tx\t{v[1]+"/"+v[2] if v[3]=="0|1" else v[2]+"/"+v[1]}\n' for v in rec))
        wgs.append(dict(sample=s,phased=file_spec(p),genotype=file_spec(g),snpsplit=file_spec(sp)))
        for t in tissues:
            cp=root/(s+'-'+t);cp.mkdir()
            entries=[(100,'10',9,1,2,8)]
            if s==ids[0] and t=='Heart':entries.extend([(100,'20',7,3,1,9),(200,'20',6,4,1,9)])
            if s==ids[1] and t=='Heart':entries.append((300,'30',9,1,2,8))
            psdict=sorted({v[1] for v in entries},reverse=True);common.write(cp/'phase_dictionary.json',psdict)
            entries.sort(key=lambda v:(v[0],psdict.index(v[1])))
            rawentries=entries+([(200,'20',1,0,0,1)] if s==ids[1] and t=='Heart' else [])
            rawentries.sort(key=lambda v:(v[0],v[1]))
            raw=root/(s+'-'+t+'.raw.tsv.gz')
            pilot.table(raw,[dict(sample=s,tissue=t,chrom=1,cpg_pos1=v[0],ps=v[1]) for v in rawentries],['sample','tissue','chrom','cpg_pos1','ps'])
            for c in range(1,19):
                a=np.array([(1,v[0],psdict.index(v[1]),*v[2:],.125,1.,1.,1.,i+1) for i,v in enumerate(entries)] if c==1 else [],dtype=CACHE_DTYPE)
                np.save(cp/f'chr{c}.npy',a,allow_pickle=False)
            common.write(cp/'summary.json',dict(rows=len(entries)))
            cachefiles={p.name:sha(p) for p in cp.iterdir() if p.is_file()}
            common.write(cp/'validation.json',dict(status='PASS_PREPARATION_STAGE',errors=[],identity='fixture',sha256=cachefiles))
            spec=dict(sample=s,tissue=t,raw=dict(**file_spec(raw),rows=len(rawentries)),annotated=dict(rows=len(entries)),
                cache=dict(path=str(cp.resolve()),gate=str((cp/'validation.json').resolve()),gate_sha256=sha(cp/'validation.json'),sha256=cachefiles))
            units.append(spec)
            for pos,ps,a,b,c,d in entries:
                state='NO_ELIGIBLE_SITE_ANCHOR' if pos==300 else 'DIFFERENT_EXACT_PS' if pos==100 and ps=='20' else 'ALIGNED'
                gt='.' if pos==300 else '0|1' if s==ids[0] else '1|0'
                alt,ref=(c/(c+d),a/(a+b)) if gt=='0|1' else (a/(a+b),c/(c+d))
                expected.append(dict(sample=s,tissue=t,chrom=1,cpg_pos1=pos,ps=ps,H1_methylated=a,H1_unmethylated=b,H2_methylated=c,H2_unmethylated=d,
                    fisher_two_sided_p=.125,BH_q_all30=1.,anchor_GT=gt,anchor_PS='.' if pos==300 else '10' if pos==100 else '20',orientation_status=state,
                    ALT_methylation=alt if state=='ALIGNED' else '.',REF_methylation=ref if state=='ALIGNED' else '.',delta_ALT_REF=alt-ref if state=='ALIGNED' else '.'))
    anchors=[dict(chrom=1,cpg_pos1=p,anchor_pos1=ap,anchor_ref='A' if ap!='.' else '.',anchor_alt='G' if ap!='.' else '.',n_testable_animals=nt,n_available_animals=na)
             for p,ap,nt,na in [(100,90,2,2),(200,290,1,2),(300,'.',0,0)]]
    ap=root/'anchors.tsv.gz';pp=root/'profiles.tsv.gz';pilot.table(ap,anchors,list(anchors[0]));pilot.table(pp,expected,list(expected[0]))
    return dict(fixture=True,ids=ids,tissues=tissues,chromosomes=[1],units=units,wgs=wgs,reference=reference,pilot_per_stratum=2,
        frozen_anchors=str(ap.resolve()),frozen_profiles=str(pp.resolve()),expected=dict(raw_rows=10,tested_rows=9,full_eligible_coordinates=3,frozen_sites=3,regression_checked=9))

class Tests(unittest.TestCase):
    def test_real_variant_oracle(self):
        """Compare against saved results from the earlier independent implementation."""
        count=0
        for r in pilot.rows(asset_path('frozen_variant_oracle.tsv.gz')):
            if r['ref'] not in BASE or r['alt'] not in BASE or r['ref']==r['alt']:continue
            p=np.array([(int(r['pos1']),BASE[r['ref']],BASE[r['alt']],gt_code(r['GT']),parse_gq(r['phased_GQ']),0 if r['PS'] in ['.','','NA'] else 1,int(r['duplicate_phased_position']))],dtype=VAR)
            g=np.array([] if r['genotype_agreement']=='MISSING_RAW_RECORD' else [(int(r['pos1']),BASE[r['ref']],BASE[r['alt']],gt_code(r['genotype_GT']),parse_gq(r['genotype_GQ']),0,int(r['duplicate_genotype_position']))],dtype=VAR)
            sp=r['snpsplit_alleles'].split('/')
            s=np.array([(int(r['pos1']),BASE[sp[0]],BASE[sp[1]],int(r['anchor_status']=='DUPLICATE_SNPSPLIT_POSITION'))] if len(sp)==2 else [],dtype=SPLIT)
            actual=join_arrays(p,g,s,np.array([BASE.get(r['reference_base'],255)]),np.array([int(r['pos1'])] if r['duplicate_genotype_position']=='1' else [],dtype='<u4'))[0]
            self.assertEqual(STATUS[int(actual['status'])],r['anchor_status'],str(r))
            self.assertEqual(['MISSING','phased_vcf','genotype_vcf_matched_fallback'][int(actual['gq_source'])],r['gq_source'])
            if r['effective_GQ']!='.':self.assertEqual(float(actual['gq']),float(r['effective_GQ']))
            count+=1
        self.assertGreater(count,1000);print('REAL_ORACLE_ROWS='+str(count))

    def test_genotype_duplicate_without_exact_key(self):
        p=np.array([(20,0,2,1,30,1,0)],dtype=VAR);s=np.array([(20,0,2,0)],dtype=SPLIT)
        result=join_arrays(p,np.empty(0,dtype=VAR),s,np.array([0]),np.array([20],dtype='<u4'))
        self.assertEqual(STATUS[int(result[0]['status'])],'DUPLICATE_GENOTYPE_POSITION')

    def test_scan_preserves_non_snv_duplicate_evidence(self):
        d=case_dir('scan');p=d/'input.vcf.gz'
        make_vcf(p,'0235',[(10,'A','G','0/1',30,1),(10,'A','T,C','1/2',30,1),(20,'A','AT','0/1',30,1),(20,'A','GT','0/1',30,1)])
        j=run.job('scan',d/'out','fixture',1024,spec=file_spec(p),role='genotype',sample='0235')
        out=Path(scan_job(j));a=np.load(out/'chr1.npy');dp=np.load(out/'chr1.duplicate_positions.npy')
        self.assertEqual(a['dup'].tolist(),[1]);self.assertEqual(dp.tolist(),[10,20])

    def test_sort_and_source_hash_rejection(self):
        d=case_dir('reject');p=d/'input.vcf.gz';make_vcf(p,'0235',[(20,'A','G','0|1',30,1),(10,'A','G','0|1',30,1)])
        j=run.job('scan',d/'out','fixture',1024,spec=file_spec(p),role='phased',sample='0235')
        with self.assertRaisesRegex(ValueError,'coordinate sorted'):scan_job(j)
        make_vcf(p,'0235',[(10,'A','G','0|1',30,1)])
        with self.assertRaisesRegex(ValueError,'hash mismatch'):scan_job(j)

    def test_gt_flip_gq_fallback_conversion_and_ps(self):
        p=np.array([(10,0,2,1,np.nan,1,0),(20,1,3,2,30,2,0),(30,0,2,1,30,0,0)],dtype=VAR)
        g=np.array([(10,0,2,3,25,0,0)],dtype=VAR);s=np.array([(10,0,2,0),(20,3,1,0),(30,0,2,0)],dtype=SPLIT)
        a=join_arrays(p,g,s,np.array([0,1,0]));self.assertEqual(a['status'].tolist(),[0,0,4])
        self.assertEqual(a['alt_hap'].tolist(),[2,1,2]);self.assertEqual(a['gq_source'].tolist(),[2,1,1]);self.assertEqual(a['conversion'].tolist(),[1,1,1])

    def test_anchor_ranking_opportunity_before_distance(self):
        a=np.array([(90,0,2,1,30,1,0,0,2,1,1),(101,0,3,1,30,1,0,0,2,1,0)],dtype=ANCHOR)
        b=np.array([(90,0,2,2,30,1,0,0,1,1,1)],dtype=ANCHOR)
        actual=pilot.choose_anchor(100,[{'x'},{'y'}],[{'x'},set()],[a,b],[{'x':1},{'y':1}])
        self.assertEqual(actual,(90*16+2,1,2,2))
        # Duplicate tissues/PS contribute one animal, never pseudo-replicates.
        self.assertEqual(pilot.choose_anchor(100,[{'x'},{'y'}],[{'x'},{'y'}],[a,b],[{'x':1},{'y':1}])[1],2)

    def test_full_fixture_resume_and_export(self):
        d=case_dir('pipeline');cfg=fixture(d);out=d/'output'
        result=run.execute(cfg,out,'fixture-release',workers=2,fixture=True)
        self.assertEqual(result['counts']['regression_checked'],9)
        loc=read(out/'result_locations.json');p=Path(loc['pilot']['chr1']);sites=list(pilot.rows(p/'pilot_sites.tsv.gz'))
        s100=next(s for s in sites if s['cpg_pos1']=='100')
        self.assertEqual((s100['N_before_Heart'],s100['N_aligned_Heart']),('2','2'))
        h=sha(p/'pilot_profiles.tsv.gz');again=run.execute(cfg,out,'fixture-release',workers=2,fixture=True)
        self.assertEqual(sha(p/'pilot_profiles.tsv.gz'),h);self.assertEqual(read(out/'result_locations.json'),loc)
        archive=run.export(out)
        with run.tarfile.open(archive) as tar:self.assertTrue(all(not n.endswith('.npy') for n in tar.getnames()))
        # A saved source cache modified after a completed run must be rejected.
        cp=Path(cfg['units'][0]['cache']['path'])/'chr1.npy';cp.write_bytes(cp.read_bytes()+b'changed')
        with self.assertRaisesRegex(ValueError,'identity mismatch'):run.execute(cfg,out,'fixture-release',workers=1,fixture=True)

    def test_checkpoint_tamper_and_identity_rejected(self):
        d=case_dir('checkpoint');out,reused=begin(d,'x');(out/'data.txt').write_text('abc');seal(d,out,'x',{})
        with self.assertRaisesRegex(ValueError,'Different checkpoint'):current(d,'y')
        (out/'data.txt').write_text('xyz')
        with self.assertRaisesRegex(ValueError,'identity mismatch'):current(d,'x')

    def test_resource_admission_and_controller_lock(self):
        d=case_dir('resources');small=dict(cpus=8,total=GIB,available=GIB,reserve=2*GIB,usable=0,hostname='fixture')
        with patch.object(common,'resources',return_value=small):
            with self.assertRaisesRegex(ValueError,'Memory reserve'):guard(d)
            with self.assertRaisesRegex(RuntimeError,'Insufficient'):parallel(str,[dict(name='x',memory_estimate=1024)],d,'fixture')
        with lock(d):
            with self.assertRaises(AlreadyRunning):
                with lock(d):pass

    def test_p_independent_selection(self):
        a=np.array([(1,10,0,9,1,1,9,.9,1,1,1,1),(1,20,0,9,1,1,9,.01,1,1,1,2)],dtype=CACHE_DTYPE)
        b=a[:1].copy();d={('s1','Heart'):a,('s2','Heart'):b};initial=pilot.independent_selection(d,['s1','s2'],1)
        a['p']=1e-100;a['a']=90000;a['q_all']=0
        changed=pilot.independent_selection(d,['s1','s2'],1)
        self.assertEqual(initial[2:],changed[2:]);self.assertEqual(initial[1].tolist(),[2,1])

def main():
    global TEST_ROOT
    parser=argparse.ArgumentParser();parser.add_argument('--record',required=True);args=parser.parse_args()
    record=Path(args.record).resolve();record.parent.mkdir(parents=True,exist_ok=True)
    # Keep fixture paths short on Windows; launcher pins TMPDIR on Linux.
    TEST_ROOT=Path(tempfile.mkdtemp(prefix='asmA_'))
    start=time.monotonic();suite=unittest.defaultTestLoader.loadTestsFromTestCase(Tests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    rh=release() if (ROOT/'release.json').exists() else 'UNPACKAGED_DEVELOPMENT'
    common.write(record,dict(status='PASS' if result.wasSuccessful() else 'FAIL',tests=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),release_sha256=rh,
        platform=platform.platform(),python=sys.version,seconds=time.monotonic()-start,work=str(TEST_ROOT)))
    raise SystemExit(0 if result.wasSuccessful() else 1)

if __name__=='__main__':main()
