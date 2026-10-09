"""Focused scientific accounting and failure-preservation tests."""
from __future__ import annotations
from collections import Counter
import csv
import gzip
import io
import itertools
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import CALLS, STATES, fingerprint, read_json, write_json, sha, signature
from profile import CORE, engine, Accumulator, assign, cpg_observations, merge_cpg, run_records, SITE_FIELDS, DETAIL_FIELDS
from validate import validate
from run import complete_region, parser, select, prepare_split, check_full_budget
import worker

def record(name='r', start=100, seq='CAAAT', xm='Z....', flag=0, mate=-1, xg='CT', cigar=None, mapq=60):
    fields = [name, str(flag), '1', str(start+1), str(mapq), cigar or str(len(seq))+'M',
              '=' if mate >= 0 else '*', str(mate+1) if mate >= 0 else '0', '0', seq, 'I'*len(seq),
              'XM:Z:'+xm, 'XG:Z:'+xg, 'XR:Z:CT']
    result = CORE.parse_sam_line('\t'.join(fields))
    result.audit_quality = fields[10]
    return result

def samline(r):
    return '\t'.join([r.qname, str(r.flag), r.rname, str(r.pos0+1), str(r.mapq),r.cigar,
        r.mate_rname,str(r.mate_pos0+1 if r.mate_pos0>=0 else 0),str(r.tlen),r.seq,r.audit_quality,
        'XM:Z:'+r.tags['XM'],'XG:Z:'+r.tags['XG'],'XR:Z:CT'])+'\n'

def classify(records, links, contexts=None):
    contexts = contexts if contexts is not None else {p:{'ELIGIBLE_PHASED_LINK'} for p in links}
    return assign(records, links, sorted(links), engine.conversion_core(links, Counter()), sorted(contexts), contexts)

class ProfileTests(unittest.TestCase):
    def test_no_snp_is_retained(self):
        a = Accumulator('0235','Heart','1',[100],[0],{}, {})
        a.consume((record(),))
        self.assertEqual(a.sites[100][STATES.index('NO_KNOWN_SNV_OVERLAP')*4],1)
        self.assertEqual(a.qc['FRAGMENT_cpg_observations'],1)

    def test_single_phased_snp_assigns(self):
        l={104:CORE.SNPInfo(104,'p','T','A')}
        self.assertEqual(classify((record(),),l)['hap'],'H1')

    def test_two_consistent_snps(self):
        l={103:CORE.SNPInfo(103,'p','A','T'),104:CORE.SNPInfo(104,'p','T','A')}
        self.assertEqual(classify((record(),),l)['state'],'ASSIGNED')

    def test_conflicting_haplotypes(self):
        l={103:CORE.SNPInfo(103,'p','A','T'),104:CORE.SNPInfo(104,'p','A','T')}
        self.assertEqual(classify((record(),),l)['state'],'CONFLICTING_HAPLOTYPES')

    def test_multiple_phase_sets(self):
        l={103:CORE.SNPInfo(103,'p','A','T'),104:CORE.SNPInfo(104,'q','T','A')}
        self.assertEqual(classify((record(),),l)['state'],'MULTIPLE_PHASE_SETS')

    def test_unphased_snp(self):
        self.assertEqual(classify((record(),),{}, {104:{'UNPHASED_HETEROZYGOUS'}})['state'],'NO_ELIGIBLE_PHASED_LINK')

    def test_unknown_base_at_link(self):
        l={104:CORE.SNPInfo(104,'p','T','A')}
        self.assertEqual(classify((record(seq='CAAAN'),),l)['state'],'NO_USABLE_LINK_BASE')

    def test_nonallelic_base(self):
        l={104:CORE.SNPInfo(104,'p','T','A')}
        self.assertEqual(classify((record(seq='CAAAC'),),l)['state'],'NO_ALLELE_SUPPORT')

    def test_conversion_rules_match_pinned_engine(self):
        for h1,h2 in itertools.permutations('ACGT',2):
            if frozenset((h1,h2)) in CORE.CONVERSION_AMBIGUOUS:
                continue
            for base,xg in itertools.product('ACGT',('CT','GA')):
                expected=engine.allele_support(base,h1,h2,xg)
                result=classify((record(seq='CAAA'+base,xg=xg),),{104:CORE.SNPInfo(104,'p',h1,h2)})
                self.assertEqual(result['hap'],'.' if expected is None else ('H1','H2')[expected])

    def test_mate_snp_conflict(self):
        l={104:CORE.SNPInfo(104,'p','T','A')}
        self.assertEqual(classify((record(),record(seq='CAAAA')),l)['state'],'MATE_SNP_CONFLICT')

    def test_conversion_before_mate_merge(self):
        l={104:CORE.SNPInfo(104,'p','C','A')}
        self.assertEqual(classify((record(seq='CAAAC'),record(seq='CAAAT')),l)['hap'],'H1')

    def test_ga_strand_reference_coordinate(self):
        self.assertEqual(cpg_observations(record(start=101,xg='GA'),[100]),{100:'METHYLATED'})

    def test_cigar_deletion_not_coverage(self):
        r=record(seq='CT',xm='Zz',cigar='1M1D1M')
        self.assertEqual(cpg_observations(r,[100,101,102]),{100:'METHYLATED',102:'UNMETHYLATED'})

    def test_soft_clip_not_coverage(self):
        r=record(seq='ACA',xm='zZ.',cigar='1S2M')
        self.assertEqual(cpg_observations(r,[100]),{100:'METHYLATED'})

    def test_uncallable_is_not_unmethylated(self):
        self.assertEqual(cpg_observations(record(xm='.....'),[100]),{100:'UNCALLABLE'})

    def test_mate_cpg_merge(self):
        self.assertEqual(merge_cpg([{1:'METHYLATED'},{1:'UNMETHYLATED'}]),{1:'MATE_CPG_CONFLICT'})
        self.assertEqual(merge_cpg([{1:'UNCALLABLE'},{1:'METHYLATED'}]),{1:'METHYLATED'})
        self.assertEqual(merge_cpg([{1:'METHYLATED'},{1:'METHYLATED'}]),{1:'METHYLATED'})

    def test_mate_only_assignment(self):
        l={110:CORE.SNPInfo(110,'p','A','T')}
        a=Accumulator('0235','Heart','1',[100],[0],l,{110:{'ELIGIBLE_PHASED_LINK'}})
        a.consume((record(flag=65,mate=110),record(start=110,seq='A',xm='.',flag=129,mate=100)))
        self.assertEqual(a.sites[100][-1],1)
        self.assertEqual(a.phases[(100,'p')][8],1)

    def test_overlapping_mates_not_double_counted(self):
        a=Accumulator('0235','Heart','1',[100],[0],{}, {})
        a.consume((record(),record()))
        self.assertEqual(a.qc['READ_cpg_observations'],2)
        self.assertEqual(a.qc['FRAGMENT_cpg_observations'],1)

    def test_fetch_boundary_is_visible(self):
        a=Accumulator('0235','Heart','1',[100],[0],{}, {},fetch=CORE.Region('1',0,200))
        a.consume((record(flag=65,mate=1000),))
        offset=(len(STATES)+STATES.index('INCOMPLETE_FETCH_CONTEXT'))*4
        self.assertEqual(a.sites[100][offset],1)

    def test_disrupted_cpg_retained_as_flag(self):
        a=Accumulator('0235','Heart','1',[100],[1],{}, {})
        a.consume((record(),))
        self.assertIn(100,a.sites)

    def test_missing_tags_fail(self):
        r=record()
        del r.tags['XM']
        with self.assertRaises(ValueError):
            cpg_observations(r,[100])

    def test_duplicate_alignment_filtered(self):
        a=Accumulator('0235','Heart','1',[100],[0],{}, {})
        run_records([samline(record(flag=1024))],a,10,10,8)
        self.assertFalse(a.sites)

    def test_record_budget_stop(self):
        a=Accumulator('0235','Heart','1',[100],[0],{}, {})
        with self.assertRaises(RuntimeError):
            run_records([samline(record())],a,0,10,8)

    def test_unsorted_stream_fails(self):
        a=Accumulator('0235','Heart','1',[100],[0],{}, {})
        with self.assertRaises(ValueError):
            run_records([samline(record(start=101)),samline(record())],a,10,10,8)

    def test_membership_and_phase_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            task={'sample':'0235','tissue':'Heart','chrom':'1','task_id':'test','start0':0,'end0':200}
            write_json(d/'task.json',task)
            with gzip.open(d/'observation_detail.tsv.gz','wt',encoding='utf-8',newline='') as f:
                w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(DETAIL_FIELDS)
                l={104:CORE.SNPInfo(104,'p','T','A')}
                a=Accumulator('0235','Heart','1',[100],[0],l,{104:{'ELIGIBLE_PHASED_LINK'}},detail=w)
                a.consume((record(),))
                a.consume((record(name='b',seq='TAAAA',xm='z....'),))
                a.consume((record(name='c',seq='TAAAN',xm='z....'),))
                a.write(d)
            gate=validate(d)
            self.assertEqual(gate['membership_rows_reconciled'],6)
            write_json(d/'timing.json',{'total_seconds':1})
            self.assertTrue(complete_region(d,'test'))
            with open(d/'entity_assignment_context.tsv','a',encoding='utf-8') as f:
                f.write('tampered\n')
            with self.assertRaises(ValueError):
                complete_region(d,'test')

    def test_same_cpg_different_ps_stays_separate(self):
        links={102:CORE.SNPInfo(102,'p','A','T'),110:CORE.SNPInfo(110,'q','T','A')}
        a=Accumulator('0235','Heart','1',[105],[0],links,{p:{'ELIGIBLE_PHASED_LINK'} for p in links})
        a.consume((record(start=102,seq='AAAC',xm='...Z'),))
        a.consume((record(start=105,seq='CAAAAT',xm='Z.....'),))
        self.assertEqual(set(a.phases),{(105,'p'),(105,'q')})

    def test_no_methylated_site_is_still_retained(self):
        a=Accumulator('0235','Heart','1',[100],[0],{}, {})
        a.consume((record(seq='TAAAT',xm='z....'),))
        self.assertIn(100,a.sites)

    def test_synthetic_sam_vcf_worker_and_replay(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            header='##fileformat=VCFv4.2\n##contig=<ID=1,length=200>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t0235\n'
            variants='1\t105\t.\tT\tA\t60\tPASS\t.\tGT:PS\t0|1:101\n1\t119\t.\tC\tT\t60\tPASS\t.\tGT:PS\t0/1:.\n'
            for role in ('phased_vcf','genotype_vcf'):
                (root/(role+'.vcf')).write_text(header+variants,encoding='utf-8')
            (root/'split.tsv').write_text('1\t1\t105\t1\tT/A\n2\t1\t119\t1\tC/T\n',encoding='utf-8')
            seq=list('A'*200)
            for p in (100,105):
                seq[p:p+2]=list('CG')
            (root/'ref.fa').write_text('>1\n'+''.join(seq)+'\n',encoding='utf-8')
            records=[record(name='a',seq='CGAATTGAAAT',xm='Z....z.....'),
                     record(name='b',seq='TGAAACGAAAT',xm='z....Z.....'),
                     record(name='c',start=105,seq='CG',xm='Z.')]
            (root/'input.sam').write_text(''.join(samline(r) for r in records),encoding='utf-8')
            paths={'bam':root/'input.sam','phased_vcf':root/'phased_vcf.vcf','genotype_vcf':root/'genotype_vcf.vcf',
                   'snpsplit_snps':root/'split.tsv','reference':root/'ref.fa'}
            task={'sample':'0235','tissue':'Heart','chrom':'1','start0':90,'end0':130,'chrom_length':200,
                  'inputs':{k:fingerprint(v) for k,v in paths.items()},'split_subset':str(root/'split.tsv'),
                  'split_subset_sha256':sha(root/'split.tsv'),'tools':{'samtools':'synthetic_samtools','bcftools':'synthetic_bcftools'},
                  'write_detail':True,'samtools_threads':0,'max_seconds':30,'max_records':100,'max_memory_gib':8,'task_id':'fixture'}
            def fake_stream(command,error_path,deadline=None):
                Path(error_path).write_text('',encoding='utf-8')
                source=Path(command[-1]) if command[0]=='synthetic_bcftools' else root/'input.sam'
                yield from (line for line in source.read_text(encoding='utf-8').splitlines(True) if not line.startswith('#'))
            def fasta_positions(**kwargs):
                region=kwargs['region']
                text=''.join(seq)
                return [p for p in range(region.start0,region.end0-1) if text[p:p+2]=='CG'], {'fixture':1}
            original_stream=CORE.stream_command
            try:
                with patch.object(worker,'stream',fake_stream),patch.object(CORE,'enumerate_shared_reference_cpgs',fasta_positions):
                    worker.run(task,root/'first')
                    worker.run(task,root/'second')
            finally:
                CORE.stream_command=original_stream
            gate=read_json(root/'first'/'validation.json')
            self.assertEqual(gate['covered_cpg_rows'],2)
            self.assertEqual(gate['exact_ps_rows'],2)
            self.assertEqual(gate['membership_rows_reconciled'],6)
            for name in ('cpg_assignment.tsv.gz','cpg_haplotype_by_ps.tsv.gz','entity_assignment_context.tsv','measurement_summary.tsv'):
                opener=gzip.open if name.endswith('.gz') else open
                with opener(root/'first'/name,'rt',encoding='utf-8') as f:
                    first=f.read()
                with opener(root/'second'/name,'rt',encoding='utf-8') as f:
                    self.assertEqual(first,f.read())
            self.assertTrue(complete_region(root/'first','fixture'))
            with self.assertRaises(FileExistsError):
                worker.run(task,root/'first')

    def test_no_coverage_validation_is_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            write_json(d/'task.json',{'sample':'0235','tissue':'Heart','chrom':'1','task_id':'empty'})
            a=Accumulator('0235','Heart','1',[100],[0],{}, {})
            a.write(d)
            result=validate(d)
            self.assertEqual(result['covered_cpg_rows'],0)
            self.assertEqual(read_json(d/'count_summary.json')['reference_cpgs_with_no_retained_alignment'],1)

    def test_partial_snp_preparation_is_preserved_on_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            source=d/'input.tsv'
            source.write_text('1\t1\t105\t1\tT/A\n',encoding='utf-8')
            unit={'sample':'0235','inputs':{'snpsplit_snps':fingerprint(source)}}
            partial=d/'prepared_snps/0235'
            partial.mkdir(parents=True)
            (partial/'preserve.txt').write_text('old partial evidence',encoding='utf-8')
            ready=prepare_split(unit,['1'],d)
            self.assertNotEqual(ready,partial)
            self.assertEqual((partial/'preserve.txt').read_text(encoding='utf-8'),'old partial evidence')
            self.assertEqual(prepare_split(unit,['1'],d),ready)

    def test_full_resource_gate_checks_all_three_budgets(self):
        request={'max_seconds':100,'max_output_bytes':100,'max_memory_gib':1}
        plan={'forecast_seconds_with_factor_2':10,'forecast_output_bytes_with_factor_2':10,
              'forecast_largest_chromosome_memory_bytes':10}
        check_full_budget(plan,request)
        for field in plan:
            changed=dict(plan);changed[field]=2**40
            with self.assertRaises(RuntimeError):
                check_full_budget(changed,request)

if __name__=='__main__':
    unittest.main()
