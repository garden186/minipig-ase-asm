"""Scientific regression tests for raw evidence and post-extraction GQ changes."""
import csv
import gzip
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from test_profile import record, samline
from profile import CORE, Accumulator
from common import CALLS, STATES, rows, read_json, write_json, fingerprint, sha
from evidence import save_catalog, reconstruct, numeric_gq, genotype_alleles
from reassign import region, write_asm_inputs
from validate import validate
import worker


def variants(links, gqs, raw_gqs=None):
    output=[]
    for p,v in links.items():
        for source in ('phased_vcf','genotype_vcf'):
            value = gqs[p] if source=='phased_vcf' or raw_gqs is None else raw_gqs[p]
            output.append(dict(source=source,chrom='1',pos1=p+1,ref=v.h1,alt=v.h2,QUAL='60',FILTER='PASS',
                GT='0|1' if source=='phased_vcf' else '0/1',GQ='.' if value is None else str(value),
                DP='30',AD='15,15',PL='.',PS=v.ps if source=='phased_vcf' else '.',eligible_link=1))
    return output


def fixture(path, links, gqs, groups, raw_gqs=None, disrupted=0, fetch=None):
    path.mkdir()
    task=dict(sample='0235',tissue='Heart',chrom='1',start0=90,end0=150,task_id='fixture')
    write_json(path/'task.json',task)
    contexts={p:{'ELIGIBLE_PHASED_LINK'} for p in links}
    save_catalog(path,variants(links,gqs,raw_gqs),links,contexts)
    with gzip.open(path/'fragment_evidence.jsonl.gz','wt',encoding='utf-8') as f:
        acc=Accumulator('0235','Heart','1',[100],[disrupted],links,contexts,evidence=f,fetch=fetch)
        for group in groups: acc.consume(group)
        acc.write(path)
    write_asm_inputs(path)
    validate(path)
    return path


class EvidenceTests(unittest.TestCase):
    def test_low_gq_conflict_removal_rescues_assignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            links={103:CORE.SNPInfo(103,'p','A','T'),104:CORE.SNPInfo(104,'p','A','T')}
            src=fixture(root/'src',links,{103:40,104:5},[(record(),)])
            result=region(src,root/'filtered',20)
            r=next(rows(root/'filtered/cpg_asm_inputs.tsv.gz'))
            self.assertEqual((r['H1_methylated'],r['H2_methylated']),('1','0'))
            self.assertEqual(next(rows(root/'filtered/fragment_assignment_transitions.tsv')),
                {'baseline_state':'CONFLICTING_HAPLOTYPES','filtered_state':'ASSIGNED','fragments':'1'})
            self.assertEqual(result['retained_links'],1)

    def test_remaining_high_gq_link_preserves_fragment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            links={103:CORE.SNPInfo(103,'p','A','T'),104:CORE.SNPInfo(104,'p','T','A')}
            src=fixture(root/'src',links,{103:40,104:5},[(record(),)])
            region(src,root/'filtered',20)
            self.assertEqual(next(rows(root/'filtered/cpg_asm_inputs.tsv.gz'))['H1_methylated'],'1')

    def test_all_links_removed_keep_total_methylation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            src=fixture(root/'src',{104:CORE.SNPInfo(104,'p','T','A')},{104:5},[(record(),)])
            region(src,root/'filtered',20)
            r=next(rows(root/'filtered/cpg_assignment.tsv.gz'))
            self.assertEqual(r['FRAGMENT_NO_ELIGIBLE_PHASED_LINK_METHYLATED'],'1')
            self.assertEqual(list(rows(root/'filtered/cpg_asm_inputs.tsv.gz')),[])

    def test_multiple_ps_can_be_resolved_by_gq_without_merging(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            src=fixture(root/'src',{103:CORE.SNPInfo(103,'p','A','T'),104:CORE.SNPInfo(104,'q','T','A')},
                        {103:40,104:5},[(record(),)])
            region(src,root/'filtered',20)
            self.assertEqual(next(rows(root/'filtered/cpg_asm_inputs.tsv.gz'))['ps'],'p')

    def test_missing_gq_is_distinct_from_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            src=fixture(root/'src',{104:CORE.SNPInfo(104,'p','T','A')},{104:None},[(record(),)])
            self.assertEqual(region(src,root/'none',None)['retained_links'],1)
            self.assertEqual(region(src,root/'zero',0)['retained_links'],0)

    def test_exact_matched_genotype_gq_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            src=fixture(root/'src',{104:CORE.SNPInfo(104,'p','T','A')},{104:None},[(record(),)],raw_gqs={104:30})
            c=read_json(src/'link_catalog.json')['links']['104']
            self.assertEqual(c['gq_source'],'genotype_vcf_matched_fallback')
            self.assertEqual(region(src,root/'filtered',30)['retained_links'],1)

    def test_nonmatching_raw_genotype_cannot_supply_gq(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); links={104:CORE.SNPInfo(104,'p','T','A')}
            v=variants(links,{104:None},{104:40});v[1]['GT']='1/1'
            c=save_catalog(root,v,links,{104:{'ELIGIBLE_PHASED_LINK'}})
            self.assertIsNone(c['104']['gq'])
            self.assertFalse(c['104']['genotype_alleles_agree'])

    def test_mate_conflict_and_exact_raw_bases_are_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            src=fixture(root/'src',{104:CORE.SNPInfo(104,'p','T','A')},{104:40},[(record(),record(seq='CAAAA'))])
            with gzip.open(src/'fragment_evidence.jsonl.gz','rt') as f: e=json.loads(next(f))
            self.assertEqual([r['snps'][0][1] for r in e['reads']],['T','A'])
            self.assertEqual(e['reads'][0]['snps'][0][2:], [40,4])
            self.assertEqual(e['baseline']['state'],'MATE_SNP_CONFLICT')

    def test_converted_bases_normalize_before_conflict(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            src=fixture(root/'src',{104:CORE.SNPInfo(104,'p','C','A')},{104:40},[(record(seq='CAAAC'),record(seq='CAAAT'))])
            self.assertEqual(read_json(src/'evidence_validation.json')['fragment_rows'],1)
            region(src,root/'filtered',20)
            self.assertEqual(next(rows(root/'filtered/cpg_asm_inputs.tsv.gz'))['H1_methylated'],'1')

    def test_mate_only_link_and_cpg_conflict_replay(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            a=record(flag=65,mate=110)
            b=record(start=110,seq='A',xm='.',flag=129,mate=100)
            src=fixture(root/'src',{110:CORE.SNPInfo(110,'p','A','T')},{110:40},[(a,b),(record(name='b'),record(name='b',xm='z....'))])
            region(src,root/'filtered',20)
            r=next(rows(root/'filtered/cpg_assignment.tsv.gz'))
            self.assertEqual(r['fragment_assigned_using_mate_only'],'1')
            self.assertEqual(r['FRAGMENT_NO_KNOWN_SNV_OVERLAP_MATE_CPG_CONFLICT'],'1')

    def test_empty_evidence_and_no_snp_fragments_are_valid(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            empty=fixture(root/'empty',{}, {},[])
            self.assertEqual(region(empty,root/'empty_replay',20)['fragment_rows'],0)
            src=fixture(root/'src',{}, {},[(record(xm='.....'),)])
            region(src,root/'filtered',20)
            self.assertEqual(next(rows(root/'filtered/cpg_assignment.tsv.gz'))['FRAGMENT_NO_KNOWN_SNV_OVERLAP_UNCALLABLE'],'1')

    def test_genotype_mask_is_not_changed_by_link_gq(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            src=fixture(root/'src',{104:CORE.SNPInfo(104,'p','T','A')},{104:40},[(record(),)],disrupted=1)
            region(src,root/'filtered',20)
            self.assertEqual(next(rows(root/'filtered/cpg_asm_inputs.tsv.gz'))['genotype_disrupted'],'1')

    def test_source_mutation_and_overwrite_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);src=fixture(root/'src',{}, {},[(record(),)])
            region(src,root/'filtered',20)
            with self.assertRaises(FileExistsError): region(src,root/'filtered',20)
            (src/'fragment_evidence.jsonl.gz').write_bytes(b'changed')
            with self.assertRaises(ValueError): region(src,root/'different',20)

    def test_numeric_gq_validation(self):
        for v in (None,'.','NaN','inf','-1','bad'): self.assertIsNone(numeric_gq(v))
        self.assertEqual(numeric_gq('0'),0)

    def test_duplicate_variant_provenance_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            links={104:CORE.SNPInfo(104,'p','T','A')};v=variants(links,{104:40})
            with self.assertRaises(ValueError): save_catalog(Path(tmp),v+v,links,{})

    def test_ga_coordinate_and_soft_clip_atoms(self):
        from evidence import atoms
        r=record(start=101,seq='ACA',xm='.Z.',xg='GA',cigar='1S2M')
        self.assertEqual(atoms(r,[100],1,True),[[100,'C',40,1,'Z']])

    def test_missing_mate_context_is_not_repaired_by_gq(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);links={104:CORE.SNPInfo(104,'p','T','A')}
            src=fixture(root/'src',links,{104:40},[(record(flag=65,mate=1000),)],fetch=CORE.Region('1',90,150))
            region(src,root/'filtered',20)
            self.assertEqual(next(rows(root/'filtered/cpg_assignment.tsv.gz'))['FRAGMENT_INCOMPLETE_FETCH_CONTEXT_METHYLATED'],'1')

    def test_phased_gq_precedence_is_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);links={104:CORE.SNPInfo(104,'p','T','A')}
            c=save_catalog(root,variants(links,{104:5},{104:40}),links,{})['104']
            self.assertEqual((c['gq'],c['gq_source'],c['genotype_gq']),(5,'phased_vcf',40))

    @unittest.skipUnless(shutil.which('samtools') and shutil.which('bcftools'), 'Real samtools/bcftools are required on the server')
    def test_actual_bam_vcf_tools_and_gq_reassignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            samtools,bcftools=shutil.which('samtools'),shutil.which('bcftools')
            sequence='A'*100+'CGAAT'+'A'*95
            (root/'ref.fa').write_text('>1\n'+sequence+'\n',encoding='utf-8')
            subprocess.run([samtools,'faidx',str(root/'ref.fa')],check=True,capture_output=True)
            header='@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:200\n'
            (root/'input.sam').write_text(header+samline(record(seq='CGAAT')),encoding='utf-8')
            subprocess.run([samtools,'sort','-T',str(root/'sort_work'),'-o',str(root/'input.bam'),str(root/'input.sam')],check=True,capture_output=True)
            subprocess.run([samtools,'index',str(root/'input.bam')],check=True,capture_output=True)
            header=('##fileformat=VCFv4.2\n##contig=<ID=1,length=200>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t0235\n')
            variant='1\t105\t.\tT\tA\t60\tPASS\t.\tGT:PS:GQ\t0|1:101:5\n'
            (root/'variants.vcf').write_text(header+variant,encoding='utf-8')
            subprocess.run([bcftools,'view','-Oz','-o',str(root/'variants.vcf.gz'),str(root/'variants.vcf')],check=True,capture_output=True)
            subprocess.run([bcftools,'index','-t',str(root/'variants.vcf.gz')],check=True,capture_output=True)
            (root/'split.tsv').write_text('1\t1\t105\t1\tT/A\n',encoding='utf-8')
            paths={'bam':root/'input.bam','reference':root/'ref.fa','phased_vcf':root/'variants.vcf.gz',
                   'genotype_vcf':root/'variants.vcf.gz','snpsplit_snps':root/'split.tsv'}
            task=dict(sample='0235',tissue='Heart',chrom='1',start0=90,end0=130,chrom_length=200,
                inputs={k:fingerprint(v) for k,v in paths.items()},split_subset=str(root/'split.tsv'),
                split_subset_sha256=sha(root/'split.tsv'),tools=dict(samtools=samtools,bcftools=bcftools),
                write_detail=True,samtools_threads=0,max_seconds=30,max_records=100,max_memory_gib=8,task_id='actual_tools')
            original=CORE.stream_command
            try: worker.run(task,root/'actual')
            finally: CORE.stream_command=original
            self.assertEqual(region(root/'actual',root/'filtered',20)['retained_links'],0)
            self.assertEqual(read_json(root/'actual/evidence_validation.json')['fragment_rows'],1)


if __name__=='__main__': unittest.main()
