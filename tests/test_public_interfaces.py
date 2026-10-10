import importlib.util,json,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
prep=load('preprocessing','preprocessing/run.py')
select=load('select_final','scripts/select_primary_tissue_pairs.py')
phase=load('phase_helpers','analysis/integration/phase_helpers.py')
sys.path.insert(0,str(ROOT/'analysis/asm_individual'))
infer=load('individual_inference','analysis/asm_individual/infer.py')
class Tests(unittest.TestCase):
    def config(self,d):return dict(work_dir=str(Path(d)/'work space'),reference_fasta=str(Path(d)/'ref.fa'),annotation_gtf=str(Path(d)/'annot.gtf'),phaser_script=str(Path(d)/'phaser.py'))
    def samples(self):return [dict(sample='S01',assay=a,tissue='Heart',read1='input R1.fastq.gz',read2='input R2.fastq.gz',optical_duplicate_distance='2500') for a in ['WGS','RNA','WGBS']]
    def test_all_preprocessing_stages(self):
        with tempfile.TemporaryDirectory() as d:
            for stage in ['reference','wgs','rna','phaser','wgbs']:
                steps=prep.plan(self.config(d),self.samples(),stage,'S01');self.assertTrue(steps);self.assertTrue(all(s['command'] for s in steps))
    def test_wgs_optical_value_is_explicit(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.samples();r[0]['optical_duplicate_distance']=''
            with self.assertRaisesRegex(ValueError,'optical_duplicate_distance'):prep.plan(self.config(d),r,'wgs','S01')
    def test_final_wgbs_stops_before_haplotype_split(self):
        with tempfile.TemporaryDirectory() as d:
            steps=prep.plan(self.config(d),self.samples(),'wgbs','S01');names=[s['name'] for s in steps]
            self.assertEqual(names[-2:],['coordinate-sort','ASM-BAM-index']);self.assertNotIn('SNPsplit',names)
    def test_phaser_wgs_counts_and_rna_exclusion(self):
        with tempfile.TemporaryDirectory() as d:
            steps=prep.plan(self.config(d),self.samples(),'phaser','S01');calls={s['name']:s['command'] for s in steps}
            self.assertIn('--haplo_count_bam_exclude 1',calls['phASER-rna']);self.assertNotIn('--haplo_count_bam_exclude',calls['phASER-wgs'])
    def test_same_ps_phase_and_direction(self):self.assertEqual(phase.self_test(),8)
    def test_final_pair_selection_without_omnibus(self):
        self.assertTrue(select.selected(dict(tissue_enrichment_q_global_BH='.04',is_primary_recurrent='1',gene_omnibus_q='1')))
        self.assertFalse(select.selected(dict(tissue_enrichment_q_global_BH='.06',is_primary_recurrent='1')))
    def test_incomplete_or_overlapping_genome_rejected(self):
        with self.assertRaisesRegex(ValueError,'Truncated'):infer.coverage([dict(chrom='1',start0=0,end0=90)],{'1':100})
        with self.assertRaisesRegex(ValueError,'Gap/overlap'):infer.coverage([dict(chrom='1',start0=0,end0=60),dict(chrom='1',start0=50,end0=100)],{'1':100})
        infer.coverage([dict(chrom='1',start0=0,end0=50),dict(chrom='1',start0=50,end0=100)],{'1':100})
    def test_probability_ordered_fisher_and_full_bh(self):
        from scipy.stats import fisher_exact
        import numpy as np
        for counts in [(8,2,1,9),(0,10,8,2),(3,0,7,0),(20,55,45,10)]:
            self.assertAlmostEqual(infer.fisher(counts)[0],fisher_exact([counts[:2],counts[2:]]).pvalue,places=12)
        np.testing.assert_allclose(infer.bh([.001,.02,.02,.8]),[.004,.02666666666666667,.02666666666666667,.8])
if __name__=='__main__':unittest.main()
