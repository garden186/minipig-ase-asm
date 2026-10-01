"""Check block integration with coverage evidence and no SNP-level ASE calls."""
import unittest

from link_engine import join


class Index:
    def __init__(self, exon_transcript='T1'):
        self.exon_transcript = exon_transcript

    def payload(self, profile):
        return {'transcripts': [dict(
            transcript_id='T1' if profile == 'promoter' else self.exon_transcript,
            gene_id='G1', is_promoter=profile == 'promoter',
            is_exon=profile == 'exon', ase_exon_union_transcript_eligible=1,
            is_ensembl_canonical=1)]}

    def annotation(self, chrom, position):
        return {'annotation_profile_id': 'exon'}


class BlockIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.cpg = dict(sample='0235', tissue='Heart', chrom='2', cpg_pos1='90', ps='100',
                        annotation_gene_id='G1', annotation_gene_name='GENE',
                        annotation_profile_id='promoter', delta_M='-0.5',
                        fisher_two_sided_p='0.001', source_BH_q_within_animal='0.01',
                        source_q_eligibility='PASS', H1_methylated='3', H1_unmethylated='7',
                        H2_methylated='8', H2_unmethylated='2', H1_depth='10', H2_depth='10',
                        **{f'nominal_ge_{d}': '1' for d in (3, 5, 8, 10)})
        self.block = dict(sample='0235', tissue='Heart', gene_id='G1', gene_contig='2',
                          haplotype_block_id='B1', measurement_id='M1',
                          wgs_phase_set='100', wgs_phase_status='CONCORDANT',
                          wgs_phase_orientation='FLIPPED', n_wgs_phase_compared_snps='1',
                          n_wgs_phase_direct_snps='0', n_wgs_phase_flipped_snps='1',
                          n_wgs_phase_mismatch_snps='0', local_phase_status='PASS',
                          v026_technical_qc_status='PASS', block_snps='2_110_A_G',
                          haplotypeA='A', haplotypeB='G', wgs_balanced_gene_snps='2_110_A_G',
                          ase_a_count='21', ase_b_count='47', ase_test_status='TESTED',
                          ase_p_exact='0.001', ase_q_bh='0.01', v026_evidence_class='CORE',
                          v026_is_core_ase='1', v026_statistical_signal='1')
        self.site = dict(sample='0235', tissue='Heart', chrom='chr2', position='110',
                         variant_id='2_110_A_G', all_parent_block_ids='B1',
                         matched_promoter_gene_ids='G1', total_informative_fragments='8')

    def run_join(self, **kwargs):
        return join([self.cpg], [self.block], [self.site], kwargs.get('idx', Index()), {})

    def test_coverage_only_input_preserves_block_effect(self):
        blocks, audit, attempts = self.run_join()
        self.assertEqual(len(blocks), 1)
        row = blocks[0]
        self.assertEqual((row['RNA_H1_count'], row['RNA_H2_count']), (47, 21))
        self.assertAlmostEqual(row['RNA_fraction_difference_H1_minus_H2'], 26 / 68)
        self.assertEqual(row['direction_relation'], 'INVERSE')
        self.assertEqual(row['same_transcript_exonic_SNPs'], '2_110_A_G')
        self.assertEqual(row['compatible_transcripts'], 'T1')
        self.assertEqual(row['ASE_block_q'], '0.01')
        self.assertFalse(any(key.startswith('ASE_site_') for key in row))
        self.assertEqual(audit[0]['n_same_PS_same_transcript_blocks'], 1)
        self.assertEqual(attempts[0]['accepted'], 1)

    def test_zero_coverage_and_different_transcript_do_not_support_block(self):
        for idx, count in [(Index(), '0'), (Index('T2'), '8')]:
            with self.subTest(count=count, transcript=idx.exon_transcript):
                self.site['total_informative_fragments'] = count
                blocks, audit, _ = self.run_join(idx=idx)
                self.assertEqual(blocks[0]['same_transcript_RNA_support'], 0)
                self.assertEqual(audit[0]['n_same_PS_same_transcript_blocks'], 0)

    def test_different_parent_does_not_supply_transcript_support(self):
        self.site['all_parent_block_ids'] = 'B2'
        blocks, _, _ = self.run_join()
        self.assertEqual(blocks[0]['same_transcript_RNA_support'], 0)

    def test_phase_mismatch_or_failed_qc_rejects_block(self):
        for field, value in [('wgs_phase_set', '101'), ('v026_technical_qc_status', 'FAIL')]:
            with self.subTest(field=field):
                original = self.block[field]
                self.block[field] = value
                blocks, _, attempts = self.run_join()
                self.assertEqual(blocks, [])
                self.assertEqual(attempts[0]['accepted'], 0)
                self.block[field] = original

    def test_unbalanced_variant_is_rejected(self):
        self.block['wgs_balanced_gene_snps'] = '2_120_C_T'
        with self.assertRaisesRegex(ValueError, 'balanced gene evidence'):
            self.run_join()


if __name__ == '__main__':
    unittest.main()
