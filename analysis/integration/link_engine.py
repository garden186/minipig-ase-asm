"""Haplotype-block join with original phase and transcript-support rules."""
from collections import defaultdict
import phase_helpers as L
require=L.require

def cpgkey(r):return (r['sample'],r['tissue'],r['chrom'].removeprefix('chr'),r['cpg_pos1'],r['ps'])

def significant_block(b):
    return b['ase_test_status']=='TESTED' and float(b['ase_q_bh'])<=.05

def ranking_fields(delta,h1,h2):
    r=L.effects(h1,h2,delta)
    r['abs_delta_M']=abs(delta)
    r['abs_RNA_fraction_difference']=abs((h1-h2)/(h1+h2)) if h1+h2 else ''
    r['effect_magnitude_product']=abs(delta)*abs((h1-h2)/(h1+h2)) if h1+h2 else ''
    return r

def join(cpgs,blocks,sites,idx,genes):
    bg=defaultdict(list);sg=defaultdict(list)
    for b in blocks:bg[b['sample'],b['tissue'],b['gene_id']].append(b)
    for s in sites:
        for gene in L.split(s['matched_promoter_gene_ids']):sg[s['sample'],s['tissue'],gene].append(s)
    blocks_out=[];audit=[];attempts=[]
    exon_cache={}
    def base(c,b):
        return dict(sample=c['sample'],tissue=c['tissue'],minimum_depth_each_haplotype=3,gene_id=c['annotation_gene_id'],gene_name=c['annotation_gene_name'],
            chrom=c['chrom'],ps=c['ps'],cpg_pos1=c['cpg_pos1'],ASM_P=c['fisher_two_sided_p'],source_ASM_q=c['source_BH_q_within_animal'],source_ASM_q_eligibility=c['source_q_eligibility'],
            ASM_H1_methylated=c['H1_methylated'],ASM_H1_unmethylated=c['H1_unmethylated'],ASM_H2_methylated=c['H2_methylated'],ASM_H2_unmethylated=c['H2_unmethylated'],
            ASM_H1_depth=c['H1_depth'],ASM_H2_depth=c['H2_depth'],ASM_delta_H1_minus_H2=c['delta_M'],
            ASE_block_id=b['haplotype_block_id'],ASE_measurement_id=b['measurement_id'],ASE_block_P=b['ase_p_exact'],ASE_block_q=b['ase_q_bh'],
            ASE_block_test_status=b['ase_test_status'],ASE_block_class=b['v026_evidence_class'],ASE_block_core=b['v026_is_core_ase'],
            ASE_block_FDR=int(significant_block(b)),ASE_block_statistical_signal=b['v026_statistical_signal'],
            ASE_gene_status=genes.get((c['sample'],c['tissue'],c['annotation_gene_id']),{}).get('gene_ase_status_v026','NO_GENE_SUMMARY'),
            phase_orientation=b['wgs_phase_orientation'],phase_compared_SNPs=b['n_wgs_phase_compared_snps'],phase_direct_SNPs=b['n_wgs_phase_direct_snps'],
            phase_flipped_SNPs=b['n_wgs_phase_flipped_snps'],phase_mismatch_SNPs=b['n_wgs_phase_mismatch_snps'],local_phase_status=b['local_phase_status'],
            phase_support='HIGH' if int(b['n_wgs_phase_compared_snps'])>=2 else 'SINGLETON',phase_evidence='SHA256_VERIFIED_V026_SOURCE_BLOCK',fresh_WGS_VCF_replay=False)
    for c in cpgs:
        key=(c['sample'],c['tissue'],c['annotation_gene_id']);relevant=bg[key]
        linked=[b for b in relevant if b['v026_technical_qc_status']=='PASS' and L.phase_ok(b,c['chrom'],c['ps'])]
        for b in relevant:
            attempts.append(dict(sample=c['sample'],tissue=c['tissue'],chrom=c['chrom'],cpg_pos1=c['cpg_pos1'],ps=c['ps'],gene_id=key[2],
                ASE_block_id=b['haplotype_block_id'],RNA_PS=b['wgs_phase_set'],RNA_chrom=b['gene_contig'],technical_qc=b['v026_technical_qc_status'],
                local_phase_status=b['local_phase_status'],phase_status=b['wgs_phase_status'],phase_orientation=b['wgs_phase_orientation'],
                n_compared=b['n_wgs_phase_compared_snps'],n_direct=b['n_wgs_phase_direct_snps'],n_flipped=b['n_wgs_phase_flipped_snps'],n_mismatch=b['n_wgs_phase_mismatch_snps'],
                exact_PS_and_orientation=int(L.phase_ok(b,c['chrom'],c['ps'])),accepted=int(b in linked)))
        promoters={t['transcript_id']:t for t in idx.payload(c['annotation_profile_id'])['transcripts']
            if t['gene_id']==key[2] and t['is_promoter'] and t['ase_exon_union_transcript_eligible']}
        support=defaultdict(set);canonical=defaultdict(set);transcripts=defaultdict(set)
        for s in sg[key]:
            if s['chrom'].removeprefix('chr')!=c['chrom'].removeprefix('chr'):continue
            ek=(s['chrom'].removeprefix('chr'),s['position'],key[2])
            if ek not in exon_cache:
                annotation=idx.annotation(s['chrom'],int(s['position']))
                exon_cache[ek]={t['transcript_id']:t for t in idx.payload(annotation['annotation_profile_id'])['transcripts']
                    if t['gene_id']==key[2] and t['is_exon'] and t['ase_exon_union_transcript_eligible']}
            exons=exon_cache[ek];common=promoters.keys() & exons.keys()
            if not common:continue
            for b in linked:
                bid=b['haplotype_block_id']
                if bid not in L.split(s['all_parent_block_ids']):continue
                require(s['variant_id'] in L.split(b['wgs_balanced_gene_snps']),'RNA SNP absent from balanced gene evidence')
                L.variant_alleles(b,s['variant_id'])
                for tid in sorted(common):
                    t=promoters[tid]
                    if int(s['total_informative_fragments'])>0:
                        support[bid].add(s['variant_id']);transcripts[bid].add(tid)
                        if int(t['is_ensembl_canonical']):canonical[bid].add(s['variant_id'])
        for b in linked:
            h1,h2=L.orient(int(b['ase_a_count']),int(b['ase_b_count']),b['wgs_phase_orientation'])
            r=base(c,b);r.update({f'nominal_ge_{d}':c[f'nominal_ge_{d}'] for d in (3,5,8,10)});r.update(ranking_fields(float(c['delta_M']),h1,h2));bid=b['haplotype_block_id']
            r.update(same_transcript_RNA_support=int(bool(support[bid])),same_canonical_transcript_RNA_support=int(bool(canonical[bid])),
                n_same_transcript_exonic_SNPs=len(support[bid]),same_transcript_exonic_SNPs=','.join(sorted(support[bid])),
                n_same_canonical_exonic_SNPs=len(canonical[bid]),compatible_transcripts=','.join(sorted(transcripts[bid])))
            blocks_out.append(r)
        audit.append(dict(sample=key[0],tissue=key[1],gene_id=key[2],gene_name=c['annotation_gene_name'],chrom=c['chrom'],ps=c['ps'],cpg_pos1=c['cpg_pos1'],
            n_source_gene_blocks=len(relevant),n_resolved_same_PS_blocks=len(linked),n_same_PS_same_transcript_blocks=sum(bool(support[b['haplotype_block_id']]) for b in linked),
            observed_RNA_phase_sets=','.join(sorted({b['wgs_phase_set'] for b in relevant if b['wgs_phase_set']})),
            n_local_phase_PASS_same_transcript_blocks=sum(bool(support[b['haplotype_block_id']]) and b['local_phase_status']=='PASS' for b in linked),
            status='NO_MATCHED_RNA_SAMPLE' if (c['sample'],c['tissue'])==('0309','Liver') else 'SAME_PS_BLOCK' if linked else 'NO_SOURCE_BLOCK' if not relevant else 'NO_RESOLVED_SAME_PS_BLOCK',
            **{f'nominal_ge_{d}':c[f'nominal_ge_{d}'] for d in (3,5,8,10)}))
    return blocks_out,audit,attempts
