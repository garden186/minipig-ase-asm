import csv,gzip
import numpy as np
STATE=['ALIGNED','NO_ELIGIBLE_SITE_ANCHOR','NO_MATCHING_VARIANT_RECORD','DIFFERENT_EXACT_PS',
       'DUPLICATE_PHASED_POSITION','NOT_HETEROZYGOUS','NOT_PHASED','MISSING_PS','REFERENCE_MISMATCH',
       'DUPLICATE_SNPSPLIT_POSITION','SNPSPLIT_MISSING','SNPSPLIT_ORDER_MISMATCH','DUPLICATE_GENOTYPE_POSITION',
       'GENOTYPE_DISCORDANT','GQ_MISSING','GQ_BELOW_20']
PREP_STATUS=['ELIGIBLE_ANCHOR','DUPLICATE_PHASED_POSITION','NOT_HETEROZYGOUS','NOT_PHASED','MISSING_PS',
       'REFERENCE_MISMATCH','DUPLICATE_SNPSPLIT_POSITION','SNPSPLIT_MISSING','SNPSPLIT_ORDER_MISMATCH',
       'DUPLICATE_GENOTYPE_POSITION','GENOTYPE_DISCORDANT','GQ_MISSING','GQ_BELOW_20']
GT=['.','0|1','1|0','0/1','1/0','0/0','1/1','OTHER','0|0','1|1']
PROFILE=np.dtype([('pos','<u4'),('sample','u1'),('tissue','u1'),('ps','<u4'),('source_row','<u4'),
    ('a','<u4'),('b','<u4'),('c','<u4'),('d','<u4'),('original_p','<f8'),('original_q','<f8'),
    ('state','u1'),('GT','u1'),('anchor_ps','<u4'),('alt_hap','u1'),('logp','<f8',(2,))])
SITE=np.dtype([('pos','<u4'),('anchor','<u8'),('anchor_testable','u1'),('anchor_available','u1'),('candidate_snps','<u4'),
    ('N_before','u1',(3,)),('N','u1',(3,))])
HIT=np.dtype([('chrom','u1'),('pos','<u4'),('N','u1'),('logp','<f8'),('logq_bh','<f8'),('logq','<f8')])

def rows(p):
    with gzip.open(p,'rt',encoding='utf-8',newline='') as f:yield from csv.DictReader(f,delimiter='\t')

def table(p,data,fields):
    with gzip.open(p,'wt',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
