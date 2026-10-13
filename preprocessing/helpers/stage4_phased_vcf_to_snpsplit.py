#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Convert phased VCF genotypes to SNPsplit allele order:
#   left allele  (before '|') -> genome1
#   right allele (after  '|') -> genome2
#
# genome1/genome2 are local haplotype labels. Their physical homolog
# identity is consistent only within the same CHROM + PS phase block.
# Do not merge G1/G2 labels across different phase sets or interpret
# them as paternal/maternal without parent-of-origin phasing.

Critical logic (verified in CpelAsm.jl line 372):
- VCF GT='0|1': REF=hap1, ALT=hap2 → SNP file: Ref=REF, SNP=ALT (no swap)
- VCF GT='1|0': REF=hap2, ALT=hap1 → SNP file: Ref=ALT, SNP=REF (SWAP)

After swap:
- 'Ref' column in SNP file = always hap1 allele
- 'SNP' column in SNP file = always hap2 allele
- SNPsplit G1 reads = hap1 reads (local to the same chromosome and phase set)
- SNPsplit G2 reads = hap2 reads (local to the same chromosome and phase set)

Output format (SNPsplit standard, tab-separated):
    SNP-ID  Chromosome  Position  Strand  Ref/SNP

Excluded variants:
- Indels (SNPsplit only handles SNPs reliably)
- Multi-allelic (>1 ALT)
- Unphased (GT='0/1' or no PS tag)
- Non-PASS filter
- Both alleles same length but >1 (MNV) — exclude for safety

Usage:
    python3 phased_vcf_to_snpsplit.py \\
        --vcf 0349.phased.wgs.honest.vcf.gz \\
        --out 0349.snpsplit_snps.txt \\
        [--chr 18]   # optional: process single chromosome for testing
"""

import argparse
import gzip
import sys
from collections import Counter


def open_vcf(path):
    """Open VCF (gzipped or plain)."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--vcf', required=True, help='Phased VCF (gzipped OK)')
    p.add_argument('--out', required=True, help='Output SNP file (SNPsplit format)')
    p.add_argument('--chr', default=None, help='Process only this chromosome (for testing)')
    p.add_argument('--require-ps', action='store_true', default=True,
                   help='Require PS tag (default: True)')
    p.add_argument('--include-indels', action='store_true',
                   help='Include indels (default: SNPs only)')
    return p.parse_args()


def main():
    args = parse_args()

    counts = Counter()
    sample_idx = None  # column index of sample in genotype field

    with open_vcf(args.vcf) as fin, open(args.out, 'w') as fout:
        for line in fin:
            line = line.rstrip('\n')

            # ----- Headers -----
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                # Sample column is index 9 (0-based) for first sample
                fields = line.split('\t')
                if len(fields) < 10:
                    sys.exit('[ERROR] No sample column in VCF')
                sample_idx = 9
                continue

            # ----- Records -----
            fields = line.split('\t')
            if len(fields) < 10:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info, fmt = fields[:9]
            sample = fields[sample_idx]

            # Optional chromosome filter
            if args.chr is not None and chrom != args.chr:
                continue

            counts['total'] += 1

            # FILTER must be PASS or '.'
            if filt not in ('PASS', '.'):
                counts['filter_fail'] += 1
                continue

            # Multi-allelic check
            if ',' in alt:
                counts['multiallelic'] += 1
                continue

            # SNP-only check (skip indels)
            if not args.include_indels:
                if len(ref) != 1 or len(alt) != 1:
                    counts['indel_or_mnv'] += 1
                    continue

            # Parse FORMAT and sample fields
            fmt_keys = fmt.split(':')
            sample_vals = sample.split(':')
            if 'GT' not in fmt_keys:
                counts['no_gt'] += 1
                continue

            gt_idx = fmt_keys.index('GT')
            gt = sample_vals[gt_idx]

            # Phasing check
            if '|' not in gt:
                counts['unphased'] += 1
                continue

            # PS tag check (WhatsHap phase set)
            if args.require_ps and 'PS' not in fmt_keys:
                counts['no_ps'] += 1
                continue

            # Heterozygous check
            alleles = gt.split('|')
            if len(alleles) != 2 or alleles[0] == alleles[1]:
                counts['homozygous'] += 1
                continue
            if set(alleles) != {'0', '1'}:
                # e.g., 0|2 or 1|2 — biallelic only after multi-allelic filter
                # so this should not happen, but just in case
                counts['unexpected_gt'] += 1
                continue

            # ===== KEY LOGIC: phasing-aware swap =====
            # gt='0|1' → REF=hap1, ALT=hap2 → write 'REF/ALT' (no swap)
            # gt='1|0' → REF=hap2, ALT=hap1 → write 'ALT/REF' (SWAP)
            if gt == '0|1':
                ref_out, snp_out = ref, alt
                counts['gt_0_1'] += 1
            elif gt == '1|0':
                ref_out, snp_out = alt, ref  # swap!
                counts['gt_1_0_swapped'] += 1
            else:
                # shouldn't reach here after checks above
                counts['unexpected_gt'] += 1
                continue

            # SNP-ID: use VCF ID if present, else chr_pos
            if vid == '.' or vid == '':
                snp_id = f'{chrom}_{pos}'
            else:
                snp_id = vid

            # Write SNPsplit format
            # SNP-ID  Chromosome  Position  Strand  Ref/SNP
            # Strand always 1 (SNPsplit ignores it per user guide)
            fout.write(f'{snp_id}\t{chrom}\t{pos}\t1\t{ref_out}/{snp_out}\n')
            counts['written'] += 1

    # ----- Report -----
    sys.stderr.write('\n=== Phased VCF → SNPsplit SNP file ===\n')
    sys.stderr.write(f'Input:  {args.vcf}\n')
    sys.stderr.write(f'Output: {args.out}\n')
    if args.chr:
        sys.stderr.write(f'Chromosome filter: {args.chr}\n')
    sys.stderr.write(f'\nVariant statistics:\n')
    sys.stderr.write(f'  Total records:        {counts["total"]:>12,}\n')
    sys.stderr.write(f'  Written (SNPs only):  {counts["written"]:>12,}\n')
    sys.stderr.write(f'    GT=0|1 (no swap):   {counts["gt_0_1"]:>12,}\n')
    sys.stderr.write(f'    GT=1|0 (SWAPPED):   {counts["gt_1_0_swapped"]:>12,}\n')
    sys.stderr.write(f'\nExcluded:\n')
    sys.stderr.write(f'  Filter fail:          {counts["filter_fail"]:>12,}\n')
    sys.stderr.write(f'  Multi-allelic:        {counts["multiallelic"]:>12,}\n')
    sys.stderr.write(f'  Indel/MNV:            {counts["indel_or_mnv"]:>12,}\n')
    sys.stderr.write(f'  No GT field:          {counts["no_gt"]:>12,}\n')
    sys.stderr.write(f'  Unphased:             {counts["unphased"]:>12,}\n')
    sys.stderr.write(f'  No PS tag:            {counts["no_ps"]:>12,}\n')
    sys.stderr.write(f'  Homozygous:           {counts["homozygous"]:>12,}\n')
    sys.stderr.write(f'  Unexpected GT:        {counts["unexpected_gt"]:>12,}\n')

    # Sanity check
    if counts['written'] > 0:
        ratio_1_0 = counts['gt_1_0_swapped'] / counts['written']
        sys.stderr.write(f'\n1|0 fraction: {ratio_1_0:.3f}')
        if 0.40 < ratio_1_0 < 0.60:
            sys.stderr.write('  [OK: balanced as expected]\n')
        else:
            sys.stderr.write('  [WARN: unusual ratio — investigate]\n')


if __name__ == '__main__':
    main()
