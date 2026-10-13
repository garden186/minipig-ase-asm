#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate N-masked reference FASTA by replacing heterozygous SNP positions
with 'N' for SNPsplit-compatible alignment.

Input:
- Reference FASTA (e.g., Sus_scrofa.Sscrofa11.1.dna.toplevel.fa)
- Phased VCF (or any VCF with het sites)

Output:
- N-masked FASTA: same as reference, but with 'N' at every het SNP position

Logic:
- For each VCF record passing filters (PASS, biallelic, SNP, het GT 0|1 or 1|0):
  set FASTA[chrom][pos-1] = 'N'
- Indels are NOT masked (would require length changes)
- Position is 1-based in VCF, converted to 0-based for sequence indexing

Notes on case preservation:
- Reference FASTA may have soft-masked (lowercase) regions for repeats
- We preserve case where possible: 'N' replaces both upper and lowercase bases
- Output uses 'N' (uppercase) consistently at het sites

Usage:
    python3 generate_nmasked_reference.py \\
        --fasta /path/to/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \\
        --vcf  /path/to/0349.phased.wgs.honest.vcf.gz \\
        --out  /path/to/Sus_scrofa.Sscrofa11.1.Nmasked.fa \\
        [--chr 18]   # optional: process single chromosome for testing

Verification:
- Output has same length per chromosome as input (no insertions/deletions)
- N count in output = N count in input + heterozygous SNPs masked
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict


def open_vcf(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--fasta', required=True, help='Reference FASTA (uncompressed)')
    p.add_argument('--vcf',   required=True, help='Phased VCF (gzipped OK)')
    p.add_argument('--out',   required=True, help='Output N-masked FASTA')
    p.add_argument('--chr',   default=None, help='Process only this chromosome (for testing)')
    return p.parse_args()


def read_fasta(path, chr_filter=None):
    """
    Read FASTA into {chrom: bytearray(sequence)}.
    Use bytearray for fast in-place mutation.
    """
    sys.stderr.write(f'[INFO] Reading FASTA: {path}\n')
    seqs = {}
    headers = {}  # preserve full header line for output
    cur_chrom = None
    cur_chunks = []

    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Save previous
                if cur_chrom is not None:
                    if chr_filter is None or cur_chrom == chr_filter:
                        seqs[cur_chrom] = bytearray(''.join(cur_chunks), 'ascii')
                # New chrom
                full_header = line[1:].rstrip()
                cur_chrom = full_header.split()[0]  # first token is chrom name
                headers[cur_chrom] = full_header
                cur_chunks = []
            else:
                cur_chunks.append(line.rstrip())
        # Save last
        if cur_chrom is not None:
            if chr_filter is None or cur_chrom == chr_filter:
                seqs[cur_chrom] = bytearray(''.join(cur_chunks), 'ascii')

    sys.stderr.write(f'[INFO] Loaded {len(seqs)} chromosomes\n')
    for chrom, seq in sorted(seqs.items()):
        sys.stderr.write(f'         {chrom}: {len(seq):,} bp\n')

    return seqs, headers


def mask_heterozygous_sites(seqs, vcf_path, chr_filter=None):
    """
    Walk through VCF and set 'N' at every het SNP position.
    Returns: dict of per-chromosome mask counts.
    """
    sys.stderr.write(f'[INFO] Reading VCF: {vcf_path}\n')

    mask_counts = defaultdict(int)
    skipped = defaultdict(int)
    sample_idx = None

    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                sample_idx = 9
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 10:
                continue

            chrom, pos, _, ref, alt, _, filt, _, fmt = fields[:9]
            sample = fields[sample_idx]

            # Optional chromosome filter
            if chr_filter is not None and chrom != chr_filter:
                continue

            # Skip if chromosome not in FASTA
            if chrom not in seqs:
                skipped['chrom_not_in_fasta'] += 1
                continue

            # Filter
            if filt not in ('PASS', '.'):
                skipped['filter_fail'] += 1
                continue

            # Biallelic
            if ',' in alt:
                skipped['multiallelic'] += 1
                continue

            # SNP only (skip indels — masking would shift positions)
            if len(ref) != 1 or len(alt) != 1:
                skipped['indel'] += 1
                continue

            # Parse GT
            fmt_keys = fmt.split(':')
            sample_vals = sample.split(':')
            if 'GT' not in fmt_keys:
                skipped['no_gt'] += 1
                continue
            gt = sample_vals[fmt_keys.index('GT')]

            # Heterozygous?
            sep = '|' if '|' in gt else ('/' if '/' in gt else None)
            if sep is None:
                skipped['unparseable_gt'] += 1
                continue
            alleles = gt.split(sep)
            if len(alleles) != 2 or alleles[0] == alleles[1]:
                skipped['homozygous'] += 1
                continue
            if set(alleles) != {'0', '1'}:
                skipped['unexpected_gt'] += 1
                continue

            # Mask: convert pos to 0-based index
            pos0 = int(pos) - 1
            if pos0 < 0 or pos0 >= len(seqs[chrom]):
                skipped['out_of_bounds'] += 1
                continue

            # In-place mask
            seqs[chrom][pos0] = ord('N')
            mask_counts[chrom] += 1

    return mask_counts, skipped


def write_fasta(seqs, headers, out_path, line_width=60):
    """Write FASTA with line wrapping (default 60 chars)."""
    sys.stderr.write(f'[INFO] Writing N-masked FASTA: {out_path}\n')
    with open(out_path, 'w') as f:
        for chrom in sorted(seqs.keys()):
            f.write(f'>{headers.get(chrom, chrom)}\n')
            seq = seqs[chrom].decode('ascii')
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width])
                f.write('\n')


def main():
    args = parse_args()

    # 1) Load FASTA
    seqs, headers = read_fasta(args.fasta, chr_filter=args.chr)
    if not seqs:
        sys.exit('[ERROR] No sequences loaded from FASTA')

    # 2) Count Ns before masking (for verification)
    n_before = {chrom: seq.count(ord('N')) + seq.count(ord('n'))
                for chrom, seq in seqs.items()}

    # 3) Mask heterozygous sites
    mask_counts, skipped = mask_heterozygous_sites(seqs, args.vcf, chr_filter=args.chr)

    # 4) Count Ns after masking
    n_after = {chrom: seq.count(ord('N')) + seq.count(ord('n'))
               for chrom, seq in seqs.items()}

    # 5) Write output
    write_fasta(seqs, headers, args.out)

    # 6) Report
    sys.stderr.write('\n=== N-masked reference generation complete ===\n')
    sys.stderr.write(f'Input FASTA: {args.fasta}\n')
    sys.stderr.write(f'Input VCF:   {args.vcf}\n')
    sys.stderr.write(f'Output:      {args.out}\n')
    if args.chr:
        sys.stderr.write(f'Chr filter:  {args.chr}\n')

    sys.stderr.write('\nPer-chromosome masking:\n')
    sys.stderr.write(f'  {"chrom":<10} {"length":>15} {"N before":>15} {"masked":>15} {"N after":>15}\n')
    total_masked = 0
    for chrom in sorted(seqs.keys()):
        masked = mask_counts.get(chrom, 0)
        total_masked += masked
        sys.stderr.write(
            f'  {chrom:<10} {len(seqs[chrom]):>15,} '
            f'{n_before[chrom]:>15,} {masked:>15,} {n_after[chrom]:>15,}\n'
        )
    sys.stderr.write(f'  {"TOTAL":<10} '
                     f'{sum(len(s) for s in seqs.values()):>15,} '
                     f'{sum(n_before.values()):>15,} '
                     f'{total_masked:>15,} '
                     f'{sum(n_after.values()):>15,}\n')

    # Verification: N_after - N_before should equal total_masked
    delta_n = sum(n_after.values()) - sum(n_before.values())
    sys.stderr.write(f'\nVerification: N delta = {delta_n:,} (should equal masked = {total_masked:,})\n')
    if delta_n == total_masked:
        sys.stderr.write('  [OK] N count matches masking count\n')
    else:
        sys.stderr.write(f'  [WARN] Mismatch by {total_masked - delta_n:,} '
                         '(may indicate ref base was already N at some positions)\n')

    if skipped:
        sys.stderr.write('\nSkipped variants:\n')
        for reason, count in sorted(skipped.items()):
            sys.stderr.write(f'  {reason:<25} {count:>12,}\n')


if __name__ == '__main__':
    main()
