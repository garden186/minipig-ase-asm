#!/usr/bin/env python3
"""Portable execution of the study's final NGS preprocessing commands (Linux).

The default prints commands. --execute runs them in a new stage directory.
FASTQ/reference data are external; sample identifiers never imply raw filenames.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def command(*args):
    return shlex.join([str(x) for x in args])


def read_samples(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    seen = set()
    for r in rows:
        for name in ('sample', 'assay', 'tissue', 'read1', 'read2'):
            if not r.get(name): raise ValueError('Missing sample-sheet field: ' + name)
        if r['assay'] not in ('WGS', 'RNA', 'WGBS'): raise ValueError('Assay must be WGS, RNA, or WGBS')
        if not all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', r[k]) for k in ('sample', 'tissue')):
            raise ValueError('Use simple sample/tissue identifiers; paths may contain spaces')
        key = r['sample'], r['assay'], r['tissue']
        if key in seen: raise ValueError('Duplicate sample/assay/tissue')
        seen.add(key)
    return rows


def plan(config, rows, stage, sample):
    work = Path(config['work_dir']).expanduser().resolve()
    fasta = Path(config['reference_fasta']).expanduser().resolve()
    gtf = Path(config['annotation_gtf']).expanduser().resolve()
    ref = work/'reference'/'reference.fa'
    star = work/'reference'/'star'
    tools = config.get('tools', {})
    tool = lambda name: tools.get(name, name)
    threads = int(config.get('threads', 16))
    if threads < 1: raise ValueError('threads must be positive')
    steps = []
    def add(name, cmd): steps.append({'name': name, 'command': cmd})
    def mkdir(p): add('mkdir', command('mkdir', '-p', p))
    wgs = work/'wgs'/str(sample)
    vcf = wgs/'phased.vcf.gz'
    if stage == 'reference':
        mkdir(ref.parent); mkdir(star)
        # Keep input reference immutable; generated indices belong to work_dir.
        add('copy-reference', command('cp', fasta, ref))
        add('fasta-index', command(tool('samtools'), 'faidx', ref))
        add('sequence-dictionary', command(tool('gatk'), 'CreateSequenceDictionary', '-R', ref, '-O', ref.with_suffix('.dict')))
        add('bwa-index', command(tool('bwa'), 'index', ref))
        add('STAR-index', command(tool('STAR'), '--runMode', 'genomeGenerate', '--runThreadN', threads,
            '--genomeDir', star, '--genomeFastaFiles', ref, '--sjdbGTFfile', gtf, '--sjdbOverhang', 150))
        return steps
    if not sample: raise ValueError('--sample is required for this stage')
    selected = [r for r in rows if r['sample'] == sample]
    if not selected: raise ValueError('Sample is absent from sample sheet')
    if stage == 'wgs':
        raw = [r for r in selected if r['assay'] == 'WGS']
        if len(raw) != 1: raise ValueError('Exactly one WGS library per sample is required')
        r = raw[0]; mkdir(wgs); mkdir(wgs/'tmp')
        r1, r2 = wgs/'trim_R1.fq.gz', wgs/'trim_R2.fq.gz'
        add('fastp', command(tool('fastp'), '--in1', r['read1'], '--in2', r['read2'], '--out1', r1, '--out2', r2,
            '--adapter_sequence', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', '--adapter_sequence_r2', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
            '--cut_right', '--cut_right_window_size', 4, '--cut_right_mean_quality', 15, '--qualified_quality_phred', 15,
            '--unqualified_percent_limit', 40, '--n_base_limit', 5, '--length_required', 50, '--trim_poly_g', '--poly_g_min_len', 10,
            '--thread', 8, '--json', wgs/'fastp.json', '--html', wgs/'fastp.html'))
        add('BWA-sort', command(tool('bwa'), 'mem', '-t', threads, '-K', 100000000, '-Y', '-R',
            f'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:{sample}', ref, r1, r2) + ' | ' +
            command(tool('samtools'), 'sort', '-@', 8, '-m', '4G', '-T', wgs/'tmp/sort', '-o', wgs/'aligned.bam', '-'))
        # The original WGS code inferred this from the instrument prefix.
        # The explicit sample-sheet value avoids guessing on relocated inputs.
        pixel_value = r.get('optical_duplicate_distance') or config.get('wgs_optical_duplicate_distance')
        if pixel_value is None or pixel_value == '':
            raise ValueError('Set optical_duplicate_distance for this WGS library (historical instrument-dependent value: 100 or 2500)')
        pixel = int(pixel_value)
        if pixel <= 0: raise ValueError('optical_duplicate_distance must be positive')
        add('mark-duplicates', command(tool('gatk'), '--java-options', f'-Xmx16g -Djava.io.tmpdir={wgs / "tmp"}',
            'MarkDuplicates', '-I', wgs/'aligned.bam', '-O', wgs/'marked.bam', '-M', wgs/'duplicates.txt',
            '--REMOVE_DUPLICATES', 'false', '--VALIDATION_STRINGENCY', 'LENIENT', '--ASSUME_SORT_ORDER', 'coordinate',
            '--OPTICAL_DUPLICATE_PIXEL_DISTANCE', pixel, '--TAGGING_POLICY', 'OpticalOnly'))
        add('marked-index', command(tool('samtools'), 'index', wgs/'marked.bam'))
        add('coverage-QC', command(tool('mosdepth'), '-t', 4, '--by', 1000, '--no-per-base', '--fast-mode',
            '--mapq', 20, '--flag', 3844, wgs/'coverage', wgs/'marked.bam'))
        add('flagstat', command(tool('samtools'), 'flagstat', '-@', 4, wgs/'marked.bam')+' > '+shlex.quote(str(wgs/'flagstat.txt')))
        add('phase-ready', command(tool('samtools'), 'view', '-@', threads, '-b', '-f', 2, '-F', 3844, '-q', 30,
            '-o', wgs/'phase_ready.bam', wgs/'marked.bam'))
        add('phase-ready-index', command(tool('samtools'), 'index', wgs/'phase_ready.bam'))
        mount = str(work)
        image = config.get('deepvariant_image', 'google/deepvariant:1.10.0')
        add('DeepVariant', command(tool('docker'), 'run', '--rm', '-v', mount+':/work', image,
            '/opt/deepvariant/bin/run_deepvariant', '--model_type=WGS', '--ref=/work/reference/reference.fa',
            f'--reads=/work/wgs/{sample}/marked.bam', f'--output_vcf=/work/wgs/{sample}/variants.vcf.gz',
            f'--output_gvcf=/work/wgs/{sample}/variants.g.vcf.gz', '--num_shards='+str(threads), '--sample_name='+sample,
            '--regions='+' '.join(map(str, range(1,19))), f'--intermediate_results_dir=/work/wgs/{sample}/dv_intermediate'))
        bcf = tool('bcftools')
        add('het-SNVs', command(bcf, 'view', '-r', ','.join(map(str,range(1,19))), wgs/'variants.vcf.gz', '-Ou')+' | '+
            command(bcf, 'view', '-m2', '-M2', '-Ou')+' | '+command(bcf, 'norm', '-f', ref, '-c', 'w', '-Ou')+' | '+
            command(bcf, 'view', '-v', 'snps', '-f', 'PASS', '-Ou')+' | '+
            command(bcf, 'view', '-g', 'het', '-i', 'FORMAT/DP>=10', '-Oz', '-o', wgs/'heterozygous.vcf.gz'))
        add('het-index', command(bcf, 'index', '-t', wgs/'heterozygous.vcf.gz'))
        add('WhatsHap', command(tool('whatshap'), 'phase', '-o', vcf, '--reference', ref, '--tag', 'PS', '--only-snvs',
            '--sample', sample, '--mapping-quality', 30, '--output-read-list', wgs/'phased_reads.txt', wgs/'heterozygous.vcf.gz', wgs/'phase_ready.bam'))
        add('phase-index', command(bcf, 'index', '-t', vcf))
        add('phase-QC', command(tool('whatshap'), 'stats', '--tsv', wgs/'phase_stats.tsv', '--block-list', wgs/'phase_blocks.tsv', vcf))
    elif stage == 'rna':
        raw = sorted([r for r in selected if r['assay']=='RNA'], key=lambda r:r['tissue'])
        if not raw: raise ValueError('No RNA libraries')
        het = work/'rna'/sample/'heterozygous.vcf'; mkdir(het.parent)
        add('uncompress-VCF', command(tool('bcftools'), 'view', '-Ov', '-o', het, wgs/'heterozygous.vcf.gz'))
        junctions=[]
        for r in raw:
            d=work/'rna'/sample/r['tissue']; mkdir(d); mkdir(d/'pass1')
            add('RNA-trim', command(tool('trim_galore'), '--paired', '--quality', 20, '--length', 20, '--cores', 4,
                '--basename', 'reads', '--output_dir', d, r['read1'], r['read2']))
            base=[tool('STAR'),'--runMode','alignReads','--runThreadN',threads,'--genomeDir',star,
                '--readFilesIn',d/'reads_val_1.fq.gz',d/'reads_val_2.fq.gz','--readFilesCommand','zcat',
                '--outFilterType','BySJout','--outFilterMultimapNmax',20,'--alignSJoverhangMin',8,
                '--alignSJDBoverhangMin',1,'--alignIntronMin',20,'--alignIntronMax',1000000,'--alignMatesGapMax',1000000]
            add('STAR-pass1', command(*base,'--outSAMtype','None','--outFileNamePrefix',str(d/'pass1')+'/'))
            junctions.append(d/'pass1/SJ.out.tab')
        for r in raw:
            d=work/'rna'/sample/r['tissue']; mkdir(d/'tmp')
            base=[tool('STAR'),'--runMode','alignReads','--runThreadN',threads,'--genomeDir',star,
                '--readFilesIn',d/'reads_val_1.fq.gz',d/'reads_val_2.fq.gz','--readFilesCommand','zcat',
                '--outFilterType','BySJout','--outFilterMultimapNmax',20,'--alignSJoverhangMin',8,
                '--alignSJDBoverhangMin',1,'--alignIntronMin',20,'--alignIntronMax',1000000,'--alignMatesGapMax',1000000]
            add('STAR-pass2-WASP', command(*base,'--sjdbFileChrStartEnd',*junctions,'--sjdbOverhang',150,
                '--outSAMtype','BAM','SortedByCoordinate','--outSAMattributes','NH','HI','AS','nM','NM','MD','vA','vG','vW',
                '--waspOutputMode','SAMtag','--varVCFfile',het,'--outFilterMismatchNmax',999,'--outFilterMismatchNoverReadLmax',0.04,
                '--outSAMattrRGline','ID:'+sample+'_'+r['tissue'],'SM:'+sample,'PL:ILLUMINA','LB:'+sample+'_'+r['tissue'],
                '--outFileNamePrefix',str(d)+'/'))
            add('RNA-mark-duplicates', command(tool('gatk'),'--java-options',f'-Xmx16g -Djava.io.tmpdir={d / "tmp"}',
                'MarkDuplicates','-I',d/'Aligned.sortedByCoord.out.bam','-O',d/'marked.bam','-M',d/'duplicates.txt',
                '--REMOVE_DUPLICATES','false','--VALIDATION_STRINGENCY','LENIENT','--ASSUME_SORT_ORDER','coordinate',
                '--OPTICAL_DUPLICATE_PIXEL_DISTANCE',2500,'--TAGGING_POLICY','OpticalOnly'))
            add('RNA-filter-WASP', command(tool('samtools'),'view','-h',d/'marked.bam')+' | '+
                command('awk','/^@/ || !/vW:i:[2-7]/')+' | '+command(tool('samtools'),'view','-bh','-F',2304,'-@',4,'-')+
                ' > '+shlex.quote(str(d/'phase_ready.bam')))
            add('RNA-index',command(tool('samtools'),'index',d/'phase_ready.bam'))
    elif stage == 'phaser':
        tissues=sorted(r['tissue'] for r in selected if r['assay']=='RNA')
        if not tissues: raise ValueError('No RNA libraries')
        script=config['phaser_script']; d=work/'phaser'/sample; mkdir(d);mkdir(d/'tmp_rna');mkdir(d/'tmp_wgs')
        bams=[wgs/'phase_ready.bam']+[work/'rna'/sample/t/'phase_ready.bam' for t in tissues]
        for kind,b,exclude in [('rna',bams,True),('wgs',bams[:1],False)]:
            argv=[tool('python'),script,'--vcf',vcf,'--bam',','.join(map(str,b)),'--sample',sample,
                '--mapq',','.join(['30']+(['255']*len(tissues) if exclude else [])),
                '--isize',','.join(['1000']+(['1000000']*len(tissues) if exclude else [])),
                '--paired_end',','.join(['1']*len(b)),'--baseq',10,'--threads',int(config.get('phaser_threads',2)),
                '--cc_threshold',0.05,'--gw_phase_method',0,'--gw_phase_vcf',1,'--unphased_vars',1,'--pass_only',1,
                '--o',d/kind,'--temp_dir',d/('tmp_'+kind)]
            if exclude: argv+=['--haplo_count_bam_exclude',1]
            add('phASER-'+kind,command(*argv))
    elif stage == 'wgbs':
        raw=[r for r in selected if r['assay']=='WGBS']
        if not raw: raise ValueError('No WGBS libraries')
        d=work/'wgbs'/sample; genome=d/'genome';mkdir(genome)
        add('allele-list',command(tool('python'),HERE/'helpers/stage4_phased_vcf_to_snpsplit.py','--vcf',vcf,'--out',d/'alleles.tsv'))
        add('N-mask',command(tool('python'),HERE/'helpers/stage4_generate_nmasked_reference.py','--fasta',ref,'--vcf',vcf,'--out',genome/'masked.fa'))
        add('Bismark-index',command(tool('bismark_genome_preparation'),'--bowtie2','--parallel',4,genome))
        for r in raw:
            q=d/r['tissue'];mkdir(q)
            add('WGBS-trim',command(tool('trim_galore'),'--paired','--cores',4,'--clip_R1',10,'--clip_R2',15,
                '--three_prime_clip_R1',10,'--three_prime_clip_R2',10,'--quality',20,'--basename','reads','--output_dir',q,r['read1'],r['read2']))
            add('Bismark',command(tool('bismark'),'--genome',genome,'--bowtie2','--score_min','L,0,-0.6','-X',700,'--dovetail',
                '--parallel',int(config.get('bismark_parallel',4)),'--output_dir',q,'-1',q/'reads_val_1.fq.gz','-2',q/'reads_val_2.fq.gz'))
            bam=q/'reads_val_1_bismark_bt2_pe.bam';dedup=q/'reads_val_1_bismark_bt2_pe.deduplicated.bam'
            add('Bismark-dedup',command(tool('deduplicate_bismark'),'--paired','--bam','--output_dir',q,bam))
            add('coordinate-sort',command(tool('samtools'),'sort','-@',threads,'-o',q/'asm_input.bam',dedup))
            add('ASM-BAM-index',command(tool('samtools'),'index',q/'asm_input.bam'))
    else: raise ValueError('Unknown stage')
    return steps


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['reference','wgs','rna','phaser','wgbs'])
    p.add_argument('--config',required=True,type=Path)
    p.add_argument('--samples',required=True,type=Path)
    p.add_argument('--sample')
    p.add_argument('--execute',action='store_true')
    args=p.parse_args();cfg=json.loads(args.config.read_text(encoding='utf-8-sig'))
    # Config paths resolve beside the config; FASTQ paths beside the sample sheet.
    for key in ['work_dir','reference_fasta','annotation_gtf','phaser_script']:
        if key in cfg:
            value=Path(cfg[key]).expanduser()
            cfg[key]=str((args.config.parent/value if not value.is_absolute() else value).resolve())
    rows=read_samples(args.samples)
    for row in rows:
        for key in ['read1','read2']:
            value=Path(row[key]).expanduser()
            row[key]=str((args.samples.parent/value if not value.is_absolute() else value).resolve())
    steps=plan(cfg,rows,args.stage,args.sample)
    for s in steps: print(s['command'])
    if not args.execute:return
    if sys.platform!='linux':raise SystemExit('Execution requires Linux; command review works on other systems.')
    work=Path(cfg['work_dir']).expanduser().resolve();log=work/'run_records'/(args.stage+'_'+str(args.sample or 'all'))
    log.mkdir(parents=True,exist_ok=False)
    (log/'plan.json').write_text(json.dumps({'config':cfg,'samples':rows,'steps':steps},indent=2)+'\n')
    for i,s in enumerate(steps):
        with (log/f'{i:03d}_{s["name"]}.log').open('w') as h:
            subprocess.run(['bash','-e','-o','pipefail','-c',s['command']],stdout=h,stderr=subprocess.STDOUT,check=True)
    (log/'COMPLETE').write_text('All commands exited successfully.\n')


if __name__=='__main__':main()
