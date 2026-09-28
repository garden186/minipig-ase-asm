"""Portable GATK ASEReadCounter launch with original fragment-counting options."""
import argparse,csv,json,shlex,subprocess,sys
from pathlib import Path
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True,type=Path);p.add_argument('--reference',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--gatk',default='gatk');p.add_argument('--execute',action='store_true');a=p.parse_args()
    with a.manifest.open(encoding='utf-8-sig',newline='') as h:rows=list(csv.DictReader(h,delimiter='\t'))
    if not rows:raise ValueError('Empty manifest')
    seen=set();jobs=[];manifest=[]
    for r in rows:
        for k in ['sample','tissue','bam','vcf']:
            if not r.get(k):raise ValueError('Missing '+k)
        import re
        if not all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',r[k]) for k in ['sample','tissue']):raise ValueError('Invalid sample/tissue identifier')
        name=r['sample']+'_'+r['tissue']
        if name in seen:raise ValueError('Duplicate unit')
        seen.add(name);tmp=a.output/'tmp'/name
        argv=[a.gatk,'--java-options',f'-Xms1g -Xmx4g -XX:+UseSerialGC -XX:ActiveProcessorCount=1 -Djava.io.tmpdir={tmp}',
            'ASEReadCounter','-R',str(a.reference),'-I',r['bam'],'-V',r['vcf']]
        for c in range(1,19):argv+=['--intervals',str(c)]
        argv+=['-O',str(a.output/'counts'/(name+'.ase_readcounter.tsv')),'--tmp-dir',str(tmp),'--output-format','TABLE',
            '--min-base-quality','10','--min-mapping-quality','255','--min-depth-of-non-filtered-base','1','--max-depth-per-sample','0',
            '--count-overlap-reads-handling','COUNT_FRAGMENTS_REQUIRE_SAME_BASE','--read-filter','ProperlyPairedReadFilter',
            '--read-filter','NotSupplementaryAlignmentReadFilter','--read-filter','FragmentLengthReadFilter','--max-fragment-length','1000000']
        jobs.append((name,tmp,argv));manifest.append(dict(sample=r['sample'],tissue_raw=r.get('tissue_raw',r['tissue']),tissue=r['tissue'],bam=r['bam'],vcf=r['vcf'],job_id=name))
        print(shlex.join(argv))
    if not a.execute:return
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'counts').mkdir();(a.output/'logs').mkdir()
    with (a.output/'input_manifest.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(manifest[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(manifest)
    (a.output/'commands.json').write_text(json.dumps([x[2] for x in jobs],indent=2)+'\n')
    for name,tmp,argv in jobs:
        tmp.mkdir(parents=True)
        with (a.output/'logs'/(name+'.log')).open('w') as h:subprocess.run(argv,stdout=h,stderr=subprocess.STDOUT,check=True)
    (a.output/'COMPLETE').write_text('All GATK commands exited successfully; validate downstream count tables.\n')
if __name__=='__main__':main()
