"""Build a reusable interval/profile index with frozen historical annotation rules."""
from bisect import bisect_right
from collections import defaultdict
from functools import lru_cache
import json
from pathlib import Path
import sqlite3
import sys
import time

from common import ROOT,VERSION,POLICY,GTF_SHA256,Guard,digest,sha,read_json,write_json,writer,rows
sys.path.insert(0,str(ROOT/'vendor'))
from legacy_v033_annotation import load_gtf,annotate_gene,Window,DEFAULT_CANONICAL_TAG,promoter_interval

GENE_REMOVE={'chrom','k','window_id','first_cpg_index','last_cpg_index','cpg_positions_1based','start0','end0','span_bp'}
TX_REMOVE={'chrom','k','window_id','cpg_positions_1based','start0','end0'}
EMPTY=dict(status='NO_DIRECT_PROTEIN_CODING_RELATION',gene_id='.',gene_name='.',primary_feature='NO_DIRECT_PROTEIN_CODING_RELATION',gene_count=0,transcript_count=0,canonical_promoter=0,any_transcript_promoter=0)

def window(chrom,pos0):
    return Window(chrom=chrom,k=1,window_id='annotation_anchor',first_cpg_index=0,last_cpg_index=0,
                  cpg_positions_1based=str(pos0+1),start0=pos0,end0=pos0+1,span_bp=1)

def profile(chrom,pos0,genes):
    gr=[];tr=[]
    for gene in sorted(genes,key=lambda g:g.gene_id):
        g,tx=annotate_gene(window(chrom,pos0),gene,DEFAULT_CANONICAL_TAG,2000,500)
        if g is None:continue
        gr.append({k:v for k,v in g.items() if k not in GENE_REMOVE})
        for t in tx:
            row={k:v for k,v in t.items() if k not in TX_REMOVE}
            row['ase_exon_union_transcript_eligible']=int(row['transcript_biotype'] in ('','NA','protein_coding'))
            for label,col in [('promoter','promoter_overlap_bp'),('five_prime_utr','five_prime_utr_overlap_bp'),
                ('three_prime_utr','three_prime_utr_overlap_bp'),('CDS','cds_overlap_bp'),('exon','exon_overlap_bp'),('intron','intron_overlap_bp')]:
                row['is_'+label]=int(row[col]>0)
            tr.append(row)
    summary=dict(EMPTY)
    if gr:
        chosen=[g for g in gr if g['include_primary_gene_feature_summary']==1]
        if chosen:
            rank={f:i for i,f in enumerate(POLICY['primary_order'])}
            best=min(chosen,key=lambda g:(rank[g['primary_feature']],-g['gene_body_overlap_bp'],g['gene_id']))
            summary.update(status=best['gene_relationship'],gene_id=best['gene_id'],gene_name=best['gene_name'],primary_feature=best['primary_feature'])
        else:
            status='NO_CANONICAL_PRIMARY_RELATION' if all(g['canonical_status']=='NO_ENSEMBL_CANONICAL' for g in gr) else 'ALTERNATIVE_ONLY_NO_PRIMARY_RELATION'
            summary.update(status=status,primary_feature=status)
        summary.update(gene_count=len(gr),transcript_count=len(tr),
            canonical_promoter=int(any(g['canonical_promoter_overlap'] for g in gr)),
            any_transcript_promoter=int(any(g['any_transcript_promoter_overlap'] for g in gr)))
    return dict(summary=summary,genes=gr,transcripts=tr)

def boundaries(genes):
    points=set();starts=defaultdict(list);ends=defaultdict(list)
    for g in genes:
        starts[g.relation_start0].append(g);ends[g.relation_end0].append(g)
        points.update((g.relation_start0,g.relation_end0,g.start0,g.end0))
        for t in g.transcripts.values():
            p=promoter_interval(t,2000,500);points.update((p.start0,p.end0,t.start0,t.end0))
            for intervals in (t.exons,t.cds,t.five_prime_utr,t.three_prime_utr):
                for v in intervals:points.update((v.start0,v.end0))
    return sorted(points),starts,ends

def build(gtf,out,expected_hash=GTF_SHA256,chromosomes=None,max_gib=8,min_free_gib=20):
    out=Path(out);gtf=Path(gtf);started=time.monotonic()
    if sha(gtf)!=expected_hash:raise ValueError('Reference GTF checksum differs')
    out.mkdir(parents=True,exist_ok=False);guard=Guard(out,max_gib,min_free_gib)
    genes,qc=load_gtf(gtf,set(chromosomes or [str(i) for i in range(1,19)]),DEFAULT_CANONICAL_TAG,2000,500)
    contract=dict(version=VERSION,policy=POLICY,gtf_sha256=expected_hash,chromosomes=sorted(genes,key=int),
        historical_engine_sha256=sha(ROOT/'vendor'/'legacy_v033_annotation.py'))
    db=sqlite3.connect(out/'profiles.sqlite');db.execute('CREATE TABLE profiles (id TEXT PRIMARY KEY, payload TEXT NOT NULL)')
    sf,sw=writer(out/'segments.tsv.gz',['chrom','start0','end0','annotation_profile_id'])
    seen=set();n=0;tx_n=0;gene_n=0;gh=th=None
    try:
        with sf:
            for chrom in sorted(genes,key=int):
                points,starts,ends=boundaries(genes[chrom]);active={}
                for left,right in zip(points,points[1:]):
                    for g in ends[left]:active.pop(g.gene_id,None)
                    for g in starts[left]:active[g.gene_id]=g
                    if not active:continue
                    payload=profile(chrom,left,active.values())
                    if not payload['genes']:continue
                    pid=digest(dict(contract=contract,chrom=chrom,payload=payload))
                    sw.writerow(dict(chrom=chrom,start0=left,end0=right,annotation_profile_id=pid));n+=1
                    if pid not in seen:
                        seen.add(pid);db.execute('INSERT INTO profiles VALUES (?,?)',(pid,json.dumps(payload,sort_keys=True,separators=(',',':'))))
                        for group in ('genes','transcripts'):
                            for row in payload[group]:
                                if group=='genes':
                                    if gh is None:gh,gw=writer(out/'gene_profiles.tsv.gz',['annotation_profile_id']+list(row))
                                    gw.writerow(dict(annotation_profile_id=pid,**row));gene_n+=1
                                else:
                                    if th is None:th,tw=writer(out/'transcript_profiles.tsv.gz',['annotation_profile_id']+list(row))
                                    tw.writerow(dict(annotation_profile_id=pid,**row));tx_n+=1
                    if n%2000==0:guard();db.commit()
                print('REFERENCE_CHROMOSOME='+chrom+' SEGMENTS='+str(n),flush=True)
        db.commit()
    finally:
        db.close()
        if gh:gh.close()
        if th:th.close()
    guard()
    for name in ('gene_profiles.tsv.gz','transcript_profiles.tsv.gz'):
        if not (out/name).exists():
            h,_=writer(out/name,['annotation_profile_id']);h.close()
    meta=dict(status='PASS_REFERENCE_ANNOTATION_INDEX',errors=[],contract=contract,gtf_path=str(gtf.resolve()),gtf_qc=qc,
        profiles=len(seen),segments=n,gene_relations=gene_n,transcript_relations=tx_n,elapsed_seconds=time.monotonic()-started,
        peak_rss_bytes=guard.peak,sha256={name:sha(out/name) for name in ('profiles.sqlite','segments.tsv.gz','gene_profiles.tsv.gz','transcript_profiles.tsv.gz')})
    write_json(out/'validation.json',meta)
    return meta

class Index:
    def __init__(self,directory,verify=True):
        self.directory=Path(directory);self.meta=read_json(self.directory/'validation.json')
        if self.meta.get('status')!='PASS_REFERENCE_ANNOTATION_INDEX' or self.meta.get('errors')!=[]:raise ValueError('Reference index is incomplete')
        if (self.meta['contract']['policy']!=POLICY or self.meta['contract'].get('version')!=VERSION or
                self.meta['contract'].get('historical_engine_sha256')!=sha(ROOT/'vendor'/'legacy_v033_annotation.py')):
            raise ValueError('Reference annotation policy/engine differs')
        if verify:
            for name,h in self.meta['sha256'].items():
                if sha(self.directory/name)!=h:raise ValueError('Reference index checksum differs: '+name)
        self.starts=defaultdict(list);self.ends=defaultdict(list);self.ids=defaultdict(list)
        for r in rows(self.directory/'segments.tsv.gz'):
            c=r['chrom'];a=int(r['start0']);b=int(r['end0'])
            if a>=b or (self.ends[c] and a<self.ends[c][-1]):raise ValueError('Overlapping reference segments')
            self.starts[c].append(a);self.ends[c].append(b);self.ids[c].append(r['annotation_profile_id'])
        self.db=sqlite3.connect('file:'+self.directory.joinpath('profiles.sqlite').resolve().as_posix()+'?mode=ro',uri=True)
    @lru_cache(maxsize=4096)
    def payload(self,pid):
        if pid=='.':return dict(summary=dict(EMPTY),genes=[],transcripts=[])
        v=self.db.execute('SELECT payload FROM profiles WHERE id=?',(pid,)).fetchone()
        if v is None:raise ValueError('Reference profile is missing')
        return json.loads(v[0])
    def locate(self,chrom,pos0):
        c=chrom.removeprefix('chr');j=bisect_right(self.starts[c],pos0)-1
        return self.ids[c][j] if j>=0 and pos0<self.ends[c][j] else '.'
    def annotation(self,chrom,pos1):
        if chrom.removeprefix('chr') not in self.meta['contract']['chromosomes'] or pos1<1:
            raise ValueError('Coordinate is outside the reference index scope')
        pid=self.locate(chrom,pos1-1);s=self.payload(pid)['summary']
        return dict(annotation_profile_id=pid,annotation_status=s['status'],annotation_gene_id=s['gene_id'],
            annotation_gene_name=s['gene_name'],annotation_primary_feature=s['primary_feature'],annotation_gene_count=str(s['gene_count']),
            annotation_transcript_count=str(s['transcript_count']),annotation_canonical_promoter=str(s['canonical_promoter']),
            annotation_any_transcript_promoter=str(s['any_transcript_promoter']),
            annotation_dyad_crosses_segment_boundary=str(int(self.locate(chrom,pos1)!=pid)))
    def close(self):self.db.close();self.payload.cache_clear()
