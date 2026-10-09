"""Return discoveries, full-family summaries and matched inference/effect counts."""
import collections,gzip,math,shutil
import numpy as np
from common import *
from data import *
from schema import STATE,GT
from stats import display

HIT_FIELDS=['pair','r','chrom','cpg_pos1','N_paired','N_first_all','N_second_all',
    'first_ALT_high_N','first_REF_high_N','first_tie_N','second_ALT_high_N','second_REF_high_N','second_tie_N',
    'D_positive_N','D_negative_N','D_zero_N',
    'gene_id','gene_name','gene_assignment','anchor','mean_D_pp','median_D_pp','min_D_pp','max_D_pp',
    'median_first_all_pp','median_second_all_pp','median_first_paired_pp','median_second_paired_pp',
    'singleton_animals','OR_PC_P','log_OR_PC_P','BH_q','log_BH_q','BY_q','log_BY_q','BH05','BH10','BY05','M',
    'first_existing_R_BH05_r_directions','second_existing_R_BH05_r_directions']
PROFILE_FIELDS=['chrom','cpg_pos1','sample','tissue','source_row','source_ps_index','anchor_ps_index','anchor_GT','orientation_status',
    'H1_M','H1_U','H2_M','H2_U','ALT_M','ALT_U','REF_M','REF_U','delta_ALT_REF','anchor','N_aligned_tissue']
PAIR_FIELDS=['pair','chrom','cpg_pos1','sample','paired','first_delta','second_delta','D_pp','D_sign',
    'OR_P','log_OR_P','full_states']

def anchor_text(chrom,k):
    k=int(k)
    return f"{chrom}:{k//16}:{'ACGT'[(k%16)//4]}:{'ACGT'[k%4]}" if k else '.'

def existing_recurrence():
    d=collections.defaultdict(list)
    for r in rows(asset_path('assets/directional_recurrent_BH10.tsv.gz')):
        if r['BH05']=='1':d[(int(r['chrom']),int(r['cpg_pos1']),r['tissue'])].append(f"{r['direction']}:r{r['r']}:q{r['direction_adjusted_q']}")
    return d

def hits_record(a,z,meta,genes,recurrence):
    p=meta['pair_index'];first,second=[TISSUES[t] for t in PAIRS[p]];r=meta['r'];c=int(a['chrom']);pos=int(a['pos']);g=genes[int(z['gene'])]
    x=dict(pair=PAIR_NAMES[p],r=r,chrom=c,cpg_pos1=pos,N_paired=int(z['N']),N_first_all=int(z['N_first']),N_second_all=int(z['N_second']),
        gene_id=g['gene_id'],gene_name=g['gene_name'],gene_assignment='frozen_primary_protein_coding_gene; no nearest-gene assignment',anchor=anchor_text(c,z['anchor']),
        singleton_animals=int(z['singleton_animals']),M=meta['M'])
    for prefix,field in [('first','sign_counts_first'),('second','sign_counts_second')]:
        for k,value in zip(['ALT_high_N','REF_high_N','tie_N'],z[field]):x[prefix+'_'+k]=int(value)
    for k,value in zip(['D_positive_N','D_negative_N','D_zero_N'],z['D_sign_counts']):x[k]=int(value)
    for key in ['mean_D','median_D','min_D','max_D','median_first_all','median_second_all','median_first_paired','median_second_paired']:x[key+'_pp']=float(z[key])*100
    for name,field,logname in [('OR_PC_P','logp','log_OR_PC_P'),('BH_q','logq_BH','log_BH_q'),('BY_q','logq_BY','log_BY_q')]:
        x[name]=float(display(a[field]));x[logname]=float(a[field])
    x.update(BH05=int(a['logq_BH']<=math.log(.05)),BH10=int(a['logq_BH']<=math.log(.10)),BY05=int(a['logq_BY']<=math.log(.05)),
        first_existing_R_BH05_r_directions=';'.join(recurrence.get((c,pos,first),[])) or '.',
        second_existing_R_BH05_r_directions=';'.join(recurrence.get((c,pos,second),[])) or '.')
    return x

def select(job):
    out,reused=begin(job['parent'],job['identity'])
    if reused:return str(out)
    genes=read(asset_path('assets/genes.json'));recurrence=existing_recurrence();selection=collections.defaultdict(set)
    hits=[];top=[];summaries=[];cache={}
    def getz(tile):
        if tile not in cache:cache[tile]=np.load(Path(job['tiles'][tile])/'results.npy',mmap_mode='r',allow_pickle=False)
        return cache[tile]
    for path in job['families']:
        s=read(Path(path)/'summary.json');a=np.load(Path(path)/'family.npy',mmap_mode='r',allow_pickle=False);summaries.append(s)
        chosen=np.flatnonzero(a['logq_BH']<=math.log(.10))
        for target,indices in [(hits,chosen),(top,range(min(25,len(a))))]:
            for i in indices:
                record=a[i];tile=int(record['tile']);row=int(record['row']);z=getz(tile)[row]
                require(z['pos']==record['pos'] and int(z['pair'])==s['pair_index'],'Family-to-tile row mismatch')
                target.append(hits_record(record,z,s,genes,recurrence));selection[tile].add(row)
    table(out/'tissue_contrasts_BH10.tsv.gz',hits,HIT_FIELDS);table(out/'top25_per_family.tsv.gz',top,HIT_FIELDS)
    simple=[{k:s[k] for k in ['pair','r','M','BH05','BH10','BY05','min_P_log','min_BH_log']} for s in summaries]
    table(out/'family_summary.tsv',simple,list(simple[0]))
    write(out/'return_selection.json',{str(k):sorted(v) for k,v in selection.items()})
    pair_counts={}
    for p,name in enumerate(PAIR_NAMES):
        before=np.zeros(11,dtype=np.int64);after=before.copy();reasons=collections.Counter();mixed=ties=0
        for path in job['tiles']:
            s=read(Path(path)/'summary.json')['per_pair'][name]
            before+=s['available_by_N'];after+=s['tested_by_N'];reasons.update(s['reasons'])
            mixed+=s['retained_with_mixed_signs'];ties+=s['retained_with_ties']
        require(np.array_equal(before[2:],after[2:]),'Full paired N inventory lost')
        pair_counts[name]=dict(available_by_N=before.tolist(),tested_by_N=after.tolist(),reasons=dict(reasons),retained_with_mixed_signs=mixed,retained_with_ties=ties)
    filterrows=[dict(pair=name,N=N,available=s['available_by_N'][N],tested=s['tested_by_N'][N]) for name,s in pair_counts.items() for N in range(11)]
    table(out/'paired_N_summary.tsv',filterrows,list(filterrows[0]))
    plot_groups=collections.defaultdict(list)
    for x in hits:
        if x['BH05']:plot_groups[(x['pair'],x['chrom'],x['cpg_pos1'])].append(x)
    plot=[]
    for key,values in sorted(plot_groups.items()):
        maximum=max(values,key=lambda x:x['r']);rec={k:maximum[k] for k in ['pair','chrom','cpg_pos1','gene_id','gene_name','anchor','N_paired','median_D_pp','mean_D_pp','min_D_pp','max_D_pp']}
        rec.update(max_significant_r=maximum['r'],significant_r=','.join(str(x['r']) for x in sorted(values,key=lambda x:x['r'])),BH_q_at_max_r=maximum['BH_q'],BY_q_at_max_r=maximum['BY_q'])
        plot.append(rec)
    fields=['pair','chrom','cpg_pos1','gene_id','gene_name','anchor','N_paired','median_D_pp','mean_D_pp','min_D_pp','max_D_pp','max_significant_r','significant_r','BH_q_at_max_r','BY_q_at_max_r']
    table(out/'pairwise_plot_source_BH05.tsv',plot,fields)
    return seal(job['parent'],out,job['identity'],dict(families=27,BH10_hypotheses=len(hits),BH05_hypotheses=sum(x['BH05'] for x in hits),
        unique_BH05_CpGs=len({(x['chrom'],x['cpg_pos1']) for x in hits if x['BH05']}),
        BY05_hypotheses=sum(x['BY05'] for x in hits),returned_pair_sites=sum(map(len,selection.values())),
        pair_availability=pair_counts,unique_BH05_by_pair={p:len({x['cpg_pos1']*100+x['chrom'] for x in plot if x['pair']==p}) for p in PAIR_NAMES},
        note='Unique/maximum-r summaries are descriptive, not union-level FDR. Top25 includes non-significant rows. No observed-sign filter.'))

def export_tile(job):
    out,reused=begin(job['parent'],job['identity'])
    if reused:return str(out)
    source=Path(job['tile']);s=read(source/'summary.json');chrom=s['chrom'];original=Path(s['source'])
    # Reverify the exact source arrays before exporting biological counts.
    for name in ['profiles.npy','sites.npy']:
        verify(original/name,job['spec']['checks'][str(original/name).replace('\\','/')])
    sites=np.load(original/'sites.npy',mmap_mode='r',allow_pickle=False);pr=np.load(original/'profiles.npy',mmap_mode='r',allow_pickle=False)
    z=np.load(source/'results.npy',mmap_mode='r',allow_pickle=False);ind=np.load(source/'individual.npy',mmap_mode='r',allow_pickle=False)
    chosen=z[job['rows']];positions=set(map(int,chosen['pos']));profiles=[];aligned={}
    for r in pr[np.isin(pr['pos'],list(positions))]:
        pos=int(r['pos']);a=sites[np.searchsorted(sites['pos'],pos)];animal=int(r['sample']);t=int(r['tissue'])
        hm,hu,rm,ru=[int(r[k]) for k in 'abcd']
        x=dict(chrom=chrom,cpg_pos1=pos,sample=IDS[animal],tissue=TISSUES[t],source_row=int(r['source_row']),source_ps_index=int(r['ps']),
            anchor_ps_index=int(r['anchor_ps']),anchor_GT=GT[int(r['GT'])],orientation_status=STATE[int(r['state'])],H1_M=hm,H1_U=hu,H2_M=rm,H2_U=ru,
            ALT_M='.',ALT_U='.',REF_M='.',REF_U='.',delta_ALT_REF='.',anchor=anchor_text(chrom,a['anchor']),N_aligned_tissue=int(a['N'][t]))
        if r['state']==0:
            if r['alt_hap']==2:hm,hu,rm,ru=rm,ru,hm,hu
            delta=hm/(hm+hu)-rm/(rm+ru);x.update(ALT_M=hm,ALT_U=hu,REF_M=rm,REF_U=ru,delta_ALT_REF=delta)
            key=(pos,animal,t);require(key not in aligned,'Duplicated export animal/tissue');aligned[key]=delta
        profiles.append(x)
    paired=[]
    for index in job['rows']:
        result=z[index];record=ind[index];pos=int(result['pos']);p=int(result['pair']);u,v=PAIRS[p]
        for animal in range(10):
            d1=aligned.get((pos,animal,u));d2=aligned.get((pos,animal,v));both=d1 is not None and d2 is not None
            x=dict(pair=PAIR_NAMES[p],chrom=chrom,cpg_pos1=pos,sample=IDS[animal],paired=int(both),first_delta=d1 if d1 is not None else '.',second_delta=d2 if d2 is not None else '.',
                D_pp=100*float(record['D'][animal]) if both else '.',D_sign=int(record['D_sign'][animal]) if both else '.',OR_P='.',log_OR_P='.',full_states='.')
            if both:
                require(np.isfinite(record['logp'][animal]) and abs(float(record['D'][animal])-(d1-d2))<1e-12,'Export differs from inference animals or effects')
                x.update(OR_P=float(display(record['logp'][animal])),log_OR_P=float(record['logp'][animal]),full_states=int(record['full_states'][animal]))
            else:require(not np.isfinite(record['logp'][animal]),'Unavailable animal was included in inference')
            paired.append(x)
    table(out/'all_eligible_profiles.tsv.gz',profiles,PROFILE_FIELDS);table(out/'paired_individual_results.tsv.gz',paired,PAIR_FIELDS)
    return seal(job['parent'],out,job['identity'],dict(chrom=chrom,returned_CpGs=len(positions),profiles=len(profiles),paired_summary_rows=len(paired)))

def merge_exports(paths,output):
    for name,fields in [('all_eligible_profiles.tsv.gz',PROFILE_FIELDS),('paired_individual_results.tsv.gz',PAIR_FIELDS)]:
        def stream():
            for path in paths:yield from rows(Path(path)/name)
        table(Path(output)/name,stream(),fields)
