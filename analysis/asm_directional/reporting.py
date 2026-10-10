"""All selected coordinates, all eligible profiles, readable counts and diagnostics."""
import csv,gzip,math,collections,itertools
import numpy as np
from common import *
from pipeline import table,STATE,GT,rows
from stats import display,LOG_MIN

def assemble(job):
    out,reused=begin(job['parent'],job['identity'])
    if reused:return str(out)
    config=job['config'];families=[read(Path(p)/'summary.json') for p in job['families']]
    selected=collections.defaultdict(set);hits=[]
    for path in job['families']:
        for r in rows(Path(path)/'hits_BH10.tsv.gz'):
            hits.append(r);selected[int(r['chrom'])].add(int(r['cpg_pos1']))
    fields=['chrom','cpg_pos1','tissue','direction','r','N','log_PC_P','PC_P','log_BH_q','BH_q','log_direction_adjusted_q','direction_adjusted_q','BH05','display_underflow','family_size']
    table(out/'directional_recurrent_BH10.tsv.gz',hits,fields)
    source_ps={(u['sample'],u['tissue']):read(Path(u['cache']['path'])/'phase_dictionary.json') for u in config['source_inputs']['units']}
    anchor_ps={}
    for s in IDS:
        sp=config['locations']['phase_index'][f'{s}-chr1'];anchor_ps[s]=read(read(Path(sp)/'summary.json')['phase_dictionary'])
    candidates=[];candidate_sites=[];audit_profiles=[];source_wanted=collections.defaultdict(set);whole_N=collections.Counter();whole_states=collections.Counter()
    for path in job['tiles']:
        guard(out,fixture=job.get('fixture',False));p=Path(path);summary=read(p/'summary.json');c=summary['chrom'];sites=np.load(p/'sites.npy',mmap_mode='r');pr=np.load(p/'profiles.npy',mmap_mode='r')
        whole_states.update(summary['states'])
        for t,h in summary['N_histogram'].items():whole_N.update({f'{t}:{i}':int(x) for i,x in enumerate(h)})
        keep=set(selected[c]);site_map={int(a['pos']):a for a in sites if int(a['pos']) in keep}
        for pos,a in site_map.items():
            k=int(a['anchor']);row=dict(chrom=c,cpg_pos1=pos,anchor_pos1=k//16,anchor_ref='ACGT'[(k%16)//4],anchor_alt='ACGT'[k%4],anchor_testable_animals=int(a['anchor_testable']),anchor_available_animals=int(a['anchor_available']))
            for ti,t in enumerate(TISSUES):row['N_before_'+t]=int(a['N_before'][ti]);row['N_aligned_'+t]=int(a['N'][ti])
            candidate_sites.append(row)
        # The original pilot is returned even if no hypotheses pass BH.
        oldpilot=Path(config['locations']['pilot'][f'chr{c}']);pilotpos={int(r['cpg_pos1']) for r in rows(oldpilot/'pilot_sites.tsv.gz')}
        wanted=keep|pilotpos
        for r in pr[np.isin(pr['pos'],list(wanted))]:
            s=IDS[int(r['sample'])];t=TISSUES[int(r['tissue'])];pos=int(r['pos']);ix=int(np.searchsorted(sites['pos'],pos));a=sites[ix];k=int(a['anchor'])
            A,B,C,D=[int(r[x]) for x in ['a','b','c','d']]
            row=dict(sample=s,tissue=t,chrom=c,cpg_pos1=pos,ps=source_ps[(s,t)][int(r['ps'])],source_row=int(r['source_row']),
                H1_methylated=A,H1_unmethylated=B,H2_methylated=C,H2_unmethylated=D,original_fisher_two_sided_p=float(r['original_p']),original_global30_q=float(r['original_q']),
                anchor_pos1=k//16 if k else '.',anchor_ref='ACGT'[(k%16)//4] if k else '.',anchor_alt='ACGT'[k%4] if k else '.',
                anchor_GT=GT[int(r['GT'])],anchor_PS=anchor_ps[s][int(r['anchor_ps'])] if r['anchor_ps'] else '.',orientation_status=STATE[int(r['state'])],
                ALT_methylated='.',ALT_unmethylated='.',REF_methylated='.',REF_unmethylated='.',ALT_methylation='.',REF_methylation='.',delta_ALT_REF='.',
                logP_ALT_gt_REF='.',logP_REF_gt_ALT='.',P_ALT_gt_REF='.',P_REF_gt_ALT='.',aligned_N=int(a['N'][int(r['tissue'])]))
            if r['state']==0:
                if r['alt_hap']==2:A,B,C,D=C,D,A,B
                row.update(ALT_methylated=A,ALT_unmethylated=B,REF_methylated=C,REF_unmethylated=D,ALT_methylation=A/(A+B),REF_methylation=C/(C+D),delta_ALT_REF=A/(A+B)-C/(C+D),
                    logP_ALT_gt_REF=float(r['logp'][0]),logP_REF_gt_ALT=float(r['logp'][1]),P_ALT_gt_REF=float(display(r['logp'][0])),P_REF_gt_ALT=float(display(r['logp'][1])))
            if pos in keep:candidates.append(row);source_wanted[(s,t,c)].add(int(r['source_row']))
            if pos in pilotpos:audit_profiles.append(row)
    require(len(audit_profiles)==config['expected_pilot_profiles'],'Incomplete original pilot export')
    profile_fields=list((candidates or audit_profiles)[0]) if candidates or audit_profiles else ['sample','tissue','chrom','cpg_pos1']
    table(out/'candidate_all_eligible_profiles.tsv.gz',candidates,profile_fields)
    table(out/'audit_pilot_profiles.tsv.gz',audit_profiles,profile_fields)
    sitefields=list(candidate_sites[0]) if candidate_sites else ['chrom','cpg_pos1','anchor_pos1','anchor_ref','anchor_alt']
    table(out/'candidate_sites.tsv.gz',candidate_sites,sitefields)
    # Preserve all eligible animals when describing sign/size, not only nominal hits.
    grouped=collections.defaultdict(list)
    for r in candidates:grouped[(r['chrom'],r['cpg_pos1'],r['tissue'])].append(r)
    effects=[]
    for (c,pos,t),rr in sorted(grouped.items()):
        good=[r for r in rr if r['orientation_status']=='ALIGNED'];delta=np.array([r['delta_ALT_REF'] for r in good],dtype=float)
        x=dict(chrom=c,cpg_pos1=pos,tissue=t,N_aligned=len(good),eligible_profiles=len(rr),aligned_animals=';'.join(r['sample'] for r in good),
            ALT_high_animals=';'.join(r['sample'] for r in good if r['ALT_methylated']*r['REF_unmethylated']>r['REF_methylated']*r['ALT_unmethylated']),
            REF_high_animals=';'.join(r['sample'] for r in good if r['ALT_methylated']*r['REF_unmethylated']<r['REF_methylated']*r['ALT_unmethylated']),
            tied_animals=';'.join(r['sample'] for r in good if r['ALT_methylated']*r['REF_unmethylated']==r['REF_methylated']*r['ALT_unmethylated']),
            nominal_ALT_P05_animals=';'.join(r['sample'] for r in good if r['logP_ALT_gt_REF']<math.log(.05)),
            nominal_REF_P05_animals=';'.join(r['sample'] for r in good if r['logP_REF_gt_ALT']<math.log(.05)),
            median_delta=float(np.median(delta)) if len(delta) else '.',minimum_delta=float(delta.min()) if len(delta) else '.',maximum_delta=float(delta.max()) if len(delta) else '.')
        effects.append(x)
    table(out/'candidate_effect_summary.tsv.gz',effects,list(effects[0]) if effects else ['chrom','cpg_pos1','tissue'])
    table(out/'family_summary.tsv.gz',families,list(families[0]))
    write(out/'selected_source_rows.json',{s+'|'+t+'|'+str(c):sorted(v) for (s,t,c),v in source_wanted.items()})
    both=[];hit_lookup=collections.defaultdict(set)
    for r in hits:
        if int(r['BH05']):hit_lookup[(r['chrom'],r['cpg_pos1'],r['tissue'],r['r'])].add(r['direction'])
    for k,v in hit_lookup.items():
        if len(v)==2:both.append(dict(zip(['chrom','cpg_pos1','tissue','r'],k)))
    table(out/'both_directions_BH05.tsv.gz',both,['chrom','cpg_pos1','tissue','r'])
    return seal(job['parent'],out,job['identity'],dict(families=len(families),hypotheses=sum(s['hypotheses'] for s in families),BH05_rows=sum(s['BH05'] for s in families),BH10_rows=len(hits),
        unique_BH05_coordinates=len({(r['chrom'],r['cpg_pos1']) for r in hits if int(r['BH05'])}),unique_BH10_coordinates=len(candidate_sites),candidate_profiles=len(candidates),
        pilot_profiles=len(audit_profiles),both_direction_BH05_hypotheses=len(both),whole_N_histogram=dict(whole_N),whole_orientation_states=dict(whole_states),
        note='Unique-coordinate union is descriptive. No union/max-r FDR claim. Same-tissue directional results do not filter original tissue OR analysis.'))

def annotation_job(job):
    out,reused=begin(job['parent'],job['identity'])
    if reused:return str(out)
    spec=job['spec'];cp=Path(spec['cache']['path']);verify(cp/'validation.json',spec['cache']['gate_sha256']);gate=read(cp/'validation.json')
    output=[];seen=set();fields=None;source_snapshots={}
    for c,want in job['wanted'].items():
        name=f'chr{c}.source.tsv.gz';require(name in gate['sha256'],'Missing authenticated annotation sidecar');source_snapshots[str(cp/name)]=verify(cp/name,gate['sha256'][name])
        wanted=set(want)
        with gzip.open(cp/name,'rt',encoding='utf-8',newline='') as f:
            reader=csv.DictReader(f,delimiter='\t')
            require(fields is None or fields==reader.fieldnames,'Annotation schema differs across chromosomes');fields=reader.fieldnames
            for r in reader:
                if int(r['source_row']) in wanted:
                    require(r['sample']==spec['sample'] and r['tissue']==spec['tissue'],'Annotation unit mismatch')
                    output.append(r);seen.add((int(c),int(r['source_row'])))
    for path,before in source_snapshots.items():require(snapshot(path)==before,'Annotation source changed while reading: '+path)
    write(out/'source_snapshots.json',source_snapshots)
    require(len(seen)==sum(len(w) for w in job['wanted'].values())==len(output),'Missing/duplicate selected source rows')
    table(out/'selected_original_annotations.tsv.gz',output,fields or ['sample','tissue','chrom','cpg_pos1'])
    return seal(job['parent'],out,job['identity'],dict(sample=spec['sample'],tissue=spec['tissue'],rows=len(output),
        note='Original annotation/source fields preserved; no nearest-gene assignment or reannotation.'))
