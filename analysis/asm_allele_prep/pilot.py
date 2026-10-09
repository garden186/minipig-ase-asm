"""P-independent opportunity pilot plus a separate frozen-coordinate regression."""
from __future__ import annotations
import csv, gzip, itertools, time
import numpy as np
from common import *
from indexing import keys, STATUS, GT, BASE

BASES='ACGT'
POPCOUNT=np.array([i.bit_count() for i in range(1024)],dtype='u1')

def rows(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt',encoding='utf-8',newline='') as f:yield from csv.DictReader(f,delimiter='\t')

def table(path,data,fields):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'wt',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)

def at(a,pos):
    lo,hi=np.searchsorted(a['pos'],pos,side='left'),np.searchsorted(a['pos'],pos,side='right')
    return a[lo:hi]

def choose_anchor(pos,available,tested,phase_arrays,phase_maps):
    """Only exact PS membership and animal opportunity enter this ranking."""
    keyparts=[];abits=[];tbits=[]
    for i in range(len(phase_arrays)):
        a=phase_arrays[i];pm=phase_maps[i]
        for ps in sorted(available[i]):
            code=pm.get(ps)
            if code is None:continue
            lo,hi=np.searchsorted(a['ps'],code,side='left'),np.searchsorted(a['ps'],code,side='right')
            k=keys(a[lo:hi])
            if not len(k):continue
            keyparts.append(k);abits.append(np.full(len(k),1<<i,dtype='<u2'))
            tbits.append(np.full(len(k),(1<<i) if ps in tested[i] else 0,dtype='<u2'))
    if not keyparts:return None,0,0,0
    k,ix=np.unique(np.concatenate(keyparts),return_inverse=True)
    am=np.zeros(len(k),dtype='<u2');tm=am.copy()
    np.bitwise_or.at(am,ix,np.concatenate(abits));np.bitwise_or.at(tm,ix,np.concatenate(tbits))
    na=POPCOUNT[am];nt=POPCOUNT[tm];dist=np.abs((k//16).astype(np.int64)-pos)
    order=np.lexsort((k,dist,-na.astype(int),-nt.astype(int)))
    j=int(order[0]);return int(k[j]),int(nt[j]),int(na[j]),len(k)

def independent_selection(unit_arrays,ids,per_stratum):
    animal_pos=[];multi=[]
    for s in ids:
        groups=[a['pos'] for (sid,t),a in unit_arrays.items() if sid==s]
        animal_pos.append(np.unique(np.concatenate(groups)))
    pos,n=np.unique(np.concatenate(animal_pos),return_counts=True)
    chosen=set();strata=[]
    for count in range(1,len(ids)+1):
        p=pos[n==count];indices=np.unique(np.rint(np.linspace(0,len(p)-1,min(per_stratum,len(p)))).astype(int)) if len(p) else []
        chosen.update(map(int,p[indices]))
        strata.append(dict(N_union=count,coordinates=len(p),selected=len(indices)))
    for a in unit_arrays.values():
        if len(a)>1:multi.extend(map(int,a['pos'][1:][a['pos'][1:]==a['pos'][:-1]]))
    return pos,n,chosen,set(multi),strata

def pilot_job(job):
    out,reuse=begin(job['parent'],job['identity'])
    if reuse:return str(out)
    start=time.monotonic();guard(out);c=job['chrom'];ids=job['ids'];tissues=job['tissues']
    unit={};avail={};ps={};cacheps={}
    for spec in job['units']:
        k=(spec['sample'],spec['tissue']);up=Path(job['unit_paths']['-'.join(k)])
        unit[k]=np.load(Path(spec['cache']['path'])/f'chr{c}.npy',mmap_mode='r',allow_pickle=False)
        avail[k]=np.load(up/f'chr{c}.available.npy',mmap_mode='r',allow_pickle=False)
        ps[k]=read(up/'phase_dictionary.json');cacheps[k]=read(Path(spec['cache']['path'])/'phase_dictionary.json')
    arrays=[];maps=[];allvars=[];allkeys=[];dicts=[]
    for s in ids:
        p=Path(job['join_paths'][f'{s}-chr{c}'])
        arrays.append(np.load(p/'eligible_by_ps.npy',mmap_mode='r',allow_pickle=False))
        v=np.load(p/'variants.npy',mmap_mode='r',allow_pickle=False);allvars.append(v);allkeys.append(keys(v))
        d=read(read(p/'summary.json')['phase_dictionary']);dicts.append(d);maps.append({v:i for i,v in enumerate(d)})
    positions,ns,independent,multi,strata=independent_selection(unit,ids,job['per_stratum'])
    np.save(out/'eligible_positions.npy',positions,allow_pickle=False)
    np.save(out/'eligible_union_N.npy',ns.astype('u1'),allow_pickle=False)
    full_tissue_N={}
    for t in tissues:
        pp,nn=np.unique(np.concatenate([np.unique(a['pos']) for (s,tt),a in unit.items() if tt==t]),return_counts=True)
        np.save(out/f'eligible_positions_{t}.npy',pp,allow_pickle=False)
        np.save(out/f'eligible_N_{t}.npy',nn.astype('u1'),allow_pickle=False)
        full_tissue_N[t]=dict(coordinates=len(pp),N_histogram=np.bincount(nn,minlength=len(ids)+1).tolist())
    del pp,nn
    frozen={int(r['cpg_pos1']):r for r in rows(job['frozen_anchors']) if int(r['chrom'])==c}
    oldrows={(r['sample'],r['tissue'],int(r['cpg_pos1']),r['ps']):r for r in rows(job['frozen_profiles']) if int(r['chrom'])==c}
    chosen=sorted(independent|multi|set(frozen));sites=[];profiles=[];mismatches=[];reg_checked=0
    seen_old=set();reason_counts={};n_distributions={}
    for z,pos in enumerate(chosen):
        if z%100==0:guard(out);print(f'PILOT=chr{c} SITES={z}/{len(chosen)}',flush=True)
        available=[set() for _ in ids];tested=[set() for _ in ids];byunit={}
        for (s,t),a in unit.items():
            i=ids.index(s);rr=at(a,pos);byunit[(s,t)]=rr
            tested[i].update(cacheps[(s,t)][int(r['ps'])] for r in rr)
            available[i].update(ps[(s,t)][int(r['ps'])] for r in at(avail[(s,t)],pos))
        anchor,nt,na,ncand=choose_anchor(pos,available,tested,arrays,maps)
        tags=';'.join(x for flag,x in [(pos in independent,'independent'),(pos in frozen,'frozen'),(pos in multi,'multiple_PS')] if flag)
        ap=anchor//16 if anchor is not None else '.';ref=BASES[(anchor%16)//4] if anchor is not None else '.'
        alt=BASES[anchor%4] if anchor is not None else '.'
        site=dict(chrom=c,cpg_pos1=pos,groups=tags,anchor_pos1=ap,anchor_ref=ref,anchor_alt=alt,
                  anchor_testable_animals=nt,anchor_available_animals=na,candidates=ncand)
        before={t:set() for t in tissues};after={t:set() for t in tissues};aligned_ps={t:{} for t in tissues}
        if pos in frozen:
            old=frozen[pos]
            expected=(old['anchor_pos1'],old['anchor_ref'],old['anchor_alt'],int(old['n_testable_animals']),int(old['n_available_animals']))
            actual=(str(ap),ref,alt,nt,na)
            if actual!=expected:mismatches.append(dict(type='anchor',chrom=c,pos=pos,expected=expected,actual=actual))
        for (s,t),rr in byunit.items():
            i=ids.index(s);v=None
            if anchor is not None:
                q=np.searchsorted(allkeys[i],anchor)
                if q<len(allkeys[i]) and int(allkeys[i][q])==anchor:v=allvars[i][q]
            for r in rr:
                phase=cacheps[(s,t)][int(r['ps'])];before[t].add(s)
                state='NO_ELIGIBLE_SITE_ANCHOR' if anchor is None else 'NO_MATCHING_VARIANT_RECORD' if v is None else STATUS[int(v['status'])] if v['status'] else 'DIFFERENT_EXACT_PS' if dicts[i][int(v['ps'])]!=phase else 'ALIGNED'
                aa,bb,cc,dd=(int(r[k]) for k in ['a','b','c','d'])
                row=dict(sample=s,tissue=t,chrom=c,cpg_pos1=pos,ps=phase,groups=tags,H1_methylated=aa,H1_unmethylated=bb,
                    H2_methylated=cc,H2_unmethylated=dd,original_fisher_two_sided_p=float(r['p']),original_q_all30=float(r['q_all']),
                    anchor_pos1=ap,anchor_ref=ref,anchor_alt=alt,anchor_GT=GT[int(v['gt'])] if v is not None else '.',
                    anchor_PS=dicts[i][int(v['ps'])] if v is not None else '.',orientation_status=state,ALT_methylated='.',ALT_unmethylated='.',
                    REF_methylated='.',REF_unmethylated='.',ALT_methylation='.',REF_methylation='.',delta_ALT_REF='.')
                if state=='ALIGNED':
                    A,B,C,D=(aa,bb,cc,dd) if v['alt_hap']==1 else (cc,dd,aa,bb)
                    row.update(ALT_methylated=A,ALT_unmethylated=B,REF_methylated=C,REF_unmethylated=D,
                        ALT_methylation=A/(A+B),REF_methylation=C/(C+D),delta_ALT_REF=A/(A+B)-C/(C+D))
                    after[t].add(s);require(s not in aligned_ps[t],'Multiple aligned PS rows for one animal/tissue/site')
                    aligned_ps[t][s]=phase
                for tag in tags.split(';'):reason_counts[(tag,t,state)]=reason_counts.get((tag,t,state),0)+1
                key=(s,t,pos,phase)
                if pos in frozen:
                    old=oldrows.get(key)
                    if old is None:mismatches.append(dict(type='unexpected_eligible_profile',key=key))
                    else:
                        seen_old.add(key);reg_checked+=1;bad=[]
                        for field in ['H1_methylated','H1_unmethylated','H2_methylated','H2_unmethylated','orientation_status','anchor_GT','anchor_PS']:
                            if str(row[field])!=old[field]:bad.append(field)
                        for newfield,oldfield in [('original_fisher_two_sided_p','fisher_two_sided_p'),('original_q_all30','BH_q_all30'),('ALT_methylation','ALT_methylation'),('REF_methylation','REF_methylation'),('delta_ALT_REF','delta_ALT_REF')]:
                            x,y=row[newfield],old[oldfield]
                            tolerance=0 if newfield.startswith('original_') else 1e-15
                            if (str(x)=='.')!=(y=='.') or (str(x)!='.' and not np.isclose(float(x),float(y),rtol=2e-12,atol=tolerance)):bad.append(newfield)
                        if bad:mismatches.append(dict(type='profile',key=key,fields=bad))
                profiles.append(row)
        for t in tissues:
            site['N_before_'+t]=len(before[t]);site['N_aligned_'+t]=len(after[t])
            for tag in tags.split(';'):
                k=(tag,t,len(before[t]),len(after[t]));n_distributions[k]=n_distributions.get(k,0)+1
        for t,u in itertools.combinations(tissues,2):
            shared=after[t]&after[u]
            site['N_aligned_'+t+'_'+u]=len(shared)
            site['N_aligned_exactPS_'+t+'_'+u]=sum(aligned_ps[t][s]==aligned_ps[u][s] for s in shared)
        sites.append(site)
    for key in oldrows.keys()-seen_old:mismatches.append(dict(type='missing_eligible_profile',key=key))
    table(out/'pilot_sites.tsv.gz',sites,list(sites[0]) if sites else ['chrom','cpg_pos1'])
    table(out/'pilot_profiles.tsv.gz',profiles,list(profiles[0]) if profiles else ['sample','tissue','chrom','cpg_pos1'])
    write(out/'regression.json',dict(checked=reg_checked,expected=len(oldrows),errors=mismatches))
    write(out/'pilot_summary.json',dict(strata=strata,full_tissue_eligible_N=full_tissue_N,
        reason_counts=[dict(group=k[0],tissue=k[1],state=k[2],rows=v) for k,v in sorted(reason_counts.items())],
        N_distribution=[dict(group=k[0],tissue=k[1],N_before=k[2],N_aligned=k[3],sites=v) for k,v in sorted(n_distributions.items())]))
    require(not mismatches,'Frozen regression mismatch; inspect '+str(out/'regression.json'))
    return seal(job['parent'],out,job['identity'],dict(chrom=c,full_eligible_coordinates=len(positions),
        independent_sites=len(independent),frozen_sites=len(frozen),multiple_PS_sites=len(multi),pilot_sites=len(sites),
        pilot_eligible_rows=len(profiles),regression_checked=reg_checked,regression_errors=0,seconds=time.monotonic()-start,
        scope='Full eligible opportunity inventory; allele alignment only on pilot and frozen regression coordinates. No new P/q.'))
