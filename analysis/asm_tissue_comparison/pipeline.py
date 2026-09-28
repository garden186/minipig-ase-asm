"""Analyze every aligned paired CpG using full-support exact OR and Fisher-PC."""
import collections,math,time
from pathlib import Path
import numpy as np
from scipy.special import digamma
from common import *
from data import *
from exact import pair_test,paired_difference
from stats import partial_conjunction,bh_log
from schema import STATE,GT


def analyze(counts,sites,genes,progress=lambda x:None):
    filt,signs,present=filters(counts,sites['pos'])
    delta=deltas(counts)
    selected=[(p,j) for p in range(3) for j in np.flatnonzero(filt['reason'][:,p]==0)]
    out=np.zeros(len(selected),dtype=RESULT);individual=np.zeros(len(selected),dtype=INDIVIDUAL)
    out['log_pc']=np.nan
    for field in ['logp','D']:individual[field]=np.nan
    for row,(p,j) in enumerate(selected):
        u,v=PAIRS[p];paired=present[j,:,u]&present[j,:,v];animals=np.flatnonzero(paired)
        z=out[row];z['pos']=sites['pos'][j];z['pair']=p;z['N']=len(animals)
        z['N_first']=filt['N'][j,u];z['N_second']=filt['N'][j,v]
        z['sign_counts_first']=[filt[key][j,u] for key in ['positive','negative','ties']]
        z['sign_counts_second']=[filt[key][j,v] for key in ['positive','negative','ties']]
        z['anchor']=sites['anchor'][j];z['gene']=genes[j]
        z['median_first_all']=np.median(delta[j,present[j,:,u],u]);z['median_second_all']=np.median(delta[j,present[j,:,v],v])
        z['median_first_paired']=np.median(delta[j,animals,u]);z['median_second_paired']=np.median(delta[j,animals,v])
        for animal in animals:
            one=tuple(map(int,counts[j,animal,u]));two=tuple(map(int,counts[j,animal,v]))
            lp,nfull=pair_test(one,two);difference,dsign=paired_difference(one,two)
            individual[row]['logp'][animal]=lp;individual[row]['full_states'][animal]=nfull
            individual[row]['D'][animal]=difference;individual[row]['D_sign'][animal]=dsign
        d=individual[row]['D'][animals];ds=individual[row]['D_sign'][animals]
        z['mean_D']=np.mean(d);z['median_D']=np.median(d);z['min_D']=np.min(d);z['max_D']=np.max(d)
        z['D_sign_counts']=[np.count_nonzero(ds>0),np.count_nonzero(ds<0),np.count_nonzero(ds==0)]
        z['singleton_animals']=np.count_nonzero(individual[row]['full_states'][animals]==1)
        if row%1000==0:progress(dict(completed=row,total=len(selected),cache=pair_test.cache_info()._asdict()))
    for N in range(2,11):
        ix=np.flatnonzero(out['N']==N)
        if not len(ix):continue
        for source,dest in [('logp','log_pc')]:
            arr=individual[source][ix];values=arr[np.isfinite(arr)].reshape(len(ix),N)
            for r in range(2,N+1):out[dest][ix,r-2]=partial_conjunction(values,r)
    return filt,out,individual


def tile_job(job):
    spec=job['spec'];snaps={p:verify(p,h) for p,h in spec['checks'].items()}
    out,reused=begin(job['parent'],job['identity'])
    if reused:return str(out)
    started=time.monotonic();guard(out,fixture=job.get('fixture',False))
    src=Path(spec['path']);sites=np.load(src/'sites.npy',mmap_mode='r',allow_pickle=False)
    pr=np.load(src/'profiles.npy',mmap_mode='r',allow_pickle=False)
    require(len(sites)==spec['sites'] and len(pr)==spec['profiles'],'Tile input count changed')
    require(np.count_nonzero(pr['state']==0)==spec['aligned'],'Tile aligned count changed')
    counts,phase,links=orient(pr,sites)
    regression_rows=[r for r in rows(asset_path('assets/real_profiles.tsv.gz')) if int(r['chrom'])==spec['chrom'] and spec['start']<=int(r['cpg_pos1'])<spec['end']]
    regression_checks=0
    if not job.get('fixture',False):
        wanted={int(r['cpg_pos1']) for r in regression_rows}
        observed={(int(z['pos']),IDS[int(z['sample'])],TISSUES[int(z['tissue'])],int(z['source_row'])):z for z in pr[np.isin(pr['pos'],list(wanted))]}
        for expected in regression_rows:
            key=(int(expected['cpg_pos1']),expected['sample'],expected['tissue'],int(expected['source_row']))
            require(key in observed,'Frozen real source row missing');z=observed[key]
            require([int(z[k]) for k in 'abcd']==[int(expected[k]) for k in ['H1_methylated','H1_unmethylated','H2_methylated','H2_unmethylated']],'Frozen integer counts changed')
            require(STATE[int(z['state'])]==expected['orientation_status'] and GT[int(z['GT'])]==expected['anchor_GT'],'Frozen orientation changed')
            if z['state']==0:
                j=np.searchsorted(sites['pos'],z['pos']);i=int(z['sample']);t=int(z['tissue'])
                require(counts[j,i,t].tolist()==[int(expected[k]) for k in ['ALT_methylated','ALT_unmethylated','REF_methylated','REF_unmethylated']],'Frozen ALT/REF counts differ')
            regression_checks+=1
    segments=np.load(asset_path(f"assets/chr{spec['chrom']}.primary.npy"),allow_pickle=False)
    genes=locate(sites['pos'],segments)
    def progress(d):
        guard(out,fixture=job.get('fixture',False));write(out/'worker_progress.json',d)
    filt,results,individual=analyze(counts,sites,genes,progress)
    # Counts stay in the authenticated source; result row points to its site by pos.
    for n,a in [('availability.npy',filt),('results.npy',results),('individual.npy',individual)]:np.save(out/n,a,allow_pickle=False)
    per_pair={}
    for p,name in enumerate(PAIR_NAMES):
        u,v=PAIRS[p];eligible=filt['reason'][:,p]==0
        per_pair[name]=dict(available_by_N=np.bincount(filt['paired_N'][:,p],minlength=11).tolist(),
            tested_by_N=np.bincount(results['N'][results['pair']==p],minlength=11).tolist(),
            reasons={reason:int(np.count_nonzero(filt['reason'][:,p]==code)) for code,reason in enumerate(REASONS)},
            retained_with_mixed_signs=int(np.count_nonzero(eligible&(((filt['positive'][:,u]>0)&(filt['negative'][:,u]>0))|((filt['positive'][:,v]>0)&(filt['negative'][:,v]>0))))),
            retained_with_ties=int(np.count_nonzero(eligible&((filt['ties'][:,u]>0)|(filt['ties'][:,v]>0)))))
        require(per_pair[name]['available_by_N'][2:]==per_pair[name]['tested_by_N'][2:],'Eligible rows lost from analysis')
    for p,before in snaps.items():require(snapshot(p)==before,'Input changed during tile computation: '+p)
    write(out/'source_snapshots.json',snaps)
    return seal(job['parent'],out,job['identity'],dict(chrom=spec['chrom'],start=spec['start'],end=spec['end'],source=spec['path'],
        input_sites=len(sites),input_profiles=len(pr),aligned_profiles=spec['aligned'],passed_pair_rows=len(results),
        per_pair=per_pair,seconds=time.monotonic()-started,sign_filter=False,regression_checks=regression_checks,
        singleton_animal_tests=int(results['singleton_animals'].sum()),individual_tests=int(results['N'].sum())))


def family_job(job):
    out,reused=begin(job['parent'],job['identity'])
    if reused:return str(out)
    guard(out,fixture=job.get('fixture',False));p=job['pair'];r=job['r'];parts=[];fullM=0
    for tile,path in enumerate(job['tiles']):
        source=Path(path);s=read(source/'summary.json');f=np.load(source/'availability.npy',mmap_mode='r',allow_pickle=False)
        z=np.load(source/'results.npy',mmap_mode='r',allow_pickle=False)
        fullM+=int(np.count_nonzero(f['paired_N'][:,p]>=r))
        ix=np.flatnonzero((z['pair']==p)&(z['N']>=r));a=np.zeros(len(ix),dtype=FAMILY)
        a['chrom']=s['chrom'];a['pos']=z['pos'][ix];a['tile']=tile;a['row']=ix;a['N']=z['N'][ix];a['logp']=z['log_pc'][ix,r-2];parts.append(a)
    records=np.concatenate(parts);M=len(records)
    if M:
        require(M==fullM,'Complete aligned eligible family mismatch')
        records['logq_BH']=bh_log(records['logp'])
        harmonic=sum(1./i for i in range(1,fullM+1)) if fullM<1000 else float(digamma(fullM+1)+np.euler_gamma)
        records['logq_BY']=np.minimum(0.,records['logq_BH']+math.log(harmonic))
        records=records[np.lexsort((records['pos'],records['chrom'],records['logp']))]
    require(M==fullM,'Empty result cannot discard an eligible family')
    np.save(out/'family.npy',records,allow_pickle=False)
    return seal(job['parent'],out,job['identity'],dict(pair=PAIR_NAMES[p],pair_index=p,r=r,M=M,
        BH05=int(np.count_nonzero(records['logq_BH']<=math.log(.05))),BH10=int(np.count_nonzero(records['logq_BH']<=math.log(.10))),
        BY05=int(np.count_nonzero(records['logq_BY']<=math.log(.05))),
        min_P_log=float(records['logp'].min()) if M else None,min_BH_log=float(records['logq_BH'].min()) if M else None,
        caution='BH has cross-CpG assumptions. BY uses the same complete family as a dependence sensitivity; neither gives union/max-r FDR.'))
