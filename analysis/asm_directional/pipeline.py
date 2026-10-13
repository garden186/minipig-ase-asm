"""Read authenticated stage-1 indexes; write tiled alignment and directional PC."""
import csv,gzip,collections,itertools,math,time
import numpy as np
from common import *
from anchors import assign,variant_keys
from stats import directional_logp,partial_conjunction,bh_log,adjusted_direction_logq,display,LOG_MIN

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

def slice_positions(a,start,end):return a[np.searchsorted(a['pos'],start):np.searchsorted(a['pos'],end)]

def opportunity_pairs(a,positions,mapping):
    if not len(a):return np.empty(0,dtype='<u8')
    j=np.searchsorted(positions,a['pos']);good=j<len(positions);good[good]&=positions[j[good]]==a['pos'][good]
    phases=mapping[a['ps']];good&=phases>0
    return np.unique((j[good].astype(np.uint64)<<32)|phases[good].astype(np.uint64))

def subset_expected(config,c,start,end):
    p=Path(config['locations']['pilot'][f'chr{c}'])
    sites={int(r['cpg_pos1']):r for r in rows(p/'pilot_sites.tsv.gz') if start<=int(r['cpg_pos1'])<end}
    prof={(int(r['cpg_pos1']),r['sample'],r['tissue'],r['ps']):r for r in rows(p/'pilot_profiles.tsv.gz') if start<=int(r['cpg_pos1'])<end}
    return sites,prof

def align_tile(job):
    out,reused=begin(job['parent'],job['identity'])
    if reused:return str(out)
    started=time.monotonic();check_resources=lambda:guard(out,fixture=job.get('fixture',False));check_resources();config=job['config'];loc=config['locations'];c=job['chrom'];start=job['start'];end=job['end']
    positions=np.load(Path(loc['pilot'][f'chr{c}'])/'eligible_positions.npy',mmap_mode='r')
    positions=positions[np.searchsorted(positions,start):np.searchsorted(positions,end)].copy()
    n=len(positions);units={};cache_ps={};dicts=[];maps=[];phase=[];allvars=[]
    for s in IDS:
        vp=Path(loc['phase_index'][f'{s}-chr{c}']);summ=read(vp/'summary.json');d=read(summ['phase_dictionary'])
        dicts.append(d);maps.append({s:i for i,s in enumerate(d)});phase.append(np.load(vp/'eligible_by_ps.npy',mmap_mode='r'))
        allvars.append(np.load(vp/'variants.npy',mmap_mode='r'))
    available=[[] for _ in IDS];tested=[[] for _ in IDS]
    for spec in config['source_inputs']['units']:
        s,t=spec['sample'],spec['tissue'];i=IDS.index(s);unit=s+'-'+t;up=Path(loc['units'][unit]);cp=Path(spec['cache']['path'])
        d=read(up/'phase_dictionary.json');mapping=np.array([maps[i].get(x,0) for x in d],dtype='<u4')
        ar=slice_positions(np.load(up/f'chr{c}.available.npy',mmap_mode='r'),start,end)
        available[i].append(opportunity_pairs(ar,positions,mapping))
        d=read(cp/'phase_dictionary.json');cache_ps[(s,t)]=d;mapping=np.array([maps[i].get(x,0) for x in d],dtype='<u4')
        ar=slice_positions(np.load(cp/f'chr{c}.npy',mmap_mode='r'),start,end)
        units[(s,t)]=(ar,mapping);tested[i].append(opportunity_pairs(ar,positions,mapping))
    av=[np.unique(np.concatenate(x)) for x in available];te=[np.unique(np.concatenate(x)) for x in tested]
    ak,nt,na,nc,perf=assign(positions,av,te,phase,check_resources)
    sites=np.zeros(n,dtype=SITE);sites['pos']=positions;sites['anchor']=ak;sites['anchor_testable']=nt;sites['anchor_available']=na;sites['candidate_snps']=nc
    profiles=np.empty(sum(len(v[0]) for v in units.values()),dtype=PROFILE);cursor=0;states=collections.Counter();missing=0
    matrices=np.full((n,3,10,2),np.nan,dtype='<f8');animal_variants=[]
    for i,v in enumerate(allvars):
        keys=variant_keys(v);q=np.searchsorted(keys,ak);found=q<len(keys);found[found]&=keys[q[found]]==ak[found]
        chosen=np.zeros(n,dtype=v.dtype);chosen[found]=v[q[found]];animal_variants.append((chosen,found))
    for (s,t),(ar,mapping) in units.items():
        check_resources();i=IDS.index(s);ti=TISSUES.index(t);j=np.searchsorted(positions,ar['pos']);v,found=animal_variants[i];v=v[j];f=found[j]
        length=len(ar);pr=profiles[cursor:cursor+length];cursor+=length
        pr['pos']=ar['pos'];pr['sample']=i;pr['tissue']=ti;pr['ps']=ar['ps'];pr['source_row']=ar['source_row']
        for field in ['a','b','c','d']:pr[field]=ar[field]
        pr['original_p']=ar['p'];pr['original_q']=ar['q_all'];pr['GT']=v['gt'];pr['anchor_ps']=v['ps'];pr['alt_hap']=v['alt_hap'];pr['logp']=np.nan
        st=np.full(length,STATE.index('ALIGNED'),dtype='u1')
        for code,name in enumerate(PREP_STATUS[1:],1):st[f&(v['status']==code)]=STATE.index(name)
        st[f&(v['status']==0)&(v['ps']!=mapping[ar['ps']])]=STATE.index('DIFFERENT_EXACT_PS')
        st[~f]=STATE.index('NO_MATCHING_VARIANT_RECORD');st[ak[j]==0]=STATE.index('NO_ELIGIBLE_SITE_ANCHOR');pr['state']=st
        good=st==0;gj=j[good];require(len(np.unique(gj))==len(gj),'Multiple aligned PS per animal/tissue/CpG')
        sites['N_before'][np.unique(j),ti]+=1;sites['N'][gj,ti]+=1
        counts=np.column_stack([ar[k][good] for k in ['a','b','c','d']]);swap=v['alt_hap'][good]==2
        counts[swap]=counts[swap][:,[2,3,0,1]];lp,num=directional_logp(counts);missing+=num
        pr['logp'][good]=lp;matrices[gj,ti,i,:]=lp
        states.update({STATE[int(k)]:int(value) for k,value in zip(*np.unique(st,return_counts=True))})
    pc=np.full((n,3,2,9),np.nan,dtype='<f8')
    for ti in range(3):
        for N in range(2,11):
            ix=np.flatnonzero(sites['N'][:,ti]==N)
            if not len(ix):continue
            values=matrices[ix,ti];mask=np.isfinite(values[:,:,0]);require(np.all(mask.sum(axis=1)==N),'N differs from actual two-direction input')
            values=values[mask].reshape(len(ix),N,2)
            for direction in range(2):
                for r in range(2,N+1):pc[ix,ti,direction,r-2]=partial_conjunction(values[:,:,direction],r)
    require(np.array_equal(np.any(np.isfinite(matrices[:,:,:,0]),axis=1).sum(axis=1),nt),'Anchor union N differs from actual aligned profiles')
    # Full table/anchor regression across every returned stage-1 pilot coordinate.
    es,ep=subset_expected(config,c,start,end);checked=0;two_sided_samples=[]
    from scipy.stats import fisher_exact
    two_cache={}
    for pos,r in es.items():
        ix=int(np.searchsorted(positions,pos));require(ix<n and positions[ix]==pos,'Pilot coordinate disappeared');site=sites[ix];key=int(site['anchor'])
        expected=0 if r['anchor_pos1']=='.' else int(r['anchor_pos1'])*16+'ACGT'.index(r['anchor_ref'])*4+'ACGT'.index(r['anchor_alt'])
        require(key==expected and int(site['anchor_testable'])==int(r['anchor_testable_animals']) and int(site['anchor_available'])==int(r['anchor_available_animals']) and int(site['candidate_snps'])==int(r['candidates']),'Optimized anchor differs from pilot')
        for ti,t in enumerate(TISSUES):require(int(site['N_before'][ti])==int(r['N_before_'+t]) and int(site['N'][ti])==int(r['N_aligned_'+t]),'Pilot tissue N differs')
    for pr in profiles[np.isin(profiles['pos'],list(es))]:
        s=IDS[int(pr['sample'])];t=TISSUES[int(pr['tissue'])];ps=cache_ps[(s,t)][int(pr['ps'])];key=(int(pr['pos']),s,t,ps)
        require(key in ep,'Unexpected pilot profile');r=ep.pop(key);checked+=1
        require(STATE[int(pr['state'])]==r['orientation_status'],'Pilot orientation differs')
        require([int(pr[k]) for k in ['a','b','c','d']]==[int(r[k]) for k in ['H1_methylated','H1_unmethylated','H2_methylated','H2_unmethylated']],'Pilot original counts differ')
        require(float(pr['original_p'])==float(r['original_fisher_two_sided_p']) and float(pr['original_q'])==float(r['original_q_all30']),'Pilot original P/q differs')
        if pr['state']==0:
            require(dicts[int(pr['sample'])][int(pr['anchor_ps'])]==ps,'Pilot exact PS mismatch')
            counts=tuple(int(pr[k]) for k in ['a','b','c','d']);oriented=counts if pr['alt_hap']==1 else counts[2:]+counts[:2]
            if counts not in two_cache:two_cache[counts]=float(fisher_exact(np.array(counts).reshape(2,2),alternative='two-sided').pvalue)
            require(math.isclose(two_cache[counts],float(pr['original_p']),rel_tol=2e-8,abs_tol=0),'Original two-sided Fisher replay failure')
            if oriented not in two_cache:two_cache[oriented]=float(fisher_exact(np.array(oriented).reshape(2,2),alternative='two-sided').pvalue)
            require(math.isclose(two_cache[oriented],float(pr['original_p']),rel_tol=2e-8,abs_tol=0),'Allele swap changed two-sided Fisher')
            require(float(r['ALT_methylation'])==oriented[0]/(oriented[0]+oriented[1]) and float(r['REF_methylation'])==oriented[2]/(oriented[2]+oriented[3]),'Pilot fractions differ')
    require(not ep,'Missing pilot profiles')
    order=np.lexsort((profiles['ps'],profiles['sample'],profiles['tissue'],profiles['pos']));profiles=profiles[order]
    np.save(out/'sites.npy',sites,allow_pickle=False);np.save(out/'profiles.npy',profiles,allow_pickle=False);np.save(out/'pc_logp.npy',pc,allow_pickle=False)
    return seal(job['parent'],out,job['identity'],dict(chrom=c,start=start,end=end,sites=n,profiles=len(profiles),states=dict(states),
        N_histogram={t:np.bincount(sites['N'][:,i],minlength=11).tolist() for i,t in enumerate(TISSUES)},
        family_sizes={t:{str(r):int((sites['N'][:,i]>=r).sum()) for r in range(2,11)} for i,t in enumerate(TISSUES)},
        pilot_sites_checked=len(es),pilot_profiles_checked=checked,two_sided_unique_tables_checked=len(two_cache),log_tail_fallback=missing,anchor_performance=perf,seconds=time.monotonic()-started))

def family_job(job):
    out,reused=begin(job['parent'],job['identity'])
    if reused:return str(out)
    start=time.monotonic();guard(out,fixture=job.get('fixture',False));ti=job['tissue'];direction=job['direction'];r=job['r'];chunks=[]
    for path in job['tiles']:
        p=Path(path);s=np.load(p/'sites.npy',mmap_mode='r');pc=np.load(p/'pc_logp.npy',mmap_mode='r');ix=np.flatnonzero(s['N'][:,ti]>=r)
        ar=np.empty(len(ix),dtype=HIT);ar['chrom']=read(p/'summary.json')['chrom'];ar['pos']=s['pos'][ix];ar['N']=s['N'][ix,ti];ar['logp']=pc[ix,ti,direction,r-2];chunks.append(ar)
    full=np.concatenate(chunks) if chunks else np.empty(0,dtype=HIT)
    order=np.lexsort((full['pos'],full['chrom']));full=full[order]
    keys=full['chrom'].astype(np.uint64)*2**32+full['pos'];require(not len(keys) or np.all(keys[1:]>keys[:-1]),'Duplicate full-family hypothesis')
    full['logq_bh']=bh_log(full['logp']);full['logq']=adjusted_direction_logq(full['logq_bh'])
    np.save(out/'full_family.npy',full,allow_pickle=False)
    hit=full[full['logq']<=math.log(.1)];np.save(out/'hits_BH10.npy',hit,allow_pickle=False)
    keep=np.argsort(full['logp'],kind='stable')[:25];np.save(out/'top25.npy',full[keep],allow_pickle=False)
    records=[]
    for a in hit:
        records.append(dict(chrom=int(a['chrom']),cpg_pos1=int(a['pos']),tissue=TISSUES[ti],direction=['ALT_gt_REF','REF_gt_ALT'][direction],r=r,N=int(a['N']),
            log_PC_P=float(a['logp']),PC_P=float(display(a['logp'])),log_BH_q=float(a['logq_bh']),BH_q=float(display(a['logq_bh'])),
            log_direction_adjusted_q=float(a['logq']),direction_adjusted_q=float(display(a['logq'])),BH05=int(a['logq']<=math.log(.05)),
            display_underflow=int(a['logp']<LOG_MIN),family_size=len(full)))
    fields=['chrom','cpg_pos1','tissue','direction','r','N','log_PC_P','PC_P','log_BH_q','BH_q','log_direction_adjusted_q','direction_adjusted_q','BH05','display_underflow','family_size']
    table(out/'hits_BH10.tsv.gz',records,fields)
    return seal(job['parent'],out,job['identity'],dict(tissue=TISSUES[ti],direction=['ALT_gt_REF','REF_gt_ALT'][direction],r=r,hypotheses=len(full),BH05=int((full['logq']<=math.log(.05)).sum()),
        BH10=len(hit),minimum_logp=float(full['logp'].min()) if len(full) else None,minimum_adjusted_q=float(display(full['logq'].min())) if len(full) else None,
        display_underflow=int((full['logp']<LOG_MIN).sum()),seconds=time.monotonic()-start,
        correction='Separate directional full-background BH, then min(1,2*qBH). No P-min selection, subset BH or union-across-r guarantee.'))
