"""Preserve observations needed for exact SNP-GQ reassignment without BAM rereading."""
from __future__ import annotations
import bisect
from collections import Counter
import gzip
import json
import math
from pathlib import Path
from common import read_json, write_json, write_table, rows, sha, CALLS

VARIANT_FIELDS = ['source','chrom','pos1','ref','alt','QUAL','FILTER','GT','GQ','DP','AD','PL','PS','eligible_link']


def variant_row(role, fields, sample_col, links):
    fmt = dict(zip(fields[8].split(':'),fields[sample_col].split(':')))
    return dict(source=role,chrom=fields[0],pos1=int(fields[1]),ref=fields[3],alt=fields[4],
                QUAL=fields[5],FILTER=fields[6],**{k:fmt.get(k,'.') for k in ('GT','GQ','DP','AD','PL','PS')},
                eligible_link=int(int(fields[1])-1 in links))


def numeric_gq(value):
    if value in (None,'','.'): return None
    try: n = float(value)
    except (TypeError,ValueError): return None
    return n if math.isfinite(n) and n >= 0 else None


def genotype_alleles(row):
    alleles = [row['ref'].upper()] + row['alt'].upper().split(',')
    try:
        gt = row['GT'].replace('|','/').split('/')
        if len(gt) != 2 or any(not x.isdigit() for x in gt): return None
        return tuple(alleles[int(x)] for x in gt)
    except IndexError: return None


def save_catalog(directory, variants, links, contexts):
    variants = sorted(variants,key=lambda r:(r['source'],r['pos1'],r['ref'],r['alt']))
    index = {}
    phased_by_pos = {}
    for r in variants:
        key = (r['source'],r['pos1'],r['ref'],r['alt'])
        if key in index: raise ValueError('Duplicate VCF allele key in evidence catalog')
        index[key] = r
        if r['source']=='phased_vcf': phased_by_pos.setdefault(r['pos1'],[]).append(r)
    catalog = {}
    for p,info in sorted(links.items()):
        candidates = [r for r in phased_by_pos.get(p+1,[])
                      if genotype_alleles(r)==(info.h1,info.h2) and r['PS']==info.ps]
        if len(candidates) != 1: raise ValueError('Eligible link has no unique phased VCF provenance')
        phase = candidates[0]
        raw = index.get(('genotype_vcf',p+1,phase['ref'],phase['alt']))
        agrees = bool(raw and genotype_alleles(raw) and sorted(genotype_alleles(raw))==sorted((info.h1,info.h2)))
        pgq, rgq = numeric_gq(phase['GQ']), numeric_gq(raw['GQ']) if raw else None
        source = 'phased_vcf' if pgq is not None else 'genotype_vcf_matched_fallback' if agrees and rgq is not None else 'MISSING'
        catalog[str(p)] = dict(pos1=p+1,ps=info.ps,h1=info.h1,h2=info.h2,
            gq=pgq if pgq is not None else rgq if agrees else None,gq_source=source,
            phased_gq=pgq,genotype_gq=rgq,genotype_alleles_agree=agrees,
            phase_record=phase,genotype_record=raw)
    write_table(directory/'variant_records.tsv.gz',VARIANT_FIELDS,([r[k] for k in VARIANT_FIELDS] for r in variants))
    write_json(directory/'link_catalog.json',{'schema_version':'0.2.0','links':catalog,
        'contexts':{str(p):sorted(v) for p,v in sorted(contexts.items())},
        'gq_policy':'Use numeric phased GQ; only when missing, use an exact REF/ALT and genotype-matched DeepVariant GQ.',
        'limits':'GQ is not phase confidence. All existing link eligibility and CpG disruption rules are unchanged.'})
    return catalog


def atoms(record, positions, shift=0, methylation=False):
    from profile import CORE
    output = []
    quality = getattr(record,'audit_quality','*')
    for qstart,start,length in CORE.aligned_blocks(record):
        for p in positions[bisect.bisect_left(positions,start-shift):bisect.bisect_left(positions,start+length-shift)]:
            q = qstart+p+shift-start
            base = record.seq[q].upper()
            bq = ord(quality[q])-33 if quality != '*' and q < len(quality) else None
            item = [p,base,bq,q]
            if methylation: item.append(record.tags['XM'][q])
            output.append(item)
    return output


def write_fragment(handle, serial, records, positions, known_positions, result):
    reads = []
    for r in records:
        reads.append(dict(flag=r.flag,start0=r.pos0,mate_start0=r.mate_pos0,mate_chrom=r.mate_rname,
            mapq=r.mapq,cigar=r.cigar,xg=r.tags['XG'],xr=r.tags.get('XR'),
            snps=atoms(r,known_positions),cpgs=atoms(r,positions,int(r.tags['XG']=='GA'),True)))
    item = dict(schema_version='0.2.0',fragment_id=serial,qname=records[0].qname,reads=reads,
        baseline={k:result[k] for k in ('state','ps','hap')},
        incomplete_fetch_context=result['state']=='INCOMPLETE_FETCH_CONTEXT')
    handle.write(json.dumps(item,separators=(',',':'),ensure_ascii=True)+'\n')


def selected_links(catalog, min_gq):
    from profile import CORE
    if min_gq is not None and (not math.isfinite(min_gq) or min_gq < 0):
        raise ValueError('Minimum GQ must be finite and nonnegative')
    return {int(p):CORE.SNPInfo(int(p),v['ps'],v['h1'],v['h2']) for p,v in catalog['links'].items()
            if min_gq is None or v['gq'] is not None and v['gq'] >= min_gq}


def classify(reads, links, contexts, incomplete=False):
    """Reconstruct legacy conversion, mate-conflict and exact-PS rules from atoms."""
    from profile import CORE, engine
    known = {a[0] for r in reads for a in r['snps']}
    eligible = known.intersection(links)
    calls, conflicts = {}, set()
    for r in reads:
        for p,base,_bq,_q in r['snps']:
            if p not in links or base not in 'ACGT': continue
            info = links[p]
            support = engine.allele_support(base,info.h1,info.h2,r['xg'])
            normalized = (info.h1,info.h2)[support] if support is not None else base
            CORE.merge_call(calls,conflicts,p,normalized)
    states,votes = CORE.infer_ps_states(CORE.Fragment('replay',{},calls),links,1)
    if incomplete: return 'INCOMPLETE_FETCH_CONTEXT','.','.'
    if len(votes)>1: return 'MULTIPLE_PHASE_SETS','.','.'
    if len(states)==1:
        ps,hap = next(iter(states.items()))
        return 'ASSIGNED',ps,hap
    state = ('CONFLICTING_HAPLOTYPES' if votes else 'MATE_SNP_CONFLICT' if conflicts else
             'NO_KNOWN_SNV_OVERLAP' if not known else 'NO_ELIGIBLE_PHASED_LINK' if not eligible else
             'NO_USABLE_LINK_BASE' if not calls else 'NO_ALLELE_SUPPORT')
    return state,'.','.'


def methylation(read):
    return {a[0]:('METHYLATED' if a[4]=='Z' else 'UNMETHYLATED' if a[4]=='z' else 'UNCALLABLE') for a in read['cpgs']}


def reconstruct(directory, min_gq=None, verify_baseline=False):
    from profile import merge_cpg
    directory = Path(directory)
    catalog = read_json(directory/'link_catalog.json')
    links = selected_links(catalog,min_gq)
    sites, phases, transitions = Counter(), Counter(), Counter()
    rescued = Counter()
    fragments = 0
    with gzip.open(directory/'fragment_evidence.jsonl.gz','rt',encoding='utf-8') as f:
        for line in f:
            e = json.loads(line)
            fragments += 1
            if e['schema_version'] != '0.2.0' or e['fragment_id'] != fragments:
                raise ValueError('Evidence schema or sequential fragment identity differs')
            observations = [methylation(r) for r in e['reads']]
            read_states = [classify([r],links,catalog['contexts'])[0] for r in e['reads']]
            for unit,rs,obs,incomplete in (
                [('READ',[r],o,False) for r,o in zip(e['reads'],observations)] +
                [('FRAGMENT',e['reads'],merge_cpg(observations),e['incomplete_fetch_context'])]):
                if not obs: continue
                state,ps,hap = classify(rs,links,catalog['contexts'],incomplete)
                if unit=='FRAGMENT':
                    old = e['baseline']
                    if verify_baseline and (state,ps,hap)!=(old['state'],old['ps'],old['hap']):
                        raise ValueError('Saved raw SNP observations do not reproduce baseline assignment')
                    transitions[old['state'],state] += 1
                for pos,call in obs.items():
                    sites[pos,unit,state,call] += 1
                    if state=='ASSIGNED': phases[pos,ps,unit,hap,call] += 1
                    if unit=='FRAGMENT' and state=='ASSIGNED' and all(
                        s!='ASSIGNED' for s,o in zip(read_states,observations) if pos in o):
                        rescued[pos] += 1
    return sites,phases,transitions,fragments,len(links),rescued


def verify_evidence(directory):
    """Check all read/fragment CpG and haplotype counts against raw saved atoms."""
    directory = Path(directory)
    sites,phases,transitions,nfragments,nlinks,rescued = reconstruct(directory,verify_baseline=True)
    expected_sites, expected_phases = Counter(), Counter()
    from common import STATES
    for r in rows(directory/'cpg_assignment.tsv.gz'):
        if rescued[int(r['cpg_pos1'])-1] != int(r['fragment_assigned_using_mate_only']):
            raise ValueError('Mate-only assignment replay differs')
        for u in ('READ','FRAGMENT'):
            for s in STATES:
                for c in CALLS: expected_sites[int(r['cpg_pos1'])-1,u,s,c] += int(r[f'{u}_{s}_{c}'])
    for r in rows(directory/'cpg_haplotype_by_ps.tsv.gz'):
        for u in ('READ','FRAGMENT'):
            for h in ('H1','H2'):
                for c in CALLS: expected_phases[int(r['cpg_pos1'])-1,r['ps'],u,h,c] += int(r[f'{u}_{h}_{c}'])
    if sites != expected_sites or phases != expected_phases:
        raise ValueError('Evidence replay does not reconstruct all baseline CpG and phase counts')
    result = dict(status='PASS_UNFILTERED_EVIDENCE_REPLAY',errors=[],fragment_rows=nfragments,eligible_links=nlinks,
        source_sha256={p:sha(directory/p) for p in ('fragment_evidence.jsonl.gz','link_catalog.json','variant_records.tsv.gz')},
        genotype_gq_filter=None,raw_BAM_reread=False,
        limits='Saved evidence and baseline count equivalence; not biological assignment or statistical validation.')
    write_json(directory/'evidence_validation.json',result)
    return result
