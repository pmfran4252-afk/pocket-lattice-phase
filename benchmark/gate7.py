"""GATE 7 -- answer-key ligand density quality.

Reports BOTH the best and worst instance of the ligand in each entry. The
contacts used as the answer key were collected from every instance on the target
chains, so the WORST instance is the honest figure; taking the best flatters it.

  EXCLUDE  worst RSCC < 0.80, or no validation score (unverifiable)
  FLAG     occupancy < 0.70, completeness < 0.90, or the depositor did not mark
           the ligand as the subject of investigation
"""
import json, urllib.request
keys=json.load(open('chain_map.json'))
def instances(hid, lg):
    try:
        e=json.load(urllib.request.urlopen(f"https://data.rcsb.org/rest/v1/core/entry/{hid}",timeout=20))
        nps=(e.get('rcsb_entry_container_identifiers') or {}).get('non_polymer_entity_ids') or []
        res=(e.get('rcsb_entry_info') or {}).get('resolution_combined',[None])[0]
    except Exception: return None,[]
    out=[]
    for ne in nps:
        try:
            ent=json.load(urllib.request.urlopen(
              f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{hid}/{ne}",timeout=20))
        except Exception: continue
        ci=ent.get('rcsb_nonpolymer_entity_container_identifiers') or {}
        if ci.get('nonpolymer_comp_id')!=lg: continue
        for asym in ci.get('asym_ids') or []:
            try:
                inst=json.load(urllib.request.urlopen(
                  f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity_instance/{hid}/{asym}",timeout=20))
            except Exception: continue
            for s in inst.get('rcsb_nonpolymer_instance_validation_score') or []:
                if s.get('RSCC') is None: continue
                out.append(dict(asym=asym,rscc=s['RSCC'],rsr=s.get('RSR'),
                                occ=s.get('average_occupancy'),comp=s.get('completeness'),
                                subj=s.get('is_subject_of_investigation')))
    return res,out
print(f"  {'target':<9}{'lig':>8}{'res':>6}{'n':>3}{'worst':>8}{'best':>7}{'occ':>6}{'comp':>6} subj  GATE 7")
rows=[]
for sym,v in keys.items():
    res,ins=instances(v['holo'],v['lig'])
    if not ins:
        print(f"  {sym:<9}{v['lig']:>8}{(res or 0):>6.2f}{0:>3}{'--':>8}{'--':>7}{'--':>6}{'--':>6}  -   EXCLUDE (no score)")
        rows.append(dict(symbol=sym,worst=None,verdict='EXCLUDE')); continue
    w=min(i['rscc'] for i in ins); b=max(i['rscc'] for i in ins)
    wi=min(ins,key=lambda i:i['rscc'])
    flags=[]
    if (wi['occ'] or 1)<0.70: flags.append('low occupancy')
    if (wi['comp'] or 1)<0.90: flags.append('incomplete')
    if wi['subj']=='N': flags.append('not subject of investigation')
    verdict='EXCLUDE' if w<0.80 else ('PASS+flag' if flags else 'PASS')
    print(f"  {sym:<9}{v['lig']:>8}{(res or 0):>6.2f}{len(ins):>3}{w:>8.3f}{b:>7.3f}"
          f"{(wi['occ'] if wi['occ'] is not None else -1):>6.2f}"
          f"{(wi['comp'] if wi['comp'] is not None else -1):>6.2f}"
          f"  {wi['subj'] or '-'}   {verdict}"
          + (f"  [{'; '.join(flags)}]" if flags else ""))
    rows.append(dict(symbol=sym,worst=w,best=b,verdict=verdict,flags=flags))
json.dump(rows,open('gate7.json','w'),indent=1)
ex=[r['symbol'] for r in rows if r['verdict']=='EXCLUDE']
fl=[r['symbol'] for r in rows if r['verdict']=='PASS+flag']
print(f"\nEXCLUDED by gate 7 ({len(ex)}): {', '.join(ex)}")
print(f"FLAGGED   by gate 7 ({len(fl)}): {', '.join(fl)}")
