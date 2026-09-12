"""Resolve which chains in each holo entry actually ARE the target protein,
from the RCSB polymer-entity records. The answer key must come only from those
chains -- contacts on a partner protein are not the target's binding site."""
import json, os, urllib.request
occ=json.load(open('apo_occupancy.json'))
v2={t['symbol']:t for t in json.load(open('crc_v2_cohort.json'))}
v3={t['symbol']:t for t in json.load(open('crc_v3_cohort.json'))}
coh={**v2,**v3}
out={}
for sym,v in occ.items():
    acc=coh[sym]['acc']
    hid,lg=v['site'].split(':')
    try:
        e=json.load(urllib.request.urlopen(f"https://data.rcsb.org/rest/v1/core/entry/{hid}",timeout=25))
        ids=e['rcsb_entry_container_identifiers']['polymer_entity_ids']
    except Exception:
        print(f"  {sym}: entry lookup failed"); continue
    good=[]
    for pe in ids:
        try:
            d=json.load(urllib.request.urlopen(
                f"https://data.rcsb.org/rest/v1/core/polymer_entity/{hid}/{pe}",timeout=20))
        except Exception: continue
        ci=d.get('rcsb_polymer_entity_container_identifiers') or {}
        accs={x.get('database_accession') for x in (ci.get('reference_sequence_identifiers') or [])}
        L=(d.get('entity_poly') or {}).get('rcsb_sample_sequence_length',0) or 0
        if acc in accs and L>=60:
            good += list(ci.get('auth_asym_ids',[]))
    out[sym]=dict(acc=acc,holo=hid,lig=lg,target_chains=sorted(set(good)))
    print(f"  {sym:<9} {hid}:{lg:<7} acc {acc:<8} target chains {sorted(set(good)) or 'NONE FOUND'}",flush=True)
json.dump(out,open('chain_map.json','w'),indent=1)
none=[s for s,v in out.items() if not v['target_chains']]
print(f"\nentries where no chain maps to the target accession: {len(none)}  {none}")
