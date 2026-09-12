"""Recompute the cofactor gate from ATOM RECORDS, not the metadata field.
How many targets in either cohort have apo members in inconsistent states?"""
import json, os, urllib.request
COF={'GDP','GTP','GNP','GSP','GCP','G2P','GBS','ATP','ADP','AMP','ANP','AGS','ACP','APC',
     'NAD','NAI','NAP','NDP','FAD','FMN','SAM','SAH','COA','ACO','PLP','HEM','HEC','U5P','UMP'}
os.makedirs('cofcheck',exist_ok=True)
def het(pid):
    for d in ('cofcheck','crcv2','crcv3','gtp_check','kras'):
        f=f"{d}/{pid}.cif"
        if os.path.exists(f) and os.path.getsize(f)>500: break
    else:
        f=f"cofcheck/{pid}.cif"
        try:
            with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.cif",timeout=60) as r:
                open(f,'wb').write(r.read())
        except Exception: return None
    s=set()
    for line in open(f,errors='ignore'):
        if line.startswith('HETATM'):
            p=line.split()
            if len(p)>5 and p[5] in COF: s.add(p[5])
    return s
rows=[]
for tag,fn in (("v2","crc_v2_cohort.json"),("v3","crc_v3_cohort.json")):
    if not os.path.exists(fn): continue
    for t in json.load(open(fn)):
        apo_states={}
        for apo_id,dep,res,title in t['apo'][:12]:
            h=het(apo_id)
            if h is None: continue
            apo_states[apo_id]=frozenset(h)
        if not apo_states: continue
        distinct=set(apo_states.values())
        holo_states=set()
        for hid,drugs,dep,res in t['holo'][:6]:
            h=het(hid)
            if h is not None: holo_states.add(frozenset(h))
        bad = len(distinct)>1
        rows.append((tag,t['symbol'],len(apo_states),len(distinct),
                     sorted({x for s in distinct for x in s}),
                     sorted({x for s in holo_states for x in s}), bad))
        print(f"  {tag} {t['symbol']:<9} apo {len(apo_states):>2} states {len(distinct)}"
              f"  apo cofactors {sorted({x for s in distinct for x in s}) or '[]'}"
              f"  holo {sorted({x for s in holo_states for x in s}) or '[]'}"
              f"{'   <-- MIXED APO STATES' if bad else ''}",flush=True)
json.dump([{'cohort':r[0],'symbol':r[1],'n_apo':r[2],'n_states':r[3],
            'apo_cof':r[4],'holo_cof':r[5],'mixed':r[6]} for r in rows],
          open('cofactor_audit.json','w'),indent=1)
mixed=[r for r in rows if r[6]]
print(f"\ntargets with apo members in MIXED cofactor states: {len(mixed)} of {len(rows)}")
for r in mixed: print(f"   {r[0]} {r[1]}")
