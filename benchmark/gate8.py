"""GATE 8 (post-hoc, labelled as such) -- cross-crystal reproducibility of the
answer-key site.

Density fit says whether a ligand matches the map in ITS OWN crystal. It cannot
say whether the site is real. The stronger test: does the same ligand land on the
same residues in INDEPENDENT crystals of the same protein?

HRAS is why. 2CE2/XY2 fits at RSCC 0.452, which is suggestive. But XY2 in 2EVW
and 2CL6 shares ONE residue with it -- Jaccard 0.10 and 0.11 -- which is
conclusive, and independent of the density metric entirely.

THRESHOLD FIXED BEFORE LOOKING: median Jaccard >= 0.40 against siblings.
Sensitivity at 0.30 and 0.50 reported. Targets with no sibling crystal cannot be
tested and are reported separately, NOT as failures.
"""
import json, os, urllib.request, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
P=MMCIFParser(QUIET=True)
S = os.environ.get("STRUCTURE_CACHE",
                   os.path.join(os.path.dirname(os.path.dirname(
                       os.path.abspath(__file__))), "structures"))
D=f"{S}/g8"; os.makedirs(D,exist_ok=True)
def get(pid):
    for d in ('crcv2','crcv3','kras','cofcheck','gtp_check','hrascheck','g8'):
        f=f"{S}/{d}/{pid}.cif"
        if os.path.exists(f) and os.path.getsize(f)>500: return f
    f=f"{D}/{pid}.cif"
    try:
        with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.cif",timeout=60) as r:
            open(f,'wb').write(r.read())
    except Exception: return None
    return f if os.path.getsize(f)>500 else None
def target_chains(pid, acc):
    """chains in pid that carry the target accession"""
    try:
        e=json.load(urllib.request.urlopen(f"https://data.rcsb.org/rest/v1/core/entry/{pid}",timeout=20))
        out=[]
        for ne in e['rcsb_entry_container_identifiers']['polymer_entity_ids']:
            d=json.load(urllib.request.urlopen(
              f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{ne}",timeout=20))
            ci=d.get('rcsb_polymer_entity_container_identifiers') or {}
            accs={x.get('database_accession') for x in (ci.get('reference_sequence_identifiers') or [])}
            L=(d.get('entity_poly') or {}).get('rcsb_sample_sequence_length',0) or 0
            if acc in accs and L>=60: out += list(ci.get('auth_asym_ids',[]))
        return sorted(set(out))
    except Exception: return []
def contacts(pid, lig, chains, cut=4.5):
    f=get(pid)
    if not f: return set()
    try: st=P.get_structure(pid,f)[0]
    except Exception: return set()
    la=[a for c in st for r in c if r.get_resname()==lig for a in r]
    if not la: return set()
    L=np.array([a.coord for a in la],float); h=set()
    for ch in st:
        if chains and ch.id not in chains: continue
        for r in ch:
            if not is_aa(r,standard=True): continue
            A=np.array([a.coord for a in r],float)
            if np.sqrt(((A[:,None,:]-L[None,:,:])**2).sum(2)).min()<=cut: h.add(r.id[1])
    return h
def siblings(acc, lig, exclude):
    """other entries of this accession containing this ligand"""
    b={"query":{"type":"group","logical_operator":"and","nodes":[
        {"type":"terminal","service":"text","parameters":{
          "attribute":"rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
          "operator":"exact_match","value":acc}},
        {"type":"terminal","service":"text","parameters":{
          "attribute":"rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
          "operator":"exact_match","value":lig}}]},
       "return_type":"entry","request_options":{"paginate":{"start":0,"rows":50}}}
    r=urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query",
      data=json.dumps(b).encode(),headers={"Content-Type":"application/json"})
    try:
        ids=[h['identifier'] for h in json.loads(urllib.request.urlopen(r,timeout=45).read()).get('result_set',[])]
    except Exception: return []
    return [i for i in ids if i!=exclude]

cm=json.load(open(f"{S}/chain_map.json"))
g7={r['symbol']:r for r in json.load(open(f"{S}/gate7.json"))}
surv=[s for s,v in g7.items() if v['verdict']!='EXCLUDE']
print(f"testing the {len(surv)} targets that survived gate 7, plus HRAS as the known positive control\n")
print(f"  {'target':<9}{'lig':>8}{'sibs':>6}{'median J':>10}{'min J':>8}{'>=0.40':>8}  siblings")
rows=[]
for sym in surv+['HRAS']:
    v=cm[sym]; acc,hid,lg=v['acc'],v['holo'],v['lig']
    key=contacts(hid,lg,set(v['target_chains']))
    if not key:
        print(f"  {sym:<9}{lg:>8}{'--':>6}{'--':>10}{'--':>8}{'--':>8}  key site empty"); continue
    sibs=siblings(acc,lg,hid)
    js=[]
    for sp in sibs[:8]:
        ch=target_chains(sp,acc)
        c=contacts(sp,lg,set(ch) if ch else None)
        if len(c)<3: continue
        j=len(key&c)/len(key|c)
        js.append((sp,j,len(c)))
    if not js:
        print(f"  {sym:<9}{lg:>8}{len(sibs):>6}{'--':>10}{'--':>8}{'n/a':>8}  no testable sibling")
        rows.append(dict(symbol=sym,n=0,med=None,mn=None,pass40=None)); continue
    vals=[j for _,j,_ in js]
    med=float(np.median(vals)); mn=float(min(vals))
    ok='yes' if med>=0.40 else 'NO'
    print(f"  {sym:<9}{lg:>8}{len(js):>6}{med:>10.2f}{mn:>8.2f}{ok:>8}  "
          +", ".join(f"{p}:{j:.2f}" for p,j,_ in js[:4]))
    rows.append(dict(symbol=sym,n=len(js),med=med,mn=mn,pass40=med>=0.40))
json.dump(rows,open(f"{S}/gate8.json","w"),indent=1)
t=[r for r in rows if r['med'] is not None]
print(f"\ntestable: {len(t)} of {len(rows)}")
for thr in (0.30,0.40,0.50):
    f=[r['symbol'] for r in t if r['med']<thr]
    print(f"  threshold {thr:.2f}: {len(f)} fail  {f}")
