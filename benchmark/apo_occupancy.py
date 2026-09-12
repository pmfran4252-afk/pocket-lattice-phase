"""Refined gate. A pocket is OCCUPIED if a substantial fraction of its residues
are in contact with a non-water heteroatom -- not merely if one heteroatom is
nearby. GDP sits adjacent to KRAS switch-II without occupying it; an ATP
analogue in a kinase ATP site occupies the drug site outright.

Threshold fixed before looking: > 50% of site residues contacted at 4.5 A.
"""
import json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
P=MMCIFParser(QUIET=True)
def path(pid):
    for d in ('crcv2','crcv3','cofcheck','gtp_check','kras'):
        f=f"{d}/{pid}.cif"
        if os.path.exists(f) and os.path.getsize(f)>500: return f
    return None
def contacts(pid,lig,cut=4.5):
    f=path(pid)
    if not f: return set()
    try: st=P.get_structure(pid,f)[0]
    except Exception: return set()
    la=[a for c in st for r in c if r.get_resname()==lig for a in r]
    if not la: return set()
    L=np.array([a.coord for a in la],float); h=set()
    for ch in st:
        for r in ch:
            if not is_aa(r,standard=True): continue
            A=np.array([a.coord for a in r],float)
            if np.sqrt(((A[:,None,:]-L[None,:,:])**2).sum(2)).min()<=cut: h.add(r.id[1])
    return h
def occupancy(pid, site, cut=4.5):
    """Fraction of site residues contacting a non-water heteroatom."""
    f=path(pid)
    if not f: return None
    try: st=P.get_structure(pid,f)[0]
    except Exception: return None
    H=[]
    for ch in st:
        for r in ch:
            if not is_aa(r,standard=True) and r.get_resname() not in ('HOH','DOD'):
                H += [a.coord for a in r]
    if not H: return 0.0
    H=np.array(H,float); n=0; tot=0
    for ch in st:
        for r in ch:
            if is_aa(r,standard=True) and r.id[1] in site:
                tot+=1
                A=np.array([a.coord for a in r],float)
                if np.sqrt(((A[:,None,:]-H[None,:,:])**2).sum(2)).min()<=cut: n+=1
    return (n/tot) if tot else None
prox=json.load(open('apo_proximity.json'))
out={}
print(f"  {'target':<9}{'site':>5}{'clean':>7}{'occ':>6}   median occupancy of the occupied apo")
for sym,v in prox.items():
    hid,lg=v['site'].split(':')
    site=contacts(hid,lg)
    if not site: continue
    fr=[]
    for pid in v['clean']+v['dirty']:
        o=occupancy(pid,site)
        if o is not None: fr.append((pid,o))
    clean=[p for p,o in fr if o<=0.50]; occ=[p for p,o in fr if o>0.50]
    med=np.median([o for _,o in fr if o>0.50]) if occ else 0.0
    out[sym]=dict(cohort=v['cohort'],site=v['site'],n_site=len(site),
                  clean=clean,occupied=occ,s1=v['s1'],recall=v['recall'])
    print(f"  {sym:<9}{len(site):>5}{len(clean):>7}{len(occ):>6}   {med:.0%}"
          f"   {'S@1 HIT' if v['s1'] else 'miss'}"
          f"{'   <-- <3 clean' if len(clean)<3 else ''}")
json.dump(out,open('apo_occupancy.json','w'),indent=1)
bad=[k for k,v in out.items() if len(v['clean'])<3]
print(f"\nfewer than 3 unoccupied apo members: {len(bad)} of {len(out)}")
print("  "+", ".join(bad) if bad else "  none")
ok=[v for k,v in out.items() if len(v['clean'])>=3]
print(f"\nSURVIVING targets: {len(ok)}   S@1 {sum(1 for v in ok if v['s1'])}/{len(ok)}")
