"""Benchmark runner.

Requires a local cache of mmCIF structures; set STRUCTURE_CACHE to point at it.
The stored outputs of this script are in ../results/ and can be checked without
re-running it, using ../verify/verify_numbers.py.

GATE 10 RE-RUN of the definitive colorectal benchmark, with the fixed detector.

Every number in this programme was produced by a LIGSITE whose voxel lattice is
anchored to the bounding-box minimum -- one extreme atom per axis -- so an
8-atom difference elsewhere in a protein re-phases the grid and flips the call.
See ../detector-fix/RESULT.md.

This is the only cohort with real answer keys, so it is the only place the fix
can be measured rather than argued. Three detectors, everything else identical:

    atlas       the original, lo = coords.min(0) - 5.0
    snapped     lo = floor((min - margin)/spacing)*spacing
    consensus   phase-averaged over 4 sub-voxel offsets

The per-target random-among-clusters baseline is recomputed for each, because
the fixed detector finds more pockets and a higher random rate is exactly how
that cost must be charged.

Original header follows.

DEFINITIVE colorectal run. Six gates, every one traceable to a defect found.

 1 TITLE GUARD      an entry whose title says complex/bound to/inhibitor/ligand
                    is not apo, whatever the metadata reports.
 2 COFACTOR STATE   read from ATOM records, never from
                    nonpolymer_bound_components.
 3 TARGET CHAIN     the answer key comes ONLY from chains the RCSB
                    polymer-entity records map to the target's UniProt
                    accession. Contacts on a partner protein are not the
                    target's site -- YAP1's key was TEAD's myristate pocket.
 4 SITE SIZE        >= 5 residues after chain restriction. CTNNB1 falls to 1.
 5 IDENTITY         an apo member must carry >= 80% of the site residues with
                    MATCHING amino-acid identity, on a single chain.
 6 OCCUPANCY        > 50% of site residues contacting a non-water heteroatom
                    means the pocket was held open; that member is not apo.

>= 3 members must survive 5 and 6 for a target to be scored at all.
"""
import json, os, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
_HERE=os.path.dirname(os.path.abspath(__file__))
_ROOT=os.path.dirname(_HERE)
sys.path.insert(0,os.path.join(_ROOT,'lib'))
sys.path.insert(0,os.path.join(_ROOT,'detector'))
import allatom as AAM
from enspock import cluster
from stable_pockets import pockets as SP, pockets_consensus as SPC
from probe_pockets import pockets_probe as PP
RMIN=float(os.environ.get('RMIN','2.0'))
MODE=os.environ.get('DETECTOR','atlas')
def _detect(coords,radii):
    if MODE=='atlas':     return AAM.pockets(coords,radii)
    if MODE=='snapped':   return SP(coords,radii,snap=True)
    if MODE=='consensus': return SPC(coords,radii)
    if MODE=='probe':     return PP(coords,radii,r_min=RMIN)
    raise SystemExit(f'unknown DETECTOR {MODE}')
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
P=MMCIFParser(QUIET=True)
S=os.environ.get("STRUCTURE_CACHE", os.path.join(_ROOT,"structures"))
def path(pid):
    for d in ('crcv2','crcv3','kras','cofcheck','gtp_check'):
        f=f"{S}/{d}/{pid}.cif"
        if os.path.exists(f) and os.path.getsize(f)>500: return f
    return None
def site_of(hid,lg,chains,cut=4.5):
    """Answer key: ligand contacts restricted to the TARGET chains only."""
    f=path(hid)
    if not f: return {}
    st=P.get_structure(hid,f)[0]
    la=[a for c in st for r in c if r.get_resname()==lg for a in r]
    if not la: return {}
    L=np.array([a.coord for a in la],float); out={}
    for ch in st:
        if ch.id not in chains: continue
        for r in ch:
            if not is_aa(r,standard=True): continue
            A=np.array([a.coord for a in r],float)
            if np.sqrt(((A[:,None,:]-L[None,:,:])**2).sum(2)).min()<=cut: out[r.id[1]]=r.get_resname()
    return out
def member_ok(pid, site):
    """Gates 5 and 6: identity-matched single chain, site not occupied."""
    f=path(pid)
    if not f: return None
    try: st=P.get_structure(pid,f)[0]
    except Exception: return None
    H=[a.coord for c in st for r in c if not is_aa(r,standard=True)
       and r.get_resname() not in ('HOH','DOD') for a in r]
    H=np.array(H,float) if H else None
    best=None
    for ch in st:
        d={r.id[1]:r for r in ch if is_aa(r,standard=True)}
        common=[n for n in site if n in d]
        if len(common)<0.80*len(site): continue
        agree=sum(1 for n in common if d[n].get_resname()==site[n])/len(common)
        if agree<0.80: continue
        if H is not None:
            occ=0
            for n in common:
                A=np.array([a.coord for a in d[n]],float)
                if np.sqrt(((A[:,None,:]-H[None,:,:])**2).sum(2)).min()<=4.5: occ+=1
            if occ/len(common)>0.50: continue
        if best is None or agree>best[1]: best=(ch.id,agree)
    return best
def pockets(pid, refset):
    AAM.DIR=os.path.dirname(path(pid))
    ch=AAM.parse_allatom(pid)
    if not ch: return None
    b=(-1,None)
    for cid,res in ch.items():
        if len(res)<50: continue
        ov=len(set(res)&refset)
        if ov>b[0]: b=(ov,cid)
    if b[1] is None: return None
    res=ch[b[1]]; coords=[];radii=[];ca={}
    for rn,(aa3,ats) in res.items():
        if 'CA' in ats: ca[rn]=ats['CA'][0]
        for an,(xyz,el) in ats.items():
            coords.append(xyz); radii.append(AAM.VDW.get(el,1.70))
    if len(coords)<100 or not ca: return None
    coords=np.array(coords,float); radii=np.array(radii,float)
    keys=np.array(sorted(ca)); CA=np.array([ca[j] for j in keys])
    out=[]
    for p in _detect(coords,radii):
        dd=np.linalg.norm(CA[:,None,:]-p['vox'][None,:,:],axis=2).min(axis=1)
        mem=sorted(int(keys[i]) for i in range(len(keys)) if dd[i]<8.0)
        if len(mem)>=4: out.append((np.array(mem),float(p['vol'])))
    return out
def zz(x):
    x=np.asarray(x,float); s=x.std(); return (x-x.mean())/s if s>0 else x*0

cm=json.load(open(f"{S}/chain_map.json"))
v2={t['symbol']:t for t in json.load(open(f"{S}/crc_v2_cohort.json"))}
v3={t['symbol']:t for t in json.load(open(f"{S}/crc_v3_cohort.json"))}
coh={**v2,**v3}
KIN={'EGFR','FGFR1','FGFR2','FGFR4','KIT','PDGFRA','RAF1','AKT1','BRAF','ERBB2','FGFR3',
     'ABL1','PIK3CA','MAPK1','BTK','LCK','SRC','BUB1','ERBB4'}
out=[]
for sym,v in cm.items():
    if not v['target_chains']:
        print(f"  {sym:<9} EXCLUDED: no chain maps to {v['acc']} in {v['holo']}",flush=True); continue
    site=site_of(v['holo'],v['lig'],set(v['target_chains']))
    if len(site)<5:
        print(f"  {sym:<9} EXCLUDED: site is {len(site)} residues after chain restriction",flush=True); continue
    members=[]
    for apo_id,dep,res,title in coh[sym]['apo'][:12]:
        if member_ok(apo_id,site): members.append(apo_id)
    if len(members)<3:
        print(f"  {sym:<9} EXCLUDED: {len(members)} members pass identity+occupancy",flush=True); continue
    ref=set()
    for p in members:
        AAM.DIR=os.path.dirname(path(p)); c=AAM.parse_allatom(p)
        if c: ref|=set().union(*[set(x) for x in c.values()])
    allp,used=[],0
    for mi,p in enumerate(members):
        try: pk=pockets(p,ref)
        except Exception: pk=None
        if not pk: continue
        used+=1
        for m,vv in pk: allp.append((mi,m,vv))
    if used<3 or not allp:
        print(f"  {sym:<9} EXCLUDED: {used} members yielded pockets",flush=True); continue
    cls=cluster(allp,0)
    sset=set(site)
    for c in cls:
        c['freq']=len(c['members'])/used; c['maxvol']=float(np.max(c['vols']))
        c['rec']=len(c['res']&sset)/len(sset)
    fv=np.array([c['freq'] for c in cls])*zz(np.log1p([c['maxvol'] for c in cls]))
    o=np.argsort(-fv)
    rank=next((i+1 for i,k in enumerate(o) if cls[k]['rec']>=0.25),None)
    rnd=float(np.mean([c['rec']>=0.25 for c in cls]))
    r=dict(symbol=sym,acc=v['acc'],holo=v['holo'],lig=v['lig'],site=len(site),
           members=used,clusters=len(cls),rank=rank,s1=bool(rank==1),
           s3=bool(rank is not None and rank<=3),rnd=rnd,
           recall=float(max(c['rec'] for c in cls)),kinase=sym in KIN,
           score=coh[sym]['score'])
    out.append(r)
    print(f"  {sym:<9} site {len(site):>2} members {used:>2} clusters {len(cls):>2}"
          f"  rank {str(rank):>4}  random {rnd:.0%}  recall {r['recall']:.2f}"
          f"  {'S@1' if r['s1'] else '   '}  {'kinase' if r['kinase'] else ''}",flush=True)
tag=MODE if MODE!='probe' else f'probe{RMIN:g}'
json.dump(out,open(f"gate10_{tag}.json","w"),indent=1)
s1=np.array([r['s1'] for r in out],float); rn=np.array([r['rnd'] for r in out],float)
kin=np.array([r['kinase'] for r in out])
rng=np.random.default_rng(20260910)
m=s1-rn; b=np.array([rng.choice(m,len(m),True).mean() for _ in range(20000)])
lo,hi=np.percentile(b,[2.5,97.5])
print(f"\nSCORED {len(out)}   S@1 {int(s1.sum())}/{len(out)} = {s1.mean():.0%}")
print(f"  paired margin {m.mean():+.1%}  95% CI [{lo:+.1%}, {hi:+.1%}]")
for tag,mask in (("kinase",kin),("non-kinase",~kin)):
    if mask.sum()==0: continue
    a=s1[mask]; q=rn[mask]; mm=a-q
    bb=np.array([rng.choice(mm,len(mm),True).mean() for _ in range(20000)])
    l2,h2=np.percentile(bb,[2.5,97.5])
    print(f"  {tag:<12} n={mask.sum():>2}  S@1 {int(a.sum())}/{mask.sum()}"
          f"  margin {mm.mean():+.0%} [{l2:+.0%}, {h2:+.0%}]")
