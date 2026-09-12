#!/usr/bin/env python3
"""Reproduces the central result of the paper in one command.

Downloads two PDB entries, shows that their coordinates across the site of
interest are identical, and shows that the original detector nevertheless calls
the site in one and not the other -- then shows both fixes removing the
discrepancy.

    python reproduce_failing_pair.py

Needs numpy and biopython, and network access on first run.
"""
import os, sys, urllib.request
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "detector"))
from structures import parse_structure, atoms_and_radii
from stable_pockets import pockets, pockets_consensus

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "structures")
PAIR = ("8Y0G", "8Z61")
SITE = {543,544,545,546,547,548,566,567,568,570,571,598,600,601,602,603,604,605,
        606,607,608,610,611}

def fetch(pid):
    os.makedirs(CACHE, exist_ok=True)
    f = os.path.join(CACHE, f"{pid}.cif")
    if not (os.path.exists(f) and os.path.getsize(f) > 500):
        print(f"  downloading {pid} ...", flush=True)
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid}.cif", f)
    return f

def load(pid):
    ch = parse_structure(fetch(pid))
    cid = max(ch, key=lambda c: len(ch[c]))
    co, ra, ca = atoms_and_radii(ch[cid])
    return np.array(co, float), np.array(ra, float), ca

def best_jaccard(pk, ca):
    keys = np.array(sorted(ca)); X = np.array([ca[k] for k in keys]); best = 0.0
    for q in pk:
        d = np.linalg.norm(X[:, None, :] - q['vox'][None, :, :], axis=2).min(axis=1)
        m = set(int(keys[i]) for i in range(len(keys)) if d[i] < 8.0)
        if len(m) < 4:
            continue
        best = max(best, len(m & SITE) / len(m | SITE))
    return best

print(__doc__.strip().splitlines()[0]); print("=" * 72)
A, B = load(PAIR[0]), load(PAIR[1])
print(f"\n1. THE TWO ENTRIES ARE THE SAME PROTEIN, AND THE SITE IS IDENTICAL\n")
print(f"   {PAIR[0]}  {len(A[0])} atoms        {PAIR[1]}  {len(B[0])} atoms"
      f"        difference {len(B[0]) - len(A[0])}")
common = sorted(set(A[2]) & set(B[2]))
w = np.sqrt((((np.array([A[2][k] for k in common]) -
               np.array([B[2][k] for k in common])) ** 2).sum(1)).mean())
sres = [r for r in sorted(SITE) if r in A[2] and r in B[2]]
s = np.sqrt((((np.array([A[2][r] for r in sres]) -
               np.array([B[2][r] for r in sres])) ** 2).sum(1)).mean())
print(f"   whole-chain Ca RMSD        {w:.3f} A   (no superposition applied)")
print(f"   Ca RMSD across the site    {s:.3f} A   <-- identical")

print(f"\n2. THE GRID ORIGIN DIFFERS ANYWAY, BECAUSE ONE ATOM MOVED\n")
for pid, (co, ra, ca) in zip(PAIR, (A, B)):
    print(f"   {pid}  origin = {np.round(co.min(0) - 5.0, 3)}")
d = (B[0].min(0) - 5.0) - (A[0].min(0) - 5.0)
print(f"   shift = {np.round(d, 3)} A  ->  {100 * abs(d).max():.0f}% of the 1.0 A voxel,"
      f" applied to the whole protein")

print(f"\n3. MOST OF THE CAVITY MAP IS STABLE. ONE CAVITY IS NOT.\n")
def cavities(co, ra, ca, **kw):
    keys = np.array(sorted(ca)); X = np.array([ca[k] for k in keys]); out = []
    for q in pockets(co, ra, snap=False, **kw):
        d = np.linalg.norm(X[:, None, :] - q['vox'][None, :, :], axis=2).min(axis=1)
        m = set(int(keys[i]) for i in range(len(keys)) if d[i] < 8.0)
        out.append((q['nvox'], m))
    return out
def J(a, b): return len(a & b) / len(a | b) if a | b else 0.0
print(f"   {'entry':8}{'voxels':>8}{'residues':>10}{'J vs site':>11}")
for pid, (co, ra, ca) in zip(PAIR, (A, B)):
    for nv, m in cavities(co, ra, ca):
        print(f"   {pid:8}{nv:>8}{len(m):>10}{J(m, SITE):>11.2f}")
    print()
print("   The minimum-size floor is 30 voxels. Removing it:\n")
print(f"   {'entry':8}{'voxels':>8}{'J vs site':>11}")
for pid, (co, ra, ca) in zip(PAIR, (A, B)):
    best = (0, 0.0)
    for nv, m in cavities(co, ra, ca, min_vox=1):
        j = J(m, SITE)
        if j > best[1]: best = (nv, j)
    print(f"   {pid:8}{best[0]:>8}{best[1]:>11.2f}")
print(f"\n   The cavity is present in BOTH. It measures 27 voxels in one and 30 in")
print(f"   the other, against a floor of 30. A three-voxel difference, produced")
print(f"   entirely by the grid phase, decides whether the site is reported.")

print(f"\n4. BOTH FIXES REMOVE THE DISCREPANCY\n")
def best_jaccard(pk, ca):
    keys = np.array(sorted(ca)); X = np.array([ca[k] for k in keys]); best = 0.0
    for q in pk:
        d = np.linalg.norm(X[:, None, :] - q['vox'][None, :, :], axis=2).min(axis=1)
        m = set(int(keys[i]) for i in range(len(keys)) if d[i] < 8.0)
        if len(m) < 4: continue
        best = max(best, len(m & SITE) / len(m | SITE))
    return best
print(f"   {'entry':8}{'original':>12}{'snapped':>11}{'phase-averaged':>17}")
for pid, (co, ra, ca) in zip(PAIR, (A, B)):
    print(f"   {pid:8}"
          f"{best_jaccard(pockets(co, ra, snap=False), ca):>12.2f}"
          f"{best_jaccard(pockets(co, ra, snap=True), ca):>11.2f}"
          f"{best_jaccard(pockets_consensus(co, ra), ca):>17.2f}")
print(f"\n   Jaccards are against the registered 23-residue site.\n")
