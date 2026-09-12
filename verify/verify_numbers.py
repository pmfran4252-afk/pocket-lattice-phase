#!/usr/bin/env python3
"""Re-derive every headline number in the manuscript from stored results.

Deliberately does NOT re-run any pipeline and needs no structures, no network and
no detector. It reads results/*.json and recomputes the statistics with
independent code and a different bootstrap seed from the one used in the paper,
so an arithmetic error in a pipeline's own summary would surface as a
disagreement here.

    python verify/verify_numbers.py
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")

def load(name):
    with open(os.path.join(RESULTS, name)) as f:
        return json.load(f)

atlas = {r['symbol']: r for r in load("crc_definitive.json")}
fixed = {r['symbol']: r for r in load("gate10_consensus.json")}
snap  = {r['symbol']: r for r in load("gate10_snapped.json")}
g7    = {r['symbol']: r for r in load("gate7.json")}
g8    = {r['symbol']: r for r in load("gate8.json")}
common = sorted(set(atlas) & set(fixed) & set(snap))

def rate(d, keep, kinase=None):
    rows = [d[s] for s in common if keep(s)]
    if kinase is not None:
        rows = [r for r in rows if r['kinase'] == kinase]
    return int(sum(r['s1'] for r in rows)), len(rows)

def margin(d, keep, kinase=None, seed=2026, draws=40000):
    rows = [d[s] for s in common if keep(s)]
    if kinase is not None:
        rows = [r for r in rows if r['kinase'] == kinase]
    m = np.array([r['s1'] - r['rnd'] for r in rows], float)
    rng = np.random.default_rng(seed)
    b = np.array([rng.choice(m, len(m), True).mean() for _ in range(draws)])
    lo, hi = np.percentile(b, [2.5, 97.5])
    return m.mean(), lo, hi

everything = lambda s: True
g7ok = lambda s: g7.get(s, {}).get('verdict', 'PASS') != 'EXCLUDE'
g8ok = lambda s: g8.get(s, {}).get('pass40') is not False
gated = lambda s: g7ok(s) and g8ok(s)

print(__doc__.strip().splitlines()[0])
print("=" * 72)
print(f"\n{len(common)} targets common to all three runs\n")
print(f"  {'detector':<12}{'overall':>14}{'non-kinase':>14}{'kinase':>12}")
for tag, d in (("original", atlas), ("snapped", snap), ("phase-avg", fixed)):
    a, b = rate(d, everything)
    c, e = rate(d, everything, kinase=False)
    f, g = rate(d, everything, kinase=True)
    print(f"  {tag:<12}{f'{a}/{b} = {a/b:.0%}':>14}{f'{c}/{e} = {c/e:.0%}':>14}"
          f"{f'{f}/{g}':>12}")

print(f"\nafter answer-key quality filtering (gates 7 and 8)\n")
print(f"  {'detector':<12}{'non-kinase':>14}{'margin over random (fresh seed)':>36}")
for tag, d in (("original", atlas), ("phase-avg", fixed)):
    c, e = rate(d, gated, kinase=False)
    m, lo, hi = margin(d, gated, kinase=False)
    print(f"  {tag:<12}{f'{c}/{e} = {c/e:.0%}':>14}"
          f"{f'{m:+.1%}  95% CI [{lo:+.1%}, {hi:+.1%}]':>36}")

resc = [s for s in common if not atlas[s]['s1'] and fixed[s]['s1']]
lost = [s for s in common if atlas[s]['s1'] and not fixed[s]['s1']]
print(f"\nrescued by the fix : {sorted(resc)}")
print(f"lost to the fix    : {sorted(lost) if lost else 'none'}")
print(f"\n  {'target':<9}{'original':>10}{'snapped':>10}{'phase-avg':>12}   (rank)")
for s in sorted(resc):
    print(f"  {s:<9}{str(atlas[s]['rank']):>10}{str(snap[s]['rank']):>10}"
          f"{str(fixed[s]['rank']):>12}")

EXPECT = {"original non-kinase": (rate(atlas, everything, kinase=False), (7, 13)),
          "snapped non-kinase":  (rate(snap,  everything, kinase=False), (8, 13)),
          "phase-avg non-kinase":(rate(fixed, everything, kinase=False), (10, 13)),
          "original gated":      (rate(atlas, gated, kinase=False), (5, 10)),
          "phase-avg gated":     (rate(fixed, gated, kinase=False), (8, 10)),
          "nothing lost":        (lost, [])}
print(f"\n{'=' * 72}\nagainst the values stated in the manuscript:\n")
fails = 0
for k, (got, want) in EXPECT.items():
    good = got == want
    fails += not good
    print(f"  {'OK  ' if good else 'FAIL'}  {k:<24} {got}   expected {want}")
print(f"\n{'all figures reproduce' if not fails else str(fails) + ' DISAGREEMENTS'}")
sys.exit(1 if fails else 0)
