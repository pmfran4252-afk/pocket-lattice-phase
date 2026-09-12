"""Cross-member cavity clustering — replaces an unpublished dependency.

Greedy single-pass agglomeration on residue-set Jaccard overlap. Each cavity is
compared against the clusters formed so far; it joins the best-overlapping one if
that overlap reaches `jac`, and otherwise starts a new cluster.

The pass is ORDER-DEPENDENT by construction, and the original behaviour is
preserved exactly, including that dependence: cavities must be presented in the
same order (member index ascending, then the detector's own volume ranking) for
results to match. This is a property of the published method, not an
implementation detail, and it is why `allp` is built in a fixed order upstream.
"""
from __future__ import annotations


def cluster(allp, jac=0.5):
    """allp: iterable of (member_index, residue_iterable, volume).

    Returns a list of dicts: res (set), members (set of member indices),
    vols (list), counts (int).
    """
    clusters = []
    for mi, mem, vol in allp:
        s = set(int(x) for x in mem)
        best, bj = None, 0.0
        for c in clusters:
            j = len(s & c['res']) / max(len(s | c['res']), 1)
            if j > bj:
                best, bj = c, j
        if best is not None and bj >= jac:
            best['res'] |= s
            best['members'].add(mi)
            best['vols'].append(vol)
            best['counts'] += 1
        else:
            clusters.append(dict(res=set(s), members={mi}, vols=[vol], counts=1))
    return clusters
