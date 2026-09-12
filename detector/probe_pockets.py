"""Probe-occupancy pocket detection. A SEPARATE module by design.

`stable_pockets.py` is frozen by hash in
`benchmarks/ctnnb1-prereg/preregistration-v0.1.json`. Appending to it would
change its sha256 and break that artifact, so nothing here touches it. This
module is additive and unvalidated until it is scored against the 20-target
benchmark alongside the detector already in use.

WHY
---
LIGSITE asks whether a voxel is buried along >= 5 of 7 ray directions. A shallow
surface groove satisfies that, which is why this programme keeps returning
extended grooves and then caveating them: MYD88 45 residues over 31% of its
domain, NT5C2 94 residues, APC and BCL2 both flat interfaces. A 1.4 A water
probe fits in a surface dimple; a drug-sized probe does not.

It also exposed a second problem downstream. Pockets are matched between
structures by the set of residues within 8 A of any voxel, so a pocket's identity
tracks voxel extent rather than where the cavity is. On MYD88 that put three
views of one site either side of a hard Jaccard 0.50 merge threshold and
collapsed 4/4 to 2/4.

A FIRST VERSION OF THIS FILE WAS WRONG
--------------------------------------
Testing `r_min <= r_fit <= r_max` alone returned a single blob covering 75-84%
of the domain, volume 29,000-68,000 A^3: the entire solvent-accessible shell. A
sphere of 3-6 A fits perfectly well just outside a convex protein surface, so
`r_max` does not distinguish a cavity from the outside. fpocket escapes this
because an alpha sphere is a Voronoi vertex and must touch four atoms; a lattice
has no such constraint.

So the burial test is kept and the probe test is added to it, rather than
replacing it. Burial says "inside a cavity"; the probe says "big enough for a
drug". Neither is sufficient alone.

THE CRITERION
-------------
For each empty voxel centre, the largest sphere that fits without clashing is

    r_fit = min over atoms of ( |v - c| - vdw(c) )

A voxel is pocket-accessible when ALL THREE hold:

    not occupied        outside the solvent-excluded volume
    psp >= psp_min      buried along >= psp_min of 7 ray directions, the
                        LIGSITE test -- this is what says "cavity" rather than
                        "outside"
    r_fit >= r_min      a drug-sized sphere fits here, not merely a water --
                        this is the new part, and what rejects shallow grooves

`r_min` is the whole point. LIGSITE's burial test is satisfied by a surface
groove 1.4 A deep; requiring a 3 A sphere to fit is what a ligand actually
needs. The idea is fpocket's (Le Guilloux, Schmidtke, Tuffery 2009), evaluated
on the shared lattice rather than on a Voronoi tessellation.

Phase averaging is kept from `stable_pockets.pockets_consensus`, because the
lattice phase is an arbitrary choice whichever criterion is applied to it, and
averaging it out is what made detection stable under sub-resolution noise.

Each returned pocket carries `r_fit_max` -- the largest probe that fits anywhere
inside it -- which is a physical handle the voxel count does not give, and a
candidate matching key that does not depend on voxel extent.
"""
from __future__ import annotations
import numpy as np

PROBE = 1.4          # water, for the solvent-excluded surface
R_MIN_DEFAULT = 3.0  # drug-sized probe
R_MAX_DEFAULT = 6.0  # above this the sphere is in bulk solvent

from .stable_pockets import _cc, _shift, _blocked, DIRS7  # reuse, do not reimplement


def _r_fit(coords, radii, lo, dims, spacing, r_max):
    """Largest clash-free sphere radius at each voxel centre.

    Only atoms within (r_max + radius + spacing) of a voxel can constrain its
    minimum, so each atom is applied over a local box. Voxels no atom reaches
    keep +inf and are bulk solvent, which the r_max test then rejects.
    """
    out = np.full(dims, np.inf, dtype=np.float32)
    for c, r in zip(coords, radii):
        reach = r_max + r + spacing
        i0 = np.maximum(np.floor((c - reach - lo) / spacing).astype(int), 0)
        i1 = np.minimum(np.ceil((c + reach - lo) / spacing).astype(int) + 1, dims)
        if (i1 <= i0).any():
            continue
        gx = (np.arange(i0[0], i1[0]) * spacing + lo[0] + 0.5 * spacing - c[0]) ** 2
        gy = (np.arange(i0[1], i1[1]) * spacing + lo[1] + 0.5 * spacing - c[1]) ** 2
        gz = (np.arange(i0[2], i1[2]) * spacing + lo[2] + 0.5 * spacing - c[2]) ** 2
        d = np.sqrt(gx[:, None, None] + gy[None, :, None] + gz[None, None, :]) - r
        sub = out[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]]
        np.minimum(sub, d.astype(np.float32), out=sub)
    return out


def pockets_probe(coords, radii, spacing=1.0, r_min=R_MIN_DEFAULT,
                  r_max=R_MAX_DEFAULT, psp_min=5, min_vox=20, margin=5.0,
                  offsets=4, keep=0.5):
    """Probe-occupancy pockets, phase-averaged.

    `min_vox` is lower than LIGSITE's 30 because the r_min test has already
    removed shallow volume: 20 voxels of drug-accessible space is a real
    cavity, where 20 voxels of water-accessible space is not.
    """
    coords = np.asarray(coords, float)
    radii = np.asarray(radii, float)
    base = np.floor((coords.min(0) - margin) / spacing) * spacing
    shifts = [(np.array([i, i, i], float) / offsets) * spacing for i in range(offsets)]

    votes = {}
    for s in shifts:
        lo = base - s
        hi = coords.max(0) + margin
        dims = np.maximum(np.ceil((hi - lo) / spacing).astype(int), 3)
        occ = np.zeros(dims, dtype=bool)
        for c, r in zip(coords, radii + PROBE):
            i0 = np.maximum(np.floor((c - r - lo) / spacing).astype(int), 0)
            i1 = np.minimum(np.ceil((c + r - lo) / spacing).astype(int) + 1, dims)
            if (i1 <= i0).any():
                continue
            gx = (np.arange(i0[0], i1[0]) * spacing + lo[0] - c[0]) ** 2
            gy = (np.arange(i0[1], i1[1]) * spacing + lo[1] - c[1]) ** 2
            gz = (np.arange(i0[2], i1[2]) * spacing + lo[2] - c[2]) ** 2
            d2 = gx[:, None, None] + gy[None, :, None] + gz[None, None, :]
            occ[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] |= d2 <= r * r
        psp = np.zeros(dims, dtype=np.int8)
        for d in DIRS7:
            psp += (_blocked(occ, d) & _blocked(occ, tuple(-x for x in d)) & ~occ).astype(np.int8)
        rf = _r_fit(coords, radii, lo, dims, spacing, r_max)
        mask = (~occ) & (psp >= psp_min) & (rf >= r_min) & (rf <= r_max)
        if not mask.any():
            continue
        for v in np.argwhere(mask):
            w = v * spacing + lo + 0.5 * spacing
            key = tuple(np.floor(w / spacing).astype(int))
            prev = votes.get(key)
            r = float(rf[v[0], v[1], v[2]])
            if prev is None:
                votes[key] = [1, r]
            else:
                prev[0] += 1
                prev[1] = max(prev[1], r)

    thr = max(1, int(np.ceil(keep * len(shifts))))
    cells = np.array([k for k, (n, r) in votes.items() if n >= thr], dtype=int)
    if len(cells) == 0:
        return []
    origin = cells.min(0)
    idx = cells - origin
    grid = np.zeros(idx.max(0) + 1, dtype=bool)
    grid[idx[:, 0], idx[:, 1], idx[:, 2]] = True
    lab, sizes = _cc(grid)
    out = []
    for rank, li in enumerate(np.argsort(sizes)[::-1]):
        if sizes[li] < min_vox:
            break
        vx = np.argwhere(lab == li + 1)
        world = (vx + origin) * spacing + 0.5 * spacing
        rs = [votes[tuple(c)][1] for c in (vx + origin) if tuple(c) in votes]
        out.append(dict(rank=rank, nvox=int(sizes[li]),
                        vol=float(sizes[li] * spacing ** 3),
                        centre=world.mean(0), vox=world,
                        r_fit_max=float(max(rs)) if rs else 0.0,
                        r_fit_mean=float(np.mean(rs)) if rs else 0.0))
    return out


# STATUS, 2026-09-12
# -----------------
# UNVALIDATED. Do not use for any reported result until the r_min sweep against
# the 20-target benchmark is complete and `pockets_consensus` is beaten on the
# same answer keys.
#
# What is known so far:
#
#   r_min = 3.0 with the burial test is too strict. MYD88 returns ZERO pockets
#   and CTNNB1 returns one 12-residue pocket overlapping the registered site at
#   2 of 23. A drug-sized sphere is 6 A across, and a cryptic pocket is by
#   definition still closed in the apo structure, so there is a real tension
#   between "admits a drug-sized probe in apo" and "is cryptic". The sweep is
#   what settles where that line sits, not argument.
#
#   Measuring "largest sphere that fits inside a site" is itself confounded the
#   same way the first version of this detector was: without a burial
#   constraint the answer is a sphere sitting in bulk solvent beside the site,
#   not enclosed in it. All four sites tested returned 4.3 A that way --
#   including HRAS's always-open ligand site, which should have been
#   distinguishable. Any future measurement of probe size must carry the burial
#   mask with it.
