"""LIGSITE pocket detection with a lattice-locked grid origin.

THE BUG
-------
Atlas's devcache LIGSITE (benchmarks/cryptic_dispersion/devcache/allatom.py)
anchors the voxel lattice to the atomic bounding box:

    lo = coords.min(0) - 5.0

`coords.min(0)` is set by a SINGLE extreme atom. Change the atom list anywhere
in the protein and lo moves by a fraction of a voxel, which re-phases the whole
lattice against every atom at once: occupancy flips at the margins, PSP counts
change, connected components merge or split, and the min_vox floor is crossed
in either direction.

Observed on beta-catenin. 8Y0G and 8Z61 are the same crystal form at the same
resolution, 508 identical residues, 0.00 A RMSD across the site of interest,
differing by 8 atoms elsewhere in the chain. One of those atoms sits 0.142 A
lower in y, so:

    8Y0G  lo = [ -7.339  -53.839  -15.666]
    8Z61  lo = [ -7.339  -53.981  -15.666]

a 0.142 A shift -- 14% of a voxel -- applied to the entire protein. The site is
detected at Jaccard 0.75 in one and 0.12 in the other.

Atlas's production grid (qid_modalities/protein/structure/grid.py) canonicalises
the frame to principal axes, which fixes the ROTATION dependence it documents.
It does not fix this: build_voxel_grid still takes its origin from the
bounding-box minimum of the canonical coordinates, so a single extreme atom
still sets the lattice phase.

THE FIX
-------
Snap the origin down to a multiple of the spacing:

    lo = np.floor((coords.min(0) - margin) / spacing) * spacing

Now the lattice phase is a property of the absolute coordinate frame, not of
which atom happens to be extreme. A small change in the bounding box either
leaves the origin alone or moves it by exactly one full voxel -- and a whole-
voxel shift is a relabelling of indices, not a re-phasing: every atom keeps its
position relative to the nearest voxel centre, so occupancy is unchanged.

Everything else -- probe radius, 7 directions, psp_min, min_vox, connected
components -- is identical to the original, so any change in output is
attributable to the grid phase and to nothing else.
"""
from __future__ import annotations
import numpy as np

PROBE = 1.4
DIRS7 = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]


def _shift(a, d):
    out = np.zeros_like(a)
    src, dst = [], []
    for ax in range(3):
        k = d[ax]
        if k == 0:
            src.append(slice(None)); dst.append(slice(None))
        elif k > 0:
            src.append(slice(0, -k)); dst.append(slice(k, None))
        else:
            src.append(slice(-k, None)); dst.append(slice(0, k))
    out[tuple(dst)] = a[tuple(src)]
    return out


def _blocked(occ, d):
    acc = np.zeros_like(occ)
    cur = occ.copy()
    for _ in range(max(occ.shape)):
        cur = _shift(cur, d)
        if not cur.any():
            break
        acc |= cur
    return acc


def _cc(mask):
    """6-connected components, iterative label propagation."""
    lab = np.zeros(mask.shape, dtype=np.int32)
    idx = np.argwhere(mask)
    if len(idx) == 0:
        return lab, np.array([], dtype=int)
    lab[mask] = np.arange(1, mask.sum() + 1)
    while True:
        prev = lab.copy()
        for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            nb = _shift(lab, d)
            lab = np.where(mask & (nb > lab), nb, lab)
        if np.array_equal(prev, lab):
            break
    vals = lab[mask]
    uniq, inv = np.unique(vals, return_inverse=True)
    out = np.zeros(mask.shape, dtype=np.int32)
    out[mask] = inv + 1
    return out, np.bincount(inv)


def pockets(coords, radii, spacing=1.0, psp_min=5, min_vox=30, margin=5.0,
            snap=True):
    """LIGSITE pockets. `snap=False` reproduces the original, unstable behaviour."""
    coords = np.asarray(coords, float)
    radii = np.asarray(radii, float)
    raw_lo = coords.min(0) - margin
    lo = np.floor(raw_lo / spacing) * spacing if snap else raw_lo
    hi = coords.max(0) + margin
    dims = np.maximum(np.ceil((hi - lo) / spacing).astype(int), 3)
    if dims.prod() > 9_000_000:
        spacing = float(np.cbrt(dims.prod() / 9_000_000)) * spacing
        lo = np.floor(raw_lo / spacing) * spacing if snap else raw_lo
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

    lab, sizes = _cc((~occ) & (psp >= psp_min))
    out = []
    for rank, li in enumerate(np.argsort(sizes)[::-1]):
        if sizes[li] < min_vox:
            break
        vx = np.argwhere(lab == li + 1)
        out.append(dict(rank=rank, nvox=int(sizes[li]),
                        vol=float(sizes[li] * spacing ** 3),
                        centre=vx.mean(0) * spacing + lo,
                        vox=vx * spacing + lo))
    return out


def pockets_consensus(coords, radii, spacing=1.0, psp_min=5, min_vox=30,
                      margin=5.0, offsets=4, keep=0.5):
    """Phase-averaged LIGSITE.

    Snapping the origin removes the atom-list dependence, but a second
    instability remains and is inherent to the method: `psp_min` and `min_vox`
    are hard thresholds on a 1 A lattice, so a pocket sitting near either one
    flips on coordinate noise of a tenth of an Angstrom. The grid phase is an
    arbitrary choice in the first place -- there is no reason the lattice should
    sit where it happens to sit.

    So run the detection at several sub-voxel phases and keep the voxels that
    are called a pocket in at least `keep` of them. The arbitrary choice is
    averaged out rather than snapped to one value, and a voxel that only
    qualifies at one particular phase no longer decides the outcome.

    Returns the same dicts as `pockets`, with `nphase` added: the mean fraction
    of phases in which that pocket's voxels were called.
    """
    coords = np.asarray(coords, float)
    radii = np.asarray(radii, float)
    base = np.floor((coords.min(0) - margin) / spacing) * spacing
    shifts = [(np.array([i, j, k], float) / offsets) * spacing
              for i in range(offsets) for j in range(offsets) for k in range(offsets)]
    # sample the phase cube on its diagonal rather than exhaustively: offsets^3
    # grids is needless, and the diagonal spans the full range in every axis
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
        mask = (~occ) & (psp >= psp_min)
        for v in np.argwhere(mask):
            # map the voxel centre back to an absolute cell, so votes from
            # different phases are counted in one common frame
            w = v * spacing + lo + 0.5 * spacing
            key = tuple(np.floor(w / spacing).astype(int))
            votes[key] = votes.get(key, 0) + 1

    thr = max(1, int(np.ceil(keep * len(shifts))))
    cells = np.array([k for k, n in votes.items() if n >= thr], dtype=int)
    if len(cells) == 0:
        return []
    frac = np.array([votes[tuple(c)] / len(shifts) for c in cells])
    origin = cells.min(0)
    idx = cells - origin
    dims = idx.max(0) + 1
    mask = np.zeros(dims, dtype=bool)
    mask[idx[:, 0], idx[:, 1], idx[:, 2]] = True
    lab, sizes = _cc(mask)
    out = []
    for rank, li in enumerate(np.argsort(sizes)[::-1]):
        if sizes[li] < min_vox:
            break
        vx = np.argwhere(lab == li + 1)
        world = (vx + origin) * spacing + 0.5 * spacing
        sel = np.array([votes.get(tuple(c), 0) for c in (vx + origin)]) / len(shifts)
        out.append(dict(rank=rank, nvox=int(sizes[li]),
                        vol=float(sizes[li] * spacing ** 3),
                        centre=world.mean(0), vox=world,
                        nphase=float(sel.mean())))
    return out
