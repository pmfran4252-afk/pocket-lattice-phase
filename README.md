# A lattice-phase instability in grid-based pocket detection

Code and data for the manuscript in [`MANUSCRIPT.md`](MANUSCRIPT.md).

Grid-based cavity detectors rasterise a protein onto a voxel lattice. In a widely
derived LIGSITE implementation the lattice origin is `coords.min(0) − margin` — an
extremum, fixed by a single atom. Change that atom and the lattice re-phases against
every atom in the protein.

## Reproduce the central result in one command

```bash
pip install -r requirements.txt
python reproduce_failing_pair.py
```

Downloads two PDB entries and shows the failure directly. Expected output:

```
   whole-chain Ca RMSD        0.104 A
   Ca RMSD across the site    0.000 A   <-- identical

   8Y0G  origin = [ -7.339 -53.839 -15.666]
   8Z61  origin = [ -7.339 -53.981 -15.666]
   shift = 0.142 A  ->  14% of the 1.0 A voxel

   entry     voxels  J vs site
   8Y0G          27       0.87
   8Z61          30       0.78
```

Two entries of the same protein, same crystal form, same resolution, whose
coordinates across the site of interest are **identical to two decimal places**. The
cavity is present in both. It measures 27 voxels in one and 30 in the other against a
floor of 30, so one reports the site and the other discards it. **A three-voxel
difference, produced entirely by the grid phase, decides the result.**

Most of the cavity map is stable — three cavities correspond one-to-one at Jaccard
0.92, 0.96 and 1.00. The instability is not diffuse; it is a hard integer threshold
sitting where the phase noise is.

## Verify the benchmark without running anything

```bash
python verify/verify_numbers.py
```

Recomputes every headline figure from `results/*.json` with independent code and a
different bootstrap seed. No structures, no network, seconds to run. This is the
check a sceptical reader should run first.

## Layout

```
detector/     stable_pockets.py  original, snapped and phase-averaged variants
              probe_pockets.py   an unvalidated alternative criterion, used for
                                 nothing reported in the paper
lib/          structures.py      standalone mmCIF reader
              clustering.py      cross-member cavity clustering
benchmark/    the cohort runner and the ten answer-key filters
results/      per-target results as written by each run
registration/ a prospective site prediction, frozen by SHA-256 before any other
              method was applied to that surface
verify/       re-derives the paper's numbers from results/
```

`lib/` replaces two files that previously lived in an unpublished repository. Both
replacements were verified to be drop-in exact: the reader is identical to the
original on 40 structures (chains, residues, atom names, coordinates, elements), and
the clustering is identical on 3,000 randomised inputs.

## What this repository does not claim

The registered prediction in `registration/` is **unvalidated**. No ligand has been
modelled at the site and no experiment has been performed. It carries four
falsification conditions specified in advance, and makes no claim that a ligand bound
there would have any biological effect.

The benchmark is twenty targets, of which ten survive answer-key quality filtering.
The headline figure rests on n = 10 with a correspondingly wide interval. It
demonstrates that the fixes matter; it is not a precise estimate of accuracy.

## Licence

MIT, see [`LICENSE`](LICENSE). `stable_pockets.py` is frozen by hash in the
registration artifact and should not be edited in place.
