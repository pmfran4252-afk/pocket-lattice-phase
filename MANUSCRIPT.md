# A lattice-phase instability in grid-based pocket detection, and what fixing it reveals about benchmark hygiene

**Draft v0.1 — 2026-09-12. Not submitted. All numbers verified 2026-09-12 by
independent re-derivation from stored per-target results (23/23 checks).**

---

## Abstract

Grid-based cavity detection underlies much of structural drug discovery. We report
that a widely-derived LIGSITE implementation anchors its voxel lattice to the atomic
bounding box, so the lattice phase is set by a single extreme atom. Adding eight atoms
anywhere in a protein shifts the origin by a fraction of a voxel and re-phases the
grid against every atom simultaneously. We demonstrate this with two PDB entries of
the same protein, same crystal form, same resolution, whose coordinates across the
site of interest are identical to two decimal places. Most of the cavity map is
reproduced between them, but one cavity measures 27 voxels in one entry and 30 in the
other against a 30-voxel floor, so the site is reported in one and discarded in the
other. The instability is not diffuse: it is a hard integer threshold sitting where
the phase noise is.

Snapping the origin to the lattice removes the atom-list dependence; averaging over
sub-voxel phases additionally removes a threshold brittleness that snapping leaves
behind. On a 20-target benchmark with answer keys taken from published co-crystals,
the two fixes together raise success-at-rank-1 on non-kinase targets from 7/13 to
10/13, and from 5/10 to 8/10 after answer-key quality filtering, with no target
regressing.

Building that benchmark required ten filters, each added after a specific defect
produced a wrong number. Several concern metadata rather than geometry: an RCSB
summary field silently under-reports bound ligands, a second under-reports engineered
substitutions, and answer keys assembled without chain restriction can describe a
partner protein's pocket rather than the target's. We give the failure each filter
prevents and the cost of omitting it.

We report four predictions this framework caused us to withdraw, and one it did not,
which is registered prospectively with its residue list, code and inputs frozen by
hash.

---

## 1. Introduction

A drug acts by occupying a cavity. Finding the cavity is the first computational step
in most structure-based campaigns, and grid-based geometric detectors — LIGSITE
(Hendlich et al. 1997) and its many descendants — remain in wide use because they are
fast, parameter-light and require no training data.

These methods rasterise a protein onto a voxel lattice, classify empty voxels by
burial, and cluster the buried ones. The lattice is a measurement instrument. This
paper is about a property of that instrument that we believe is under-appreciated:
**where the lattice is placed is usually decided by the coordinates themselves, and
that decision is not stable.**

We found this while attempting to nominate a site prospectively, and the discovery
inverted the project: results we had believed for nine days had to be re-derived, four
were withdrawn, and the benchmark figure we had been quoting was substantially
pessimistic rather than optimistic.

---

## 2. Results

### 2.1 A failing pair

PDB entries 8Y0G and 8Z61 are both human β-catenin, space group C222₁, 2.50 Å, from
the same depositor, with 508 modelled residues each. Across the 23 residues that
concerned us they are identical:

    whole-chain Cα RMSD            0.104 Å   (no superposition applied)
    Cα RMSD across the 23 residues 0.000 Å

They differ by eight atoms elsewhere in the chain. One of those atoms lies 0.142 Å
lower in *y* than anything in 8Y0G, and the detector's grid origin is
`coords.min(0) − margin`:

    8Y0G   origin = [ −7.339,  −53.839,  −15.666 ]
    8Z61   origin = [ −7.339,  −53.981,  −15.666 ]

a shift of 0.142 Å, 14% of the 1.0 Å voxel, applied to the whole protein at once.

**The consequence is narrow and specific, and worth stating precisely.** Most of the
cavity map is stable. The detector returns three cavities for 8Y0G and four for
8Z61, and three of them correspond one-to-one at Jaccard 0.92, 0.96 and 1.00. What
differs is a single cavity:

    entry   rank  voxels  residues  Jaccard vs the site of interest
    8Y0G       0      75        24              0.09
    8Y0G       1      67        23              0.00
    8Y0G       2      66        26              0.02
    8Z61       0      76        25              0.09
    8Z61       1      72        23              0.00
    8Z61       2      69        26              0.02
    8Z61       3      30        18              0.78   <-- present only here

Removing the minimum-size floor shows the cavity is present in both structures and
was never absent from either:

    8Y0G    27 voxels   Jaccard 0.87
    8Z61    30 voxels   Jaccard 0.78

**The floor is `min_vox = 30`.** A three-voxel difference, produced entirely by the
grid phase, decides whether the site is reported at all. One structure reports it;
the other, at identical coordinates, does not.

### 2.2 The mechanism is a hard threshold sitting in the phase noise

Two things combine. `coords.min(0)` is an extremum, so it is fixed by exactly one atom
per axis, and any edit to the atom list that changes that atom re-phases the lattice
against every atom — the affected site need not be near the atom responsible. And the
detector's decisions are hard integer cutoffs on that lattice: burial along ≥5 of 7
directions, connected components of ≥30 voxels. A cavity whose size lands near the
floor is decided by the phase rather than by the protein.

This predicts that the instability should be concentrated in marginal cavities and
absent from large ones, which is what §2.1 shows, and that it should be worst for the
small, shallow sites that cryptic-pocket detection is specifically aimed at.

We quantified the resulting instability by perturbing coordinates at σ = 0.10 Å, well
below the coordinate error at 2.2–2.8 Å resolution, and asking how often the site was
still reported:

    structure   role in the ensemble    recovered at σ = 0.10 Å
    8Z0U        supports the cluster           78%
    2Z6H        does NOT support it            75%
    9I8K        supports                       62%
    8Z61        supports                       60%
    8Y0G        does NOT support it            42%

**There is no separation between supporters and non-supporters.** A reproducibility
statistic of the form "the site appears in N of M structures" is therefore not N
successes and M−N failures, and any null model conditioning on those per-structure
calls inherits the problem.

### 2.3 Two fixes, addressing different things

**Origin snapping.** Setting `origin = floor((min − margin)/spacing) · spacing` makes
the phase a property of the absolute coordinate frame. A small change to the bounding
box now either leaves the origin alone or moves it by exactly one whole voxel, which
is a relabelling of indices rather than a re-phasing: every atom retains its position
relative to the nearest voxel centre. Applied to the failing pair, both entries adopt
origin [−8, −54, −16], both report the cavity, and both return Jaccard 0.83 against
the site (0.91 with phase averaging). All reported Jaccards in §2.1–2.3 are against
the registered 23-residue list of §2.7.

**Phase averaging.** Snapping fixes the atom-list dependence but not a second problem
that is intrinsic to the method: burial and size thresholds are hard cutoffs on a 1 Å
lattice, so a cavity near either flips on a tenth of an Ångström. Since the phase is
arbitrary to begin with, we run the detection at four sub-voxel offsets and retain
voxels called in a majority.

    recovery at σ = 0.10 Å      original   snapped   phase-averaged
    8Y0G                          30%        55%         100%
    8Z61                          50%        35%          95%
    2Z6H                          85%        95%         100%

Both fixes are necessary. On the benchmark, BCL2 moves from rank 9 to rank 13 under
snapping alone and only reaches rank 1 under phase averaging.

**Verification of faithfulness.** Our reimplementation with snapping disabled
reproduces the original implementation voxel-for-voxel on test structures (8Y0G:
[75, 67, 66]; 8Z61: [76, 72, 69, 30]), so the changes reported here are attributable
to the fixes and not to a rewrite.

### 2.4 Benchmark impact

Twenty targets with answer keys derived from published co-crystals, never shown to the
detector. Success is the top-ranked cluster overlapping the key. Per-target random
baselines are recomputed for each detector, so a detector that finds more cavities is
charged for the additional chances.

    detector          overall        non-kinase      non-kinase margin over random
    original        14/20   70%      7/13   54%      +31.9%  95% CI [+3.8%, +59.3%]
    snapped         15/20   75%      8/13   62%      +42.4%  95% CI [+16.3%, +67.1%]
    phase-averaged  17/20   85%     10/13   77%      +57.5%  95% CI [+32.0%, +80.0%]

Kinases are reported separately and remain 7/7 throughout: locating a kinase ATP site,
the largest cavity on the protein, does not test cryptic-pocket detection. APC moves
from rank 8 to 1, BCL2 from 9 to 1, EP300 from 2 to 1, and AKT1 re-enters the scored
cohort at rank 1 having previously been excluded for yielding too few cavities.
**No target regressed.**

After answer-key quality filtering (§2.5, gates 7–8), non-kinase success is 8/10 =
80%, margin +60.1%, 95% CI [+32.5%, +83.2%], against 5/10 = 50% and +29.4% before the
fix.

### 2.5 Ten filters, each from a defect that produced a wrong number

Half of what determines a benchmark result here is not the predictor. It is whether
the answer key describes what one thinks it describes.

| # | Filter | The defect it prevents |
|---|---|---|
| 1 | Title guard | An entry titled "…and ligand" whose metadata field lists none. `"complexed with"` does not contain `"complex with"`; both spellings are needed. |
| 2 | Cofactor state from ATOM records | `rcsb_entry_info.nonpolymer_bound_components` under-reports. It omitted nucleotides across every small GTPase tested, and admitted TYMS 1HVY ("*complexed with dUMP and Raltitrexed*") and 2ONB as apo. |
| 3 | Target-chain restriction | Contacts collected from all chains describe a partner's pocket. One target's key was another protein's myristate site; three co-crystals did not contain the target protein at all. Only 15 of 29 targets are chain A alone. |
| 4 | Minimum site size after restriction | One target falls to a single residue once restricted. |
| 5 | Sequence identity of apo members | Numbering offsets of +10 between entries of the same protein are common. |
| 6 | Site occupancy | A "holo structure with the ligand deleted at parse time" is not apo. A first version using minimum distance was wrong and rejected a physiological cofactor resting beside the site; fraction-of-site occupancy is correct. |
| 7 | Ligand density quality | Worst-instance RSCC ≥ 0.80. Excluded one rank-1 hit whose key fits its map at 0.452 with 50% occupancy. |
| 8 | Cross-crystal ligand reproducibility | Does the same ligand contact the same residues in independent crystals? Median Jaccard ≥ 0.40. |
| 9 | Cluster support and method mixing | A cluster must appear in ≥2 structures. NMR ensembles pack ~6% less densely than crystals and return cavities 4–7× larger; mixing them with X-ray members without checking lets an NMR artefact outrank a real site. |
| 10 | Stability | Perturb at σ = 0.10 Å and confirm the call is stable before quoting any reproducibility statistic (§2.2). |

Two further observations belong with these.

**An RCSB summary field under-reports engineered substitutions.**
`rcsb_polymer_entity.pdbx_mutation` reports no mutation for a structure carrying five
substitutions and an eight-residue deletion, because the depositor annotated them as
`conflict` rather than `engineered mutation`. Across 135 entries the field disagreed
with a direct sequence diff twelve times, eleven in that direction. Construct
independence must be measured by alignment against the reference sequence.

**A specificity statistic that rejects its own positive control is not evidence.** We
briefly withdrew a prediction on the basis of a perturb-and-compare-to-random-patches
test, then found the test rejects a confirmed, clinically drugged pocket at the same
value. With a 46-residue patch on a 166-residue protein the statistic is dominated by
patch size relative to protein size. We report it only as a relative comparator
against a known positive.

### 2.6 Four withdrawn predictions

Each was withdrawn by a filter added after an earlier failure.

**APOBEC3B.** A convergence between two paralogues at Jaccard 0.51 proved to rest on a
single NMR entry, reproducing in none of seven crystals (median Jaccard 0.038). Filter
9. Re-tested after the detector fix — because the same pattern proved to be false
negatives elsewhere — and the withdrawal held (median 0.050). The fixed detector
instead returns the catalytic active site, containing all four zinc ligands, in three
independent crystals: a positive control, and a demonstration that the original
detector had been ranking an NMR packing artefact above the real site.

**MYD88.** A cluster reported in 4/4 structures at p = 0.0078 falls to 2/4 under the
fixed detector, where random patches reach 2/4 86% of the time. The site itself is
present in all four structures; what fails is the merge, with one structure scoring
Jaccard 0.43 and 0.42 against two others under a 0.50 clustering threshold. **The
residue-level observation stands and the statistic does not** — an instructive case of
an arbitrary threshold downstream of the detector deciding a published number.

**TYMS.** Cleared the reproducibility null but failed specificity against the positive
control by sixfold, and was not recovered at all in one structure under perturbation.

**A β-catenin support statistic.** Reported as 5/8 at p = 0.0008, withdrawn under
filter 10, and re-derived as 8/8 after the detector fix — the three apparent negatives
were artefacts of the lattice phase.

### 2.7 One prediction, registered

The site that survived is on β-catenin, 23 residues on the convex outer face of the
armadillo superhelix. It is registered with its residue list, the detector's SHA-256
and all eight input-structure hashes frozen before any other method was applied to that
surface.

    reproducibility        8 of 8 wild-type apo structures, p = 0.0000
    crystal contacts       5% of the site buried by symmetry mates vs 12% background
    partner interfaces     22 partners with solved complexes; Jaccard 0.03 against
                           their union; zero site residues within 5 Å of any partner
                           across 13 superposed complexes
    groove geometry        projection on the groove axis: groove +7.2 Å, protein mean
                           −1.2 Å, site −28.3 Å, consistent across all 8 structures
    published ligands      Jaccard 0.00 against the small-molecule contact set
    positive control       HRAS, same pipeline and run, co-crystals withheld: 12/12
                           structures, Jaccard 0.63 against its own ligand contacts

Four falsification conditions are specified in advance. We make **no functional
claim**: nothing is known to occur on that surface, and a ligand bound there may have
no biological effect.

---

## 3. Discussion

The instability we report is a property of how the lattice is positioned, not of
LIGSITE's burial criterion. Any method that derives a grid origin from an extremum of
the input coordinates is exposed. We note that a current production implementation
canonicalises the molecular frame to principal axes — correctly solving the *rotation*
dependence, which its documentation describes — while still deriving the origin from
the bounding-box minimum of the canonicalised coordinates, and recording that choice
as a known property rather than a defect.

The more general point concerns benchmark construction. Six of our ten filters address
the answer key rather than the predictor, and omitting them inflates or deflates
results in both directions: filter 3 alone corrupted three rows and admitted three
targets whose co-crystals do not contain the target protein. A cryptic-pocket
benchmark assembled from PDB metadata without these checks is measuring something
other than what it reports.

### Limitations

**The cohort is small.** Twenty targets, of which thirteen are non-kinase and ten
survive answer-key filtering. The headline 8/10 rests on ten targets and the interval
is correspondingly wide. This is a demonstration that the fixes matter, not a
precise estimate of accuracy.

**Filter 8 is a spot check here.** Most co-crystals are deposited once per compound,
so cross-crystal ligand reproducibility was testable on only a minority of rows. It is
powerful prospectively and weak retrospectively.

**We have not established that phase averaging is optimal.** Four offsets and a
majority vote were chosen a priori and not tuned. A probe-occupancy criterion that
would additionally reject shallow surface grooves is under evaluation in a separate
module and is not used for any result reported here.

**The registered prediction is unvalidated.** No ligand has been modelled at the site
and no experiment has been performed.

---

## 4. Methods

### 4.1 Cavity detection

Protein heavy atoms are rasterised onto a cubic lattice of 1.0 Å spacing with a 5.0 Å
margin. Each atom occupies voxels within its van der Waals radius plus a 1.4 Å probe,
so the occupied region is the solvent-excluded rather than the van der Waals volume.
Empty voxels are scored by protein–solvent–protein burial along seven directions —
the three Cartesian axes and four cube diagonals — and a voxel is a cavity voxel when
it is buried along at least five of the seven. Cavity voxels are grouped into
6-connected components and components of at least 30 voxels are returned, ranked by
volume.

**Origin placement.** The original implementation sets the origin to
`coords.min(0) − margin`. The snapped variant sets it to
`floor((coords.min(0) − margin)/spacing) · spacing`. The phase-averaged variant
evaluates the lattice at four offsets, `−(i/4)·spacing` for `i = 0…3` applied equally
on all three axes, maps each resulting cavity voxel back to a common absolute cell,
and retains cells called in at least half of the offsets before connected-component
analysis. Offsets and the majority threshold were fixed a priori and not tuned.

**Equivalence control.** The reimplementation used here reproduces the original
voxel-for-voxel when snapping is disabled, verified on test structures by comparing
component sizes. All differences reported are therefore attributable to origin
placement and phase averaging.

### 4.2 Mapping cavities to residues, and clustering across an ensemble

A cavity is expressed as the set of residues whose Cα lies within 8.0 Å of any of its
voxels; cavities touching fewer than four residues are discarded. Across the members
of an apo ensemble, cavities are merged into clusters by greedy agglomeration on
residue-set Jaccard overlap with a merge threshold of 0.50. A cluster's *support* is
the number of distinct structures contributing to it.

Clusters are ranked by `frequency × z(log(1 + maximum volume))`, where frequency is
support divided by the number of contributing structures and the z-score is taken
across the clusters of that target.

### 4.3 Benchmark cohort and success criterion

Targets are proteins with at least three apo structures surviving the filters in §4.4
and at least one co-crystal supplying an answer key. The answer key is the set of
residues within 4.5 Å of the bound ligand, restricted to chains that RCSB
polymer-entity records map to the target's UniProt accession.

A cluster *hits* when it recovers at least 25% of the key. Success-at-rank-1 is
scored when the top-ranked cluster hits. The per-target random baseline is the
fraction of that target's clusters that would hit, so a detector returning more
clusters is charged for the additional opportunities; the reported margin is
success-at-rank-1 minus that baseline, averaged across targets. Confidence intervals
are 95% percentile bootstrap over targets, 20,000 resamples. Kinases are reported
separately throughout.

### 4.4 The ten filters

Filters 1–6 admit or reject structures; 7–8 admit or reject answer keys; 9–10 govern
what may be quoted from a result.

1. **Title guard.** Entries whose title contains any of `complex with`, `complexed`,
   `bound to`, `inhibitor`, `covalent`, `compound`, `ligand`, `in complex` are not apo.
2. **Bound components.** Read from `non_polymer_entity_ids`, never from
   `rcsb_entry_info.nonpolymer_bound_components`. A component counts as drug-like at
   ≥12 heavy atoms excluding an explicit cofactor and buffer list.
3. **Target chains.** Answer-key contacts are collected only from chains whose
   polymer entity carries the target's UniProt accession.
4. **Site size.** At least five residues after chain restriction.
5. **Identity.** An apo member must carry ≥80% of the key's residue positions on a
   single chain with ≥80% amino-acid agreement.
6. **Occupancy.** A member is rejected when more than 50% of key residues contact a
   non-water heteroatom within 4.5 Å.
7. **Ligand density.** Worst-instance real-space correlation coefficient ≥0.80, from
   the wwPDB validation pipeline.
8. **Cross-crystal reproducibility.** Median Jaccard ≥0.40 between the ligand contact
   sets of independent crystals of the same ligand, where more than one exists.
9. **Support and method.** A cluster must appear in at least two structures. Mean
   heavy-atom neighbour count within 6.0 Å is compared across ensemble members before
   NMR and X-ray structures are combined.
10. **Stability.** Coordinates are perturbed with isotropic Gaussian noise at
    σ = 0.10 Å and the detection repeated; a call whose recovery rate does not separate
    from that of size-matched random patches is not quoted.

**Construct independence.** Where apo members must be independent, substitutions are
counted by aligning each entity's canonical sequence to the UniProt reference with
`difflib`, retaining matching blocks of ≥10 residues and summing
`min(reference advance, entity advance)` across each gap. Leading and trailing
mismatch is truncation or an expression tag and is not counted; internal insertions
are not counted. `rcsb_polymer_entity.pdbx_mutation` is not used.

### 4.5 Null models

**Reproducibility null.** Each detected cavity is replaced by a spatial patch of the
same residue count drawn from the same structure — the *n* residues nearest a
uniformly chosen seed residue by Cα distance — and the ensemble is re-clustered
identically. The reported *p* is the fraction of draws in which some cluster reaches
at least the observed support. Draw counts are stated per experiment.

**Specificity comparison.** Under the same σ = 0.10 Å perturbation, the recovery rate
of the candidate site is compared with that of size-matched random patches from the
same structure. This is reported only relative to a known-positive control, for the
reason given in §2.5.

### 4.6 Geometry

**Crystal contacts.** Symmetry mates are generated from the `REMARK 290` matrices and
unit cell, over all lattice translations in {−1, 0, 1}³, retaining copies whose
bounding box approaches the reference. A residue is buried when any of its atoms lies
within 4.5 Å of a mate, and site burial is compared with the whole-chain rate.

**Superposition.** Complexes are superposed onto the apo reference by Kabsch
alignment on common Cα positions, with the rotation forced right-handed; the partner
chains are carried by the same transform and closest approach is measured per residue.

**Groove projection.** For each structure, the vector from the molecular centroid to
the centroid of the known-interface residues defines the groove direction. Each
residue's Cα is projected onto that unit vector. No superposition and no partner
coordinates are used.

### 4.7 Registration

The prediction, the detector source and every input structure are hashed with SHA-256.
The registration artifact carries the residue list, endpoints, falsification conditions
and stated limitations, and is itself hashed over its canonical JSON serialisation with
the hash field removed. Detector source and artifact are committed before any further
method is applied to the nominated surface.

## 5. Data and code availability

All per-target results, the ten filters, both detector implementations, the
registration artifact and its hashes are in the project repository. Every figure in
§2 is re-derivable from stored JSON without re-running any pipeline; the verification
script that does so is included.

## 6. Author contributions and competing interests

*(To be completed.)*
