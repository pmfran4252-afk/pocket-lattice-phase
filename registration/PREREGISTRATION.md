# Pre-registration: a ligandable site on the convex outer face of β-catenin

    protocol_id   ctnnb1-outer-face-v0.1
    frozen        2026-09-11
    target        CTNNB1 / catenin beta-1, UniProt P35222
    status        prediction registered; no ligand has been modelled, no
                  functional evidence claimed

Written before any docking, fragment work, or literature search on this surface.
Nothing downstream of this file may alter the residue list.

## The prediction

**A small molecule can be made to bind a pocket lined by these 23 residues,**
in the numbering of PDB 8Y0G:

    543 544 545 546 547 548 · 566 567 568 570 571
    598 600 601 602 603 604 605 606 607 608 · 610 611

    A543 H544 Q545 D546 T547 Q548 · M566 E567 E568 V570 E571
    L598 V600 Q601 L602 L603 Y604 S605 P606 I607 I608 · I610 Q611

One connected component at 8 Å, 18 Å across, radius of gyration 7.7 Å, 4% of the
armadillo-repeat domain.

**The site is on the convex outer face of the armadillo superhelix, not in the
concave groove where every known β-catenin partner and every clinical agent
binds.** That is the substance of the claim, and it is what makes it worth
registering rather than a restatement of the Wnt field.

## What was done, and by what code

    8 wild-type apo X-ray structures, 1.00-2.80 A, no substitution against
    P35222 over the span each covers:

        8Y0G  2Z6H  8Y14  8Z0U  8Z10  8Z5J  8Z61  9I8K

    detector   src/loki/stable_pockets.py :: pockets_consensus
               sha256 e693d692846e9fb0d3787aa987d0c3fedfe3ef32a53802ab4b326abce07aa6e9
               LIGSITE with the voxel origin snapped to the lattice and averaged
               over 4 sub-voxel phases. The unfixed detector, which anchored the
               grid to a single extreme atom, is documented in
               ../detector-fix/RESULT.md.

    support    the cluster appears in 8 of 8 structures
    null       p = 0.0000, 2000 draws, each real pocket replaced by a spatial
               patch of the same size from the same structure and re-clustered

## The evidence this rests on, each by a separate path

    reproducibility      8 of 8 apo structures, p = 0.0000
    not crystal packing  5% of the site buried by symmetry mates vs 12%
                         background (REMARK 290, 4.5 A)
    not an oligomer      every apo deposit is single-chain
    not a partner site   22 distinct partners with solved complexes -- TCF3,
                         TCF7L2, APC, ICAT, Axin-1, Axin, BCL9, E-cadherin,
                         beta-TrCP x2, GSK-3beta, MDM2, PAK4, TAX1BP3, cyclic
                         peptides, FOG-001's Helicon polypeptide. Union of all
                         22 interfaces is 179 residues; Jaccard 0.03; 19 of the
                         24 residues in the earlier list appear at NO interface
    outer face, test 1   13 complexes superposed at 0.96-2.74 A RMSD; ZERO site
                         residues within 5 A of any partner in any complex;
                         median closest approach 14.5 A against 5.8 A for the
                         groove
    outer face, test 2   projection on the groove axis, no superposition and no
                         partner structures: groove +7.2 A, protein average
                         -1.2 A, this site -28.3 A, consistent across all 8
    not the ligand site  Jaccard 0.00 against the small-molecule contact set on
                         full-length beta-catenin (8 residues at 302-380, from
                         9I8W). Keys built only from entries whose CTNNB1 chain
                         is the full protein -- three of six candidate
                         co-crystals contain it only as a 33-residue
                         phosphopeptide bound to beta-TrCP
    positive control     HRAS, same pipeline, same run: 12 of 12 structures,
                         Jaccard 0.63 against its own ligand contacts,
                         p = 0.0000, co-crystals never shown to the caller

## Declared in advance

  * **One prediction.** The 23 residues above are final. No re-ranking, no
    re-clustering, no threshold change, whatever any subsequent result shows.
  * The detector is frozen by hash. A later version does not get to re-derive
    this list.
  * Under the previous, broken detector the same site came out as 24 residues
    at 5 of 8 support. Jaccard between the two lists is 0.880 — 608 added, 569
    and 599 dropped. Both lists are recorded so neither can be quietly
    preferred later.
  * The specificity statistic that briefly withdrew this site is not used as a
    gate, because it rejects HRAS's confirmed pocket at the same p = 0.100.
    It is reported as a comparison against that control and nothing more.

## How this is falsified

**CONFIRMED if** a structure or a localised biophysical measurement shows a
small molecule or fragment bound with a contact set reaching **Jaccard ≥ 0.30**
against the 23 residues. Any source: crystallography, cryo-EM, NMR chemical-shift
perturbation, HDX-MS, or a covalent adduct mapped by MS.

**FALSIFIED if any of the following:**

    F1  a fragment or compound screen that demonstrably samples this surface
        -- >= 500 fragments, readout able to localise on the ARM domain --
        returns hits elsewhere on beta-catenin and none reaching Jaccard 0.30
        here
    F2  >= 3 further wild-type apo beta-catenin X-ray structures are deposited
        and the cluster fails to reach 50% support in the enlarged ensemble,
        scored by the hashed detector
    F3  the surface is shown to be the binding interface of a protein partner
        not among the 22 tested, at Jaccard >= 0.30
    F4  the site is shown to be an artefact of the construct shared by the six
        C222(1) depositions -- for instance if a structure from a different
        expression construct and space group lacks it while the hashed detector
        is applied unchanged

**NOT falsification:** absence of a published ligand. No one has looked here.

## What this cannot establish, stated before the fact

**Function.** Novelty is not mechanism. Every β-catenin partner binds the
concave groove because that is where the biology is, and nothing is known to
happen on the convex face. A ligand at this site might do nothing. That risk is
different in kind from the ones this analysis retired — those were about whether
the site is real.

**Competitive position.** β-catenin is not an uncrowded target: 51 clinical
studies, with FOG-001, tegavivint, E7386, PRI-724/OP-724 and ALN-BCAT live, and
nothing approved. Every direct agent binds the groove. This site is
differentiated from all of them, which is the argument for it — and it will
matter less if any of them succeeds.

**Construct breadth.** Six of the eight apo structures are one depositor in
C222(1). Two pairs share identical unit cells. The ensemble spans three space
groups and two depositors, not eight independent experiments. F4 exists for
this reason.

## Provenance

    ../detector-fix/RESULT.md          the lattice bug and the fix
    ../crc-definitive/GATE10.md        the fix measured on answer keys,
                                       non-kinase 54% -> 77%
    ../oncogene-screen/RESULT.md       how this target was selected
    ../oncogene-screen/GATE10.md       this re-run, and why TYMS was dropped
    ../detector-stability/RESULT.md    Gate 10, and what it withdrew
    ../null-calibration/RESULT.md      why p in 0.06-0.12 is failure
