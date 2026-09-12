"""Standalone mmCIF reading — replaces an unpublished dependency.

The benchmark originally called a parser that lives in a separate, unlicensed
repository. Nothing outside that machine could run it. This module reproduces
that parser's contract exactly, on top of Biopython's MMCIF2Dict, so the package
is self-contained and a reviewer can execute it.

The contract, preserved bit for bit:

    model 1 only                     pdbx_PDB_model_num == '1'
    polymer atoms only               group_PDB == 'ATOM'
    an accepted residue list         the 20 standard plus MSE PTR SEP TPO CSO
                                     CME CAS LLP -- tested BEFORE canonicalisation
    one altloc                       label_alt_id in {'.', '?', 'A'}
    no hydrogens                     type_symbol not in {H, D}
    author numbering throughout      auth_asym_id / auth_seq_id / auth_atom_id
    modified residues canonicalised  MSE -> MET and so on, AFTER the accept test

Returns {chain: {resnum: (resname, {atom_name: ((x, y, z), element)})}}.

Verified against the original on the benchmark structures: identical chains,
residues, atom names, coordinates and elements.
"""
from __future__ import annotations
import os
from collections import defaultdict
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

VDW = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'SE': 1.90, 'P': 1.80}
VDW_DEFAULT = 1.70
PROBE = 1.40

ACCEPTED = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU',
            'LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL',
            'MSE','PTR','SEP','TPO','CSO','CME','CAS','LLP'}
CANON = {'MSE': 'MET', 'PTR': 'TYR', 'SEP': 'SER', 'TPO': 'THR', 'CSO': 'CYS',
         'CME': 'CYS', 'CAS': 'CYS', 'LLP': 'LYS'}


def parse_structure(path):
    """Parse one mmCIF file. Returns None when the file is absent or too small."""
    if not os.path.exists(path) or os.path.getsize(path) < 500:
        return None
    d = MMCIF2Dict(path)
    n = len(d.get('_atom_site.group_PDB', []))
    if n == 0:
        return None

    def col(name, default=None):
        v = d.get('_atom_site.' + name)
        if v is None:
            return [default] * n
        return v

    grp   = col('group_PDB')
    model = col('pdbx_PDB_model_num', '1')
    comp  = col('auth_comp_id')
    alt   = col('label_alt_id', '.')
    elem  = col('type_symbol', '')
    xs, ys, zs = col('Cartn_x'), col('Cartn_y'), col('Cartn_z')
    seq   = col('auth_seq_id')
    asym  = col('auth_asym_id')
    atom  = col('auth_atom_id')

    chains = defaultdict(dict)
    for i in range(n):
        if (model[i] or '1') != '1':
            continue
        if grp[i] != 'ATOM':
            continue
        if comp[i] not in ACCEPTED:
            continue
        if (alt[i] or '.') not in ('.', '?', 'A'):
            continue
        el = (elem[i] or '').upper()
        if el in ('H', 'D'):
            continue
        try:
            xyz = (float(xs[i]), float(ys[i]), float(zs[i]))
            rn = int(seq[i])
        except (TypeError, ValueError):
            continue
        aa3 = CANON.get(comp[i], comp[i])
        rec = chains[asym[i]].setdefault(rn, (aa3, {}))
        rec[1][atom[i]] = (xyz, el)
    return dict(chains)


def atoms_and_radii(residues):
    """Flatten one chain's residue dict into (coords, radii, ca_by_resnum)."""
    coords, radii, ca = [], [], {}
    for rn, (aa3, ats) in residues.items():
        if 'CA' in ats:
            ca[rn] = ats['CA'][0]
        for _, (xyz, el) in ats.items():
            coords.append(xyz)
            radii.append(VDW.get(el, VDW_DEFAULT))
    return coords, radii, ca
