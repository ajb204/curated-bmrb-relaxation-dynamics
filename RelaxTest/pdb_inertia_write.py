#!/usr/bin/env python3
"""
pdb_inertia_write.py

Compute principal axes from a PDB, rotate coordinates into that internal frame,
and write a new PDB with rotated coordinates.

Usage (CLI):
    python pdb_inertia_write.py input.pdb output_rotated.pdb --mass-weighted
    python pdb_inertia_write.py input.pdb output_cov_rotated.pdb --covariance
"""

import numpy as np
from typing import Tuple, List, Optional
import argparse
import sys

# Minimal atomic mass table
_ATOMIC_MASSES = {
    "H": 1.00784, "C": 12.00, "N": 14.00307, "O": 15.999, "S": 32.06,
    "P": 30.97376, "FE": 55.845, "MG": 24.305, "ZN": 65.38, "CL": 35.45,
}

def _element_from_pdb_line(line: str) -> str:
    el = line[76:78].strip()
    if el:
        return el.upper()
    atom_name = line[12:16].strip()
    letters = "".join([c for c in atom_name if c.isalpha()])
    return letters[:2].upper() if letters else ""

def load_pdb_coords(pdb_path: str) -> Tuple[np.ndarray, List[str], List[int], List[str]]:
    """
    Parse PDB file returning:
      coords: (N,3) numpy array
      elements: list of element tokens (len N)
      atom_line_indices: list of the line indices in file that correspond to atoms (for replacement)
      lines: full file lines (list of strings)
    """
    coords = []
    elements = []
    atom_line_indices = []
    with open(pdb_path, "r") as fh:
        lines = fh.readlines()

    for i, line in enumerate(lines):
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                continue
            coords.append([x, y, z])
            elements.append(_element_from_pdb_line(line))
            atom_line_indices.append(i)

    return np.array(coords, dtype=float), elements, atom_line_indices, lines

def coords_from_array(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    if arr.ndim == 1 and arr.size == 3:
        return arr.reshape(1, 3)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError("coords must be shape (N,3)")
    return arr

def _masses_from_elements(elements: List[str]) -> np.ndarray:
    masses = []
    for el in elements:
        if not el:
            masses.append(1.0); continue
        el_up = el.strip().upper()
        if el_up in _ATOMIC_MASSES:
            masses.append(_ATOMIC_MASSES[el_up])
        else:
            if len(el_up) > 1 and el_up[0] in _ATOMIC_MASSES:
                masses.append(_ATOMIC_MASSES[el_up[0]])
            else:
                masses.append(1.0)
    return np.array(masses, dtype=float)

def center_of_mass(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> np.ndarray:
    coords = coords_from_array(coords)
    if masses is None:
        return coords.mean(axis=0)
    masses = np.asarray(masses, dtype=float)
    total = masses.sum()
    if total == 0:
        raise ValueError("sum of masses is zero")
    return (masses[:, None] * coords).sum(axis=0) / total

def inertia_tensor(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> np.ndarray:
    coords = coords_from_array(coords)
    if masses is None:
        masses = np.ones(coords.shape[0], dtype=float)
    else:
        masses = np.asarray(masses, dtype=float)
    com = center_of_mass(coords, masses)
    rel = coords - com
    x = rel[:, 0]; y = rel[:, 1]; z = rel[:, 2]; m = masses
    I_xx = np.sum(m * (y**2 + z**2))
    I_yy = np.sum(m * (x**2 + z**2))
    I_zz = np.sum(m * (x**2 + y**2))
    I_xy = -np.sum(m * x * y)
    I_xz = -np.sum(m * x * z)
    I_yz = -np.sum(m * y * z)
    I = np.array([[I_xx, I_xy, I_xz],
                  [I_xy, I_yy, I_yz],
                  [I_xz, I_yz, I_zz]], dtype=float)
    return I

def compute_principal_axes(coords: np.ndarray,
                           elements: Optional[List[str]] = None,
                           mass_weighted: bool = False
                           ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    coords = coords_from_array(coords)
    if mass_weighted and elements is not None:
        masses = _masses_from_elements([e.upper() for e in elements])
    else:
        masses = None
    com = center_of_mass(coords, masses)
    I = inertia_tensor(coords, masses)
    eigvals, eigvecs = np.linalg.eigh(I)
    idx = np.argsort(eigvals)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    if np.linalg.det(eigvecs) < 0:
        eigvecs[:, 0] = -eigvecs[:, 0]
    return I, eigvals, eigvecs, com

def compute_covariance_axes(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute covariance matrix of positions (shape spread), get eigenvectors.
    Returns covariance matrix, eigenvalues (descending), eigenvectors (columns), com
    """
    coords = coords_from_array(coords)
    com = coords.mean(axis=0)
    rel = coords - com
    cov = np.cov(rel.T, bias=True)  # population covariance
    eigvals, eigvecs = np.linalg.eigh(cov)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    if np.linalg.det(eigvecs) < 0:
        eigvecs[:, 0] = -eigvecs[:, 0]
    return cov, eigvals, eigvecs, com

def rotate_to_internal_frame(coords: np.ndarray, com: np.ndarray, eigvecs: np.ndarray) -> np.ndarray:
    coords = coords_from_array(coords)
    rel = coords - com
    coords_int = (eigvecs.T @ rel.T).T
    return coords_int

def _format_pdb_atom_line(line: str, x: float, y: float, z: float) -> str:
    """
    Replace only columns 31–54 (0-based 30:54) with new XYZ coordinates.
    Preserves everything else exactly as in the original PDB.
    """
    # Ensure line is long enough
    orig = line.rstrip("\n")
    if len(orig) < 54:
        orig = orig.ljust(54)

    new_xyz = f"{x:8.3f}{y:8.3f}{z:8.3f}"
    new_line = orig[:30] + new_xyz + orig[54:]
    # ensure trailing newline
    return new_line + ("\n" if not new_line.endswith("\n") else "")

def write_pdb_with_coords(lines: List[str], atom_line_indices: List[int],
                          new_coords: np.ndarray, out_path: str) -> None:
    """
    Write a PDB file to out_path replacing the coordinates in the atom lines
    with new_coords (N,3). Assumes atom_line_indices and new_coords correspond.
    """
    if len(atom_line_indices) != new_coords.shape[0]:
        raise ValueError("number of atom lines and new_coords rows must match")
    new_lines = list(lines)  # copy
    for idx, (line_idx, coord) in enumerate(zip(atom_line_indices, new_coords)):
        orig_line = lines[line_idx]
        x, y, z = float(coord[0]), float(coord[1]), float(coord[2])
        new_line = _format_pdb_atom_line(orig_line, x, y, z)
        new_lines[line_idx] = new_line
    with open(out_path, "w") as fh:
        fh.writelines(new_lines)

# ---------- CLI ----------
def main_cli():
    parser = argparse.ArgumentParser(
        description="Rotate PDB coordinates into principal/internal axis frame and write rotated PDB."
    )
    parser.add_argument("pdb_in", help="input PDB file")
    parser.add_argument("pdb_out", help="output rotated PDB file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--mass-weighted", action="store_true",
                       help="compute inertia tensor mass-weighted by element masses")
    group.add_argument("--covariance", action="store_true",
                       help="use covariance of positions (shape) instead of inertia tensor")
    parser.add_argument("--preview", action="store_true", help="print summary info")
    args = parser.parse_args()

    coords, elements, atom_idxs, lines = load_pdb_coords(args.pdb_in)
    if coords.size == 0:
        print("No ATOM/HETATM coordinates found in PDB.", file=sys.stderr)
        sys.exit(1)

    if args.covariance:
        cov, eigvals, eigvecs, com = compute_covariance_axes(coords)
        method = "covariance"
    else:
        I, eigvals, eigvecs, com = compute_principal_axes(coords, elements, mass_weighted=args.mass_weighted)
        method = "inertia (mass-weighted)" if args.mass_weighted else "inertia (unit masses)"

    coords_internal = rotate_to_internal_frame(coords, com, eigvecs)
    write_pdb_with_coords(lines, atom_idxs, coords_internal, args.pdb_out)

    if args.preview:
        print(f"Method: {method}")
        print("Center (used):", com)
        print("Eigenvalues:", eigvals)
        print("Eigenvectors (columns):\n", eigvecs)
        print("Wrote rotated PDB to:", args.pdb_out)

if __name__ == "__main__":
    main_cli()
