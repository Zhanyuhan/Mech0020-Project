#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced PLY->STL converter (non-destructive)

Functionality Overview
- Read original STL (triangular mesh) and use only its topology faces.
- Read PLY point cloud (original / PC±2σ).
- Rigid ICP (SVD-Umeyama) to roughly align PLY to STL vertex distribution (align after normalization for better robustness).
- Use nearest neighbor to map each STL vertex to the nearest point in PLY point cloud, obtaining "deformed STL vertices".
- Optional Laplacian (umbrella) smoothing on "deformed STL vertices", trying to keep boundaries fixed.
- **Will not overwrite original files**, output to new directory.

Directory Convention (consistent with existing code/output)
- PLY directory: <ply_dir>/anatomical_<ID>/
    - anatomical_<ID>_original.ply
    - anatomical_<ID>_PC1_minus2sigma.ply
    - anatomical_<ID>_PC1_plus2sigma.ply
    - anatomical_<ID>_PC2_minus2sigma.ply
    - anatomical_<ID>_PC2_plus2sigma.ply
    - anatomical_<ID>_PC3_minus2sigma.ply
    - anatomical_<ID>_PC3_plus2sigma.ply
- STL directory: <stl_dir>/<ID>.stl

Example Usage
    python ply_to_stl_converter_smooth.py \\
        --ply_dir "C:/Users/Administrator/Desktop/MyScalismoProject/data/Final" \\
        --stl_dir "C:/Users/Administrator/Desktop/MyScalismoProject/data/ArotaModels" \\
        --output_dir "C:/Users/Administrator/Desktop/MyScalismoProject/data/GeneratedSTL_Smoothed" \\
        --icp_iters 10 --smooth_iters 10 --smooth_lambda 0.3 --all


Dependencies
- numpy
- scipy (recommended, provides cKDTree; if unavailable, fallback to O(N^2) nearest neighbor, slower with large data)
"""

import os
import sys
import argparse
import struct
import numpy as np
from pathlib import Path

try:
    from scipy.spatial import cKDTree
except Exception:
    cKDTree = None


# ------------------------------- I/O ---------------------------------

class PLYReader:
    @staticmethod
    def read_ply(file_path: str) -> np.ndarray:
        """Read ASCII or binary PLY, return Nx3 float32 array (only x y z)"""
        verts = []
        with open(file_path, 'rb') as f:
            line = f.readline().decode('ascii', errors='ignore').strip()
            if line.lower() != 'ply':
                raise ValueError(f"Not a PLY file: {file_path}")
            fmt = f.readline().decode('ascii', errors='ignore').strip()
            is_binary = 'binary' in fmt

            vertex_count = 0
            prop_count = 0
            while True:
                line = f.readline().decode('ascii', errors='ignore').strip()
                if not line:
                    raise ValueError("Incomplete PLY header")
                if line.startswith('element vertex'):
                    vertex_count = int(line.split()[-1])
                elif line.startswith('property'):
                    prop_count += 1
                elif line == 'end_header':
                    break

            if is_binary:
                # Assume float32 properties & at least first three are xyz
                for _ in range(vertex_count):
                    raw = f.read(4 * prop_count)
                    xyz = struct.unpack('<' + 'f' * prop_count, raw)[:3]
                    verts.append([xyz[0], xyz[1], xyz[2]])
            else:
                for _ in range(vertex_count):
                    parts = f.readline().decode('ascii', errors='ignore').strip().split()
                    x, y, z = map(float, parts[:3])
                    verts.append([x, y, z])

        return np.asarray(verts, dtype=np.float32)


class STLProcessor:
    @staticmethod
    def read_stl(file_path: str):
        """Read binary STL, return (V, F): V vertices (deduplicated) float32, F triangle face indices int32"""
        with open(file_path, 'rb') as f:
            f.read(80)  # header
            tri_n = struct.unpack('<I', f.read(4))[0]
            verts = []
            faces = []
            vmap = {}
            vid = 0
            for _ in range(tri_n):
                f.read(12)  # normal
                tri = []
                for _ in range(3):
                    x = struct.unpack('<f', f.read(4))[0]
                    y = struct.unpack('<f', f.read(4))[0]
                    z = struct.unpack('<f', f.read(4))[0]
                    key = (round(x, 6), round(y, 6), round(z, 6))
                    if key not in vmap:
                        vmap[key] = vid
                        verts.append([x, y, z])
                        vid += 1
                    tri.append(vmap[key])
                faces.append(tri)
                f.read(2)  # attr
        return np.asarray(verts, dtype=np.float32), np.asarray(faces, dtype=np.int32)

    @staticmethod
    def write_stl(file_path: str, V: np.ndarray, F: np.ndarray):
        """Write binary STL. Output triangles according to F, recalculate normals using right-hand rule."""
        with open(file_path, 'wb') as f:
            hdr_text = 'Generated (smoothed) STL'
            header = hdr_text.encode('ascii') + b' ' * (80 - len(hdr_text))
            f.write(header)
            f.write(struct.pack('<I', F.shape[0]))
            for tri in F:
                p0, p1, p2 = V[tri[0]], V[tri[1]], V[tri[2]]
                n = np.cross(p1 - p0, p2 - p0)
                ln = np.linalg.norm(n)
                if ln > 1e-12:
                    n = n / ln
                else:
                    n = np.array([0.0, 0.0, 0.0], dtype=np.float32)
                f.write(struct.pack('<fff', n[0], n[1], n[2]))
                for idx in tri:
                    v = V[idx]
                    f.write(struct.pack('<fff', float(v[0]), float(v[1]), float(v[2])))
                f.write(struct.pack('<H', 0))


# -------------------------- Geometry Utils ---------------------------

def build_kdtree(pts: np.ndarray):
    if cKDTree is None:
        # O(N^2) fallback version (slow with large data)
        class Dummy:
            def __init__(self, data): self.data = data
            def query(self, q):
                q = np.atleast_2d(q)
                idxs = []
                dists = []
                for p in q:
                    d = np.linalg.norm(self.data - p, axis=1)
                    i = int(np.argmin(d))
                    idxs.append(i); dists.append(float(d[i]))
                return np.asarray(dists), np.asarray(idxs, dtype=np.int64)
        return Dummy(pts)
    return cKDTree(pts)


def center_scale(pts: np.ndarray):
    c = pts.mean(axis=0, keepdims=True)
    centered = pts - c
    s = float(np.linalg.norm(centered, axis=1).max())
    if s < 1e-12:
        s = 1.0
    return centered / s, c, s


def umeyama_align(src: np.ndarray, dst: np.ndarray):
    """
    Estimate rigid + scale registration: dst ≈ s * R @ src + t
    Return (R, t, s)
    """
    src_m = src.mean(axis=0)
    dst_m = dst.mean(axis=0)
    X = src - src_m
    Y = dst - dst_m
    C = X.T @ Y / src.shape[0]
    U, S, Vt = np.linalg.svd(C)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    var = (X ** 2).sum() / src.shape[0]
    s = (S.sum() / var) if var > 1e-12 else 1.0
    t = dst_m - s * (R @ src_m)
    return R, t, float(s)


def icp_rigid(src_pts: np.ndarray, dst_pts: np.ndarray, iters: int = 10):
    """
    Iterative Closest Point (rigid + scale), align in normalized space.
    Return (aligned_src_pts, (R_total, t_total, s_total))
    """
    src = src_pts.copy()
    R_total = np.eye(3, dtype=np.float32)
    t_total = np.zeros(3, dtype=np.float32)
    s_total = 1.0
    tree = build_kdtree(dst_pts)
    for _ in range(max(0, iters)):
        _, idx = tree.query(src)
        corr = dst_pts[idx]
        R, t, s = umeyama_align(src, corr)
        src = (s * (R @ src.T).T) + t
        # Accumulate (approximate)
        R_total = R @ R_total
        s_total = s_total * s
        t_total = (R @ t_total) + t
    return src, (R_total, t_total, s_total)


def map_stl_vertices_from_ply(stl_V: np.ndarray, ply_V: np.ndarray, icp_iters: int = 10):
    """
    Align PLY to STL (normalized ICP), then use nearest neighbor to map STL vertices to PLY original coordinate system.
    Return (mapped_V, mean_unit_dist)
    """
    stl_u, stl_c, stl_s = center_scale(stl_V)
    ply_u, ply_c, ply_s = center_scale(ply_V)

    ply_u_aligned, _ = icp_rigid(ply_u, stl_u, iters=icp_iters)

    # Perform nearest neighbor in "STL unit space" (more stable), get PLY corresponding indices, then use original ply_V coordinates
    tree = build_kdtree(ply_u_aligned)
    dists, idx = tree.query(stl_u)
    mapped = ply_V[idx]  # Original coordinates
    return mapped.astype(np.float32), float(np.mean(dists))


# --------------------------- Smoothing -------------------------------

def laplacian_smooth(V: np.ndarray, F: np.ndarray, iters=10, lam=0.3, fix_boundary=True):
    """
    Simple umbrella Laplacian smoothing:
    - iters: number of iterations
    - lam: step size (0~1, smaller is more stable)
    - fix_boundary: try to keep boundary vertices fixed
    """
    V = V.copy()
    n = V.shape[0]
    neighbors = [[] for _ in range(n)]
    edge_count = {}
    for tri in F:
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            neighbors[a].append(b)
            neighbors[b].append(a)
            key = tuple(sorted((int(a), int(b))))
            edge_count[key] = edge_count.get(key, 0) + 1

    boundary = np.zeros(n, dtype=bool)
    if fix_boundary:
        for (a, b), c in edge_count.items():
            if c == 1:
                boundary[a] = True
                boundary[b] = True

    for _ in range(max(0, iters)):
        V_new = V.copy()
        for i in range(n):
            if fix_boundary and boundary[i]:
                continue
            nbrs = neighbors[i]
            if not nbrs:
                continue
            avg = V[nbrs].mean(axis=0)
            V_new[i] = V[i] + lam * (avg - V[i])
        V = V_new
    return V


# --------------------------- Pipeline --------------------------------

VARIANTS = [
    ("original",         "anatomical_{id}_original.ply"),
    ("PC1_minus2sigma",  "anatomical_{id}_PC1_minus2sigma.ply"),
    ("PC1_plus2sigma",   "anatomical_{id}_PC1_plus2sigma.ply"),
    ("PC2_minus2sigma",  "anatomical_{id}_PC2_minus2sigma.ply"),
    ("PC2_plus2sigma",   "anatomical_{id}_PC2_plus2sigma.ply"),
    ("PC3_minus2sigma",  "anatomical_{id}_PC3_minus2sigma.ply"),
    ("PC3_plus2sigma",   "anatomical_{id}_PC3_plus2sigma.ply"),
]


def process_one_model(model_id: str, ply_dir: Path, stl_dir: Path, out_dir: Path,
                      icp_iters=10, smooth_iters=10, smooth_lambda=0.3, fix_boundary=True) -> bool:
    print(f"\n>>> Processing model: {model_id}")
    ply_base = ply_dir / f"anatomical_{model_id}"
    stl_path = stl_dir / f"{model_id}.stl"

    if not ply_base.exists():
        print(f"  ! Missing PLY directory: {ply_base}")
        return False
    if not stl_path.exists():
        print(f"  ! Missing original STL: {stl_path}")
        return False

    # Read STL
    stl_V, stl_F = STLProcessor.read_stl(str(stl_path))
    model_out_dir = out_dir / f"anatomical_{model_id}"
    model_out_dir.mkdir(parents=True, exist_ok=True)

    ok_any = False
    for tag, name_tpl in VARIANTS:
        ply_path = ply_base / name_tpl.format(id=model_id)
        if not ply_path.exists():
            # Missing variant is normal, skip
            continue

        print(f"  - Processing variant: {tag}")
        try:
            ply_V = PLYReader.read_ply(str(ply_path))
            mapped_V, mean_u = map_stl_vertices_from_ply(stl_V, ply_V, icp_iters=icp_iters)
            print(f"    Mapping complete (unit space average distance ~ {mean_u:.6f})")

            # Optional smoothing
            if smooth_iters > 0 and smooth_lambda > 0:
                mapped_V = laplacian_smooth(mapped_V, stl_F,
                                            iters=int(smooth_iters),
                                            lam=float(smooth_lambda),
                                            fix_boundary=bool(fix_boundary))
                print(f"    Smoothing complete (iters={smooth_iters}, λ={smooth_lambda})")

            out_path = model_out_dir / f"anatomical_{model_id}_{tag}.stl"
            STLProcessor.write_stl(str(out_path), mapped_V, stl_F)
            print(f"    ✅ Output written: {out_path}")
            ok_any = True

        except Exception as e:
            print(f"    ❌ Variant {tag} failed: {e}")

    if not ok_any:
        print("  ! No output generated (all variants may be missing or parsing failed)")
    return ok_any
def infer_model_ids(ply_dir: Path, stl_dir: Path):
    ids = set()
    # Use STL filenames as reference
    for stl in stl_dir.glob("*.stl"):
        ids.add(stl.stem)
    # Intersection: must have corresponding PLY directory
    ids2 = []
    for mid in sorted(ids):
        if (ply_dir / f"anatomical_{mid}").exists():
            ids2.append(mid)
    return ids2


def main():
    ap = argparse.ArgumentParser(description="PLY point cloud -> STL (using original STL topology) mapping with optional smoothing (non-destructive)")
    ap.add_argument('--ply_dir', required=True, help='PLY "Final" root directory')
    ap.add_argument('--stl_dir', required=True, help='Original STL root directory')
    ap.add_argument('--output_dir', required=True, help='Output directory (new, does not overwrite original files)')

    ap.add_argument('--models', default='', help='Only process these IDs, comma-separated, e.g.: 001,ABC,case3')
    ap.add_argument('--all', action='store_true', help='Auto-scan and process all (with .stl in stl_dir and corresponding anatomical_<ID> in ply_dir)')

    ap.add_argument('--icp_iters', type=int, default=10, help='ICP iteration count (rigid registration)')
    ap.add_argument('--smooth_iters', type=int, default=10, help='Laplacian smoothing iteration count (0 means no smoothing)')
    ap.add_argument('--smooth_lambda', type=float, default=0.3, help='Laplacian smoothing step size λ (recommended 0.1~0.4)')
    ap.add_argument('--no_fix_boundary', action='store_true', help='Do not fix boundary vertices (default is to fix them)')

    args = ap.parse_args()

    ply_dir = Path(args.ply_dir)
    stl_dir = Path(args.stl_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not ply_dir.exists():
        print(f"PLY root directory does not exist: {ply_dir}")
        sys.exit(1)
    if not stl_dir.exists():
        print(f"STL root directory does not exist: {stl_dir}")
        sys.exit(1)

    # Select IDs to process
    model_ids = []
    if args.models.strip():
        model_ids = [s.strip() for s in args.models.split(',') if s.strip()]
    elif args.all:
        model_ids = infer_model_ids(ply_dir, stl_dir)
    else:
        print("Please use --models to specify IDs, or add --all to process all.")
        sys.exit(1)

    if not model_ids:
        print("No model IDs to process.")
        sys.exit(0)

    print(f"Will process {len(model_ids)} models: {model_ids}")

    ok_count = 0
    for mid in model_ids:
        if process_one_model(mid, ply_dir, stl_dir, out_dir,
                             icp_iters=args.icp_iters,
                             smooth_iters=args.smooth_iters,
                             smooth_lambda=args.smooth_lambda,
                             fix_boundary=not args.no_fix_boundary):
            ok_count += 1

    print(f"\nComplete: {ok_count}/{len(model_ids)} successful. Output directory: {out_dir}")


if __name__ == "__main__":
    main()
