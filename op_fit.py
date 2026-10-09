"""Outline Patch engine - runs in a Python with numpy and scipy.

    python op_fit.py job.json out.json

Lays iso-sections across a disk-like region of a body, so Fusion can loft ONE
smooth face through them. The Form a patch is cut from is many faces (one per
T-spline region); thickened, they give a rim of dozens of tiny edges on which no
fillet builds, and Fusion's Merge Faces refuses them. A single lofted face
through points on the Form has no such seams.

The region is flattened conformally (LSCM, Levy et al. 2002), so the sections
are straight lines in a domain that does not care how far the patch wraps
around the bone - the femur guide's normals span 108 degrees, which rules out a
height field. Section points come back on the TESSELLATION; Fusion then projects
each one exactly onto the body before lofting.

Job (lengths in mm):
    coords, idx      the region's triangles
    inside           per triangle, 1 inside the outline, 0 outside
    outline          the projected outline, closed polyline
    margin_mm        how far past the outline the sheet reaches
    section_mm       spacing between sections
    point_mm         spacing of points along a section
"""

import json
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import cKDTree


def weld(coords, idx, tol=1e-4):
    v = np.asarray(coords, float).reshape(-1, 3)
    f = np.asarray(idx, np.int64).reshape(-1, 3)
    key = np.round(v / tol).astype(np.int64)
    _u, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    v = v[first]
    f = inv[f]
    good = (f[:, 0] != f[:, 1]) & (f[:, 1] != f[:, 2]) & (f[:, 0] != f[:, 2])
    return v, f, good


def components(nv, f):
    """Connected components of faces through shared edges."""
    e = np.sort(np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    fid = np.tile(np.arange(len(f)), 3)
    order = np.lexsort((e[:, 1], e[:, 0]))
    e, fid = e[order], fid[order]
    same = np.all(e[1:] == e[:-1], axis=1)
    a, b = fid[:-1][same], fid[1:][same]
    g = sp.coo_matrix((np.ones(len(a)), (a, b)), shape=(len(f), len(f)))
    n, lab = sp.csgraph.connected_components(g, directed=False)
    return lab


def boundary_loops(f):
    e = np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    s = np.sort(e, axis=1)
    _u, inv, cnt = np.unique(s, axis=0, return_inverse=True, return_counts=True)
    be = e[cnt[inv.ravel()] == 1]
    nxt = {int(a): int(b) for a, b in be}
    loops, seen = [], set()
    for start in list(nxt):
        if start in seen:
            continue
        loop, v = [], start
        while v not in seen and v in nxt:
            seen.add(v)
            loop.append(v)
            v = nxt[v]
        loops.append(loop)
    return loops


def lscm(v, f, pin):
    """Least-squares conformal map. Two pinned vertices fix the similarity."""
    nv = len(v)
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    e1 = b - a
    x = e1 / np.linalg.norm(e1, axis=1, keepdims=True)
    n = np.cross(e1, c - a)
    area2 = np.linalg.norm(n, axis=1)
    n = n / area2[:, None]
    y = np.cross(n, x)
    # local 2D coordinates of the three corners
    p = np.zeros((len(f), 3, 2))
    p[:, 1, 0] = np.linalg.norm(e1, axis=1)
    p[:, 2, 0] = ((c - a) * x).sum(1)
    p[:, 2, 1] = ((c - a) * y).sum(1)
    w = 1.0 / np.sqrt(area2)
    rows, cols, vals = [], [], []
    for k in range(3):
        i, j = (k + 1) % 3, (k + 2) % 3
        # W_k = (x_j - x_i) + i (y_j - y_i), complex weights of corner k
        wr = (p[:, j, 0] - p[:, i, 0]) * w
        wi = (p[:, j, 1] - p[:, i, 1]) * w
        t = np.arange(len(f))
        vk = f[:, k]
        # real part rows 0..F-1, imaginary rows F..2F-1; unknowns u (0..nv-1), v (nv..)
        rows += [t, t, t + len(f), t + len(f)]
        cols += [vk, vk + nv, vk, vk + nv]
        vals += [wr, -wi, wi, wr]
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(2 * len(f), 2 * nv))
    fixed = np.array([pin[0], pin[0] + nv, pin[1], pin[1] + nv])
    fixval = np.array([0.0, 0.0, 1.0, 0.0])
    free = np.setdiff1d(np.arange(2 * nv), fixed)
    rhs = -A[:, fixed] @ fixval
    Af = A[:, free]
    sol = spla.spsolve((Af.T @ Af).tocsc(), Af.T @ rhs)
    uv = np.zeros(2 * nv)
    uv[free] = sol
    uv[fixed] = fixval
    return np.c_[uv[:nv], uv[nv:]]


def main():
    job = json.load(open(sys.argv[1], encoding='utf-8'))
    v, f, good = weld(job['coords'], job['idx'])
    inside = np.asarray(job['inside'], bool)[good]
    f = f[good]
    outline = np.asarray(job['outline'], float).reshape(-1, 3)
    margin = float(job.get('margin_mm', 4.0))
    log = []

    # crop to the inside plus a margin past the outline, one connected piece
    dense = np.vstack([outline[i] + (outline[(i + 1) % len(outline)] - outline[i]) * t
                       for i in range(len(outline)) for t in np.linspace(0, 1, 5, endpoint=False)])
    cen = v[f].mean(1)
    dist = cKDTree(dense).query(cen)[0]
    keep = inside | (dist < margin)
    f = f[keep]
    inside = inside[keep]
    lab = components(len(v), f)
    best = np.bincount(lab[inside]).argmax()
    f = f[lab == best]
    inside_kept = inside[lab == best]
    used = np.unique(f)
    remap = -np.ones(len(v), np.int64)
    remap[used] = np.arange(len(used))
    v = v[used]
    f = remap[f]
    loops = boundary_loops(f)
    loops.sort(key=len, reverse=True)
    log.append('region: %d vertices, %d triangles, %d boundary loop(s)'
               % (len(v), len(f), len(loops)))
    if len(loops) != 1:
        log.append('  more than one boundary: the region has holes; using the outer one')

    # flatten; pin the two boundary vertices farthest apart
    outer = np.array(loops[0])
    bp = v[outer]
    i0 = int(np.argmax(np.linalg.norm(bp - bp.mean(0), axis=1)))
    i1 = int(np.argmax(np.linalg.norm(bp - bp[i0], axis=1)))
    uv = lscm(v, f, (int(outer[i0]), int(outer[i1])))
    # scale the domain to millimetres (median edge ratio) and align u with its long axis
    e = np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    s3 = np.linalg.norm(v[e[:, 0]] - v[e[:, 1]], axis=1)
    s2 = np.linalg.norm(uv[e[:, 0]] - uv[e[:, 1]], axis=1)
    uv *= np.median(s3 / np.maximum(s2, 1e-12))
    uv -= uv.mean(0)
    _w, vec = np.linalg.eigh(np.cov(uv.T))
    uv = uv @ vec[:, ::-1]
    # flips: a conformal map of a disk has none; count them to be sure
    a2 = ((uv[f[:, 1]] - uv[f[:, 0]])[:, 0] * (uv[f[:, 2]] - uv[f[:, 0]])[:, 1] -
          (uv[f[:, 1]] - uv[f[:, 0]])[:, 1] * (uv[f[:, 2]] - uv[f[:, 0]])[:, 0])
    sign = np.sign(np.median(a2))
    flipped = int((a2 * sign <= 0).sum())
    log.append('flattened: domain %.1f x %.1f mm, %d flipped triangles'
               % (np.ptp(uv[:, 0]), np.ptp(uv[:, 1]), flipped))

    # Sections span ONE rectangle of the domain: the inside's extent plus the
    # margin. Clipping each line to the region instead made neighbouring
    # sections jump by up to 14 mm wherever the region is not convex, and the
    # loft folded. Positions come from a smooth tensor B-spline fitted to the
    # whole region, which also carries the sheet a little past the data.
    from scipy.interpolate import BSpline

    # Which domain points lie over the region. Not matplotlib's TriFinder: it
    # rejects a triangulation with any fold, and sliver faces along an old
    # imprint left 6 folded triangles of 5770 in the femur's flattening.
    cuv = uv[f].mean(1)
    ctree = cKDTree(cuv)

    def over_region(q):
        _d, near = ctree.query(q, k=min(12, len(f)))
        hit = np.zeros(len(q), bool)
        for col in range(near.shape[1]):
            t = near[:, col]
            A, B, C = uv[f[t, 0]], uv[f[t, 1]], uv[f[t, 2]]
            den = (B[:, 1] - C[:, 1]) * (A[:, 0] - C[:, 0]) + (C[:, 0] - B[:, 0]) * (A[:, 1] - C[:, 1])
            den = np.where(np.abs(den) < 1e-18, 1e-18, den)
            l1 = ((B[:, 1] - C[:, 1]) * (q[:, 0] - C[:, 0]) + (C[:, 0] - B[:, 0]) * (q[:, 1] - C[:, 1])) / den
            l2 = ((C[:, 1] - A[:, 1]) * (q[:, 0] - C[:, 0]) + (A[:, 0] - C[:, 0]) * (q[:, 1] - C[:, 1])) / den
            hit |= (l1 >= -1e-9) & (l2 >= -1e-9) & (1 - l1 - l2 >= -1e-9)
        return hit
    ins_v = np.unique(f[inside_kept])
    lo = uv[ins_v].min(0) - margin
    hi = uv[ins_v].max(0) + margin
    hc = float(job.get('control_mm', 3.0))
    k = 3

    def knots(a, b):
        n = max(1, int(np.ceil((b - a) / hc)))
        inner = np.linspace(a, b, n + 1)
        return np.r_[[a] * k, inner, [b] * k]
    tu, tv = knots(lo[0], hi[0]), knots(lo[1], hi[1])
    nu, nv = len(tu) - k - 1, len(tv) - k - 1
    sel = np.all((uv >= lo) & (uv <= hi), axis=1)
    Bu = BSpline.design_matrix(uv[sel, 0], tu, k).toarray()
    Bv = BSpline.design_matrix(uv[sel, 1], tv, k).toarray()
    A = (Bu[:, :, None] * Bv[:, None, :]).reshape(len(Bu), nu * nv)

    def d2(n):
        D = np.zeros((max(0, n - 2), n))
        for r in range(n - 2):
            D[r, r:r + 3] = (1.0, -2.0, 1.0)
        return D
    R = np.vstack([np.kron(d2(nu), np.eye(nv)), np.kron(np.eye(nu), d2(nv))])
    lam = float(job.get('fair', 1e-3)) * len(A) / max(1, len(R))
    M = np.vstack([A, np.sqrt(lam) * R])
    rhs = np.vstack([v[sel], np.zeros((len(R), 3))])
    ctrl, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    fit_err = np.linalg.norm(A @ ctrl - v[sel], axis=1)
    log.append('surface fit: %d x %d control points, data within %.3f mm max, %.3f rms'
               % (nu, nv, fit_err.max(), np.sqrt((fit_err ** 2).mean())))

    ds = float(job.get('section_mm', 2.0))
    dp = float(job.get('point_mm', 2.0))
    us = np.linspace(lo[0], hi[0], max(4, int(round((hi[0] - lo[0]) / ds)) + 1))
    vs = np.linspace(lo[1], hi[1], max(4, int(round((hi[1] - lo[1]) / dp)) + 1))
    Bv_s = BSpline.design_matrix(vs, tv, k).toarray()
    sections, on_body = [], []
    for u in us:
        Bu_s = BSpline.design_matrix(np.array([u]), tu, k).toarray()
        Ak = (Bu_s[:, :, None] * Bv_s[:, None, :]).reshape(len(vs), nu * nv)
        P = Ak @ ctrl
        inside_pts = over_region(np.c_[np.full(len(vs), u), vs])
        sections.append(P.tolist())
        on_body.append(inside_pts.astype(int).tolist())
    log.append('%d sections of %d points over %.1f x %.1f mm of domain'
               % (len(sections), len(vs), hi[0] - lo[0], hi[1] - lo[1]))
    json.dump({'ok': True, 'sections': sections, 'on_body': on_body, 'log': log,
               'flipped': flipped},
              open(sys.argv[2], 'w', encoding='utf-8'))


if __name__ == '__main__':
    try:
        main()
    except Exception:
        import traceback
        json.dump({'ok': False, 'error': traceback.format_exc()},
                  open(sys.argv[2], 'w', encoding='utf-8'))
