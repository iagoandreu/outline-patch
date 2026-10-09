"""Outline Patch - the work. Re-imported on every command event and dev job.

WHAT IT DOES
    You draw a closed outline on (or near) a curved body - a 3D spline on a Form,
    fit points snapped to an organic surface, anything. The tool returns the part of
    the body's surface that lies inside the outline, with the body's exact curvature,
    and thickens it into a new body: 3-matic's "attach curve + uniform offset".

WHY NOT IMPRINT THE CURVE
    Split Face with a curve projects the curve onto the faces. A curve that does not
    lie exactly on the surface - and a hand-drawn 3D spline never does - leaves an
    edge that drifts or jumps between faces. Here the curve only steers a BAND:
      1. the outline is sampled every 0.5 mm and each sample projected to its
         exact closest point c(s) on the body (Surface: tessellation for the start,
         Gauss-Newton on the face for the foot), with the outward normal n(s);
      2. the normals are smoothed along the outline;
      3. a spline is fitted through c(s) and its control points are offset by
         -h and +h along the smoothed normal: a ruled band between the two
         offsets, whose mid-line is that spline;
      4. the faces near the outline are split with the band - an exact
         surface/surface intersection - and the side inside the loop is lifted
         out as a surface patch and thickened.
    The cut lands at the closest-point projection of the outline, whatever the
    outline's distance from the body.

All geometry is in Fusion's internal centimetres; reports are in millimetres.
"""

import adsk.core
import adsk.fusion

import json
import math
import os
import sys
import time
import traceback
import types

ADDIN_DIR = os.path.dirname(os.path.realpath(__file__))
CM = 10.0                                   # cm -> mm
MM = 0.1                                    # mm -> cm

CMD_ID = 'OutlinePatch_Create'
CMD_NAME = 'Outline Patch'
RES_DIR = os.path.join(ADDIN_DIR, 'resources', 'OutlinePatch')
SETTINGS_PATH = os.path.join(ADDIN_DIR, 'last_settings.json')
LOG_PATH = os.path.join(ADDIN_DIR, 'OutlinePatch.log')

PANELS = ['SurfaceCreatePanel', 'SolidCreatePanel']

KEEP_SMALL = 'Smaller side of the outline'
KEEP_LARGE = 'Larger side of the outline'
KEEPS = [KEEP_SMALL, KEEP_LARGE]

DIR_OUT = 'Outward (along the surface normal)'
DIR_IN = 'Inward'
DIR_NONE = 'No thickness, surface patch only'
DIRECTIONS = [DIR_OUT, DIR_IN, DIR_NONE]

DEFAULTS = {
    'thickness_mm': 3.7,
    'direction': DIR_OUT,
    'keep': KEEP_SMALL,
    'smooth_mm': 2.0,         # normal smoothing length along the outline
    'curtain_mm': 1.5,        # half-height of the cutting band
    'spacing_mm': 0.5,        # outline sampling
    'rail_mm': 1.0,           # spacing of the band's centre-line points
    'single_face': True,      # loft ONE face over the patch, so the rim can be filleted
    'rim_mm': 2.0,            # fillet on the outer rim of the thickened body; 0 = none
    'min_bend_mm': 0.0,       # ease convex outline corners tighter than this; 0 = untouched
    'sheet_margin_mm': 4.0,   # how far past the outline the single face reaches
    'section_mm': 4.0,        # spacing of the loft sections
    'keep_patch': True,
}

MAX_SAMPLES = 900
MAX_RAIL_POINTS = 600


def log(msg):
    try:
        with open(LOG_PATH, 'a', encoding='utf-8') as fh:
            fh.write('%s  %s\n' % (time.strftime('%H:%M:%S'), msg))
    except Exception:
        pass


def short_error(tb):
    lines = [l for l in tb.strip().splitlines() if l.strip()]
    return lines[-1] if lines else tb


def state():
    st = sys.modules.get('_op_state')
    if st is None:
        st = types.ModuleType('_op_state')
        st.last_report = ''
        sys.modules['_op_state'] = st
    return st


def load_settings():
    o = dict(DEFAULTS)
    try:
        with open(SETTINGS_PATH, encoding='utf-8') as fh:
            o.update(json.load(fh))
    except Exception:
        pass
    return o


def save_settings(o):
    try:
        with open(SETTINGS_PATH, 'w', encoding='utf-8') as fh:
            json.dump({k: o[k] for k in DEFAULTS if k in o}, fh, indent=1)
    except Exception:
        pass


def is_parametric(design):
    return design.designType == adsk.fusion.DesignTypes.ParametricDesignType


# --------------------------------------------------------------------------- #
# small vector helpers (Fusion's Python has no numpy)
# --------------------------------------------------------------------------- #

def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a, b, s=1.0):
    return (a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2])


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a):
    return math.sqrt(_dot(a, a))


def _unit(a):
    n = _norm(a)
    return (a[0] / n, a[1] / n, a[2] / n) if n > 0 else a


def _pt(a):
    return adsk.core.Point3D.create(a[0], a[1], a[2])


def _tup(p):
    return (p.x, p.y, p.z)


# --------------------------------------------------------------------------- #
# the outline
# --------------------------------------------------------------------------- #

def curve_points(curve, tol_cm=0.001):
    ev = curve.worldGeometry.evaluator
    ok, p0, p1 = ev.getParameterExtents()
    ok, pts = ev.getStrokes(p0, p1, tol_cm)
    return [_tup(p) for p in pts]


def chain_outline(curves, gap_cm=0.005):
    """Order the selected sketch curves end to end into ONE closed polyline (cm).

    A single closed spline is its own loop. Several curves are chained greedily by
    their nearest free ends, reversing as needed; the loop must close within gap_cm.
    """
    polys = [curve_points(c) for c in curves if not c.isConstruction]
    if not polys:
        raise ValueError('the outline has no non-construction curves')
    if len(polys) == 1:
        p = polys[0]
        if _norm(_sub(p[0], p[-1])) > gap_cm:
            raise ValueError('the outline is open: its ends are %.3f mm apart'
                             % (_norm(_sub(p[0], p[-1])) * CM))
        return p[:-1]
    loop = list(polys.pop(0))
    while polys:
        tail = loop[-1]
        best = None
        for i, p in enumerate(polys):
            for rev in (False, True):
                d = _norm(_sub(tail, p[-1] if rev else p[0]))
                if best is None or d < best[0]:
                    best = (d, i, rev)
        d, i, rev = best
        if d > gap_cm:
            raise ValueError('the outline curves do not join: a gap of %.3f mm'
                             % (d * CM))
        p = polys.pop(i)
        if rev:
            p = p[::-1]
        loop.extend(p[1:])
    if _norm(_sub(loop[0], loop[-1])) > gap_cm:
        raise ValueError('the outline does not close: its ends are %.3f mm apart'
                         % (_norm(_sub(loop[0], loop[-1])) * CM))
    return loop[:-1]


def resample_closed(pts, spacing_cm, cap):
    """Uniform arc-length resampling of a closed polyline."""
    segs = [_norm(_sub(pts[(i + 1) % len(pts)], pts[i])) for i in range(len(pts))]
    total = sum(segs)
    n = max(24, min(cap, int(round(total / spacing_cm))))
    step = total / n
    out, i, acc = [], 0, 0.0
    for k in range(n):
        target = k * step
        while acc + segs[i] < target and i < len(segs) - 1:
            acc += segs[i]
            i += 1
        t = (target - acc) / segs[i] if segs[i] > 0 else 0.0
        a, b = pts[i], pts[(i + 1) % len(pts)]
        out.append(_add(a, _sub(b, a), t))
    return out, total, step


# --------------------------------------------------------------------------- #
# closest point and outward normal on the body
# --------------------------------------------------------------------------- #

FIT_ENGINE = os.path.join(ADDIN_DIR, 'op_fit.py')
LAST_FIT = []          # job and result files of the last engine run, for inspection
NO_WINDOW = 0x08000000 if os.name == 'nt' else 0


def python_cmd():
    """An interpreter with numpy, scipy and matplotlib (Fusion's has none):
    OUTLINEPATCH_PYTHON, then python_path.txt beside the add-in, then py -3.10."""
    import subprocess
    cands = []
    if os.environ.get('OUTLINEPATCH_PYTHON'):
        cands.append([os.environ['OUTLINEPATCH_PYTHON']])
    try:
        with open(os.path.join(ADDIN_DIR, 'python_path.txt'), encoding='utf-8') as fh:
            cands.append([fh.read().strip()])
    except Exception:
        pass
    cands.append([os.path.expandvars(r'%LOCALAPPDATA%\Programs\Python\Python310\python.exe')])
    cands.append(['py', '-3.10'])
    for c in cands:
        try:
            r = subprocess.run(c + ['-c', 'import numpy, scipy, matplotlib'],
                               capture_output=True, timeout=60, creationflags=NO_WINDOW)
            if r.returncode == 0:
                return c
        except Exception:
            pass
    return None


def run_fit(job):
    import subprocess
    import tempfile
    exe = python_cmd()
    if not exe:
        raise RuntimeError('the single-face option needs a Python with numpy, scipy and '
                           'matplotlib; put its path in python_path.txt beside the add-in')
    d = tempfile.mkdtemp(prefix='outline_patch_')
    jp, op = os.path.join(d, 'job.json'), os.path.join(d, 'out.json')
    with open(jp, 'w', encoding='utf-8') as fh:
        json.dump(job, fh)
    r = subprocess.run(exe + [FIT_ENGINE, jp, op], capture_output=True, text=True,
                       timeout=600, creationflags=NO_WINDOW)
    if not os.path.exists(op):
        raise RuntimeError('the fit engine produced nothing:\n%s' % r.stderr[-600:])
    with open(op, encoding='utf-8') as fh:
        res = json.load(fh)
    LAST_FIT[:] = [jp, op]
    if not res.get('ok'):
        raise RuntimeError('the fit engine failed:\n%s' % res.get('error', '')[-600:])
    return res


def tessellate_faces(faces, tol_cm=0.005):
    """Triangles of the given faces, in mm: coords, idx and the face index of each."""
    coords, idx, owner = [], [], []
    for k, f in enumerate(faces):
        calc = f.meshManager.createMeshCalculator()
        calc.surfaceTolerance = tol_cm
        m = calc.calculate()
        base = len(coords) // 3
        coords.extend(v * CM for v in m.nodeCoordinatesAsDouble)
        ni = m.nodeIndices
        idx.extend(base + i for i in ni)
        owner.extend([k] * (len(ni) // 3))
    return coords, idx, owner


def single_face_sheet(comp, body, side, outside_ring, cps, o, notes, created):
    """ONE lofted face over the patch region, its points exactly on the body.

    Returns (sheet body, its deviation from the body in mm). The loft and its
    sketch are appended to `created`, to be removed once the patch is lifted out.
    """
    faces = list(side) + list(outside_ring)
    coords, idx, owner = tessellate_faces(faces)
    inside = [1 if k < len(side) else 0 for k in owner]
    res = run_fit({'coords': coords, 'idx': idx, 'inside': inside,
                   'outline': [[c * CM for c in q] for q in cps],
                   'margin_mm': float(o.get('sheet_margin_mm', 4.0)),
                   'section_mm': float(o.get('section_mm', 2.0)),
                   'point_mm': float(o.get('section_mm', 2.0))})
    notes.extend('fit: ' + l for l in res.get('log', []))
    secs = [[(q[0] * MM, q[1] * MM, q[2] * MM) for q in sec] for sec in res['sections']]
    flags = res['on_body']
    flat = [q for sec, fl in zip(secs, flags) for q, b in zip(sec, fl) if b]
    proj = Surface(comp, body, flat, 0.05)
    # points over the region go exactly onto the body; the few past it keep the
    # fitted surface's smooth extension (they are cut away)
    exact = [[proj.closest(q)[1] if b else q for q, b in zip(sec, fl)]
             for sec, fl in zip(secs, flags)]

    sk = comp.sketches.add(comp.xYConstructionPlane)
    sk.name = 'Outline Patch sections'
    created.append(sk)
    sk.isComputeDeferred = True
    splines = []
    for sec in exact:
        coll = adsk.core.ObjectCollection.create()
        for q in sec:
            coll.add(sk.modelToSketchSpace(_pt(q)))
        splines.append(sk.sketchCurves.sketchFittedSplines.add(coll))
    sk.isComputeDeferred = False
    lin = comp.features.loftFeatures.createInput(
        adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    lin.isSolid = False
    for spl in splines:
        lin.loftSections.add(comp.features.createPath(spl, False))
    loft = comp.features.loftFeatures.add(lin)
    created.append(loft)
    sheet = loft.bodies.item(0)
    sheet.name = 'Outline Patch sheet'

    return sheet, None



def deviation_from(comp, face_body, body, n_max=3000):
    """How far a body's faces stray from another body, sampled on its tessellation."""
    calc = face_body.meshManager.createMeshCalculator()
    calc.surfaceTolerance = 0.002
    m = calc.calculate()
    xyz = m.nodeCoordinatesAsDouble
    pts = [(xyz[i], xyz[i + 1], xyz[i + 2]) for i in range(0, len(xyz), 3)]
    pts = pts[::max(1, len(pts) // n_max)]
    proj = Surface(comp, body, pts, 0.05)
    devs = []
    for q in pts:
        try:
            devs.append(proj.closest(q)[0] * CM)
        except ValueError:
            pass
    devs.sort()
    if not devs:
        return None
    return {'max_mm': round(devs[-1], 4), 'p99_mm': round(devs[int(0.99 * (len(devs) - 1))], 4),
            'mean_mm': round(sum(devs) / len(devs), 4), 'points': len(devs)}


def _closest_on_triangle(p, a, b, c):
    """Closest point to p on triangle abc (Ericson, Real-Time Collision Detection)."""
    ab, ac, ap = _sub(b, a), _sub(c, a), _sub(p, a)
    d1, d2 = _dot(ab, ap), _dot(ac, ap)
    if d1 <= 0 and d2 <= 0:
        return a
    bp = _sub(p, b)
    d3, d4 = _dot(ab, bp), _dot(ac, bp)
    if d3 >= 0 and d4 <= d3:
        return b
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        return _add(a, ab, d1 / (d1 - d3))
    cp = _sub(p, c)
    d5, d6 = _dot(ab, cp), _dot(ac, cp)
    if d6 >= 0 and d5 <= d6:
        return c
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        return _add(a, ac, d2 / (d2 - d6))
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        return _add(b, _sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6)))
    den = 1.0 / (va + vb + vc)
    return _add(_add(a, ab, vb * den), ac, vc * den)


class Surface:
    """Exact closest point + outward normal on a BRep body.

    MeasureManager.measureMinimumDistance with a temporary Point3D is NOT usable:
    against a sphere it returned distance 0 and the query point as the foot for a
    point 0.4 mm off. SurfaceEvaluator.getParameterAtPoint is not a projection
    either - from 0.5 mm off a Form's NURBS face it converged 94 mm away. So:
      1. tessellate only the faces near the outline, at 0.05 mm,
      2. nearest triangle through a hash grid, giving a foot within ~0.05 mm on
         a known face,
      3. Gauss-Newton on that face's exact surface in (u, v) to the true foot.
    getParameterAtPoint is only ever asked about a point ~on the surface.
    """

    TESS_TOL = 0.005          # cm; 0.001 on a 30 mm ball once killed Fusion
    SLACK = 0.01              # cm; faces this much farther than the nearest are tried too

    def __init__(self, comp, body, outline, step):
        self.comp = comp
        self.body = body
        self.solid = body.isSolid
        self.sign = {}
        self.faces_hit = {}
        self.searches = 0
        faces = {}
        stride = max(1, int(round(0.2 / step)))          # a probe every ~2 mm
        for p in outline[::stride]:
            found = comp.findBRepUsingPoint(
                _pt(p), adsk.fusion.BRepEntityTypes.BRepFaceEntityType, 0.5, False)
            for k in range(found.count):
                f = adsk.fusion.BRepFace.cast(found.item(k))
                if f and f.body == self.body:
                    faces[f.tempId] = f
        if not faces:
            raise ValueError('the outline is more than 5 mm from the body')
        self.tris = []
        for f in faces.values():
            calc = f.meshManager.createMeshCalculator()
            calc.surfaceTolerance = self.TESS_TOL
            m = calc.calculate()
            xyz = m.nodeCoordinatesAsDouble
            idx = m.nodeIndices
            P = [(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]) for i in range(len(xyz) // 3)]
            for t in range(0, len(idx), 3):
                self.tris.append((P[idx[t]], P[idx[t + 1]], P[idx[t + 2]], f))
        self.n_faces = len(faces)
        self.cell = 0.1
        self.grid = {}
        c = self.cell
        for ti, (a, b, cc, _f) in enumerate(self.tris):
            lo = [int(math.floor(min(a[k], b[k], cc[k]) / c)) for k in range(3)]
            hi = [int(math.floor(max(a[k], b[k], cc[k]) / c)) for k in range(3)]
            for x in range(lo[0], hi[0] + 1):
                for y in range(lo[1], hi[1] + 1):
                    for z in range(lo[2], hi[2] + 1):
                        self.grid.setdefault((x, y, z), []).append(ti)

    def _nearest_triangle(self, p):
        c = self.cell
        k0 = [int(math.floor(p[k] / c)) for k in range(3)]
        best, seen, per_face = None, set(), {}
        for ring in range(0, 51):
            for x in range(k0[0] - ring, k0[0] + ring + 1):
                for y in range(k0[1] - ring, k0[1] + ring + 1):
                    for z in range(k0[2] - ring, k0[2] + ring + 1):
                        if max(abs(x - k0[0]), abs(y - k0[1]), abs(z - k0[2])) != ring:
                            continue
                        for ti in self.grid.get((x, y, z), ()):
                            if ti in seen:
                                continue
                            seen.add(ti)
                            a, b, cc, f = self.tris[ti]
                            q = _closest_on_triangle(p, a, b, cc)
                            d = _norm(_sub(p, q))
                            cur = per_face.get(f.tempId)
                            if cur is None or d < cur[0]:
                                per_face[f.tempId] = (d, q, f)
                            if best is None or d < best:
                                best = d
            # every triangle closer than `ring` cells has been seen by now
            if best is not None and best + self.SLACK <= ring * c:
                break
        if best is None:
            return []
        return sorted([v for v in per_face.values() if v[0] <= best + self.SLACK],
                      key=lambda v: v[0])

    def _refine(self, face, p, q0):
        """Gauss-Newton foot of p on the face's exact surface, starting at q0."""
        ev = face.evaluator
        ok, uv = ev.getParameterAtPoint(_pt(q0))
        if not ok:
            return q0, None
        for _ in range(8):
            ok, S = ev.getPointAtParameter(uv)
            ok2, du, dv = ev.getFirstDerivative(uv)
            if not (ok and ok2):
                break
            S, du, dv = _tup(S), _tup(du), _tup(dv)
            r = _sub(p, S)
            a, b, c = _dot(du, du), _dot(du, dv), _dot(dv, dv)
            g1, g2 = _dot(du, r), _dot(dv, r)
            det = a * c - b * b
            if abs(det) < 1e-18:
                break
            su, sv = (c * g1 - b * g2) / det, (a * g2 - b * g1) / det
            nxt = adsk.core.Point2D.create(uv.x + su, uv.y + sv)
            if not ev.isParameterOnFace(nxt):
                break
            uv = nxt
            if abs(su) * math.sqrt(a) + abs(sv) * math.sqrt(c) < 1e-7:
                break
        ok, S = ev.getPointAtParameter(uv)
        return (_tup(S) if ok else q0), uv

    def _face_sign(self, face, q, n):
        key = face.tempId
        if key not in self.sign:
            s = -1.0 if face.isParamReversed else 1.0
            if self.solid:
                # measured, not assumed: step out along s*n and ask the solid
                probe = _add(q, n, 0.01 * s)
                pc = self.body.pointContainment(_pt(probe))
                if pc == adsk.fusion.PointContainment.PointInsidePointContainment:
                    s = -s
            self.sign[key] = s
        return self.sign[key]

    def closest(self, p):
        self.searches += 1
        cands = self._nearest_triangle(p)
        if not cands:
            raise ValueError('an outline point is more than 5 mm from the body')
        # Refine on each nearby face and keep the TRUE foot: p - q parallel to the
        # normal. The refined foot can be farther from p than the tessellation
        # foot (from inside a convex body the chords lie in between), so distance
        # alone cannot pick it; near a face edge the step leaves the first face and
        # the foot belongs to its neighbour.
        best = None
        for _d0, q0, face in cands:
            q, uv = self._refine(face, p, q0)
            ev = face.evaluator
            if uv is None:
                ok, uv = ev.getParameterAtPoint(_pt(q))
            ok, n = ev.getNormalAtParameter(uv)
            if not ok:
                continue
            n = _unit(_tup(n))
            r = _sub(p, q)
            lr = _norm(r)
            ang = 0.0 if lr < 1e-5 else math.degrees(
                math.acos(min(1.0, abs(_dot(r, n)) / lr)))
            score = (ang > 1.0, lr if ang <= 1.0 else ang)
            if best is None or score < best[0]:
                best = (score, face, q, n, ang)
            if ang <= 1.0 and len(cands) == 1:
                break
        if best is None:
            raise ValueError('could not evaluate the surface normal at an outline point')
        _score, face, q, n, ang = best
        self.worst_angle = max(getattr(self, 'worst_angle', 0.0), ang)
        if ang > 1.0:
            self.loose_feet = getattr(self, 'loose_feet', 0) + 1
        s = self._face_sign(face, q, n)
        self.faces_hit[face.tempId] = face
        return _norm(_sub(p, q)), q, (n[0] * s, n[1] * s, n[2] * s)


def bend_curvature(cps, nrm, step):
    """Signed curvature of the outline IN the surface (geodesic), per sample, 1/cm.
    Positive bends towards the inside of the loop (a convex corner of the patch),
    negative away from it (a notch)."""
    n = len(cps)
    t = [_unit(_sub(cps[(i + 1) % n], cps[i - 1])) for i in range(n)]
    kg = []
    for i in range(n):
        dt = _sub(t[(i + 1) % n], t[i - 1])
        b = _unit(_cross(nrm[i], t[i]))       # in the surface, to the left of travel
        kg.append(_dot(dt, b) / (2.0 * step))
    # which side is inside: the loop's total turning is positive walking with the
    # inside on the left
    if sum(kg) < 0:
        kg = [-k for k in kg]
    return kg


def ease_bends(surf, cps, nrm, step, r_min, max_iter=400):
    """Relax only the CONVEX corners tighter than r_min, in the surface, keeping every
    point on the body. Returns (points, normals, largest displacement in cm).

    A rim fillet's spine runs r inside the rim: R - r round a convex corner of
    radius R, which goes negative when R < r and the fillet cannot build;
    R + r round a notch, which is always fine. Measured on the femur guide:
    tightest convex corner 1.66 mm, rims of 1.8 mm built and 1.9 mm did not;
    easing the notches instead changed nothing."""
    pts, nrs = list(cps), list(nrm)
    orig = list(cps)
    n = len(pts)
    for _ in range(max_iter):
        kg = bend_curvature(pts, nrs, step)
        tight = [i for i in range(n) if kg[i] > 1.0 / r_min]
        if not tight:
            break
        move = set()
        for i in tight:
            for j in range(-3, 4):
                move.add((i + j) % n)
        for i in move:
            lap = _sub(_add(pts[i - 1], pts[(i + 1) % n]), _add(pts[i], pts[i]))
            lap = _sub(lap, _add((0, 0, 0), nrs[i], _dot(lap, nrs[i])))   # tangential only
            _d, q, nn = surf.closest(_add(pts[i], lap, 0.25))
            pts[i], nrs[i] = q, nn
    moved = max(_norm(_sub(a, b)) for a, b in zip(pts, orig))
    return pts, nrs, moved


def smooth_normals(normals, step_cm, sigma_cm):
    """Gaussian smoothing of unit normals along a closed outline, renormalised."""
    n = len(normals)
    if sigma_cm <= 0:
        return list(normals)
    w = max(1, int(math.ceil(3.0 * sigma_cm / step_cm)))
    w = min(w, n // 2)
    kern = [math.exp(-0.5 * (k * step_cm / sigma_cm) ** 2) for k in range(-w, w + 1)]
    out = []
    for i in range(n):
        acc = (0.0, 0.0, 0.0)
        for j, kw in zip(range(-w, w + 1), kern):
            acc = _add(acc, normals[(i + j) % n], kw)
        out.append(_unit(acc))
    return out


# --------------------------------------------------------------------------- #
# the cutting band
# --------------------------------------------------------------------------- #

def _closed_spline(sk, pts):
    coll = adsk.core.ObjectCollection.create()
    for p in pts:
        coll.add(sk.modelToSketchSpace(_pt(p)))
    spl = sk.sketchCurves.sketchFittedSplines.add(coll)
    spl.isClosed = True
    return spl


def _offset_nurbs(g, normals_at, h):
    """Two copies of the NURBS curve g with every control point moved by -h / +h
    along the normal near it. Same degree, knots and weights: equal parameters
    correspond exactly, and their mid-curve IS g (a B-spline is affine in its
    control points)."""
    ok, ctrl, deg, knots, rational, weights, periodic = g.getData()
    inner, outer = [], []
    for c in ctrl:
        n = normals_at(_tup(c))
        inner.append(_pt(_add(_tup(c), n, -h)))
        outer.append(_pt(_add(_tup(c), n, h)))
    if rational:
        mk = lambda pts: adsk.core.NurbsCurve3D.createRational(pts, deg, knots, weights, periodic)
    else:
        mk = lambda pts: adsk.core.NurbsCurve3D.createNonRational(pts, deg, knots, periodic)
    return mk(inner), mk(outer)


def build_curtain(comp, centre, normals_at, h, notes):
    """The cutting band: a ruled surface between the centre spline offset -h and +h
    along the smoothed normals. Returns (temp band, timeline items to remove).

    The rails are NOT two separately fitted splines: each would be parametrised
    by its own normalised chord length, and since offsetting along the normal of
    a curved surface changes the length, equal parameters drift apart along the
    loop. On the femur the band then crossed the surface up to 0.2 mm off the
    outline. Offsetting the CONTROL POINTS of one spline keeps the knots, so the
    band's mid-line is the centre spline itself, through every projected point.
    """
    sk = comp.sketches.add(comp.xYConstructionPlane)
    sk.name = 'Outline Patch centre'
    try:
        sk.isComputeDeferred = True
        spl = _closed_spline(sk, centre)
        sk.isComputeDeferred = False
        a, b = _offset_nurbs(spl.worldGeometry, normals_at, h)
        tbm = adsk.fusion.TemporaryBRepManager.get()
        wa, _ea = tbm.createWireFromCurves([a])
        wb, _eb = tbm.createWireFromCurves([b])
        if wa is None or wb is None:
            raise RuntimeError('createWireFromCurves returned nothing')
        band = tbm.createRuledSurface(wa.wires.item(0), wb.wires.item(0))
        if band is None:
            raise RuntimeError('createRuledSurface returned nothing')
    except Exception:
        _delete_quietly([sk], notes)
        raise RuntimeError('could not build the cutting band: %s'
                           % short_error(traceback.format_exc()))
    notes.append('band: ruled between control-point offsets of one spline')
    return band, [sk]


# --------------------------------------------------------------------------- #
# the pipeline
# --------------------------------------------------------------------------- #

class _Grid:
    """Point-to-polyline distance through a uniform hash grid (cm)."""

    def __init__(self, pts, cell):
        self.pts = pts
        self.cell = cell
        self.g = {}
        n = len(pts)
        for i in range(n):
            a, b = pts[i], pts[(i + 1) % n]
            for p in (a, b, _add(a, _sub(b, a), 0.5)):
                self.g.setdefault(self._key(p), set()).add(i)

    def _key(self, p):
        c = self.cell
        return (int(math.floor(p[0] / c)), int(math.floor(p[1] / c)), int(math.floor(p[2] / c)))

    def dist(self, p):
        kx, ky, kz = self._key(p)
        best = None
        n = len(self.pts)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for i in self.g.get((kx + dx, ky + dy, kz + dz), ()):
                        a, b = self.pts[i], self.pts[(i + 1) % n]
                        ab = _sub(b, a)
                        L2 = _dot(ab, ab)
                        t = 0.0 if L2 == 0 else max(0.0, min(1.0, _dot(_sub(p, a), ab) / L2))
                        d = _norm(_sub(p, _add(a, ab, t)))
                        if best is None or d < best:
                            best = d
        return best if best is not None else float('inf')


def _edge_points(edge, k=5):
    ev = edge.evaluator
    ok, p0, p1 = ev.getParameterExtents()
    out = []
    for i in range(k):
        ok, p = ev.getPointAtParameter(p0 + (p1 - p0) * (i + 0.5) / k)
        if ok:
            out.append(_tup(p))
    return out


def _dist_to_faces(faces, p):
    """Distance from a point NEAR a face to it, through the face's evaluator: only
    well-posed for points close to the surface, which is the only use here."""
    best = float('inf')
    for f in faces:
        ev = f.evaluator
        ok, uv = ev.getParameterAtPoint(_pt(p))
        if not ok:
            continue
        ok, q = ev.getPointAtParameter(uv)
        if ok:
            best = min(best, _norm(_sub(p, _tup(q))))
    return best


def _regions(body, boundary):
    """Union-find the faces across every edge that is NOT on the cut."""
    parent = {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    faces = list(body.faces)
    for f in faces:
        parent[f.tempId] = f.tempId
    for e in body.edges:
        if e.tempId in boundary:
            continue
        fs = [f.tempId for f in e.faces]
        for x in fs[1:]:
            ra, rb = find(fs[0]), find(x)
            if ra != rb:
                parent[ra] = rb
    groups = {}
    for f in faces:
        groups.setdefault(find(f.tempId), []).append(f)
    return list(groups.values())


def _collection(items):
    coll = adsk.core.ObjectCollection.create()
    for it in items:
        coll.add(it)
    return coll


def _delete_quietly(items, notes):
    for it in items:
        try:
            if it is not None and it.isValid:
                it.deleteMe()
        except Exception:
            notes.append('cleanup: %s' % short_error(traceback.format_exc()))


def make_patch(app, curves, body, o):
    """The whole job. Returns a report dict; raises with a readable message."""
    t0 = time.time()
    design = adsk.fusion.Design.cast(app.activeProduct)
    parametric = is_parametric(design)
    if not adsk.fusion.BRepBody.cast(body):
        raise ValueError('pick a solid or surface body; a mesh body has to be '
                         'converted first (Mesh to Part)')
    if body.assemblyContext is not None:
        raise ValueError('pick the body in its own component (activate the '
                         'component that owns it)')
    comp = body.parentComponent
    notes = []

    # 1. the outline, resampled
    raw = chain_outline(curves)
    spacing = float(o['spacing_mm']) * MM
    pts, length, step = resample_closed(raw, spacing, MAX_SAMPLES)

    tm = {}
    t1 = time.time()
    # 2. closest point and outward normal at every sample
    surf = Surface(comp, body, pts, step)
    tm['tessellate'] = round(time.time() - t1, 2)
    t1 = time.time()
    off, cps, nrm = [], [], []
    for p in pts:
        d, q, n = surf.closest(p)
        off.append(d)
        cps.append(q)
        nrm.append(n)
    tm['project'] = round(time.time() - t1, 2)
    kg0 = bend_curvature(cps, nrm, step)
    conv = [k for k in kg0 if k > 0]
    bend_before = (1.0 / max(conv)) * CM if conv else None
    eased = 0.0
    if float(o.get('min_bend_mm', 0.0)) > 0:
        cps, nrm, eased = ease_bends(surf, cps, nrm, step, float(o['min_bend_mm']) * MM)
    # 3. normals smoothed along the outline: the "interpolate the path" step
    nrm_s = smooth_normals(nrm, step, float(o['smooth_mm']) * MM)
    turn = max(math.degrees(math.acos(max(-1.0, min(1.0, _dot(a, b)))))
               for a, b in zip(nrm, nrm_s))

    # 4. the band's centre line: a spline through the projected points
    h = float(o['curtain_mm']) * MM
    # Between its points the centre is a spline through the projected outline.
    # Measured on a 354 mm femur outline: 2 mm spacing left the cut 0.058 mm off
    # the projected outline in the tightest turn (radius ~2 mm), 1 mm 0.017 mm,
    # which is the 0.5 mm sampling's own chord error; the band then took 7 s.
    stride = max(1, int(round(float(o.get('rail_mm', 1.0)) * MM / step)),
                 int(math.ceil(len(cps) / float(MAX_RAIL_POINTS))))
    centre = cps[::stride]

    def normals_at(q):
        # the smoothed normal of the projected sample nearest to q
        best, bi = None, 0
        for i0 in range(len(cps)):
            d = _norm(_sub(q, cps[i0]))
            if best is None or d < best:
                best, bi = d, i0
        return nrm_s[bi]

    if o.get('dump'):
        with open(o['dump'], 'w', encoding='utf-8') as fh:
            json.dump({'samples': pts, 'closest': cps, 'normals': nrm,
                       'normals_smoothed': nrm_s, 'centre': centre,
                       'offsets': off}, fh)

    created = []          # timeline entries to remove at the end
    try:
        t1 = time.time()
        band, cleanup = build_curtain(comp, centre, normals_at, h, notes)
        tm['band'] = round(time.time() - t1, 2)
        t1 = time.time()
        created.extend(cleanup)

        # 5. a working copy of the body (never the user's) and the band, in a base feature
        tbm = adsk.fusion.TemporaryBRepManager.get()
        work_tmp = tbm.copy(body)
        # every face within 0.2 mm of the projected outline, not only the faces
        # the samples landed on: thin faces along the outline (an old imprint, a
        # sliver of the Form) are crossed between samples and left gaps in the cut
        near_faces = dict(surf.faces_hit)
        for c in cps:
            found = comp.findBRepUsingPoint(
                _pt(c), adsk.fusion.BRepEntityTypes.BRepFaceEntityType, 0.02, False)
            for k in range(found.count):
                f = adsk.fusion.BRepFace.cast(found.item(k))
                if f and f.body == body:
                    near_faces[f.tempId] = f
        hit_pts = [f.pointOnFace for f in near_faces.values()]
        base = comp.features.baseFeatures.add() if parametric else None
        if base:
            base.startEdit()
        try:
            work = comp.bRepBodies.add(work_tmp, base) if base else comp.bRepBodies.add(work_tmp)
            tool = comp.bRepBodies.add(band, base) if base else comp.bRepBodies.add(band)
        finally:
            if base:
                base.finishEdit()
        if base:
            created.insert(0, base)
            work = base.bodies.item(0)
            tool = base.bodies.item(1)
        work.name = 'Outline Patch work'

        # 6. split only the faces the outline crosses, with the band
        targets = adsk.core.ObjectCollection.create()
        for f in work.faces:
            q = f.pointOnFace
            if any(q.distanceTo(h0) < 1e-7 for h0 in hit_pts):
                targets.add(f)
        if targets.count == 0:
            raise RuntimeError('none of the copied faces matched the faces the '
                               'outline lies on')
        sfi = comp.features.splitFaceFeatures.createInput(targets, tool, False)
        split = comp.features.splitFaceFeatures.add(sfi)
        tm['split'] = round(time.time() - t1, 2)
        t1 = time.time()
        created.append(split)
        work = split.bodies.item(0) if split.bodies.count else work

        # 7. the cut = the edges lying on the band. Two signals, either one
        # enough: the band's own face (getParameterAtPoint is unreliable at a few
        # points even this close - edges 0.016 mm from the outline read 0.1-1 mm
        # off the band) and the projected outline polyline. An old edge that
        # coincides with the cut (an earlier imprint of a similar curve) is a cut
        # edge as much as a new one; a sliver between two such chains is dropped
        # by the area rule below.
        grid = _Grid(cps, max(4.0 * step, 0.2))
        band_faces = list(tool.faces)
        boundary, bmax = set(), 0.0
        for e in work.edges:
            pts_e = _edge_points(e)
            if not pts_e:
                continue
            dline = max(grid.dist(q) for q in pts_e)
            if dline > 0.03:
                continue
            dband = max(_dist_to_faces(band_faces, q) for q in pts_e)
            if dband < 0.002 or dline < 0.002:
                boundary.add(e.tempId)
                bmax = max(bmax, dline)
        if not boundary:
            raise RuntimeError('the band did not cut the body; raise the band '
                               'height or check the outline lies on the body')
        groups = _regions(work, boundary)
        area = lambda g: sum(f.area for f in g)
        total = sum(area(g) for g in groups)
        big = sorted([g for g in groups if area(g) > 0.01 * total], key=area)
        if len(big) < 2:
            # where the open ends of the cut are says why it did not close
            ends = {}
            for e in work.edges:
                if e.tempId in boundary:
                    for v in (e.startVertex, e.endVertex):
                        ends[v.tempId] = ends.get(v.tempId, 0) + 1
            loose = [vid for vid, k in ends.items() if k % 2]
            raise RuntimeError('the cut does not close into a loop on the body '
                               '(%d cut edges, %d loose ends, regions mm2: %s)'
                               % (len(boundary), len(loose),
                                  ', '.join('%.1f' % (area(g) * 100) for g in
                                            sorted(groups, key=area)[-4:])))
        side = big[0] if o['keep'] == KEEP_SMALL else big[-1]
        side_ids = set(f.tempId for f in side)
        keep_pts = [f.pointOnFace for f in side]

        sheet_dev = None
        want_single = (o.get('single_face', True) and float(o.get('rim_mm', 0)) > 0
                       and o['direction'] != DIR_NONE)
        if want_single and o['keep'] == KEEP_SMALL:
            # one smooth face instead of the Form's many: the faces across the cut
            # from the patch give the sheet its margin
            ring = {}
            for e in work.edges:
                if e.tempId in boundary:
                    for f in e.faces:
                        if f.tempId not in side_ids:
                            ring[f.tempId] = f
            sheet, sheet_dev = single_face_sheet(comp, body, side, ring.values(), cps,
                                                 o, notes, created)
            tm['single_face'] = round(time.time() - t1, 2)
            t1 = time.time()
            # On ONE face, with a curve that lies on it, imprinting is reliable,
            # and it gives the patch a single closed edge. Cut by the band, the
            # same sheet came out with a 109-edge boundary, a 109-face side wall
            # after the thicken, and a rim no fillet would build on.
            centre_spl = cleanup[0].sketchCurves.item(0)
            try:
                sfi2 = comp.features.splitFaceFeatures.createInput(
                    _collection(list(sheet.faces)), centre_spl, False)
                sfi2.setClosestPointSplitType()
                split2 = comp.features.splitFaceFeatures.add(sfi2)
                notes.append('single face cut by the centre spline')
            except Exception:
                notes.append('centre-spline cut refused (%s), cutting with the band'
                             % short_error(traceback.format_exc()))
                sfi2 = comp.features.splitFaceFeatures.createInput(
                    _collection(list(sheet.faces)), tool, False)
                split2 = comp.features.splitFaceFeatures.add(sfi2)
            created.append(split2)
            sheet = split2.bodies.item(0) if split2.bodies.count else sheet
            # the patch is the face that reaches no open edge of the sheet
            inner = [f for f in sheet.faces
                     if all(e.faces.count == 2 for e in f.edges)]
            if len(inner) != 1:
                raise RuntimeError('cutting the single face left %d inner faces'
                                   % len(inner))
            work = sheet
            side = inner
            side_ids = set(f.tempId for f in side)
            keep_pts = [f.pointOnFace for f in side]

        tm['cut_and_regions'] = round(time.time() - t1, 2)
        t1 = time.time()
        # 8. lift the kept faces out as a surface patch
        patch_tmp = tbm.copy(work)
        drop = [f for f in patch_tmp.faces
                if not any(f.pointOnFace.distanceTo(q) < 1e-7 for q in keep_pts)]
        if len(list(patch_tmp.faces)) - len(drop) != len(side_ids):
            raise RuntimeError('matched %d of %d patch faces in the copy'
                               % (len(list(patch_tmp.faces)) - len(drop), len(side_ids)))
        if drop and not tbm.deleteFaces(drop, True):
            raise RuntimeError('could not trim the patch out of the copy')
    finally:
        # nothing of the working state stays in the design, whatever happened
        for it in reversed(created):
            _delete_quietly([it], notes)
        created = []

    tm['extract_and_cleanup'] = round(time.time() - t1, 2)
    t1 = time.time()
    pbase = comp.features.baseFeatures.add() if parametric else None
    if pbase:
        pbase.startEdit()
    try:
        patch = comp.bRepBodies.add(patch_tmp, pbase) if pbase else comp.bRepBodies.add(patch_tmp)
    finally:
        if pbase:
            pbase.finishEdit()
    if pbase:
        pbase.name = 'Outline Patch'
        patch = pbase.bodies.item(0)
    stem = curves[0].parentSketch.name if curves[0].parentSketch else 'Outline'
    patch.name = '%s patch' % stem

    rep = {'outline_mm': round(length * CM, 2), 'samples': len(pts),
           'spacing_mm': round(step * CM, 3),
           'outline_to_body_max_mm': round(max(off) * CM, 3),
           'outline_to_body_mean_mm': round(sum(off) / len(off) * CM, 3),
           'normal_smoothing_turn_max_deg': round(turn, 2),
           'tightest_convex_bend_mm': round(bend_before, 2) if bend_before else None,
           'outline_eased_mm': round(eased * CM, 3),
           'faces_crossed': len(surf.faces_hit), 'faces_split': len(near_faces),
           'faces_tessellated': surf.n_faces,
           'foot_angle_worst_deg': round(getattr(surf, 'worst_angle', 0.0), 3),
           'cut_edges': len(boundary),
           'cut_to_outline_max_mm': round(bmax * CM, 4),
           'patch': patch.name, 'patch_faces': patch.faces.count,
           'patch_area_mm2': round(patch.area * 100.0, 2),
           'patch_edges': patch.edges.count,
           'single_face_dev': (deviation_from(comp, patch, body)
                               if sheet_dev is not None or (want_single and patch.faces.count == 1)
                               else None)}

    tm['patch_body'] = round(time.time() - t1, 2)
    t1 = time.time()
    # 9. thicken into a new body
    if o['direction'] != DIR_NONE and float(o['thickness_mm']) > 0:
        t = float(o['thickness_mm']) * MM * (1.0 if o['direction'] == DIR_OUT else -1.0)
        # a lofted face has its own normal: make "outward" the body's outward
        pf = patch.faces.item(0)
        q = _tup(pf.pointOnFace)
        ok, uv = pf.evaluator.getParameterAtPoint(_pt(q))
        ok, npf = pf.evaluator.getNormalAtParameter(uv)
        npf = _unit(_tup(npf))
        if pf.isParamReversed:
            npf = (-npf[0], -npf[1], -npf[2])
        _d, _q, nbody = Surface(comp, body, [q], 0.05).closest(q)
        if _dot(npf, nbody) < 0:
            t = -t
        faces = adsk.core.ObjectCollection.create()
        for f in patch.faces:
            faces.add(f)
        before = set(b.entityToken for b in comp.bRepBodies)
        tin = comp.features.thickenFeatures.createInput(
            faces, adsk.core.ValueInput.createByReal(t), False,
            adsk.fusion.FeatureOperations.NewBodyFeatureOperation, False)
        th = comp.features.thickenFeatures.add(tin)
        th.name = 'Outline Patch thicken'
        solid = None
        for b in th.bodies:
            if b.isSolid and b.entityToken not in before:
                solid = b
        if solid is None:
            raise RuntimeError('the thicken produced no new solid')
        solid.name = '%s %.2f mm' % (stem, float(o['thickness_mm']))
        rep.update({'body': solid.name, 'volume_mm3': round(solid.volume * 1000.0, 2),
                    'shells': solid.shells.count, 'solid_faces': solid.faces.count})
        if not o.get('keep_patch', True):
            patch.isLightBulbOn = False
        r_rim = float(o.get('rim_mm', 0.0))
        if r_rim > 0:
            t2 = time.time()
            edges, _counts = outer_perimeter(solid, patch, abs(t))
            ok, msg = _try_fillet(comp, edges, r_rim * MM, True)
            rep['rim'] = {'radius_mm': r_rim, 'edges': len(edges), 'ok': ok,
                          'msg': '' if ok else msg[:200]}
            if not ok and bend_before and bend_before < r_rim * 1.6:
                rep['rim']['hint'] = (
                    'the outline has convex corners down to R %.2f mm; a rim needs '
                    'clearly more room than its radius. Use a smaller rim, or ease the '
                    'corners (Advanced, about %.1f mm for a %.1f mm rim) - it moves the '
                    'outline locally and the report says by how much.'
                    % (bend_before, 1.6 * r_rim, r_rim))
            if ok:
                rep['volume_mm3'] = round(solid.volume * 1000.0, 2)
            tm['rim'] = round(time.time() - t2, 2)
    tm['thicken'] = round(time.time() - t1, 2)
    rep['timing_s'] = tm
    rep['loose_feet'] = getattr(surf, 'loose_feet', 0)
    rep['notes'] = notes
    rep['secs'] = round(time.time() - t0, 2)
    return rep


def format_report(r):
    lines = [
        'Outline: %.1f mm, %d samples every %.2f mm'
        % (r['outline_mm'], r['samples'], r['spacing_mm']),
        'Outline sits off the body by %.3f mm max, %.3f mm mean'
        % (r['outline_to_body_max_mm'], r['outline_to_body_mean_mm']),
        'Cut edge follows the projected outline to %.4f mm (%d edges)'
        % (r['cut_to_outline_max_mm'], r['cut_edges']),
        'Tightest convex corner of the outline: R %.2f mm%s' % (
            r.get('tightest_convex_bend_mm') or 0,
            ('; bends eased, outline moved %.3f mm max' % r['outline_eased_mm'])
            if r.get('outline_eased_mm') else ''),
        'Patch: %d faces, %.1f mm2' % (r['patch_faces'], r['patch_area_mm2']),
    ]
    if r.get('single_face_dev'):
        lines.append('Single face follows the body to %.4f mm max, %.4f mean'
                     % (r['single_face_dev']['max_mm'], r['single_face_dev']['mean_mm']))
    if 'body' in r:
        lines.append('Body "%s": %.1f mm3, %d shell(s)'
                     % (r['body'], r['volume_mm3'], r['shells']))
    if r.get('rim'):
        lines.append('Rim %.1f mm: %s' % (r['rim']['radius_mm'],
                     'rounded' if r['rim']['ok'] else 'NOT rounded'))
        if r['rim'].get('hint'):
            lines.append('  ' + r['rim']['hint'])
    lines.append('%.1f s' % r['secs'])
    return '\n'.join(lines)


# --------------------------------------------------------------------------- #
# the dialog
# --------------------------------------------------------------------------- #

def on_created(app, cmd):
    o = load_settings()
    inputs = cmd.commandInputs
    cmd.setDialogInitialSize(380, 420)
    cmd.okButtonText = 'Create'
    cmd.isExecutedWhenPreEmpted = False

    sel = inputs.addSelectionInput('outline', 'Outline',
                                   'The closed outline: sketch curves, 3D splines welcome')
    sel.addSelectionFilter('SketchCurves')
    sel.setSelectionLimits(1, 0)
    sel.tooltip = ('A closed loop of sketch curves drawn on or near the body. It does '
                   'not have to lie exactly on the surface: it only steers the cut.')
    bsel = inputs.addSelectionInput('body', 'Body', 'The body the patch lies on')
    bsel.addSelectionFilter('SolidBodies')
    bsel.addSelectionFilter('SurfaceBodies')
    bsel.setSelectionLimits(1, 1)

    dd = inputs.addDropDownCommandInput('keep', 'Keep',
                                        adsk.core.DropDownStyles.TextListDropDownStyle)
    for label in KEEPS:
        dd.listItems.add(label, label == o['keep'])
    dd = inputs.addDropDownCommandInput('direction', 'Thicken',
                                        adsk.core.DropDownStyles.TextListDropDownStyle)
    for label in DIRECTIONS:
        dd.listItems.add(label, label == o['direction'])
    inputs.addValueInput('thickness', 'Thickness', 'mm',
                         adsk.core.ValueInput.createByReal(o['thickness_mm'] * MM))
    inputs.addBoolValueInput('keep_patch', 'Keep the surface patch visible', True, '',
                             bool(o['keep_patch']))
    v = inputs.addValueInput('rim', 'Round the outer rim', 'mm',
                             adsk.core.ValueInput.createByReal(o['rim_mm'] * MM))
    v.tooltip = ('Fillet radius on the outer rim of the thickened body; 0 for none.\n\n'
                 'With a rim, the patch is rebuilt as ONE smooth face through points on '
                 'the body (within ~0.1 mm of it), because the many faces of a Form give '
                 'a rim no fillet builds on. A rim needs convex corners of the outline '
                 'clearly rounder than its radius.')
    v = inputs.addValueInput('min_bend', 'Ease corners tighter than', 'mm',
                             adsk.core.ValueInput.createByReal(o['min_bend_mm'] * MM))
    v.tooltip = ('0 leaves the outline exactly as drawn. Otherwise convex corners of '
                 'the outline tighter than this radius are relaxed, in the surface, '
                 'until they are not - only there. The report gives the largest move. '
                 'On the femur guide, 3.2 mm let a 2 mm rim build and moved the outline '
                 'by at most 0.5 mm at seven spots.')

    adv = inputs.addGroupCommandInput('advanced', 'Advanced')
    adv.isExpanded = False
    a = adv.children
    v = a.addValueInput('smooth', 'Normal smoothing', 'mm',
                        adsk.core.ValueInput.createByReal(o['smooth_mm'] * MM))
    v.tooltip = ('The surface normal is taken at every point of the outline and '
                 'averaged along it over about this length, so the cutting band '
                 'turns smoothly. 0 uses the raw normals.')
    v = a.addValueInput('curtain', 'Band half-height', 'mm',
                        adsk.core.ValueInput.createByReal(o['curtain_mm'] * MM))
    v.tooltip = ('How far the cutting band reaches either side of the surface. '
                 'Keep it below the smallest concave radius the outline crosses.')
    v = a.addValueInput('spacing', 'Outline sampling', 'mm',
                        adsk.core.ValueInput.createByReal(o['spacing_mm'] * MM))
    v.tooltip = 'How often the outline is projected onto the body.'
    v = a.addValueInput('rail', 'Band point spacing', 'mm',
                        adsk.core.ValueInput.createByReal(o['rail_mm'] * MM))
    v.tooltip = ('The cutting band runs through a projected point this often and is '
                 'a spline in between. Smaller follows tight turns more closely and '
                 'takes longer: 1 mm keeps the cut within ~0.02 mm of the projected '
                 'outline on a femur.')
    inputs.addTextBoxCommandInput('info', '', state().last_report or
                                  'Select the outline, then the body.', 6, True)


def read_options(inputs):
    return {
        'keep': inputs.itemById('keep').selectedItem.name,
        'direction': inputs.itemById('direction').selectedItem.name,
        'thickness_mm': inputs.itemById('thickness').value * CM,
        'keep_patch': inputs.itemById('keep_patch').value,
        'smooth_mm': inputs.itemById('smooth').value * CM,
        'curtain_mm': inputs.itemById('curtain').value * CM,
        'spacing_mm': inputs.itemById('spacing').value * CM,
        'rail_mm': inputs.itemById('rail').value * CM,
        'rim_mm': inputs.itemById('rim').value * CM,
        'min_bend_mm': inputs.itemById('min_bend').value * CM,
    }


def on_input_changed(app, inputs, changed_id):
    if changed_id == 'direction':
        solid = inputs.itemById('direction').selectedItem.name != DIR_NONE
        inputs.itemById('thickness').isEnabled = solid
        inputs.itemById('rim').isEnabled = solid


def on_validate(inputs):
    return (inputs.itemById('outline').selectionCount >= 1 and
            inputs.itemById('body').selectionCount == 1 and
            inputs.itemById('spacing').value > 0 and
            inputs.itemById('rail').value > 0 and
            inputs.itemById('rim').value >= 0 and
            inputs.itemById('min_bend').value >= 0 and
            inputs.itemById('curtain').value > 0)


def on_execute(app, ui, inputs):
    o = read_options(inputs)
    save_settings(o)
    sel = inputs.itemById('outline')
    curves = [adsk.fusion.SketchCurve.cast(sel.selection(i).entity)
              for i in range(sel.selectionCount)]
    body = adsk.fusion.BRepBody.cast(inputs.itemById('body').selection(0).entity)
    try:
        r = make_patch(app, curves, body, o)
    except ValueError as exc:
        ui.messageBox(str(exc), CMD_NAME)
        return
    state().last_report = format_report(r)
    log(json.dumps(r))
    if r.get('rim') and not r['rim']['ok']:
        ui.messageBox(format_report(r), CMD_NAME)


def on_destroy(app):
    pass


def clear_graphics(app):
    pass


# --------------------------------------------------------------------------- #
# buttons
# --------------------------------------------------------------------------- #

TOOLTIP = ('Cut the patch of a body that lies inside a closed outline, exactly, and '
           'thicken it into a new body.\n\n'
           'Draw the outline on or near the surface - a 3D spline is fine. At every '
           'point the local surface normal is found and smoothed along the outline, '
           'and a band along those normals cuts the surface, so the patch keeps the '
           'body\'s exact curvature and follows the outline\'s projection.')


def install(app, ui, created_handler_factory):
    old = ui.commandDefinitions.itemById(CMD_ID)
    if old:
        old.deleteMe()
    folder = RES_DIR if os.path.isdir(RES_DIR) else ''
    cd = ui.commandDefinitions.addButtonDefinition(CMD_ID, CMD_NAME, TOOLTIP, folder)
    cd.commandCreated.add(created_handler_factory())
    placed = []
    for pid in PANELS:
        panel = ui.allToolbarPanels.itemById(pid)
        if panel:
            old = panel.controls.itemById(CMD_ID)
            if old:
                old.deleteMe()
            panel.controls.addCommand(cd)
            placed.append(pid)
    return ', '.join(placed) or None


def uninstall(ui):
    for pid in PANELS:
        panel = ui.allToolbarPanels.itemById(pid)
        if panel:
            c = panel.controls.itemById(CMD_ID)
            if c:
                c.deleteMe()
    cd = ui.commandDefinitions.itemById(CMD_ID)
    if cd:
        cd.deleteMe()


# --------------------------------------------------------------------------- #
# dev bridge - how this gets tested inside a live Fusion
# --------------------------------------------------------------------------- #

def _design(app):
    return adsk.fusion.Design.cast(app.activeProduct)


def _body_by_name(app, name):
    for comp in _design(app).allComponents:
        for b in comp.bRepBodies:
            if b.name == name:
                return b
    raise ValueError('no body named %r' % name)


def _sketch_by_name(app, name):
    for comp in _design(app).allComponents:
        for sk in comp.sketches:
            if sk.name == name:
                return sk
    raise ValueError('no sketch named %r' % name)


def dev_ping(app, job):
    d = _design(app)
    return {'last_fit': LAST_FIT,'doc': app.activeDocument.name if app.activeDocument else None,
            'docs': [app.documents.item(i).name for i in range(app.documents.count)],
            'parametric': is_parametric(d) if d else None}


def dev_activate(app, job):
    for i in range(app.documents.count):
        doc = app.documents.item(i)
        if doc.name == job['name'] or doc.name.startswith(job['name']):
            doc.activate()
            return {'active': app.activeDocument.name}
    raise ValueError('no open document named %r' % job['name'])


def dev_run(app, job):
    """Run the pipeline on a named sketch (curves `indices`, default all) and body."""
    sk = _sketch_by_name(app, job['sketch'])
    idx = job.get('indices')
    curves = [sk.sketchCurves.item(i) for i in (idx if idx is not None
                                                else range(sk.sketchCurves.count))]
    o = dict(DEFAULTS)
    o.update(job.get('options') or {})
    r = make_patch(app, curves, _body_by_name(app, job['body']), o)
    r['text'] = format_report(r)
    return r


def dev_make_test(app, job):
    """A fresh document: a solid sphere-like blob and a wavy 3D outline hovering
    around it - drawn deliberately OFF the surface to prove the cut does not care."""
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    d = _design(app)
    if not is_parametric(d):
        d.designType = adsk.fusion.DesignTypes.ParametricDesignType
    root = d.rootComponent
    R = float(job.get('radius_mm', 20.0)) * MM
    tbm = adsk.fusion.TemporaryBRepManager.get()
    sph = tbm.createSphere(adsk.core.Point3D.create(0, 0, 0), R)
    base = root.features.baseFeatures.add()
    base.startEdit()
    root.bRepBodies.add(sph, base)
    base.finishEdit()
    body = base.bodies.item(0)
    body.name = 'TestBall'
    sk = root.sketches.add(root.xYConstructionPlane)
    sk.name = 'TestOutline'
    coll = adsk.core.ObjectCollection.create()
    n = 24
    wobble = float(job.get('wobble_mm', 0.4)) * MM
    for k in range(n):
        t = 2 * math.pi * k / n
        th = 0.55 + 0.12 * math.sin(3 * t)          # polar angle of the loop
        r = R + wobble * math.sin(5 * t)            # off the surface on purpose
        coll.add(adsk.core.Point3D.create(r * math.sin(th) * math.cos(t),
                                          r * math.sin(th) * math.sin(t),
                                          r * math.cos(th)))
    spl = sk.sketchCurves.sketchFittedSplines.add(coll)
    spl.isClosed = True
    return {'doc': doc.name, 'body': body.name, 'sketch': sk.name}


def dev_close(app, job):
    doc = app.activeDocument
    name = doc.name
    doc.close(False)
    return {'closed': name, 'now': app.activeDocument.name if app.activeDocument else None}


def dev_timeline(app, job):
    d = _design(app)
    tl = d.timeline
    names = []
    for i in range(max(0, tl.count - int(job.get('last', 12))), tl.count):
        try:
            names.append(tl.item(i).entity.name)
        except Exception:
            names.append('?')
    bodies = [{'name': b.name, 'solid': b.isSolid, 'visible': b.isVisible,
               'faces': b.faces.count} for b in d.rootComponent.bRepBodies
              if job.get('all') or 'Outline' in b.name or 'patch' in b.name
              or b.name in job.get('names', [])]
    return {'timeline_tail': names, 'count': tl.count, 'bodies': bodies}


def dev_measure(app, job):
    """Tessellate a body (mm) to a JSON file, for offline measurement."""
    b = _body_by_name(app, job['body'])
    calc = b.meshManager.createMeshCalculator()
    calc.surfaceTolerance = float(job.get('tolerance_mm', 0.01)) * MM
    m = calc.calculate()
    with open(job['path'], 'w', encoding='utf-8') as fh:
        json.dump({'coords_mm': [v * CM for v in m.nodeCoordinatesAsDouble],
                   'idx': list(m.nodeIndices)}, fh)
    return {'path': job['path'], 'tris': m.triangleCount}


def dev_remove_last(app, job):
    """Delete the most recent Outline Patch output (thicken + its base feature)."""
    d = _design(app)
    tl = d.timeline
    removed = []
    for want in ('Outline Patch thicken', 'Outline Patch'):
        for i in range(tl.count - 1, -1, -1):
            ent = tl.item(i).entity
            if ent is not None and getattr(ent, 'name', None) == want:
                bodies = [b.name for b in getattr(ent, 'bodies', [])]
                ent.deleteMe()
                removed.append({'feature': want, 'bodies': bodies})
                break
    return {'removed': removed, 'timeline_count': tl.count}


def dev_make_from_files(app, job):
    """A fresh document with a Form from a TSM file and a closed fitted spline
    through the points of a JSON file (mm) - the real case, on a clean body."""
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    d = _design(app)
    if not is_parametric(d):
        d.designType = adsk.fusion.DesignTypes.ParametricDesignType
    root = d.rootComponent
    ff = root.features.formFeatures.add()
    ff.startEdit()
    ff.tSplineBodies.addByTSMDescription(open(job['tsm'], encoding='utf-8').read())
    ff.finishEdit()
    body = ff.bodies.item(0)
    body.name = 'Form'
    pts = json.load(open(job['points'], encoding='utf-8'))['points_mm']
    step = max(1, int(job.get('every', 4)))
    sk = root.sketches.add(root.xYConstructionPlane)
    sk.name = 'Outline'
    coll = adsk.core.ObjectCollection.create()
    for q in pts[:-1][::step]:
        coll.add(adsk.core.Point3D.create(q[0] * MM, q[1] * MM, q[2] * MM))
    spl = sk.sketchCurves.sketchFittedSplines.add(coll)
    spl.isClosed = True
    return {'doc': doc.name, 'body': body.name, 'faces': body.faces.count,
            'fit_points': coll.count}


def outer_perimeter(solid, patch, thickness_cm):
    """The edges between the offset (top) faces of a thickened patch and its side
    walls. Faces are told apart by where they sit relative to the patch: the top
    ones a thickness away with the same normal, the walls across it."""
    comp = solid.parentComponent
    face_pts = [_tup(f.pointOnFace) for f in solid.faces]
    surf = Surface(comp, patch, face_pts, 0.05)
    kind = {}
    for f, p in zip(solid.faces, face_pts):
        ev = f.evaluator
        ok, uv = ev.getParameterAtPoint(_pt(p))
        ok, n = ev.getNormalAtParameter(uv)
        n = _unit(_tup(n))
        sgn = -1.0 if f.isParamReversed else 1.0
        probe = _add(p, n, 0.01 * sgn)
        if solid.pointContainment(_pt(probe)) == \
                adsk.fusion.PointContainment.PointInsidePointContainment:
            sgn = -sgn
        n = (n[0] * sgn, n[1] * sgn, n[2] * sgn)
        d, _q, npatch = surf.closest(p)
        dot = _dot(n, npatch)
        if abs(dot) < 0.5:
            kind[f.tempId] = 'side'
        elif d > 0.5 * thickness_cm:
            kind[f.tempId] = 'top'
        else:
            kind[f.tempId] = 'bottom'
    edges = []
    for e in solid.edges:
        ks = sorted(kind.get(f.tempId) for f in e.faces)
        if ks == ['side', 'top']:
            edges.append(e)
    counts = {}
    for k in kind.values():
        counts[k] = counts.get(k, 0) + 1
    return edges, counts


def _order_chain(edges):
    """Edges of a closed loop in walking order."""
    left = list(edges)
    out = [left.pop(0)]
    v = out[0].endVertex
    while left:
        nxt = None
        for i, e in enumerate(left):
            if e.startVertex == v or e.endVertex == v:
                nxt = i
                break
        if nxt is None:
            out.extend(left)
            break
        e = left.pop(nxt)
        v = e.endVertex if e.startVertex == v else e.startVertex
        out.append(e)
    return out


def _try_fillet(comp, edges, r_cm, keep, chain=True):
    coll = adsk.core.ObjectCollection.create()
    for e in edges:
        coll.add(e)
    fin = comp.features.filletFeatures.createInput()
    fin.edgeSetInputs.addConstantRadiusEdgeSet(
        coll, adsk.core.ValueInput.createByReal(r_cm), chain)
    try:
        f = comp.features.filletFeatures.add(fin)
    except Exception:
        return False, short_error(traceback.format_exc())
    if not keep:
        f.deleteMe()
    return True, f.name if keep else ''


def dev_fillet(app, job):
    """Fillet the outer perimeter of a thickened patch; on failure, find where."""
    solid = _body_by_name(app, job['body'])
    patch = _body_by_name(app, job['patch'])
    comp = solid.parentComponent
    edges, counts = outer_perimeter(solid, patch, float(job.get('thickness_mm', 3.7)) * MM)
    lens = sorted(e.length * CM for e in edges)
    out = {'faces': counts, 'perimeter_edges': len(edges),
           'perimeter_mm': round(sum(lens), 2),
           'shortest_mm': [round(x, 4) for x in lens[:5]]}
    r = float(job.get('radius_mm', 2.0)) * MM
    for rr in job.get('also_mm', []):
        okr, _m = _try_fillet(comp, edges, rr * MM, False)
        out['whole_%gmm' % rr] = okr
    ok, msg = _try_fillet(comp, edges, r, bool(job.get('keep')))
    out['whole'] = {'ok': ok, 'msg': msg}
    if not ok and job.get('locate', True):
        # walk the chain in pieces to see where a fillet of this radius breaks
        chain = _order_chain(edges)
        k = int(job.get('pieces', 12))
        per = max(1, len(chain) // k)
        pieces = []
        for i in range(0, len(chain), per):
            part = chain[i:i + per]
            okp, msgp = _try_fillet(comp, part, r, False, chain=False)
            mid = _edge_points(part[len(part) // 2], 1)[0]
            pieces.append({'edges': len(part), 'ok': okp,
                           'at_mm': [round(v * CM, 1) for v in mid],
                           'msg': '' if okp else msgp[:160]})
        out['pieces'] = pieces
    return out


DEV_COMMANDS = {
    'fillet': dev_fillet,
    'make_from_files': dev_make_from_files,
    'remove_last': dev_remove_last,
    'ping': dev_ping,
    'activate': dev_activate,
    'run': dev_run,
    'make_test': dev_make_test,
    'close': dev_close,
    'timeline': dev_timeline,
    'measure': dev_measure,
}
