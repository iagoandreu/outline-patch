"""Offline tests for the parts of Outline Patch that need no Fusion.

    py -3.10 tests/test_geometry.py

adsk is stubbed (see tests/stub), so op_core imports outside Fusion. Every check
has an answer known in advance.
"""

import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'stub'))
sys.path.insert(0, os.path.dirname(HERE))

import op_core as O                                          # noqa: E402

FAILED = []


def check(name, cond, detail=''):
    print('  %s %s%s' % ('ok  ' if cond else 'FAIL', name,
                         '' if cond else '   %s' % (detail,)))
    if not cond:
        FAILED.append(name)


def circle(n, r=1.0, z=0.0):
    return [(r * math.cos(2 * math.pi * k / n), r * math.sin(2 * math.pi * k / n), z)
            for k in range(n)]


print('resampling a closed polyline')
pts, total, step = O.resample_closed(circle(400, 2.0), 0.05, 10000)
check('a circle of radius 2 measures 2 pi', abs(total - 4 * math.pi) < 1e-3, total)
check('samples are evenly spaced', abs(step - total / len(pts)) < 1e-12)
gaps = [O._norm(O._sub(pts[(i + 1) % len(pts)], pts[i])) for i in range(len(pts))]
check('consecutive gaps equal the step', max(abs(g - step) for g in gaps) < 2e-4,
      max(abs(g - step) for g in gaps))
pts2, _t, _s = O.resample_closed(circle(400, 2.0), 0.0001, 300)
check('the sample cap holds', len(pts2) == 300, len(pts2))

print('smoothing normals along the outline')
random.seed(1)
base = [(0.0, 0.0, 1.0)] * 200
noisy = [O._unit((x + random.uniform(-0.2, 0.2), y + random.uniform(-0.2, 0.2), z))
         for x, y, z in base]
sm = O.smooth_normals(noisy, 0.05, 0.2)
ang = lambda ns: [math.degrees(math.acos(min(1.0, O._dot(n, (0, 0, 1))))) for n in ns]
rms = lambda v: math.sqrt(sum(x * x for x in v) / len(v))
err_before, err_after = rms(ang(noisy)), rms(ang(sm))
# a Gaussian of sigma 0.2 over samples 0.05 apart averages ~14 of them: the rms
# of independent noise falls by about sqrt(14) = 3.7
check('smoothing pulls noisy normals towards the true one', err_after < err_before / 3.0,
      (err_before, err_after))
check('smoothed normals stay unit length', max(abs(O._norm(n) - 1) for n in sm) < 1e-12)
ring = [(math.cos(t), math.sin(t), 0.0) for t in (2 * math.pi * k / 360 for k in range(360))]
sm = O.smooth_normals(ring, 1.0, 0.0)
check('zero smoothing returns the normals untouched', sm == ring)
sm = O.smooth_normals(ring, 2 * math.pi / 360, 0.05)
# a uniformly turning field is symmetric about every sample: only the length shrinks
turn = max(math.degrees(math.acos(min(1.0, O._dot(a, b)))) for a, b in zip(ring, sm))
check('a uniformly turning field keeps its directions', turn < 1e-4, turn)

print('closest point on a triangle')
a, b, c = (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)
cases = [((0.2, 0.2, 0.7), (0.2, 0.2, 0.0)),        # face interior
         ((-1.0, -1.0, 0.3), (0.0, 0.0, 0.0)),      # vertex a
         ((2.0, -0.5, 0.0), (1.0, 0.0, 0.0)),       # vertex b
         ((0.5, -1.0, 0.0), (0.5, 0.0, 0.0)),       # edge ab
         ((-1.0, 0.5, 1.0), (0.0, 0.5, 0.0)),       # edge ac
         ((1.0, 1.0, 0.0), (0.5, 0.5, 0.0))]        # edge bc
for p, want in cases:
    got = O._closest_on_triangle(p, a, b, c)
    check('closest to %s is %s' % (p, want), O._norm(O._sub(got, want)) < 1e-12, got)
random.seed(2)
worst = 0.0
for _ in range(2000):
    p = tuple(random.uniform(-2, 2) for _ in range(3))
    got = O._closest_on_triangle(p, a, b, c)
    # brute force over a fine barycentric grid
    best = min(O._norm(O._sub(p, (u / 60.0, v / 60.0, 0.0)))
               for u in range(61) for v in range(61 - u))
    worst = max(worst, O._norm(O._sub(p, got)) - best)
check('never farther than a 1/60 barycentric grid', worst < 1e-12, worst)

print('point-to-outline distance through the grid')
loop = circle(720, 3.0)
g = O._Grid(loop, 0.2)
d_on = max(g.dist(q) for q in circle(97, 3.0))
check('points on the circle are on the polyline (to its chord error)', d_on < 3.0 * (1 - math.cos(math.pi / 720)) + 1e-9, d_on)
check('a point 0.05 off reads 0.05', abs(g.dist((3.05, 0.0, 0.0)) - 0.05) < 1e-9)
check('a point beyond the grid reach reads infinity', g.dist((10.0, 0.0, 0.0)) == float('inf'))

print('')
if FAILED:
    print('%d FAILED: %s' % (len(FAILED), ', '.join(FAILED)))
    sys.exit(1)
print('all geometry tests passed')
