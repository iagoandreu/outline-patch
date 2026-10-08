# Outline Patch

A Fusion 360 add-in that takes a closed outline drawn on or near a curved body and
returns the part of the body's surface inside it, with the body's exact curvature,
thickened into a new body. It is 3-matic's *attach curve + uniform offset*, for
surgical guides and baseplates drawn on a bone Form.

**Surface › Create › Outline Patch** (also on Solid › Create).

| Input | |
|---|---|
| Outline | a closed loop of sketch curves; a 3D spline is fine, and it does not have to lie on the surface |
| Body | a solid or surface body (a scan has to become a body first, e.g. with Mesh to Part) |
| Keep | the smaller or the larger side of the outline |
| Thicken | outward along the surface normal, inward, or no thickness (surface only) |
| Thickness | 3.7 mm by default |

The result is two timeline entries: **Outline Patch**, a base feature holding the surface
patch, and **Outline Patch thicken**. Nothing else is left behind.

## Why not just split the face with the curve

Split Face with a curve *projects* the curve onto the faces. A hand-drawn 3D spline
never lies exactly on the surface, and where it does not, the imprinted edge drifts or
jumps between faces. Here the curve never cuts anything. It only steers a band:

1. The outline is sampled every 0.5 mm and every sample is projected to its exact
   closest point on the body, with the outward normal there.
2. The normals are smoothed along the outline (2 mm by default), so the band turns smoothly.
3. One spline is fitted through the projected points, and its control points are
   offset by ±1.5 mm along the smoothed normals. The ruled surface between the two
   offsets is the band. Because both offsets share one set of knots, the band's
   mid-line *is* that spline.
4. The faces near the outline are split with the band (an exact surface/surface
   intersection), and the side inside the loop is lifted out and thickened.

So the cut lands at the closest-point projection of the outline, whatever the
outline's distance from the body.

## Measured

On a 499 mm artificial femur (Form of 6044 quads, 415 faces) with a 354.5 mm
outline drawn on the original STL (0.46 mm max, 0.05 mm mean off the Form):

| | |
|---|---|
| patch boundary vs projected outline | 0.021 mm max, 0.004 mm mean |
| patch | one disk, 2861 mm² |
| baseplate | one closed shell, 11 346 mm³ |
| outer face vs the original scan | 3.70 mm median, 3.58–3.78 p1–p99 |
| time | 16 s on a clean Form |

On a sphere of radius 20 mm with an outline 0.4 mm off the surface the projected points
match the exact answer to 1e-6 mm, and a 3.7 mm outward shell on a 378.2 mm² cap comes
out at 1673.8 mm³ against the analytic 1673.6 mm³.

## Traps found on the way

- `MeasureManager.measureMinimumDistance` with a temporary `Point3D` returned distance
  0 and the query point itself, for a point 0.4 mm off a sphere. It is not used.
- `SurfaceEvaluator.getParameterAtPoint` is not a projection. From 0.5 mm off a Form's
  NURBS face it converged 94 mm away. It is only ever asked about points already on the
  surface (the projection is done on a tessellation, then refined by Gauss-Newton in
  (u, v)).
- Two separately fitted rail splines drift apart in parameter, because offsetting along
  the normal changes the length. The ruled band then crossed the surface 0.2 mm off the
  outline. Offsetting the control points of a single spline fixes it.
- An old imprint of a similar curve on the body leaves edges within 0.1 mm of the cut
  and thin faces along it. Every face within 0.2 mm of the outline is split, coincident
  old edges count as cut, and slivers under 1 % of the area are ignored when the side
  is chosen. It also makes the split slow (38 s instead of 16 s).

## Install

```
powershell -ExecutionPolicy Bypass -File install.ps1
```

then in Fusion: Utilities › Add-Ins › Scripts and Add-Ins › Add-Ins › Outline Patch ›
Run (tick *Run on Startup*). `op_core.py` reloads on every click; only
`OutlinePatch.py` needs a Stop/Run.

## Tests

```
py -3.10 tests/test_geometry.py
```

Offline, with `adsk` stubbed: resampling, normal smoothing, the closest point on a
triangle and the outline distance grid. The pipeline itself is tested live through the
dev bridge (`devbridge.on` beside the add-in), with `make_test` (a sphere and an outline
0.4 mm off it) and `make_from_files` (a Form from a TSM file and an outline from a JSON
file).
