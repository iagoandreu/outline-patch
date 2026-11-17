# Outline Patch

A Fusion 360 add-in that takes a closed outline drawn on or near a curved body and
returns the part of the body's surface inside it, with the body's exact curvature,
thickened into a new body (an attached curve plus a uniform offset), for
surgical guides and baseplates drawn on a bone Form.

## How to use

### Where the button is

Design workspace. The **Outline Patch** button is placed in two panels:

- **Surface** tab > **Create** panel
- **Solid** tab > **Create** panel

### Before you start

- **A body**: a solid or surface body, for example the bone as a Form or a solid.
  Mesh bodies cannot be picked; convert a scan first (for example with Mesh to
  Part).
- The body must be picked in its own component. If it sits in a sub-component,
  activate that component first; otherwise the tool stops with "pick the body in
  its own component".
- **An outline**: a closed loop of sketch curves drawn on or near the body. A
  single closed 3D spline is fine, and so are several curves that join end to end
  (gaps up to 0.05 mm). It does not have to lie exactly on the surface: every
  point is projected to the closest point on the body. Construction curves are
  ignored.
- **External Python - only for the rounded rim.** When **Round the outer rim** is
  above 0, the body is thickened, and Keep is the smaller side, the patch is
  rebuilt as one smooth face. That step runs in a separate Python 3 with **numpy**
  and **scipy** (`pip install numpy scipy`), because Fusion's Python has neither.
  The add-in uses the first interpreter that can `import numpy, scipy`, in this
  order:
  1. the path in the environment variable **`OUTLINEPATCH_PYTHON`** (set before
     Fusion starts);
  2. the path written in **`python_path.txt`** beside the installed add-in
     (`%APPDATA%\Autodesk\Autodesk Fusion 360\API\AddIns\OutlinePatch\python_path.txt`),
     one line, the full path to `python.exe`;
  3. `%LOCALAPPDATA%\Programs\Python\Python310\python.exe`;
  4. `py -3.10`.

  If none works, Create stops with "the single-face option needs a Python with
  numpy and scipy". Without a rim (Round the outer rim = 0), no external Python is
  needed.

### Step by step

1. Click **Outline Patch**.
2. **Outline**: select the sketch curves of the closed outline.
3. **Body**: select the body the patch lies on.
4. Set **Keep**, **Thicken**, **Thickness** and the rim (table below).
5. Click **Create**. A run takes seconds without a rim and can take a couple of
   minutes with one.

| Input | Default | What it does |
|---|---|---|
| **Outline** | none | The closed outline: one or more sketch curves. |
| **Body** | none | One solid or surface body. |
| **Keep** | Smaller side of the outline | **Smaller side of the outline** or **Larger side of the outline**: which of the two regions the outline splits the surface into becomes the patch. |
| **Thicken** | Outward (along the surface normal) | **Outward (along the surface normal)**, **Inward**, or **No thickness, surface patch only**. With the last one, Thickness and the rim are disabled. |
| **Thickness** | 3.7 mm | Thickness of the new body. |
| **Keep the surface patch visible** | on | Off: the surface patch is hidden after the thicken (it stays in the design). |
| **Round the outer rim** | 2 mm | Fillet radius on the outer rim of the thickened body. 0 = no fillet. Above 0 the patch is rebuilt as one smooth face first (needs the external Python). |
| **Ease corners tighter than** | 0 mm | 0 leaves the outline exactly as drawn. Above 0, convex corners of the outline tighter than this radius are relaxed in the surface until they are not, and only there. The report gives the largest move. |

**Advanced** group (collapsed):

| Input | Default | What it does |
|---|---|---|
| **Normal smoothing** | 2 mm | The surface normal at each outline point is averaged along the outline over about this length, so the cutting band turns smoothly. 0 uses the raw normals. |
| **Band half-height** | 1.5 mm | How far the cutting band reaches either side of the surface. Keep it below the smallest concave radius the outline crosses. |
| **Outline sampling** | 0.5 mm | How often the outline is projected onto the body. |
| **Band point spacing** | 1 mm | The band passes through a projected point this often and is a spline in between. Smaller follows tight turns more closely and takes longer. |

All values are remembered for the next time (`last_settings.json` in the add-in
folder).

### What you get

In the component that owns the body:

- a surface body named `<outline sketch name> patch` - the part of the body's
  surface inside the outline. In a parametric design it sits in a base feature
  named **Outline Patch**;
- unless Thicken is "No thickness": a thicken feature **Outline Patch thicken** and
  a new solid body named `<outline sketch name> <thickness> mm`, for example
  `Sketch3 3.70 mm`;
- with a rim that built: a fillet feature **Outline Patch rim**.

The original body and the outline are not changed. The working copies, bands and
helper sketches are removed at the end.

### Checking the result

The report appears in the dialog's text box the next time you open **Outline
Patch** (and as a message box straight away if the rim could not be rounded).
It is also written to `OutlinePatch.log` in the add-in folder. It gives:

| Line | What to look at |
|---|---|
| Outline: length, samples, spacing | - |
| Outline sits off the body by ... max, ... mean | How far your drawn outline was from the surface. Large values mean the projection moved it a lot. |
| Cut edge follows the projected outline to ... mm | Should be a few hundredths of a mm. |
| Tightest convex corner of the outline: R ... mm | A rim needs clearly more room than its radius. If you eased corners, how far the outline moved. |
| Patch: faces, area in mm2 | One face when the smooth single face was used. |
| Single face follows the body to ... max, ... mean | Only with a rim: how far the rebuilt face is from the original body. |
| Body "...": volume in mm3, shells | Should be 1 shell. |
| Rim ... mm: rounded / NOT rounded | If not rounded, a hint with a smaller rim or an ease value follows. |

Before you trust the part, also check:

- the cut edge in the viewport against your outline, especially at tight turns;
- the thickening direction (outward away from the bone, or inward);
- the single-face deviation when a rim was used, since that face is a fitted
  surface, not the original one.

### Where it fits

Use it once the bone is a Form or solid body and the guide or plate outline is
drawn on it (by hand, or from Intersection Sketch or Bone Curve). It gives the
bone-contacting base of a guide or plate with the bone's exact curvature, ready for
drill tubes, screw holes and export for printing.

### Known limits

- The outline must be closed. An open outline, or curves with gaps over 0.05 mm,
  are refused.
- The outline is sampled at no more than 900 points; on a very long outline the
  real spacing becomes larger than **Outline sampling**.
- A rim needs room at the outline's convex corners: the fillet fails when a corner
  radius is close to the rim radius. Use a smaller rim or **Ease corners tighter
  than** (it moves the outline locally).
- The smooth single face, and so a reliable rim on a Form, is only built when Keep
  is the smaller side.
- The single face is fitted through points projected onto the body. It is close to
  the body (about 0.1 mm on a femur test), not identical.
- With a rim, a run on a large bone takes about two minutes; without, seconds.
- An old imprint of a similar curve on the body makes the split slower.

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

## A rounded rim

**Round the outer rim** (2 mm by default) fillets the edge between the top face and the
side wall. On the Form's own faces no fillet builds at any radius (0.3 mm included):
43 patch faces thicken into 43 top faces, 85 side walls and an 85-edge rim with edges of
0.05 mm, and Fusion's Merge Faces refuses to join T-spline faces. So with a rim the patch
is rebuilt as **one smooth face**: the region is flattened conformally (LSCM, in
`op_fit.py`, run by an external Python with numpy and scipy), sections are
laid every 4 mm across one rectangle of that domain, every point over the region is
projected exactly onto the body, and one face is lofted through them. It is cut with the
band's centre spline, which lies on it, so the patch has a single closed edge. On the femur
that face is within 0.097 mm max, 0.010 mm mean of the Form (6 mm sections: 0.29 mm, and
the rim failed).

A rim also needs room at the outline's **convex corners**: the fillet's spine runs a
radius inside the rim, R - r round a corner of radius R. The GUIDE outline has convex
corners down to R 1.66 mm; rims up to 1.8 mm built, 1.9 and 2.0 mm did not (2.1 mm did,
a degenerate coincidence). **Ease corners tighter than** relaxes only those corners, in
the surface, and reports the largest move: 3.2 mm let the 2 mm rim build and moved the
outline by at most 0.50 mm, at seven spots 2.5-5 mm long; the other 650 samples did not
move. It is off by default, because it changes the outline.

With a rim the run takes about two minutes on the femur (Thicken of the lofted face
alone ~90 s); without one, ~16 s.

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

The rim's single face needs a Python with numpy and scipy; the add-in looks for
`OUTLINEPATCH_PYTHON`, then `python_path.txt` beside it, then
`%LOCALAPPDATA%\Programs\Python\Python310\python.exe`, then `py -3.10`
(see *Before you start*).

## Tests

```
py -3.10 tests/test_geometry.py
```

Offline, with `adsk` stubbed: resampling, normal smoothing, the closest point on a
triangle and the outline distance grid. The pipeline itself is tested live through the
dev bridge (`devbridge.on` beside the add-in), with `make_test` (a sphere and an outline
0.4 mm off it) and `make_from_files` (a Form from a TSM file and an outline from a JSON
file).
