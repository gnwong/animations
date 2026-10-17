GRMHD MAGNETIC-FIELD ANIMATION STARTER

Purpose
  Reproduce the accepted square scientific illustration from an AthenaK GRMHD
  .bin snapshot, then carry the same recipe to a sequence of snapshots.
  Start by reading QUALITY_CRITERIA.txt and HANDOFF.txt.

What is here
  settings.json       Complete accepted square recipe; no input paths.
  pipeline.py         CLI for sample / prepare / render / run / batch / doctor.
  prepare.py          Bin resampling, display-density surfaces, native B traces.
  renderer.py         Lighting, native-field tubes, full 3D occlusion, bloom.
  common.py           Config validation, binary header inspection, provenance.
  keyframes.py        Time-based keyframe selection and settings interpolation.
  tests.py            Binary-layout, interpolation, and stable-ID tests.
  run.sh              Local runtime/cache paths and chosen Python interpreter.
  requirements.txt    Versions tested here; wongutils comes from your environment.
  reference/          Accepted square PNG, source/reader provenance, validation.

No snapshot, resampled cache, large VTK, VisIt session, concept art, or discarded
experiment is needed in this directory. Copy this directory as one unit.
wongutils must already be importable in the chosen Python environment.

Environment: Linux, Python 3.11, OpenGL rendering
  Activate an environment that imports wongutils.grmhd.athenak.AthenaKSnapshot.
  If necessary, install the other dependencies with:
    python -m pip install -r requirements.txt
  These pins reproduce the locally tested environment, not a requirement to
  replace a working remote environment. Verify another environment as below.

  Inspect imports/versions:
    bash run.sh doctor
  Verify an actual translucent OpenGL render and depth peeling:
    bash run.sh doctor --render
  Run CPU-only workflow tests:
    TMPDIR="$PWD/runtime/tmp" PYTHONDONTWRITEBYTECODE=1 python tests.py
  (Run doctor first so runtime/tmp exists, or use the system TMPDIR.)

Headless Linux
  The tested VTK 9.3 build uses X11. SetOffScreenRendering alone does not make
  that build independent of an X server. If a GPU X server is available, use
  its DISPLAY. Otherwise install Xvfb, xauth, and Mesa/GLX through the host's
  package manager, then run, for example:
    xvfb-run -a -s '-screen 0 1920x1080x24' bash run.sh doctor --render
    xvfb-run -a -s '-screen 0 1920x1080x24' bash run.sh run /data/snapshot.bin --work work/example --out frames/example.png
  Software GL may be slower. Other EGL/OSMesa/VTK builds must pass doctor and
  a reference comparison before using them for production. The renderer fails
  if depth peeling is unavailable, rather than silently changing the depth.
  VisIt is not required: the accepted image uses standalone VTK rendering.

Render one snapshot
  From this directory, with your Python environment activated:
    bash run.sh run /data/torus.mhd_w_bcc.01000.bin --work work/01000 --out frames/01000.png

  Or run the stages separately:
    bash run.sh sample /data/torus.mhd_w_bcc.01000.bin --work work/01000
    bash run.sh prepare --work work/01000
    bash run.sh render --work work/01000 --out frames/01000.png

  Preview at smaller IMAGE resolution without changing the physical mesh:
    bash run.sh render --work work/01000 --out previews/01000.png --width 512 --height 512
  Add --passes to render/run/batch to save disk, raw beauty, and horizon-mask
  PNGs alongside the final image. The default saves the final PNG + JSON only.

  --settings PATH accepts a partial JSON override of settings.json. Unknown
  settings are rejected. Width/height CLI overrides affect the image only.
  Raw density/B are not modified. Spin and simulation time come from the bin
  header, including the horizon geometry; they are not fixed to the reference.

Input contract and coordinates
  Full 3D AthenaK binary output version 1.1, read by installed wongutils.
  Required native variables: dens, bcc1, bcc2, bcc3. Coordinates/components are
  Cartesian Kerr-Schild in code length units rg. A different GRMHD format or
  coordinate basis needs an explicit reader/basis adapter; do not merely
  rename variables. Required header groups include mesh, meshblock, coord/a.
  The sampling box must lie inside the simulation domain. Ghost zones are
  populated before trilinear interpolation. Samples are float32 (z,y,x,4).
  Use completed snapshot files; this starter does not watch for ongoing writes.

Accepted processing recipe
  401^3 sinh-stretched grid: x/y +/-150 rg, z +/-240 rg, stretch=4.
  Density display gain .3; Gaussian widths .65/2.2/3.4 in grid-index units.
  Middle blend R6-30, stronger outer blend R12-40; five density contours.
  Windowed-sinc smoothing, radially blended, then decimation and smooth normals.
  B is NOT Gaussian-smoothed. RK45 uses every second sampled grid node.
  64 seeds at z=+/-2.8, radii .35/.7/1.15/1.7, eight azimuths each; the
  pattern rotates with simulation time while its topology/identities stay fixed.
  128 possible seed/direction branches; 96 chosen by fixed permutation.
  The baseline palette is gold inside, warmer gray/silver outside, blue-white
  field tubes. Outer surface opacity is reduced at grazing angles and large R.
  Camera inclination=52 deg, azimuth=200 deg, roll=-10 deg, scale=82.
  A final 180-degree screen rotation makes the approaching negative-z jet go up.
  Full 3D geometry occlusion precedes conventional bloom and sharpening.
  Output defaults to 1600 x 1600 with 2x supersampling and area downsampling.
  Bloom/sharpening follow projected physical scale when size/aspect/zoom change.
  These are illustrative colors/opacity, not physical temperature or GRRT
  emissivity. This is a coordinate-space simulation illustration, not a lensed
  observer image. No AI generation or simulated polarization vectors are added.

Moving magnetic-field footpoints (enabled by default)
  The original splash_rendering script used Omega_F=Omega_H/2 and phase=Omega_F*t.
  The portable version carries that rotation law onto the accepted 64-seed
  ring layout, rather than restoring the old spherical seed cloud. In G=M=c=1
  code units, with spin a and r_H=1+sqrt(1-a^2):
    Omega_H = a / (2*r_H)
    Omega_F = field.seed_omega_fraction * Omega_H
    phi(t) = phi_initial + field.seed_phase
             + Omega_F * (snapshot_header_time - field.seed_reference_time)
  seed_motion="rotate", seed_omega_fraction=.5, seed_reference_time=10450,
  seed_phase=0 radians are the defaults. Thus the accepted snapshot at t=10450
  retains its exact original seed positions/image. Positive/negative spin
  determines rotation direction; spin zero gives no rotation. Use an appropriate
  reference_time/phase for a new run. Set reference_time=0 to use the original
  unanchored Omega_F*t convention, or seed_motion="fixed" to disable rotation.
  Inspect phase/frequency without loading the meshblock payload:
    bash run.sh seed-info /data/snapshot.bin
    bash run.sh seed-info /data/snapshot.bin --points
  Actual header time is used, never filename index or movie-frame number.
  New native-field integral curves are traced from the rotated starting points
  in each snapshot. Neither rendered curves nor images are rotated to imitate
  this effect. Seed IDs and the chosen 96 directions remain attached to the
  same pattern positions. Prepared/final metadata record phase, frequency,
  time, spin, and all actual seed coordinates; cache keys include the motion.
  This is a corotating seeding prescription, not fluid-velocity advection of
  Lagrangian plasma footpoints or measured field-line angular velocity.
  The rotation policy is global/analytic, not a per-keyframe interpolated phase.
  Sparse snapshot cadence can alias the rotation; review neighboring frames.

Artifacts, caching, and resources
  Each work directory holds samples.npy/json, disk.vtk, fields.vtk, prepared.json.
  Work directories belong to individual snapshots. Never run two workers on
  the same work directory. One worker is used by default.
  Cache signatures include source path/stat/header fingerprint, sampling,
  reader module hashes, geometry/trace settings, and VTK version. Changing a
  rendering-only setting reuses geometry; changing density_gain rebuilds it.
  Source payloads are not cryptographically hashed in full on every run.
  Completed batch renders are reused only if settings, source fingerprint,
  versions, and final PNG hash still match.
  After editing sampler/geometry/rendering implementation code, use fresh work
  directories and output paths. Cache signatures cover settings/dependencies,
  not arbitrary edits to this package's source. Reference provenance records
  the shipped source hashes for comparison in a new context.
  The default sample cache alone is about 1 GB per snapshot. Native loading,
  ghost zones, smoothing, surfaces, and rendering need additional RAM/VRAM.
  Allow roughly 32 GB RAM per worker initially and measure your dataset before
  increasing concurrency. Keep caches for the current review set, not blindly
  for an entire long sequence. Source files are always opened read only.

100-keyframe workflow
  Use snapshots from ONE simulation, not a mixture of runs/spins. The original
  00800 and 01000 comparison files were different spins/runs; do not treat them
  as successive animation frames.

  1. Select approximately uniformly spaced simulation times:
    python keyframes.py plan --glob '/data/run/*.bin' --count 100 --out keyframes.json
     This reads headers only, preserves endpoints, and rejects mixed spins or
     duplicate times. Fewer than 100 inputs requires a smaller --count.

  2. Materialize settings for those snapshots and render review images:
    python keyframes.py schedule --plan keyframes.json --keyframes-only --out keyframe_schedule.json
    bash run.sh batch --schedule keyframe_schedule.json --work-root work/keyframes --out-dir previews/keyframes --width 512 --height 512

  3. Review all keyframes against QUALITY_CRITERIA.txt. Prefer one global
     setting correction over many local patches. Edit per-keyframe overrides
     and notes in keyframes.json; mark reviewed=true after inspection. Render
     difficult keyframes at full resolution before accepting them.
     Example overrides:
       "overrides": {"render.scale": 86, "render.warm_glow": 0.55}
     Each anchor starts from base_settings. An omitted override means the BASE
     value at that anchor, not the previous anchor's value. Edit the base file
     for persistent global changes. Regenerate the schedule after edits.

  4. Interpolate continuous settings at every available raw snapshot time:
    python keyframes.py schedule --plan keyframes.json --glob '/data/run/*.bin' --require-reviewed --out schedule.json
    bash run.sh batch --schedule schedule.json --work-root work/sequence --out-dir frames

     PCHIP interpolation preserves scalar shapes; azimuth/roll are unwrapped
     along the shortest angular path first. No extrapolation is performed.
     Every resulting configuration is validated, including coupled radius
     constraints. Structural changes (grid/seed layout, rotation policy,
     contours, trace count)
     cannot be per-keyframe overrides. They require a new global recipe/review.
     Footpoint phase is evaluated analytically from each snapshot's header time,
     independently of how many keyframes are selected. Regenerate schedules
     created with the earlier fixed-footpoint version of this package.

  5. Review adjacent frames and interval midpoints. Additional keyframes near
     transients may be needed even if the initial 100 look good. The code
     interpolates RENDER SETTINGS, not GRMHD states, field trajectories, or
     bitmap images. There must be actual simulation outputs between keyframes.

  Batch supports --start / --stop for disjoint job ranges; numbering is taken
  from the schedule. It resumes matching completed outputs. For an optional
  viewing movie, with ffmpeg installed separately:
    ffmpeg -framerate 30 -i frames/%06d.png -c:v libx264 -crf 18 -pix_fmt yuv420p preview.mp4
  That frame rate assumes the scheduled snapshots have the intended playback
  cadence. Uneven simulation-time spacing needs an explicit timing decision.

Runtime paths
  run.sh directs temporary/cache files below this package's runtime/ directory.
  Override GRMHD_RENDER_RUNTIME for a node-local scratch directory if desired.
  Set PYTHON to an interpreter path if the active python is not the right one.
  The package has no workstation-specific paths or symlinks to the old project.
