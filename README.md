# SplineMaker

SplineMaker is an interactive C++/OpenGL tool for reconstructing smooth three-dimensional particle trajectories from discrete time-step files. It groups measurements by particle id, fits independent cubic B-splines for `x(t)`, `y(t)`, and `z(t)`, and evaluates analytic velocity and acceleration from spline derivatives.

## Input Data

Place a data directory named `DA_ppp_0_005` next to the executable or run the executable from a working directory that contains that folder.

Each `*.dat` file is treated as one time step. Files are sorted alphabetically, so their names must preserve chronological order. Each file must contain three header lines followed by rows in this format:

```text
x y z particle_id
```

Rows with the same `particle_id` across files are assembled into one trajectory:

```text
P_i = { (x_k, y_k, z_k, t_k) }
```

where `t_k` is the file index in the sorted sequence.

## Trajectory Model

For every coordinate component, SplineMaker fits a cubic uniform B-spline

```text
r(t) = sum_j c_j B_3(t - j)
```

independently for:

```text
x(t), y(t), z(t)
```

The cubic basis `B_3` has compact support on `[-2, 2]`, so only four neighboring control points contribute to each evaluation. The interpolating solve is a tridiagonal linear system for the spline control points `c_j`; it is solved with Thomas elimination.

Velocity and acceleration are obtained analytically:

```text
v(t) = dr/dt
a(t) = d^2r/dt^2
```

The code evaluates these by differentiating the spline basis rather than by finite differences, which makes the derived quantities smoother than direct frame-to-frame differencing.

## Optional Smoothing Optimization

The `V` command refines spline control points by minimizing a smoothness-penalized objective:

```text
J(c) = sum_i |S(t_i; c) - p_i|^2
       + lambda * sum_j |c_{j-1} - 3c_j + 3c_{j+1} - c_{j+2}|^2
```

The first term keeps the spline near measured particle positions. The second term penalizes the discrete third derivative of the control polygon, suppressing high-frequency noise while preserving the large-scale trajectory.

## Field Reconstruction

The project also contains experimental grid-field reconstruction code. It maps scattered particle velocities and accelerations onto a Cartesian grid using local weighted least squares. Around each query point, it fits a quadratic scalar model

```text
f(x, y, z) ~= f_0 + grad(f) . r + 1/2 r^T H r
```

with Gaussian distance weights. The fitted coefficients provide field values, gradients, and Hessian terms. The pressure branch forms a Poisson-style reconstruction from acceleration divergence:

```text
Delta p = -rho div(a)
```

with boundary variants of the same local least-squares fit. This branch is useful for research experiments but is not the primary trajectory-viewing workflow.

## Visualization

The OpenGL window shows measured particle points, raw polylines, fitted splines, optional acceleration vectors, axes, status text, controls, and a color scale. Spline colors encode velocity magnitude by default; when acceleration vectors are enabled, the scale switches to acceleration magnitude.

## Controls

| Key / Mouse | Action |
|---|---|
| Mouse drag | Rotate the scene |
| `W` / `S` | Move camera forward / backward |
| `A` / `D` | Strafe camera left / right |
| `Q` / `E` | Move camera up / down |
| `1` | Toggle measured particle points |
| `2` | Toggle raw trajectory polylines |
| `3` | Toggle fitted splines |
| `4` | Toggle velocity layer placeholder |
| `5` | Toggle acceleration vectors |
| `9` / `0` | Previous / next time step |
| `[` / `]` | Decrease / increase color and vector scale |
| `I` / `O` | Move grid slice in X |
| `K` / `L` | Move grid slice in Y |
| `,` / `.` | Move grid slice in Z |
| `B` | Smooth grid fields |
| `M` | Run pressure Poisson iterations |
| `V` | Optimize all splines |
| `Z` | Calculate error metrics |
| `X` | Save reconstructed tracks, velocities, and accelerations |
| Space | Toggle continuous redraw |

## Build

The project uses qmake, C++17, OpenGL, GLU, and GLUT.

```text
qmake SplineMaker.pro
nmake
```

On Windows with the bundled 32-bit `glut32.lib`, build with an x86 MSVC developer environment. The runtime `glut32.dll` must be available next to the executable or through the system DLL search path.

## Source Layout

| File | Purpose |
|---|---|
| `main.cpp` | Application state, OpenGL rendering, interaction, field experiments |
| `fileloader.*` | Data-file discovery and particle-track loading |
| `splineinterpolator.*` | Cubic B-spline fitting, evaluation, velocity, acceleration, smoothing optimization |
| `leastsquaressolver.*` | Local weighted least-squares derivative and Poisson helpers |
| `tools.*` | Timing, color mapping, finite-difference helper coefficients |
| `globals.*` | Shared constants and viewer state |

## Notes

SplineMaker assumes monotonically ordered input snapshots and at least four samples for a trajectory to receive a valid cubic spline. Shorter tracks are loaded but skipped by spline rendering and optimization.
