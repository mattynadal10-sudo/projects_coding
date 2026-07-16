# Rope Simulation Projects (`main.cpp` and `main_ext.cpp`)

This repository contains two closely related Open Inventor / SoQt programs:

- `main.cpp`: baseline 3-D rope simulation.
- `main_ext.cpp`: extended rope simulation with a holonomic distance constraint, right-side wall collision, and distance logging.

Both programs model a rope as a chain of `N = 13` particles connected by linear spring-damper elements and render the rope as both discrete spheres and a smooth `SoNurbsCurve`.

---

## 1. Dynamics formulation

### 1.1 Common particle model

Each rope node is treated as a point mass with state

- position `q_i in R^3`
- velocity `v_i in R^3`

for `i = 0, ..., 12`.

The first particle is fixed, and the remaining particles are free.

For each free particle, the equation of motion is the sum of:

- gravity
- air drag
- wind force
- spring-damper forces from adjacent particles

In code form, the total force on particle `i` is

`f_i = f_gravity + f_drag + f_wind + sum(f_spring-damper)`

with

- gravity: `f_gravity = (0, -m g, 0)`
- drag: `f_drag = -c_air v_i`
- spring between adjacent particles `i` and `j`:
  - `d = q_j - q_i`
  - `L = ||d||`
  - `dir = d / L`
  - spring force magnitude: `k (L - L0)`
  - damper magnitude: `b ((v_j - v_i) · dir)`
  - total spring-damper force: `(k (L - L0) + b ((v_j - v_i) · dir)) dir`

where

- `m` is the particle mass
- `g` is gravitational acceleration
- `k` is spring stiffness
- `b` is the spring damping coefficient
- `L0` is the nominal segment rest length.

### 1.2 Numerical integrator

Both files use a **semi-implicit / symplectic Euler** method.

For each time step `dt`:

1. assemble all forces
2. compute acceleration from Newton’s second law
3. update velocity first
4. update position using the new velocity

That is,

- `a_i = f_i / m`
- `v_i^(n+1) = v_i^n + dt a_i`
- `q_i^(n+1) = q_i^n + dt v_i^(n+1)`

This integrator is explicit, inexpensive, and more stable for mechanical systems than standard forward Euler.

### 1.3 Visualization

Both programs visualize the rope in two ways:

- each particle is rendered as a small `SoSphere`
- the full rope is rendered as a smooth cubic `SoNurbsCurve`

The NURBS curve is built from the particle positions with repeated endpoints to obtain a clamped smooth rope-like shape.

---

## 2. `main.cpp`: baseline model

### 2.1 Baseline dynamics

`main.cpp` implements the base mass-spring-damper rope with:

- one fixed endpoint
- gravity
- linear air drag
- time-dependent wind
- nearest-neighbor springs and dampers

The initial condition is chosen as a **static hanging approximation**: each segment is initialized with an extension proportional to the weight of the masses below it. This reduces the initial transient and starts the rope close to an equilibrium hanging shape.

### 2.2 Wind-force model in `main.cpp`

Wind is toggled with the keyboard and is spatially and temporally varying. The model applies a force with:

- oscillatory `+x` component
- oscillatory `-z` component
- phase variation along the rope index

In the code:

- `fx = 40 + 18 sin(2.2 t + phase)`
- `fz = -28 - 12 cos(1.7 t + 0.4 phase)`
- `phase = 0.2 i`

This produces a nonuniform excitation along the rope, creating visible 3-D motion.

### 2.3 User interaction in `main.cpp`

- `W`: toggle wind on/off
- `R`: reset the rope to the initial hanging configuration
- `Esc`: quit

---

## 3. `main_ext.cpp`: extended model

`main_ext.cpp` extends the baseline model with three additional features:

1. a holonomic distance constraint between particles 5 and 6
2. a right-side wall collision model
3. distance logging to `distance.txt`

### 3.1 Holonomic distance constraint

The extended file enforces a fixed distance between particles

- `C_A = 4`
- `C_B = 5`

which correspond to particles 5 and 6 in 1-based indexing.

The constraint is

`C(q) = ||q_A - q_B||^2 - L_c^2 = 0`

where `L_c` is initialized from the starting distance between those particles.

The code uses a **single Lagrange multiplier** `lambda` to impose this constraint at the velocity level. The scalar multiplier is solved from the constraint Jacobian and a Baumgarte stabilization term:

- predict unconstrained velocities
- evaluate `J v*`
- compute a scalar `lambda`
- apply the corresponding corrective velocity impulse

A position projection step is then applied to reduce residual constraint drift after integration.

This gives:

- a constrained rope segment length
- reduced drift over time
- more stable long-run behavior than a purely penalty-based approach.

### 3.2 Wall collision model

The extended simulation places a textured vertical wall on the **right-hand side** at `x = WALL_X`.

Collision handling proceeds by:

1. checking whether a particle has penetrated the wall plane beyond its allowed radius
2. reflecting its position back across the wall limit
3. reversing the normal component of velocity with restitution
4. damping the tangential components slightly
5. re-projecting the distance constraint afterward

The implemented collision parameters are:

- restitution coefficient `WALL_RESTITUTION`
- tangential damping `WALL_TANGENTIAL_DAMP`
- collision gap `WALL_EPS`
- separating kick `WALL_SEPARATING_SPEED`

This produces a visible rebound instead of a particle becoming glued to the wall.

### 3.3 Wind-force model in `main_ext.cpp`

In the extended file, the wind is simplified to a strong deterministic push in the positive `x` direction:

- `fx = 65 + 5 i`
- `fy = 0`
- `fz = 0`

This model is intentionally chosen to drive the rope toward the wall and generate collisions.

### 3.4 Distance logging

After each simulation step, the code appends the current constrained distance

`||q_5 - q_6||`

to `distance.txt`.

This makes it possible to evaluate how well the constraint is maintained during the simulation.

### 3.5 User interaction in `main_ext.cpp`

- `W`: toggle wind on/off
- `R`: reset rope
- `Esc`: quit

---

## 4. Optional / extra features implemented

Across the two files, the following optional enhancements are present:

- smooth NURBS rope rendering instead of only polyline segments
- near-equilibrium initialization for the baseline rope
- Qt event filter for reliable keyboard handling in the extended version
- textured wall in `main_ext.cpp`
- persistent logging of constrained distance in `main_ext.cpp`
- Baumgarte-stabilized constraint solve plus position projection in `main_ext.cpp`

---

## 5. Compiler, libraries, and build requirements

### Compiler

The tested build setup is:

- `clang++` with C++17

### Required libraries

- Coin3D / Open Inventor
- SoQt
- Qt (via SoQt)
- `pkg-config`

On the user machine shown during development, SoQt flags were obtained with:

```bash
/opt/homebrew/bin/pkg-config --cflags --libs SoQt
```

If needed, Eigen can also be included on the command line, although these two source files do not directly use it.

---

## 6. Exact compile commands

### 6.1 Compile `main.cpp`

```bash
clang++ -std=c++17 -g -O0 main.cpp -o main $(pkg-config --cflags --libs SoQt)
```

On the Homebrew/macOS setup used during development, the explicit form is:

```bash
clang++ -std=c++17 -g -O0 \
  /Users/mattemasce/Documents/coin_3d_setup_folder/project_1/main.cpp \
  -o /Users/mattemasce/Documents/coin_3d_setup_folder/project_1/main \
  $(/opt/homebrew/bin/pkg-config --cflags --libs SoQt)
```

### 6.2 Compile `main_ext.cpp`

```bash
clang++ -std=c++17 -g -O0 main_ext.cpp -o main_ext $(pkg-config --cflags --libs SoQt)
```

On the Homebrew/macOS setup used during development, the explicit form is:

```bash
clang++ -std=c++17 -g -O0 \
  /Users/mattemasce/Documents/coin_3d_setup_folder/project_1/main_ext.cpp \
  -o /Users/mattemasce/Documents/coin_3d_setup_folder/project_1/main_ext \
  $(/opt/homebrew/bin/pkg-config --cflags --libs SoQt)
```

If `pkg-config` is not on the shell path, replace `pkg-config` with the full path, for example:

```bash
/opt/homebrew/bin/pkg-config --cflags --libs SoQt
```

---

## 7. Run instructions

### Run the baseline program

```bash
./main
```

### Run the extended program

```bash
./main_ext
```

Notes:

- Click once inside the render window if keyboard events do not register immediately.
- In `main_ext.cpp`, the simulation appends measurements to `distance.txt` in the current working directory.
- The extended version also expects a texture file named `wall.jpg` in the working directory if the wall texture is to appear correctly.

---

## 8. Summary

- `main.cpp` is the base 13-particle mass-spring-damper rope with gravity, drag, wind, NURBS rendering, and interactive wind toggle.
- `main_ext.cpp` adds a Lagrange-multiplier-based distance constraint, wall collision, and distance logging, while preserving the same real-time visualization structure.

