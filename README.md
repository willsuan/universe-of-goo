# Universe of Goo

Physics simulation milestones for a computer graphics course, inspired by
*World of Goo*. Each milestone is a self-contained C++ / Eigen / Polyscope
project with its own build, tests, and write-up.

## Layout

```
universe-of-goo/
├── milestone1/        # Mass-spring simulation
│   ├── CMakeLists.txt
│   ├── README.md      # Full feature list, build + run + test instructions
│   ├── cmake/
│   ├── deps/polyscope/
│   ├── src/           # main.cpp, SimParameters.h, SceneObjects.h
│   ├── tests/         # 73 unit tests, 2229 assertions
│   └── goo1.pdf       # Spec
└── milestone2/        # Rigid + flexible rods, constraint methods
    ├── CMakeLists.txt
    ├── README.md      # Full feature list, build + run + test instructions
    ├── cmake/
    ├── deps/polyscope/
    ├── src/           # main.cpp
    └── tests/         # Non-GUI constraint tests
```

## Milestone I — Mass-Spring Simulation

2D particles connected by Hookean springs, with gravity, viscous damping,
and a floor. Four integrators (Explicit Euler, Velocity Verlet, Implicit
Euler, Implicit Midpoint), analytic force Jacobians for the implicit solvers,
Cauchy-strain spring snapping, and a saw tool that cuts springs and removes
particles with correct index remapping.

73 unit tests cover geometry, mass matrix, gravity, spring force + Hessian
(finite-difference verified), damping, floor penalty, all four integrators,
fixed particles, snapping, saw collisions, and energy conservation.

See [`milestone1/README.md`](milestone1/README.md) for build and run
instructions.

## Milestone II — Rigid and Flexible Rods

Adds rigid-rod connectors implemented three different ways:

- **Penalty method** — rod constraint potential + force added to the system
- **Step-and-project** — unconstrained step followed by Newton solve over
  `(q, lambda)` that projects back onto the constraint manifold
- **Lagrange multipliers** — Newton solve on `lambda` each step

Flexible rods are modeled by splitting the rod into inert internal particles
plus unsnappable springs, with bending hinges between consecutive segments.
Spring parameters scale with rest length (`ks/L`, mass `rho*L`); hinge
stiffness is `2*kb/(L1+L2)`. Elastic bending force is toggleable.

See [`milestone2/README.md`](milestone2/README.md) for build and run
instructions, along with the full feature checklist.

## Build (either milestone)

```bash
cd milestone1   # or milestone2
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
./bin/goo1                # or ./bin/goo2
./bin/goo1_tests          # or ./bin/goo2 --run-tests
```
