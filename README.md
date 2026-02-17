# Universe of Goo

A 2D mass-spring physics simulation inspired by World of Goo, built with C++, Eigen, and Polyscope.

## Features

- **Particle system** with mass, position, velocity, and fixed/free DOFs
- **Spring connectors** with automatic creation based on proximity, Hookean force model (stiffness = k_b / rest length), and Cauchy strain-based snapping
- **Force computation** including gravity, spring forces, viscous damping, and floor penalty forces, with analytic force Jacobian (Hessian) for implicit solvers
- **Four numerical integrators**: Explicit Euler, Velocity Verlet, Implicit Euler (Newton's method), and Implicit Midpoint (Newton's method)
- **Saw tool** for cutting springs and deleting particles, with proper index remapping
- **Interactive GUI** via Polyscope with click-to-add particles, saw placement, and real-time parameter tuning

## Build

```
mkdir build
cd build
cmake ..
make -j4
```

For Release mode (compiler optimizations):
```
cmake -DCMAKE_BUILD_TYPE=Release ..
```

## Run

```
./build/bin/goo1
```

**Controls:**
- Left-click to add particles (springs auto-connect to nearby particles)
- Switch to saw mode in the GUI to place saws that cut springs/particles
- Toggle gravity, springs, damping, floor in the sidebar
- Switch integrators and adjust time step, stiffness, etc. in real time

## Tests

73 unit tests covering all physics functions (2229 assertions):

```
cd build
cmake --build . --target goo1_tests
./bin/goo1_tests
```

Test categories: geometry helpers, particle/spring creation, configuration vectors, mass matrix, gravity, spring forces (including finite-difference verification of force and Hessian), damping, floor penalty, all four integrators, fixed particles, spring snapping, saw collision with index remapping, out-of-bounds deletion, energy conservation, and integrator stability.

## Project Structure

```
src/
  main.cpp          - All simulation logic and Polyscope GUI
  SimParameters.h   - Simulation parameters (time step, stiffness, etc.)
  SceneObjects.h    - Particle, Spring, Saw, Connector data structures
tests/
  test_main.cpp     - 73 unit tests with custom minimal test framework
deps/
  polyscope/        - Polyscope visualization library
cmake/
  libigl.cmake      - libigl/Eigen dependency fetching
```

