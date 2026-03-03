# Universe of Goo Milestone II

## Team
- Name: Will Suan
- Teammate: none

## Build Instructions
Use the normal CMake flow from the project root:

```bash
mkdir -p build
cd build
cmake ..
make -j4
```

For optimized builds:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```

## Run Instructions
From `build/`:

```bash
./bin/goo2
```

## Test Instructions
This project includes built-in non-GUI tests so constraint logic can be checked from terminal only:

```bash
./bin/goo2 --run-tests
```

The test suite covers:
- flexible rod creation (inert particles, spring masses, and bending hinges)
- rigid rod step-and-project projection quality
- saw deletion behavior for flexible rods and stencil cleanup
- constrained Lagrange multiplier smoke test (finite states after integration)

## Milestone II Feature Checklist
- Rigid rod creation mode is supported through `Connector Type` dropdown.
- Rod rendering is visually distinct from springs.
- Rods are deleted correctly by saws, endpoint deletion, and scene reset.
- Penalty method is implemented with rod constraint potential and force.
- Step-and-project is implemented with Newton solve over `(q, lambda)`.
- Lagrange multiplier method is implemented with Newton solve on lambda.
- Rigid rod lambdas are packed/unpacked in global configuration vectors.
- Flexible rods are implemented by splitting into inert internal particles and unsnappable springs.
- Flexible rod spring parameters use:
  - rest length `L`
  - stiffness `ks / L`
  - mass `rho * L`
- Bending hinges are created between consecutive rod segments with stiffness `2*kb/(L1+L2)`.
- Elastic bending force is implemented and toggleable through `Bending Enabled`.
- Bending hinges are pruned automatically whenever springs/particles are deleted by saws or culling.

## Aesthetic Improvements Beyond Baseline
- Updated color palette for clearer visual grouping:
  - snappable springs, unsnappable rope springs, rigid rods, and bending-active springs each have unique color styling
  - inert particles are visually different from regular and fixed particles
- Particle color now scales with effective mass, so mass changes from rod springs are easier to read visually.
- Connector widths are tuned by type to make rods feel more structural and rope segments easier to track.
- Floor and object palette was adjusted for better contrast and less visual clutter.

## Known Limitations / Notes
- Very short flexible rod segments can still destabilize due to stiff bending forces (as expected from the spec note).
- The Lagrange multiplier mode uses Newton with zero initialization each frame; if the system gets extremely stiff, it may need smaller timestep for best behavior.
