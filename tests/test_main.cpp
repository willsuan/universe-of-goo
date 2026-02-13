// Universe of Goo — Milestone I Test Suite
//
// This file tests all implemented physics functions by including main.cpp
// with the real main() renamed via a preprocessor macro, then defining our
// own main() that runs a battery of unit tests.
//
// Test categories:
//   1. Geometry helpers      (point-to-line, point-to-segment distance)
//   2. Particle/spring creation (addParticle, spring connectivity)
//   3. Configuration vectors (buildConfiguration, unbuildConfiguration)
//   4. Mass matrix           (computeMassInverse)
//   5. Gravity force         (computeForceAndHessian with gravity only)
//   6. Spring force          (computeForceAndHessian with springs only)
//   7. Damping force         (computeForceAndHessian with damping only)
//   8. Floor force           (computeForceAndHessian with floor only)
//   9. Explicit Euler        (numericalIntegration)
//  10. Velocity Verlet       (numericalIntegration)
//  11. Implicit Euler        (numericalIntegration)
//  12. Implicit Midpoint     (numericalIntegration)
//  13. Spring snapping       (pruneOverstrainedSprings)
//  14. Saw collision         (deleteSawedObjects)
//  15. Fixed particle DOFs   (fixed particles don't move)
//  16. Energy conservation   (implicit midpoint preserves energy)
//  17. Reset simulation      (initSimulation clears everything)

// Rename the real main() so we can define our own test runner
#define main original_main
#include "../src/main.cpp"
#undef main

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <functional>

// Minimal test framework

static int g_tests_passed = 0;
static int g_tests_failed = 0;
static int g_assertions_passed = 0;
static int g_assertions_failed = 0;
static std::string g_current_test;
static std::vector<std::string> g_failures;

#define ASSERT_TRUE(expr)                                                     \
    do {                                                                      \
        if (!(expr)) {                                                        \
            std::ostringstream ss;                                            \
            ss << "  FAIL: " << #expr << "  (" << __FILE__ << ":" << __LINE__ << ")"; \
            std::cerr << ss.str() << std::endl;                              \
            g_assertions_failed++;                                            \
            g_failures.push_back(g_current_test + " -> " + ss.str());        \
        } else {                                                              \
            g_assertions_passed++;                                            \
        }                                                                     \
    } while (0)

#define ASSERT_FALSE(expr) ASSERT_TRUE(!(expr))

#define ASSERT_EQ(a, b)                                                       \
    do {                                                                      \
        if ((a) != (b)) {                                                     \
            std::ostringstream ss;                                            \
            ss << "  FAIL: " << #a << " == " << #b << "  (got " << (a)       \
               << " vs " << (b) << ")  (" << __FILE__ << ":" << __LINE__ << ")"; \
            std::cerr << ss.str() << std::endl;                              \
            g_assertions_failed++;                                            \
            g_failures.push_back(g_current_test + " -> " + ss.str());        \
        } else {                                                              \
            g_assertions_passed++;                                            \
        }                                                                     \
    } while (0)

#define ASSERT_NEAR(a, b, tol)                                                \
    do {                                                                      \
        double _a = (a), _b = (b), _t = (tol);                               \
        if (std::abs(_a - _b) > _t) {                                        \
            std::ostringstream ss;                                            \
            ss << "  FAIL: |" << #a << " - " << #b << "| <= " << _t          \
               << "  (got |" << _a << " - " << _b << "| = "                  \
               << std::abs(_a - _b) << ")  (" << __FILE__ << ":" << __LINE__ << ")"; \
            std::cerr << ss.str() << std::endl;                              \
            g_assertions_failed++;                                            \
            g_failures.push_back(g_current_test + " -> " + ss.str());        \
        } else {                                                              \
            g_assertions_passed++;                                            \
        }                                                                     \
    } while (0)

#define ASSERT_GT(a, b)                                                       \
    do {                                                                      \
        double _a = (a), _b = (b);                                            \
        if (!(_a > _b)) {                                                     \
            std::ostringstream ss;                                            \
            ss << "  FAIL: " << #a << " > " << #b << "  (got " << _a         \
               << " vs " << _b << ")  (" << __FILE__ << ":" << __LINE__ << ")"; \
            std::cerr << ss.str() << std::endl;                              \
            g_assertions_failed++;                                            \
            g_failures.push_back(g_current_test + " -> " + ss.str());        \
        } else {                                                              \
            g_assertions_passed++;                                            \
        }                                                                     \
    } while (0)

#define ASSERT_LT(a, b) ASSERT_GT(b, a)

void runTest(const std::string& name, std::function<void()> fn)
{
    g_current_test = name;
    int failed_before = g_assertions_failed;
    std::cout << "[ RUN      ] " << name << std::endl;

    // Reset simulation state before each test
    initSimulation();
    params_ = SimParameters();  // Reset to defaults

    try {
        fn();
    } catch (const std::exception& e) {
        std::cerr << "  EXCEPTION: " << e.what() << std::endl;
        g_assertions_failed++;
        g_failures.push_back(name + " -> EXCEPTION: " + e.what());
    }

    if (g_assertions_failed == failed_before) {
        std::cout << "[       OK ] " << name << std::endl;
        g_tests_passed++;
    } else {
        std::cout << "[  FAILED  ] " << name << std::endl;
        g_tests_failed++;
    }
}

// Helper: reset everything and return to a clean slate
void resetAll()
{
    initSimulation();
    params_ = SimParameters();
}

// TEST 1: Point-to-line distance
void test_distancePointToLine_perpendicular()
{
    // Point (1, 1) to line along x-axis (y=0): distance should be 1
    Eigen::Vector2d p(1.0, 1.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(2.0, 0.0);
    ASSERT_NEAR(distancePointToLine(p, q1, q2), 1.0, 1e-10);
}

void test_distancePointToLine_onLine()
{
    // Point on the line: distance should be 0
    Eigen::Vector2d p(1.0, 0.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(2.0, 0.0);
    ASSERT_NEAR(distancePointToLine(p, q1, q2), 0.0, 1e-10);
}

void test_distancePointToLine_diagonal()
{
    // Point (0, 1) to the line y = x: distance = 1/sqrt(2)
    Eigen::Vector2d p(0.0, 1.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(1.0, 1.0);
    ASSERT_NEAR(distancePointToLine(p, q1, q2), 1.0 / std::sqrt(2.0), 1e-10);
}

void test_distancePointToLine_negativeProjection()
{
    // Point behind q1: infinite line still gives perpendicular distance
    Eigen::Vector2d p(-1.0, 1.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(2.0, 0.0);
    ASSERT_NEAR(distancePointToLine(p, q1, q2), 1.0, 1e-10);
}

void test_distancePointToLine_degenerateSegment()
{
    // Degenerate line (q1 == q2): distance is just point-to-point
    Eigen::Vector2d p(3.0, 4.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(0.0, 0.0);
    ASSERT_NEAR(distancePointToLine(p, q1, q2), 5.0, 1e-10);
}

// TEST 2: Point-to-segment distance
void test_distancePointToSegment_perpendicular()
{
    // Same as line test when projection falls on segment
    Eigen::Vector2d p(1.0, 1.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(2.0, 0.0);
    ASSERT_NEAR(distancePointToSegment(p, q1, q2), 1.0, 1e-10);
}

void test_distancePointToSegment_clampToQ1()
{
    // Point projects behind q1 -> closest point is q1
    Eigen::Vector2d p(-1.0, 1.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(2.0, 0.0);
    double expected = std::sqrt(2.0);  // dist to q1
    ASSERT_NEAR(distancePointToSegment(p, q1, q2), expected, 1e-10);
}

void test_distancePointToSegment_clampToQ2()
{
    // Point projects beyond q2 -> closest point is q2
    Eigen::Vector2d p(3.0, 1.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(2.0, 0.0);
    double expected = std::sqrt(2.0);  // dist to q2
    ASSERT_NEAR(distancePointToSegment(p, q1, q2), expected, 1e-10);
}

void test_distancePointToSegment_endpoint()
{
    // Point is exactly at an endpoint
    Eigen::Vector2d p(0.0, 0.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(1.0, 0.0);
    ASSERT_NEAR(distancePointToSegment(p, q1, q2), 0.0, 1e-10);
}

void test_distancePointToSegment_midpoint()
{
    // Point directly above midpoint of segment
    Eigen::Vector2d p(0.5, 2.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(1.0, 0.0);
    ASSERT_NEAR(distancePointToSegment(p, q1, q2), 2.0, 1e-10);
}

// TEST 3: Particle and spring creation
void test_addParticle_single()
{
    // Adding a single particle: no springs created
    addParticle(0.0, 0.0);
    ASSERT_EQ((int)particles_.size(), 1);
    ASSERT_EQ((int)connectors_.size(), 0);
    ASSERT_NEAR(particles_[0].pos[0], 0.0, 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], 0.0, 1e-10);
    ASSERT_NEAR(particles_[0].vel[0], 0.0, 1e-10);
    ASSERT_NEAR(particles_[0].vel[1], 0.0, 1e-10);
}

void test_addParticle_twoClose()
{
    // Two particles within maxSpringDist: one spring created
    params_.maxSpringDist = 0.5;
    addParticle(0.0, 0.0);
    addParticle(0.1, 0.0);
    ASSERT_EQ((int)particles_.size(), 2);
    ASSERT_EQ((int)connectors_.size(), 1);

    Spring* s = dynamic_cast<Spring*>(connectors_[0]);
    ASSERT_TRUE(s != nullptr);
    ASSERT_NEAR(s->restlen, 0.1, 1e-10);
}

void test_addParticle_twoFar()
{
    // Two particles beyond maxSpringDist: no spring created
    params_.maxSpringDist = 0.05;
    addParticle(0.0, 0.0);
    addParticle(0.1, 0.0);
    ASSERT_EQ((int)particles_.size(), 2);
    ASSERT_EQ((int)connectors_.size(), 0);
}

void test_addParticle_springStiffness()
{
    // Spring stiffness should be springStiffness / restLength
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 200.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    ASSERT_EQ((int)connectors_.size(), 1);

    Spring* s = dynamic_cast<Spring*>(connectors_[0]);
    ASSERT_TRUE(s != nullptr);
    ASSERT_NEAR(s->stiffness, 200.0 / 0.5, 1e-10);  // k = kb / L
    ASSERT_NEAR(s->restlen, 0.5, 1e-10);
}

void test_addParticle_multipleConnections()
{
    // Three particles in a triangle: adding 3rd should create 2 springs (to 1st and 2nd)
    params_.maxSpringDist = 0.5;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);   // 1 spring to particle 0
    addParticle(0.1, 0.15);  // within 0.25 of both -> 2 springs
    ASSERT_EQ((int)particles_.size(), 3);
    ASSERT_EQ((int)connectors_.size(), 3);  // 0-1, 2-0, 2-1
}

void test_addParticle_fixedParticle()
{
    // Fixed particles should have infinite mass and fixed flag
    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    ASSERT_TRUE(particles_[0].fixed);
    ASSERT_TRUE(std::isinf(particles_[0].mass));
}

void test_addParticle_massPreserved()
{
    // Each particle should have the mass set at creation time
    params_.particleMass = 3.5;
    addParticle(0.0, 0.0);
    params_.particleMass = 7.0;
    addParticle(1.0, 0.0);
    ASSERT_NEAR(particles_[0].mass, 3.5, 1e-10);
    ASSERT_NEAR(particles_[1].mass, 7.0, 1e-10);
}

// TEST 4: Configuration vector building
void test_buildConfiguration_basic()
{
    addParticle(1.0, 2.0);
    addParticle(3.0, 4.0);
    particles_[0].vel = Eigen::Vector2d(0.5, 0.6);
    particles_[1].vel = Eigen::Vector2d(0.7, 0.8);

    Eigen::VectorXd q, qprev, qdot;
    buildConfiguration(q, qprev, qdot);

    ASSERT_EQ(q.size(), 4);
    ASSERT_NEAR(q[0], 1.0, 1e-10);
    ASSERT_NEAR(q[1], 2.0, 1e-10);
    ASSERT_NEAR(q[2], 3.0, 1e-10);
    ASSERT_NEAR(q[3], 4.0, 1e-10);
    ASSERT_NEAR(qdot[0], 0.5, 1e-10);
    ASSERT_NEAR(qdot[1], 0.6, 1e-10);
    ASSERT_NEAR(qdot[2], 0.7, 1e-10);
    ASSERT_NEAR(qdot[3], 0.8, 1e-10);
}

void test_buildConfiguration_prevpos()
{
    addParticle(1.0, 2.0);
    particles_[0].prevpos = Eigen::Vector2d(0.9, 1.9);

    Eigen::VectorXd q, qprev, qdot;
    buildConfiguration(q, qprev, qdot);

    ASSERT_NEAR(qprev[0], 0.9, 1e-10);
    ASSERT_NEAR(qprev[1], 1.9, 1e-10);
}

void test_unbuildConfiguration_basic()
{
    addParticle(0.0, 0.0);
    addParticle(0.0, 0.0);

    Eigen::VectorXd q(4), qdot(4);
    q << 1.0, 2.0, 3.0, 4.0;
    qdot << 5.0, 6.0, 7.0, 8.0;
    unbuildConfiguration(q, qdot);

    ASSERT_NEAR(particles_[0].pos[0], 1.0, 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], 2.0, 1e-10);
    ASSERT_NEAR(particles_[1].pos[0], 3.0, 1e-10);
    ASSERT_NEAR(particles_[1].pos[1], 4.0, 1e-10);
    ASSERT_NEAR(particles_[0].vel[0], 5.0, 1e-10);
    ASSERT_NEAR(particles_[0].vel[1], 6.0, 1e-10);
}

void test_unbuildConfiguration_savesPrevpos()
{
    addParticle(1.0, 2.0);

    Eigen::VectorXd q(2), qdot(2);
    q << 1.5, 2.5;
    qdot << 0.0, 0.0;
    unbuildConfiguration(q, qdot);

    // prevpos should be the OLD position (1.0, 2.0)
    ASSERT_NEAR(particles_[0].prevpos[0], 1.0, 1e-10);
    ASSERT_NEAR(particles_[0].prevpos[1], 2.0, 1e-10);
    // pos should be the NEW position (1.5, 2.5)
    ASSERT_NEAR(particles_[0].pos[0], 1.5, 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], 2.5, 1e-10);
}

void test_buildUnbuild_roundtrip()
{
    // Build and unbuild should be inverses of each other
    addParticle(1.1, 2.2);
    addParticle(3.3, 4.4);
    particles_[0].vel = Eigen::Vector2d(5.5, 6.6);
    particles_[1].vel = Eigen::Vector2d(7.7, 8.8);

    Eigen::VectorXd q, qprev, qdot;
    buildConfiguration(q, qprev, qdot);

    // Modify and unbuild
    unbuildConfiguration(q, qdot);

    ASSERT_NEAR(particles_[0].pos[0], 1.1, 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], 2.2, 1e-10);
    ASSERT_NEAR(particles_[1].vel[0], 7.7, 1e-10);
    ASSERT_NEAR(particles_[1].vel[1], 8.8, 1e-10);
}

// TEST 5: Mass matrix
void test_massMatrix_singleParticle()
{
    params_.particleMass = 2.0;
    addParticle(0.0, 0.0);

    Eigen::SparseMatrix<double> Minv;
    computeMassInverse(Minv);

    ASSERT_EQ(Minv.rows(), 2);
    ASSERT_EQ(Minv.cols(), 2);
    ASSERT_NEAR(Minv.coeff(0, 0), 0.5, 1e-10);
    ASSERT_NEAR(Minv.coeff(1, 1), 0.5, 1e-10);
    ASSERT_NEAR(Minv.coeff(0, 1), 0.0, 1e-10);
}

void test_massMatrix_fixedParticle()
{
    params_.particleFixed = true;
    addParticle(0.0, 0.0);

    Eigen::SparseMatrix<double> Minv;
    computeMassInverse(Minv);

    // Fixed particle: M^{-1} = 0
    ASSERT_NEAR(Minv.coeff(0, 0), 0.0, 1e-10);
    ASSERT_NEAR(Minv.coeff(1, 1), 0.0, 1e-10);
}

void test_massMatrix_mixedParticles()
{
    params_.particleMass = 4.0;
    params_.particleFixed = false;
    params_.maxSpringDist = 0.01;  // prevent springs from adding mass
    addParticle(0.0, 0.0);

    params_.particleFixed = true;
    addParticle(1.0, 0.0);

    Eigen::SparseMatrix<double> Minv;
    computeMassInverse(Minv);

    ASSERT_EQ(Minv.rows(), 4);
    ASSERT_NEAR(Minv.coeff(0, 0), 0.25, 1e-10);  // 1/4
    ASSERT_NEAR(Minv.coeff(1, 1), 0.25, 1e-10);
    ASSERT_NEAR(Minv.coeff(2, 2), 0.0, 1e-10);   // fixed
    ASSERT_NEAR(Minv.coeff(3, 3), 0.0, 1e-10);
}

void test_massMatrix_isDiagonal()
{
    params_.particleMass = 1.0;
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);
    addParticle(1.0, 0.0);
    addParticle(2.0, 0.0);

    Eigen::SparseMatrix<double> Minv;
    computeMassInverse(Minv);

    // Check all off-diagonal entries are zero
    for (int i = 0; i < Minv.rows(); i++)
        for (int j = 0; j < Minv.cols(); j++)
            if (i != j)
                ASSERT_NEAR(Minv.coeff(i, j), 0.0, 1e-10);
}

// TEST 6: Gravity force
void test_gravity_singleParticle()
{
    params_.gravityEnabled = true;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.particleMass = 2.0;
    params_.gravityG = -9.8;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 0.0);

    Eigen::VectorXd q(2), qprev(2), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0;
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    // F_x = 0, F_y = m * g = 2.0 * (-9.8) = -19.6
    ASSERT_NEAR(F[0], 0.0, 1e-10);
    ASSERT_NEAR(F[1], 2.0 * (-9.8), 1e-10);
}

void test_gravity_fixedParticle()
{
    params_.gravityEnabled = true;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.particleFixed = true;

    addParticle(0.0, 0.0);

    Eigen::VectorXd q(2), qprev(2), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0;
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    // Fixed particle: force should be zero
    ASSERT_NEAR(F[0], 0.0, 1e-10);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
}

void test_gravity_zeroHessian()
{
    params_.gravityEnabled = true;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 0.0);

    Eigen::VectorXd q(2), qprev(2), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.5;
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    // Gravity Hessian should be zero (constant force)
    ASSERT_EQ(H.nonZeros(), 0);
}

void test_gravity_disabled()
{
    params_.gravityEnabled = false;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 0.0);

    Eigen::VectorXd q(2), qprev(2), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0;
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    ASSERT_NEAR(F[0], 0.0, 1e-10);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
}

// TEST 7: Spring force
void test_springForce_atRest()
{
    // Two particles at rest length: spring force should be zero
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    Eigen::VectorXd q(4), qprev(4), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0, 0.5, 0.0;
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    // At rest length, force should be zero
    ASSERT_NEAR(F[0], 0.0, 1e-10);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
    ASSERT_NEAR(F[2], 0.0, 1e-10);
    ASSERT_NEAR(F[3], 0.0, 1e-10);
}

void test_springForce_stretched()
{
    // Stretch the spring: force should pull particles together
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);  // rest length = 0.5, stiffness = 100/0.5 = 200

    Spring* s = dynamic_cast<Spring*>(connectors_[0]);
    double k = s->stiffness;
    double L = s->restlen;

    // Place particles further apart than rest length
    Eigen::VectorXd q(4), qprev(4), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0, 1.0, 0.0;  // distance = 1.0, rest = 0.5
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    // Force on p1: k * (r - L) * dhat = 200 * (1.0 - 0.5) * (1,0) = (100, 0)
    ASSERT_NEAR(F[0], 100.0, 1e-8);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
    // Force on p2: equal and opposite
    ASSERT_NEAR(F[2], -100.0, 1e-8);
    ASSERT_NEAR(F[3], 0.0, 1e-10);
}

void test_springForce_compressed()
{
    // Compress the spring: force should push particles apart
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);  // rest length = 0.5, k = 200

    Eigen::VectorXd q(4), qprev(4), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0, 0.2, 0.0;  // distance = 0.2, rest = 0.5
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    // Force on p1: k * (r - L) * dhat = 200 * (0.2 - 0.5) * (1,0) = (-60, 0)
    // (pushed away from p2, i.e., in negative x direction)
    ASSERT_NEAR(F[0], -60.0, 1e-8);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
    ASSERT_NEAR(F[2], 60.0, 1e-8);
    ASSERT_NEAR(F[3], 0.0, 1e-10);
}

void test_springForce_newtonThirdLaw()
{
    // Spring forces should sum to zero (Newton's third law)
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.3, 0.4);  // rest length = 0.5

    Eigen::VectorXd q(4), qprev(4), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0, 0.6, 0.8;  // distance = 1.0
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    // Total force should be zero
    ASSERT_NEAR(F[0] + F[2], 0.0, 1e-10);
    ASSERT_NEAR(F[1] + F[3], 0.0, 1e-10);
}

void test_springForce_hessianSymmetric()
{
    // The Hessian matrix should be symmetric
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    Eigen::VectorXd q(4), qprev(4), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0, 0.8, 0.3;
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    Eigen::MatrixXd Hdense = Eigen::MatrixXd(H);
    for (int i = 0; i < Hdense.rows(); i++)
        for (int j = 0; j < Hdense.cols(); j++)
            ASSERT_NEAR(Hdense(i, j), Hdense(j, i), 1e-10);
}

void test_springForce_finiteDifferenceCheck()
{
    // Verify the force by finite-difference of potential energy
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    Spring* s = dynamic_cast<Spring*>(connectors_[0]);
    double k = s->stiffness;
    double L = s->restlen;

    Eigen::VectorXd q(4), qprev(4);
    q << 0.1, 0.2, 0.7, 0.4;
    qprev = q;

    Eigen::VectorXd F;
    Eigen::SparseMatrix<double> H;
    computeForceAndHessian(q, qprev, F, H);

    // Compute potential energy at q: V = (k/2) * (||p2-p1|| - L)^2
    auto computeV = [&](const Eigen::VectorXd& qq) -> double {
        Eigen::Vector2d p1(qq[0], qq[1]);
        Eigen::Vector2d p2(qq[2], qq[3]);
        double r = (p2 - p1).norm();
        return 0.5 * k * (r - L) * (r - L);
    };

    // Finite-difference gradient: F_i ≈ -(V(q + eps*e_i) - V(q - eps*e_i)) / (2*eps)
    double eps = 1e-6;
    for (int i = 0; i < 4; i++)
    {
        Eigen::VectorXd qp = q, qm = q;
        qp[i] += eps;
        qm[i] -= eps;
        double Ffd = -(computeV(qp) - computeV(qm)) / (2.0 * eps);
        ASSERT_NEAR(F[i], Ffd, 1e-4);
    }
}

void test_springForce_hessianFiniteDifferenceCheck()
{
    // Verify Hessian H = dF/dq by finite-difference of F
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    Eigen::VectorXd q(4), qprev(4);
    q << 0.1, 0.2, 0.7, 0.4;
    qprev = q;

    Eigen::VectorXd F;
    Eigen::SparseMatrix<double> H;
    computeForceAndHessian(q, qprev, F, H);
    Eigen::MatrixXd Hdense = Eigen::MatrixXd(H);

    // Finite-difference: H_{ij} ≈ (F_i(q + eps*e_j) - F_i(q - eps*e_j)) / (2*eps)
    double eps = 1e-6;
    for (int j = 0; j < 4; j++)
    {
        Eigen::VectorXd qp = q, qm = q;
        qp[j] += eps;
        qm[j] -= eps;

        Eigen::VectorXd Fp, Fm;
        Eigen::SparseMatrix<double> Hp, Hm;
        computeForceAndHessian(qp, qprev, Fp, Hp);
        computeForceAndHessian(qm, qprev, Fm, Hm);

        for (int i = 0; i < 4; i++)
        {
            double Hfd = (Fp[i] - Fm[i]) / (2.0 * eps);
            ASSERT_NEAR(Hdense(i, j), Hfd, 1e-3);
        }
    }
}

// TEST 8: Damping force
void test_dampingForce_zeroVelocity()
{
    // If both particles have zero velocity (q == qprev), damping force = 0
    params_.gravityEnabled = false;
    params_.springsEnabled = false;
    params_.dampingEnabled = true;
    params_.floorEnabled = false;
    params_.maxSpringDist = 1.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    Eigen::VectorXd q(4), qprev(4), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0, 0.5, 0.0;
    qprev = q;  // zero finite-difference velocity

    computeForceAndHessian(q, qprev, F, H);

    ASSERT_NEAR(F[0], 0.0, 1e-10);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
    ASSERT_NEAR(F[2], 0.0, 1e-10);
    ASSERT_NEAR(F[3], 0.0, 1e-10);
}

void test_dampingForce_opposesRelativeVelocity()
{
    // Particles moving apart -> damping force pulls them together
    params_.gravityEnabled = false;
    params_.springsEnabled = false;
    params_.dampingEnabled = true;
    params_.floorEnabled = false;
    params_.dampingStiffness = 2.0;
    params_.timeStep = 0.01;
    params_.maxSpringDist = 1.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    double h = params_.timeStep;
    double kdamp = params_.dampingStiffness;

    Eigen::VectorXd q(4), qprev(4), F;
    Eigen::SparseMatrix<double> H;
    // p1 moved left, p2 moved right (particles moving apart)
    q << -0.01, 0.0, 0.51, 0.0;
    qprev << 0.0, 0.0, 0.5, 0.0;

    computeForceAndHessian(q, qprev, F, H);

    // v1 = (-0.01 - 0.0)/h = -1.0, v2 = (0.51 - 0.5)/h = 1.0
    // F1 = kdamp * (v2 - v1) = 2 * (1 - (-1)) = 4.0 (pulls p1 toward p2)
    // F2 = kdamp * (v1 - v2) = 2 * (-1 - 1) = -4.0 (pulls p2 toward p1)
    ASSERT_NEAR(F[0], 4.0, 1e-8);
    ASSERT_NEAR(F[2], -4.0, 1e-8);
}

void test_dampingForce_newtonThirdLaw()
{
    params_.gravityEnabled = false;
    params_.springsEnabled = false;
    params_.dampingEnabled = true;
    params_.floorEnabled = false;
    params_.maxSpringDist = 1.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    Eigen::VectorXd q(4), qprev(4), F;
    Eigen::SparseMatrix<double> H;
    q << 0.01, 0.02, 0.52, 0.03;
    qprev << 0.0, 0.0, 0.5, 0.0;

    computeForceAndHessian(q, qprev, F, H);

    // Total force should sum to zero
    ASSERT_NEAR(F[0] + F[2], 0.0, 1e-10);
    ASSERT_NEAR(F[1] + F[3], 0.0, 1e-10);
}

// TEST 9: Floor force
void test_floorForce_aboveFloor()
{
    // Particle above floor: no floor force
    params_.gravityEnabled = false;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = true;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 0.0);  // y = 0, well above floor at y = -0.5

    Eigen::VectorXd q(2), qprev(2), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0;
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    ASSERT_NEAR(F[0], 0.0, 1e-10);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
}

void test_floorForce_belowFloor()
{
    // Particle below floor: upward restoring force
    params_.gravityEnabled = false;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = true;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, -0.6);

    Eigen::VectorXd q(2), qprev(2), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, -0.6;
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    // Floor at y = -0.5, particle at y = -0.6, penetration = 0.1
    // F_y = kfloor * 0.1 = 10000 * 0.1 = 1000 (upward)
    ASSERT_GT(F[1], 0.0);
    ASSERT_NEAR(F[1], 1000.0, 1e-6);
    ASSERT_NEAR(F[0], 0.0, 1e-10);  // no x force
}

void test_floorForce_atExactFloor()
{
    // Particle exactly at floor: no force (boundary case)
    params_.gravityEnabled = false;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = true;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, -0.5);

    Eigen::VectorXd q(2), qprev(2), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, -0.5;
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    // Exactly at floor: y == floorY, so condition (y < floorY) is false -> no force
    ASSERT_NEAR(F[1], 0.0, 1e-10);
}

void test_floorForce_preventsHeavyFallthrough()
{
    // Simulate a heavy particle (100 kg) with gravity and floor
    // It should NOT fall through y = -1.0
    params_.gravityEnabled = true;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = true;
    params_.particleMass = 100.0;
    params_.maxSpringDist = 0.01;
    params_.integrator = SimParameters::TI_IMPLICIT_EULER;
    params_.timeStep = 0.001;

    addParticle(0.0, 0.0);

    // Simulate 500 steps
    for (int step = 0; step < 500; step++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;
    }

    // Particle should not have fallen below y = -1.0
    ASSERT_GT(particles_[0].pos[1], -1.0);
}

// TEST 10: Explicit Euler integrator
void test_explicitEuler_freefall()
{
    // Single particle in gravity, no floor, explicit Euler
    params_.gravityEnabled = true;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_EXPLICIT_EULER;
    params_.particleMass = 1.0;
    params_.gravityG = -10.0;
    params_.timeStep = 0.01;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 0.0);

    // One step: q_new = q + h*v = (0,0), v_new = v + h*g = (0, -0.1)
    Eigen::VectorXd q, qprev, v;
    buildConfiguration(q, qprev, v);
    numericalIntegration(q, qprev, v);

    ASSERT_NEAR(q[0], 0.0, 1e-10);
    ASSERT_NEAR(q[1], 0.0, 1e-10);      // position hasn't changed yet (v was 0)
    ASSERT_NEAR(v[0], 0.0, 1e-10);
    ASSERT_NEAR(v[1], -0.1, 1e-10);     // velocity updated: 0 + 0.01 * (-10) = -0.1
}

void test_explicitEuler_multiStep()
{
    // Multi-step freefall check
    params_.gravityEnabled = true;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_EXPLICIT_EULER;
    params_.particleMass = 1.0;
    params_.gravityG = -10.0;
    params_.timeStep = 0.001;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 1.0);

    int nsteps = 100;
    for (int i = 0; i < nsteps; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;
    }

    // After 100 steps of 0.001 = 0.1 seconds
    // Approximate: y ≈ 1.0 - 0.5 * 10 * 0.1^2 = 1.0 - 0.05 = 0.95
    // Explicit Euler is first-order, so some error is expected
    ASSERT_NEAR(particles_[0].pos[1], 0.95, 0.01);
}

// TEST 11: Velocity Verlet integrator
void test_velocityVerlet_freefall()
{
    params_.gravityEnabled = true;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_VELOCITY_VERLET;
    params_.particleMass = 1.0;
    params_.gravityG = -10.0;
    params_.timeStep = 0.001;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 1.0);

    int nsteps = 100;
    for (int i = 0; i < nsteps; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;
    }

    // After 0.1 seconds: y ≈ 1.0 - 0.5 * 10 * 0.01 = 0.95
    ASSERT_NEAR(particles_[0].pos[1], 0.95, 0.01);
}

// TEST 12: Implicit Euler integrator
void test_implicitEuler_freefall()
{
    params_.gravityEnabled = true;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_IMPLICIT_EULER;
    params_.particleMass = 1.0;
    params_.gravityG = -10.0;
    params_.timeStep = 0.001;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 1.0);

    int nsteps = 100;
    for (int i = 0; i < nsteps; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;
    }

    ASSERT_NEAR(particles_[0].pos[1], 0.95, 0.01);
}

void test_implicitEuler_springOscillation()
{
    // Two particles connected by a spring, no gravity
    // The system should oscillate without exploding
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_IMPLICIT_EULER;
    params_.particleMass = 1.0;
    params_.springStiffness = 100.0;
    params_.maxSpringDist = 2.0;
    params_.timeStep = 0.01;

    // Fixed anchor at origin, free particle at (1, 0)
    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(1.0, 0.0);  // rest length = 1.0

    // Stretch the free particle
    particles_[1].pos[0] = 1.5;

    // Simulate 1000 steps and verify no explosion
    for (int i = 0; i < 1000; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;

        // Should stay bounded
        ASSERT_LT(std::abs(particles_[1].pos[0]), 5.0);
        ASSERT_LT(std::abs(particles_[1].pos[1]), 5.0);
    }
}

// TEST 13: Implicit Midpoint integrator
void test_implicitMidpoint_freefall()
{
    params_.gravityEnabled = true;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_IMPLICIT_MIDPOINT;
    params_.particleMass = 1.0;
    params_.gravityG = -10.0;
    params_.timeStep = 0.001;
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 1.0);

    int nsteps = 100;
    for (int i = 0; i < nsteps; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;
    }

    // Implicit midpoint is second-order so there's a small offset vs the
    // first-order approximation 0.95; use a wider tolerance.
    ASSERT_NEAR(particles_[0].pos[1], 0.95, 0.03);
}

void test_implicitMidpoint_energyConservation()
{
    // Implicit midpoint is symplectic: total energy should be approximately conserved
    // for a spring system without damping
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_IMPLICIT_MIDPOINT;
    params_.particleMass = 1.0;
    params_.springStiffness = 50.0;
    params_.maxSpringDist = 2.0;
    params_.timeStep = 0.005;
    params_.maxSpringStrain = 100.0;  // prevent snapping

    // Fixed anchor, free particle stretched
    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(1.0, 0.0);

    particles_[1].pos[0] = 1.3;  // stretch by 0.3

    Spring* s = dynamic_cast<Spring*>(connectors_[0]);
    double k = s->stiffness;
    double L = s->restlen;

    // Compute initial energy: E = 0.5*m*v^2 + 0.5*k*(r-L)^2
    auto computeEnergy = [&]() -> double {
        double r = (particles_[1].pos - particles_[0].pos).norm();
        double KE = 0.5 * particles_[1].mass * particles_[1].vel.squaredNorm();
        double PE = 0.5 * k * (r - L) * (r - L);
        return KE + PE;
    };

    double E0 = computeEnergy();

    for (int i = 0; i < 500; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;
    }

    double Ef = computeEnergy();

    // Energy should be conserved within ~20% for implicit midpoint (some numerical
    // dissipation is expected with finite step size and Newton tolerance)
    ASSERT_NEAR(E0, Ef, 0.2 * E0);
}

// TEST 14: Fixed particles
void test_fixedParticles_dontMove()
{
    params_.gravityEnabled = true;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_EXPLICIT_EULER;
    params_.timeStep = 0.01;
    params_.maxSpringDist = 0.01;

    params_.particleFixed = true;
    addParticle(0.5, 0.5);

    Eigen::Vector2d initialPos = particles_[0].pos;

    for (int i = 0; i < 100; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;
    }

    ASSERT_NEAR(particles_[0].pos[0], initialPos[0], 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], initialPos[1], 1e-10);
    ASSERT_NEAR(particles_[0].vel[0], 0.0, 1e-10);
    ASSERT_NEAR(particles_[0].vel[1], 0.0, 1e-10);
}

void test_fixedParticles_withSpring()
{
    // Fixed particle connected to free particle via spring
    // Only the free particle should move
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_EXPLICIT_EULER;
    params_.timeStep = 0.001;
    params_.springStiffness = 100.0;
    params_.maxSpringDist = 2.0;

    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(1.0, 0.0);

    // Stretch free particle
    particles_[1].pos[0] = 1.5;

    for (int i = 0; i < 10; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;
    }

    // Fixed particle hasn't moved
    ASSERT_NEAR(particles_[0].pos[0], 0.0, 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], 0.0, 1e-10);
    // Free particle has moved (pulled back toward rest length)
    ASSERT_LT(particles_[1].pos[0], 1.5);
}

// TEST 15: Spring snapping
void test_springSnapping_noSnap()
{
    params_.maxSpringStrain = 0.5;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);  // rest length = 0.5

    // Move to stretch by 20% (below 50% threshold)
    particles_[1].pos[0] = 0.6;

    ASSERT_EQ((int)connectors_.size(), 1);
    pruneOverstrainedSprings();
    ASSERT_EQ((int)connectors_.size(), 1);  // spring survives
}

void test_springSnapping_snap()
{
    params_.maxSpringStrain = 0.5;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);  // rest length = 0.5

    // Move to stretch by 60% (above 50% threshold)
    particles_[1].pos[0] = 0.8;  // strain = (0.8-0.5)/0.5 = 0.6

    ASSERT_EQ((int)connectors_.size(), 1);
    pruneOverstrainedSprings();
    ASSERT_EQ((int)connectors_.size(), 0);  // spring snapped
}

void test_springSnapping_exactThreshold()
{
    params_.maxSpringStrain = 0.5;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    // Exactly at threshold: strain = 0.5 (should NOT snap, needs to exceed)
    particles_[1].pos[0] = 0.75;  // strain = (0.75-0.5)/0.5 = 0.5

    pruneOverstrainedSprings();
    ASSERT_EQ((int)connectors_.size(), 1);
}

void test_springSnapping_selective()
{
    // Multiple springs, only some should snap
    params_.maxSpringStrain = 0.5;
    params_.maxSpringDist = 2.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    addParticle(0.2, 0.0);

    int initialSprings = connectors_.size();
    ASSERT_GT(initialSprings, 0);

    // Stretch particle 1 far away (should snap springs to it)
    particles_[1].pos[0] = 5.0;

    pruneOverstrainedSprings();

    // Some springs should have snapped
    ASSERT_LT((int)connectors_.size(), initialSprings);
}

// TEST 16: Saw collision
void test_sawCollision_particleInSaw()
{
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 0.0);
    addParticle(1.0, 0.0);
    saws_.push_back(Saw(Eigen::Vector2d(0.0, 0.0), 0.1));  // saw right on particle 0

    ASSERT_EQ((int)particles_.size(), 2);
    deleteSawedObjects();
    ASSERT_EQ((int)particles_.size(), 1);  // particle 0 deleted
    ASSERT_NEAR(particles_[0].pos[0], 1.0, 1e-10);  // particle 1 survives
}

void test_sawCollision_particleOutsideSaw()
{
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 0.0);
    saws_.push_back(Saw(Eigen::Vector2d(1.0, 0.0), 0.1));  // far from particle

    deleteSawedObjects();
    ASSERT_EQ((int)particles_.size(), 1);  // particle survives
}

void test_sawCollision_springCut()
{
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(1.0, 0.0);  // spring between them

    ASSERT_EQ((int)connectors_.size(), 1);

    // Saw at midpoint of spring, big enough to touch it
    saws_.push_back(Saw(Eigen::Vector2d(0.5, 0.0), 0.1));

    deleteSawedObjects();
    ASSERT_EQ((int)connectors_.size(), 0);  // spring cut
    ASSERT_EQ((int)particles_.size(), 2);   // particles survive
}

void test_sawCollision_springNotCut()
{
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(1.0, 0.0);

    // Saw far from the spring
    saws_.push_back(Saw(Eigen::Vector2d(0.5, 2.0), 0.1));

    deleteSawedObjects();
    ASSERT_EQ((int)connectors_.size(), 1);  // spring survives
}

void test_sawCollision_particleDeleteCascadesSpring()
{
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 100.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    ASSERT_EQ((int)connectors_.size(), 1);

    // Kill particle 0 with saw -> spring should also die
    saws_.push_back(Saw(Eigen::Vector2d(0.0, 0.0), 0.1));

    deleteSawedObjects();
    ASSERT_EQ((int)particles_.size(), 1);
    ASSERT_EQ((int)connectors_.size(), 0);
}

void test_sawCollision_indexRemapping()
{
    // Delete a particle in the middle and verify spring indices are remapped
    params_.maxSpringDist = 0.01;  // no springs at first

    addParticle(0.0, 0.0);   // particle 0
    addParticle(0.0, 0.5);   // particle 1 (will be killed by saw)
    addParticle(1.0, 0.0);   // particle 2

    // Manually add a spring between particle 0 and particle 2
    // (along x-axis, far from the saw at (0, 0.5))
    connectors_.push_back(new Spring(0, 2, 0.0, 100.0, 1.0, true));

    // Kill particle 1 with a saw
    saws_.push_back(Saw(Eigen::Vector2d(0.0, 0.5), 0.1));

    deleteSawedObjects();

    ASSERT_EQ((int)particles_.size(), 2);
    ASSERT_EQ((int)connectors_.size(), 1);
    // After remapping: old particle 0 -> new 0, old particle 2 -> new 1
    ASSERT_EQ(connectors_[0]->p1, 0);
    ASSERT_EQ(connectors_[0]->p2, 1);
}

void test_sawCollision_multipleSaws()
{
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 0.0);
    addParticle(1.0, 0.0);
    addParticle(2.0, 0.0);

    saws_.push_back(Saw(Eigen::Vector2d(0.0, 0.0), 0.1));
    saws_.push_back(Saw(Eigen::Vector2d(2.0, 0.0), 0.1));

    deleteSawedObjects();

    ASSERT_EQ((int)particles_.size(), 1);
    ASSERT_NEAR(particles_[0].pos[0], 1.0, 1e-10);
}

// TEST 17: Out-of-bounds particle deletion
void test_outOfBounds_deletion()
{
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 0.0);
    addParticle(6.0, 0.0);  // out of bounds (> 5.0)

    deleteSawedObjects();

    ASSERT_EQ((int)particles_.size(), 1);
    ASSERT_NEAR(particles_[0].pos[0], 0.0, 1e-10);
}

void test_outOfBounds_yAxis()
{
    params_.maxSpringDist = 0.01;

    addParticle(0.0, 0.0);
    addParticle(0.0, -6.0);  // out of bounds

    deleteSawedObjects();

    ASSERT_EQ((int)particles_.size(), 1);
}

// TEST 18: Reset simulation
void test_resetSimulation()
{
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    saws_.push_back(Saw(Eigen::Vector2d(2.0, 0.0), 0.1));

    ASSERT_GT((int)particles_.size(), 0);
    ASSERT_GT((int)connectors_.size(), 0);
    ASSERT_GT((int)saws_.size(), 0);

    initSimulation();

    ASSERT_EQ((int)particles_.size(), 0);
    ASSERT_EQ((int)connectors_.size(), 0);
    ASSERT_EQ((int)saws_.size(), 0);
    ASSERT_NEAR(time_, 0.0, 1e-10);
}

// TEST 19: Combined forces
void test_combinedForces_gravityAndSpring()
{
    // Verify gravity + spring forces add correctly
    params_.gravityEnabled = true;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.maxSpringDist = 2.0;
    params_.springStiffness = 100.0;
    params_.particleMass = 1.0;
    params_.gravityG = -10.0;

    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(1.0, 0.0);

    Eigen::VectorXd q(4), qprev(4), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0, 1.5, 0.0;  // stretched
    qprev = q;

    computeForceAndHessian(q, qprev, F, H);

    // Free particle (index 1) should feel:
    // Spring: k*(r-L)*(-dhat) = (100/1)*(1.5-1.0)*(-1,0) = (-50, 0)
    // Gravity: (0, -10)
    // Total: (-50, -10)
    ASSERT_NEAR(F[2], -50.0, 1e-6);
    ASSERT_NEAR(F[3], -10.0, 1e-6);
}

// TEST 20: Integrator stability comparison
void test_explicitEuler_instability()
{
    // Explicit Euler with large time step should show energy growth (instability)
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_EXPLICIT_EULER;
    params_.particleMass = 1.0;
    params_.springStiffness = 1000.0;  // stiff spring
    params_.maxSpringDist = 2.0;
    params_.timeStep = 0.01;  // large for stiff spring
    params_.maxSpringStrain = 1000.0;  // prevent snapping

    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(1.0, 0.0);
    particles_[1].pos[0] = 1.3;

    double initialPos = particles_[1].pos[0];

    for (int i = 0; i < 100; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;
    }

    // For a stiff spring with explicit Euler and large dt, amplitude should grow
    // (This is the expected instability of explicit Euler)
    double finalDist = std::abs(particles_[1].pos[0]);
    ASSERT_GT(finalDist, 1.3);  // energy has grown
}

void test_implicitEuler_stability()
{
    // Same stiff spring but with implicit Euler: should remain stable
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_IMPLICIT_EULER;
    params_.particleMass = 1.0;
    params_.springStiffness = 1000.0;
    params_.maxSpringDist = 2.0;
    params_.timeStep = 0.01;
    params_.maxSpringStrain = 1000.0;

    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(1.0, 0.0);
    particles_[1].pos[0] = 1.3;

    for (int i = 0; i < 100; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;
    }

    // Implicit Euler should stay bounded
    ASSERT_LT(std::abs(particles_[1].pos[0]), 5.0);
    ASSERT_LT(std::abs(particles_[1].pos[1]), 5.0);
}

// TEST 21: Empty simulation edge cases
void test_emptySimulation_noParticles()
{
    // Calling simulateOneStep with no particles should not crash
    Eigen::VectorXd q, qprev, v;
    buildConfiguration(q, qprev, v);
    ASSERT_EQ(q.size(), 0);
}

void test_emptySimulation_sawNoParticles()
{
    // Saws with no particles should not crash
    saws_.push_back(Saw(Eigen::Vector2d(0.0, 0.0), 0.5));
    deleteSawedObjects();
    ASSERT_EQ((int)particles_.size(), 0);
}

// TEST 22: Damping effect
void test_dampingReducesOscillation()
{
    // With damping, oscillation amplitude should decrease over time
    params_.gravityEnabled = false;
    params_.springsEnabled = true;
    params_.dampingEnabled = true;
    params_.floorEnabled = false;
    params_.integrator = SimParameters::TI_IMPLICIT_EULER;
    params_.particleMass = 1.0;
    params_.springStiffness = 100.0;
    params_.dampingStiffness = 5.0;
    params_.maxSpringDist = 2.0;
    params_.timeStep = 0.001;
    params_.maxSpringStrain = 100.0;

    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(1.0, 0.0);
    particles_[1].pos[0] = 1.5;  // stretch by 0.5

    double maxDeviation_early = 0;
    double maxDeviation_late = 0;

    for (int i = 0; i < 2000; i++)
    {
        Eigen::VectorXd q, qprev, v;
        buildConfiguration(q, qprev, v);
        numericalIntegration(q, qprev, v);
        unbuildConfiguration(q, v);
        time_ += params_.timeStep;

        double deviation = std::abs(particles_[1].pos[0] - 1.0);
        if (i < 500)
            maxDeviation_early = std::max(maxDeviation_early, deviation);
        if (i >= 1500)
            maxDeviation_late = std::max(maxDeviation_late, deviation);
    }

    // Late oscillations should be smaller than early ones
    ASSERT_LT(maxDeviation_late, maxDeviation_early);
}

// MAIN: Run all tests
int main(int argc, char** argv)
{
    std::cout << " Universe of Goo — Milestone I Test Suite" << std::endl;
    std::cout << std::endl;

    //Geometry Helpers
    std::cout << "Geometry Helpers" << std::endl;
    runTest("PointToLine_Perpendicular",      test_distancePointToLine_perpendicular);
    runTest("PointToLine_OnLine",             test_distancePointToLine_onLine);
    runTest("PointToLine_Diagonal",           test_distancePointToLine_diagonal);
    runTest("PointToLine_NegativeProjection",  test_distancePointToLine_negativeProjection);
    runTest("PointToLine_Degenerate",          test_distancePointToLine_degenerateSegment);
    runTest("PointToSegment_Perpendicular",    test_distancePointToSegment_perpendicular);
    runTest("PointToSegment_ClampToQ1",        test_distancePointToSegment_clampToQ1);
    runTest("PointToSegment_ClampToQ2",        test_distancePointToSegment_clampToQ2);
    runTest("PointToSegment_Endpoint",         test_distancePointToSegment_endpoint);
    runTest("PointToSegment_Midpoint",         test_distancePointToSegment_midpoint);

    //Particle/Spring Creation
    std::cout << std::endl << "Particle/Spring Creation" << std::endl;
    runTest("AddParticle_Single",              test_addParticle_single);
    runTest("AddParticle_TwoClose",            test_addParticle_twoClose);
    runTest("AddParticle_TwoFar",              test_addParticle_twoFar);
    runTest("AddParticle_SpringStiffness",     test_addParticle_springStiffness);
    runTest("AddParticle_MultipleConnections", test_addParticle_multipleConnections);
    runTest("AddParticle_FixedParticle",       test_addParticle_fixedParticle);
    runTest("AddParticle_MassPreserved",       test_addParticle_massPreserved);

    //Configuration Vectors
    std::cout << std::endl << "Configuration Vectors" << std::endl;
    runTest("BuildConfig_Basic",               test_buildConfiguration_basic);
    runTest("BuildConfig_PrevPos",             test_buildConfiguration_prevpos);
    runTest("UnbuildConfig_Basic",             test_unbuildConfiguration_basic);
    runTest("UnbuildConfig_SavesPrevPos",      test_unbuildConfiguration_savesPrevpos);
    runTest("BuildUnbuild_Roundtrip",          test_buildUnbuild_roundtrip);

    //Mass Matrix
    std::cout << std::endl << "Mass Matrix" << std::endl;
    runTest("MassMatrix_SingleParticle",       test_massMatrix_singleParticle);
    runTest("MassMatrix_FixedParticle",        test_massMatrix_fixedParticle);
    runTest("MassMatrix_MixedParticles",       test_massMatrix_mixedParticles);
    runTest("MassMatrix_IsDiagonal",           test_massMatrix_isDiagonal);

    //Gravity Force
    std::cout << std::endl << "Gravity Force" << std::endl;
    runTest("Gravity_SingleParticle",          test_gravity_singleParticle);
    runTest("Gravity_FixedParticle",           test_gravity_fixedParticle);
    runTest("Gravity_ZeroHessian",             test_gravity_zeroHessian);
    runTest("Gravity_Disabled",                test_gravity_disabled);

    //Spring Force
    std::cout << std::endl << "Spring Force" << std::endl;
    runTest("SpringForce_AtRest",              test_springForce_atRest);
    runTest("SpringForce_Stretched",           test_springForce_stretched);
    runTest("SpringForce_Compressed",          test_springForce_compressed);
    runTest("SpringForce_NewtonThirdLaw",      test_springForce_newtonThirdLaw);
    runTest("SpringForce_HessianSymmetric",    test_springForce_hessianSymmetric);
    runTest("SpringForce_FiniteDiffCheck",     test_springForce_finiteDifferenceCheck);
    runTest("SpringForce_HessianFiniteDiff",   test_springForce_hessianFiniteDifferenceCheck);

    //Damping Force
    std::cout << std::endl << "Damping Force" << std::endl;
    runTest("Damping_ZeroVelocity",            test_dampingForce_zeroVelocity);
    runTest("Damping_OpposesRelVelocity",      test_dampingForce_opposesRelativeVelocity);
    runTest("Damping_NewtonThirdLaw",          test_dampingForce_newtonThirdLaw);

    //Floor Force
    std::cout << std::endl << "Floor Force" << std::endl;
    runTest("Floor_AboveFloor",                test_floorForce_aboveFloor);
    runTest("Floor_BelowFloor",                test_floorForce_belowFloor);
    runTest("Floor_AtExactFloor",              test_floorForce_atExactFloor);
    runTest("Floor_PreventsHeavyFallthrough",  test_floorForce_preventsHeavyFallthrough);

    //Integrators
    std::cout << std::endl << "Explicit Euler" << std::endl;
    runTest("ExplicitEuler_Freefall",          test_explicitEuler_freefall);
    runTest("ExplicitEuler_MultiStep",         test_explicitEuler_multiStep);

    std::cout << std::endl << "Velocity Verlet" << std::endl;
    runTest("VelocityVerlet_Freefall",         test_velocityVerlet_freefall);

    std::cout << std::endl << "Implicit Euler" << std::endl;
    runTest("ImplicitEuler_Freefall",          test_implicitEuler_freefall);
    runTest("ImplicitEuler_SpringOscillation", test_implicitEuler_springOscillation);

    std::cout << std::endl << "Implicit Midpoint" << std::endl;
    runTest("ImplicitMidpoint_Freefall",       test_implicitMidpoint_freefall);
    runTest("ImplicitMidpoint_EnergyConserv",  test_implicitMidpoint_energyConservation);

    //Fixed Particles
    std::cout << std::endl << "Fixed Particles" << std::endl;
    runTest("FixedParticles_DontMove",         test_fixedParticles_dontMove);
    runTest("FixedParticles_WithSpring",       test_fixedParticles_withSpring);

    //Spring Snapping
    std::cout << std::endl << "Spring Snapping" << std::endl;
    runTest("SpringSnap_NoSnap",               test_springSnapping_noSnap);
    runTest("SpringSnap_Snap",                 test_springSnapping_snap);
    runTest("SpringSnap_ExactThreshold",       test_springSnapping_exactThreshold);
    runTest("SpringSnap_Selective",            test_springSnapping_selective);

    //Saw Collision
    std::cout << std::endl << "Saw Collision" << std::endl;
    runTest("Saw_ParticleInSaw",               test_sawCollision_particleInSaw);
    runTest("Saw_ParticleOutsideSaw",          test_sawCollision_particleOutsideSaw);
    runTest("Saw_SpringCut",                   test_sawCollision_springCut);
    runTest("Saw_SpringNotCut",                test_sawCollision_springNotCut);
    runTest("Saw_ParticleDeleteCascade",       test_sawCollision_particleDeleteCascadesSpring);
    runTest("Saw_IndexRemapping",              test_sawCollision_indexRemapping);
    runTest("Saw_MultipleSaws",                test_sawCollision_multipleSaws);

    //Out of Bounds
    std::cout << std::endl << "Out of Bounds" << std::endl;
    runTest("OutOfBounds_Deletion",            test_outOfBounds_deletion);
    runTest("OutOfBounds_YAxis",               test_outOfBounds_yAxis);

    //Reset
    std::cout << std::endl << "Reset Simulation" << std::endl;
    runTest("ResetSimulation",                 test_resetSimulation);

    //Combined Forces
    std::cout << std::endl << "Combined Forces" << std::endl;
    runTest("Combined_GravityAndSpring",       test_combinedForces_gravityAndSpring);

    //Stability
    std::cout << std::endl << "Integrator Stability" << std::endl;
    runTest("ExplicitEuler_Instability",       test_explicitEuler_instability);
    runTest("ImplicitEuler_Stability",         test_implicitEuler_stability);

    //Edge Cases
    std::cout << std::endl << "Edge Cases" << std::endl;
    runTest("Empty_NoParticles",               test_emptySimulation_noParticles);
    runTest("Empty_SawNoParticles",            test_emptySimulation_sawNoParticles);

    //Damping Behavior
    std::cout << std::endl << "Damping Behavior" << std::endl;
    runTest("Damping_ReducesOscillation",      test_dampingReducesOscillation);

    // Summary
    std::cout << std::endl;
    std::cout << " RESULTS: " << g_tests_passed << " passed, "
              << g_tests_failed << " failed  ("
              << g_assertions_passed << " assertions passed, "
              << g_assertions_failed << " assertions failed)" << std::endl;

    if (!g_failures.empty())
    {
        std::cout << std::endl << "FAILURES:" << std::endl;
        for (size_t i = 0; i < g_failures.size(); i++)
            std::cout << "  " << g_failures[i] << std::endl;
    }

    return g_tests_failed > 0 ? 1 : 0;
}
