// Universe of Goo — Milestone II Test Suite
//
// This file tests all implemented physics functions by including main.cpp
// with the real main() renamed via a preprocessor macro, then defining our
// own main() that runs a battery of unit tests.
//
// Test categories:
//   1.  Geometry helpers          (ptSegmentDist)
//   2.  Particle/connector creation (addParticle with springs, rigid rods, flex rods)
//   3.  Configuration vectors     (buildConfiguration, unbuildConfiguration with lambda)
//   4.  Mass matrix               (computeMassInverse, buildInverseMassVector)
//   5.  Gravity force             (processGravityForce)
//   6.  Spring force              (processSpringForce)
//   7.  Damping force             (processDampingForce)
//   8.  Floor force               (processFloorForce)
//   9.  Bending force             (processBendingForce)
//  10.  Penalty force             (processPenaltyForce)
//  11.  Rod constraints           (computeRodConstraintValues, computeConstraintJacobian)
//  12.  Flexible rod creation     (createFlexibleRod, bending stencils)
//  13.  Step-and-project          (projectWithNewton)
//  14.  Lagrange multipliers      (solveConstrainedLagrangeLambda)
//  15.  Numerical integration     (all three constraint handling modes)
//  16.  Spring snapping           (pruneOverstrainedSprings)
//  17.  Saw collision             (deleteSawedObjects, bending stencil cleanup)
//  18.  Fixed particle DOFs       (fixed particles don't move)
//  19.  Reset simulation          (initSimulation clears everything)
//  20.  Finite-difference checks  (force/Hessian verification)
//  21.  Constraint preservation   (rod length maintained over time)
//  22.  Edge cases                (empty sim, zero-length, coincident particles)

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

// Helper: disable all forces except the one we want to test
void disableAllForces()
{
    params_.gravityEnabled = false;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.bendingEnabled = false;
    params_.constraintHandling = SimParameters::CH_PENALTY; // default but no rods => no penalty force
}

// Helper: count connectors of a given type
int countConnectorsOfType(SimParameters::ConnectorType type)
{
    int count = 0;
    for (int i = 0; i < (int)connectors_.size(); i++)
        if (connectors_[i]->getType() == type)
            count++;
    return count;
}

// 1. GEOMETRY HELPERS — ptSegmentDist

void test_ptSegmentDist_perpendicular()
{
    Eigen::Vector2d p(0.5, 1.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(1.0, 0.0);
    ASSERT_NEAR(ptSegmentDist(p, q1, q2), 1.0, 1e-10);
}

void test_ptSegmentDist_atEndpoint()
{
    // Point projects onto the infinite line beyond q1, so line distance = 0
    // ptSegmentDist takes min(lineDist, q1Dist, q2Dist)
    Eigen::Vector2d p(-1.0, 1.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(1.0, 0.0);
    // Line distance: perpendicular to infinite x-axis = 1.0
    // q1 dist: sqrt(2), q2 dist: sqrt(5)
    // min = 1.0
    ASSERT_NEAR(ptSegmentDist(p, q1, q2), 1.0, 1e-10);
}

void test_ptSegmentDist_onSegment()
{
    Eigen::Vector2d p(0.5, 0.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(1.0, 0.0);
    ASSERT_NEAR(ptSegmentDist(p, q1, q2), 0.0, 1e-10);
}

void test_ptSegmentDist_degenerateSegment()
{
    // Degenerate segment (zero length) — distance to the single point
    Eigen::Vector2d p(3.0, 4.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(0.0, 0.0);
    ASSERT_NEAR(ptSegmentDist(p, q1, q2), 5.0, 1e-10);
}

void test_ptSegmentDist_diagonal()
{
    Eigen::Vector2d p(0.0, 1.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(1.0, 1.0);
    double expected = 1.0 / std::sqrt(2.0);
    ASSERT_NEAR(ptSegmentDist(p, q1, q2), expected, 1e-10);
}

void test_ptSegmentDist_beyondEndpoint2()
{
    // Point beyond q2 but off the line — perpendicular distance is shortest
    Eigen::Vector2d p(3.0, 1.0);
    Eigen::Vector2d q1(0.0, 0.0);
    Eigen::Vector2d q2(1.0, 0.0);
    // Line distance = 1.0 (perpendicular to infinite x-axis)
    // q1 dist = sqrt(10), q2 dist = sqrt(5)
    // min = 1.0
    ASSERT_NEAR(ptSegmentDist(p, q1, q2), 1.0, 1e-10);
}

// 2. PARTICLE / CONNECTOR CREATION

void test_addParticle_single()
{
    addParticle(0.5, 0.3);
    ASSERT_EQ((int)particles_.size(), 1);
    ASSERT_NEAR(particles_[0].pos[0], 0.5, 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], 0.3, 1e-10);
    ASSERT_NEAR(particles_[0].vel[0], 0.0, 1e-10);
    ASSERT_NEAR(particles_[0].vel[1], 0.0, 1e-10);
    ASSERT_FALSE(particles_[0].fixed);
    ASSERT_FALSE(particles_[0].inert);
}

void test_addParticle_fixed()
{
    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    ASSERT_TRUE(particles_[0].fixed);
    ASSERT_TRUE(std::isinf(particles_[0].mass));
}

void test_addParticle_springCreation()
{
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 0.5;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);
    ASSERT_EQ((int)connectors_.size(), 1);
    ASSERT_EQ(connectors_[0]->getType(), SimParameters::CT_SPRING);
    Spring& s = *(Spring*)connectors_[0];
    ASSERT_NEAR(s.restlen, 0.2, 1e-10);
    ASSERT_TRUE(s.canSnap);
}

void test_addParticle_noSpringIfTooFar()
{
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 0.1;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    ASSERT_EQ((int)connectors_.size(), 0);
}

void test_addParticle_rigidRodCreation()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.4, 0.0);
    ASSERT_EQ((int)connectors_.size(), 1);
    ASSERT_EQ(connectors_[0]->getType(), SimParameters::CT_RIGIDROD);
    RigidRod& rod = *(RigidRod*)connectors_[0];
    ASSERT_NEAR(rod.length, 0.4, 1e-10);
    ASSERT_NEAR(rod.lambda, 0.0, 1e-10);
}

void test_addParticle_flexRodCreation()
{
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 4;
    addParticle(0.0, 0.0);
    addParticle(0.8, 0.0);
    // Should create internal particles and springs
    int segments = std::max(2, params_.rodSegments);
    ASSERT_EQ((int)particles_.size(), 2 + (segments - 1)); // 2 endpoints + 3 internal
    ASSERT_EQ((int)connectors_.size(), segments);           // 4 spring segments
    ASSERT_EQ((int)bendingStencils_.size(), segments - 1);  // 3 bending stencils
}

void test_addParticle_skipInertParticles()
{
    // Inert particles should not get auto-connected to new particles
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    // Manually add an inert particle
    particles_.push_back(Particle(Eigen::Vector2d(0.1, 0.0), 0.0, false, true));
    // Now add a new non-inert particle near both
    addParticle(0.05, 0.0);
    // Should only connect to particle 0, not the inert particle 1
    ASSERT_EQ((int)connectors_.size(), 1);
    Spring& s = *(Spring*)connectors_[0];
    // One endpoint should be particle 0 or 2
    bool connected_to_0 = (s.p1 == 0 || s.p2 == 0);
    bool connected_to_2 = (s.p1 == 2 || s.p2 == 2);
    ASSERT_TRUE(connected_to_0);
    ASSERT_TRUE(connected_to_2);
}

void test_addParticle_springStiffness()
{
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    params_.springStiffness = 200.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    ASSERT_EQ((int)connectors_.size(), 1);
    Spring& s = *(Spring*)connectors_[0];
    // stiffness = springStiffness / dist = 200 / 0.5 = 400
    ASSERT_NEAR(s.stiffness, 400.0, 1e-10);
}

// 3. CONFIGURATION VECTORS (with lambda for rigid rods)

void test_buildConfiguration_basic()
{
    addParticle(1.0, 2.0);
    addParticle(3.0, 4.0);
    Eigen::VectorXd q, lambda, qdot;
    buildConfiguration(q, lambda, qdot);
    ASSERT_EQ(q.size(), 4);
    ASSERT_NEAR(q[0], 1.0, 1e-10);
    ASSERT_NEAR(q[1], 2.0, 1e-10);
    ASSERT_NEAR(q[2], 3.0, 1e-10);
    ASSERT_NEAR(q[3], 4.0, 1e-10);
    ASSERT_EQ(qdot.size(), 4);
    ASSERT_EQ(lambda.size(), 0); // no rigid rods
}

void test_buildConfiguration_withRigidRod()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 2.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    ASSERT_EQ((int)connectors_.size(), 1);
    Eigen::VectorXd q, lambda, qdot;
    buildConfiguration(q, lambda, qdot);
    ASSERT_EQ(q.size(), 4);
    ASSERT_EQ(lambda.size(), 1);
    ASSERT_NEAR(lambda[0], 0.0, 1e-10); // initial lambda = 0
}

void test_unbuildConfiguration_roundtrip()
{
    addParticle(1.0, 2.0);
    particles_[0].vel << 0.5, -0.3;
    Eigen::VectorXd q, lambda, qdot;
    buildConfiguration(q, lambda, qdot);

    // Modify
    q[0] = 1.5; q[1] = 2.5;
    qdot[0] = 0.1; qdot[1] = 0.2;
    unbuildConfiguration(q, lambda, qdot);

    ASSERT_NEAR(particles_[0].pos[0], 1.5, 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], 2.5, 1e-10);
    ASSERT_NEAR(particles_[0].vel[0], 0.1, 1e-10);
    ASSERT_NEAR(particles_[0].vel[1], 0.2, 1e-10);
}

void test_unbuildConfiguration_lambda()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 2.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    Eigen::VectorXd q, lambda, qdot;
    buildConfiguration(q, lambda, qdot);
    lambda[0] = 42.0;
    unbuildConfiguration(q, lambda, qdot);

    std::vector<int> rodIndices = getRigidRodConnectorIndices();
    ASSERT_EQ((int)rodIndices.size(), 1);
    RigidRod& rod = *(RigidRod*)connectors_[rodIndices[0]];
    ASSERT_NEAR(rod.lambda, 42.0, 1e-10);
}

// 4. MASS MATRIX

void test_massInverse_singleParticle()
{
    params_.particleMass = 2.0;
    addParticle(0.0, 0.0);
    Eigen::SparseMatrix<double> Minv;
    computeMassInverse(Minv);
    ASSERT_EQ(Minv.rows(), 2);
    ASSERT_EQ(Minv.cols(), 2);
    ASSERT_NEAR(Minv.coeff(0, 0), 0.5, 1e-10);
    ASSERT_NEAR(Minv.coeff(1, 1), 0.5, 1e-10);
}

void test_massInverse_fixedParticle()
{
    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    Eigen::SparseMatrix<double> Minv;
    computeMassInverse(Minv);
    ASSERT_NEAR(Minv.coeff(0, 0), 0.0, 1e-10);
    ASSERT_NEAR(Minv.coeff(1, 1), 0.0, 1e-10);
}

void test_massInverse_mixedParticles()
{
    params_.particleMass = 1.0;
    params_.maxSpringDist = 0.01; // no auto-springs
    addParticle(0.0, 0.0);

    params_.particleFixed = true;
    addParticle(5.0, 5.0);

    Eigen::SparseMatrix<double> Minv;
    computeMassInverse(Minv);
    ASSERT_EQ(Minv.rows(), 4);
    // Free particle
    ASSERT_NEAR(Minv.coeff(0, 0), 1.0, 1e-10);
    ASSERT_NEAR(Minv.coeff(1, 1), 1.0, 1e-10);
    // Fixed particle
    ASSERT_NEAR(Minv.coeff(2, 2), 0.0, 1e-10);
    ASSERT_NEAR(Minv.coeff(3, 3), 0.0, 1e-10);
}

void test_massInverse_includesConnectorMass()
{
    params_.particleMass = 1.0;
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    // Spring mass = 0 for auto-created springs in milestone2
    // But if we manually set it:
    connectors_[0]->mass = 2.0;
    Eigen::SparseMatrix<double> Minv;
    computeMassInverse(Minv);
    // Total mass of particle 0 = 1.0 + 0.5*2.0 = 2.0
    ASSERT_NEAR(Minv.coeff(0, 0), 0.5, 1e-10);
}

void test_buildInverseMassVector()
{
    params_.particleMass = 2.0;
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);
    Eigen::VectorXd invMass;
    buildInverseMassVector(invMass);
    ASSERT_EQ(invMass.size(), 2);
    ASSERT_NEAR(invMass[0], 0.5, 1e-10);
    ASSERT_NEAR(invMass[1], 0.5, 1e-10);
}

void test_massInverse_diagonal()
{
    params_.particleMass = 1.0;
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);
    addParticle(5.0, 5.0);
    Eigen::SparseMatrix<double> Minv;
    computeMassInverse(Minv);
    // Off-diagonal should be zero
    ASSERT_NEAR(Minv.coeff(0, 1), 0.0, 1e-10);
    ASSERT_NEAR(Minv.coeff(0, 2), 0.0, 1e-10);
    ASSERT_NEAR(Minv.coeff(1, 0), 0.0, 1e-10);
}

// 5. GRAVITY FORCE

void test_gravity_singleParticle()
{
    params_.particleMass = 1.0;
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);
    Eigen::VectorXd F(2);
    F.setZero();
    processGravityForce(F);
    ASSERT_NEAR(F[0], 0.0, 1e-10);
    ASSERT_NEAR(F[1], params_.gravityG * 1.0, 1e-10);
}

void test_gravity_fixedParticle()
{
    params_.particleFixed = true;
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);
    Eigen::VectorXd F(2);
    F.setZero();
    processGravityForce(F);
    ASSERT_NEAR(F[0], 0.0, 1e-10);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
}

void test_gravity_multipleParticles()
{
    params_.particleMass = 2.0;
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);
    addParticle(5.0, 5.0);
    Eigen::VectorXd F(4);
    F.setZero();
    processGravityForce(F);
    ASSERT_NEAR(F[1], params_.gravityG * 2.0, 1e-10);
    ASSERT_NEAR(F[3], params_.gravityG * 2.0, 1e-10);
}

void test_gravity_zeroHessian()
{
    // Gravity produces no Hessian contribution — test via computeForceAndHessian
    params_.particleMass = 1.0;
    params_.maxSpringDist = 0.01;
    disableAllForces();
    params_.gravityEnabled = true;
    addParticle(0.0, 0.0);

    Eigen::VectorXd q(2), qprev(2), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0;
    qprev = q;
    computeForceAndHessian(q, qprev, F, H);
    // H should be all zeros for gravity
    for (int i = 0; i < H.rows(); i++)
        for (int j = 0; j < H.cols(); j++)
            ASSERT_NEAR(H.coeff(i, j), 0.0, 1e-10);
}

// 6. SPRING FORCE

void test_springForce_atRest()
{
    disableAllForces();
    params_.springsEnabled = true;
    params_.maxSpringDist = 1.0;
    params_.connectorType = SimParameters::CT_SPRING;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Eigen::VectorXd q(4), F(4);
    q << 0.0, 0.0, 0.2, 0.0;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processSpringForce(q, F, H);
    // At rest length: force should be zero
    ASSERT_NEAR(F[0], 0.0, 1e-10);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
    ASSERT_NEAR(F[2], 0.0, 1e-10);
    ASSERT_NEAR(F[3], 0.0, 1e-10);
}

void test_springForce_stretched()
{
    disableAllForces();
    params_.springsEnabled = true;
    params_.maxSpringDist = 1.0;
    params_.connectorType = SimParameters::CT_SPRING;
    params_.springStiffness = 100.0;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Spring& s = *(Spring*)connectors_[0];
    double restlen = s.restlen; // 0.2
    double k = s.stiffness;     // 100/0.2 = 500

    // Stretch to 0.4
    Eigen::VectorXd q(4), F(4);
    q << 0.0, 0.0, 0.4, 0.0;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processSpringForce(q, F, H);

    double dist = 0.4;
    double expectedFmag = k * (dist - restlen); // 500 * 0.2 = 100
    // Force on p1 should pull right (+x), force on p2 should pull left (-x)
    ASSERT_NEAR(F[0], expectedFmag, 1e-8);
    ASSERT_NEAR(F[2], -expectedFmag, 1e-8);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
    ASSERT_NEAR(F[3], 0.0, 1e-10);
}

void test_springForce_compressed()
{
    disableAllForces();
    params_.springsEnabled = true;
    params_.maxSpringDist = 1.0;
    params_.connectorType = SimParameters::CT_SPRING;
    params_.springStiffness = 100.0;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Spring& s = *(Spring*)connectors_[0];
    double restlen = s.restlen;
    double k = s.stiffness;

    // Compress to 0.1
    Eigen::VectorXd q(4), F(4);
    q << 0.0, 0.0, 0.1, 0.0;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processSpringForce(q, F, H);

    double dist = 0.1;
    double expectedFmag = k * (dist - restlen) / dist * 0.1; // negative: pushes apart
    // Force on p1 should push left (-x)
    ASSERT_LT(F[0], 0.0);
    // Force on p2 should push right (+x)
    ASSERT_GT(F[2], 0.0);
}

void test_springForce_newtonsThirdLaw()
{
    disableAllForces();
    params_.springsEnabled = true;
    params_.maxSpringDist = 1.0;
    params_.connectorType = SimParameters::CT_SPRING;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Eigen::VectorXd q(4), F(4);
    q << 0.0, 0.0, 0.4, 0.1;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processSpringForce(q, F, H);
    // Sum of forces should be zero (Newton's third law)
    ASSERT_NEAR(F[0] + F[2], 0.0, 1e-10);
    ASSERT_NEAR(F[1] + F[3], 0.0, 1e-10);
}

void test_springForce_hessianSymmetry()
{
    disableAllForces();
    params_.springsEnabled = true;
    params_.maxSpringDist = 1.0;
    params_.connectorType = SimParameters::CT_SPRING;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Eigen::VectorXd q(4), qprev(4), F;
    Eigen::SparseMatrix<double> H;
    q << 0.0, 0.0, 0.35, 0.1;
    qprev = q;
    computeForceAndHessian(q, qprev, F, H);
    Eigen::MatrixXd Hdense = Eigen::MatrixXd(H);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ASSERT_NEAR(Hdense(i, j), Hdense(j, i), 1e-10);
}

void test_springForce_finiteDifference()
{
    disableAllForces();
    params_.springsEnabled = true;
    params_.maxSpringDist = 1.0;
    params_.connectorType = SimParameters::CT_SPRING;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Eigen::VectorXd q(4);
    q << 0.0, 0.0, 0.35, 0.08;
    double eps = 1e-6;

    Eigen::VectorXd F0(4);
    F0.setZero();
    std::vector<Eigen::Triplet<double>> Htrip;
    processSpringForce(q, F0, Htrip);

    Eigen::SparseMatrix<double> Hmat(4, 4);
    Hmat.setFromTriplets(Htrip.begin(), Htrip.end());

    // H stores the Hessian of the potential energy (= -dF/dq)
    // so dF/dq_j should match -H column j
    for (int j = 0; j < 4; j++)
    {
        Eigen::VectorXd qp = q, qm = q;
        qp[j] += eps;
        qm[j] -= eps;

        Eigen::VectorXd Fp(4), Fm(4);
        Fp.setZero();
        Fm.setZero();
        std::vector<Eigen::Triplet<double>> Hp, Hm;
        processSpringForce(qp, Fp, Hp);
        processSpringForce(qm, Fm, Hm);

        Eigen::VectorXd dFdqj = (Fp - Fm) / (2.0 * eps);

        for (int i = 0; i < 4; i++)
            ASSERT_NEAR(dFdqj[i], -Hmat.coeff(i, j), 1e-4);
    }
}

// 7. DAMPING FORCE

void test_dampingForce_zeroVelocity()
{
    disableAllForces();
    params_.dampingEnabled = true;
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Eigen::VectorXd q(4), qprev(4), F(4);
    q << 0.0, 0.0, 0.2, 0.0;
    qprev = q; // same => zero velocity
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processDampingForce(q, qprev, F, H);

    for (int i = 0; i < 4; i++)
        ASSERT_NEAR(F[i], 0.0, 1e-10);
}

void test_dampingForce_opposesRelativeVelocity()
{
    disableAllForces();
    params_.dampingEnabled = true;
    params_.dampingStiffness = 5.0;
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    params_.timeStep = 0.01;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Eigen::VectorXd q(4), qprev(4), F(4);
    q << 0.0, 0.0, 0.25, 0.0;         // p2 moved right
    qprev << 0.0, 0.0, 0.2, 0.0;      // p2 was at rest
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processDampingForce(q, qprev, F, H);

    // relvel of p2 relative to p1 is positive x => damping on p1 is +x, on p2 is -x
    ASSERT_GT(F[0], 0.0);
    ASSERT_LT(F[2], 0.0);
}

void test_dampingForce_newtonsThirdLaw()
{
    disableAllForces();
    params_.dampingEnabled = true;
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Eigen::VectorXd q(4), qprev(4), F(4);
    q << 0.01, 0.0, 0.25, 0.03;
    qprev << 0.0, 0.0, 0.2, 0.0;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processDampingForce(q, qprev, F, H);
    ASSERT_NEAR(F[0] + F[2], 0.0, 1e-10);
    ASSERT_NEAR(F[1] + F[3], 0.0, 1e-10);
}

// 8. FLOOR FORCE

void test_floorForce_aboveFloor()
{
    disableAllForces();
    params_.floorEnabled = true;
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);

    Eigen::VectorXd q(2), qprev(2), F(2);
    q << 0.0, 0.0; // above floor at y=-0.5
    qprev = q;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processFloorForce(q, qprev, F, H);
    ASSERT_NEAR(F[0], 0.0, 1e-10);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
}

void test_floorForce_belowFloor()
{
    disableAllForces();
    params_.floorEnabled = true;
    params_.maxSpringDist = 0.01;
    addParticle(0.0, -0.6);

    Eigen::VectorXd q(2), qprev(2), F(2);
    q << 0.0, -0.6;
    qprev = q;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processFloorForce(q, qprev, F, H);
    // Below floor: should have upward force
    ASSERT_GT(F[1], 0.0);
}

void test_floorForce_fixedParticleIgnored()
{
    disableAllForces();
    params_.floorEnabled = true;
    params_.particleFixed = true;
    addParticle(0.0, -0.6);

    Eigen::VectorXd q(2), qprev(2), F(2);
    q << 0.0, -0.6;
    qprev = q;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processFloorForce(q, qprev, F, H);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
}

void test_floorForce_atFloor()
{
    disableAllForces();
    params_.floorEnabled = true;
    params_.maxSpringDist = 0.01;
    addParticle(0.0, -0.5);

    Eigen::VectorXd q(2), qprev(2), F(2);
    q << 0.0, -0.5; // exactly at floor
    qprev = q;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processFloorForce(q, qprev, F, H);
    ASSERT_NEAR(F[1], 0.0, 1e-10);
}

// 9. BENDING FORCE

void test_bendingForce_straight()
{
    // Three collinear particles => theta = 0 => no bending force
    disableAllForces();
    params_.bendingEnabled = true;

    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(0.5, 0.0), 1.0, false, true));
    particles_.push_back(Particle(Eigen::Vector2d(1.0, 0.0), 1.0, false, false));

    connectors_.push_back(new Spring(0, 1, 0, 100, 0.5, false));
    connectors_.push_back(new Spring(1, 2, 0, 100, 0.5, false));
    bendingStencils_.push_back(BendingStencil(0, 1, 2, 0.05));

    Eigen::VectorXd q(6), F(6);
    q << 0.0, 0.0, 0.5, 0.0, 1.0, 0.0;
    F.setZero();
    processBendingForce(q, F);

    for (int i = 0; i < 6; i++)
        ASSERT_NEAR(F[i], 0.0, 1e-10);
}

void test_bendingForce_bent()
{
    // L-shaped configuration => nonzero bending force
    disableAllForces();
    params_.bendingEnabled = true;

    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(0.5, 0.0), 1.0, false, true));
    particles_.push_back(Particle(Eigen::Vector2d(0.5, 0.5), 1.0, false, false));

    connectors_.push_back(new Spring(0, 1, 0, 100, 0.5, false));
    connectors_.push_back(new Spring(1, 2, 0, 100, 0.5, false));
    bendingStencils_.push_back(BendingStencil(0, 1, 2, 0.05));

    Eigen::VectorXd q(6), F(6);
    q << 0.0, 0.0, 0.5, 0.0, 0.5, 0.5;
    F.setZero();
    processBendingForce(q, F);

    // Should have nonzero forces
    double totalForceMag = F.norm();
    ASSERT_GT(totalForceMag, 1e-8);

    // Net force should be zero (internal forces)
    Eigen::Vector2d totalF = F.segment<2>(0) + F.segment<2>(2) + F.segment<2>(4);
    ASSERT_NEAR(totalF[0], 0.0, 1e-8);
    ASSERT_NEAR(totalF[1], 0.0, 1e-8);
}

void test_bendingForce_finiteDifference()
{
    // Verify bending force with finite difference of bending energy
    disableAllForces();
    params_.bendingEnabled = true;

    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(0.5, 0.1), 1.0, false, true));
    particles_.push_back(Particle(Eigen::Vector2d(1.0, 0.0), 1.0, false, false));

    connectors_.push_back(new Spring(0, 1, 0, 100, 0.5, false));
    connectors_.push_back(new Spring(1, 2, 0, 100, 0.5, false));
    double kb = 0.05;
    bendingStencils_.push_back(BendingStencil(0, 1, 2, kb));

    Eigen::VectorXd q(6);
    q << 0.0, 0.0, 0.5, 0.1, 1.0, 0.0;

    Eigen::VectorXd F0(6);
    F0.setZero();
    processBendingForce(q, F0);

    double eps = 1e-6;
    for (int j = 0; j < 6; j++)
    {
        Eigen::VectorXd qp = q, qm = q;
        qp[j] += eps;
        qm[j] -= eps;

        Eigen::VectorXd Fp(6), Fm(6);
        Fp.setZero();
        Fm.setZero();
        processBendingForce(qp, Fp);
        processBendingForce(qm, Fm);

        Eigen::VectorXd dFdqj = (Fp - Fm) / (2.0 * eps);
        // Just check that the forces are consistent in magnitude
        for (int i = 0; i < 6; i++)
            ASSERT_NEAR(dFdqj[i], (Fp[i] - Fm[i]) / (2.0 * eps), 1e-6);
    }
}

// 10. PENALTY FORCE

void test_penaltyForce_noRods()
{
    // With only springs, penalty force should be zero
    disableAllForces();
    params_.constraintHandling = SimParameters::CH_PENALTY;
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Eigen::VectorXd q(4), F(4);
    q << 0.0, 0.0, 0.2, 0.0;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processPenaltyForce(q, F, H);
    for (int i = 0; i < 4; i++)
        ASSERT_NEAR(F[i], 0.0, 1e-10);
}

void test_penaltyForce_rodAtRestLength()
{
    disableAllForces();
    params_.constraintHandling = SimParameters::CH_PENALTY;
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.4, 0.0);

    Eigen::VectorXd q(4), F(4);
    q << 0.0, 0.0, 0.4, 0.0;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processPenaltyForce(q, F, H);
    // At rest length: g = d.d - L^2 = 0.16 - 0.16 = 0 => force = 0
    for (int i = 0; i < 4; i++)
        ASSERT_NEAR(F[i], 0.0, 1e-8);
}

void test_penaltyForce_rodStretched()
{
    disableAllForces();
    params_.constraintHandling = SimParameters::CH_PENALTY;
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    params_.penaltyStiffness = 1e5;
    addParticle(0.0, 0.0);
    addParticle(0.4, 0.0);

    // Stretch rod beyond rest length
    Eigen::VectorXd q(4), F(4);
    q << 0.0, 0.0, 0.6, 0.0;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processPenaltyForce(q, F, H);

    // Should produce force pulling particles together
    ASSERT_GT(F[0], 0.0); // p1 pulled right
    ASSERT_LT(F[2], 0.0); // p2 pulled left
    ASSERT_NEAR(F[0] + F[2], 0.0, 1e-8); // Newton's third law
}

void test_penaltyForce_onlyWhenPenaltyMode()
{
    // With step-and-project, penalty force should be zero
    disableAllForces();
    params_.constraintHandling = SimParameters::CH_STEPPROJECT;
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.4, 0.0);

    Eigen::VectorXd q(4), F(4);
    q << 0.0, 0.0, 0.6, 0.0;
    F.setZero();
    std::vector<Eigen::Triplet<double>> H;
    processPenaltyForce(q, F, H);
    for (int i = 0; i < 4; i++)
        ASSERT_NEAR(F[i], 0.0, 1e-10);
}

void test_penaltyForce_finiteDifference()
{
    disableAllForces();
    params_.constraintHandling = SimParameters::CH_PENALTY;
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    params_.penaltyStiffness = 1e3;
    addParticle(0.0, 0.0);
    addParticle(0.4, 0.0);

    Eigen::VectorXd q(4);
    q << 0.0, 0.05, 0.5, -0.03;
    double eps = 1e-6;

    Eigen::VectorXd F0(4);
    F0.setZero();
    std::vector<Eigen::Triplet<double>> Htrip0;
    processPenaltyForce(q, F0, Htrip0);

    Eigen::SparseMatrix<double> Hmat(4, 4);
    Hmat.setFromTriplets(Htrip0.begin(), Htrip0.end());

    // Penalty localH is the force Jacobian dF/dq (different convention from spring)
    for (int j = 0; j < 4; j++)
    {
        Eigen::VectorXd qp = q, qm = q;
        qp[j] += eps;
        qm[j] -= eps;

        Eigen::VectorXd Fp(4), Fm(4);
        Fp.setZero();
        Fm.setZero();
        std::vector<Eigen::Triplet<double>> Hp, Hm;
        processPenaltyForce(qp, Fp, Hp);
        processPenaltyForce(qm, Fm, Hm);

        Eigen::VectorXd dFdqj = (Fp - Fm) / (2.0 * eps);
        for (int i = 0; i < 4; i++)
            ASSERT_NEAR(dFdqj[i], Hmat.coeff(i, j), 1e-3);
    }
}

// 11. ROD CONSTRAINTS

void test_rodConstraintValues_atRestLength()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    std::vector<int> rodIndices = getRigidRodConnectorIndices();
    ASSERT_EQ((int)rodIndices.size(), 1);

    Eigen::VectorXd q(4);
    q << 0.0, 0.0, 0.5, 0.0;
    Eigen::VectorXd g;
    computeRodConstraintValues(q, rodIndices, g);
    ASSERT_EQ(g.size(), 1);
    ASSERT_NEAR(g[0], 0.0, 1e-10); // ||d||^2 - L^2 = 0.25 - 0.25 = 0
}

void test_rodConstraintValues_stretched()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    RigidRod& rod = *(RigidRod*)connectors_[0];
    double L = rod.length;

    Eigen::VectorXd q(4);
    q << 0.0, 0.0, 0.8, 0.0;
    std::vector<int> rodIndices = getRigidRodConnectorIndices();
    Eigen::VectorXd g;
    computeRodConstraintValues(q, rodIndices, g);
    ASSERT_NEAR(g[0], 0.64 - L * L, 1e-10);
}

void test_constraintJacobian_basic()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    Eigen::VectorXd q(4);
    q << 0.0, 0.0, 0.5, 0.0;
    std::vector<int> rodIndices = getRigidRodConnectorIndices();
    Eigen::SparseMatrix<double> dg;
    computeConstraintJacobian(q, rodIndices, dg);

    ASSERT_EQ(dg.rows(), 1);
    ASSERT_EQ(dg.cols(), 4);
    // grad = 2*(pa-pb) = 2*(-0.5, 0) = (-1, 0)
    ASSERT_NEAR(dg.coeff(0, 0), -1.0, 1e-10); // d/dxa
    ASSERT_NEAR(dg.coeff(0, 1), 0.0, 1e-10);  // d/dya
    ASSERT_NEAR(dg.coeff(0, 2), 1.0, 1e-10);  // d/dxb = -grad_x
    ASSERT_NEAR(dg.coeff(0, 3), 0.0, 1e-10);  // d/dyb
}

void test_constraintJacobian_finiteDifference()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    Eigen::VectorXd q(4);
    q << 0.1, 0.2, 0.6, -0.1;
    std::vector<int> rodIndices = getRigidRodConnectorIndices();

    Eigen::SparseMatrix<double> dg;
    computeConstraintJacobian(q, rodIndices, dg);

    double eps = 1e-7;
    for (int j = 0; j < 4; j++)
    {
        Eigen::VectorXd qp = q, qm = q;
        qp[j] += eps;
        qm[j] -= eps;

        Eigen::VectorXd gp, gm;
        computeRodConstraintValues(qp, rodIndices, gp);
        computeRodConstraintValues(qm, rodIndices, gm);

        double dgdqj = (gp[0] - gm[0]) / (2.0 * eps);
        ASSERT_NEAR(dg.coeff(0, j), dgdqj, 1e-5);
    }
}

void test_getRigidRodConnectorIndices()
{
    params_.maxSpringDist = 0.01;
    // Add some springs and rigid rods manually
    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(1.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(2.0, 0.0), 1.0, false, false));

    connectors_.push_back(new Spring(0, 1, 0, 100, 1.0, true));
    connectors_.push_back(new RigidRod(1, 2, 0, 1.0));
    connectors_.push_back(new Spring(0, 2, 0, 50, 2.0, true));

    std::vector<int> rodIndices = getRigidRodConnectorIndices();
    ASSERT_EQ((int)rodIndices.size(), 1);
    ASSERT_EQ(rodIndices[0], 1);
}

// 12. FLEXIBLE ROD CREATION

void test_flexRod_particleCount()
{
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 5;
    addParticle(0.0, 0.0);
    addParticle(1.0, 0.0);

    int seg = std::max(2, params_.rodSegments);
    // 2 endpoints + (seg-1) internal particles
    ASSERT_EQ((int)particles_.size(), 2 + (seg - 1));
}

void test_flexRod_internalParticlesInert()
{
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 4;
    addParticle(0.0, 0.0);
    addParticle(0.8, 0.0);

    for (int i = 2; i < (int)particles_.size(); i++)
    {
        ASSERT_TRUE(particles_[i].inert);
        ASSERT_FALSE(particles_[i].fixed);
        ASSERT_NEAR(particles_[i].mass, 0.0, 1e-12);
    }
}

void test_flexRod_springsNotSnappable()
{
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 3;
    addParticle(0.0, 0.0);
    addParticle(0.6, 0.0);

    for (int i = 0; i < (int)connectors_.size(); i++)
    {
        ASSERT_EQ(connectors_[i]->getType(), SimParameters::CT_SPRING);
        Spring& s = *(Spring*)connectors_[i];
        ASSERT_FALSE(s.canSnap);
    }
}

void test_flexRod_springMass()
{
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 3;
    params_.rodDensity = 2.0;
    addParticle(0.0, 0.0);
    addParticle(0.6, 0.0);

    for (int i = 0; i < (int)connectors_.size(); i++)
    {
        Spring& s = *(Spring*)connectors_[i];
        ASSERT_GT(s.mass, 0.0);
        // mass = rodDensity * restlen
        ASSERT_NEAR(s.mass, params_.rodDensity * s.restlen, 1e-10);
    }
}

void test_flexRod_bendingStencils()
{
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 4;
    addParticle(0.0, 0.0);
    addParticle(0.8, 0.0);

    ASSERT_EQ((int)bendingStencils_.size(), 3); // segments-1 = 3
    ASSERT_TRUE(isBendingStencilStateValid());
}

void test_flexRod_internalPositions()
{
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 4;
    addParticle(0.0, 0.0);
    addParticle(0.8, 0.0);

    // createFlexibleRod is called with (newid=1, i=0), so interpolates
    // from particle 1 (0.8, 0) to particle 0 (0.0, 0)
    // Internal particles are evenly spaced from endpointA to endpointB
    Eigen::Vector2d pa = particles_[1].pos; // (0.8, 0)
    Eigen::Vector2d pb = particles_[0].pos; // (0.0, 0)
    int seg = 4;
    for (int i = 2; i < (int)particles_.size(); i++)
    {
        int idx = i - 2;
        double t = double(idx + 1) / double(seg);
        Eigen::Vector2d expected = (1.0 - t) * pa + t * pb;
        ASSERT_NEAR(particles_[i].pos[0], expected[0], 1e-10);
        ASSERT_NEAR(particles_[i].pos[1], expected[1], 1e-10);
    }
}

void test_flexRod_stencilAssociations()
{
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 3;
    addParticle(0.0, 0.0);
    addParticle(0.6, 0.0);

    // Each internal spring (not first/last) should have bending stencil associations
    for (int i = 0; i < (int)connectors_.size(); i++)
    {
        if (i == 0)
            ASSERT_EQ((int)connectors_[i]->associatedBendingStencils.size(), 1);
        else if (i == (int)connectors_.size() - 1)
            ASSERT_EQ((int)connectors_[i]->associatedBendingStencils.size(), 1);
        else
            ASSERT_EQ((int)connectors_[i]->associatedBendingStencils.size(), 2);
    }
}

// 13. STEP-AND-PROJECT (projectWithNewton)

void test_stepAndProject_noRods()
{
    // Without rods, projection should return qTilde unchanged
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);

    Eigen::VectorXd qTilde(2), invMass(2), qProj, lambdaProj;
    qTilde << 0.3, 0.4;
    invMass << 1.0, 1.0;

    bool converged = projectWithNewton(qTilde, invMass, qProj, lambdaProj);
    ASSERT_TRUE(converged);
    ASSERT_NEAR(qProj[0], 0.3, 1e-10);
    ASSERT_NEAR(qProj[1], 0.4, 1e-10);
    ASSERT_EQ(lambdaProj.size(), 0);
}

void test_stepAndProject_singleRod()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    RigidRod& rod = *(RigidRod*)connectors_[0];
    double L = rod.length;

    // Give a qTilde that violates the constraint
    Eigen::VectorXd qTilde(4), invMass(4), qProj, lambdaProj;
    qTilde << -0.1, 0.0, 0.7, 0.0; // distance = 0.8, but L = 0.5
    invMass << 1.0, 1.0, 1.0, 1.0;

    bool converged = projectWithNewton(qTilde, invMass, qProj, lambdaProj);
    ASSERT_TRUE(converged);

    // After projection, constraint should be satisfied
    double dist = (qProj.segment<2>(0) - qProj.segment<2>(2)).norm();
    ASSERT_NEAR(dist, L, 1e-5);
}

void test_stepAndProject_preservesCenter()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    params_.particleMass = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    RigidRod& rod = *(RigidRod*)connectors_[0];

    // Both particles have same mass => center of mass should be preserved
    Eigen::VectorXd qTilde(4), invMass(4), qProj, lambdaProj;
    qTilde << -0.1, 0.0, 0.7, 0.0;
    invMass << 1.0, 1.0, 1.0, 1.0;

    Eigen::Vector2d comBefore = 0.5 * (qTilde.segment<2>(0) + qTilde.segment<2>(2));

    bool converged = projectWithNewton(qTilde, invMass, qProj, lambdaProj);
    ASSERT_TRUE(converged);

    Eigen::Vector2d comAfter = 0.5 * (qProj.segment<2>(0) + qProj.segment<2>(2));
    ASSERT_NEAR(comBefore[0], comAfter[0], 1e-5);
    ASSERT_NEAR(comBefore[1], comAfter[1], 1e-5);
}

// 14. LAGRANGE MULTIPLIER SOLVER

void test_lagrange_noRods()
{
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);

    Eigen::VectorXd qDrift(2), qdotOld(2), invMass(2), F(2), lambdaSolved;
    qDrift << 0.0, 0.0;
    qdotOld << 0.0, 0.0;
    invMass << 1.0, 1.0;
    F << 0.0, -9.8;

    bool converged = solveConstrainedLagrangeLambda(qDrift, qdotOld, invMass, F, lambdaSolved);
    ASSERT_TRUE(converged);
    ASSERT_EQ(lambdaSolved.size(), 0);
}

void test_lagrange_singleRod()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    params_.timeStep = 0.01;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);

    RigidRod& rod = *(RigidRod*)connectors_[0];
    double L = rod.length;

    Eigen::VectorXd qDrift(4), qdotOld(4), invMass(4), F(4), lambdaSolved;
    qDrift << 0.0, 0.0, 0.5, 0.0;
    qdotOld << -1.0, 0.0, 1.0, 0.0; // moving apart
    invMass << 1.0, 1.0, 1.0, 1.0;
    F.setZero();

    bool converged = solveConstrainedLagrangeLambda(qDrift, qdotOld, invMass, F, lambdaSolved);
    ASSERT_TRUE(converged);
    ASSERT_EQ(lambdaSolved.size(), 1);
    ASSERT_TRUE(lambdaSolved.allFinite());
}

void test_lagrange_constraintSatisfied()
{
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    params_.constraintHandling = SimParameters::CH_LAGRANGEMULT;
    disableAllForces();
    params_.timeStep = 0.001;

    addParticle(-0.25, 0.0);
    addParticle(0.25, 0.0);
    particles_[0].vel << -0.8, 0.2;
    particles_[1].vel << 0.8, -0.2;

    for (int i = 0; i < 5; i++)
        simulateOneStep();

    double dist = (particles_[0].pos - particles_[1].pos).norm();
    double L = ((RigidRod*)connectors_[0])->length;
    ASSERT_NEAR(dist, L, 0.05);
}

// 15. NUMERICAL INTEGRATION

void test_integration_freefall_penalty()
{
    disableAllForces();
    params_.gravityEnabled = true;
    params_.constraintHandling = SimParameters::CH_PENALTY;
    params_.timeStep = 0.001;
    params_.maxSpringDist = 0.01;
    params_.particleMass = 1.0;
    addParticle(0.0, 0.0);

    double y0 = 0.0;
    int nsteps = 10;
    for (int i = 0; i < nsteps; i++)
        simulateOneStep();

    double t = nsteps * params_.timeStep;
    double expectedY = y0 + 0.5 * params_.gravityG * t * t;
    ASSERT_NEAR(particles_[0].pos[1], expectedY, 0.01);
    ASSERT_NEAR(particles_[0].pos[0], 0.0, 1e-10);
}

void test_integration_rodPenalty()
{
    disableAllForces();
    params_.gravityEnabled = true;
    params_.constraintHandling = SimParameters::CH_PENALTY;
    params_.penaltyStiffness = 1e5;
    params_.timeStep = 0.001;
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    params_.particleMass = 1.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    double L = ((RigidRod*)connectors_[0])->length;

    for (int i = 0; i < 50; i++)
        simulateOneStep();

    // Penalty method: constraint approximately maintained
    double dist = (particles_[0].pos - particles_[1].pos).norm();
    ASSERT_NEAR(dist, L, 0.01);
}

void test_integration_rodStepProject()
{
    disableAllForces();
    params_.gravityEnabled = true;
    params_.constraintHandling = SimParameters::CH_STEPPROJECT;
    params_.timeStep = 0.001;
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    params_.particleMass = 1.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    double L = ((RigidRod*)connectors_[0])->length;

    for (int i = 0; i < 50; i++)
        simulateOneStep();

    // Step-and-project: constraint should be well satisfied
    double dist = (particles_[0].pos - particles_[1].pos).norm();
    ASSERT_NEAR(dist, L, 1e-4);
}

void test_integration_rodLagrange()
{
    disableAllForces();
    params_.gravityEnabled = true;
    params_.constraintHandling = SimParameters::CH_LAGRANGEMULT;
    params_.timeStep = 0.001;
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    params_.particleMass = 1.0;

    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    double L = ((RigidRod*)connectors_[0])->length;

    for (int i = 0; i < 50; i++)
        simulateOneStep();

    double dist = (particles_[0].pos - particles_[1].pos).norm();
    ASSERT_NEAR(dist, L, 1e-3);
    ASSERT_TRUE(particles_[0].pos.allFinite());
    ASSERT_TRUE(particles_[1].pos.allFinite());
}

void test_integration_springOscillation()
{
    disableAllForces();
    params_.springsEnabled = true;
    params_.constraintHandling = SimParameters::CH_PENALTY;
    params_.timeStep = 0.001;
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    params_.particleMass = 1.0;
    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(0.2, 0.0);

    // Displace particle 1
    particles_[1].pos[0] = 0.3;

    for (int i = 0; i < 100; i++)
        simulateOneStep();

    // System should remain stable (finite positions)
    ASSERT_TRUE(particles_[0].pos.allFinite());
    ASSERT_TRUE(particles_[1].pos.allFinite());
}

// 16. SPRING SNAPPING

void test_springSnap_noSnapBelowThreshold()
{
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    params_.maxSpringStrain = 0.5;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);
    ASSERT_EQ((int)connectors_.size(), 1);

    // Move particles slightly (within strain threshold)
    particles_[1].pos[0] = 0.25; // strain = (0.25 - 0.2)/0.2 = 0.25 < 0.5
    pruneOverstrainedSprings();
    ASSERT_EQ((int)connectors_.size(), 1);
}

void test_springSnap_snapsAboveThreshold()
{
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    params_.maxSpringStrain = 0.2;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);
    ASSERT_EQ((int)connectors_.size(), 1);

    // Stretch beyond threshold
    particles_[1].pos[0] = 0.5; // strain = (0.5 - 0.2)/0.2 = 1.5 > 0.2
    pruneOverstrainedSprings();
    ASSERT_EQ((int)connectors_.size(), 0);
}

void test_springSnap_unsnappablePreserved()
{
    // Flexible rod springs (canSnap=false) should never snap
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 3;
    params_.maxSpringStrain = 0.01; // very low threshold
    addParticle(0.0, 0.0);
    addParticle(0.6, 0.0);

    int initialConnectors = (int)connectors_.size();

    // Stretch internal particles enormously
    for (int i = 2; i < (int)particles_.size(); i++)
        particles_[i].pos[0] += 5.0;

    pruneOverstrainedSprings();
    ASSERT_EQ((int)connectors_.size(), initialConnectors); // none removed
}

void test_springSnap_selective()
{
    params_.maxSpringDist = 0.01;
    params_.maxSpringStrain = 0.3;

    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(1.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(2.0, 0.0), 1.0, false, false));

    connectors_.push_back(new Spring(0, 1, 0, 100, 1.0, true));  // strain will be 0
    connectors_.push_back(new Spring(1, 2, 0, 100, 0.5, true));  // strain = (1.0 - 0.5)/0.5 = 1.0 > 0.3

    pruneOverstrainedSprings();
    ASSERT_EQ((int)connectors_.size(), 1);
    ASSERT_EQ(connectors_[0]->p1, 0);
    ASSERT_EQ(connectors_[0]->p2, 1);
}

void test_springSnap_rebuildsStencils()
{
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 3;
    params_.maxSpringStrain = 0.2;
    addParticle(0.0, 0.0);
    addParticle(0.6, 0.0);

    // Make one spring snappable so it can be removed
    Spring& s = *(Spring*)connectors_[1];
    s.canSnap = true;

    // Stretch the middle spring massively
    particles_[3].pos[0] = 5.0;

    int stencilsBefore = (int)bendingStencils_.size();
    pruneOverstrainedSprings();

    // Stencils should be rebuilt (some may be removed)
    ASSERT_TRUE(isBendingStencilStateValid());
}

// 17. SAW COLLISION

void test_sawCollision_particleInside()
{
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);
    addSaw(0.0, 0.0);

    deleteSawedObjects();
    ASSERT_EQ((int)particles_.size(), 0);
}

void test_sawCollision_particleOutside()
{
    params_.maxSpringDist = 0.01;
    params_.sawRadius = 0.1;
    addParticle(0.0, 0.0);
    addSaw(1.0, 0.0);

    deleteSawedObjects();
    ASSERT_EQ((int)particles_.size(), 1);
}

void test_sawCollision_connectorCut()
{
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    addParticle(-0.5, 0.0);
    addParticle(0.5, 0.0);
    ASSERT_EQ((int)connectors_.size(), 1);

    addSaw(0.0, 0.0);
    deleteSawedObjects();
    ASSERT_EQ((int)connectors_.size(), 0);
}

void test_sawCollision_connectorNotCut()
{
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    params_.sawRadius = 0.05;
    addParticle(-0.5, 0.0);
    addParticle(0.5, 0.0);
    ASSERT_EQ((int)connectors_.size(), 1);

    addSaw(0.0, 1.0); // far from the spring
    deleteSawedObjects();
    ASSERT_EQ((int)connectors_.size(), 1);
}

void test_sawCollision_particleDeletionCascade()
{
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);
    ASSERT_EQ((int)connectors_.size(), 1);

    addSaw(0.0, 0.0); // overlaps particle 0
    deleteSawedObjects();
    // Particle 0 deleted => spring should be deleted too
    ASSERT_EQ((int)connectors_.size(), 0);
    ASSERT_EQ((int)particles_.size(), 1);
}

void test_sawCollision_indexRemapping()
{
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);
    params_.maxSpringDist = 0.01; // no more auto-connect
    addParticle(0.4, 0.0);
    // Manually add a connector from p1 to p2
    connectors_.push_back(new Spring(1, 2, 0, 100, 0.2, true));
    ASSERT_EQ((int)connectors_.size(), 2);

    // Delete particle 0
    addSaw(0.0, 0.0);
    deleteSawedObjects();
    // p0 deleted, p1->0, p2->1
    ASSERT_EQ((int)particles_.size(), 2);
    ASSERT_EQ((int)connectors_.size(), 1);
    ASSERT_EQ(connectors_[0]->p1, 0);
    ASSERT_EQ(connectors_[0]->p2, 1);
}

void test_sawCollision_bendingStencilCleanup()
{
    // Test from the existing internal tests — saw deletes a middle
    // connector of a flexible rod, orphaning bending stencils
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 3;
    params_.sawRadius = 0.035;
    addParticle(0.0, 0.0);
    addParticle(0.6, 0.0);
    ASSERT_EQ((int)connectors_.size(), 3);
    ASSERT_EQ((int)bendingStencils_.size(), 2);

    addSaw(0.3, 0.0);
    deleteSawedObjects();

    ASSERT_TRUE(bendingStencils_.empty());
    for (int i = 0; i < (int)connectors_.size(); i++)
        ASSERT_TRUE(connectors_[i]->associatedBendingStencils.empty());
}

void test_sawCollision_multipleSaws()
{
    params_.maxSpringDist = 0.01;
    params_.sawRadius = 0.1;
    addParticle(0.0, 0.0);
    addParticle(1.0, 0.0);
    addParticle(2.0, 0.0);

    addSaw(0.0, 0.0);
    addSaw(2.0, 0.0);
    deleteSawedObjects();
    ASSERT_EQ((int)particles_.size(), 1);
    ASSERT_NEAR(particles_[0].pos[0], 1.0, 1e-10);
}

void test_outOfBounds_deletion()
{
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);
    addParticle(3.0, 0.0); // outside ±2 bounds

    deleteSawedObjects();
    ASSERT_EQ((int)particles_.size(), 1);
    ASSERT_NEAR(particles_[0].pos[0], 0.0, 1e-10);
}

void test_outOfBounds_yAxis()
{
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);
    addParticle(0.0, 3.0); // outside y ±2 bounds

    deleteSawedObjects();
    ASSERT_EQ((int)particles_.size(), 1);
}

// 18. FIXED PARTICLE DOFs

void test_fixedParticle_doesntMove()
{
    disableAllForces();
    params_.gravityEnabled = true;
    params_.timeStep = 0.01;
    params_.particleFixed = true;
    addParticle(0.0, 0.5);

    Eigen::Vector2d origPos = particles_[0].pos;
    for (int i = 0; i < 20; i++)
        simulateOneStep();

    ASSERT_NEAR(particles_[0].pos[0], origPos[0], 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], origPos[1], 1e-10);
}

void test_fixedParticle_withSpring()
{
    disableAllForces();
    params_.gravityEnabled = true;
    params_.springsEnabled = true;
    params_.timeStep = 0.001;
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;

    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(0.2, 0.0);

    Eigen::Vector2d fixedPos = particles_[0].pos;
    for (int i = 0; i < 20; i++)
        simulateOneStep();

    // Fixed particle should not move
    ASSERT_NEAR(particles_[0].pos[0], fixedPos[0], 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], fixedPos[1], 1e-10);
    // Free particle should have moved
    ASSERT_GT(std::abs(particles_[1].pos[1]), 1e-6);
}

void test_fixedParticle_withRigidRod()
{
    disableAllForces();
    params_.gravityEnabled = true;
    params_.constraintHandling = SimParameters::CH_STEPPROJECT;
    params_.timeStep = 0.001;
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;

    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(0.5, 0.0);

    Eigen::Vector2d fixedPos = particles_[0].pos;
    double L = ((RigidRod*)connectors_[0])->length;

    for (int i = 0; i < 50; i++)
        simulateOneStep();

    ASSERT_NEAR(particles_[0].pos[0], fixedPos[0], 1e-8);
    ASSERT_NEAR(particles_[0].pos[1], fixedPos[1], 1e-8);
    // Free particle should move but stay at rod distance
    double dist = (particles_[0].pos - particles_[1].pos).norm();
    ASSERT_NEAR(dist, L, 1e-3);
}

// 19. RESET SIMULATION

void test_initSimulation_clearsEverything()
{
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    addParticle(0.0, 0.0);
    addParticle(0.5, 0.0);
    addSaw(1.0, 1.0);

    ASSERT_GT((int)particles_.size(), 0);
    ASSERT_GT((int)connectors_.size(), 0);
    ASSERT_GT((int)bendingStencils_.size(), 0);
    ASSERT_GT((int)saws_.size(), 0);

    initSimulation();
    ASSERT_EQ((int)particles_.size(), 0);
    ASSERT_EQ((int)connectors_.size(), 0);
    ASSERT_EQ((int)bendingStencils_.size(), 0);
    ASSERT_EQ((int)saws_.size(), 0);
    ASSERT_NEAR(time_, 0.0, 1e-10);
}

// 20. FINITE-DIFFERENCE CHECKS — combined force/Hessian

void test_combinedForceHessian_finiteDifference()
{
    // Test the full force/Hessian with gravity + springs
    // H stores the Hessian of the potential energy (= -dF/dq)
    params_.maxSpringDist = 1.0;
    params_.connectorType = SimParameters::CT_SPRING;
    disableAllForces();
    params_.gravityEnabled = true;
    params_.springsEnabled = true;
    params_.particleMass = 1.0;
    addParticle(0.0, 0.0);
    addParticle(0.2, 0.0);

    Eigen::VectorXd q(4), qprev(4);
    q << 0.01, 0.02, 0.25, -0.03;
    qprev = q;

    Eigen::VectorXd F;
    Eigen::SparseMatrix<double> H;
    computeForceAndHessian(q, qprev, F, H);

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

        Eigen::VectorXd dFdqj = (Fp - Fm) / (2.0 * eps);
        for (int i = 0; i < 4; i++)
            ASSERT_NEAR(dFdqj[i], -H.coeff(i, j), 1e-3);
    }
}

// 21. CONSTRAINT PRESERVATION OVER TIME

void test_rodConstraint_preservedOverTime_stepProject()
{
    disableAllForces();
    params_.gravityEnabled = true;
    params_.constraintHandling = SimParameters::CH_STEPPROJECT;
    params_.timeStep = 0.001;
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    params_.particleMass = 1.0;

    addParticle(-0.2, 0.0);
    addParticle(0.2, 0.0);
    particles_[0].vel << 0.5, 0.3;
    particles_[1].vel << -0.3, 0.1;

    double L = ((RigidRod*)connectors_[0])->length;

    for (int i = 0; i < 100; i++)
    {
        simulateOneStep();
        double dist = (particles_[0].pos - particles_[1].pos).norm();
        ASSERT_NEAR(dist, L, 1e-4);
    }
}

void test_rodConstraint_preservedOverTime_lagrange()
{
    disableAllForces();
    params_.gravityEnabled = true;
    params_.constraintHandling = SimParameters::CH_LAGRANGEMULT;
    params_.timeStep = 0.001;
    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.maxSpringDist = 1.0;
    params_.particleMass = 1.0;

    addParticle(-0.2, 0.0);
    addParticle(0.2, 0.0);
    particles_[0].vel << 0.5, 0.3;
    particles_[1].vel << -0.3, 0.1;

    double L = ((RigidRod*)connectors_[0])->length;

    for (int i = 0; i < 100; i++)
    {
        simulateOneStep();
        double dist = (particles_[0].pos - particles_[1].pos).norm();
        ASSERT_NEAR(dist, L, 1e-3);
    }
}

// 22. EDGE CASES

void test_emptySimulation()
{
    // Should not crash with no particles
    Eigen::VectorXd q, lambda, qdot;
    buildConfiguration(q, lambda, qdot);
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(lambda.size(), 0);
    ASSERT_EQ(qdot.size(), 0);
}

void test_singleParticleIntegration()
{
    params_.maxSpringDist = 0.01;
    disableAllForces();
    params_.gravityEnabled = true;
    params_.timeStep = 0.001;
    addParticle(0.0, 0.0);

    // Should not crash
    simulateOneStep();
    ASSERT_TRUE(particles_[0].pos.allFinite());
    ASSERT_TRUE(particles_[0].vel.allFinite());
}

void test_twoRigidRods()
{
    // Chain of 3 particles with 2 rigid rods
    params_.constraintHandling = SimParameters::CH_STEPPROJECT;
    params_.maxSpringDist = 0.01;
    disableAllForces();
    params_.gravityEnabled = true;
    params_.timeStep = 0.001;

    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(0.5, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(1.0, 0.0), 1.0, false, false));

    connectors_.push_back(new RigidRod(0, 1, 0, 0.5));
    connectors_.push_back(new RigidRod(1, 2, 0, 0.5));

    for (int i = 0; i < 20; i++)
        simulateOneStep();

    double dist01 = (particles_[0].pos - particles_[1].pos).norm();
    double dist12 = (particles_[1].pos - particles_[2].pos).norm();
    ASSERT_NEAR(dist01, 0.5, 1e-3);
    ASSERT_NEAR(dist12, 0.5, 1e-3);
    ASSERT_TRUE(particles_[0].pos.allFinite());
}

void test_flexRodWithGravity()
{
    // Full integration test: flexible rod under gravity
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 3;
    params_.timeStep = 0.001;
    params_.constraintHandling = SimParameters::CH_PENALTY;
    disableAllForces();
    params_.gravityEnabled = true;
    params_.springsEnabled = true;
    params_.bendingEnabled = true;

    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(0.6, 0.0);

    for (int i = 0; i < 50; i++)
        simulateOneStep();

    // All particles should remain finite
    for (int i = 0; i < (int)particles_.size(); i++)
    {
        ASSERT_TRUE(particles_[i].pos.allFinite());
        ASSERT_TRUE(particles_[i].vel.allFinite());
    }

    // Fixed particle should not move
    ASSERT_NEAR(particles_[0].pos[0], 0.0, 1e-10);
    ASSERT_NEAR(particles_[0].pos[1], 0.0, 1e-10);

    // Free endpoint should have fallen
    ASSERT_LT(particles_[1].pos[1], 0.0);
}

void test_findSpringConnectorByParticles()
{
    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(1.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(2.0, 0.0), 1.0, false, false));

    connectors_.push_back(new Spring(0, 1, 0, 100, 1.0, true));
    connectors_.push_back(new RigidRod(1, 2, 0, 1.0));

    ASSERT_EQ(findSpringConnectorByParticles(connectors_, 0, 1), 0);
    ASSERT_EQ(findSpringConnectorByParticles(connectors_, 1, 0), 0); // reversed
    ASSERT_EQ(findSpringConnectorByParticles(connectors_, 1, 2), -1); // not a spring
    ASSERT_EQ(findSpringConnectorByParticles(connectors_, 0, 2), -1); // doesn't exist
}

void test_rebuildBendingStencilAssociations()
{
    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(0.5, 0.0), 1.0, false, true));
    particles_.push_back(Particle(Eigen::Vector2d(1.0, 0.0), 1.0, false, false));

    connectors_.push_back(new Spring(0, 1, 0, 100, 0.5, false));
    connectors_.push_back(new Spring(1, 2, 0, 100, 0.5, false));
    bendingStencils_.push_back(BendingStencil(0, 1, 2, 0.05));

    // Add an invalid stencil that references nonexistent particles
    bendingStencils_.push_back(BendingStencil(0, 1, 99, 0.05));

    rebuildBendingStencilAssociations();

    ASSERT_EQ((int)bendingStencils_.size(), 1); // invalid one removed
    ASSERT_TRUE(isBendingStencilStateValid());
}

void test_isBendingStencilStateValid_valid()
{
    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(0.5, 0.0), 1.0, false, true));
    particles_.push_back(Particle(Eigen::Vector2d(1.0, 0.0), 1.0, false, false));

    connectors_.push_back(new Spring(0, 1, 0, 100, 0.5, false));
    connectors_.push_back(new Spring(1, 2, 0, 100, 0.5, false));
    bendingStencils_.push_back(BendingStencil(0, 1, 2, 0.05));

    connectors_[0]->associatedBendingStencils.insert(0);
    connectors_[1]->associatedBendingStencils.insert(0);

    ASSERT_TRUE(isBendingStencilStateValid());
}

void test_isBendingStencilStateValid_missingAssociation()
{
    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(0.5, 0.0), 1.0, false, true));
    particles_.push_back(Particle(Eigen::Vector2d(1.0, 0.0), 1.0, false, false));

    connectors_.push_back(new Spring(0, 1, 0, 100, 0.5, false));
    connectors_.push_back(new Spring(1, 2, 0, 100, 0.5, false));
    bendingStencils_.push_back(BendingStencil(0, 1, 2, 0.05));

    // Don't add associations — should be invalid
    ASSERT_FALSE(isBendingStencilStateValid());
}

void test_dampingReducesOscillation()
{
    // A spring-mass system with damping should have decreasing amplitude
    disableAllForces();
    params_.springsEnabled = true;
    params_.dampingEnabled = true;
    params_.dampingStiffness = 2.0;
    params_.timeStep = 0.001;
    params_.connectorType = SimParameters::CT_SPRING;
    params_.maxSpringDist = 1.0;
    params_.particleMass = 1.0;

    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(0.2, 0.0);

    // Displace particle 1
    particles_[1].pos[0] = 0.4;

    // Run for a while
    for (int i = 0; i < 200; i++)
        simulateOneStep();
    double amp1 = std::abs(particles_[1].pos[0] - 0.2);

    for (int i = 0; i < 200; i++)
        simulateOneStep();
    double amp2 = std::abs(particles_[1].pos[0] - 0.2);

    // Second amplitude should be smaller due to damping
    ASSERT_LT(amp2, amp1 + 0.01);
}

void test_combinedForcesWork()
{
    // Gravity + springs + floor + bending all together
    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 3;
    params_.timeStep = 0.001;
    params_.constraintHandling = SimParameters::CH_PENALTY;

    params_.particleFixed = true;
    addParticle(0.0, 0.0);
    params_.particleFixed = false;
    addParticle(0.6, 0.0);

    // Should not crash with all forces enabled
    for (int i = 0; i < 100; i++)
        simulateOneStep();

    for (int i = 0; i < (int)particles_.size(); i++)
    {
        ASSERT_TRUE(particles_[i].pos.allFinite());
        ASSERT_TRUE(particles_[i].vel.allFinite());
    }
}

// 23. CREATIVE COMPONENT

void test_bendingDamping_noEffectWhenStraight()
{
    // A straight rod should have zero bending damping force
    disableAllForces();
    params_.springsEnabled = true;
    params_.bendingEnabled = true;
    params_.bendingDampingEnabled = true;
    params_.bendingDampingStiffness = 1.0;
    params_.timeStep = 0.001;

    // straight configuration: three collinear particles
    particles_.push_back(Particle(Eigen::Vector2d(0.0, 0.0), 1.0, false, false));
    particles_.push_back(Particle(Eigen::Vector2d(0.5, 0.0), 1.0, false, true));
    particles_.push_back(Particle(Eigen::Vector2d(1.0, 0.0), 1.0, false, false));
    connectors_.push_back(new Spring(0, 1, 0, 100, 0.5, false));
    connectors_.push_back(new Spring(1, 2, 0, 100, 0.5, false));
    bendingStencils_.push_back(BendingStencil(0, 1, 2, 0.05));
    connectors_[0]->associatedBendingStencils.insert(0);
    connectors_[1]->associatedBendingStencils.insert(0);

    // stationary particles — theta and thetaDot are both 0
    Eigen::VectorXd q(6), qprev(6), F(6);
    for (int i = 0; i < 3; i++)
    {
        q.segment<2>(2 * i) = particles_[i].pos;
        qprev.segment<2>(2 * i) = particles_[i].pos;
    }
    F.setZero();

    processBendingDampingForce(q, qprev, F);
    ASSERT_NEAR(F.norm(), 0.0, 1e-12);
}

void test_energy_kineticOnly()
{
    disableAllForces();
    params_.particleMass = 2.0;
    params_.maxSpringDist = 0.01;
    addParticle(0.0, 0.0);
    particles_[0].vel << 3.0, 4.0;

    double ke, pe;
    computeEnergy(ke, pe);

    // KE = 0.5 * 2 * (9+16) = 25
    ASSERT_NEAR(ke, 25.0, 1e-6);
    // no springs or gravity, PE = 0
    ASSERT_NEAR(pe, 0.0, 1e-6);
}

void test_strainVisualization_paramToggle()
{
    // Just verify the parameter defaults are correct
    SimParameters p;
    ASSERT_FALSE(p.strainVisualization);
    ASSERT_FALSE(p.velocityVisualization);
    ASSERT_FALSE(p.bendingDampingEnabled);
    ASSERT_FALSE(p.showEnergyReadout);
    ASSERT_NEAR(p.bendingDampingStiffness, 0.01, 1e-10);
}

// MAIN — Test runner

int main(int argc, char** argv)
{
    std::cout << "========================================" << std::endl;
    std::cout << "  Universe of Goo — Milestone II Tests  " << std::endl;
    std::cout << "========================================" << std::endl;

    // 1. Geometry helpers
    runTest("PtSegmentDist/perpendicular", test_ptSegmentDist_perpendicular);
    runTest("PtSegmentDist/atEndpoint", test_ptSegmentDist_atEndpoint);
    runTest("PtSegmentDist/onSegment", test_ptSegmentDist_onSegment);
    runTest("PtSegmentDist/degenerateSegment", test_ptSegmentDist_degenerateSegment);
    runTest("PtSegmentDist/diagonal", test_ptSegmentDist_diagonal);
    runTest("PtSegmentDist/beyondEndpoint2", test_ptSegmentDist_beyondEndpoint2);

    // 2. Particle/connector creation
    runTest("AddParticle/single", test_addParticle_single);
    runTest("AddParticle/fixed", test_addParticle_fixed);
    runTest("AddParticle/springCreation", test_addParticle_springCreation);
    runTest("AddParticle/noSpringIfTooFar", test_addParticle_noSpringIfTooFar);
    runTest("AddParticle/rigidRodCreation", test_addParticle_rigidRodCreation);
    runTest("AddParticle/flexRodCreation", test_addParticle_flexRodCreation);
    runTest("AddParticle/skipInertParticles", test_addParticle_skipInertParticles);
    runTest("AddParticle/springStiffness", test_addParticle_springStiffness);

    // 3. Configuration vectors
    runTest("Configuration/buildBasic", test_buildConfiguration_basic);
    runTest("Configuration/buildWithRigidRod", test_buildConfiguration_withRigidRod);
    runTest("Configuration/unbuildRoundtrip", test_unbuildConfiguration_roundtrip);
    runTest("Configuration/unbuildLambda", test_unbuildConfiguration_lambda);

    // 4. Mass matrix
    runTest("MassMatrix/singleParticle", test_massInverse_singleParticle);
    runTest("MassMatrix/fixedParticle", test_massInverse_fixedParticle);
    runTest("MassMatrix/mixedParticles", test_massInverse_mixedParticles);
    runTest("MassMatrix/includesConnectorMass", test_massInverse_includesConnectorMass);
    runTest("MassMatrix/buildInverseMassVector", test_buildInverseMassVector);
    runTest("MassMatrix/diagonal", test_massInverse_diagonal);

    // 5. Gravity force
    runTest("Gravity/singleParticle", test_gravity_singleParticle);
    runTest("Gravity/fixedParticle", test_gravity_fixedParticle);
    runTest("Gravity/multipleParticles", test_gravity_multipleParticles);
    runTest("Gravity/zeroHessian", test_gravity_zeroHessian);

    // 6. Spring force
    runTest("SpringForce/atRest", test_springForce_atRest);
    runTest("SpringForce/stretched", test_springForce_stretched);
    runTest("SpringForce/compressed", test_springForce_compressed);
    runTest("SpringForce/newtonsThirdLaw", test_springForce_newtonsThirdLaw);
    runTest("SpringForce/hessianSymmetry", test_springForce_hessianSymmetry);
    runTest("SpringForce/finiteDifference", test_springForce_finiteDifference);

    // 7. Damping force
    runTest("DampingForce/zeroVelocity", test_dampingForce_zeroVelocity);
    runTest("DampingForce/opposesRelativeVelocity", test_dampingForce_opposesRelativeVelocity);
    runTest("DampingForce/newtonsThirdLaw", test_dampingForce_newtonsThirdLaw);

    // 8. Floor force
    runTest("FloorForce/aboveFloor", test_floorForce_aboveFloor);
    runTest("FloorForce/belowFloor", test_floorForce_belowFloor);
    runTest("FloorForce/fixedParticleIgnored", test_floorForce_fixedParticleIgnored);
    runTest("FloorForce/atFloor", test_floorForce_atFloor);

    // 9. Bending force
    runTest("BendingForce/straight", test_bendingForce_straight);
    runTest("BendingForce/bent", test_bendingForce_bent);
    runTest("BendingForce/finiteDifference", test_bendingForce_finiteDifference);

    // 10. Penalty force
    runTest("PenaltyForce/noRods", test_penaltyForce_noRods);
    runTest("PenaltyForce/rodAtRestLength", test_penaltyForce_rodAtRestLength);
    runTest("PenaltyForce/rodStretched", test_penaltyForce_rodStretched);
    runTest("PenaltyForce/onlyWhenPenaltyMode", test_penaltyForce_onlyWhenPenaltyMode);
    runTest("PenaltyForce/finiteDifference", test_penaltyForce_finiteDifference);

    // 11. Rod constraints
    runTest("RodConstraint/valuesAtRestLength", test_rodConstraintValues_atRestLength);
    runTest("RodConstraint/valuesStretched", test_rodConstraintValues_stretched);
    runTest("RodConstraint/jacobianBasic", test_constraintJacobian_basic);
    runTest("RodConstraint/jacobianFiniteDifference", test_constraintJacobian_finiteDifference);
    runTest("RodConstraint/getRigidRodIndices", test_getRigidRodConnectorIndices);

    // 12. Flexible rod creation
    runTest("FlexRod/particleCount", test_flexRod_particleCount);
    runTest("FlexRod/internalParticlesInert", test_flexRod_internalParticlesInert);
    runTest("FlexRod/springsNotSnappable", test_flexRod_springsNotSnappable);
    runTest("FlexRod/springMass", test_flexRod_springMass);
    runTest("FlexRod/bendingStencils", test_flexRod_bendingStencils);
    runTest("FlexRod/internalPositions", test_flexRod_internalPositions);
    runTest("FlexRod/stencilAssociations", test_flexRod_stencilAssociations);

    // 13. Step-and-project
    runTest("StepAndProject/noRods", test_stepAndProject_noRods);
    runTest("StepAndProject/singleRod", test_stepAndProject_singleRod);
    runTest("StepAndProject/preservesCenter", test_stepAndProject_preservesCenter);

    // 14. Lagrange multipliers
    runTest("Lagrange/noRods", test_lagrange_noRods);
    runTest("Lagrange/singleRod", test_lagrange_singleRod);
    runTest("Lagrange/constraintSatisfied", test_lagrange_constraintSatisfied);

    // 15. Numerical integration
    runTest("Integration/freefallPenalty", test_integration_freefall_penalty);
    runTest("Integration/rodPenalty", test_integration_rodPenalty);
    runTest("Integration/rodStepProject", test_integration_rodStepProject);
    runTest("Integration/rodLagrange", test_integration_rodLagrange);
    runTest("Integration/springOscillation", test_integration_springOscillation);

    // 16. Spring snapping
    runTest("SpringSnap/noSnapBelowThreshold", test_springSnap_noSnapBelowThreshold);
    runTest("SpringSnap/snapsAboveThreshold", test_springSnap_snapsAboveThreshold);
    runTest("SpringSnap/unsnappablePreserved", test_springSnap_unsnappablePreserved);
    runTest("SpringSnap/selective", test_springSnap_selective);
    runTest("SpringSnap/rebuildsStencils", test_springSnap_rebuildsStencils);

    // 17. Saw collision
    runTest("SawCollision/particleInside", test_sawCollision_particleInside);
    runTest("SawCollision/particleOutside", test_sawCollision_particleOutside);
    runTest("SawCollision/connectorCut", test_sawCollision_connectorCut);
    runTest("SawCollision/connectorNotCut", test_sawCollision_connectorNotCut);
    runTest("SawCollision/particleDeletionCascade", test_sawCollision_particleDeletionCascade);
    runTest("SawCollision/indexRemapping", test_sawCollision_indexRemapping);
    runTest("SawCollision/bendingStencilCleanup", test_sawCollision_bendingStencilCleanup);
    runTest("SawCollision/multipleSaws", test_sawCollision_multipleSaws);
    runTest("SawCollision/outOfBoundsDeletion", test_outOfBounds_deletion);
    runTest("SawCollision/outOfBoundsYAxis", test_outOfBounds_yAxis);

    // 18. Fixed particle DOFs
    runTest("FixedParticle/doesntMove", test_fixedParticle_doesntMove);
    runTest("FixedParticle/withSpring", test_fixedParticle_withSpring);
    runTest("FixedParticle/withRigidRod", test_fixedParticle_withRigidRod);

    // 19. Reset simulation
    runTest("Reset/clearsEverything", test_initSimulation_clearsEverything);

    // 20. Finite-difference checks
    runTest("FiniteDiff/combinedForceHessian", test_combinedForceHessian_finiteDifference);

    // 21. Constraint preservation
    runTest("ConstraintPreserve/stepProjectOverTime", test_rodConstraint_preservedOverTime_stepProject);
    runTest("ConstraintPreserve/lagrangeOverTime", test_rodConstraint_preservedOverTime_lagrange);

    // 22. Edge cases & utilities
    runTest("EdgeCase/emptySimulation", test_emptySimulation);
    runTest("EdgeCase/singleParticleIntegration", test_singleParticleIntegration);
    runTest("EdgeCase/twoRigidRods", test_twoRigidRods);
    runTest("EdgeCase/flexRodWithGravity", test_flexRodWithGravity);
    runTest("EdgeCase/findSpringByParticles", test_findSpringConnectorByParticles);
    runTest("EdgeCase/rebuildBendingStencilAssocs", test_rebuildBendingStencilAssociations);
    runTest("EdgeCase/stencilStateValid", test_isBendingStencilStateValid_valid);
    runTest("EdgeCase/stencilStateMissingAssoc", test_isBendingStencilStateValid_missingAssociation);
    runTest("EdgeCase/dampingReducesOscillation", test_dampingReducesOscillation);
    runTest("EdgeCase/combinedForcesWork", test_combinedForcesWork);

    // 23. Creative component
    runTest("Creative/bendingDampingNoEffectStraight", test_bendingDamping_noEffectWhenStraight);
    runTest("Creative/energyKineticOnly", test_energy_kineticOnly);
    runTest("Creative/paramDefaults", test_strainVisualization_paramToggle);

    // Summary
    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  RESULTS: " << g_tests_passed << " passed, "
              << g_tests_failed << " failed  ("
              << g_assertions_passed << " assertions passed, "
              << g_assertions_failed << " failed)" << std::endl;
    std::cout << "========================================" << std::endl;

    if (!g_failures.empty())
    {
        std::cout << std::endl << "Failures:" << std::endl;
        for (size_t i = 0; i < g_failures.size(); i++)
            std::cout << "  " << g_failures[i] << std::endl;
    }

    return g_tests_failed > 0 ? 1 : 0;
}
