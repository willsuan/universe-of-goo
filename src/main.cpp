#include "polyscope/polyscope.h"

#include "polyscope/messages.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

#include <iostream>
#include <unordered_set>
#include <utility>
#include <deque>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "SimParameters.h"
#include "SceneObjects.h"

bool running_;
SimParameters params_;
double time_;
std::vector<Particle, Eigen::aligned_allocator<Particle> > particles_;
std::vector<Connector*> connectors_;
std::vector<Saw> saws_;

Eigen::MatrixXd renderQ;
Eigen::MatrixXi renderF;
Eigen::MatrixXd renderC;

struct MouseClick
{
    double x;
    double y;
    SimParameters::ClickMode mode;
};

std::deque<MouseClick> mouseClicks_;

double getTotalParticleMass(int idx)
{
    double mass = particles_[idx].mass;
    for (std::vector<Connector*>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
    {
        if ((*it)->p1 == idx || (*it)->p2 == idx)
            mass += 0.5 * (*it)->mass;
    }
    return mass;
}


void initSimulation()
{
    time_ = 0;
    particles_.clear();
    for (std::vector<Connector*>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
        delete* it;
    connectors_.clear();
    saws_.clear();
}

void updateRenderGeometry()
{
    double baseradius = 0.02;
    double pulsefactor = 0.1;
    double pulsespeed = 50.0;

    int sawteeth = 20;
    double sawdepth = 0.1;
    double sawangspeed = 10.0;

    double baselinewidth = 0.005;

    int numcirclewedges = 20;

    // this is terrible. But, easiest to get up and running

    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector3d> vertexColors;
    std::vector<Eigen::Vector3i> faces;

    int idx = 0;

    double eps = 1e-4;


    if (params_.floorEnabled)
    {
        for (int i = 0; i < 6; i++)
        {
            vertexColors.push_back(Eigen::Vector3d(0.3, 1.0, 0.3));
        }

        verts.push_back(Eigen::Vector3d(-2, -0.5, eps));
        verts.push_back(Eigen::Vector3d(2, -0.5, eps));
        verts.push_back(Eigen::Vector3d(-2, -1, eps));

        faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));

        verts.push_back(Eigen::Vector3d(-2, -1, eps));
        verts.push_back(Eigen::Vector3d(2, -0.5, eps));
        verts.push_back(Eigen::Vector3d(2, -1, eps));
        faces.push_back(Eigen::Vector3i(idx + 3, idx + 4, idx + 5));
        idx += 6;
    }


    for (std::vector<Connector*>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
    {
        Eigen::Vector3d color;
        if ((*it)->associatedBendingStencils.empty())
            color << 0.0, 0.0, 1.0;
        else
            color << 0.75, 0.5, 0.75;
        Eigen::Vector2d sourcepos = particles_[(*it)->p1].pos;
        Eigen::Vector2d destpos = particles_[(*it)->p2].pos;

        Eigen::Vector2d vec = destpos - sourcepos;
        Eigen::Vector2d perp(-vec[1], vec[0]);
        perp /= perp.norm();

        double dist = (sourcepos - destpos).norm();

        double width = baselinewidth / (1.0 + 20.0 * dist * dist);

        for (int i = 0; i < 4; i++)
            vertexColors.push_back(color);

        verts.push_back(Eigen::Vector3d(sourcepos[0] + width * perp[0], sourcepos[1] + width * perp[1], -eps));
        verts.push_back(Eigen::Vector3d(sourcepos[0] - width * perp[0], sourcepos[1] - width * perp[1], -eps));
        verts.push_back(Eigen::Vector3d(destpos[0] + width * perp[0], destpos[1] + width * perp[1], -eps));
        verts.push_back(Eigen::Vector3d(destpos[0] - width * perp[0], destpos[1] - width * perp[1], -eps));

        faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));
        faces.push_back(Eigen::Vector3i(idx + 2, idx + 1, idx + 3));
        idx += 4;
    }

    int nparticles = particles_.size();

    for (int i = 0; i < nparticles; i++)
    {
        double radius = baseradius * sqrt(getTotalParticleMass(i));
        radius *= (1.0 + pulsefactor * sin(pulsespeed * time_));

        Eigen::Vector3d color(0, 0, 0);

        if (particles_[i].fixed)
        {
            radius = baseradius;
            color << 1.0, 0, 0;
        }

        for (int j = 0; j < numcirclewedges + 2; j++)
        {
            vertexColors.push_back(color);
        }


        verts.push_back(Eigen::Vector3d(particles_[i].pos[0], particles_[i].pos[1], 0));

        const double PI = 3.1415926535898;
        for (int j = 0; j <= numcirclewedges; j++)
        {
            verts.push_back(Eigen::Vector3d(particles_[i].pos[0] + radius * cos(2 * PI * j / numcirclewedges),
                particles_[i].pos[1] + radius * sin(2 * PI * j / numcirclewedges), 0));
        }

        for (int j = 0; j <= numcirclewedges; j++)
        {
            faces.push_back(Eigen::Vector3i(idx, idx + j + 1, idx + 1 + ((j + 1) % (numcirclewedges + 1))));
        }

        idx += numcirclewedges + 2;
    }

    for (std::vector<Saw>::iterator it = saws_.begin(); it != saws_.end(); ++it)
    {
        double outerradius = it->radius;
        double innerradius = (1.0 - sawdepth) * outerradius;

        Eigen::Vector3d color(0.5, 0.5, 0.5);

        int spokes = 2 * sawteeth;
        for (int j = 0; j < spokes + 2; j++)
        {
            vertexColors.push_back(color);
        }

        verts.push_back(Eigen::Vector3d(it->pos[0], it->pos[1], 0));

        const double PI = 3.1415926535898;
        for (int i = 0; i <= spokes; i++)
        {
            double radius = (i % 2 == 0) ? innerradius : outerradius;
            verts.push_back(Eigen::Vector3d(it->pos[0] + radius * cos(2 * PI * i / spokes + sawangspeed * time_),
                it->pos[1] + radius * sin(2 * PI * i / spokes + sawangspeed * time_), 0));
        }

        for (int j = 0; j <= spokes; j++)
        {
            faces.push_back(Eigen::Vector3i(idx, idx + j + 1, idx + 1 + ((j + 1) % (spokes + 1))));
        }

        idx += spokes + 2;
    }

    renderQ.resize(verts.size(), 3);
    renderC.resize(vertexColors.size(), 3);
    for (int i = 0; i < verts.size(); i++)
    {
        renderQ.row(i) = verts[i];
        renderC.row(i) = vertexColors[i];
    }
    renderF.resize(faces.size(), 3);
    for (int i = 0; i < faces.size(); i++)
        renderF.row(i) = faces[i];
}

void addParticle(double x, double y)
{
    Eigen::Vector2d newpos(x, y);
    double mass = params_.particleMass;
    if (params_.particleFixed)
        mass = std::numeric_limits<double>::infinity();

    int newid = particles_.size();
    particles_.push_back(Particle(newpos, mass, params_.particleFixed, false));

    // Connect new particle to all existing particles within maxSpringDist.
    // Each spring's rest length equals the current distance between the two
    // particles it connects. Spring stiffness is kb / L, where kb =
    // springStiffness and L = rest length, so longer springs are weaker.
    // The spring mass is set to 0 (massless springs; mass is tracked on
    // particles only).
    for (int i = 0; i < newid; i++)
    {
        if (particles_[i].inert)
            continue;
        double dist = (particles_[i].pos - newpos).norm();
        if (dist <= params_.maxSpringDist)
        {
            double restlen = dist;
            double stiffness = (restlen > 1e-10) ? params_.springStiffness / restlen : params_.springStiffness;
            connectors_.push_back(new Spring(newid, i, 0.0, stiffness, restlen, true));
        }
    }
}

void addSaw(double x, double y)
{
    saws_.push_back(Saw(Eigen::Vector2d(x, y), params_.sawRadius));
}

// Packs particle positions, previous positions, and velocities into flat
// vectors of dimension 2*n (where n = number of particles). The layout is:
//   q     = [x0, y0, x1, y1, ..., x_{n-1}, y_{n-1}]
//   qprev = same layout, but using each particle's previous-step position
//   qdot  = same layout, but using each particle's velocity
// Fixed particles ARE included; their DOFs will simply never change because
// forces on them are zeroed out and their velocities stay at zero.
void buildConfiguration(Eigen::VectorXd& q, Eigen::VectorXd& qprev, Eigen::VectorXd& qdot)
{
    int n = particles_.size();
    q.resize(2 * n);
    qprev.resize(2 * n);
    qdot.resize(2 * n);
    for (int i = 0; i < n; i++)
    {
        q[2 * i]     = particles_[i].pos[0];
        q[2 * i + 1] = particles_[i].pos[1];
        qprev[2 * i]     = particles_[i].prevpos[0];
        qprev[2 * i + 1] = particles_[i].prevpos[1];
        qdot[2 * i]     = particles_[i].vel[0];
        qdot[2 * i + 1] = particles_[i].vel[1];
    }
}

// Reverses buildConfiguration: writes the flat configuration vectors back
// into the per-particle pos and vel fields so the renderer can display the
// updated state. Also saves the current position as prevpos for the damping
// force computation at the next time step.
void unbuildConfiguration(const Eigen::VectorXd& q, const Eigen::VectorXd& qdot)
{
    int n = particles_.size();
    for (int i = 0; i < n; i++)
    {
        particles_[i].prevpos = particles_[i].pos;
        particles_[i].pos[0] = q[2 * i];
        particles_[i].pos[1] = q[2 * i + 1];
        particles_[i].vel[0] = qdot[2 * i];
        particles_[i].vel[1] = qdot[2 * i + 1];
    }
}

// Inverse mass matrix (sparse). M is a 2n x 2n diagonal matrix where each
// particle i contributes M[2i,2i] = M[2i+1,2i+1] = m_i. For fixed particles,
// mass = infinity, so M^{-1} entry = 0, which correctly prevents any force
// from accelerating them. We build M^{-1} directly using triplets.
void computeMassInverse(Eigen::SparseMatrix<double>& Minv)
{
    int n = particles_.size();
    Minv.resize(2 * n, 2 * n);
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(2 * n);
    for (int i = 0; i < n; i++)
    {
        if (!particles_[i].fixed)
        {
            double invmass = 1.0 / getTotalParticleMass(i);
            trips.push_back(Eigen::Triplet<double>(2 * i,     2 * i,     invmass));
            trips.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, invmass));
        }
    }
    Minv.setFromTriplets(trips.begin(), trips.end());
}


// Computes the total configurational force F = -dV/dq and the force Jacobian
// H = dF/dq, assembled from all enabled force types. Each force type
// contributes additively to F and H. We use Eigen triplets to build H as a
// sparse matrix, since most entries are zero (each force only couples a small
// number of DOFs).
//
// Force types:
//   1. Gravity       -- constant downward force on each non-fixed particle
//   2. Springs       -- Hookean elastic springs with stiffness kb/L
//   3. Damping       -- viscous damping approximated via position differences
//   4. Floor penalty -- one-sided quadratic penalty force at y = -0.5
void computeForceAndHessian(const Eigen::VectorXd& q, const Eigen::VectorXd& qprev, Eigen::VectorXd& F, Eigen::SparseMatrix<double>& H)
{
    int n = particles_.size();
    int ndofs = 2 * n;
    F.resize(ndofs);
    F.setZero();
    std::vector<Eigen::Triplet<double>> Htrips;

    // 1. GRAVITY: V = sum_i (-m_i * g * y_i), F_y_i = m_i * g, Hessian = 0
    if (params_.gravityEnabled)
    {
        for (int i = 0; i < n; i++)
        {
            if (!particles_[i].fixed)
                F[2 * i + 1] += getTotalParticleMass(i) * params_.gravityG;
        }
    }

    // 2. SPRINGS: V = (k/2)(||p2-p1|| - L)^2, k = stiffness, L = rest length.
    // Let d = p2-p1, r = ||d||. Force on p1: f1 = k(r-L)*d/r.
    // Force Jacobian block: let B = k*[(1-L/r)*I + (L/r^3)*d*d^T].
    // Then dF1/dp1 = -B, dF1/dp2 = +B, dF2/dp1 = +B, dF2/dp2 = -B.
    if (params_.springsEnabled)
    {
        for (std::vector<Connector*>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
        {
            Spring* s = dynamic_cast<Spring*>(*it);
            if (!s) continue;

            int i1 = s->p1;
            int i2 = s->p2;
            double k = s->stiffness;
            double L = s->restlen;

            Eigen::Vector2d p1(q[2 * i1], q[2 * i1 + 1]);
            Eigen::Vector2d p2(q[2 * i2], q[2 * i2 + 1]);
            Eigen::Vector2d d = p2 - p1;
            double r = d.norm();

            if (r < 1e-10) continue;

            Eigen::Vector2d dhat = d / r;
            double forcemag = k * (r - L);
            Eigen::Vector2d f1 = forcemag * dhat;
            Eigen::Vector2d f2 = -forcemag * dhat;

            F[2 * i1]     += f1[0];
            F[2 * i1 + 1] += f1[1];
            F[2 * i2]     += f2[0];
            F[2 * i2 + 1] += f2[1];

            Eigen::Matrix2d I2 = Eigen::Matrix2d::Identity();
            Eigen::Matrix2d ddT = d * d.transpose();
            Eigen::Matrix2d B = k * ((1.0 - L / r) * I2 + (L / (r * r * r)) * ddT);

            for (int a = 0; a < 2; a++)
            {
                for (int b = 0; b < 2; b++)
                {
                    Htrips.push_back(Eigen::Triplet<double>(2 * i1 + a, 2 * i1 + b, -B(a, b)));
                    Htrips.push_back(Eigen::Triplet<double>(2 * i1 + a, 2 * i2 + b,  B(a, b)));
                    Htrips.push_back(Eigen::Triplet<double>(2 * i2 + a, 2 * i1 + b,  B(a, b)));
                    Htrips.push_back(Eigen::Triplet<double>(2 * i2 + a, 2 * i2 + b, -B(a, b)));
                }
            }
        }
    }

    // 3. VISCOUS DAMPING: F1 = kdamp*(v2-v1), F2 = kdamp*(v1-v2), where
    // v_i = (q_i - q_i_prev) / h (finite-difference velocity approximation).
    // Hessian: dF1/dq1 = -(kdamp/h)*I, dF1/dq2 = +(kdamp/h)*I, etc.
    if (params_.dampingEnabled)
    {
        double h = params_.timeStep;
        double kdamp = params_.dampingStiffness;

        for (std::vector<Connector*>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
        {
            Spring* s = dynamic_cast<Spring*>(*it);
            if (!s) continue;

            int i1 = s->p1;
            int i2 = s->p2;

            Eigen::Vector2d v1 = (Eigen::Vector2d(q[2 * i1], q[2 * i1 + 1]) -
                                  Eigen::Vector2d(qprev[2 * i1], qprev[2 * i1 + 1])) / h;
            Eigen::Vector2d v2 = (Eigen::Vector2d(q[2 * i2], q[2 * i2 + 1]) -
                                  Eigen::Vector2d(qprev[2 * i2], qprev[2 * i2 + 1])) / h;

            Eigen::Vector2d f1 = kdamp * (v2 - v1);
            Eigen::Vector2d f2 = kdamp * (v1 - v2);

            F[2 * i1]     += f1[0];
            F[2 * i1 + 1] += f1[1];
            F[2 * i2]     += f2[0];
            F[2 * i2 + 1] += f2[1];

            double coeff = kdamp / h;
            for (int a = 0; a < 2; a++)
            {
                Htrips.push_back(Eigen::Triplet<double>(2 * i1 + a, 2 * i1 + a, -coeff));
                Htrips.push_back(Eigen::Triplet<double>(2 * i1 + a, 2 * i2 + a,  coeff));
                Htrips.push_back(Eigen::Triplet<double>(2 * i2 + a, 2 * i1 + a,  coeff));
                Htrips.push_back(Eigen::Triplet<double>(2 * i2 + a, 2 * i2 + a, -coeff));
            }
        }
    }

    // 4. FLOOR PENALTY: one-sided quadratic penalty at y = -0.5.
    // V = (kfloor/2)*(y+0.5)^2 if y < -0.5, else 0.
    // F_y = kfloor*(-0.5-y), dF_y/dy = -kfloor.
    // kfloor = 10000 prevents 100 kg particles from reaching y = -1.0.
    if (params_.floorEnabled)
    {
        double floorY = -0.5;
        double kfloor = 10000.0;

        for (int i = 0; i < n; i++)
        {
            if (particles_[i].fixed) continue;

            double yi = q[2 * i + 1];
            if (yi < floorY)
            {
                double penetration = yi - floorY;
                F[2 * i + 1] += -kfloor * penetration;
                Htrips.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, -kfloor));
            }
        }
    }

    // Zero out forces on fixed particles (safety measure on top of M^{-1}=0)
    for (int i = 0; i < n; i++)
    {
        if (particles_[i].fixed)
        {
            F[2 * i] = 0;
            F[2 * i + 1] = 0;
        }
    }

    H.resize(ndofs, ndofs);
    H.setFromTriplets(Htrips.begin(), Htrips.end());
}

// Advances the simulation by one time step of size h = params_.timeStep,
// using the integrator selected in params_.integrator. Four methods are
// supported: Explicit Euler, Velocity Verlet, Implicit Euler, and Implicit
// Midpoint. The latter two require Newton's method to solve nonlinear systems.
void numericalIntegration(Eigen::VectorXd& q, Eigen::VectorXd& qprev, Eigen::VectorXd& qdot)
{
    double h = params_.timeStep;
    int ndofs = q.size();

    Eigen::SparseMatrix<double> Minv;
    computeMassInverse(Minv);

    switch (params_.integrator)
    {
    // Explicit Euler: q_{i+1} = q_i + h*qdot_i, qdot_{i+1} = qdot_i + h*M^{-1}*F(q_i)
    // First-order, conditionally stable. No linear system solve needed.
    case SimParameters::TI_EXPLICIT_EULER:
    {
        Eigen::VectorXd F;
        Eigen::SparseMatrix<double> H;
        computeForceAndHessian(q, qprev, F, H);

        Eigen::VectorXd qnew = q + h * qdot;
        Eigen::VectorXd qdotnew = qdot + h * (Minv * F);

        q = qnew;
        qdot = qdotnew;
        break;
    }

    // Velocity Verlet: position update first, then force at new position.
    // Symplectic, better energy conservation than explicit Euler.
    case SimParameters::TI_VELOCITY_VERLET:
    {
        Eigen::VectorXd qnew = q + h * qdot;

        Eigen::VectorXd F;
        Eigen::SparseMatrix<double> H;
        computeForceAndHessian(qnew, q, F, H);

        Eigen::VectorXd qdotnew = qdot + h * (Minv * F);

        q = qnew;
        qdot = qdotnew;
        break;
    }

    // Implicit Euler: q_{i+1} = q_i + h*qdot_{i+1}, qdot_{i+1} = qdot_i + h*M^{-1}*F(q_{i+1}).
    // Unconditionally stable. Solved via Newton's method on the residual
    // g(q) = q - q_old - h*v_old - h^2*M^{-1}*F(q), with Jacobian
    // dg/dq = I - h^2*M^{-1}*dF/dq.
    case SimParameters::TI_IMPLICIT_EULER:
    {
        Eigen::VectorXd qguess = q + h * qdot;

        for (int iter = 0; iter < params_.NewtonMaxIters; iter++)
        {
            Eigen::VectorXd F;
            Eigen::SparseMatrix<double> dF;
            computeForceAndHessian(qguess, qprev, F, dF);

            Eigen::VectorXd residual = qguess - q - h * qdot - h * h * (Minv * F);

            if (residual.norm() < params_.NewtonTolerance)
                break;

            Eigen::SparseMatrix<double> I(ndofs, ndofs);
            I.setIdentity();
            Eigen::SparseMatrix<double> A = I - h * h * (Minv * dF);

            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            Eigen::VectorXd dq = solver.solve(-residual);

            qguess += dq;
        }

        qdot = (qguess - q) / h;
        q = qguess;
        break;
    }

    // Implicit Midpoint: q_{i+1} = q_i + h*(qdot_{i+1}+qdot_i)/2,
    // qdot_{i+1} = qdot_i + h*M^{-1}*F((q_{i+1}+q_i)/2).
    // Second-order symplectic integrator. Newton's method Jacobian:
    // dg/dq_{i+1} = I - (h^2/4)*M^{-1}*dF(qmid).
    case SimParameters::TI_IMPLICIT_MIDPOINT:
    {
        Eigen::VectorXd qguess = q + h * qdot;

        for (int iter = 0; iter < params_.NewtonMaxIters; iter++)
        {
            Eigen::VectorXd qmid = 0.5 * (qguess + q);
            Eigen::VectorXd qprevmid = 0.5 * (qprev + q);

            Eigen::VectorXd F;
            Eigen::SparseMatrix<double> dF;
            computeForceAndHessian(qmid, qprevmid, F, dF);

            Eigen::VectorXd residual = qguess - q - h * qdot - 0.5 * h * h * (Minv * F);

            if (residual.norm() < params_.NewtonTolerance)
                break;

            Eigen::SparseMatrix<double> I(ndofs, ndofs);
            I.setIdentity();
            Eigen::SparseMatrix<double> A = I - 0.25 * h * h * (Minv * dF);

            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            Eigen::VectorXd dq = solver.solve(-residual);

            qguess += dq;
        }

        qdot = (qguess - q) / h;
        q = qguess;
        break;
    }
    }

    // Zero out velocity on fixed particles (safety on top of M^{-1}=0)
    for (int i = 0; i < (int)particles_.size(); i++)
    {
        if (particles_[i].fixed)
        {
            qdot[2 * i] = 0;
            qdot[2 * i + 1] = 0;
            q[2 * i] = particles_[i].pos[0];
            q[2 * i + 1] = particles_[i].pos[1];
        }
    }
}


// Computes the shortest distance from point p to the infinite line through
// q1 and q2. The closest point c = q1 + t*(q2-q1) satisfies (c-p).(q2-q1)=0,
// giving t = (p-q1).(q2-q1) / ||q2-q1||^2.
double distancePointToLine(const Eigen::Vector2d& p, const Eigen::Vector2d& q1, const Eigen::Vector2d& q2)
{
    Eigen::Vector2d lineDir = q2 - q1;
    double lineLenSq = lineDir.squaredNorm();
    if (lineLenSq < 1e-20) return (p - q1).norm();

    double t = (p - q1).dot(lineDir) / lineLenSq;
    Eigen::Vector2d closest = q1 + t * lineDir;
    return (p - closest).norm();
}

// Computes the shortest distance from point p to the finite line segment
// connecting q1 to q2. Same as distancePointToLine but with t clamped to [0,1].
double distancePointToSegment(const Eigen::Vector2d& p, const Eigen::Vector2d& q1, const Eigen::Vector2d& q2)
{
    Eigen::Vector2d lineDir = q2 - q1;
    double lineLenSq = lineDir.squaredNorm();
    if (lineLenSq < 1e-20) return (p - q1).norm();

    double t = (p - q1).dot(lineDir) / lineLenSq;
    t = std::max(0.0, std::min(1.0, t));
    Eigen::Vector2d closest = q1 + t * lineDir;
    return (p - closest).norm();
}

// Detects and removes all particles and springs that are touching a saw or
// have drifted out of bounds. Also handles index remapping after particle
// deletion so spring endpoint indices remain valid.
void deleteSawedObjects()
{
    if (saws_.empty() && particles_.empty()) return;

    int n = particles_.size();
    std::vector<bool> particleDead(n, false);

    for (int i = 0; i < n; i++)
    {
        // Out-of-bounds check
        if (std::abs(particles_[i].pos[0]) > 5.0 || std::abs(particles_[i].pos[1]) > 5.0)
        {
            particleDead[i] = true;
            continue;
        }

        // Saw-particle collision: particle center inside saw disk
        for (std::vector<Saw>::iterator sit = saws_.begin(); sit != saws_.end(); ++sit)
        {
            if ((particles_[i].pos - sit->pos).norm() < sit->radius)
            {
                particleDead[i] = true;
                break;
            }
        }
    }

    std::vector<bool> connectorDead(connectors_.size(), false);

    for (int i = 0; i < (int)connectors_.size(); i++)
    {
        if (particleDead[connectors_[i]->p1] || particleDead[connectors_[i]->p2])
        {
            connectorDead[i] = true;
            continue;
        }

        Eigen::Vector2d p1 = particles_[connectors_[i]->p1].pos;
        Eigen::Vector2d p2 = particles_[connectors_[i]->p2].pos;

        for (std::vector<Saw>::iterator sit = saws_.begin(); sit != saws_.end(); ++sit)
        {
            double dist = distancePointToSegment(sit->pos, p1, p2);
            if (dist < sit->radius)
            {
                connectorDead[i] = true;
                break;
            }
        }
    }

    for (int i = (int)connectors_.size() - 1; i >= 0; i--)
    {
        if (connectorDead[i])
        {
            delete connectors_[i];
            connectors_.erase(connectors_.begin() + i);
        }
    }

    std::vector<int> indexMap(n, -1);
    int newIdx = 0;
    for (int i = 0; i < n; i++)
    {
        if (!particleDead[i])
            indexMap[i] = newIdx++;
    }

    for (int i = n - 1; i >= 0; i--)
    {
        if (particleDead[i])
            particles_.erase(particles_.begin() + i);
    }

    for (std::vector<Connector*>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
    {
        (*it)->p1 = indexMap[(*it)->p1];
        (*it)->p2 = indexMap[(*it)->p2];
    }
}

// After each time step, compute Cauchy strain epsilon = (L_current - L_rest) / L_rest
// for each spring. Delete any spring whose strain exceeds maxSpringStrain.
void pruneOverstrainedSprings()
{
    for (int i = (int)connectors_.size() - 1; i >= 0; i--)
    {
        Spring* s = dynamic_cast<Spring*>(connectors_[i]);
        if (!s || !s->canSnap) continue;

        Eigen::Vector2d p1 = particles_[s->p1].pos;
        Eigen::Vector2d p2 = particles_[s->p2].pos;
        double currentLen = (p2 - p1).norm();
        double strain = (currentLen - s->restlen) / s->restlen;

        if (strain > params_.maxSpringStrain)
        {
            delete connectors_[i];
            connectors_.erase(connectors_.begin() + i);
        }
    }
}

bool simulateOneStep()
{
    // Create configurational vectors
    Eigen::VectorXd q, qprev, v;
    buildConfiguration(q, qprev, v);
    // Use them for one step of time integration
    numericalIntegration(q, qprev, v);
    // Unpack the DOFs back into the particles for rendering
    unbuildConfiguration(q, v);

    // Cleanup: delete sawed objects and snapped springs
    pruneOverstrainedSprings();
    deleteSawedObjects();

    // Time advances
    time_ += params_.timeStep;
    return false;
}

void callback()
{
    ImGui::SetNextWindowSize(ImVec2(500., 0.));
    ImGui::Begin("UI", nullptr);

    if (ImGui::CollapsingHeader("Simulation Control", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Run/Pause Sim", ImVec2(-1, 0)))
        {
            running_ = !running_;
        }
        if (ImGui::Button("Reset Sim", ImVec2(-1, 0)))
        {
            running_ = false;
            initSimulation();
        }
    }
    if (ImGui::CollapsingHeader("UI Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Combo("Click Adds", (int*)&params_.clickMode, "Particles\0Saws\0\0");
    }
    if (ImGui::CollapsingHeader("Simulation Options"))
    {
        ImGui::InputDouble("Timestep", &params_.timeStep);
        ImGui::Combo("Integrator", (int*)&params_.integrator, "Explicit Euler\0Implicit Euler\0Implicit Midpoint\0Velocity Verlet\0\0");
        ImGui::InputDouble("Newton Tolerance", &params_.NewtonTolerance);
        ImGui::InputInt("Newton Max Iters", &params_.NewtonMaxIters);
    }
    if (ImGui::CollapsingHeader("Forces"))
    {
        ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
        ImGui::InputDouble("  Gravity g", &params_.gravityG);
        ImGui::Checkbox("Springs Enabled", &params_.springsEnabled);
        ImGui::InputDouble("  Max Strain", &params_.maxSpringStrain);
        ImGui::Checkbox("Damping Enabled", &params_.dampingEnabled);
        ImGui::InputDouble("  Viscosity", &params_.dampingStiffness);
        ImGui::Checkbox("Floor Enabled", &params_.floorEnabled);
    }


    if (ImGui::CollapsingHeader("New Particles"))
    {
        ImGui::Checkbox("Is Fixed", &params_.particleFixed);
        ImGui::InputDouble("Mass", &params_.particleMass);
    }

    if (ImGui::CollapsingHeader("New Saws"))
    {
        ImGui::InputDouble("Radius", &params_.sawRadius);
    }

    if (ImGui::CollapsingHeader("New Springs"))
    {
        ImGui::InputDouble("Max Spring Dist", &params_.maxSpringDist);
        ImGui::InputDouble("Base Stiffness", &params_.springStiffness);
    }

    ImGuiIO& io = ImGui::GetIO();
    io.DisplayFramebufferScale = ImVec2(1, 1);
    // this now only works on macs with retina displays - maybe funky on older macbook airs?
    // more robust solution is documented here: https://github.com/ocornut/imgui/issues/5081
    #if defined(__APPLE__)
        io.DisplayFramebufferScale = ImVec2(2,2);
    #endif

    if (io.MouseClicked[0] && !io.WantCaptureMouse) {
        MouseClick mc;
        glm::vec2 screenCoords{ io.MousePos.x * io.DisplayFramebufferScale.x, io.MousePos.y * io.DisplayFramebufferScale.y};

        glm::mat4 proj = polyscope::view::getCameraPerspectiveMatrix();

        glm::vec4 ndc{ -1.0f + 2.0f * screenCoords.x / (polyscope::view::bufferWidth  ) , 1.0f - 2.0f * screenCoords.y / (polyscope::view::bufferHeight ), 0, 1 };
        glm::vec4 camera = glm::inverse(proj) * ndc;
        mc.x = camera[0];
        mc.y = camera[1];
        mc.mode = params_.clickMode;
        mouseClicks_.push_back(mc);
    }

    ImGui::End();
}

int main(int argc, char **argv)
{
  polyscope::view::setWindowSize(1600, 800);
  polyscope::view::setWindowResizable(false);
  polyscope::view::style = polyscope::view::NavigateStyle::Planar;
  polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
  polyscope::options::buildGui = false;
  polyscope::options::openImGuiWindowForUserCallback = false;


  polyscope::options::autocenterStructures = false;
  polyscope::options::autoscaleStructures = false;

  initSimulation();

  polyscope::init();

  polyscope::options::automaticallyComputeSceneExtents = false;
  polyscope::state::lengthScale = 1.;
  polyscope::state::boundingBox =
      std::tuple<glm::vec3, glm::vec3>{ {-2., -1., -1.}, {2., 1., 1.} };

  polyscope::state::userCallback = callback;

  while (!polyscope::render::engine->windowRequestsClose())
  {
      if (running_)
          simulateOneStep();
      updateRenderGeometry();
      auto * surf = polyscope::registerSurfaceMesh("UI", renderQ, renderF);
      surf->setTransparency(0.9);
      auto * color = surf->addVertexColorQuantity("Colors", renderC);
      color->setEnabled(true);

      polyscope::frameTick();
      while (!mouseClicks_.empty())
      {
          MouseClick mc = mouseClicks_.front();
          mouseClicks_.pop_front();
          switch (mc.mode)
          {
          case SimParameters::ClickMode::CM_ADDPARTICLE:
          {
              addParticle(mc.x, mc.y);
              break;
          }
          case SimParameters::ClickMode::CM_ADDSAW:
          {
              addSaw(mc.x, mc.y);
              break;
          }
          }
      }
  }

  return 0;
}
