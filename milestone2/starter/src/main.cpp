#include "polyscope/polyscope.h"

#include "polyscope/messages.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_set>
#include <utility>
#include <deque>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/Dense>

#include "SimParameters.h"
#include "SceneObjects.h"

bool running_;
SimParameters params_;
double time_;
std::vector<Particle, Eigen::aligned_allocator<Particle> > particles_;
std::vector<Connector*> connectors_;
std::vector<Saw> saws_;
std::vector<BendingStencil> bendingStencils_;

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

double getInverseMassForParticle(int idx)
{
    double mass = getTotalParticleMass(idx);
    if (!std::isfinite(mass) || mass <= 0.0)
        return 0.0;
    return 1.0 / mass;
}

void buildInverseMassVector(Eigen::VectorXd& invMass)
{
    invMass.resize(2 * particles_.size());
    for (int i = 0; i < (int)particles_.size(); i++)
    {
        double invm = getInverseMassForParticle(i);
        invMass[2 * i] = invm;
        invMass[2 * i + 1] = invm;
    }
}

std::vector<int> getRigidRodConnectorIndices()
{
    std::vector<int> rodIndices;
    for (int i = 0; i < (int)connectors_.size(); i++)
    {
        if (connectors_[i]->getType() == SimParameters::CT_RIGIDROD)
            rodIndices.push_back(i);
    }
    return rodIndices;
}

int findSpringConnectorByParticles(const std::vector<Connector*>& connectorList, int pA, int pB)
{
    for (int i = 0; i < (int)connectorList.size(); i++)
    {
        if (connectorList[i]->getType() != SimParameters::CT_SPRING)
            continue;

        int c1 = connectorList[i]->p1;
        int c2 = connectorList[i]->p2;
        if ((c1 == pA && c2 == pB) || (c1 == pB && c2 == pA))
            return i;
    }
    return -1;
}

void rebuildBendingStencilAssociations()
{
    // quick cleanup pass so every stencil points to two real springs after deletes
    for (int i = 0; i < (int)connectors_.size(); i++)
        connectors_[i]->associatedBendingStencils.clear();

    std::vector<BendingStencil> validStencils;
    validStencils.reserve(bendingStencils_.size());
    for (int i = 0; i < (int)bendingStencils_.size(); i++)
    {
        const BendingStencil& stencil = bendingStencils_[i];
        if (stencil.p1 < 0 || stencil.p1 >= (int)particles_.size() ||
            stencil.p2 < 0 || stencil.p2 >= (int)particles_.size() ||
            stencil.p3 < 0 || stencil.p3 >= (int)particles_.size())
            continue;

        int springA = findSpringConnectorByParticles(connectors_, stencil.p1, stencil.p2);
        int springB = findSpringConnectorByParticles(connectors_, stencil.p2, stencil.p3);
        if (springA < 0 || springB < 0 || springA == springB)
            continue;

        validStencils.push_back(stencil);
    }

    bendingStencils_.swap(validStencils);

    for (int i = 0; i < (int)bendingStencils_.size(); i++)
    {
        const BendingStencil& stencil = bendingStencils_[i];
        int springA = findSpringConnectorByParticles(connectors_, stencil.p1, stencil.p2);
        int springB = findSpringConnectorByParticles(connectors_, stencil.p2, stencil.p3);
        if (springA >= 0)
            connectors_[springA]->associatedBendingStencils.insert(i);
        if (springB >= 0)
            connectors_[springB]->associatedBendingStencils.insert(i);
    }
}

void computeRodConstraintValues(const Eigen::VectorXd& q, const std::vector<int>& rodIndices, Eigen::VectorXd& g)
{
    g.resize(rodIndices.size());
    for (int i = 0; i < (int)rodIndices.size(); i++)
    {
        RigidRod& rod = *(RigidRod*)connectors_[rodIndices[i]];
        Eigen::Vector2d pa = q.segment<2>(2 * rod.p1);
        Eigen::Vector2d pb = q.segment<2>(2 * rod.p2);
        Eigen::Vector2d d = pa - pb;
        g[i] = d.dot(d) - rod.length * rod.length;
    }
}

void computeConstraintJacobian(const Eigen::VectorXd& q, const std::vector<int>& rodIndices, Eigen::SparseMatrix<double>& dg)
{
    int n = q.size();
    int m = rodIndices.size();
    dg.resize(m, n);
    dg.setZero();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(4 * m);
    for (int i = 0; i < m; i++)
    {
        RigidRod& rod = *(RigidRod*)connectors_[rodIndices[i]];
        Eigen::Vector2d pa = q.segment<2>(2 * rod.p1);
        Eigen::Vector2d pb = q.segment<2>(2 * rod.p2);
        Eigen::Vector2d grad = 2.0 * (pa - pb);

        triplets.push_back(Eigen::Triplet<double>(i, 2 * rod.p1, grad[0]));
        triplets.push_back(Eigen::Triplet<double>(i, 2 * rod.p1 + 1, grad[1]));
        triplets.push_back(Eigen::Triplet<double>(i, 2 * rod.p2, -grad[0]));
        triplets.push_back(Eigen::Triplet<double>(i, 2 * rod.p2 + 1, -grad[1]));
    }

    dg.setFromTriplets(triplets.begin(), triplets.end());
}

void createFlexibleRod(int endpointA, int endpointB)
{
    // we split the rod into equal chunks, then each chunk is just an unsnappable spring
    int segmentCount = std::max(2, params_.rodSegments);
    Eigen::Vector2d pa = particles_[endpointA].pos;
    Eigen::Vector2d pb = particles_[endpointB].pos;

    std::vector<int> chain;
    chain.reserve(segmentCount + 1);
    chain.push_back(endpointA);

    for (int i = 1; i < segmentCount; i++)
    {
        double t = double(i) / double(segmentCount);
        Eigen::Vector2d mid = (1.0 - t) * pa + t * pb;
        particles_.push_back(Particle(mid, 0.0, false, true));
        chain.push_back((int)particles_.size() - 1);
    }
    chain.push_back(endpointB);

    std::vector<int> springConnectorIds;
    springConnectorIds.reserve(segmentCount);
    for (int i = 0; i < segmentCount; i++)
    {
        int p1 = chain[i];
        int p2 = chain[i + 1];
        Eigen::Vector2d diff = particles_[p2].pos - particles_[p1].pos;
        double restLen = diff.norm();
        if (restLen <= 1e-12)
            continue;

        double springMass = params_.rodDensity * restLen;
        double springStiffness = params_.rodStretchingStiffness / restLen;
        connectors_.push_back(new Spring(p1, p2, springMass, springStiffness, restLen, false));
        springConnectorIds.push_back((int)connectors_.size() - 1);
    }

    for (int i = 0; i + 1 < (int)springConnectorIds.size(); i++)
    {
        Spring& s1 = *(Spring*)connectors_[springConnectorIds[i]];
        Spring& s2 = *(Spring*)connectors_[springConnectorIds[i + 1]];
        double denom = s1.restlen + s2.restlen;
        if (denom <= 1e-12)
            continue;

        double kb = 2.0 * params_.rodBendingStiffness / denom;
        bendingStencils_.push_back(BendingStencil(chain[i], chain[i + 1], chain[i + 2], kb));
        int stencilIdx = (int)bendingStencils_.size() - 1;
        connectors_[springConnectorIds[i]]->associatedBendingStencils.insert(stencilIdx);
        connectors_[springConnectorIds[i + 1]]->associatedBendingStencils.insert(stencilIdx);
    }
}

void processPenaltyForce(const Eigen::VectorXd& q, Eigen::VectorXd& F, std::vector<Eigen::Triplet<double>>& H)
{
    if (params_.constraintHandling != SimParameters::CH_PENALTY)
        return;

    std::vector<int> rodIndices = getRigidRodConnectorIndices();
    for (int i = 0; i < (int)rodIndices.size(); i++)
    {
        RigidRod& rod = *(RigidRod*)connectors_[rodIndices[i]];
        Eigen::Vector2d pa = q.segment<2>(2 * rod.p1);
        Eigen::Vector2d pb = q.segment<2>(2 * rod.p2);
        Eigen::Vector2d d = pa - pb;

        double g = d.dot(d) - rod.length * rod.length;
        Eigen::Vector2d force = -4.0 * params_.penaltyStiffness * g * d;
        F.segment<2>(2 * rod.p1) += force;
        F.segment<2>(2 * rod.p2) -= force;

        Eigen::Matrix2d I = Eigen::Matrix2d::Identity();
        Eigen::Matrix2d localH = -8.0 * params_.penaltyStiffness * (d * d.transpose()) - 4.0 * params_.penaltyStiffness * g * I;
        for (int r = 0; r < 2; r++)
            for (int c = 0; c < 2; c++)
            {
                H.push_back(Eigen::Triplet<double>(2 * rod.p1 + r, 2 * rod.p1 + c, localH(r, c)));
                H.push_back(Eigen::Triplet<double>(2 * rod.p2 + r, 2 * rod.p2 + c, localH(r, c)));
                H.push_back(Eigen::Triplet<double>(2 * rod.p1 + r, 2 * rod.p2 + c, -localH(r, c)));
                H.push_back(Eigen::Triplet<double>(2 * rod.p2 + r, 2 * rod.p1 + c, -localH(r, c)));
            }
    }
}

bool projectWithNewton(const Eigen::VectorXd& qTilde, const Eigen::VectorXd& invMass, Eigen::VectorXd& qProjected, Eigen::VectorXd& lambdaProjected)
{
    std::vector<int> rodIndices = getRigidRodConnectorIndices();
    int n = qTilde.size();
    int m = rodIndices.size();
    if (m == 0)
    {
        qProjected = qTilde;
        lambdaProjected.resize(0);
        return true;
    }

    qProjected = qTilde;
    lambdaProjected = Eigen::VectorXd::Zero(m);

    // this is the step-and-project Newton solve from the spec
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    for (int iter = 0; iter < params_.NewtonMaxIters; iter++)
    {
        Eigen::VectorXd g;
        Eigen::SparseMatrix<double> dg;
        computeRodConstraintValues(qProjected, rodIndices, g);
        computeConstraintJacobian(qProjected, rodIndices, dg);

        Eigen::VectorXd top = qProjected - qTilde;
        Eigen::VectorXd coupling = dg.transpose() * lambdaProjected;
        for (int i = 0; i < n; i++)
            top[i] += invMass[i] * coupling[i];

        Eigen::VectorXd f(n + m);
        f.head(n) = top;
        f.tail(m) = g;
        if (f.norm() < params_.NewtonTolerance)
            return true;

        std::vector<Eigen::Triplet<double>> dfTriplets;
        dfTriplets.reserve(n + 32 * m);
        for (int i = 0; i < n; i++)
            dfTriplets.push_back(Eigen::Triplet<double>(i, i, 1.0));

        for (int i = 0; i < m; i++)
        {
            if (std::abs(lambdaProjected[i]) < 1e-14)
                continue;
            RigidRod& rod = *(RigidRod*)connectors_[rodIndices[i]];
            double weight = 2.0 * lambdaProjected[i];

            int a0 = 2 * rod.p1;
            int a1 = a0 + 1;
            int b0 = 2 * rod.p2;
            int b1 = b0 + 1;

            dfTriplets.push_back(Eigen::Triplet<double>(a0, a0, weight * invMass[a0]));
            dfTriplets.push_back(Eigen::Triplet<double>(a0, b0, -weight * invMass[a0]));
            dfTriplets.push_back(Eigen::Triplet<double>(b0, a0, -weight * invMass[b0]));
            dfTriplets.push_back(Eigen::Triplet<double>(b0, b0, weight * invMass[b0]));

            dfTriplets.push_back(Eigen::Triplet<double>(a1, a1, weight * invMass[a1]));
            dfTriplets.push_back(Eigen::Triplet<double>(a1, b1, -weight * invMass[a1]));
            dfTriplets.push_back(Eigen::Triplet<double>(b1, a1, -weight * invMass[b1]));
            dfTriplets.push_back(Eigen::Triplet<double>(b1, b1, weight * invMass[b1]));
        }

        for (int outer = 0; outer < dg.outerSize(); outer++)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(dg, outer); it; ++it)
            {
                int row = it.row();
                int col = it.col();
                double val = it.value();
                dfTriplets.push_back(Eigen::Triplet<double>(col, n + row, invMass[col] * val));
                dfTriplets.push_back(Eigen::Triplet<double>(n + row, col, val));
            }
        }

        Eigen::SparseMatrix<double> df(n + m, n + m);
        df.setFromTriplets(dfTriplets.begin(), dfTriplets.end());

        solver.analyzePattern(df);
        solver.factorize(df);
        if (solver.info() != Eigen::Success)
            return false;

        Eigen::VectorXd delta = solver.solve(-f);
        if (solver.info() != Eigen::Success)
            return false;

        qProjected += delta.head(n);
        lambdaProjected += delta.tail(m);
    }

    Eigen::VectorXd finalG;
    computeRodConstraintValues(qProjected, rodIndices, finalG);
    return finalG.norm() < params_.NewtonTolerance;
}

bool solveConstrainedLagrangeLambda(const Eigen::VectorXd& qDrift, const Eigen::VectorXd& qdotOld, const Eigen::VectorXd& invMass, const Eigen::VectorXd& F, Eigen::VectorXd& lambdaSolved)
{
    std::vector<int> rodIndices = getRigidRodConnectorIndices();
    int n = qDrift.size();
    int m = rodIndices.size();
    if (m == 0)
    {
        lambdaSolved.resize(0);
        return true;
    }

    double h = params_.timeStep;
    Eigen::SparseMatrix<double> dgDrift;
    computeConstraintJacobian(qDrift, rodIndices, dgDrift);

    Eigen::MatrixXd B = Eigen::MatrixXd(dgDrift.transpose());
    for (int i = 0; i < n; i++)
        B.row(i) *= h * h * invMass[i];

    Eigen::VectorXd base = qDrift + h * qdotOld;
    for (int i = 0; i < n; i++)
        base[i] += h * h * invMass[i] * F[i];

    // Newton on lambda only, since qDrift is already known for this step
    lambdaSolved = Eigen::VectorXd::Zero(m);
    for (int iter = 0; iter < params_.NewtonMaxIters; iter++)
    {
        Eigen::VectorXd evalQ = base + B * lambdaSolved;
        Eigen::VectorXd f;
        computeRodConstraintValues(evalQ, rodIndices, f);
        if (f.norm() < params_.NewtonTolerance)
            return true;

        Eigen::SparseMatrix<double> dgEval;
        computeConstraintJacobian(evalQ, rodIndices, dgEval);
        Eigen::MatrixXd df = Eigen::MatrixXd(dgEval * B);

        Eigen::VectorXd delta = df.colPivHouseholderQr().solve(-f);
        if (!delta.allFinite())
            return false;
        lambdaSolved += delta;
    }

    Eigen::VectorXd finalEval = base + B * lambdaSolved;
    Eigen::VectorXd finalF;
    computeRodConstraintValues(finalEval, rodIndices, finalF);
    return finalF.norm() < params_.NewtonTolerance;
}


void initSimulation()
{
    time_ = 0;
    particles_.clear();
    for (std::vector<Connector*>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
        delete* it;
    connectors_.clear();
    saws_.clear();
    bendingStencils_.clear();
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
            vertexColors.push_back(Eigen::Vector3d(0.78, 0.92, 0.78));
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
        switch ((*it)->getType())
        {
        case SimParameters::CT_SPRING:
        {
            Spring& spring = *(Spring*)*it;
            Eigen::Vector3d color;
            if ((*it)->associatedBendingStencils.empty())
            {
                if (spring.canSnap)
                    color << 0.16, 0.36, 0.96;
                else
                    color << 0.10, 0.62, 0.72;
            }
            else
                color << 0.96, 0.62, 0.32;
            Eigen::Vector2d sourcepos = particles_[(*it)->p1].pos;
            Eigen::Vector2d destpos = particles_[(*it)->p2].pos;

            Eigen::Vector2d vec = destpos - sourcepos;
            double dist = vec.norm();
            if (dist < 1e-12)
                continue;
            Eigen::Vector2d perp(-vec[1], vec[0]);
            perp /= dist;

            double styleScale = spring.canSnap ? 0.90 : 1.15;
            double width = styleScale * baselinewidth / (1.0 + 16.0 * dist * dist);

            for (int i = 0; i < 4; i++)
                vertexColors.push_back(color);

            verts.push_back(Eigen::Vector3d(sourcepos[0] + width * perp[0], sourcepos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(sourcepos[0] - width * perp[0], sourcepos[1] - width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] + width * perp[0], destpos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] - width * perp[0], destpos[1] - width * perp[1], -eps));

            faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));
            faces.push_back(Eigen::Vector3i(idx + 2, idx + 1, idx + 3));
            idx += 4;

            break;
        }
        case SimParameters::CT_RIGIDROD:
        {
            Eigen::Vector3d color;
            if ((*it)->associatedBendingStencils.empty())
                color << 0.96, 0.40, 0.78;
            else
                color << 0.98, 0.90, 0.30;

            Eigen::Vector2d sourcepos = particles_[(*it)->p1].pos;
            Eigen::Vector2d destpos = particles_[(*it)->p2].pos;
            Eigen::Vector2d vec = destpos - sourcepos;
            double dist = vec.norm();
            if (dist < 1e-12)
                continue;
            Eigen::Vector2d perp(-vec[1], vec[0]);
            perp /= dist;

            double width = 1.45 * baselinewidth;

            for (int i = 0; i < 4; i++)
                vertexColors.push_back(color);

            verts.push_back(Eigen::Vector3d(sourcepos[0] + width * perp[0], sourcepos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(sourcepos[0] - width * perp[0], sourcepos[1] - width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] + width * perp[0], destpos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] - width * perp[0], destpos[1] - width * perp[1], -eps));

            faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));
            faces.push_back(Eigen::Vector3i(idx + 2, idx + 1, idx + 3));
            idx += 4;

            break;
        }
        default:
            break;
        }
    }

    int nparticles = particles_.size();

    for (int i = 0; i < nparticles; i++)
    {
        double totalMass = getTotalParticleMass(i);
        double renderMass = std::max(0.05, std::min(8.0, totalMass));
        double radius = baseradius * sqrt(renderMass);
        radius *= (1.0 + pulsefactor * sin(pulsespeed * time_));

        Eigen::Vector3d color(0, 0, 0);

        if (particles_[i].fixed)
        {
            radius = baseradius;
            color << 0.96, 0.22, 0.22;
        }
        else if (particles_[i].inert)
        {
            radius = 0.75 * baseradius;
            color << 0.20, 0.42, 0.42;
        }
        else
        {
            double warm = std::max(0.0, std::min(1.0, (renderMass - 0.5) / 3.0));
            color << 0.05 + 0.58 * warm, 0.08 + 0.36 * warm, 0.14;
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

    int numparticles = particles_.size() - 1;

    for (int i = 0; i < numparticles; i++)
    {
        if (particles_[i].inert)
            continue;
        Eigen::Vector2d pos = particles_[i].pos;
        double dist = (pos - newpos).norm();
        if (dist <= params_.maxSpringDist && dist > 1e-12)
        {
            if (params_.connectorType == SimParameters::CT_SPRING)
            {
                connectors_.push_back(new Spring(newid, i, 0, params_.springStiffness / dist, dist, true));
            }
            else if (params_.connectorType == SimParameters::CT_RIGIDROD)
            {
                connectors_.push_back(new RigidRod(newid, i, 0, dist));
            }
            else if (params_.connectorType == SimParameters::CT_FLEXROD)
            {
                createFlexibleRod(newid, i);
            }
        }
    }
    
}

void addSaw(double x, double y)
{
    saws_.push_back(Saw(Eigen::Vector2d(x, y), params_.sawRadius));
}


void buildConfiguration(Eigen::VectorXd& q, Eigen::VectorXd& lambda, Eigen::VectorXd& qdot)
{
    int ndofs = 2 * particles_.size();
    q.resize(ndofs);
    qdot.resize(ndofs);

    for (int i = 0; i < (int)particles_.size(); i++)
    {
        q.segment<2>(2 * i) = particles_[i].pos;
        qdot.segment<2>(2 * i) = particles_[i].vel;
    }

    std::vector<int> rodIndices = getRigidRodConnectorIndices();
    lambda.resize(rodIndices.size());
    for (int i = 0; i < (int)rodIndices.size(); i++)
    {
        lambda[i] = ((RigidRod*)connectors_[rodIndices[i]])->lambda;
    }
}

void unbuildConfiguration(const Eigen::VectorXd& q, const Eigen::VectorXd &lambda, const Eigen::VectorXd& qdot)
{
    int ndofs = q.size();
    assert(ndofs == int(2 * particles_.size()));

    for (int i = 0; i < ndofs / 2; i++)
    {
        particles_[i].pos = q.segment<2>(2 * i);
        particles_[i].vel = qdot.segment<2>(2 * i);
    }

    std::vector<int> rodIndices = getRigidRodConnectorIndices();
    assert(lambda.size() == (int)rodIndices.size());
    for (int i = 0; i < (int)rodIndices.size(); i++)
    {
        ((RigidRod*)connectors_[rodIndices[i]])->lambda = lambda[i];
    }
}

void computeMassInverse(Eigen::SparseMatrix<double>& Minv)
{
    int ndofs = 2 * int(particles_.size());

    Minv.resize(ndofs, ndofs);
    Minv.setZero();

    std::vector<Eigen::Triplet<double> > Minvcoeffs;
    for (int i = 0; i < ndofs / 2; i++)
    {
        double invm = getInverseMassForParticle(i);
        Minvcoeffs.push_back(Eigen::Triplet<double>(2 * i, 2 * i, invm));
        Minvcoeffs.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, invm));
    }

    Minv.setFromTriplets(Minvcoeffs.begin(), Minvcoeffs.end());
}

void processGravityForce(Eigen::VectorXd& F)
{
    int nparticles = (int)particles_.size();
    for (int i = 0; i < nparticles; i++)
    {
        if (!particles_[i].fixed)
        {
            F[2 * i + 1] += params_.gravityG * getTotalParticleMass(i);
        }
    }
}

void processSpringForce(const Eigen::VectorXd& q, Eigen::VectorXd& F, std::vector<Eigen::Triplet<double> >& H)
{
    int nsprings = (int)connectors_.size();

    for (int i = 0; i < nsprings; i++)
    {
        if (connectors_[i]->getType() != SimParameters::CT_SPRING)
            continue;
        Spring& s = *(Spring*)connectors_[i];
        Eigen::Vector2d p1 = q.segment<2>(2 * s.p1);
        Eigen::Vector2d p2 = q.segment<2>(2 * s.p2);
        double dist = (p2 - p1).norm();
        if (dist < 1e-12)
            continue;
        Eigen::Vector2d localF = s.stiffness * (dist - s.restlen) / dist * (p2 - p1);
        F.segment<2>(2 * s.p1) += localF;
        F.segment<2>(2 * s.p2) -= localF;

        Eigen::Matrix2d I;
        I << 1, 0, 0, 1;
        Eigen::Matrix2d localH = s.stiffness * (1.0 - s.restlen / dist) * I;
        localH += s.stiffness * s.restlen * (p2 - p1) * (p2 - p1).transpose() / dist / dist / dist;

        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
            {
                H.push_back(Eigen::Triplet<double>(2 * s.p1 + j, 2 * s.p1 + k, localH.coeff(j, k)));
                H.push_back(Eigen::Triplet<double>(2 * s.p2 + j, 2 * s.p2 + k, localH.coeff(j, k)));
                H.push_back(Eigen::Triplet<double>(2 * s.p1 + j, 2 * s.p2 + k, -localH.coeff(j, k)));
                H.push_back(Eigen::Triplet<double>(2 * s.p2 + j, 2 * s.p1 + k, -localH.coeff(j, k)));
            }
    }
}

void processDampingForce(const Eigen::VectorXd& q, const Eigen::VectorXd& qprev, Eigen::VectorXd& F, std::vector<Eigen::Triplet<double> >& H)
{
    int nsprings = (int)connectors_.size();

    for (int i = 0; i < nsprings; i++)
    {
        if (connectors_[i]->getType() != SimParameters::CT_SPRING)
            continue;
        Spring& s = *(Spring*)connectors_[i];
        Eigen::Vector2d p1 = q.segment<2>(2 * s.p1);
        Eigen::Vector2d p2 = q.segment<2>(2 * s.p2);
        Eigen::Vector2d p1prev = qprev.segment<2>(2 * s.p1);
        Eigen::Vector2d p2prev = qprev.segment<2>(2 * s.p2);

        Eigen::Vector2d relvel = (p2 - p2prev) / params_.timeStep - (p1 - p1prev) / params_.timeStep;
        Eigen::Vector2d localF = params_.dampingStiffness * relvel;
        F.segment<2>(2 * s.p1) += localF;
        F.segment<2>(2 * s.p2) -= localF;

        Eigen::Matrix2d I;
        I << 1, 0, 0, 1;
        Eigen::Matrix2d localH = params_.dampingStiffness * I / params_.timeStep;

        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
            {
                H.push_back(Eigen::Triplet<double>(2 * s.p1 + j, 2 * s.p1 + k, localH.coeff(j, k)));
                H.push_back(Eigen::Triplet<double>(2 * s.p2 + j, 2 * s.p2 + k, localH.coeff(j, k)));
                H.push_back(Eigen::Triplet<double>(2 * s.p1 + j, 2 * s.p2 + k, -localH.coeff(j, k)));
                H.push_back(Eigen::Triplet<double>(2 * s.p2 + j, 2 * s.p1 + k, -localH.coeff(j, k)));
            }
    }
}

void processFloorForce(const Eigen::VectorXd& q, const Eigen::VectorXd& qprev, Eigen::VectorXd& F, std::vector<Eigen::Triplet<double> >& H)
{
    int nparticles = particles_.size();

    double basestiffness = 10000;
    double basedrag = 1000.0;

    for (int i = 0; i < nparticles; i++)
    {
        if (q[2 * i + 1] < -0.5 && !particles_[i].fixed)
        {
            double vel = (q[2 * i + 1] - qprev[2 * i + 1]) / params_.timeStep;
            double dist = -0.5 - q[2 * i + 1];

            F[2 * i + 1] += basestiffness * dist - basedrag * dist * vel;

            H.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, basestiffness
                - 0.5 * basedrag / params_.timeStep
                + basedrag * qprev[2 * i + 1] / params_.timeStep
                - 2.0 * basedrag * q[2 * i + 1] / params_.timeStep));
        }
    }
}

void processBendingForce(const Eigen::VectorXd& q, Eigen::VectorXd& F)
{
    for (int i = 0; i < (int)bendingStencils_.size(); i++)
    {
        const BendingStencil& stencil = bendingStencils_[i];
        Eigen::Vector2d pi = q.segment<2>(2 * stencil.p1);
        Eigen::Vector2d pj = q.segment<2>(2 * stencil.p2);
        Eigen::Vector2d pk = q.segment<2>(2 * stencil.p3);

        Eigen::Vector2d e1 = pj - pi;
        Eigen::Vector2d e2 = pk - pj;
        double l1 = e1.norm();
        double l2 = e2.norm();
        if (l1 < 1e-12 || l2 < 1e-12)
            continue;

        double cross = e1[0] * e2[1] - e1[1] * e2[0];
        double dot = e1.dot(e2);
        double theta = 2.0 * atan2(cross, l1 * l2 + dot);

        Eigen::Vector2d rot1(e1[1], -e1[0]);
        Eigen::Vector2d rot2(e2[1], -e2[0]);
        Eigen::Vector2d Fi = stencil.kb * theta * rot1 / (l1 * l1);
        Eigen::Vector2d Fk = stencil.kb * theta * rot2 / (l2 * l2);
        Eigen::Vector2d Fj = -Fi - Fk;

        F.segment<2>(2 * stencil.p1) += Fi;
        F.segment<2>(2 * stencil.p2) += Fj;
        F.segment<2>(2 * stencil.p3) += Fk;
    }
}

void computeForceAndHessian(const Eigen::VectorXd& q, const Eigen::VectorXd& qprev, Eigen::VectorXd& F, Eigen::SparseMatrix<double>& H)
{
    F.resize(q.size());
    F.setZero();
    H.resize(q.size(), q.size());
    H.setZero();

    std::vector<Eigen::Triplet<double> > Hcoeffs;
    if (params_.gravityEnabled)
        processGravityForce(F);
    if (params_.springsEnabled)
        processSpringForce(q, F, Hcoeffs);
    if (params_.dampingEnabled)
        processDampingForce(q, qprev, F, Hcoeffs);
    if (params_.floorEnabled)
        processFloorForce(q, qprev, F, Hcoeffs);
    if (params_.bendingEnabled)
        processBendingForce(q, F);
    processPenaltyForce(q, F, Hcoeffs);

    H.setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());
}


void numericalIntegration(Eigen::VectorXd& q, Eigen::VectorXd& lambda, Eigen::VectorXd& qdot)
{
    Eigen::VectorXd F;
    Eigen::SparseMatrix<double> H;
    Eigen::SparseMatrix<double> Minv;
    Eigen::VectorXd invMass;

    computeMassInverse(Minv);
    buildInverseMassVector(invMass);

    const double h = params_.timeStep;
    Eigen::VectorXd oldq = q;
    Eigen::VectorXd oldqdot = qdot;

    if (params_.constraintHandling == SimParameters::CH_LAGRANGEMULT)
    {
        q = oldq + h * oldqdot;
        computeForceAndHessian(q, oldq, F, H);

        std::vector<int> rodIndices = getRigidRodConnectorIndices();
        if (rodIndices.empty())
        {
            qdot = oldqdot + h * Minv * F;
            lambda.resize(0);
            return;
        }

        Eigen::VectorXd momentum = oldqdot;
        for (int i = 0; i < momentum.size(); i++)
        {
            if (invMass[i] > 0.0)
                momentum[i] /= invMass[i];
            else
                momentum[i] = 0.0;
        }

        Eigen::VectorXd solvedLambda;
        bool converged = solveConstrainedLagrangeLambda(q, oldqdot, invMass, F, solvedLambda);
        if (!converged)
            solvedLambda = Eigen::VectorXd::Zero(rodIndices.size());
        Eigen::SparseMatrix<double> dg;
        computeConstraintJacobian(q, rodIndices, dg);
        momentum += h * F + h * dg.transpose() * solvedLambda;

        qdot.resize(momentum.size());
        for (int i = 0; i < momentum.size(); i++)
            qdot[i] = invMass[i] * momentum[i];

        lambda = solvedLambda;
        return;
    }

    q = oldq + h * oldqdot;
    computeForceAndHessian(q, oldq, F, H);
    qdot = oldqdot + h * Minv * F;

    if (params_.constraintHandling == SimParameters::CH_STEPPROJECT)
    {
        Eigen::VectorXd qTilde = q;
        Eigen::VectorXd qProjected, lambdaProjected;
        bool converged = projectWithNewton(qTilde, invMass, qProjected, lambdaProjected);
        if (converged)
        {
            q = qProjected;
            lambda = lambdaProjected;
            qdot += (qProjected - qTilde) / h;
        }
        else
        {
            lambda = Eigen::VectorXd::Zero(getRigidRodConnectorIndices().size());
        }
    }
}

double ptSegmentDist(const Eigen::Vector2d& p, const Eigen::Vector2d& q1, const Eigen::Vector2d& q2)
{
    Eigen::Vector2d seg = q2 - q1;
    double seglen2 = seg.dot(seg);
    if (seglen2 < 1e-12)
        return (p - q1).norm();
    double t = (p - q1).dot(seg) / seglen2;
    double linedistsq = (q1 + t * (q2 - q1) - p).squaredNorm();
    double q1dist = (p - q1).squaredNorm();
    double q2dist = (p - q2).squaredNorm();
    double mindistsq = std::min(linedistsq, std::min(q1dist, q2dist));
    return sqrt(mindistsq);
}

void detectSawedConnectors(std::set<int>& connectorsToDelete)
{
    for (int i = 0; i < (int)connectors_.size(); i++)
    {
        Eigen::Vector2d pos1 = particles_[connectors_[i]->p1].pos;
        Eigen::Vector2d pos2 = particles_[connectors_[i]->p2].pos;
        double maxx = std::max(pos1[0], pos2[0]);
        double minx = std::min(pos1[0], pos2[0]);
        double maxy = std::max(pos1[1], pos2[1]);
        double miny = std::min(pos1[1], pos2[1]);
        for (std::vector<Saw>::iterator saw = saws_.begin(); saw != saws_.end(); ++saw)
        {
            Eigen::Vector2d sawpos = saw->pos;
            double sawr = saw->radius;

            if (sawpos[0] - sawr > maxx || sawpos[0] + sawr < minx || sawpos[1] - sawr > maxy || sawpos[1] + sawr < miny)
                continue;

            double sawspringdist = ptSegmentDist(sawpos, pos1, pos2);
            if (sawspringdist <= sawr)
            {
                connectorsToDelete.insert(i);
                break;
            }
        }
    }
}

void detectSawedParticles(std::set<int>& particlesToDelete)
{
    for (int i = 0; i < (int)particles_.size(); i++)
    {
        Eigen::Vector2d partpos = particles_[i].pos;

        if (!(fabs(partpos[0]) < 2 && fabs(partpos[1]) < 2))
        {
            particlesToDelete.insert(i);
            continue;
        }

        for (std::vector<Saw>::iterator it = saws_.begin(); it != saws_.end(); ++it)
        {
            Eigen::Vector2d sawpos = it->pos;
            double sqdist = (sawpos - partpos).squaredNorm();
            if (sqdist < it->radius * it->radius)
            {
                particlesToDelete.insert(i);
                break;
            }
        }
    }
}

void deleteSawedObjects()
{
    std::set<int> particlestodelete;
    std::set<int> connectorstodelete;
    detectSawedParticles(particlestodelete);
    detectSawedConnectors(connectorstodelete);

    if (!particlestodelete.empty())
    {
        for (int i = 0; i < (int)connectors_.size(); i++)
        {
            if (particlestodelete.count(connectors_[i]->p1) || particlestodelete.count(connectors_[i]->p2))
                connectorstodelete.insert(i);
        }
    }

    if (particlestodelete.empty() && connectorstodelete.empty())
        return;

    std::vector<Particle, Eigen::aligned_allocator<Particle>> newparticles;
    std::vector<Connector*> newconnectors;
    std::vector<int> remainingparticlemap(particles_.size(), -1);

    for (int i = 0; i < (int)particles_.size(); i++)
    {
        if (particlestodelete.count(i))
            continue;
        remainingparticlemap[i] = (int)newparticles.size();
        newparticles.push_back(particles_[i]);
    }

    for (int i = 0; i < (int)connectors_.size(); i++)
    {
        Connector* conn = connectors_[i];
        bool remove = connectorstodelete.count(i) > 0;
        if (!remove)
        {
            if (conn->p1 < 0 || conn->p2 < 0 || conn->p1 >= (int)remainingparticlemap.size() || conn->p2 >= (int)remainingparticlemap.size())
                remove = true;
            else if (remainingparticlemap[conn->p1] < 0 || remainingparticlemap[conn->p2] < 0)
                remove = true;
        }

        if (remove)
        {
            delete conn;
            continue;
        }

        conn->p1 = remainingparticlemap[conn->p1];
        conn->p2 = remainingparticlemap[conn->p2];
        newconnectors.push_back(conn);
    }

    particles_ = newparticles;
    connectors_ = newconnectors;
    rebuildBendingStencilAssociations();
}

void pruneOverstrainedSprings()
{   
    int nsprings = connectors_.size();

    std::vector<int> toremove;
    for (int i = 0; i < nsprings; i++)
    {
        if (connectors_[i]->getType() != SimParameters::CT_SPRING)
            continue;
        Spring& s = *(Spring*)connectors_[i];
        if (s.canSnap)
        {
            Eigen::Vector2d srcpos = particles_[s.p1].pos;
            Eigen::Vector2d dstpos = particles_[s.p2].pos;
            double dist = (dstpos - srcpos).norm();
            if (s.restlen <= 1e-12)
                continue;

            double strain = (dist - s.restlen) / s.restlen;
            if (strain > params_.maxSpringStrain)
                toremove.push_back(i);
        }
    }

    for (std::vector<int>::reverse_iterator it = toremove.rbegin(); it != toremove.rend(); ++it)
    {
        delete connectors_[*it];
        connectors_.erase(connectors_.begin() + *it);
    }

    if (!toremove.empty())
        rebuildBendingStencilAssociations();
}

bool simulateOneStep()
{
    // Create configurational vectors
    Eigen::VectorXd q, lambda, qdot;
    buildConfiguration(q, lambda, qdot);
    // Use them for one step of time integration
    numericalIntegration(q, lambda, qdot);
    // Unpack the DOFs back into the particles for rendering
    unbuildConfiguration(q, lambda, qdot);

    // Cleanup: delete sawed objects and snapped springs
    pruneOverstrainedSprings();
    deleteSawedObjects();
    
    // Time advances
    time_ += params_.timeStep;
    return false;
}

bool isBendingStencilStateValid()
{
    for (int i = 0; i < (int)bendingStencils_.size(); i++)
    {
        const BendingStencil& stencil = bendingStencils_[i];
        int springA = findSpringConnectorByParticles(connectors_, stencil.p1, stencil.p2);
        int springB = findSpringConnectorByParticles(connectors_, stencil.p2, stencil.p3);
        if (springA < 0 || springB < 0 || springA == springB)
            return false;
        if (connectors_[springA]->associatedBendingStencils.count(i) == 0)
            return false;
        if (connectors_[springB]->associatedBendingStencils.count(i) == 0)
            return false;
    }
    return true;
}

bool testFlexibleRodCreation()
{
    params_ = SimParameters();
    initSimulation();

    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 5;
    params_.rodDensity = 1.5;
    params_.rodStretchingStiffness = 120.0;
    params_.rodBendingStiffness = 0.1;

    addParticle(0.0, 0.0);
    addParticle(0.6, 0.0);

    int s = std::max(2, params_.rodSegments);
    if ((int)particles_.size() != 2 + (s - 1))
        return false;
    if ((int)connectors_.size() != s)
        return false;
    if ((int)bendingStencils_.size() != s - 1)
        return false;

    for (int i = 2; i < (int)particles_.size(); i++)
    {
        if (!particles_[i].inert || particles_[i].fixed || std::abs(particles_[i].mass) > 1e-12)
            return false;
    }

    for (int i = 0; i < (int)connectors_.size(); i++)
    {
        if (connectors_[i]->getType() != SimParameters::CT_SPRING)
            return false;
        Spring& spring = *(Spring*)connectors_[i];
        if (spring.canSnap)
            return false;
        if (spring.mass <= 0.0)
            return false;
    }

    return isBendingStencilStateValid();
}

bool testStepAndProjectConstraint()
{
    params_ = SimParameters();
    initSimulation();

    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.constraintHandling = SimParameters::CH_STEPPROJECT;
    params_.maxSpringDist = 2.0;
    params_.gravityEnabled = false;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.bendingEnabled = false;
    params_.timeStep = 0.01;

    addParticle(-0.2, 0.0);
    addParticle(0.2, 0.0);
    if (connectors_.size() != 1 || connectors_[0]->getType() != SimParameters::CT_RIGIDROD)
        return false;

    particles_[0].vel << -1.0, 0.15;
    particles_[1].vel << 1.0, -0.15;

    simulateOneStep();

    double dist = (particles_[0].pos - particles_[1].pos).norm();
    double rest = ((RigidRod*)connectors_[0])->length;
    return std::abs(dist - rest) < 1e-5;
}

bool testSawDeletionPrunesBendingStencils()
{
    params_ = SimParameters();
    initSimulation();

    params_.connectorType = SimParameters::CT_FLEXROD;
    params_.maxSpringDist = 2.0;
    params_.rodSegments = 3;
    params_.sawRadius = 0.035;

    addParticle(0.0, 0.0);
    addParticle(0.6, 0.0);
    if (connectors_.size() != 3 || bendingStencils_.size() != 2)
        return false;

    addSaw(0.3, 0.0);
    deleteSawedObjects();

    if (!bendingStencils_.empty())
        return false;
    if (connectors_.size() != 2)
        return false;
    for (int i = 0; i < (int)connectors_.size(); i++)
    {
        if (connectors_[i]->getType() != SimParameters::CT_SPRING)
            return false;
        if (!connectors_[i]->associatedBendingStencils.empty())
            return false;
    }
    return true;
}

bool testLagrangeMultiplierSmoke()
{
    params_ = SimParameters();
    initSimulation();

    params_.connectorType = SimParameters::CT_RIGIDROD;
    params_.constraintHandling = SimParameters::CH_LAGRANGEMULT;
    params_.maxSpringDist = 2.0;
    params_.gravityEnabled = false;
    params_.springsEnabled = false;
    params_.dampingEnabled = false;
    params_.floorEnabled = false;
    params_.bendingEnabled = false;
    params_.timeStep = 0.01;

    addParticle(-0.25, 0.0);
    addParticle(0.25, 0.0);
    particles_[0].vel << -0.8, 0.0;
    particles_[1].vel << 0.8, 0.0;

    simulateOneStep();

    for (int i = 0; i < (int)particles_.size(); i++)
    {
        if (!particles_[i].pos.allFinite() || !particles_[i].vel.allFinite())
            return false;
    }
    return true;
}

bool runInternalTests()
{
    struct NamedTest
    {
        std::string name;
        bool (*testFn)();
    };

    std::vector<NamedTest> tests;
    tests.push_back({ "Flexible rod creation", testFlexibleRodCreation });
    tests.push_back({ "Step-and-project rod constraint", testStepAndProjectConstraint });
    tests.push_back({ "Saw cleanup of rope stencils", testSawDeletionPrunesBendingStencils });
    tests.push_back({ "Lagrange multiplier smoke test", testLagrangeMultiplierSmoke });

    bool allPassed = true;
    for (int i = 0; i < (int)tests.size(); i++)
    {
        bool passed = tests[i].testFn();
        std::cout << (passed ? "[PASS] " : "[FAIL] ") << tests[i].name << std::endl;
        allPassed = allPassed && passed;
    }
    if (allPassed)
        std::cout << "All tests passed." << std::endl;
    else
        std::cout << "At least one test failed." << std::endl;
    return allPassed;
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
        ImGui::Combo("Connector Type", (int*)&params_.connectorType, "Springs\0Rigid Rods\0Flexible Rods\0\0");
    }
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Combo("Constraint Handling", (int*)&params_.constraintHandling, "Penalty Method\0Step and Project\0Lagrange Multipliers\0\0");

        ImGui::InputDouble("Timestep", &params_.timeStep);
        ImGui::InputDouble("Newton Tolerance", &params_.NewtonTolerance);
        ImGui::InputInt("Newton Max Iters", &params_.NewtonMaxIters);
        ImGui::InputDouble("Penalty Stiffness", &params_.penaltyStiffness);
    }
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
        ImGui::InputDouble("  Gravity g", &params_.gravityG);
        ImGui::Checkbox("Springs Enabled", &params_.springsEnabled);
        ImGui::InputDouble("  Max Strain", &params_.maxSpringStrain);
        ImGui::Checkbox("Damping Enabled", &params_.dampingEnabled);
        ImGui::InputDouble("  Viscosity", &params_.dampingStiffness);
        ImGui::Checkbox("Floor Enabled", &params_.floorEnabled);
        ImGui::Checkbox("Bending Enabled", &params_.bendingEnabled);
        //viewer.imgui->addWindow(Eigen::Vector2i(1000, 0), "New Objects");
    }


    if (ImGui::CollapsingHeader("New Particles", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Is Fixed", &params_.particleFixed);
        ImGui::InputDouble("Mass", &params_.particleMass);
    }

    if (ImGui::CollapsingHeader("New Saws", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputDouble("Radius", &params_.sawRadius);
    }

    if (ImGui::CollapsingHeader("New Springs", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputDouble("Max Spring Dist", &params_.maxSpringDist);
        ImGui::InputDouble("Base Stiffness", &params_.springStiffness);
    }


    if (ImGui::CollapsingHeader("New Rods", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputInt("Num Segments", &params_.rodSegments);
        ImGui::InputDouble("Density", &params_.rodDensity);
        ImGui::InputDouble("Stretching Stiffness", &params_.rodStretchingStiffness);
        ImGui::InputDouble("Bending Stiffness", &params_.rodBendingStiffness);
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
  if (argc > 1 && std::string(argv[1]) == "--run-tests")
  {
      return runInternalTests() ? 0 : 1;
  }

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

