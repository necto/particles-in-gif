#include "common.h"
#include "potential.h"

double r0, eps, sigma, sigma6, m;
Point G;

double quad(double arg)
{
    double sq = arg*arg;
    return sq*sq;
}

Point computeForce(int i, const Shreds& particles)
{
    Point grad = {0., 0.};
    for (int j = 0; j < N; ++j)
    {
        Point delta = particles[i].r - particles[j].r;
        double distsq = delta.lensq();
        double dist8 = quad(distsq);

        if (distsq < r0*r0 && distsq > epsilon)
            grad += delta*(24*eps*sigma6/dist8*(2*sigma6/dist8 - 1));
    }
    return (grad - G)*(1/m);
}

Point computeEnergy(const Shreds& particles)
{
    double U = 0;
    double K = 0;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            Point delta = particles[i].r - particles[j].r;
            double distsq = delta.lensq();
            double dist6 = distsq*distsq*distsq;
            if (distsq < r0*r0)
                U += 4*eps*sigma6/dist6*(sigma6/dist6 - 1);
        }
        K += particles[i].v.lensq()*m/2;
    }
    return {U, K};
}

void readUParams(const Config& cfg)
{
    r0 = getProperty<double>("r0", cfg);
    eps = getProperty<double>("epsilon", cfg);
    sigma = getProperty<double>("sigma", cfg);
    sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
    m = getProperty<double>("m", cfg);
    double gg = getProperty<double>("G", cfg);
    G.x = 0;
    G.y = gg;
}
