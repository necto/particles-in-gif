#include "common.h"
#include "potential.h"

double A, rmin, rmax, a, b, r0;

double ro(double dist)
{
    return (rmax - dist)/(rmax - rmin)*std::exp(-b*(dist - r0)/r0);
}

double dro(double dist)
{
    return -std::exp(-b*(dist-r0)/r0)*(1+(rmax-dist)*b/r0)/(rmax-rmin);
}

double dphi(double dist)
{
    return std::exp(a*(dist-r0)/r0)*a*a*(dist-r0)/r0/r0;
}

double dF(double x)
{
    return 2*x*(1-std::log(x*x));
}

Point computeForce(int i, const Shreds& particles)
{
    double sro = 0;
    Point sdro = {0., 0.};
    Point sdphi = {0., 0.};
    for (int j = 0; j < N; ++j)
    {
        Point delta = particles[j].r - particles[i].r;
        double distance = delta.len();
        if (i != j &&
            distance > epsilon)
        {
            sro += ro(distance);
            sdro += delta*(1/distance)*dro(distance);
            sdphi += delta*(1/distance)*dphi(distance);
        }
    }
    return sdro*A*dF(sro) + sdphi;
}
Point computeEnergy(const Shreds& particles)
{
    return {0., 0.};
}
void readUParams(const Config& cfg)
{
    A = getProperty<double>("A", cfg);
    rmin = getProperty<double>("rmin", cfg);
    rmax = getProperty<double>("rmax", cfg);
    a = getProperty<double>("a", cfg);
    b = getProperty<double>("b", cfg);
    r0 = getProperty<double>("r0", cfg);
}
