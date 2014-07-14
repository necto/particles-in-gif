#include "common.h"
#include "potential.h"

double A, rmin, rmax;

Point computeForceComponent(Point i, Point j)
{
    Point delta = j - i;
    double distance = delta.len();
    if (distance < epsilon) return {0.,0.};
    return delta*(A*(2*std::log(distance) - 1)) +
        delta*(-std::exp(-(distance - rmin)/(rmax-rmin))/distance/(rmax-rmin));
}

Point computeForce(int i, const Shreds& particles)
{
    Point a = {0., 0.};
    for (int j = 0; j < N; ++j)
        if (j != i)
        {
            a += computeForceComponent(particles[i].r, particles[j].r);
        }
    return a;
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
}
