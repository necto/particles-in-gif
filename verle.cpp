#include "common.h"
#include "integrator.h"
#include "potential.h"

const double deviationFactor = 1;

Shreds tryMove(const Shreds& particles, double dt)
{
    Shreds ret;
    ret.resize(N);
    for (int i = 0; i < N; ++i)
    {
        ret[i].a = computeForce(i, particles);

        ret[i].r = particles[i].r + particles[i].v*dt
            + particles[i].a*(dt*dt/2);
        ret[i].v = particles[i].v + (particles[i].a + ret[i].a)*(dt/2);
    }
    return ret;
}


