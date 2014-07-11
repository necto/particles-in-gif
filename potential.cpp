
#if POTENTIAL == 1
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
double computeEnergy(const Shreds& particles)
{
    return {0., 0.};
}

#elif POTENTIAL == 2
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
double computeEnergy(const Shreds& particles)
{
    return {0., 0.};
}

#elif POTENTIAL == 3
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

#else
#error unknown potential (only #1 and #2 are coded)
#endif
