
double computeDeviation(const Shreds& p1, const Shreds& p2)
{
    double dev = 0;
    for (int i = 0; i < N; ++i)
    {
        dev += (p1[i].r - p2[i].r).len();
    }
#if INTEGRATOR == 1
    return dev/15;
#elif INTEGRATOR == 2
    return dev;
#else
#error unknown POTENTIAL
#endif//POTENTIAL
}

#if INTEGRATOR == 1
Shreds tryMove(const Shreds& particles, double dt)
{
    Shreds k[4];
    for (int i = 0; i < 4; ++i)
        k[i].resize(N);
    Shreds ret;
    ret.resize(N);
    for (int i = 0; i < N; ++i)
    {
        k[0][i] = particles[i];
        k[0][i].a = computeForce(i, particles);

        ret[i].r = particles[i].r + k[0][i].v*(dt/2);
        ret[i].v = particles[i].v + k[0][i].a*(dt/2);
        ret[i].a = {0., 0.};
    }
    for (int i = 0; i < N; ++i)
    {
        k[1][i] = ret[i];
        k[1][i].a = computeForce(i, ret);

        ret[i].r = particles[i].r + k[1][i].v*(dt/2);
        ret[i].v = particles[i].v + k[1][i].a*(dt/2);
        ret[i].a = {0., 0.};
    }
    for (int i = 0; i < N; ++i)
    {
        k[2][i] = ret[i];
        k[2][i].a = computeForce(i, ret);

        ret[i].r = particles[i].r + k[2][i].v*dt;
        ret[i].v = particles[i].v + k[2][i].a*dt;
        ret[i].a = {0., 0.};
    }
    for (int i = 0; i < N; ++i)
    {
        k[3][i] = ret[i];
        k[3][i].a = computeForce(i, ret);

        ret[i].r = particles[i].r + (k[0][i].v + k[1][i].v*2 +
                                     k[2][i].v*2 + k[3][i].v)*(dt/6);
        ret[i].v = particles[i].v + (k[0][i].a + k[1][i].a*2 +
                                     k[2][i].a*2 + k[3][i].a)*(dt/6);
        ret[i].a = {0., 0.};
    }
    return ret;
}
#elif INTEGRATOR == 2
Shreds tryMove(const Shreds& particles, double dt)
{
    Shreds ret;
    ret.resize(N);
    for (int i = 0; i < N; ++i)
    {
        ret[i].a = computeForce(i, ret);

        ret[i].r = particles[i].r + particles[i].v*dt
            + particles[i].a*(dt*dt/2);
        ret[i].v = particles[i].v + (particles[i].a + ret[i].a)*(dt/2);
    }
    return ret;
}
#else
#error Integrator != 1,2 is undefined
#endif//INTEGRATOR

void KeepInBox(Shreds* parts)
{
    if (keepInBox)
    {
        for (int i = 0; i < N; ++i)
        {
            if ((*parts)[i].r.x < 0)
            {
                (*parts)[i].r.x *= -1;
                if ((*parts)[i].v.x < 0)
                    (*parts)[i].v.x *= -1;
            }
            if ((*parts)[i].r.y < 0)
            {
                (*parts)[i].r.y *= -1;
                if ((*parts)[i].v.y < 0)
                    (*parts)[i].v.y *= -1;
            }
            if ((*parts)[i].r.x > box.x)
            {
                (*parts)[i].r.x = 2*box.x - (*parts)[i].r.x;
                if ((*parts)[i].v.x > 0)
                    (*parts)[i].v.x *= -1;
            }
            if ((*parts)[i].r.y > box.y)
            {
                (*parts)[i].r.y = 2*box.y - (*parts)[i].r.y;
                if ((*parts)[i].v.y > 0)
                    (*parts)[i].v.y *= -1;
            }
        }
    }
}

Scene moveParticles(const Scene& scene, double step, double* h)
{
    Shreds p1 = tryMove(tryMove(scene.particles, *h), *h);
    Shreds p2 = tryMove(scene.particles, *h*2);
    double dt = *h;
    while (dt < step)
    {
        p1 = tryMove(tryMove(p1, *h), *h);
        p2 = tryMove(p2, *h*2);
        dt += *h*2;
    }
    double dev = computeDeviation(p1, p2);
    if (dev < maxDev)
    {
        if (*h*2 <= step)
            *h = *h*2;
        KeepInBox(&p2);
        return Scene({p2, scene.time + dt, *h, dev});
    }
    cout <<"selecting step:";
    while (dev > maxDev)
    {
        *h = *h/2;
        cout <<*h <<" ";
        cout.flush();
        p2 = p1;
        p1 = tryMove(scene.particles, *h);
        dt = *h;
        while (dt < step)
        {
            p1 = tryMove(p1, *h);
            dt += *h;
        }
        dev = computeDeviation(p1, p2);
        if (*h < epsilon)
        {
            cerr <<"Insuperable diversion" <<endl;
            std::raise(SIGTERM);
        }
    }
    cout <<endl;
    KeepInBox(&p1);
    return Scene({p1, scene.time + dt, *h, dev});
}
