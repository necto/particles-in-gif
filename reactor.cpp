#include <Magick++.h>
#include <Magick++/STL.h>
#include <sstream>
#include <list>
#include "common.h"

using Magick::InitializeMagick;
using Magick::Geometry;
using Magick::Image;
using Magick::writeImages;
using Magick::DrawableCircle;
using Magick::DrawableText;
using Magick::DrawableFont;
using Magick::DrawableLine;
using Magick::DrawablePointSize;
using std::ostringstream;
using std::ifstream;
using std::list;

#ifndef POTENTIAL
#define POTENTIAL 3
#endif

#ifndef INTEGRATOR
#define INTEGRATOR 2
#endif

int N;
double A, rmin, rmax, T, a, b, r0;
double eps, sigma, sigma6, m;
Point G;
double maxDev;
int screenWidth;
int delay, endDelay;
double vscalefactor;

Point box;
bool keepInBox;

const double epsilon = 1e-10;

struct Scene
{
    Shreds particles;
    double time;
    double h;
    double deviation;
    double U;
    double K;
};

void readConfig(Config* cfg, const char* fname)
{
    try
    {
        cfg->readFile(fname);
    }
    catch(FileIOException &fioex)
    {
        std::cerr << "I/O error while reading file" <<fname <<"." << std::endl;
        std::raise(SIGTERM);
    }
    catch(ParseException &pex)
    {
        std::cerr << "Parse error at " << pex.getLine()
                  << " - " << pex.getError() << std::endl;
        std::raise(SIGTERM);
    }
}

Shreds initShredsDual(const string& fnamePrefix)
{
    string rfname = fnamePrefix + dataRPostfix;
    string vfname = fnamePrefix + dataVPostifx;
    ifstream r(rfname);
    ifstream v(vfname);
    if (r.fail())
    {
        cerr << "I/O error while opening file " <<rfname <<endl;
        std::raise(SIGTERM);
    }
    if (v.fail())
    {
        cerr << "I/O error while opening file " <<vfname <<endl;
        std::raise(SIGTERM);
    }

    Shreds ret;
    ret.resize(N);
    for (int i = 0; i < N; ++i)
    {
        if (v.eof() || r.eof())
        {
            cerr << "Not enough data (coors/vels) supplied" <<endl;
            std::raise(SIGTERM);
        }
        Particle p;
        r >>p.r.x >>p.r.y;
        v >>p.v.x >>p.v.y;
        p.a = {0., 0.};
        ret[i] = p;
    }


    if (r.fail())
    {
        cerr << "I/O error while reading coordinates from file " <<rfname <<endl;
        std::raise(SIGTERM);
    }
    if (v.fail())
    {
        cerr << "I/O error while reading velocities from  file " <<vfname <<endl;
        std::raise(SIGTERM);
    }
    r.close();
    v.close();
    return ret;
}

Shreds initShreds(const Config& cfg)
{
    string dataFilePrefix = getProperty<string>("dataFilePrefix", cfg, "");
    if (dataFilePrefix.empty())
    {
        const Setting &ps = getProperty<const Setting&>("particles", cfg);
        Shreds ret;
        for (int i = 0; i < ps.getLength(); ++i)
        {
            const Setting& r = ps[i]["r"];
            const Setting& v = ps[i]["v"];
            ret.push_back({{r["x"], r["y"]},{v["x"], v["y"]}});
        }
        return ret;
    }
    return initShredsDual(dataFilePrefix);
}

struct Box
{
    Point leftup;
    Point rightdown;

    inline double width() const
    { return rightdown.x - leftup.x; }

    inline Point limits() const
    { return rightdown - leftup; }

    inline double aspectRatio() const
    {return limits().y/limits().x;}

    inline void enlarge(double percent)
    {
        leftup = leftup - limits()*0.1;
        rightdown = rightdown + limits()*0.1;
    }
};

Box getBox(const vector<Scene>& seq)
{
    Box ret = {{seq[0].particles[0].r.x, seq[0].particles[0].r.y},
               {seq[0].particles[0].r.x,seq[0].particles[0].r.y}};
    for (auto ps = seq.begin(); ps != seq.end(); ps++)
        for (auto i = ps->particles.begin(); i != ps->particles.end(); i++)
        {
            if (i->r.x < ret.leftup.x) ret.leftup.x = i->r.x;
            if (i->r.y < ret.leftup.y) ret.leftup.y = i->r.y;
            if (i->r.x > ret.rightdown.x) ret.rightdown.x = i->r.x;
            if (i->r.y > ret.rightdown.y) ret.rightdown.y = i->r.y;
        }
    return ret;
}

void drawFrame(Image* img, const Shreds &parts, const Box& box, double scale)
{
    img->fillColor("red");
    img->strokeColor("black");
    img->strokeWidth(1);
    for (auto p = parts.begin(); p != parts.end(); p++)
    {
        Point coor = (p->r - box.leftup)*scale;
        int radius = 1;
        img->draw(DrawableCircle(coor.x - radius, coor.y - radius,
                                 coor.x + radius, coor.y + radius));

        Point vend = coor + p->v*vscalefactor*scale;
        img->draw(DrawableLine(coor.x, coor.y, vend.x, vend.y));
    }
}

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

string point2string(Point p)
{
    ostringstream oss;
    oss <<std::setprecision(3) <<"(" <<p.x <<", " <<p.y <<")";
    return oss.str();
}

string double2string(double val)
{
    ostringstream oss;
    oss <<std::setprecision(3) <<val;
    return oss.str();
}

void drawInfo(Image* img, const Box& box, double scale, double t,
              double h, double deviation, double U, double K,
              const Geometry& size)
{
    img->fillColor("black");
    img->strokeColor("black");
    img->strokeWidth(0);

    list<Magick::Drawable> text
        ({DrawablePointSize(15),
                DrawableFont("sans-serif", Magick::StyleType::NormalStyle, 300,
                             Magick::StretchType::NormalStretch),
                Magick::Drawable()});


    text.back() = DrawableText(size.width() - 80, size.height() - 20,
                               point2string(box.rightdown));
    img->draw(text);
    text.back() = DrawableText(10, 20, point2string(box.leftup));
    img->draw(text);
    text.back() = DrawableText(10, size.height() - 20,
                               point2string({box.leftup.x, box.rightdown.y}));
    img->draw(text);
    text.back() = DrawableText(10, size.height()/2,
                               double2string((box.leftup.y + box.rightdown.y)/2));
    img->draw(text);
    text.back() = DrawableText(size.width()/2, size.height() - 20,
                               double2string((box.leftup.x + box.rightdown.x)/2));
    img->draw(text);
    text.back() = DrawableText(size.width() - 110, 20,
                               "t:      " + double2string(t) + "\n"
                               "h:      " + double2string(h) + "\n"
                               "e r r: "  + double2string(deviation) + "\n"
                               "K:      " + double2string(K) + "\n"
                               "U:      " + double2string(U) + "\n"
                               "E:      " + double2string(K+U));
    img->draw(text);
}

void makeMovie(const vector<Scene>& sequence, string fname)
{
    Box box = getBox(sequence);
    box.enlarge(0.1);
    double width = screenWidth;
    Geometry size(width, box.aspectRatio()*width);
    double scale = width / box.width();

    cout <<"Rendering frames..." <<endl;
    vector<Image> v;
    int fifth = sequence.size()/5;
    for (int i = 0; i < sequence.size(); ++i)
    {
        v.push_back(Image(size, "white"));
        drawFrame(&v[i], sequence[i].particles, box, scale);
        drawInfo(&v[i], box, scale, sequence[i].time, sequence[i].h,
                 sequence[i].deviation, sequence[i].U, sequence[i].K,
                 size);
        v[i].animationDelay(delay);
        if ((i - (i/fifth)*fifth) == 0)
        {
            cout <<i/fifth*20 <<"% ";
            cout.flush();
        }
    }
    cout <<endl;
    v.back().animationDelay(endDelay);
    cout <<"Storing frames to file "<<fname<<endl;
    writeImages(v.begin(), v.end(), fname);
}

int main(int argc,char **argv)
{
    InitializeMagick(*argv);
    Config cfg;

    const char* fname = "start.cfg";
    if (argc>1) fname = argv[1];
    else
    {
        cout <<"No options supplied. Using default configuration file: " <<fname <<endl;
    }
    readConfig(&cfg, fname);
    N = getProperty<int>("N", cfg);
    A = getProperty<double>("A", cfg);
    T = getProperty<double>("T", cfg);
    rmin = getProperty<double>("rmin", cfg);
    rmax = getProperty<double>("rmax", cfg);
    maxDev = getProperty<double>("deviation", cfg);
    screenWidth = getProperty<int>("width", cfg);
    delay = getProperty<int>("delay", cfg);
    endDelay = getProperty<int>("endDelay", cfg);
    vscalefactor = getProperty<double>("vScaleFactor", cfg);
#if POTENTIAL == 2
    a = getProperty<double>("a", cfg);
    b = getProperty<double>("b", cfg);
    r0 = getProperty<double>("r0", cfg);
#elif POTENTIAL == 3
    r0 = getProperty<double>("r0", cfg);
    eps = getProperty<double>("epsilon", cfg);
    sigma = getProperty<double>("sigma", cfg);
    sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
    m = getProperty<double>("m", cfg);
    double gg = getProperty<double>("G", cfg);
    G.x = 0;
    G.y = gg;
#endif
    string outFname = getProperty<string>("rezultFile", cfg);
    keepInBox = getProperty<bool>("keepInBox", cfg, false);
    if (keepInBox)
        box = getProperty<Point>("box", cfg, Point{0., 0.});
    

    vector<Scene> sequence;
    sequence.push_back({initShreds(cfg), 0., 0., 0.});
    double step = getProperty<double>("step", cfg);
    double h = step;

    for (int i = 1; ; ++i)
    {
        sequence.push_back(moveParticles(sequence[i-1], step, &h));
        Point E = computeEnergy(sequence[i].particles);
        sequence[i].U = E.x;
        sequence[i].K = E.y;
        cout <<"computed state # " <<i <<" time: " <<sequence[i].time <<endl;
        if (sequence[i].time > T) break;
    }

    makeMovie(sequence, outFname);

    return 0;
}
