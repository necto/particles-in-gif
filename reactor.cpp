#include <Magick++.h>
#include <Magick++/STL.h>
#include "common.h"

using Magick::InitializeMagick;
using Magick::Geometry;
using Magick::Image;
using Magick::writeImages;
using Magick::DrawableCircle;
using Magick::DrawableText;

int N;
double A, rmin, rmax;

const double epsilon = 1e-10;

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

Shreds initShreds(const Config& cfg)
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

Box getBox(const vector<Shreds>& seq)
{
    Box ret = {{seq[0][0].r.x, seq[0][0].r.y},{seq[0][0].r.x,seq[0][0].r.y}};
    for (auto ps = seq.begin(); ps != seq.end(); ps++)
        for (auto i = ps->begin(); i != ps->end(); i++)
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
    }
    
}

Point computeForce(Point i, Point j)
{
    Point delta = j - i;
    double distance = delta.len();
    if (distance < epsilon) return {0.,0.};
    return delta*(A*(2*std::log(distance) - 1)) +
        delta*(-std::exp(-(distance - rmin)/(rmax-rmin))/distance/(rmax-rmin));
}

Shreds moveParticles(const Shreds& particles, double dt)
{
    Shreds k[4];
    for (int i = 0; i < 4; ++i)
        k[i].resize(N);
    Shreds ret;
    ret.resize(N);
    for (int i = 0; i < N; ++i)
    {
        k[0][i] = particles[i];
        for (int j = 0; j < N; ++j)
            if (j != i)
                k[0][i].a += computeForce(particles[i].r, particles[j].r);

        ret[i].r = particles[i].r + k[0][i].v*(dt/2);
        ret[i].v = particles[i].v + k[0][i].a*(dt/2);
        ret[i].a = {0., 0.};
    }
    for (int i = 0; i < N; ++i)
    {
        k[1][i] = ret[i];
        for (int j = 0; j < N; ++j)
            if (j != i)
                k[1][i].a += computeForce(ret[i].r, ret[j].r);

        ret[i].r = particles[i].r + k[1][i].v*(dt/2);
        ret[i].v = particles[i].v + k[1][i].a*(dt/2);
        ret[i].a = {0., 0.};
    }
    for (int i = 0; i < N; ++i)
    {
        k[2][i] = ret[i];
        for (int j = 0; j < N; ++j)
            if (j != i)
                k[2][i].a += computeForce(ret[i].r, ret[j].r);

        ret[i].r = particles[i].r + k[2][i].v*dt;
        ret[i].v = particles[i].v + k[2][i].a*dt;
        ret[i].a = {0., 0.};
    }
    for (int i = 0; i < N; ++i)
    {
        k[3][i] = ret[i];
        for (int j = 0; j < N; ++j)
            if (j != i)
                k[3][i].a += computeForce(ret[i].r, ret[j].r);

        ret[i].r = particles[i].r + (k[0][i].v + k[1][i].v*2 + k[2][i].v*2 + k[3][i].v)*(dt/6);
        ret[i].v = particles[i].v + (k[0][i].a + k[1][i].a*2 + k[2][i].a*2 + k[3][i].a)*(dt/6);
        ret[i].a = {0., 0.};
    }
    /*    for (auto p = particles.begin(); p != particles.end(); p++, i++)
    {
        k[i] = *p;
        for (int j = 0; j < N; ++j)
            if (j != i)
                k[i].a += computeForce(ret[i], ret[j]);
    }

    for (auto p = particles.begin(); p != particles.end(); p++)
    {
        Particle pa = *p;
        pa.r = pa.r + pa.v*dt;
        ret.push_back(pa);
        }*/
    return ret;
}

void makeMovie(const vector<Shreds>& sequence, string fname)
{
    Box box = getBox(sequence);
    box.enlarge(0.1);
    double width = 1000;
    Geometry size(width, box.aspectRatio()*width);
    double scale = width / box.width();

    vector<Image> v;
    for (int i = 0; i < sequence.size(); ++i)
    {
        v.push_back(Image(size, "white"));
        drawFrame(&v[i], sequence[i], box, scale);
    }

    writeImages(v.begin(), v.end(), fname);
}

int main(int argc,char **argv) 
{
    InitializeMagick(*argv);
    Config cfg;
    readConfig(&cfg, "start.cfg");
    N = getProperty<int>("N", cfg);
    A = getProperty<double>("A", cfg);
    rmin = getProperty<double>("rmin", cfg);
    rmax = getProperty<double>("rmax", cfg);
    
    vector<Shreds> sequence;
    sequence.push_back(initShreds(cfg));

    for (int i = 1; i < 100; ++i)
    {
        sequence.push_back(moveParticles(sequence[i-1], 1.0));
    }

    //v[i].draw( DrawableText(300, 300, "I love Russ пупкин"));

    makeMovie(sequence, "ttt.gif");

    return 0; 
}
