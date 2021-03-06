#include <math.h>
#include <Magick++.h>
#include <Magick++/STL.h>
#include <fstream>
#include <sstream>
#include <list>
#include <stdexcept>
#include "common.h"
#include "potential.h"
#include "integrator.h"

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
using std::cin;

double T;
int N;
double maxDev;
int screenWidth;
int delay, endDelay;
double vscalefactor;

Point box;
bool keepInBox;
bool drawEnergy = false,
    dumpPointsSeparately = true,
    stepautofit = false;

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

Point toScreen(Point r, const Box& box, double scale)
{
    Point leftdown = {box.leftup.x, box.rightdown.y};
    Point tmp = (r - leftdown)*scale;
    tmp.y *= -1;
    return tmp;
}

void drawFrame(Image* img, const Shreds &parts, const Box& box, double scale)
{
    img->fillColor("red");
    img->strokeColor("black");
    img->strokeWidth(1);
    for (auto p = parts.begin(); p != parts.end(); p++)
    {
        Point coor = toScreen(p->r, box, scale);
        int radius = 1;
        img->draw(DrawableCircle
                  (coor.x - radius, coor.y - radius,
                   coor.x + radius, coor.y + radius));

        Point vsegm = p->v*vscalefactor;
        Point vend = toScreen(p->r + vsegm, box, scale);
        img->draw(DrawableLine(coor.x, coor.y,
                               vend.x, vend.y));
    }
}

double computeDeviation(const Shreds& p1, const Shreds& p2)
{
    double dev = 0;
    for (int i = 0; i < N; ++i)
    {
        dev += (p1[i].r - p2[i].r).len();
    }
    return dev/deviationFactor;
}

void KeepInBox(Shreds* parts)
{
    if (keepInBox)
    {
        for (int i = 0; i < N; ++i)
        {
            bool moved = false;
            do
            {
                moved = false;
                if ((*parts)[i].r.x < 0)
                {
                    (*parts)[i].r.x *= -1;
                    if ((*parts)[i].v.x < 0)
                        (*parts)[i].v.x *= -1;
                    moved = true;
                }
                if ((*parts)[i].r.y < 0)
                {
                    (*parts)[i].r.y *= -1;
                    if ((*parts)[i].v.y < 0)
                        (*parts)[i].v.y *= -1;
                    moved = true;
                }
                if ((*parts)[i].r.x > box.x)
                {
                    (*parts)[i].r.x = fmod((*parts)[i].r.x, 2*box.x);
                    (*parts)[i].r.x = 2*box.x - (*parts)[i].r.x;
                    if ((*parts)[i].v.x > 0)
                        (*parts)[i].v.x *= -1;
                    moved = true;
                }
                if ((*parts)[i].r.y > box.y)
                {
                    (*parts)[i].r.y = fmod((*parts)[i].r.y, 2*box.y);
                    (*parts)[i].r.y = 2*box.y - (*parts)[i].r.y;
                    if ((*parts)[i].v.y > 0)
                        (*parts)[i].v.y *= -1;
                    moved = true;
                }
            } while (moved);
        }
    }
}

Scene moveParticles(const Scene& scene, double step, double* h)
{
    if (!stepautofit)
    {
        Shreds pp = tryMove(scene.particles, *h);
        KeepInBox(&pp);
        return Scene({pp, scene.time + *h, *h, 0.});
    }
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



string double2string(double val)
{
    ostringstream oss;
    oss <<std::setprecision(3) <<val;
    return oss.str();
}

string point2string(Point p)
{
    ostringstream oss;
    oss <<std::setprecision(2) <<"(" <<p.x <<", "
        <<p.y <<")";
    return oss.str();
}

void drawCoords(Image* img, const Box& original, double scale,
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


    text.back() = DrawableText(size.width() - 120, size.height() - 40,
                               point2string({original.rightdown.x,
                                           original.leftup.y}));
    img->draw(text);
    text.back() = DrawableText(30, 40, point2string
                               ({original.leftup.x,
                                       original.rightdown.y}));
    img->draw(text);
    text.back() = DrawableText(30, size.height() - 40,
                               point2string(original.leftup));
    img->draw(text);
    text.back() = DrawableText(30, size.height()/2,
                               double2string((original.leftup.y + original.rightdown.y)/2));
    img->draw(text);
    text.back() = DrawableText(size.width()/2, size.height() - 40,
                               double2string((original.leftup.x + original.rightdown.x)/2));
    img->draw(text);
}

void drawBox(Image* img, const Box& original,
             const Box& box, double scale,
             const Geometry& size)
{
    drawCoords(img, original, scale, size);
    img->fillColor("black");
    img->strokeColor("gray");
    img->strokeWidth(1);

    Point leftup = toScreen(original.leftup, box, scale);
    Point rightdown = toScreen(original.rightdown, box, scale);

    img->draw(DrawableLine(leftup.x, leftup.y,
                           leftup.x, rightdown.y));
    img->draw(DrawableLine(leftup.x, rightdown.y,
                           rightdown.x, rightdown.y));
    img->draw(DrawableLine(rightdown.x, rightdown.y,
                           rightdown.x, leftup.y));
    img->draw(DrawableLine(rightdown.x, leftup.y,
                           leftup.x, leftup.y));

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
    text.back() = DrawableText(size.width() - 110, 20,
                               "t:      " + double2string(t) + "\n" + 
                               (stepautofit? 
                                ("h:      " + double2string(h) + "\n"
                                 "e r r: "  + double2string(deviation) + "\n") :
                                ("#:      " + double2string(t/h) + "\n")) +
                               (drawEnergy ?
                                ("K:      " + double2string(K) + "\n"
                                 "U:      " + double2string(U) + "\n"
                                 "E:      " + double2string(K+U)):
                                std::string()));
    img->draw(text);
}

void makeMovie(const vector<Scene>& sequence, string fname, int nStep)
{
    Box originalBox = getBox(sequence);
    Box box = originalBox;
    box.enlarge(0.1);
    double width = screenWidth;
    Geometry size(width, box.aspectRatio()*width);
    double scale = width / box.width();

    cout <<"Rendering frames..." <<endl;
    vector<Image> v;
    int fifth = (sequence.size() != 0 ?
                 sequence.size()/5:
                 sequence.size()+1);
    int imgNum = 0;
    for (int i = 0; i < sequence.size(); ++i)
    {
        if (!stepautofit && (((i/nStep)*nStep - i) != 0)) continue;
        v.push_back(Image(size, "white"));
        drawBox(&v[imgNum], originalBox, box, scale, size);
        drawFrame(&v[imgNum], sequence[i].particles, box, scale);
        drawInfo(&v[imgNum], box, scale, sequence[i].time, sequence[i].h,
                 sequence[i].deviation, sequence[i].U, sequence[i].K,
                 size);
        v[imgNum].animationDelay(delay);
        if ((i - (i/fifth)*fifth) == 0)
        {
            cout <<i/fifth*20 <<"% ";
            cout.flush();
        }
        ++imgNum;
    }
    cout <<endl;
    v.back().animationDelay(endDelay);
    cout <<"Storing frames to file "<<fname<<endl;
    writeImages(v.begin(), v.end(), fname);
}

void dumpPoints(const Shreds& particles,
                std::string rfname, std::string vfname,
                std::ios_base::openmode mode)
{
    std::ofstream fr(rfname, mode);
    std::ofstream fv(vfname, mode);
    
    if (fr.bad() || fr.fail())
    {
        std::cerr <<"Error: can\'t open file " <<rfname <<std::endl;
        return;
    }
    if (fv.bad() || fv.fail())
    {
        std::cerr <<"Error: can\'t open file " <<vfname <<std::endl;
        return;
    }
    
    for (int i = 0; i < N; ++i)
    {
        fr <<particles[i].r.x <<"    "<<particles[i].r.y <<endl;
        fv <<particles[i].v.x <<"    "<<particles[i].v.y <<endl;
    }
    fr<<endl<<endl;
    fv<<endl<<endl;

    fr.close();
    fv.close();
}

void describe(int iter, const Scene& sc)
{
    cout <<"k = " <<sc.K <<"  u = " <<sc.U <<" e = " <<sc.K + sc.U<<endl;
    std::string fname = "State/" + 
        std::to_string(iter) + "_" + double2string(sc.time) + ".txt";
    std::ofstream file(fname);
    
    if (file.bad() || file.fail())
    {
        std::cerr <<"Error: can\'t open file " <<fname
                  <<" (is State dir existing?)." <<std::endl;
        return;
    }
    file <<"#moment="<<iter<<"	t=" <<sc.time <<"	k="
         <<sc.K <<"	u="<<sc.U <<"	e=" <<sc.K + sc.U<<endl;
    file <<"#n= " <<N<<endl;
    if (dumpPointsSeparately)
    {
        dumpPoints(sc.particles,
                   "State/particles.r" + std::to_string(iter) + ".txt",
                   "State/particles.v" + std::to_string(iter) + ".txt",
                   std::ios_base::out|std::ios_base::trunc);
    }
    else
    {
        file <<"#x	y	vx	vy"<<endl;
        for (int i = 0; i < N; ++i)
        {
            file <<sc.particles[i].r.x <<"    "<<sc.particles[i].r.y <<"    "
                 <<sc.particles[i].v.x <<"    "<<sc.particles[i].v.y <<endl;
        }
    }
    file.close();
}

void describeNeccessaryFrames(vector<Scene>& sequence)
{
    int stN = 0;
    do
    {
        cout <<"Input state # to get it's energy:>";
        cout.flush();
        std::string word;
        if (cin >> word)
        {
            try
            {
                stN = std::stoi(word);
            }
            catch(std::invalid_argument ex)
            {
                break;
            }
            if (stN == 0) break;
            if (0 > stN || stN >= sequence.size())
                cout <<"Nonexistent state "<< stN <<endl;
            else
                describe(stN, sequence[stN]);
        }
        else break;

    } while (true);
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
    maxDev = getProperty<double>("deviation", cfg);
    screenWidth = getProperty<int>("width", cfg);
    delay = getProperty<int>("delay", cfg);
    endDelay = getProperty<int>("endDelay", cfg);
    vscalefactor = getProperty<double>("vScaleFactor", cfg);
    readUParams(cfg);
    string outFname = getProperty<string>("rezultFile", cfg);
    keepInBox = getProperty<bool>("keepInBox", cfg, false);
    if (keepInBox)
        box = getProperty<Point>("box", cfg, Point{0., 0.});
    drawEnergy = getProperty<bool>("drawEnergy", cfg, false);
    dumpPointsSeparately = getProperty<bool>("dumpPointsSeparately", cfg, true);
    bool askForSpecificDumps = getProperty<bool>("askForSpecificDumps", cfg, true);
    bool dumpPointsOnEachFrame = getProperty<bool>("dumpPointsOnEachFrame", cfg, true);
    double t0 = getProperty<double>("t0", cfg, 0.);
    

    vector<Scene> sequence;
    sequence.push_back({initShreds(cfg), t0, 0., 0.});

    stepautofit = getProperty<bool>("stepautofit", cfg, false);
    double step = 0, h = 0;
    int NStep = 0;
    if (stepautofit)
    {
        step = getProperty<double>("step", cfg);
        h = step;
        NStep = 0;
        T = getProperty<double>("T", cfg) + t0;
    }
    else
    {
        h = getProperty<double>("h", cfg);
        NStep = getProperty<int>("Nstep", cfg);
        step = 0;
        T = getProperty<int>("NT", cfg)*h + t0;
    }


    if (dumpPointsOnEachFrame)
    {
        dumpPoints(sequence[0].particles, "r.txt", "v.txt",
                   std::ios_base::trunc|std::ios_base::out);
    }
    for (int i = 1; ; ++i)
    {
        sequence.push_back(moveParticles(sequence[i-1], step, &h));
        Point E = computeEnergy(sequence[i].particles);
        sequence[i].U = E.x;
        sequence[i].K = E.y;
        cout <<"computed state # " <<i;
        bool anotherFrame = false;
        if (stepautofit)
        {
            anotherFrame = true;
        }
        else
            if ((i/NStep)*NStep - i == 0)
            {
                cout <<"(frame #"<< i/NStep <<")";
                anotherFrame = true;
            }
        cout <<" time: " <<sequence[i].time <<endl;
        if (anotherFrame && dumpPointsOnEachFrame)
            dumpPoints(sequence[i].particles, "r.txt", "v.txt",
                       std::ios_base::app|std::ios_base::ate|
                       std::ios_base::out);
        if (sequence[i].time > T) break;
    }

    makeMovie(sequence, outFname, NStep);

    if (askForSpecificDumps)
        describeNeccessaryFrames(sequence);

    return 0;
}
