
#include <stdlib.h>
#include <csignal>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <libconfig.h++>

using libconfig::Config;
using libconfig::Setting;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using libconfig::FileIOException;
using libconfig::ParseException;
using libconfig::SettingNotFoundException;
using libconfig::SettingException;

double PI = 3.1415926;

struct Point
{
    double x, y;
    double len() const { return sqrt(x*x + y*y);}
};

struct Particle
{
    Point r, v;
};

double random(double max)
{
    return (double)rand()*max/RAND_MAX;
}

double normalRandom(double dispersion)
{
    double U = random(1);
    double V = random(1);
    //if (U > .5)
    //return sqrt(-2*log(U))*sin(2*PI*V);
    return sqrt(-2*log(U))*cos(2*PI*V)*dispersion;
}

Point makeRandom(Point limit)
{
    return {random(limit.x), random(limit.y)};
}

Point makeNormalRandomP(double dispersion)
{
    return {normalRandom(dispersion), normalRandom(dispersion)};
}

Particle makeRandom(Point limit, double vdisp)
{
    return {makeRandom(limit), makeNormalRandomP(vdisp)};
}

template<typename T>
T getProperty(const char* name, const Config& cfg)
{
    T ret;
    try
    {
        ret = cfg.lookup(name);
    }
    catch(const SettingNotFoundException &nfex)
    {
        cerr << "No '" <<name <<"' setting in configuration file." << endl;
        std::raise(SIGTERM);
    }
    catch(const SettingException &tex)
    {
        cerr << "Error on parameter " <<name <<":" <<tex.what() <<endl;
        std::raise(SIGTERM);
    }
    return ret;
}

typedef vector<Particle> Shreds;

Shreds generate(int N, Point box, double vdisp)
{
    Shreds ret;
    for (int i = 0; i < N; ++i)
    {
        ret.push_back(makeRandom(box, vdisp));
    }
    return ret;
}

void store(const Shreds& dots, Setting* arr)
{
    for (int i = 0; i < dots.size(); ++i)
    {
        Setting &p = arr->add(Setting::TypeGroup);
        Setting &r = p.add("r", Setting::TypeGroup);
        r.add("x", Setting::TypeFloat) = dots[i].r.x;
        r.add("y", Setting::TypeFloat) = dots[i].r.y;
        Setting &v = p.add("v", Setting::TypeGroup);
        v.add("x", Setting::TypeFloat) = dots[i].v.x;
        v.add("y", Setting::TypeFloat) = dots[i].v.y;
    }
}

double computeDistribution(const Shreds& dots, int* arr, int len)
{
    double vmax = 0;
    for (int i = 0; i < dots.size(); ++i)
    {
        double v = dots[i].v.len();
        if (v > vmax) vmax = v;
    }
    for (int i = 0; i < dots.size(); ++i)
    {
        double v = dots[i].v.len();
        arr[(int)(v*len/vmax + 0.5)]++;
    }
    return vmax;
}

int main(int argc, char** argv)
{
    Config cfg;
    const char* fname = "example.cfg";
    if (argc>1) fname = argv[1];
    try
    {
        cfg.readFile("example.cfg");
    }
    catch(FileIOException &fioex)
    {
        std::cerr << "I/O error while reading file." << std::endl;
        return 1;
    }
    catch(ParseException &pex)
    {
        std::cerr << "Parse error at " << pex.getLine()
                  << " - " << pex.getError() << std::endl;
        return 1;
    }

    int N = getProperty<int>("N", cfg);
    
    srand(time(NULL));

    Config out;
    Setting &root = out.getRoot();
    root.add("N", Setting::TypeInt) = N;

    Setting &coords = root.add("particles", Setting::TypeList);
   
    vector<Particle> dots = generate(N, {100., 100.}, 10.);
    const int nCells = 10;
    int dist[nCells] = {0};
    double maxV = computeDistribution(dots, dist, nCells);
    cout <<endl;
    cout <<"Generated velocity distribution:" <<endl;
    for (int i = 0; i < nCells; ++i)
    {
        int len = dist[i];
        cout <<std::setprecision(3) <<std::fixed <<std::setw(8) <<i*maxV/nCells <<" ";
        while (len--) cout <<"#";
        cout <<endl;
    }
    store(dots, &coords);

    try
    {
        out.writeFile("start.cfg");
    }
    catch(const FileIOException &fioex)
    {
        cerr << "I/O error while writing output file" <<endl;
        return 1;
    }

    return 0;
}
