
#include "common.h"

using std::ofstream;

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

void storeDual(const Shreds& dots, const string& filePrefix)
{
    string rfname = filePrefix+dataRPostfix;
    string vfname = filePrefix+dataVPostifx;
    ofstream coords(rfname);
    ofstream vels(vfname);
    if (coords.fail())
    {
        cerr << "I/O error while opening output file " <<rfname <<endl;
        std::raise(SIGTERM);
    }
    if (vels.fail())
    {
        cerr << "I/O error while opening output file " <<vfname <<endl;
        std::raise(SIGTERM);
    }

    for (int i = 0; i < dots.size(); ++i)
    {
        coords <<dots[i].r.x <<" "<<dots[i].r.y <<endl;
        vels <<dots[i].v.x <<" "<<dots[i].v.y <<endl;
    }

    if (coords.fail())
    {
        cerr << "I/O error while writing output file " <<rfname <<endl;
        std::raise(SIGTERM);
    }
    if (vels.fail())
    {
        cerr << "I/O error while writing output file " <<vfname <<endl;
        std::raise(SIGTERM);
    }
    coords.close();
    vels.close();
    cout <<"Coordinates are written to "<<rfname<<endl;
    cout <<"Velocities are wirtten to "<<vfname<<endl;
}

double computeDistribution(const Shreds& dots, int* arr, int len)
{
    double vmax = 0;
    for (int i = 0; i < len; ++i)
    {
        arr[i] = 0;
    }
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

void drawHistogramm(const Shreds& dots, int nCells)
{
    if (nCells <= 0) return;
    int* dist = new int[nCells];
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
}

void ballance(Shreds* dots)
{
    unsigned size = dots->size();
    unsigned half = size/2;
    for (unsigned i = 0; i < half; ++i)
    {
        (*dots)[i].v = (*dots)[half + i].v*-1;
    }
    if (half*2 < size)
        (*dots)[half].v = {0., 0.};
}

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

void writeFile(Config& rezult, string fname)
{
    try
    {
        rezult.writeFile(fname.c_str());
        cout <<"The generated configuration is wirtten to "<<fname <<endl;
    }
    catch(const FileIOException &fioex)
    {
        cerr << "I/O error while writing output file" <<endl;
        std::raise(SIGTERM);
    }
}

int main(int argc, char** argv)
{
    Config cfg;
    const char* fname = "example.cfg";
    if (argc>1) fname = argv[1];
    else
    {
        cout <<"No options supplied. Using default configuration file: " <<fname <<endl;
    }

    readConfig(&cfg, fname);
    int N = getProperty<int>("N", cfg);
    Point box = getProperty<Point>("box", cfg, {100., 100.});
    double vDisp = getProperty<double>("vDisp", cfg, 10.);
    string outputFname = getProperty<string>("output", cfg, "start.cfg");
    int nHistogramm = getProperty<int>("hist", cfg, 0);
    string dataFilePrefix = getProperty<string>("dataFilePrefix", cfg, "");
    bool ballanced = getProperty<bool>("ballanced", cfg, false);
    bool keepInBox = getProperty<bool>("keepInBox", cfg, false);
    
    srand(time(NULL));

    Config out;
    Setting &root = out.getRoot();
    root.add("N", Setting::TypeInt) = N;

   
    vector<Particle> dots = generate(N, box, vDisp);
    if (ballanced)
    {
        ballance(&dots);
    }
    drawHistogramm(dots, nHistogramm);
    root.add("rmin", Setting::TypeFloat) = getProperty<double>("rmin", cfg, 1.);
    root.add("rmax", Setting::TypeFloat) = getProperty<double>("rmax", cfg, 10.);
    root.add("deviation", Setting::TypeFloat) = getProperty<double>("deviation", cfg, 10.);
    root.add("A", Setting::TypeFloat) = getProperty<double>("A", cfg, 1.);
    root.add("a", Setting::TypeFloat) = getProperty<double>("a", cfg, 0.);
    root.add("b", Setting::TypeFloat) = getProperty<double>("b", cfg, 0.);
    root.add("r0", Setting::TypeFloat) = getProperty<double>("r0", cfg, 1.);
    root.add("epsilon", Setting::TypeFloat) = getProperty<double>("epsilon", cfg, 1.);
    root.add("sigma", Setting::TypeFloat) = getProperty<double>("sigma", cfg, 1.);
    root.add("m", Setting::TypeFloat) = getProperty<double>("m", cfg, 1.);
    root.add("G", Setting::TypeFloat) = getProperty<double>("G", cfg, 1.);
    root.add("T", Setting::TypeFloat) = getProperty<double>("T", cfg, 15.);
    root.add("step", Setting::TypeFloat) = getProperty<double>("step", cfg, 15.);
    root.add("width", Setting::TypeInt) = getProperty<int>("width", cfg, 1000);
    root.add("rezultFile", Setting::TypeString) = getProperty<string>("rezultFile", cfg, "sim.gif");
    root.add("delay", Setting::TypeInt) = getProperty<int>("delay", cfg, 30);
    root.add("endDelay", Setting::TypeInt) = getProperty<int>("endDelay", cfg, 300);
    root.add("vScaleFactor", Setting::TypeFloat) = getProperty<double>("vScaleFactor", cfg, 0.);
    root.add("drawEnergy", Setting::TypeBoolean) = getProperty<bool>("drawEnergy", cfg, false);
    root.add("dumpPointsSeparately", Setting::TypeBoolean) =
        getProperty<bool>("dumpPointsSeparately", cfg, true);
    root.add("askForSpecificDumps", Setting::TypeBoolean) =
        getProperty<bool>("askForSpecificDumps", cfg, true);
    root.add("dumpPointsOnEachFrame", Setting::TypeBoolean) =
        getProperty<bool>("dumpPointsOnEachFrame", cfg, true);
    root.add("t0", Setting::TypeFloat) = getProperty<double>("t0", cfg, 0.);


    if (keepInBox)
    {
        Setting &r = root.add("box", Setting::TypeGroup);
        r.add("x", Setting::TypeFloat) = box.x;
        r.add("y", Setting::TypeFloat) = box.y;
        root.add("keepInBox", Setting::TypeBoolean) = true;
    }

    if (dataFilePrefix.empty())
    {
        Setting &coords = root.add("particles", Setting::TypeList);
        store(dots, &coords);
    }
    else
    {
        root.add("dataFilePrefix", Setting::TypeString) = dataFilePrefix;
        storeDual(dots, dataFilePrefix);
    }
    writeFile(out, outputFname);

    return 0;
}

