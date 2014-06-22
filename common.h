#include <stdlib.h>
#include <cstring>
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

const double PI = 3.1415926;

struct Point
{
    double x, y;
    inline double len() const { return sqrt(x*x + y*y);}
    inline Point operator-(const Point& other) const
    {
        return {x - other.x, y - other.y};
    }
    inline Point operator+(const Point& other) const
    {
        return {x + other.x, y + other.y};
    }
    inline Point operator*(double k) const
    {
        return {x*k, y*k};
    }
    inline Point& operator +=(const Point& other)
    {
        x = other.x;
        y = other.y;
    }
};

struct Particle
{
    Point r, v;
    Point a;
};



template<typename T>
T getProperty(const char* name, const Config& cfg)
{
    try
    {
        return cfg.lookup(name);
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
}

template<>
string getProperty<string>(const char* name, const Config& cfg)
{
    return getProperty<const char*>(name, cfg);
}

template<>
Point getProperty<Point>(const char* name, const Config& cfg)
{
    Point ret;
    try
    {
        Setting& p = cfg.lookup(name);
        ret.x = p["x"];
        ret.y = p["y"];
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

template<typename T>
T getProperty(const char* name, const Config& cfg, T def)
{
    if (cfg.exists(name))
        return getProperty<T>(name, cfg);
    return def;
}

typedef vector<Particle> Shreds;
