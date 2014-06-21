#include <iostream>
#include <string>
#include <libconfig.h++>

using libconfig::Config;
using std::cout;
using std::cerr;
using std::endl;
using libconfig::FileIOException;
using libconfig::ParseException;
using libconfig::SettingNotFoundException;

int main()
{
    Config cfg;
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

    // Get the N
    try
    {
        int N = cfg.lookup("N");
        cout << "Number of particles: " << N << endl << endl;
    }
    catch(const SettingNotFoundException &nfex)
    {
        cerr << "No 'N' setting in configuration file." << endl;
    }

    return 0;
}
