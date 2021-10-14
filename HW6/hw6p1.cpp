#include <iostream>

using namespace std;

/*
 Cases a, c, and d do not work.
 Case b works!
 Case c results in a compile-time error.
*/

class Spacecraft
{
    public:
        void virtual launch()
        {
            cout << "Can't launch a generic spacecraft! Signed, the FAA." << endl;
        }
        // void launch() = 0;
};

class MMS : public Spacecraft
{
    public:
        void launch()
        {
            cout << "Launching MMS!" << endl;
        }
};

class GOESR : public Spacecraft
{
    public:
        void launch()
        {
            cout << "Launching GOES-R!" << endl;
        }
};

class JWST : public Spacecraft
{
    public:
        void launch()
        {
            cout << "Finally launching JWST!" << endl;
        }
};

void run_simulation(Spacecraft &sc)
{
    sc.launch();
}

int main()
{
    MMS mms;
    GOESR goesr;
    JWST jwst;

    run_simulation(mms);
    run_simulation(goesr);
    run_simulation(jwst);

    return 0;
}