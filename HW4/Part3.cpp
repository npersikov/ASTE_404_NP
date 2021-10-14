#include <fstream> // For file writing
#include <iostream>

using namespace std;

int main()
{
    float sum_f = 0;
    double sum_d = 0;
    double dt = 1e-7; // Not sure if this syntax works...

    std::ofstream data;
    data.open("fp.csv");
    data << "i\tsum_f\tsum_d\treal\terror_f\terror_d\n"; // Header for the output file
    //int skip = 10; // Apparently this skips every 10th particle.

    int total = 30000000;
    double real = 0;
    double error_f = 0;
    double error_d = 0;
    for(int i = 0; i < total; i++)
    {
        real = i*dt;
        sum_f += dt;
        sum_d += dt;

        if(i%10000 == 0)
        {
            error_f = (real - sum_f) / real;
            error_d = (real - sum_d) / real;
            data << i << "\t" << sum_f << "\t" << sum_d << "\t" << real << "\t" << error_f << "\t" << error_d << "\n";
        }
    }
    data.close();
}