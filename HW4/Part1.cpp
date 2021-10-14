#include <iostream>
#include <math.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main()
{

    size_t N = 100000000; // size_t is an unsigned integer type...
    float *r = new float[N]; // allocate the array

    // Allocate initial data
    for(size_t i = 0; i < N; i++)
        r[i] = i;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // Compute the magnitude...
    double sum = 0;
    for(size_t i = 0; i < N; i++)
        sum += r[i]*r[i];

    double mag = sqrt(sum);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    cout << "Mag r[" << N << "] is " << mag << endl;
    cout << "Calculation took " << time_span.count() << " seconds" << endl;

    delete[] r; //free up the memory
    return 0;
}