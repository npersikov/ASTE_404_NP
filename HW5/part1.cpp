#include <fstream> // For file writing
#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;

void doubleMatrix(double **A, size_t nr, size_t nc);

int main()
{
    // Set dimensions
    size_t nr = 2000, nc = 2000;

    // Dynamically allocate the array
    double **A = new double*[nr];
    for(size_t i = 0; i < nr; i++) A[i] = new double[nc];

    // Fill array with some values
    for(int i = 0; i < nr; i++)
        for(int j = 0; j < nc; j++)
            A[i][j] = i*j;

    // Save the time before the experiment
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // Run the matrix operations
    for(int s = 0; s < 100; s++)
        doubleMatrix(A, nr, nc);

    // Save the time after the experiment
    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    // Get the time the experiment took. Not sure what's going on to do that
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1); 

    // Display result
    cout << "Multiplication took " << time_span.count() << " seconds." << endl;

    // Clear allocated memory
    for(int i = 0; i < nr; i++)
        delete[] A[i];

    delete[] A;
}

void doubleMatrix(double **A, size_t nr, size_t nc)
{
    // for(int i = 0; i < nr; i++)
    //     for(int j = 0; j < nc; j++)
    //         A[i][j] = 2*A[i][j];

    for(int j = 0; j < nc; j++)
        for(int i = 0; i < nr; i++)
            A[i][j] = 2*A[i][j];
}