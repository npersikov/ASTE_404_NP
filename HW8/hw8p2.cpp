#include<iostream>
#include <random>
#include <math.h>
#include <fstream>

using namespace std;

std::random_device rnd_device;
std::mt19937 mt_gen(rnd_device());
std::uniform_real_distribution<double> rnd_dist;

double rnd();
double* solveNextRow(double *T_old, int width, double dt, double D, double dx);
double calcTemp(double T_left, double T_center, double T_right, double dt, double D, double dx);

int main()
{
    double dx = 1;
    int numCells = 60;
    int numTimeSteps = 15000;
    double dt = 1e-5;
    int pulseMagnitude = 300;
    int endPulseTimeStep = 3000;
    double D = 1; // Diffusion coefficient

    // Initialize 1D T array
    double *T_old = new double[numCells];
  	for (int n = 0; n < numCells; n++) // Clear data, set all temperatures to 0 K
    	T_old[n] = 0;  

    // Initialize array for saving and plotting
    double *T_all = new double[numCells*numTimeSteps];
    for (int n = 0; n < numCells*numTimeSteps; n++) // Clear data, set all temperatures to 0 K
    	T_all[n] = 0;

    bool isActivePulse = false;
    double t_pulse = 0.0;
    double i_pulse = 0.0;
    double pulseRate = 0.0;
    double currentPulseMagnitude = 0.0;
    double pulseTimeElapsed = 0.0;
    for(int ts = 0; ts < numTimeSteps; ts++)
    {
        // Resolve the current pulse or create a new one if the old one is done.
        if(ts < endPulseTimeStep)
        {
            if(!isActivePulse)
            {
                isActivePulse = true;
                t_pulse = 1e-3 * rnd();
                i_pulse = 1 + (numCells - 2) * rnd();
                pulseRate = rnd()*(100/dt);
                currentPulseMagnitude = 0.0;
                pulseTimeElapsed = 0.0;
            }
            else
            {
                if(pulseTimeElapsed >= t_pulse)
                {
                    isActivePulse = false;
                    currentPulseMagnitude = 0.0;
                }
                else if(currentPulseMagnitude < pulseMagnitude)
                {
                    currentPulseMagnitude += pulseRate*dt;
                    // Only add pulse to the rod if the pulse is on.
                    T_old[(int)i_pulse] = currentPulseMagnitude; 
                }
                pulseTimeElapsed += dt;
            }
        }

        double *T_new = solveNextRow(T_old, numCells, dt, D, dx);

        // if(T_new[30] > 5)
        //     cout << T_new[30] << endl;

        for(int cell = 0; cell < numCells; cell++)
            T_all[numCells*ts + cell] = T_old[cell];

        T_old = T_new;
        // if(T_old[30] > 5)
        //     cout << T_old[30] << endl;
    }




    // Save the file
    ofstream out("Unsteady_1D_Heat_Sim.vti");

  	out<<"<VTKFile type=\"ImageData\">\n";
  	out<<"<ImageData WholeExtent=\"0 "<<numCells-1<<" 0 "<<numTimeSteps-1<<" 0 "<<0<<"\""; 
  	out<<" Origin=\""<<0<<" "<<0<<" "<<0.0<<"\"";
  	out<<" Spacing=\""<< dx <<" " << dt*1000 <<" "<<0.0<<"\">\n";
  	out<<"<Piece Extent=\"0 "<<numCells-1<<" 0 "<<numTimeSteps-1<<" 0 "<<0<<"\">\n"; 
  	out<<"<PointData>\n";

	out<<"<DataArray Name=\"T\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int n=0;n<numCells*numTimeSteps;n++)
    {
        out<<T_all[n]<<" ";	

        // if(T_all[n] > 5)
        //     cout << T_all[n] << " K at index " << n << endl;
    }

	out<<"\n</DataArray>\n";

	out<<"</PointData>\n";
	out<<"</Piece>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";



    // Deallocate all arrays
    
    delete[] T_old;
    delete[] T_all;

    return 0;
}

/**
 * @return a random double R[0,1)
 */
double rnd()
{
	return rnd_dist(mt_gen);
}

/**
 * The logic behind solving the 1D unsteady heat equation. Currently, there are no
 * pulses yet.
 * 
 * @return an array of the next rows temperatures.
 */
double* solveNextRow(double *T_old, int width, double dt, double D, double dx)
{
    double *T_new = new double[width];
    double T_left = 0;
    double T_center = 0;
    double T_right = 0;
    for(int i = 0; i < width; i++)
    {
        if(i > 0)
            T_left = T_old[i - 1];
        else
            T_left = 0;

        T_center = T_old[i];

        if(i < width)
            T_right = T_old[i + 1];
        else
            T_right = 0;

        T_new[i] = calcTemp(T_left, T_center, T_right, dt, D, dx);
    }

    return T_new;
}

/**
 * Separating the equation that solves for future temperatures from the algrithm
 * that solves an entire row for visual clarity.
 * 
 * @return the temperature at a cell determined by the location of the last row's
 * temperatures at the same and adjacent positions.
 */
double calcTemp(double T_left, double T_center, double T_right, double dt, double D, double dx)
{
    double T = T_center + dt*D*(T_left - 2*T_center + T_right)/(dx*dx);
    // if(T > 5)
    //     cout << T << endl;
    return T;
}