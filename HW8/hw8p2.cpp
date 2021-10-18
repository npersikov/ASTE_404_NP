#include<iostream>
#include <random>
#include <math.h>
#include <fstream>

using namespace std;

std::random_device rnd_device;
std::mt19937 mt_gen(rnd_device());
std::uniform_real_distribution<double> rnd_dist;

double rnd();
double* solveGS(double *a, double *b, double *c, double **R, double *a_prev, double *b_prev, double *c_prev,
			double *T, double *T_prev, int numCells, int numTimeSteps, int it);

int main()
{
    double dx = 1;
    int numCells = 60;
    int numTimeSteps = 15000;
    double dt = 1e-5;
    double t_pulse = 1e-3 * rnd();
    int pulseMagnitude = 300;
    // Should this be a different random number? is n_i numCells?
    double i_pulse = 1 + (numCells - 2) * rnd(); // random internal point for pulse location
    double pulseRate = rnd()*(100/dt);
    int endPulseTimeStep = 3000;
    double D = 0.8; // Diffusion coefficient

    // Initialize 1D T array
    double *T = new double[numCells];

    // Clear data
  	for (int n = 0; n < numCells; n++) 
	{
    	T[n] = 0;  
        // cout<<T[n]<<endl;
	}

    // Initialize coefficient lists
    double *a = new double[numCells];
	double *b = new double[numCells];
	double *c = new double[numCells];
    double *a_prev = new double[numCells];
	double *b_prev = new double[numCells];
	double *c_prev = new double[numCells];

    // Set coefficients for k+1
    for(int i = 0; i < numCells; i++)
    {
        a[i] = -dt/2/pow(dx,2);
        b[i] = dt/pow(dx,2.0) + 1.0/D; // 1/D causes a floating point error
        cout<<a[i]<<endl;
        c[i] = a[i];
    }

    // Set coefficients for k
    for(int i = 0; i < numCells; i++)
    {
        a_prev[i] = dt/2/pow(dx,2);
        b_prev[i] = -dt/pow(dx,2.0) + 1.0/D; // 1/D causes a floating point error
        c_prev[i] = a[i];
    }

    // Allocate array of pulses
    double **R = new double*[numCells];
    for(int r = 0; r < numCells; r++)
    {
        R[r] = new double[numTimeSteps];
    }

    // Create an x by t array of the random pulses
    bool isPulse = false;
    double currentPulseTime = 0;
    double currentPulseMagnitude = 0;
    for(int t = 0; t < numTimeSteps; t++)
    {
        if(t < endPulseTimeStep)
        {
            // If there is no pulse, make one
            if(!isPulse)
            {
                t_pulse = 1e-3 * rnd();
                i_pulse = 1 + (numCells - 2) * rnd(); // random internal point for pulse location
                pulseRate = rnd()*(100/dt);
                isPulse = true;
                currentPulseTime = 0;
                currentPulseMagnitude = 0;
                // cout<<"made new pulse"<<endl;
            }
            else
            {
                // cout<<"current pulse"<<endl;
                R[(int)i_pulse][t] = currentPulseMagnitude;

                if(currentPulseMagnitude < pulseMagnitude)
                {
                    currentPulseMagnitude += pulseRate*dt;
                    // cout<<pulseRate*dt<<endl;
                }

                currentPulseTime += dt;
                // Stop the pulse if it lasts too long
                if(currentPulseTime >= t_pulse)
                {
                    // cout<<"pulse is over"<<endl;
                    isPulse = false;
                }
            } 
            // cout<<"step: "<<t<<" mag: "<<currentPulseMagnitude<<" pos: "<<(int)i_pulse<<" time: " <<currentPulseTime<< endl;
        }
    }

    // Allocate time history array
    // double **history = new double*[numCells];
    // for(int r = 0; r < numCells; r++)
    // {
    //     history[r] = new double[numTimeSteps];
    // }
    double *history = new double[numCells*numTimeSteps];
    double *T_prev = new double[numCells];
    for (int n = 0; n < numCells; n++) 
	{
    	T_prev[n] = 0;  
	}

    // Solve the system... I hope
    int histind = 0;
    for(int t = 0; t < numTimeSteps; t++)
    {   
        double *ans = solveGS(a, b, c, R, a_prev, b_prev, c_prev, T_prev, T_prev, numCells, numTimeSteps, t);
        for(int i = 0; i < numCells; i++)
        {
            history[histind] = ans[i];
            T_prev[i] = ans[i];

            if(ans[i] > 0)
                cout <<  ans[i] << endl;

            histind++;
        }
    }

    // Save the file
    ofstream out("field.vti");

  	out<<"<VTKFile type=\"ImageData\">\n";
  	out<<"<ImageData WholeExtent=\"0 "<<numTimeSteps-1<<" 0 "<<numCells-1<<" 0 "<<0<<"\""; 
  	out<<" Origin=\""<<0<<" "<<0<<" "<<0.0<<"\"";
  	out<<" Spacing=\""<<dx<<" " <<dt<<" "<<0.0<<"\">\n";
  	out<<"<Piece Extent=\"0 "<<numTimeSteps-1<<" 0 "<<numCells-1<<" 0 "<<0<<"\">\n"; 
  	out<<"<PointData>\n";

	out<<"<DataArray Name=\"T\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int n=0;n<numCells*numTimeSteps;n++) out<<history[n]<<" ";	
	out<<"\n</DataArray>\n";

	out<<"</PointData>\n";
	out<<"</Piece>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";



    // Deallocate all arrays

    for(int r = 0; r < numCells; r++)
        delete[] R[r];

    delete[] R;
    R = nullptr;


    delete[] history;

    delete[] a;
	delete[] b;
	delete[] c;

    delete[] T;

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
 * Matrix solver equation modified for 1D
 */
double* solveGS(double *a, double *b, double *c, double **R, double *a_prev, double *b_prev, double *c_prev,
			double *T, double *T_prev, int numCells, int numTimeSteps, int it)
{
    double *ans = new double[numCells];

    // cout<<T[1]<<endl;
    // cout<<T[20]<<endl;
    // cout<<T[55]<<endl;
	int nn = numCells;

    for(int n = 0; n < nn; n++)
    {
        // cout<<b[n]<<endl;
        double sum = 0;

        if(a[n] != 0)
        {
            sum += a[n]*T[n - 1];
        }
        if(c[n] != 0)
        {
            sum += c[n]*T[n + 1];
        }

        double g = 0;
        if(n > 0)
        {
            g = -a_prev[n]*T[n-1] - b_prev[n]*T[n] - c_prev[n]*T_prev[n+1] + R[n][it];
        }

        double T_star = (g - sum) / c[n];

        ans[n] = T[n] + 1.4*(T_star - T[n]);
    }

    // if(it%50 == 0)
    // {
    // 	double r2_sum = 0;

    // 	for(int n = 0; n < nn; n++)
    // 	{
    //         double g = -a[n]*T[n-1] - b[n]*T[n] - c[n]*T[n-1] + R[n][it];
    // 		double sum = 0;
    // 		if(a[n] != 0)
    // 		{
    // 			sum += a[n]*T[n - 1];
    // 		}
    // 		sum += b[n]*T[n]; // This must happen every time?

    // 		if(c[n] != 0)
    // 		{
    // 			sum += c[n]*T[n + 1];
    // 		}

    // 		double r = g - sum;
    // 		r2_sum += r*r;
    // 	}

    // 	double L2 = sqrt(r2_sum/nn);

    // 	cout << "solver iteration: " << it << ", L2 norm: " << L2 << endl;
    // 	if(L2 < 1e-6) return true; // break the loop if converged enough.
    // }

	return ans; // If you get here, that means you didn't converge. Big sad.
}
