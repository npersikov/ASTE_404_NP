/*
ASTE-404 HW5 initial file
creates and outputs an empty 2D mesh
*/

#include <iostream>			// for screen output
#include <fstream>			// for file writing
#include <math.h>			// math functions but unused	
#include <random>

using namespace std;

std::random_device rnd_device;
std::mt19937 mt_gen(rnd_device());
std::uniform_real_distribution<double> rnd_dist;

bool solveGS(double *a, double *b, double *c, double *d, double *e, double *g, 
			double *T, int ni, int nj); // GS Solver function prototype

double rnd();

int main() 
{
  	int ni = 81;   // number of nodes
  	int nj = 61;

  	double x0 = 0.1;  // origin
  	double y0 = 0.1;

  	double dx = 0.01;  // cell spacing
  	double dy = 0.02;
	// double dx = 1/(ni - 1);
	// double dy = 1/(nj - 1);

  	int nn = ni*nj;  // total number of nodes

  	// allocate 1D memory block
  	double *T = new double[ni*nj];

  	// clear data (not initialized by default)
  	for (int n=0;n<nn;n++) 
	{
    	T[n] = 0;  
	}

    // Allocate memory for the matrix coefficients
	double *a = new double[nn];
	double *b = new double[nn];
	double *c = new double[nn];
	double *d = new double[nn];
	double *e = new double[nn];
	double *g = new double[nn];

	// Clear data
	for(int n = 0; n < nn; n++)
	{
		a[n] = 0;
		b[n] = 0;
		c[n] = 0;
		d[n] = 0;
		e[n] = 0;
		g[n] = 0; 
	}

	// Set matrix values
	for(int j = 0; j < nj; j++)
	{
		for(int i = 0; i < ni; i++)
		{
			int n = j*ni + i;

			// Default code without Neumann BCs
			// if(i == 0 || j == 0)
			// {
			// 	c[n] = 1.0;
			// 	g[n] = 0.0; // Zero Kelvin on minimum i or j
			// }
			// else if(i == ni - 1 || j == nj - 1)
			// {
			// 	c[n] = 1.0;
			// 	g[n] = 100.0; // 100 Kelvin on maximum i or j
			// }

			// Implement Neumann boundary conditions like in HW2
			if(i == 0)
			{
				c[n] = -1/dx;
				d[n] = 1/dx;
				g[n] = 0.0;
			}
			else if(i == ni - 1)
			{
				c[n] = -1/dx;
				b[n] = 1/dx;
				g[n] = 100.0;
			}
			else if(j == 0)
			{
				c[n] = -1/dy;
				e[n] = 1/dy;
				g[n] = 0.0;
			}
			else if(j == nj - 1)
			{
				c[n] = -1/dy;
				a[n] = 1/dy;
				g[n] = 100.0;
			}
			else
			{
				a[n] = e[n] = 1/(dy*dy); // I didn't know you could do this
				b[n] = d[n] = 1/(dx*dx);
				c[n] = -2/(dx*dx) - 2/(dy*dy);
				g[n] = 0;
			}
		}
	}

	// Set random internal Dirichlet nodes
	for(int s = 0; s < 20; s++)
	{
		int i = (int) (1 + rnd()*(ni - 2)); // Random int within 1 and ni-2
		int j = (int) (1 + rnd()*(nj - 2));
		int n = j*ni + i;

		a[n] = b[n] = c[n] = d[n] = e[n] = g[n] = 0; // Clear all row data

		// Make a Dirichlet node
		c[n] = 1.0; // Main diagonal
		g[n] = rnd()*100; // Random value from 0 to 99
	}

	// Solve the matrix system
	solveGS(a, b, c, d, e, g, T, ni, nj);

	// Release memory allocated for the matrix
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d;
	delete[] e;
	delete[] g;

  	/* output vti file */
  	ofstream out("field.vti");

  	out<<"<VTKFile type=\"ImageData\">\n";
  	out<<"<ImageData WholeExtent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 "<<0<<"\""; 
  	out<<" Origin=\""<<x0<<" "<<y0<<" "<<0.0<<"\"";
  	out<<" Spacing=\""<<dx<<" " <<dy<<" "<<0.0<<"\">\n";
  	out<<"<Piece Extent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 "<<0<<"\">\n"; 
  	out<<"<PointData>\n";

	out<<"<DataArray Name=\"T\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int n=0;n<nn;n++) out<<T[n]<<" ";	
	out<<"\n</DataArray>\n";

	out<<"</PointData>\n";
	out<<"</Piece>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";

	// free memory
	delete[] T;
	return 0;	// normal exit
}

bool solveGS(double *a, double *b, double *c, double *d, double *e, double *g, 
			double *T, int ni, int nj)
{
	int nn = ni * nj;

	for(int it = 0; it < 10000; it++)
	{
		for(int n = 0; n < nn; n++)
		{
			double sum = 0;

			if(a[n] != 0)
			{
				sum += a[n]*T[n - ni];
			}
			if(b[n] != 0)
			{
				sum += b[n]*T[n - 1];
			}
			if(d[n] != 0)
			{
				sum += d[n]*T[n + 1];
			}
			if(e[n] != 0)
			{
				sum += e[n]*T[n + ni];
			}

			double T_star = (g[n] - sum) / c[n];

			T[n] = T[n] + 1.4*(T_star - T[n]);
		}

		if(it%50 == 0)
		{
			double r2_sum = 0;

			for(int n = 0; n < nn; n++)
			{
				double sum = 0;
				if(a[n] != 0)
				{
					sum += a[n]*T[n - ni];
				}
				if(b[n] != 0)
				{
					sum += b[n]*T[n - 1];
				}
				sum += c[n]*T[n]; // This must happen every time?

				if(d[n] != 0)
				{
					sum += d[n]*T[n + 1];
				}
				if(e[n] != 0)
				{
					sum += e[n]*T[n + ni];
				}

				double r = g[n] - sum;
				r2_sum += r*r;
			}

			double L2 = sqrt(r2_sum/nn);

			cout << "solver iteration: " << it << ", L2 norm: " << L2 << endl;
			if(L2 < 1e-6) return true; // break the loop if converged enough.
		}
	}

	return false; // If you get here, that means you didn't converge. Big sad.
}

double rnd()
{
	return rnd_dist(mt_gen);
}