#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

using namespace std;

class Matrix
{
    public:

        Matrix(int nr) : nr{nr} {}
        virtual ~Matrix() {}

        int getNr()
        {
            return nr;
        }

        virtual double dotRow(int r, const vector<double> &x) const = 0;

        friend vector<double> operator*(const Matrix &A, const vector<double> &x)
        {
            vector<double> y(A.nr);
            for(int r = 0; r < A.nr; r++)
            {
                y[r] = A.dotRow(r, x);
            }
            return y;
        }

        virtual double& operator()(int r, int c) = 0; // Used for setting and write access
        virtual double operator()(int r, int c) const = 0; // For read-only access

    protected:
        int nr = 0;
};

class DenseMatrix : public Matrix
{
    public:
        DenseMatrix(int nr) : Matrix(nr)
        {
            a = new double*[nr];
            for(int r = 0; r < nr; r++)
            {
                a[r] = new double[nr];
            }

            for(int r = 0; r < nr; r++)
                for(int c = 0; c < nr; c++)
                    a[r][c] = 0.0;
        }

        ~DenseMatrix()
        {
            for(int r = 0; r < nr; r++)
                delete[] a[r];

            delete[] a;
            a = nullptr; // Clear stored address of a
        }

        double dotRow(int r, const vector<double> &x) const
        {
            double sum = 0;
            for(int c = 0; c < nr; c++)
            {
                sum += a[r][c] * x[c];
            }
            return sum;
        }

        double& operator()(int r, int c)
        {
            return a[r][c];
        }
        double operator()(int r, int c) const
        {
            return a[r][c];
        }
    
        double **a = nullptr;
};

class IdentityMatrix : public Matrix
{
    public:
        IdentityMatrix(int nr) : Matrix(nr) {}

        double dotRow(int r, const vector<double> &x) const
        {
            return x[r];
        }

        double& operator()(int r, int c)
        {
            throw runtime_error("Unsupported Operation!");
        }
        double operator()(int r, int c) const
        {
            if(r == c)
                return 1;
            else 
                return 0;
        }
};

class FiveBandedMatrix : public Matrix
{
    public:
        FiveBandedMatrix(int nr, int d1) : Matrix(nr), d1{d1}
        {
            a = newAndClear(nr);
            b = newAndClear(nr);
            c = newAndClear(nr);
            d = newAndClear(nr);
            e = newAndClear(nr);
        } 

        ~FiveBandedMatrix()
        {
            delete[] a; a = nullptr;
            delete[] b; b = nullptr;
            delete[] c; c = nullptr;
            delete[] d; d = nullptr;
            delete[] e; e = nullptr;
        }

        double dotRow(int r, const vector<double> &x) const
        {
            double sum = 0;
            
            if(r-d1 >= 0)
                sum += a[r]*x[r-d1];
            if(r-1 >= 0) 
                sum += b[r]*x[r-1];
            sum += c[r]*x[r];
            if(r+1 < nr)
                sum += d[r]*x[r+1];
            if(r+d1 < nr) 
                sum += e[r]*x[r+d1];

            return sum;
        }

        double& operator()(int r, int c)
        {
            if(c - r == -d1) return a[r];
            else if(c - r == -1) return b[r];
            else if(c - r == 0) return this->c[r];
            else if(c - r == 1) return d[r];
            else if(c - r == d1) return e[r];
            else throw runtime_error("Unsupported operation!!");
        }

        double operator()(int r, int c) const
        {
            if(c - r == -d1) return a[r];
            else if(c - r == -1) return b[r];
            else if(c - r == 0) return this->c[r];
            else if(c - r == 1) return d[r];
            else if(c - r == d1) return e[r];
            // else throw runtime_error("Unsupported operation!");
            else return 0;
        }

    protected:
        static double *newAndClear(int nr)
        {
            double *p = new double[nr];
            for(int i = 0; i < nr; i++)
            {
                p[i] = 0;
            }
            return p;
        }

        double *a = nullptr;
        double *b = nullptr;
        double *c = nullptr;
        double *d = nullptr;
        double *e = nullptr;

        int d1;
};

vector<double> operator-(const vector<double> &a, const vector<double> &b)
{
    int n = a.size();
    vector<double> c(n);
    for(int i = 0; i < n; i++)
    {
        c[i] = a[i] - b[i];
    }
    return c;
}

double mag2(vector<double> &a)
{
    double sum = 0;
    for(double d:a)
    {
        sum += d*d;
    }
    return sum;
}

vector<double> solveGS(Matrix &A, vector<double> &g);

int main()
{
    // int ni = 51;   // number of nodes
  	// int nj = 26;

    int ni = 51;   // number of nodes
  	int nj = 26;

  	double x0 = 0.1;  // origin
  	double y0 = 0.1;

  	double dx = 0.01;  // cell spacing
  	double dy = 0.01;

  	int nn = ni*nj;  // total number of nodes

    int timesteps = 100;
    double dt = 0.01;

    double D = 0.1;
    double inlet_size = 0.04;
    double outlet_inner_radius = 0.1;
    double outlet_outer_radius = 0.2;

    vector<double> g(nn);

    FiveBandedMatrix A(nn, ni); // The spacing should be five because of the nature of this system
    FiveBandedMatrix lhs(nn, ni);

    // Fill in the A matrix and b vector. Note: dy is dr and dx is dz
    // j is indexing r, i is indexing z
    double z = 0.0;
    double r = 0.0;
    // vector<double> test(25);
    vector<double> n_vector(nj*ni);
    // double* dotprod = new double[25];
    for(int j = 0; j < nj; j++)
    {
        for(int i = 0; i < ni; i++)
        {
            int n = j*ni + i;
            z = (double)i*dx; // This is actually z
            r = (double)j*dy + dy/2; // This is actually rho/r

            n_vector[n] = 0; // All initial densities are 0 except for the boundary conditions

            // Neumann boundary conditions
            if(i == 0)
            {
                A(n,n) = 1.0/dx;
                A(n,n+1) = -1.0/dx;

                lhs(n,n) = -1.0/dx;
                lhs(n,n+1) = 1.0/dx;
            }
            else if(j == 0)
            {
                A(n,n) = 1.0/dy;
                A(n,n+ni) = -1.0/dy;

                lhs(n,n) = -1.0/dy;
                lhs(n,n+ni) = 1.0/dy;
            }
            else if(i == ni - 1)
            {
                A(n,n) = 1.0/dx;
                A(n,n-1) = -1.0/dx;

                lhs(n,n) = -1.0/dx;
                lhs(n,n-1) = 1.0/dx;
            }
            else if(j == nj - 1)
            {
                A(n,n) = 1.0/dy;
                A(n,n-ni) = -1.0/dy;

                lhs(n,n) = -1.0/dy;
                lhs(n,n-ni) = 1.0/dy;
            }
            else
            {
                A(n,n-ni) = -D*dt/(2*dx*dx);
                A(n,n-1) = -D*dt/(2*dy*dy);
                A(n,n) = 1 + D*dt/(2*r*dy) + D*dt/(dy*dy) + D*dt/(dx*dx);
                A(n,n+1) = -D*dt/(2*r*dy) - D*dt/(2*dy*dy);
                A(n,n+ni) = -D*dt/(2*dx*dx);

                lhs(n,n-ni) = D*dt/(2*dx*dx);
                lhs(n,n-1) = D*dt/(2*dy*dy);
                lhs(n,n) = 1 - D*dt/(2*r*dy) - D*dt/(dy*dy) - D*dt/(dx*dx);
                lhs(n,n+1) = D*dt/(2*r*dy) + D*dt/(2*dy*dy);
                lhs(n,n+ni) = D*dt/(2*dx*dx);

            }

            if(i == 0 && r <= inlet_size) // If at the inlet
            {
                n_vector[n] = 100.0;
            }
            if(i == ni-1 && r >= outlet_inner_radius && r <= outlet_outer_radius) // If at the outlet
            {
                n_vector[n] = 0.0;
            }

            // test[n] = 1;
        }
    }

    // for(int rows = 0; rows < 25; rows++)
    // {
    //     dotprod[rows] = A.dotRow(rows, test);
    // }

    // Print A for debugging
    // for(int j = 0; j < nn; j++)
    // {
    //     for(int i = 0; i < nn; i++)
    //     {
    //         // cout << "error is here" << endl;
    //         if(j-i == -ni || j-i == -1 || j-i == 0 || j-i == 1 || j-i == ni)
    //             cout << A(j,i) << ' ';
    //         else 
    //             cout << '0' << ' ';
    //     }
    //     cout << endl;
    // }

    for(int time_step = 0; time_step < timesteps; time_step++)
    {
        cout << "Progress: " << (double)time_step/(double)timesteps*100 << "%" << endl;
        vector<double> b(ni*nj);
        // Take the dot product of the original n vector and its corresponding matrix
        // to get the actual LHS vector b
        for(int rows = 0; rows < ni*nj; rows++)
        {
            b[rows] = lhs.dotRow(rows, n_vector);
        }

        // Solve the system for this time step and update the n vector
        n_vector = solveGS(A,b);

        // Keep applying the boundary conditions
        for(int j = 0; j < nj; j++)
        {
            for(int i = 0; i < ni; i++)
            {
                int n = j*ni + i;
                z = (double)i*dx; // This is actually z
                r = (double)j*dy + dy/2; // This is actually rho/r

                if(i == 0 && r <= inlet_size) // If at the inlet
                {
                    n_vector[n] = 100;
                }
                if(i == ni-1 && r >= outlet_inner_radius && r <= outlet_outer_radius) // If at the outlet
                {
                    n_vector[n] = 0;
                }

            }
        }
    }

    /* output vti file */
  	ofstream out("field.vti");

  	out<<"<VTKFile type=\"ImageData\">\n";
  	out<<"<ImageData WholeExtent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 "<<0<<"\""; 
  	out<<" Origin=\""<<x0<<" "<<y0<<" "<<0.0<<"\"";
  	out<<" Spacing=\""<<dx<<" " <<dy<<" "<<0.0<<"\">\n";
  	out<<"<Piece Extent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 "<<0<<"\">\n"; 
  	out<<"<PointData>\n";

	out<<"<DataArray Name=\"T\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int n=0;n<nn;n++) out<<n_vector[n]<<" ";	
	out<<"\n</DataArray>\n";

	out<<"</PointData>\n";
	out<<"</Piece>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";

	// free memory
	return 0;	// normal exit
}

vector<double> solveGS(Matrix &A, vector<double> &g)
{
    int nr = A.getNr();

    vector<double> x(nr);

    for(int it = 0; it < 10000; it++)
    {
        for(int r = 0; r < nr; r++)
        {
            double sum = A.dotRow(r,x) - A(r,r)*x[r];
            double x_star = (g[r] - sum) / A(r,r);

            x[r] += 1.4*(x_star - x[r]);
        }

        if(it % 50 == 0)
        {
            vector<double> R = A*x - g;

            double L2 = sqrt(mag2(R)/nr);

            cout << "Solver iteration: " << it << ", L2 norm: " << L2 << endl;
            if(L2 < 1e-4)
                break;
        } 
    }
    
    return std::move(x);
}