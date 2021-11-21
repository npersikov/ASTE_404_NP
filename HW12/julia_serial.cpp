#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <math.h>
#include <thread>
using namespace std;

void outputField(int **data, int ni, int nj);

// support for complex numbers
struct Complex {
	double r;
	double i;
 
	Complex(double a, double b) : r(a), i(b) { }

	double magnitude2() {return r*r + i*i;}
	
	Complex operator*(const Complex& a) {
		return Complex(r*a.r - i*a.i, i*a.r + r*a.i);
	}
	Complex operator+(const Complex& a) {
		return Complex(r + a.r, i + a.i);
	}
 };
 
   	
// computes value of julia set at [i,j]*/
int juliaValue(int i, int j, int ni, int nj)  {
	double fi = -1.0 + 2.0*i/ni;	// fi = [-1:1)
	double fj = -1.0 + 2.0*j/nj;	// fj = [-1:1)
	
	Complex c(-0.8, 0.156);	// coefficient for the image
	Complex a(fi, fj);		// pixel pos as a complex number

	int k;					// iteration counter
	for( k = 0; k < 200; k++) {
 		a = a * a + c;
		if (a.magnitude2() > 1000) break;	// check for divergence
	}
	return k;				// return 
}

void computePixels(int i, int j_0, int j_f, int ni, int nj, int** julia)
{
	for (int j=j_0;j<j_f;j++)
		julia[i][j] = juliaValue(i,j,ni,nj);
}

int main() {
	const int ni=4000;
	const int nj=4000;

	/*allocate memory for our domain*/
	int **julia=new int*[ni];
	for (int i=0;i<ni;i++) julia[i]=new int[nj];

	/*start timing*/
  	auto clock_start = chrono::high_resolution_clock::now();

	// compute pixels
	for (int i=0;i<ni;i++)
	{
		// for (int j=0;j<nj;j++)
		// 	julia[i][j] = juliaValue(i,j,ni, nj);
		thread t1(computePixels, i, 0, nj/12, ni, nj, julia);
		thread t2(computePixels, i, nj/12, nj/6, ni, nj, julia);
		thread t3(computePixels, i, nj/6, nj/4, ni, nj, julia);
		thread t4(computePixels, i, nj/4, nj/3, ni, nj, julia);
		thread t5(computePixels, i, nj/3, 5*nj/12, ni, nj, julia);
		thread t6(computePixels, i, 5*nj/12, nj/2, ni, nj, julia);
		thread t7(computePixels, i, nj/2, 7*nj/12, ni, nj, julia);
		thread t8(computePixels, i, 7*nj/12, 2*nj/3, ni, nj, julia);
		thread t9(computePixels, i, 2*nj/3, 3*nj/4, ni, nj, julia);
		thread t10(computePixels, i, 3*nj/4, 5*nj/6, ni, nj, julia);
		thread t11(computePixels, i, 5*nj/6, 11*nj/12, ni, nj, julia);
		thread t12(computePixels, i, 11*nj/12, nj, ni, nj, julia);

		t1.join();
		t2.join();
		t3.join();
		t4.join();
		t5.join();
		t6.join();
		t7.join();
		t8.join();
		t9.join();
		t10.join();
		t11.join();
		t12.join();
	}

	cout << std::thread::hardware_concurrency() << endl;

	/*capture ending time*/
  	auto clock_end = chrono::high_resolution_clock::now();

	outputField(julia,ni,nj);
 	
	std::chrono::duration<float> delta = clock_end-clock_start;
    cout << "Simulation took "<<delta.count()<< "s\n";

	// free memory
	for (int i=0;i<ni;i++) delete[] julia[i];
	delete[] julia;	

	return 0;
}

/*saves output in VTK format*/
void outputField(int **data, int ni, int nj) {
	stringstream name;
	name<<"julia.vti";

    /*open output file*/
    ofstream out(name.str());
   	if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}

	/*ImageData is vtk format for structured Cartesian meshes*/
	out<<"<VTKFile type=\"ImageData\">\n";
	out<<"<ImageData Origin=\""<<"0 0 0\" ";
	out<<"Spacing=\"1 1 1\" ";
	out<<"WholeExtent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 0\">\n";
	
	/*output data stored on nodes (point data)*/
	out<<"<PointData>\n";
	
	/*potential, scalar*/
	out<<"<DataArray Name=\"julia\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
	for (int j=0;j<nj;j++)	{
		for (int i=0;i<ni;i++) out<<data[i][j]<<" ";
		out<<"\n";
	}
	out<<"</DataArray>\n";

	/*close out tags*/
	out<<"</PointData>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
 	out.close();
}
	    
