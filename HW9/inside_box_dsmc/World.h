#ifndef _WORLD_H
#define _WORLD_H

#include <random>
#include <chrono>
#include "Vec.h"

/*object for sampling random numbers*/
class Rnd {
public:
	//constructor: set initial random seed and distribution limits
	Rnd(): mt_gen{/*std::random_device()()*/0}, rnd_dist{0,1.0} {}
	double operator() () {return rnd_dist(mt_gen);}

protected:
	std::mt19937 mt_gen;	    //random number generator
	std::uniform_real_distribution<double> rnd_dist;  //uniform distribution
};

extern Rnd rnd;		//tell the compiler that an object of type Rnd called rnd exists somewhere

/*define constants*/
namespace Const
{
	const double EPS_0 = 8.85418782e-12;  	// C/(V*m), vacuum permittivity
	const double QE = 1.602176565e-19;		// C, electron charge
	const double AMU = 1.660538921e-27;		// kg, atomic mass unit
	const double ME = 9.10938215e-31;		// kg, electron mass
	const double K = 1.380648e-23;			// J/K, Boltzmann constant
	const double PI = 3.141592653;			// pi
	const double EvToK = QE/K;				// 1eV in K ~ 11604
}

// computational domain
struct World {

	World(int ni, int nj, int nk): ni{ni}, nj{nj}, nk{nk}, nn{ni,nj,nk} {
		time_start =  std::chrono::high_resolution_clock::now();	//save starting time point
	}

	/* set class data, this-> allows us to distinguish
	 * the class member from the function argument
	 * cannot use initializer list outside the constructor	 */
	void setBoundingBox(const double3 &x0, const double3 &xd) {
		this->x0 = x0;
		this->xd = xd;
		int3 nn{ni,nj,nk};

		//set cell spacing
		dh = (xd-x0)/(nn-1);
	}
	void setTime(double dt, int max_ts) {this->dt = dt; this->max_ts=max_ts;}

	double3 getX0() {return x0;}
	double3 getXd() {return xd;}
	double3 getDh() {return dh;}

	// support for time advance
	int getTs() const {return ts;}  //returns current time step
	int getMaxTs() const {return max_ts;} // maximum number of steps
	double getDt() const {return dt;} // returns dt
	bool advanceTime() {ts++; return ts<=max_ts;}
	double getWallTime();  //returns wall time in seconds

	bool inBounds(double3 pos) {
		for (int i=0;i<3;i++)
			if (pos[i]<x0[i] || pos[i]>=xd[i]) return false;
		return true;
	}

	/*converts physical position to logical coordinate*/
	double3 XtoL(const double3 &x) const {
	  	return (x-x0)/dh;
	}

	/*converts physical position to a sequential cell index*/
	int XtoC(const double3 &x) const {
		double3 lc = XtoL(x);
		int i = (int)lc[0];
		int j = (int)lc[1];
		int k = (int)lc[2];
		return k*(ni-1)*(nj-1) + j*(ni-1) + i;
	}

	double3 X(double i, double j, double k) {
		return {x0[0]+i*dh[0], x0[1]+j*dh[1], x0[2]+k*dh[2]};
	}

	// sets sphere origin and radius
	void addSphere(const double3 &x0, double radius) {
		sphere_x0 = x0; sphere_r = radius;
	}

	/*return true if point x is inside or on the sphere*/
	bool inSphere(const double3 &x);

	/*returns the parametric position for the intersection point*/
	double lineSphereIntersect(const double3 &x1, const double3 &x2);

	/*returns random diffuse velocity sampled at a given surface point*/
	double3 sphereDiffuseVector(const double3 &x);

	const int ni, nj, nk; // number of mesh nodes
	const int3 nn;  // number of nodes stored in an int3


protected:
	double3 x0;
	double3 dh;
	double3 xd;		// diagonal corner from origin

	double dt = 0;
	int max_ts = 0;
	int ts = 0;

	double3 sphere_x0;		// sphere origin
	double sphere_r = 0;		// sphere radius

	std::chrono::time_point<std::chrono::high_resolution_clock> time_start;	//time at simulation start

};


#endif
