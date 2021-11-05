/*USC ASTE-404 example DSMC code for flow past a sphere */

#include <iostream>
#include <string>
#include <vector>
#include "Vec.h"
#include "World.h"
#include "Output.h"
#include "Source.h"
#include "Collisions.h"

using namespace std;

//create the rnd functor
Rnd rnd;


int main() {
	// set some simulation parameters
	const double nd0 = 1e21;		// desired injection number density
	const double vel0 = 10000;		// injection velocity
	const int N_sample = 1000;		// number of particles to create

	World world{41,41,61};		// 41x41x81 mesh
	world.setBoundingBox({-0.3,-0.3,0}, {0.3,0.3,1.0});   //set origin and diagonal corner

	// set time step such that particles move 0.1*dz at vel0
	double dz = world.getDh()[2];	// get cell spacing in z
	world.setTime(0.1*dz/vel0, 5000);			// time step and max number of time steps

	cout<<"Running for "<<world.getMaxTs()<<" "<<world.getDt()<<" time steps"<<endl;

	world.addSphere({0.0, 0.0, 0.25}, 0.05);

	// Lx * Ly
	double3 xd = world.getXd();
	double3 x0 = world.getX0();
	double A_source = (xd[0]-x0[0])*(xd[1]-x0[1]);

	// set macroparticle weight such that we create 100 particles per time step
	double n_dot = nd0*vel0*A_source;		// number of real particles created per second
	double n_real = n_dot*world.getDt();	// number of real particles created per time step
	double mpw = n_real/N_sample;			// set mpw such that N_create particles are generated

	// create species container for N2 neutrals
	Species neutrals("N2", 14*Const::AMU, mpw, world);  // container for neutrals

	// create a cold beam source on the z=0 face, samples N_sample particles with velocity vel0
	ColdBeamSource source(neutrals,world,N_sample,vel0);

	// initialize DSMC,
	/* collisions will be performed only once every xx time steps
	 * this "subcycling" helps reduce run time but is sometimes also helpful when colision rate is too low*/
	DSMC_MEX dsmc_mex(neutrals, world, 2);
	int num_cols = 0;

	// main loop, using a while loop instead of for
	while(world.advanceTime()) {
		source.sample();		//injects particles

		neutrals.advance();		//integrate velocity and position
		neutrals.sampleMoments();

		num_cols = dsmc_mex.perform();

		// periodically save results
		if (world.getTs()%100==0) {
			cout<<"ts="<<world.getTs()<<", np="<<neutrals.getNp()<<", num_cols_inst="<<num_cols<<endl;
			neutrals.computeGasProperties();
			Output::saveVTI(world, neutrals);
		}
	}

	/* grab ending time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;
}
