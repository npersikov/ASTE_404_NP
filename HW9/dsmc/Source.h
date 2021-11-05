#ifndef _SOURCE_H
#define _SOURCE_H

#include "Species.h"

// samples N_sample particles on the zmin face with uniform velocity
class ColdBeamSource {
public:
	ColdBeamSource(Species &sp, World &world, int N_sample, double vel0) :
		sp{sp}, world{world}, N_sample{N_sample}, vel0{vel0} {
			x0 = world.getX0()[0];
			y0 = world.getX0()[1];
			Lx = world.getXd()[0]-x0;
			Ly = world.getXd()[1]-y0;
		}

	// samples N_sample randomly positioned particles
	void sample() {
		for (int p = 0; p < N_sample; p++) {
			double3 pos { x0 + rnd() * Lx,  y0 + rnd() * Ly,  0};  // random x, y, z=0
			double3 vel {0, 0, vel0};
			sp.addParticle(pos,vel);
		}
	}
protected:
	Species &sp;
	World &world;
	int N_sample;
	double vel0;
	double x0, y0;		// origin
	double Lx, Ly;		// domain length
};


#endif
