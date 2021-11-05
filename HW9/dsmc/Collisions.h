#ifndef _DSMC_H
#define _DSMC_H

#include "World.h"
#include "Species.h"

class DSMC_MEX {
public:
	DSMC_MEX(Species &species, World &world, int skip) :
		species{species}, world{world}, skip{skip} {

		// compute parameters for VHS collisions cross-section
		//reduced mass
		mr = species.mass*species.mass/(species.mass + species.mass);
		c[0] = 4.07e-10;
		c[1] = 0.77;
		c[2]= 2*Const::K*273.15/mr;	//Bird's reference params at 273.15 K
		c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)

	}

	int perform();		// runs DSMC

protected:
	void collide(double3 &vel1, double3 &vel2, double mass1, double mass2);
	double evalSigma(double g_rel) {
		return Const::PI*c[0]*c[0]*pow(c[2]/(g_rel*g_rel),c[1]-0.5)/c[3];
	}

	Species &species;
	World &world;
	int skip;
	double sigma_vr_max = 1e-14;	//some initial value for the (sigma*vr)_max parameter
	double mr;
	double c[4];			// sigma coefficients


};

#endif
