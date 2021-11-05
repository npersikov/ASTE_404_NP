#include <vector>
#include "Collisions.h"

using namespace std;

// performs DSMC
int DSMC_MEX::perform()
{

	// only perform every "skip" time steps
	if (world.getTs()%skip!=0) return -1;

	// sort particles to cell, we only do this when performing collisions
	species.sortParticlesToCells();

	double dt = world.getDt()*skip;

	double sigma_vr_max_temp = 0;	/*reset for max computation*/
	double3 dh = world.getDh();
	double dV = dh[0]*dh[1]*dh[2];		// cell volume
	double Fn = species.mpw;	//specific weight, using Bird's notation
	int num_cols=0;		//reset collision counter

	size_t num_cells = species.cell_part_ids.size();

	/*now perform collisions*/
	for (size_t c=0;c<num_cells;c++)
	{
		// get a list of particles in this cell
		const vector<int> &parts = species.cell_part_ids[c];
		int N = parts.size();
		if (N<2) continue;		// need at least two particles for collision

		// compute number of groups according to NTC
		double ng_f = 0.5*N*N*Fn*sigma_vr_max*dt/dV;
		int ng = (int)(ng_f+0.5);	// number of groups, round

		// loop over groups, assumes we have at least two particles per cell
		for (int g=0;g<ng;g++)
		{
			int p1, p2;
			p1 = (int)(rnd()*N);	// pick array index for the firt particle

			do {
				p2 = (int)(rnd()*N); // pick array index for particle 2
			} while (p2==p1);		 // try again if p1=p2

			// set references to the actual particles
			Particle &part1 =  species.getParticle(parts[p1]);
			Particle &part2 =  species.getParticle(parts[p2]);

			// compute relative velocity
			double3 vr_vec = part1.vel -  part2.vel;
			double vr = mag(vr_vec);

			// evaluate cross section
			double sigma = evalSigma(vr);

			// eval sigma_vr
			double sigma_vr=sigma*vr;

			// update sigma_vr_max_temp
			if (sigma_vr>sigma_vr_max_temp)
				sigma_vr_max_temp=sigma_vr;

			// evaluate probability
			double P=sigma_vr/sigma_vr_max;

			/*did the collision occur?*/
			if (P>rnd())
			{
				num_cols++;
				collide(part1.vel,part2.vel,species.mass, species.mass);
			}
		}
	}

	// update the sigma_vr_max parameter
	if (num_cols){
		sigma_vr_max = sigma_vr_max_temp;
	}

	return num_cols;

}


/* collides two particles*/
void DSMC_MEX::collide(double3 &vel1, double3 &vel2, double mass1, double mass2)
{
	double3 vcm = (mass1*vel1 + mass2*vel2)/(mass1+mass2);

	/*relative velocity, magnitude remains constant through the collision*/
	double3 vr = vel1 - vel2;
	double vr_mag = mag(vr);

	/*pick two random angles, per Bird's VHS method*/
	double cos_chi = 2*rnd()-1;
	double sin_chi = sqrt(1-cos_chi*cos_chi);
	double eps = 2*Const::PI*rnd();

	/*perform rotation*/
	vr[0] = vr_mag*cos_chi;
	vr[1] = vr_mag*sin_chi*cos(eps);
	vr[2] = vr_mag*sin_chi*sin(eps);

	/*post collision velocities*/
	vel1 = vcm + mass2/(mass1+mass2)*vr;
	vel2 = vcm - mass1/(mass1+mass2)*vr;
}

