#ifndef _SPECIES_H
#define _SPECIES_H
#include <string>
#include <vector>
#include "Vec.h"
#include "Field.h"
#include "World.h"

// single particle container
struct Particle {
	Particle(const double3 &pos, const double3 &vel) : pos{pos}, vel{vel} {}

	double3 pos;	// particle position
	double3 vel;	// particle velocity
};

// container for all particles of the same chemical type
class Species {
public:
	Species(std::string name, double mass, double mpw0, World &world) :
			name(name), mass(mass), mpw(mpw0),
			den(world.nn), vel(world.nn), world(world),
			count_sum(world.nn), vel_sum(world.nn) {		}

	const std::string name;
	double mass;	 // physical mass of a single particle (kg)
	double mpw;		 // macroparticle weight, number of real particles per sim particle

	Field den;
	Field3 vel;

	void advance();  //advances particle positions and velocities through dt
	void sampleMoments();
	void computeGasProperties();	// computes density and velocity

	/*adds a new particle*/
	void addParticle(double3 pos, double3 vel);

	// samples random "thermal velocity" magnitude
	double sampleVth(double T);
	// samples random vector for isotropic thermal velocity at temperature T
	double3 sampleIsotropicVel(double T);
	double3 sampleReflectedVelocity(const double3 &pos, double v_mag);

	size_t getNp() {return particles.size();} // number of particles
	Particle& getParticle(int p) {return particles[p];}


	/* DSMC support */
	void sortParticlesToCells();
	std::vector<std::vector<int>> cell_part_ids; //cell_part_ids[c] contains a vector of all particles in cell c

protected:
	World &world;				// reference to world
	std::vector<Particle> particles;  // particle array (contiguous memory)

	Field count_sum;  // sum of particle counts
	Field3 vel_sum;
	int num_samples = 0;

};


#endif
