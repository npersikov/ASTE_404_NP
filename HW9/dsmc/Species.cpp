#include <math.h>
#include <iostream>
#include "Species.h"
#include "Field.h"
#include "Vec.h"

/*updates velocities and positions of all particles of this species*/
void Species::advance()
{

	size_t np = particles.size();
	/*loop over all particles*/
	for (size_t p=0;p<np;p++)
	{
		Particle &part  = particles[p];

		// assume no force
		double3 F{0,0,0};

		/*update velocity */
		part.vel += F;

		double part_dt = world.getDt();

		bool alive = true;

		while (part_dt>0 && alive) {
			double3 pos_old = part.pos;		// save position prior to the push
			part.pos += part.vel*part_dt;   // x_new = x_old + v*dt


			/*did this particle leave the domain?*/
			if (!world.inBounds(part.pos))
			{
				alive = false;	//kill the particle
			}
			else if (world.inSphere(part.pos)) {
				double tp = world.lineSphereIntersect(pos_old,part.pos);
				double dt_rem = (1-tp)*part_dt;
				part_dt -= dt_rem;

				//move particle *almost* to the surface
				part.pos = 	pos_old + 0.999*tp*(part.pos-pos_old);
				double v_mag1 = mag(part.vel);	//pre-impact speed
				part.vel = sampleReflectedVelocity(part.pos,v_mag1);
				continue;
			}
			else {part_dt = 0;} //use up all time step

			if (!alive) {
				particles[p] = particles[np-1]; //fill the hole
				np--;	//reduce count of valid elements
				p--;	//decrement p so this position gets checked again
			}
		}
	}

	//now delete particles[np:end]
	particles.erase(particles.begin()+np,particles.end());
}


/*returns random post-impact velocity*/
double3 Species::sampleReflectedVelocity(const double3 &pos, double v_mag1)
{
	double v_th = sampleVth(5000); //assume T_sphere = 1000K
	const double a_th = 1;		//thermal accommodation coeff
	double v_mag2 = v_mag1 + a_th*(v_th-v_mag1);
	return v_mag2*world.sphereDiffuseVector(pos); //set new velocity
}

/*adds a new particle, rewinding velocity by half dt*/
void Species::addParticle(double3 pos, double3 vel)
{
	//don't do anything (return) if pos outside domain bounds [x0,xd)
	if (!world.inBounds(pos)) return;

	double3 F{0,0,0};

	//rewind velocity back by 0.5*dt*F
    vel -=  (1/mass)*0.5*world.getDt()*F;

    //add to list
    particles.emplace_back(pos,vel);
}


/*returns random thermal velocity*/
double Species::sampleVth(double T)
{
	//thermal velocity
	double v_th = sqrt(2*Const::K*T/mass);
	//get three random velocity components
	double v1 = v_th*(rnd()+rnd()+rnd()-1.5);
	double v2 = v_th*(rnd()+rnd()+rnd()-1.5);
	double v3 = v_th*(rnd()+rnd()+rnd()-1.5);
	return 3/sqrt(2+2+2)*sqrt(v1*v1+v2*v2+v3*v3);	//magnitude
}

/*returns random isotropic velocity*/
double3 Species::sampleIsotropicVel(double T) {
	double theta = 2*Const::PI*rnd();
    double r = -1.0+2*rnd();    //pick a random direction for d[0]
    double a = sqrt(1-r*r); //scaling for unity magnitude

    double3 d;
    d[0] = r;
	d[1] = cos(theta)*a;
    d[2] = sin(theta)*a;

    double v_th = sampleVth(T);
    double3 vel = v_th*d;
    return vel;
}

//collect data for averaging
void Species::sampleMoments() {

	for (Particle &part:particles)
	{
		double3 lc = world.XtoL(part.pos);
		count_sum.scatter(lc, mpw);
		vel_sum.scatter(lc,mpw*part.vel);
	}

	num_samples++;
}

// computes density and velocity
void Species::computeGasProperties() {

	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
			for (int k=0;k<world.nk;k++) {
				den[i][j][k] = count_sum[i][j][k]/num_samples;
				if (count_sum[i][j][k]>0)
					vel[i][j][k] = vel_sum[i][j][k]/count_sum[i][j][k];
				else vel[i][j][k] = 0;

			}
	num_samples = 0;

	count_sum.clear();
	vel_sum.clear();


	//divide by node volume
	double3 dh = world.getDh();
	for (int k=0;k<world.nk;k++)
		for (int j=0;j<world.nj;j++)
			for (int i=0;i<world.ni;i++) {

				// volume of regular internal node control volume
				double dV=dh[0]*dh[1]*dh[2];

				// only partial volumes on domain boundaries, half on faces, 1/4 on edges, 1/8 on corners
				if (i==0 || i==world.ni-1) dV/=2.0;
				if (j==0 || j==world.nj-1) dV/=2.0;
				if (k==0 || k==world.nk-1) dV/=2.0;
				den[i][j][k] /= dV;
			}
}

// sorts particles into cells
void Species::sortParticlesToCells() {

	int num_cells = (world.ni-1)*(world.nj-1)*(world.nk-1);
	cell_part_ids.resize(num_cells);

	// clear past data
	for (int c=0;c<num_cells;c++)
		cell_part_ids[c].clear();

	// sort particles to cells
	for (size_t p=0;p<particles.size();p++)
	{
		Particle &part = particles[p];
		int c = world.XtoC(part.pos);
		cell_part_ids[c].push_back(p);
	}

}

