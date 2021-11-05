/*USC ASTE-404 1D PIC simulation of a plasma sheath
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <random>
using namespace std;

/*constants*/
#define EPS_0 8.85418782e-12  	// F/m, vacuum permittivity
#define K	1.38065e-23			// J/K, Boltzmann constant
#define ME 9.10938215e-31		// kg, electron mass
#define QE 1.602176565e-19		// C, elementary charge
#define AMU  1.660538921e-27	// kg, atomic mass unit
#define EV_TO_K	11604.52		// 1eV in Kelvin, QE/K

/*simulation parameters, these could come from an input file*/
#define PLASMA_DEN	1e16		// plasma density to load
#define NUM_IONS 5000000		    // number of ions
#define NUM_ELECTRONS 5000000	// number of electrons
#define DX 1e-4					// cell spacing
#define NC 100					// number of cells
#define NUM_TS	10000			// number of time steps
#define DT 1e-11				// time step size
#define ELECTRON_TEMP 3.0		// electron temperature in eV
#define ION_TEMP 1.0			// ion temperature in eV

/* Data structure to hold domain information*/
struct Domain
{
	int ni;				/*number of nodes*/
	double x0;			/*mesh origin*/
	double dx;			/*cell spacing*/
	double xl;			/*domain length*/
	double xmax;		/*domain max position*/
	
	/*data structures*/
	double *phi;		/*potential*/
	double *ef;			/*electric field*/
	double *rho;		/*charge density*/
	double *ndi;		/*ion density*/
	double *nde;		/*electron density*/
};

/* Data structure for particle storage **/
struct Particle
{
	double x;			/*position*/
	double v;			/*velocity*/
	bool alive;			/*flag to avoid removing particles*/
};

/* Data structure to hold species information*/
struct Species
{
	double mass;			/*particle mass in kg*/
	double charge;			/*particle charge in Coulomb*/
	double spwt;			/*species specific weight*/
	
	int np;					/*number of particles*/
	int np_alloc;			/*size of the allocated data array*/
	Particle *part;			/*array holding particles on the CPU*/
};

/** FUNCTION PROTOTYPES **/
double SampleVel(double v_th);
void ScatterSpecies(Species *species, double *nd);
void ComputeRho(Species *ions, Species *electrons);
bool SolvePotential(double *phi, double *rho);
bool SolvePotentialDirect(double *phi, double *rho);
void ComputeEF(double *phi, double *ef);
void PushSpecies(Species *species, double *ef);
void RewindSpecies(Species *species, double *ef);
void AddParticle(Species *species, double x, double v);	
double XtoL(double pos);
void scatter(double lc, double value, double *field);
double gather(double lc, double *field);
void WriteResults(int ts);

/* random number generator class*/
class Rnd {
public:
  std::random_device rnd_device;			// a true system-based slow RNG used to seed Mersenne Twister
  std::mt19937 mt_gen;      // random number generator
  std::uniform_real_distribution<double> rnd_dist;  //uniform dist in [0,1)
  double operator() () {  						// returns random number in [0,1)
	return rnd_dist(mt_gen);
 }
  Rnd() : mt_gen{rnd_device()} {}
};

/* GLOBAL VARIABLES */
Domain domain;
Rnd rnd;  			 // instance of random number generator

/* --------- main -------------*/
int main()
{
	int p;
	int ts;						// time step

	/*---1) initialize domain---*/
	domain.ni = NC+1;			// number of nodes
	domain.dx = DX;				// cell spacing
	domain.x0 = 0;					// origin
	domain.xl = (domain.ni-1)*domain.dx;//domain length
	domain.xmax = domain.x0+domain.xl;	//max position
	
	/*allocate data structures, remember these need to be cleared first in C++*/
	domain.phi = new double[domain.ni];		//potential
	domain.rho  = new double[domain.ni];		//charge density
	domain.ef  = new double[domain.ni];		//electric field
	domain.nde  = new double[domain.ni];		//number density of electrons
	domain.ndi	= new double[domain.ni];		//number density of ions
	
	/*save pointers so we can write phi instead of domain.phi*/
	double *phi = domain.phi;
	double *rho = domain.rho;
	double *ef = domain.ef;
	double *nde = domain.nde;
	double *ndi = domain.ndi;
	
	/*clear data, this sets 8*domain.ni bytes starting from &phi[0] to 0*/
	memset(phi,0,sizeof(double)*domain.ni);
	
	/*variables to hold species data ions*/
	Species ions;
	Species electrons;
	
	/*set ion data*/
	ions.mass = 16*AMU;		
	ions.charge = QE;
	ions.spwt = PLASMA_DEN*domain.xl/NUM_IONS;
	ions.np = 0;
	ions.np_alloc = NUM_IONS;
	ions.part = new Particle[NUM_IONS];
	
	/*set electron data*/
	electrons.mass = ME;	// electrons
	electrons.charge = -QE;
	electrons.spwt = PLASMA_DEN*domain.xl/NUM_ELECTRONS;
	electrons.np = 0;
	electrons.np_alloc = NUM_ELECTRONS;
	electrons.part = new Particle[NUM_ELECTRONS];
			
	/*load uniformly spaced ions and electrons*/
	double delta_ions = domain.xl/NUM_IONS;
	double v_thi = sqrt(2*K*ION_TEMP*EV_TO_K/ions.mass);	
	for (p=0;p<NUM_IONS;p++)
	{
		double x = domain.x0 + p*delta_ions;
		double v = SampleVel(v_thi);
		AddParticle(&ions,x,v);
	}
		
	/*now do the same for electrons*/
	double delta_electrons = domain.xl/NUM_ELECTRONS;
	double v_the = sqrt(2*K*ELECTRON_TEMP*EV_TO_K/electrons.mass);	
	for (p=0;p<NUM_ELECTRONS;p++)
	{
		double x = domain.x0 + p*delta_electrons;
		double v = SampleVel(v_the);
		AddParticle(&electrons,x,v);
	}
	
	/*compute number density*/
	ScatterSpecies(&ions,ndi);
	ScatterSpecies(&electrons,nde);
		
	/*compute charge density and solve potential*/
	ComputeRho(&ions, &electrons);
	SolvePotential(phi,rho);
	ComputeEF(phi,ef);
	
	RewindSpecies(&ions,ef);
	RewindSpecies(&electrons,ef);
	
	/*OUTPUT*/
	WriteResults(0);
		
	using namespace std::chrono;
    chrono::high_resolution_clock::time_point t1 = high_resolution_clock::now();  // grab starting clock time

 	/* MAIN LOOP*/
	for (ts = 1; ts<=NUM_TS; ts++)
	{
		/*compute number density*/
		ScatterSpecies(&ions,ndi);
		ScatterSpecies(&electrons,nde);
		
		ComputeRho(&ions, &electrons);
		SolvePotential(phi,rho);
		ComputeEF(phi,ef);
		
		/*move particles*/
		PushSpecies(&electrons,ef);
		PushSpecies(&ions,ef);
		
		/*write diagnostics*/
		if (ts%25==0)
		{
			/*max phi*/
			double max_phi = phi[0];
			for (int i=0;i<domain.ni;i++)
				if (phi[i]>max_phi) max_phi=phi[i];
			
			cout<<"TS:"<<ts<<"\t np_i:"<<ions.np<<"\t np_e:"<<electrons.np<<"\t dphi:"<<max_phi-phi[0]<<endl;
		}
		
		/*save data*/
		if (ts%1000==0) WriteResults(ts);
	}
	
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
    chrono::duration<double, std::milli> time_span = t2 - t1;

	/*free up memory*/
	delete[] phi;
	delete[] rho;
	delete[] ef;
	delete[] nde;
	delete[] ndi;
	
	/*free particles*/
	delete[] ions.part;
	delete[] electrons.part;
	
    cout<<"Time per time step: "<<time_span.count()/ts<< " ms"<<endl;

	return 0;
}

/***** HELPER FUNCTIONS *********************************************************/
  
/* samples random velocity from Maxwellian distribution using Birdsall's method*/
double SampleVel(double v_th)
{
    const int M = 12;    
    double sum = 0;
	for (int i=0;i<M;i++) sum+=rnd();

    return sqrt(0.5)*v_th*(sum-M/2.0)/sqrt(M/12.0);
}

/*scatter particles of species to the mesh*/
void ScatterSpecies(Species *species, double *den)
{
	/*initialize densities to zero*/
	memset(den,0,sizeof(double)*domain.ni);
	
	/*scatter particles to the mesh*/
	for (int p=0;p<species->np_alloc;p++)
	{
		if (!species->part[p].alive) continue;
		double lc = XtoL(species->part[p].x);
		scatter(lc, species->spwt, den);
	}
	
	/*divide by cell volume*/
	for (int i=0;i<domain.ni;i++)
		den[i]/=domain.dx;
		
	/*only half cell at boundaries*/
	den[0] *=2.0;
	den[domain.ni-1] *= 2.0;
}

/*adds new particle to the species, returns pointer to the newly added data*/
void AddParticle(Species *species, double x, double v)
{
	/*abort the simulation if we ran out of space to store this particle*/
	if (species->np > species->np_alloc-1)
	{
		printf("Too many particles!\n"); exit(-1);
	}
	
	/*store position and velocity of this particle*/
	species->part[species->np].x = x;
	species->part[species->np].v = v;
	species->part[species->np].alive = true;
		
	/*increment particle counter*/
	species->np++;
}

/*computes charge density by adding ion and electron data*/
void ComputeRho(Species *ions, Species *electrons)
{
	double *rho = domain.rho;
	
	for (int i=0;i<domain.ni;i++)
		rho[i] = ions->charge*domain.ndi[i] + electrons->charge*domain.nde[i];
}

/*Thomas algorithm for a tri-diagonal matrix*/
bool SolvePotentialDirect(double *x, double *rho)
{
	/*set coefficients, this should be pre-computed*/
	int ni = domain.ni;
	double dx2 = domain.dx*domain.dx;
	int i;
	double *a = new double[ni];
	double *b = new double[ni];
	double *c = new double[ni];
	
	/*central difference on internal nodes*/
	for (i=1;i<ni-1;i++)
	{
		a[i] = 1; b[i] = -2; c[i] = 1;
	}
	
	/*dirichlet b.c. on boundaries*/
	a[0]=0;b[0]=1;c[0]=0;
	a[ni-1]=0;b[ni-1]=1;c[ni-1]=0;
	
	/*multiply RHS*/
	for  (i=1;i<domain.ni-1;i++)
		x[i]=-rho[i]*dx2/EPS_0;
		
	x[0] = 0;
	x[ni-1] = 0;
	
	/* Modify the coefficients. */
	c[0] /= b[0];				/* Division by zero risk. */
	x[0] /= b[0];				/* Division by zero would imply a singular matrix. */
	for(i = 1; i < ni; i++)
	{
		double id = (b[i] - c[i-1] * a[i]);	/* Division by zero risk. */
		c[i] /= id;							/* Last value calculated is redundant. */
		x[i] = (x[i] - x[i-1] * a[i])/id;
	}

	/* Now back substitute. */
	for(i = ni - 2; i >= 0; i--)
		x[i] = x[i] - c[i] * x[i + 1];

	return true;		
}

/* solves potential using the Gauss Seidel Method, returns true if converged*/
bool SolvePotential(double *phi, double *rho)
{
	double L2;
	double dx2 = domain.dx*domain.dx;	/*precompute*/
	
	/*initialize boundaries*/
	phi[0]=phi[domain.ni-1]=0;
	
	/*solve potential, identical to lesson 2*/
	for (int solver_it=0;solver_it<40000;solver_it++)
	{
		/*Gauss Seidel method, phi[i-1]-2*phi[i]+phi[i+1] = -dx^2*rho[i]/eps_0*/
		for (int i=1;i<domain.ni-1;i++)
		{
			/*SOR*/
			double g = 0.5*(phi[i-1] + phi[i+1] + dx2*rho[i]/EPS_0);
			phi[i] = phi[i] + 1.4*(g-phi[i]);
		}
				
		/*check for convergence*/
		if (solver_it%25==0)
		{
			double sum = 0;
			for (int i=1;i<domain.ni-1;i++)
			{
				double R = -rho[i]/EPS_0 - (phi[i-1] - 2*phi[i] + phi[i+1])/dx2;
				sum+=R*R;
			}
			L2 = sqrt(sum)/domain.ni;
			if (L2<1e-4) {return true;}
		}
	}
	printf("Gauss-Seidel solver failed to converge, L2=%.3g!\n",L2);
	return false;	
}

/* computes electric field by differentiating potential*/
void ComputeEF(double *phi, double *ef)
{
	for (int i=1;i<domain.ni-1;i++)
		ef[i] = -(phi[i+1]-phi[i-1])/(2*domain.dx);	//central difference
	
	/*one sided difference at boundaries*/
	ef[0] = -(phi[1]-phi[0])/domain.dx;
	ef[domain.ni-1] = -(phi[domain.ni-1]-phi[domain.ni-2])/domain.dx;
}


/* moves particles of a single species*/
void PushSpecies(Species *species, double *ef)
{
	/*precompute q/m*/
	double qm = species->charge / species->mass;
	
	/*loop over particles*/
	for (int p=0;p<species->np_alloc;p++)
	{
		/*grab pointer to this particle*/
		Particle *part = &species->part[p];	
		
		if (!part->alive) continue;
			
		/*compute particle node position*/
		double lc = XtoL(part->x);
		
		/*gather electric field onto particle position*/
		double part_ef = gather(lc,ef);
		
		/*advance velocity*/
		part->v += DT*qm*part_ef;	
		
		/*advance position*/
		part->x += DT*part->v;
			
		/*remove particles leaving the domain*/
		if (part->x < domain.x0 || part->x >= domain.xmax)
		{
			/*mark particle as dead*/
			part->alive = false;
			species->np--;
			
		}
	}
}

/* rewinds particle velocities by -0.5DT*/
void RewindSpecies(Species *species, double *ef)
{
	/*precompute q/m*/
	double qm = species->charge / species->mass;
	
	for (int p=0;p<species->np_alloc;p++)
	{
		Particle *part = &species->part[p];	
			
		/*compute particle node position*/
		double lc = XtoL(part->x);
		
		/*gather electric field onto particle position*/
		double part_ef = gather(lc,ef);
		
		/*advance velocity*/
		part->v -= 0.5*DT*qm*part_ef;			
	}
}

/* converts physical coordinate to logical*/
double XtoL(double pos)
{
	double li = (pos-domain.x0)/domain.dx;
	return li;
}

/* scatters scalar value onto a field at logical coordinate lc*/
void scatter(double lc, double value, double *field)
{
	int i = (int)lc;
	double di = lc-i;
			
	field[i] += value*(1-di);
	field[i+1] += value*(di);
}

/* gathers field value at logical coordinate lc*/
double gather(double lc, double *field)
{
	int i = (int)lc;
	double di = lc-i;
				
	/*gather field value onto particle position*/
	double val = field[i]*(1-di) + field[i+1]*(di);
	return val;
}

/* writes new zone to the results file*/
void WriteResults(int ts)
{
	stringstream ss;
    ss<<"field_"<<setw(6)<<setfill('0')<<ts<<".csv";  // create file name
	ofstream out(ss.str());    // open file
    out<<"x,nd.e,nd.i,rho,phi,ef\n"<<endl;
	for (int i=0;i<domain.ni;i++)
	{
		out<<i*domain.dx<<","
		   <<domain.nde[i]<<","<<domain.ndi[i]<<","
           <<domain.rho[i]<<","<<domain.phi[i]<<","<<domain.ef[i]<<"\n";
	}			
}
