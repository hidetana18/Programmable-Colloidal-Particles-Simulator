/*--- This Code is Written by Hidenori Tanaka in 2016 ---*/
/*--- Brownian Dynamics Simulation of Self-Replicating Colloidal Clusters in Two Dimensions ---*/
/*--- Integration Method: Overdamped Langevin Equation with Forward Euler Time Step ---*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/*--- Pointers to create a file ---*/
FILE    *traj;                  // Pointer for trajectory file
FILE    *numcl;                 // Pointer for counting newly formed clusters

/*--- Control parameters ---*/
#define N       1652            // total number of particles  N = Nl*(Nl-1)+12
#define L       82.0            // Length of a side of the simulation box
#define IE0     7.8125          // Steepness of harmonic core potential
#define IE1     10.0            // Interaction between attached monomers
#define IE2     14.0            // Interaction between particles inside monomers and clusteres
#define IE3     25.0            // Interaction between particles inside clusters (Irreversible bond)
#define Nc      3               // critical number of bonds required to be formed for cluster to fall apart.
#define Dim     2               // Dimensions
#define h       0.000002        // Time step
#define t_total 2000000000      // Total time steps
#define dat     100000          // Time interval to take trajectory data
#define NeighborListTime 1000   // Time interval creating new NeighborList
#define Dif     1.0             // Diffusion coefficient
#define Rho     80.0            // probability measure (described in Miranda PNAS paper)
#define Roff    1.05            // Cut off distance of interaction potential
#define BoL     1.05            // we say bond is created when distance between two particles is less than BoL
#define NeighborRadius 2.0      // verlet radius.
int IE1note;                    // Note for interaction strength between attached monomers
int IE2note;                    // Note for interaction strength between seed cluster and attached monomer
int Lnote;                      // Note for length of a side of box
int RunNnote;                   // Note, Nth run

/*----Declare functions----*/
void    InitialSet();           // Set initial conditions
void    SetParameters();        // Set initial parameters
void    SetIM();                // Set interaction matrix
void    Integrate();            // calculate positions and velocities of next time step for every paritlces
void    Note();                 // every 'dat' time, write down informations of positions of every particles
void    NewNeighborList();      // create new verlet list
void    CalculateDistances();   // calculate distance between every possible pairs of particles
void    ModifySV();             // modify species vector when type 0 particle attaches to cluster
void    ModifyCA();             // modify Connection vector
void    CalculateEdepth();      // calculate depth of potential energy given SV,CIN,CA
void    CalculateForce();       // calculate force acting between every possible pairs
void    SumForces();            // function to sum-up all forces acting between particles
void    Renew();                // put new variables into old variables preparing for integration at next time step
void    BondBreak();            // Assgin CIN if number of bonds between monomers reach Nc
int     visited[N];             // Store particles already visited in depth searching
void    Search();               // Functinon for depth search
int     CountBond();            // Count number of bonds within a cluster
int     NewCl[N];               // 
int     NumNewCl;               //
void    RemoveDuplicate();      // Remove duplicate numbers in an array
int     SV[N];                  // species vector (vector to remember which particle is assigned to be which species)
int     CIN[N];                 // cluster identification number (if a particle is not inside a cluster, -1)
int     CA[N][3];               // Connection Array, CA[][0]=>Complementary, CA[][1,2]=>Different category
int     NumCl;                  // total number of clusters
int     t;                      // global variable to make time roop
int     a,b;                    // Variables to iterate
int     BoN;                    // Number of bonds created
int     StartP;                 // Particles to start bond searching
double  RandNormal();           // generate random normal distributed numbers with mean=0, sigma(variance)=1
double  Edep;                   // Energy depth
double  rnew[N][Dim];           // -old=> input for integration,
double  rold[N][Dim];           // -new=> output of integration)
double  fpnew[N][Dim];          // Sum of potential forces acting on single particle
double  fpold[N][Dim];          //
double  F[N][N];                //
double  f[N][N][Dim];           // Each force acting between a pair of particle
double  D[N][N][Dim];           // D[a][b][0]=rnew[a][0]-rnew[b][0]
double  Rsq[N][N];              // DX^2+DY^2
double  R[N][N] ;               // sqrt(DX^2+DY^2)
double  IM[5][5];               // Interaction Matrix
double  CMIM[5][5];             // Interaction Matrix between a particle inside cluster and an attaced monomer
double  MMIM[5][5];             // Interaction Matrix between attaced monomers
int     NeighborList[N][N];     // NeighborList[i][] represents PIN of particles inside neighbor radius
int     NofNeighbors[N];        // This remembers how many particles are within the verlet radius
int     NPIC[1000];             // Number of particles inside a cluster

/*--- Cell List ---*/
#define SideOfCell 1.2
double  CellList[30][30][4];
void    SetCellList();
int     NinCell[30][30];
void    NewCellList();






















