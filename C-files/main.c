/*--- This Code is Written by Hidenori Tanaka in 2016 ---*/
/*--- Brownian Dynamics Simulation of Self-Replicating Colloidal Clusters in Two Dimensions ---*/
/*--- Integration Method: Overdamped Langevin Equation with Forward Euler Time Step ---*/

#include "main.h"
#include <time.h>

int main( int argc, char *argv[] ){
    
    int  Input;
    char *InputIN   = argv[1];
    sscanf(InputIN, "%d", &Input);
    //IE1     = IE1prime;
    //IE1note = IE1prime;
    RunNnote= Input;
    IE1note = IE1;
    IE2note = IE2;
    Lnote   = L;

    InitialSet();
    /* time integration loop */
    t=0;
    while(t < t_total){
        // Integrate Over-damped Langevin Equation with Forward Euler time-step
        Integrate();
        
        // Output trajectory data file
        Note();
        
        // Update neighbor list
        if(t%NeighborListTime==0){
            NewNeighborList();
        }

        // Loop to calculate parwise quantities
        for(a=0;a<N;a++){
            int c;
            for(c=0;c<NofNeighbors[a];c++){
                b=NeighborList[a][c];
            
                // Given configuration of two particles calculate DX, R[a][b] between them
                CalculateDistances();
                
                // Based on R[a][b], when type 0 attaches to cluster, we assgin complementally type
                ModifySV();
                
                // Based on R[a][b] we judge connections among particles inside clusters and monomers
                ModifyCA();
                
                // Calculate Depth of Interaction to be Assigned
                CalculateEdepth();
                
                // Calculate Force acting between Particles
                CalculateForce();
            
            }
        }
        //NewCellList();
        

            
        // Sum up all forces between a particle and others to get total force acting
        SumForces();
        
        //
        Renew();

        // Count Nc and check if the bond breaking is needed
        if(t%100==0){
            BondBreak();
        }
        
        // Increment total time steps by one
        t = t + 1;
        
    }
    return 0;
}




