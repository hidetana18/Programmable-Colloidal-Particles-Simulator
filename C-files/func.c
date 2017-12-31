/*--- This Code is Written by Hidenori Tanaka in 2016 ---*/
/*--- Brownian Dynamics Simulation of Self-Replicating Colloidal Clusters in Two Dimensions ---*/
/*--- Integration Method: Overdamped Langevin Equation with Forward Euler Time Step ---*/

#include "main.h"
#include "Mt.h"
#include <time.h>

void InitialSet(){
    // Initialize random number generator with current time
    init_genrand((unsigned)time(NULL)+10000*RunNnote);

    /* set initial conditions */
    SetParameters();

    SetIM();
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void SetParameters(){
    int i;
    
    // SV, CIN, CA[][0], CA[][1], CA[][2], rold[][0], rold[][1]
    for(i=0;i<N;i++){
        SV[i]=0, CIN[i]=-1, CA[i][0]=-1, CA[i][1]=-1, CA[i][2]=-1;
        
        rold[i][0] = L*(genrand_real3()-0.5) , rold[i][1] = L*(genrand_real3()-0.5);
        
        // Make sure monomer doesn't overlap with Initial Parent Clusters
        while ( -4.5< rold[i][0] && rold[i][0] < 4.5 && -4.5 < rold[i][1] && rold[i][1] < 4.5 ){
            rold[i][0] = L*(genrand_real3()-0.5),rold[i][1] = L*(genrand_real3()-0.5);
        }
    }
    

    
   NumCl   = 4;                  // NumCl(Number of Clusters) => total number of Clusters formed.
   NPIC[0] = 4;                  // NPIC (Number of Particle inside a Cluster)
   NPIC[1] = 4;
   NPIC[2] = 4;
   NPIC[3] = 4;
    
    
   SV[0] = 1, CIN[0] = 0, CA[0][0] =-1, CA[0][1] = 1, CA[0][2] = 3, rold[0][0] = -2.0, rold[0][1] =  2.0;
   SV[1] = 2, CIN[1] = 0, CA[1][0] =-1, CA[1][1] = 0, CA[1][2] = 2, rold[1][0] = -1.0, rold[1][1] =  2.0;
   SV[2] = 1, CIN[2] = 0, CA[2][0] =-1, CA[2][1] = 1, CA[2][2] = 3, rold[2][0] = -1.0, rold[2][1] =  1.0;
   SV[3] = 2, CIN[3] = 0, CA[3][0] =-1, CA[3][1] = 0, CA[3][2] = 2, rold[3][0] = -2.0, rold[3][1] =  1.0;
    
   SV[4] = 1, CIN[4] = 1, CA[4][0] =-1, CA[4][1] = 5, CA[4][2] = 7, rold[4][0] =  1.0, rold[4][1] =  2.0;
   SV[5] = 2, CIN[5] = 1, CA[5][0] =-1, CA[5][1] = 4, CA[5][2] = 6, rold[5][0] =  2.0, rold[5][1] =  2.0;
   SV[6] = 1, CIN[6] = 1, CA[6][0] =-1, CA[6][1] = 5, CA[6][2] = 7, rold[6][0] =  2.0, rold[6][1] =  1.0;
   SV[7] = 2, CIN[7] = 1, CA[7][0] =-1, CA[7][1] = 4, CA[7][2] = 6, rold[7][0] =  1.0, rold[7][1] =  1.0;
    
   SV[8] = 3, CIN[8] = 2, CA[8][0] =-1, CA[8][1] = 9, CA[8][2] =11, rold[8][0] = -2.0, rold[8][1] = -1.0;
   SV[9] = 4, CIN[9] = 2, CA[9][0] =-1, CA[9][1] = 8, CA[9][2] =10, rold[9][0] = -1.0, rold[9][1] = -1.0;
   SV[10]= 3, CIN[10]= 2, CA[10][0]=-1, CA[10][1]= 9, CA[10][2]=11, rold[10][0]= -1.0, rold[10][1]= -2.0;
   SV[11]= 4, CIN[11]= 2, CA[11][0]=-1, CA[11][1]= 8, CA[11][2]=10, rold[11][0]= -2.0, rold[11][1]= -2.0;
    
   SV[12]= 3, CIN[12]= 3, CA[12][0]=-1, CA[12][1]=13, CA[12][2]=15, rold[12][0]=  1.0, rold[12][1]= -1.0;
   SV[13]= 4, CIN[13]= 3, CA[13][0]=-1, CA[13][1]=12, CA[13][2]=14, rold[13][0]=  2.0, rold[13][1]= -1.0;
   SV[14]= 3, CIN[14]= 3, CA[14][0]=-1, CA[14][1]=13, CA[14][2]=15, rold[14][0]=  2.0, rold[14][1]= -2.0;
   SV[15]= 4, CIN[15]= 3, CA[15][0]=-1, CA[15][1]=12, CA[15][2]=14, rold[15][0]=  1.0, rold[15][1]= -2.0;
    
    
    
    
   /*---
    
    
    NumCl   = 5;                  // NumCl(Number of Clusters) => total number of Clusters formed.
    NPIC[0] = 5;
    NPIC[1] = 5;
    NPIC[2] = 5;
    NPIC[3] = 5;
    
    
    SV[0] = 1, CIN[0] = 0, CA[0][0] =-1, CA[0][1] =-1, CA[0][2] = 1, rold[0][0] = -3.0, rold[0][1] =  2.0;
    SV[1] = 2, CIN[1] = 0, CA[1][0] =-1, CA[1][1] = 0, CA[1][2] = 2, rold[1][0] = -3.0, rold[1][1] =  1.0;
    SV[2] = 1, CIN[2] = 0, CA[2][0] =-1, CA[2][1] = 1, CA[2][2] = 3, rold[2][0] = -3.0, rold[2][1] =  0.0;
    SV[3] = 2, CIN[3] = 0, CA[3][0] =-1, CA[3][1] = 2, CA[3][2] = 4, rold[3][0] = -3.0, rold[3][1] = -1.0;
    SV[4] = 1, CIN[4] = 0, CA[4][0] =-1, CA[4][1] = 3, CA[4][2] =-1, rold[4][0] = -3.0, rold[4][1] = -2.0;
    
    SV[5] = 2, CIN[5] = 1, CA[5][0] =-1, CA[5][1] =-1, CA[5][2] = 6, rold[5][0] = -1.0, rold[5][1] =  4.0;
    SV[6] = 1, CIN[6] = 1, CA[6][0] =-1, CA[6][1] = 5, CA[6][2] = 7, rold[6][0] = -1.0, rold[6][1] =  3.0;
    SV[7] = 2, CIN[7] = 1, CA[7][0] =-1, CA[7][1] = 6, CA[7][2] = 8, rold[7][0] = -1.0, rold[7][1] =  2.0;
    SV[8] = 1, CIN[8] = 1, CA[8][0] =-1, CA[8][1] = 7, CA[8][2] = 9, rold[8][0] = -1.0, rold[8][1] =  1.0;
    SV[9] = 2, CIN[9] = 1, CA[9][0] =-1, CA[9][1] = 8, CA[9][2] =-1, rold[9][0] = -1.0, rold[9][1] =  0.0;
    
    SV[10]= 3, CIN[10]= 2, CA[10][0]=-1, CA[10][1]=-1, CA[10][2]=11, rold[10][0]=  1.0, rold[10][1]=  0.0;
    SV[11]= 4, CIN[11]= 2, CA[11][0]=-1, CA[11][1]=10, CA[11][2]=12, rold[11][0]=  1.0, rold[11][1]= -1.0;
    SV[12]= 3, CIN[12]= 2, CA[12][0]=-1, CA[12][1]=11, CA[12][2]=13, rold[12][0]=  1.0, rold[12][1]= -2.0;
    SV[13]= 4, CIN[13]= 2, CA[13][0]=-1, CA[13][1]=12, CA[13][2]=14, rold[13][0]=  1.0, rold[13][1]= -3.0;
    SV[14]= 3, CIN[14]= 2, CA[14][0]=-1, CA[14][1]=13, CA[14][2]=-1, rold[14][0]=  1.0, rold[14][1]= -4.0;
    
    SV[15]= 4, CIN[15]= 3, CA[15][0]=-1, CA[15][1]=-1, CA[15][2]=16, rold[15][0]=  3.0, rold[15][1]=  2.0;
    SV[16]= 3, CIN[16]= 3, CA[16][0]=-1, CA[16][1]=15, CA[16][2]=17, rold[16][0]=  3.0, rold[16][1]=  1.0;
    SV[17]= 4, CIN[17]= 3, CA[17][0]=-1, CA[17][1]=16, CA[17][2]=18, rold[17][0]=  3.0, rold[17][1]=  0.0;
    SV[18]= 3, CIN[18]= 3, CA[18][0]=-1, CA[18][1]=17, CA[18][2]=19, rold[18][0]=  3.0, rold[18][1]= -1.0;
    SV[19]= 4, CIN[19]= 3, CA[19][0]=-1, CA[19][1]=18, CA[19][2]=-1, rold[19][0]=  3.0, rold[19][1]= -2.0;

    
    ---*/

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void SetIM(){
 
//%%%%%%%%%%%%%%%%% Interaction Matrix %%%%%%%%%%%%%%%%%%%%//
/* Interaction Matrix between given pairs of particles */
    
//            type 0         A            B             A'           B'      //
/* 0  */    IM[0][0]= 0; IM[0][1]= 1; IM[0][2]= 1; IM[0][3]= 1; IM[0][4]= 1;
/* A  */    IM[1][0]= 1; IM[1][1]= 0; IM[1][2]= 1; IM[1][3]= 1; IM[1][4]= 0;
/* B  */    IM[2][0]= 1; IM[2][1]= 1; IM[2][2]= 0; IM[2][3]= 0; IM[2][4]= 1;
/* A' */    IM[3][0]= 1; IM[3][1]= 1; IM[3][2]= 0; IM[3][3]= 0; IM[3][4]= 1;
/* B' */    IM[4][0]= 1; IM[4][1]= 0; IM[4][2]= 1; IM[4][3]= 1; IM[4][4]= 0;
    

//%%%%%%%%%%%%%%%%% Cluster Monomer Interaction Matrix %%%%%%%%%%%%%%%%%%%%//
/* Interaction Matrix between given pairs of particles */

    
//              type 0          A              B              A'             B'    //
/* 0  */    CMIM[0][0]= 1; CMIM[0][1]= 1; CMIM[0][2]= 1; CMIM[0][3]= 1; CMIM[0][4]= 1;
/* A  */    CMIM[1][0]= 1; CMIM[1][1]= 0; CMIM[1][2]= 0; CMIM[1][3]= 1; CMIM[1][4]= 0;
/* B  */    CMIM[2][0]= 1; CMIM[2][1]= 0; CMIM[2][2]= 0; CMIM[2][3]= 0; CMIM[2][4]= 1;
/* A' */    CMIM[3][0]= 1; CMIM[3][1]= 1; CMIM[3][2]= 0; CMIM[3][3]= 0; CMIM[3][4]= 0;
/* B' */    CMIM[4][0]= 1; CMIM[4][1]= 0; CMIM[4][2]= 1; CMIM[4][3]= 0; CMIM[4][4]= 0;
    

//%%%%%%%%%%%%%%%%% Monomer Monomer Interaction Matrix %%%%%%%%%%%%%%%%%%%%//
/* Interaction Matrix between given pairs of particles */

    
    //          type 0          A              B              A'             B'    //
/* 0  */    MMIM[0][0]= 0; MMIM[0][1]= 0; MMIM[0][2]= 0; MMIM[0][3]= 0; MMIM[0][4]= 0;
/* A  */    MMIM[1][0]= 0; MMIM[1][1]= 0; MMIM[1][2]= 1; MMIM[1][3]= 0; MMIM[1][4]= 0;
/* B  */    MMIM[2][0]= 0; MMIM[2][1]= 1; MMIM[2][2]= 0; MMIM[2][3]= 0; MMIM[2][4]= 0;
/* A' */    MMIM[3][0]= 0; MMIM[3][1]= 0; MMIM[3][2]= 0; MMIM[3][3]= 0; MMIM[3][4]= 1;
/* B' */    MMIM[4][0]= 0; MMIM[4][1]= 0; MMIM[4][2]= 0; MMIM[4][3]= 1; MMIM[4][4]= 0;
    

}









//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// Integrate over-damped langevin equation to calculate position at next time step
void Integrate(){
    int i,j;
    for(i=0;i<N;i++){
        for(j=0;j<Dim;j++){
            // Over-damped Langevin equation with Forward Euler Time step
            rnew[i][j]  =  rold[i][j]  + fpold[i][j]*h + sqrt(2.0*Dif*h)*RandNormal();
            // Periodic boundary condition
            rnew[i][j]  =  rnew[i][j] - round(rnew[i][j]/L)*L;
        }
    }
}


//******************** Cell List Constructor ********************//

void NewCellList(){
    int NofCell;       // Number of Cells
    double SideofCell; // Length of side of a cell
    int i,j,k;         // Loop index
    // Vector Cell Index
    int VCI[2];
    // Scalar Cell Index
    int SCI;
    // Number of Cells in each direction
    NofCell  = round(L/Roff);
    // Length of side of an each cell
    SideofCell = L/NofCell;
    // Define the headers
    int head[NofCell*NofCell];
    int CellList[N];
    // Reset the headers
    for(i=0;i<NofCell*NofCell;i++){
        head[i]=-1;
    }
    // Scan particles to construct headers
    for(j=0;j<N;j++){
        // Calculate Vector Cell Index to which atom "j" belongs
        for(k=0;k<Dim;k++){
            VCI[k] = floor((R[j][k]+L/2)/SideofCell);
        }
        // Translate Vector Cell Index to Scaller cell index
        SCI = VCI[0]*NofCell+VCI[1];
        /* Link to the previous occupant (or EMPTY if you're the 1st) */
        CellList[i] = head[SCI];
        /* The last one goes to the header */
        head[SCI] = j;
        SCI=0;
    }
}


/*---



//%%%%%%%%%%%%%%%% Set Cell List %%%%%%%%%%%%%%//
void SetCellList(){
    
    int k,l;
    for( k=0;k<30;k++ ){
        for( l=0;l<30;l++){
            NinCell[k][l]=0;
        }
    }
    
    
    int i, CellX,CellY;
    for(i=0;i<N;i++){
        CellX = (int)(floor( rnew[i][0]/SideOfCell )+15);
        CellY = (int)(floor( rnew[i][1]/SideOfCell )+15);
        
        CellList[ CellX ][ CellY ][ NinCell[ CellX ][ CellY ] ] = i;
        NinCell[ CellX ][ CellY ] +=1 ;
    }
}



for( i=0;i<30;i++ ){
    for(j=0;j<30;j++){
        if(i+1=){
            i=
        }
    
        
        
        for (k=0; l<NinCell[i][j]-1){
            a=CellList[i][j][k];
            for(l=k+1;l<NinCell[i][j]){
            b=CellList[i][j][l]
                
        
            }
            for(l=0;l<NinCell[i][j]){
                b=CellList[i+1][j][l]
                
                
            }
            for(l=0;l<NinCell[i][j]){
                b=CellList[i][j+1][l]
                
                
            }
            for(l=0;l<NinCell[i][j]){
                b=CellList[i+1][j+1][l]
                
                
            }
            for(l=0;l<NinCell[i][j]){
                b=CellList[i+1][j-1][l]
                
                
            }
        }
    }
}

---*/




//%%%%%%%%%%%%%%%%%%%%%%% Output trajectory file %%%%%%%%%%%%%%%%%%%%%%%%%%//
void Note(){
    int i,j;
    if(t%dat == 0){
        char trajfile[60];
        sprintf(trajfile,"trajNc3N%dL%dMM%dCM%dRun%d.dat",N,Lnote,IE1note,IE2note,RunNnote);
        // Open the file
        traj     =   fopen(trajfile,"a");
        fprintf(traj, "%d\n%s\n",N,"empty");
        // Print SV, CA[][0], CA[][1], CA[][2], CIN[]
        for(i=0;i<N;i++){
            fprintf(traj,"%d\t%d\t%d\t%d\t%d\t%d\t",i, SV[i] ,CA[i][0],CA[i][1],CA[i][2],CIN[i]);
            for(j=0;j<Dim;j++){
                fprintf(traj,"%f\t", rnew[i][j] );
            }
            fprintf(traj,"\n");
        }
        fclose(traj);
    }
}


// NewNeighborList => Given rnew[PIN][Dim] of every particles, create Neighbor list [][]
void NewNeighborList(){
    int i,j,k,l;
    // Reset number of neighbors
    for(i=0;i<N;i++){
        NofNeighbors[i]=0;
    }
    
    // Loop through possible pairs "k and j"
    for(j=0;j<N-1;j++){
        for(k=j+1;k<N;k++){
            
            // Reset distance squared
            Rsq[j][k] = 0;
            
            // Calculate each component (l) of distance
            for(l=0;l<Dim;l++){
                D[j][k][l] = rnew[j][l]  -  rnew[k][l];
                D[j][k][l] = D[j][k][l]  -  round(D[j][k][l]/L)*L;
                Rsq[j][k]  = Rsq[j][k]   +  D[j][k][l]*D[j][k][l];
            }
            // Calculate distance
            R[j][k]  = sqrt(Rsq[j][k]);
            R[k][j]  = R[j][k];
            
            // If the pair is within NeighborRadius
            if(R[j][k]<NeighborRadius){
                NeighborList[j][NofNeighbors[j]]=k;
                NeighborList[k][NofNeighbors[k]]=j;
                NofNeighbors[j]=NofNeighbors[j]+1;
                NofNeighbors[k]=NofNeighbors[k]+1;
            }
        }
    }
}






//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// function to calculate distance R[a][b] between particles 'a' and 'b' and force which acts between them
void CalculateDistances(){
    int i;
    Rsq[a][b]  = 0;
    for(i=0;i<Dim;i++){
        D[a][b][i] = rnew[a][i]  - rnew[b][i];
        D[a][b][i] = D[a][b][i]  - round(D[a][b][i]/L)*L;
        // DX^2 + DY^2
        Rsq[a][b]  = Rsq[a][b]   + D[a][b][i]*D[a][b][i];
    }
    // (DX^2+DY^2)^(1/2)
    R[a][b]  = sqrt( Rsq[a][b] );
    R[b][a]  = R[a][b];
}


// Species Vector is modified only when "type 0" attaches to an "unoccupied" particle inside cluster//
void ModifySV(){
    if(R[a][b]<=BoL){
        if(SV[a]==0/*a is type 0 <=> unoccupied */ && CIN[b]!= -1 /* b is inside cluster */ && CA[b][0] == -1 /* b is unoccupied */ ){
            if(SV[b]==1){
                SV[a]=3;
            }
            else if(SV[b]==2){
                SV[a]=4;
            }
            else if(SV[b]==3){
                SV[a]=1;
            }
            else if(SV[b]==4){
                SV[a]=2;
            }
                    
                    /* make the species for the one of complementally particle*/
                                                                                                                                            }
        else if(SV[b]==0 /* b is type 0 */ && CIN[a]!= -1 /* a is inside cluster */ && CA[a][0] == -1 /* a is unoccupied */ ){
            if(SV[a]==1){
                SV[b]=3;
            }
            else if(SV[a]==2){
                SV[b]=4;
            }
            else if(SV[a]==3){
                SV[b]=1;
            }
            else if(SV[a]==4){
                SV[b]=2;
            }
        }
    }
}




//
//CA[0][0] CA[0][1] CA[0][2];
//CA[1][0] CA[1][1] CA[1][2];

//Modify Connection Array when "unoccupied" "monomer and cluster" comes togeter or "connected" particles fallen apart
// so this time (we ignoer monomer monomer interaction)
// for A-A', B-B' bonding we allow it only between particles "inside cluster" and "free monomer"

// NO CONNECTION FOR
// MONOMER-MONOMER
// MONOMER-ATTACHED MONOMER
// DIFFERENT CLUSTER

// ALLOW CONNECTION FOR
// MONOMER-CLUSTER
// ATTACHED MONOMER-ATTACHED MONOMER
// ATTACHED MONOMER-CLUSTER
// SAME CLUSTER-CLUSTER

void ModifyCA(){
    
// CONNECT
// If distance between a pair of particles is less than "Bond Length"
    if(R[a][b]<=BoL){
         
// If the pair is COMPLEMENTARY, A-A', B-B'
         if(CMIM[SV[a]][SV[b]] == 1){
             // If both particles have spot for complementary particle
             if(CA[a][0] == -1 && CA[b][0] == -1){
                 // If a is monomer and b is in cluster
                 if( CIN[a] == -1 && CIN[b] != -1 ){
                     // Make a bond
                     CA[a][0] = b, CA[b][0] = a;
                 }
                 // Or if a is monomer and b is in cluster
                 else if( CIN[a] != -1 && CIN[b] == -1 ){
                     // Make a bond
                     CA[a][0] = b, CA[b][0] = a;
                 }
             }
            // If the two particles are both from inside cluster, there will be no interaction.
            if(CIN[a] != -1 && CIN[b] != -1 && CA[a][0] == b ){
                CA[a][0] = -1, CA[b][0] = -1;
            }
        }
         
// If the pair is A-B or A'-B', we allow it only between "attached monomers"(CIN=-1) or particles in the same cluster (CIN=CIN)
         // For A-B or A'-B' interaction
         else if(MMIM[SV[a]][SV[b]]==1 && CIN[a] == CIN[b]){
            // If "a" and "b" are both attached monomers
            if((CIN[a] == -1 && CA[a][0] != -1 && CA[b][0] != -1) || CIN[a] != -1 ){
                
                if(CA[a][1] == -1 && CA[a][2] != b && CA[b][1] == -1 && CA[b][2] != a){
                 CA[a][1]=b, CA[b][1]=a;
                }

                else if(CA[a][1] == -1 && CA[a][2] != b && CA[b][2] == -1 && CA[b][1] != a ){
                 CA[a][1]=b, CA[b][2]=a;
                }

                else if(CA[a][2] == -1 && CA[a][1] != b && CA[b][1] == -1 && CA[b][2] != a){
                 CA[a][2]=b, CA[b][1]=a;
                }

                else if(CA[a][2] == -1 && CA[a][1] != b && CA[b][2] == -1 && CA[b][1] != a){
                 CA[a][2]=b, CA[b][2]=a;
                }
                
             }
         }
         
         
         
         
     }
    
    /* if connected pair is no longer in the connected distance or "complementaly bonds between two cluster"*/
    if(R[a][b]>BoL){
    
    if( CA[a][0] == b && CA[b][0] == a/*particles are connected*/){
                                                        CA[a][0] = -1, CA[b][0] = -1;
    }
    
    else if(CA[a][1]==b && CA[b][1]==a){
                                                        CA[a][1] = -1, CA[b][1] = -1;
        }
    
    else if(CA[a][1]==b && CA[b][2]==a){
                                                        CA[a][1] = -1, CA[b][2] = -1;
    }
    
    else if(CA[a][2]==b && CA[b][1]==a){
                                                        CA[a][2] = -1, CA[b][1] = -1;
    }
    
    else if(CA[a][2]==b && CA[b][2]==a){
                                                        CA[a][2] = -1, CA[b][2] = -1;
    }
    }
    
}


// Given Species Vector, Cluster Identification Number, Connection Array, calculate depth of potential between particles
/*   if "a" and "b" are connected or not => if "CA[a]=b (CA[b])=a" or not    */
/*   if particles are inside the same cluster or not => if "CIA[a]=CIA[b]" or not" */


void CalculateEdepth(){
    int i;
    // If an element of the interaction matrix is brank => Edep=0
    if(IM[SV[a]][SV[b]]==0){
        Edep = 0;
    }
    
    
    else{

        for(i=0;i<3;i++){
            // if "a" and "b" are connected,
            if(CA[a][i]==b){
                
                if(i==0 ){ // Bond between complementaly particles (We ignore Cluster-Cluster, Monomer-Monomer)
                    Edep = IE2;
                }
                else{  // Bond between different type
                    // Interaction inside cluster, irreversible bond
                    if( CIN[a] == CIN[b] && CIN[a] != -1){
                        Edep = IE3;
                    }
                    // Attached monomer - Attached monomer
                    else {
                        Edep = IE1;
                    }
                }
            }
            if(CA[a][i]==b) goto OUT;
        }
        Edep = 0;
    }


OUT:;}









// CalculateForce Calculate Forces acting between particles "a" and "b" //
void CalculateForce(){
    int i;
    // Edep == 0 <=> Hard Sphere
    if(Edep == 0){
        if(1. < R[a][b]){
            F[a][b] = 0;
        }
        // IE0=7.8125
        // Rho is steepness of the Morse potential
        else if (R[a][b]<=1 ){
            F[a][b] = -2.0*Rho*IE0*Rho*(R[a][b]-1.0) ;
        }
    }
    
    // Edep != 0 <=> Interact through potential with depth Edep
    else{
        // If the distance between particles are larger than Roff => No interaction
        if(Roff < R[a][b]){
            F[a][b] = 0;
        }
        // If the distance is shorter than cut off length => Attractive Morse potential
        else if (1< R[a][b]){
            F[a][b] = -2.0*Rho*Edep*exp(Rho*(1.0-Roff)) * ( -1 + exp(Rho*(1.0-Roff)) - exp(Rho*(1.0 - 2.0*R[a][b] +Roff)) + exp(Rho*(Roff-R[a][b])) ) ;
        }
        // If the distance is shorter than particle's diameter => Hard core, harmonic repulsion
        else if (R[a][b]<=1){
            F[a][b] = -2.0*Rho*( exp(Rho*(1.0-Roff))*(exp(Rho*(1.0-Roff)) -1)*Edep + IE0*Rho*(R[a][b]-1.0) );
        }
    }
    // Calculate X, Y components of force
    for(i=0;i<Dim;i++){
        f[a][b][i] =  F[a][b]*(D[a][b][i]/R[a][b]);
        f[b][a][i] =  -f[a][b][i];
    }
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// Function to sum up all pairwise forces
void SumForces(){
    int i,j,k,l,m,n;
    for(i=0;i<N;i++){
        for(j=0;j<Dim;j++){
            fpnew[i][j]=0.0;
        }
    }

    for(k=0;k<N;k++){
        for(l=0;l<NofNeighbors[k];l++){
            n=NeighborList[k][l];
            for(m=0;m<Dim;m++){
                fpnew[k][m]=fpnew[k][m]+f[k][n][m];
            }
        }
    }
}


//%%%%%%%%%%%%%%%%% Generate Normally Distributed Numbers with Box-Muller transform %%%%%%%%%%%%%%%%%//
// return numbers with random normal distribution //
double RandNormal(){
    double x=sqrt(-2.0*log(genrand_real3()))*sin(2.0*M_PI*genrand_real3());
    return x;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Renew(){
    int i,j;
    for(i=0;i<N;i++){
        for(j=0;j<Dim;j++){
            rold[i][j]   = rnew[i][j]  ;
            fpold[i][j]  = fpnew[i][j];
        }
    }
}








//%%%%%%%%%%%%%%%%%%%%%%%%  BondBreak function  %%%%%%%%%%%%%%%%%%%%%%%%%//
void BondBreak(){
    int i,j,l;
    // first you create "visited" array to mark particles which you've already visited
    // visited[] = 0 ( yet ), visited[] = 1 (visited)
    for(i=0;i<N;i++){
        visited[i]=0;
    }
    
    // look for attached monomer to start depth search
    int StartP;
    for(StartP=0;StartP<N;StartP++){
    // if the monomer is connected with some cluster and you haven't visited, choose this particle to be StartP
        if( CIN[StartP]==-1 && CA[StartP][0] != -1 && CIN[CA[StartP][0]] != -1 && visited[StartP]==0 ){
    // create NewCl[] array to store particles in a newly formed cluster
            for(j=0;j<N;j++){
                NewCl[j] = -1;
            }
    // number of particles in a new cluster -1
            NumNewCl = 0;
            visited[StartP]=1;
            NewCl[0]=StartP;
            Search( StartP );
            
            if(CountBond()>=Nc){
                
                for(j=0;j<NumNewCl+1;j++){ // there should be NumNewCl+1 particles within the chain
                    CIN[NewCl[j]] = NumCl; // assign "Cluster Identification Number" to the particles within the cluster
                }
                
                char filename[60];
                sprintf(filename,"NumClNc3N%dL%dMM%dCM%dRun%d.dat",N,Lnote,IE1note,IE2note,RunNnote);
                
                //printf("Number of bonds is %d\n",CountBond());
                //printf("Number of particles in the cluster is %d\n", NumNewCl + 1);
                NPIC[NumCl] = NumNewCl + 1;
                NumCl = NumCl +1;
                numcl    =   fopen(filename,"a");
                fprintf(numcl, "%d\t\t%d\t\t%d\t\t",t, NumCl,NumNewCl+1);
                
                //printf("Memeber of new cluster will be\t");
                //int i;
                //for(i=0;i<=NumNewCl;i++){
                //  printf("%d\t",NewCl[i]);
                //}
                int Parents[NumNewCl+1];
                for(l=0;l<NumNewCl+1;l++){
                    Parents[l] = CIN[CA[NewCl[l]][0]];
                }
                RemoveDuplicate( Parents);
            }
        }
    }
}















// Given "connection array", "cluster identification number", "Number of particles in the cluster"
// 1. List connected attaced monomers with NewCl
// 2. Count number of particles in the cluster
// 3. Count number of bonds within the cluster
// 4.


// Search( int FocusP ) Given a PIN of FocusP, look for another connected monomer which you haven't visited.
// NumCl[0] = StartP, NumCl[1] =
void Search( int FocusP ){
    int i;
    for(i=1;i<3;i++){
        /* if you find another "attached monomer" connected to it... */
        if( CA[FocusP][i] != -1 && visited[ CA[FocusP][i]] == 0 && CIN[ CA[FocusP][i] ] == -1 && CA[ CA[ FocusP][i] ][0] != -1 ){
            NumNewCl = NumNewCl + 1;
            NewCl[NumNewCl] = CA[FocusP][i];
            visited[ CA[FocusP][i] ] = 1;
            Search(CA[FocusP][i]);
        }
    }
}


/*--- Count number of bonds with in attached monomers ---*/
int CountBond(){
    int i,j,BonN;
    BonN = 0;
    
    /*-- Loop through all clusters exist ---*/
    for( i=0; i<=NumNewCl; i++ ){
        for(j=1;j<3;j++){
            /* counting number of bonds between "attached" monomers */
            if( NewCl[i] != -1 && CA[NewCl[i]][j] != -1 &&  CA [ CA[ NewCl[i] ][j] ][0] != -1 ){
                BonN = BonN +1;
            }
        }
    }
    return BonN/2;
}


/*--- Remove duplicate from Parents[] array ---*/
void RemoveDuplicate(int Parents[]){
    int *p;
    int i,j,k,size;
    
    size=NumNewCl+1;
    
    p=Parents;
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            if(i==j){
                continue;
            }
            else if(*(p+i)==*(p+j)){
                k=j;
                size--;
                while(k < size){
                    *(p+k)=*(p+k+1);
                    k++;
                }
                j=0;
            }
        }
    }
    //printf("\nThe array after removing duplicates is: ");
    //for(i=0;i < size;i++){
    //    printf(" %d",Parents[i]);
    //}
    //printf("\n");
    
    //printf("\nNumber of particles inside parents are: ");
    for(i=0;i < size;i++){
        //    printf(" %d",NPIC[Parents[i]]);
        fprintf(numcl, "%d",NPIC[ Parents[i] ]);
    }
    fprintf(numcl, "\n");
    fclose(numcl);
    
}













