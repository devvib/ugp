#include<bits/stdc++.h>
#include "antoine_solver.h"
using namespace std;
double P_total=1;
int no_of_comp=2;
int no_of_stages=4;
int feed_mass_flow=1;
int main(){

   int ct=1;
   int K[2][4],V[5],L[5];//make it vector
   float A1,B1,C1,T;
   while(ct--){
    int ctt=5;
     while(ctt--)
     {
        //calculation of K (y=K*x);
          double k1=K_calculator(T,A1,B1,C1,P_total);
        //calcution of S (S[1][2]=k[1][2]*V[2]/L[2]);
        double S12=S_calculator(K[1][2],V[2],L[2]);
        //make matrix A,X and B for A*X=C
        //matrix solver

                  
        
        //solving equation
         // Find the temperature using Newton's method
       double T1 = findTemperature(initialGuess, x1, A1, B1, C1, x2, A2, B2, C2, P_total, tolerance);
     }
   }
    return 0;
}