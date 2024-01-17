#include<bits/stdc++.h>
#include "antoine_solver.h"
using namespace std;
double P_total=1;
int no_of_comp=2;
int no_of_stages=4;
int feed_mass_flow=1;
int main(){

   int ct=1;
   float A1,B1,C1,T;
   while(ct--){
    int ctt=5;
     while(ctt--)
     {
        //calculation of K (y=K*x);
          double k1=k_calculator(T,A1,B1,C1,P_total);
        //calcution of S (S[1][2]=k[1][2]*V[2]/L[2]);
        //write function for S in antoine_solver or make separate header file

                  
        
        //solving equation
         // Find the temperature using Newton's method
       double T1 = findTemperature(initialGuess, x1, A1, B1, C1, x2, A2, B2, C2, P_total, tolerance);
     }
   }
    return 0;
}