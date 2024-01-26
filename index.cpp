#include <bits/stdc++.h>
#include "antoine_solver.h"
#include "linear_system_solver.h"
#include "energy_solver.h"
using namespace std;
double P_total = 1;
int no_of_comp = 2;
int no_of_stages = 4;
int feed_mass_flow = 1;
double D=0.5;
int main()
{
  vector<double>T(no_of_stages+1),A(no_of_comp+1),
  B(no_of_comp+1),C(no_of_comp+1),V(no_of_stages+1),
  L(no_of_stages+1),lf(no_of_comp+1);
  vector<vector<double>>Sij,lij,xij;
  while (true)
  {
    //what to initialise here??



   //loop2 needs T,V,L
    while (true)
    {
      // calculation of K (y=K*x);
      vector<vector<double>> K = calculate_K(A,B,C,T, P_total,no_of_stages);

      // calcution of S (S[1][2]=k[1][2]*V[2]/L[2]);
      vector<vector<double>> S = calculate_S(K, V, L,D);

      // make matrix A,X and B for A*X=C
      vector<vector<double>>l=matrix_solver(S,lf);

      // we got X all l
      // now we calculate all xij;
      vector<vector<double>>x=calculate_x(l);

      //  Find the temperature using Newton's method
      vector<double> Tnew= temp_solver(T,A,B,C,x);
      //breaking condition for T
      if(accurate(Tnew,T)){Sij=S,lij=l; T=Tnew;xij=x;break;}
    }
    // Initializing and calculating vij.
    vector<vector<double> > vij = calculate_vij(Sij,lij);
    // Initializing and calculating Yij.
    vector<vector<double> > Yij = calculate_Yij(vij);

    // Assuming that we've taken Enthalpy parameters for every component in 2D vectors named H_param and h_param.

    // Initializing and calculating Hij.

    //make H_param in A,B,C form if possible??
    vector<vector<double> > Hij = calculate_Hij(H_param,T);

    // Initializing and calculating hij.
    vector<vector<double> > hij = calculate_hij(h_param,T);

    // Initializing and calculating Hi.
    vector<double> Hi = calculate_Hi(Hij,Yij);

    // Initializing and calculating hi.
    vector<double> hi = calculate_hi(hij,xij);

    // Calculating hfi
    vector<double> hfi = calculatehfi(h_param, Feed_T);

    // Assuming distillate as 'distillate'
    vector<double> Vnew = calculate_Vnew(Hi, hi, Li, Vi, distillate, Feed, hfi);


  }
  return 0;
}