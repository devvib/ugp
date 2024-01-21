#include <bits/stdc++.h>
#include "antoine_solver.h"
#include "linear_system_solver.h"
#include "energy_solver.h"
using namespace std;
double P_total = 1;
int no_of_comp = 2;
int no_of_stages = 4;
int feed_mass_flow = 1;
int main()
{

  int ct = 1;
  int K[2][4], V[5], L[5]; // make it vector
  float A1, B1, C1, T;
  while (ct--)
  {
    int ctt = 5;
    while (ctt--)
    {
      // calculation of K (y=K*x);
      double k1 = K_calculator(T, A1, B1, C1, P_total);
      // calcution of S (S[1][2]=k[1][2]*V[2]/L[2]);
      double S12 = S_calculator(K[1][2], V[2], L[2]);
      // make matrix A,X and B for A*X=C
      vector<vector<double>> A;
      vector<double> B;

      // matrix solver

      try
      {
        vector<double> X = solveLinearSystem(A, B);
      }
      catch (const std::invalid_argument &e)
      {
        cerr << "Error: " << e.what() << std::endl;
      }
      // we got X all l

      // now we calculate all xij;

      // write function for calculation of xij;

      // solving equation
      //  Find the temperature using Newton's method
      double T1 = findTemperature(initialGuess, x1, A1, B1, C1, x2, A2, B2, C2, P_total, tolerance);
    }
    // Initializing and calculating vij.
    vector<vector<double> > vij = calculate_vij(Sij,lij);
    // Initializing and calculating Yij.
    vector<vector<double> > Yij = calculate_Yij(vij);

    // Assuming that we've taken Enthalpy parameters for every component in 2D vectors named H_param and h_param.

    // Initializing and calculating Hij.
    vector<vector<double> > Hij = calculate_Hij(H_param,T);

    // Initializing and calculating hij.
    vector<vector<double> > hij = calculate_hij(h_param,T);

    // Initializing and calculating Hi.
    vector<double> Hi = calculate_Hi(Hij,Yij);

    // Initializing and calculating hi.
    vector<double> hi = calculate_hi(hij,xij);

  }
  return 0;
}