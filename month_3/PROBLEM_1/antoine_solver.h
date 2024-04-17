//antoine_solver.h

#ifndef ANTOINE_SOLVER_H
#define ANTOINE_SOLVER_H

#include <cmath>
#include <iostream>
#include<vector>
#include<bits/stdc++.h>
using namespace std;
// Function to calculate saturation vapor pressure using Antoine equation
double calculateSaturationPressure(double T, double A, double B, double C, double val)
{
    return (exp(A - B / (T + C)))*val;
}
// Function to calculate the total vapor pressure equation
double totalPressureEquation(double T, double x1, double A1, double B1, double C1, double val1, double x2, double A2, double B2, double C2, double val2)
{
    double P1 = calculateSaturationPressure(T, A1, B1, C1, val1);
    double P2 = calculateSaturationPressure(T, A2, B2, C2, val2);
    return x1 * P1 + x2 * P2;
}

// Derivative of the total vapor pressure equation with respect to temperature
double derivativeOfTotalPressureEquation(double T, double x1, double A1, double B1, double C1, double val1, double x2, double A2, double B2, double C2, double val2)
{
   return (x1*B1*exp(A1 - B1 / (T + C1))*val1 ) / ((T+C1)*(T+C1)) + (x2*B2*exp(A2 - B2 / (T + C2))*val2 ) / ((T+C2)*(T+C2));
}

// Newton's method to solve for temperature
double findTemperature(double initialGuess, double x1, double A1, double B1, double C1, double x2, double A2, double B2, double C2, double P_total, double val1, double val2)
{
    double T_old = initialGuess;

    while (true)
    {
        // cout<<"yha pe"<<endl;
        // Calculate f(T) and f'(T)
        // val1=1;
        // val2=1;
        double f_T = totalPressureEquation(T_old, x1, A1, B1, C1, val1, x2, A2, B2, C2,val2) - P_total;
        double f_prime_T = derivativeOfTotalPressureEquation(T_old, x1, A1, B1, C1, val1, x2, A2, B2, C2, val2);

        // Update temperature
        double T_new = T_old - f_T / f_prime_T;

        if (abs(T_new - T_old) < 0.1)
        {
            break;
        }

        // Update T_old for the next iteration
        T_old = T_new;
    }

    return T_old;
}
double K_calculator(double T, double A, double B, double C,double val, double P_total)
{
    return (calculateSaturationPressure(T, A, B, C,val) / P_total);
}
double S_calculator(double K, double V, double L)
{
    return K * V / L;
}

// function that return a 2d vector of K
vector<vector<double>> calculate_K(vector<double> A, vector<double> B, vector<double> C, vector<vector<double>> gamma, vector<double> T, double P_total, int no_of_stages)
{
    int no_of_comp = A.size() - 1;
    vector<vector<double>> Kij(no_of_comp + 1, vector<double>(no_of_stages + 1, 0));

    for (int i = 1; i <= no_of_comp; i++)

        for (int j = 1; j <= no_of_stages; j++)

            Kij[i][j] = K_calculator(T[j], A[i], B[i], C[i],gamma[i][j], P_total);

    return Kij;
}
// function that return a 2d vector of S
vector<vector<double>> calculate_S(vector<vector<double>> K, vector<double> V, vector<double> L, float D)
{
    int no_of_comp = K.size() - 1;
    int no_of_stages = K[0].size() - 1;
    vector<vector<double>> S(no_of_comp + 1, vector<double>(no_of_stages + 1, 0));

    for (int i = 1; i < no_of_comp + 1; i++)

        for (int j = 1; j < no_of_stages + 1; j++)
        {

            if (j == 1)
                S[i][j] = D / L[i];
            else
                S[i][j] = (K[i][j] * V[j]) / L[j];
        }

    return S;
}
// function that return a 2d vector of Xij
vector<vector<double>> calculate_x(vector<vector<double>> l)
{
    int no_of_comp = l.size() ;
    int no_of_stages = l[0].size();
    vector<vector<double>> x(no_of_comp + 1, vector<double>(no_of_stages + 1, 0));

    for (int i = 1; i <= no_of_stages; i++)
    {
        double denominator = 0;

        for (int j = 1; j <= no_of_comp; j++)
            denominator += l[j-1][i-1];

        for (int j = 1; j <= no_of_comp; j++)
            x[j][i] = l[j-1][i-1] / denominator;
    }
    return x;
}
vector<vector<double>>calculate_Ai(vector<vector<double>>S,int comp_no){
    int i=comp_no;
    vector<vector<double>>A={
        {-(1+S[i][1]),S[i][2],0,0},
        {1,-(1+S[i][2]),S[i][3],0},
        {0,1,-(1+S[i][3]),S[i][4]},
        {0,0,1,-(1+S[i][4])},

    };
    return A;

}
vector<double> temp_solver(vector<double> T,vector<double>A,vector<double>B,vector<double>C,vector<vector<double>> gamma, vector<vector<double>>x,float P_total){
    int no_of_stages=T.size()-1;

    for(int i=1;i<=no_of_stages;i++)
     T[i]= findTemperature(T[i], x[1][i], A[1], B[1], C[1], x[2][i], A[2], B[2], C[2], P_total,gamma[1][i],gamma[2][i]);
     return T;
}

bool accurate(vector<double>Tnew,vector<double>T){
    for(int i=1;i<T.size();i++)if(abs(Tnew[i]-T[i])>0.01)return 0;
    return 1;
}


vector<vector<double>> calculate_gamma(vector<vector<double>> x, vector<double> lambda)
{
  double no_of_comp = x.size()-1 ;
  double no_of_stages = x[0].size()-1;
  vector<vector<double>> gamma(no_of_comp + 1, vector<double>(no_of_stages + 1));

  for (double i = 1; i <= no_of_stages; i++)
  {
    double val1 = ((lambda[1]/(x[1][i] + (x[2][i]*lambda[1]))) - ( lambda[2]/(x[2][i]+(x[1][i]*lambda[2]))));

    double val2 = val1;
    val2 *= x[2][i];

    double val3 = val1;
    val3 *= x[1][i];
    val3 *= -1;

    double val4 = (x[1][i]+(x[2][i]*lambda[1]));
    val4 = log(val4);
    val4 *= -1;

    double val5 = (x[2][i]+(x[1][i]*lambda[2]));
    val5 = log(val5);
    val5 *= -1;

    val4 += val2;
    val5 += val3;

    gamma[1][i] = exp(val4);
    gamma[2][i] = exp(val5);
    
  }
  return gamma;
  
}





#endif // ANTOINE_SOLVER_H
// how to use in cpp
/*
// main.cpp

#include "antoine_solver.h"

int main() {
    // Example parameters (replace with your specific values)
    double x1 = 0.5, A1 = 10.0, B1 = 1000.0, C1 = 20.0;
    double x2 = 0.5, A2 = 12.0, B2 = 1200.0, C2 = 25.0;
    double P_total = 800.0; // Replace with your total pressure
    double initialGuess = 300.0; // Replace with your initial guess
    double tolerance = 1e-6;

    // Find the temperature using Newton's method
    double temperature = findTemperature(initialGuess, x1, A1, B1, C1, x2, A2, B2, C2, P_total, tolerance);

    std::cout << "Estimated temperature: " << temperature << " degrees Celsius" << std::endl;

    return 0;
}

*/