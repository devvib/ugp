// antoine_solver.h

#ifndef ANTOINE_SOLVER_H
#define ANTOINE_SOLVER_H

#include <cmath>
#include <iostream>
using namespace std;
// Function to calculate saturation vapor pressure using Antoine equation
double calculateSaturationPressure(double T, double A, double B, double C) {
    return pow(10, (A - B / (T + C)));
}
// Function to calculate the total vapor pressure equation
double totalPressureEquation(double T, double x1, double A1, double B1, double C1, double x2, double A2, double B2, double C2) {
    double P1=calculateSaturationPressure(T,A1,B1,C1);
    double P2=calculateSaturationPressure(T,A2,B2,C2);
    return x1 * P1 + x2 * P2;
}

// Derivative of the total vapor pressure equation with respect to temperature
double derivativeOfTotalPressureEquation(double T, double x1, double A1, double B1, double C1, double x2, double A2, double B2, double C2) {
    return -x1 * B1 / ((T + C1) * log(10) * pow(10, (A1 - B1 / (T + C1)))) - x2 * B2 / ((T + C2) * log(10) * pow(10, (A2 - B2 / (T + C2))));
}

// Newton's method to solve for temperature
double findTemperature(double initialGuess, double x1, double A1, double B1, double C1, double x2, double A2, double B2, double C2, double P_total, double tolerance) {
    double T_old = initialGuess;

    while (true) {
        // Calculate f(T) and f'(T)
        double f_T = totalPressureEquation(T_old, x1, A1, B1, C1, x2, A2, B2, C2) - P_total;
        double f_prime_T = derivativeOfTotalPressureEquation(T_old, x1, A1, B1, C1, x2, A2, B2, C2);

        // Update temperature
        double T_new = T_old - f_T / f_prime_T;

        // Check for convergence
        if (abs(T_new - T_old) < tolerance) {
            break;
        }

        // Update T_old for the next iteration
        T_old = T_new;
    }

    return T_old;
}

#endif // ANTOINE_SOLVER_H
//how to use in cpp 
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
