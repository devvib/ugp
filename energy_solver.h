// energy_solver.h

#ifndef ENERGY_SOLVER_H
#define ENERGY_SOLVER_H

#include <cmath>
#include <iostream>
#include<bits/stdc++.h>
using namespace std;

// Function to calculate enthalpy for a component
double calculateEnthalpy(double A, double B, double C, double T) {
    return A + B * T + C * T * T;
}

// Function to calcualte the vapour phase composition.
inline vector<vector<double>> calculate_vij(
    const vector<vector<double>>& Sij,
    const vector<vector<double>>& lij
) {
    vector<vector<double>> vij;

    // Assuming Sij and lij have the same dimensions
    size_t numComponents = Sij.size();
    size_t numStages = Sij[0].size();  // Assuming all components have the same number of stages

    // Resize vij to match the dimensions of Sij and lij
    vij.resize(numComponents, vector<double>(numStages, 0.0));

    // Calculate vij
    for (size_t i = 1; i < numComponents; ++i) {
        for (size_t j = 2; j < numStages; ++j) {
            vij[i][j] = Sij[i][j] * lij[i][j];
        }
    }

    return vij;
}

inline vector<vector<double>> calculate_Yij(
    const vector<vector<double>>& vij
) {
    vector<vector<double>> Yij;

    // Assuming vij have the correct dimension
    size_t numComponents = vij.size();
    size_t numStages = vij[0].size();  // Assuming all components have the same number of stages

    // Resize Yij to match the dimensions of vij
    Yij.resize(numComponents, vector<double>(numStages, 0.0));

    // Calculate Yij
    for(size_t j = 2; j < numStages; ++j){
        double sumofall = 0.0;
        for(size_t i = 1; i< numComponents; ++i){
            sumofall += vij[i][j];
        }
        for(size_t i = 1; i< numComponents; ++i) {
            Yij[i][j] = vij[i][j] / sumofall;
        }

    }
    return Yij;
}


#endif  // ENERGY_SOLVER_H
