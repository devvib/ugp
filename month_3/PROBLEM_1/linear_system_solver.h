// linear_system_solver.h

#ifndef LINEAR_SYSTEM_SOLVER_H
#define LINEAR_SYSTEM_SOLVER_H

#include<algorithm>
#include "Eigen/Dense"
#include <vector>
#include <iostream>
#include "antoine_solver.h"
using namespace std;

// Function to solve the system of linear equations AX = B
vector<double> solveLinearSystem(const vector<vector<double>>& A, const vector<double>& B) {
 
      

    // Convert input vectors to Eigen matrices
     size_t numCols = A[0].size();
    Eigen::MatrixXd eigenA(A.size(), numCols);
    Eigen::VectorXd eigenB(B.size());

    // Copy data from vector<vector<double>> to Eigen::MatrixXd
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            eigenA(i, j) = A[i][j];
        }
    }

    // Copy data from vector<double> to Eigen::VectorXd
    for (size_t i = 0; i < B.size(); ++i) {
        eigenB(i) = B[i];
    }


    // Solve for X using Eigen's linear solver
    Eigen::VectorXd eigenX = eigenA.colPivHouseholderQr().solve(eigenB);

    // Convert Eigen vector to std::vector
    std::vector<double> result(eigenX.data(), eigenX.data() + eigenX.size());

    return result;
}
vector<vector<int>> transposeMatrix(const vector<vector<int>>& matrix) {
    // Get the number of rows and columns
    int rows = matrix.size();
    int cols = matrix[0].size();

    // Transpose the matrix
    vector<vector<int>> transpose(cols, vector<int>(rows));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transpose[j][i] = matrix[i][j];
        }
    }

    return transpose;
}
 vector<vector<double>> matrix_solver(vector<vector<double>>S,vector<double> lf){
        vector<vector<double>>l;
for(int i=1;i<=2;i++){
      vector<vector<double>> A=calculate_Ai(S,i);
      // vector<vector<double>> A2=calculate_Ai(S,2);
      vector<double> B={{0},{0},{-lf[i]},{0}};

      // matrix solver
      vector<double> l1;
      try
      {
         l1 = solveLinearSystem(A, B);
      }
      catch (const std::invalid_argument &e)
      {
        cerr << "Error: " << e.what() << std::endl;
      }
      l.push_back(l1);
      }
      return l;
}

#endif // LINEAR_SYSTEM_SOLVER_H
