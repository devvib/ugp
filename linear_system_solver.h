// linear_system_solver.h

#ifndef LINEAR_SYSTEM_SOLVER_H
#define LINEAR_SYSTEM_SOLVER_H

#include <Eigen/Dense>
#include <vector>

// Function to solve the system of linear equations AX = B
std::vector<double> solveLinearSystem(const std::vector<std::vector<double>>& A,
                                      const std::vector<double>& B) {
    // Convert input vectors to Eigen matrices
    Eigen::Map<const Eigen::MatrixXd> eigenA(A[0].data(), A.size(), A[0].size());
    Eigen::Map<const Eigen::VectorXd> eigenB(B.data(), B.size());

    // Solve for X using Eigen's linear solver
    Eigen::VectorXd eigenX = eigenA.colPivHouseholderQr().solve(eigenB);

    // Convert Eigen vector to std::vector
    std::vector<double> result(eigenX.data(), eigenX.data() + eigenX.size());

    return result;
}

#endif // LINEAR_SYSTEM_SOLVER_H
