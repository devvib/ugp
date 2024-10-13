# Multicomponent Distillation Simulation

## Overview
This project simulates a multicomponent distillation process using C++. It calculates temperature, vapor, and liquid flows, balancing energy and mass for a 2-component, 4-stage distillation process.

## Project Structure
- **index.cpp**: Core simulation logic.
- **antoine.h**: Manages the temperature-pressure relationship using Antoine’s equation.
- **linear_equation_solver.h**: Solves linear systems essential to the simulation.
- **energy_solver.h**: Handles energy balance computations.

## Skills Demonstrated
- **Modular Design**: Clear separation of concerns across different headers improves code maintainability.
- **Numerical Methods**: Implements solvers for energy and linear equations, showcasing problem-solving abilities.
- **C++ Best Practices**: Use of standard libraries like **Eigen** for matrix solving.

## Compilation and Usage

To compile and run the simulation:

```bash
g++ index.cpp -o index
./index

```
## UGP: Multicomponent Distillation Code Summary

### 1. Purpose:
The code simulates multicomponent distillation using rigorous methods. It involves iterative calculations to determine parameters like temperature, compositions, vapor and liquid flows, and enthalpy.

### 2. Main Process:

- **Initialization**: Initialize system properties such as pressure, number of components, stages, and feed mass flow.
- **Equilibrium Calculations**: Iteratively calculate equilibrium constants (K) and distribution coefficients (S) at each stage.
- **Temperature Calculation**: Use Newton’s method to determine stage-wise temperature.
- **Vapor and Liquid Flows**: Calculate vapor flow rates (vij) and vapor compositions (Yij).
- **Enthalpy Calculations**: Compute enthalpy for each component in both phases.
- **Mass and Energy Balances**: Ensure consistency with mass and energy balances.

---

## Multicomponent Distillation Simulation

### Overview
This project simulates a multicomponent distillation process using C++. It calculates temperature, vapor, and liquid flows, balancing energy and mass for a 2-component, 4-stage distillation process.


