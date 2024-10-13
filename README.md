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
## ugp

Multicomponent Distillation Code Summary
1. Purpose:
The code aims to simulate multicomponent distillation using rigorous methods. It involves iterative processes to determine various parameters such as temperature, compositions, vapor and liquid flows, and enthalpy during the distillation process.

2. Main Process:
Initialization:

Initialize parameters and variables representing system properties such as pressure, number of components, stages, and feed mass flow.

- Equilibrium Calculations:

Iteratively calculate equilibrium conditions, including equilibrium constants (K) and distribution coefficients (S) for each component at each stage.
- Temperature Calculation:

Use Newton's method to iteratively determine the temperature at each stage until convergence is achieved.
Vapor and Liquid Flows:

Calculate vapor flow rates (vij) and vapor phase compositions (Yij) for each component at each stage.
- Enthalpy Calculations:

Determine enthalpy parameters (Hij and hij) for each component in both liquid and vapor phases.
- Mass and Energy Balances:

Implement mass and energy balances to ensure consistency in the distillation process.
Multicomponent Distillation Simulation
Overview
This project simulates a multicomponent distillation process using C++. It calculates temperature, vapor, and liquid flows, balancing energy and mass for a 2-component, 4-stage distillation process.

