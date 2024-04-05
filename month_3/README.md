# ugp
This project is about "Rigorous methods of multicomponent distillation".
🥲Problem : We have multistage distillation column, we will pass multicompnent feed and by passing it , we want to separate out these components.
For this problem we have 2 components and 4 stages.


Multicomponent Distillation Code Summary
1. Purpose:
The code aims to simulate multicomponent distillation using rigorous methods. It involves iterative processes to determine various parameters such as temperature, compositions, vapor and liquid flows, and enthalpy during the distillation process.

2. Main Process:
Initialization:

Initialize parameters and variables representing system properties such as pressure, number of components, stages, and feed mass flow.

Equilibrium Calculations:

Iteratively calculate equilibrium conditions, including equilibrium constants (K) and distribution coefficients (S) for each component at each stage.
Temperature Calculation:

Use Newton's method to iteratively determine the temperature at each stage until convergence is achieved.
Vapor and Liquid Flows:

Calculate vapor flow rates (vij) and vapor phase compositions (Yij) for each component at each stage.
Enthalpy Calculations:

Determine enthalpy parameters (Hij and hij) for each component in both liquid and vapor phases.
Mass and Energy Balances:

Implement mass and energy balances to ensure consistency in the distillation process.
