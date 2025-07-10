# Fault-Tolerant Fractional Order Sliding Mode Control for Soft Continuum Manipulators

## Overview

This repository contains the implementation of a **Fault-Tolerant Fractional Order Sliding Mode Controller (FTFOSMC)** for soft continuum manipulators, utilizing the **Cosserat rod model** to describe their dynamic behavior. The method is aimed at achieving precise boundary control for soft manipulators under large deformations, external disturbances, and time-varying faults.

The controller leverages fractional calculus to model the memory-dependent behaviors of the manipulator, enhancing robustness and accuracy when compared to traditional controllers. Numerical simulations are provided to demonstrate the effectiveness of the proposed method under various fault conditions and model uncertainties.

## Key Features

- **Cosserat Rod Model**: Describes the large deformation behavior of the soft continuum manipulator.
- **Fractional Order Sliding Mode Controller**: Integrates the fractional order derivative for improved robustness.
- **Fault-Tolerant Control**: Designed to handle time-varying faults in the system.
- **Numerical Simulations**: Simulations are provided to demonstrate the controller's performance under different fault scenarios.
- **Data Visualization**: Includes plots of displacement, control forces, and error responses under different disturbance and fault conditions.

## Simulation Results

The numerical simulations demonstrate the control performance of the FTFOSMC under different scenarios:

1. **No disturbance**: Shows the baseline response of the manipulator.
2. **Boundary disturbance**: Response under external boundary forces.
3. **Model and boundary disturbance**: Response when both model errors and boundary disturbances are present.
4. **Time-varying fault**: Demonstrates the fault-tolerant capacity of the controller.
5. **Fractional order variation**: Evaluates the effect of different fractional order parameters on control performance.

## Repository Structure

```markdown
/Result
├── CosseratFTFOSMC\_modelerror.mat  # Model error simulation results
├── CosseratFTFOSMC.mat            # Default simulation results
├── CosseratFTFOSMC\_F001.mat       # Response with specific fault condition (F001)
├── CosseratFTFOSMC\_F005.mat       # Response with another fault condition (F005)
CantileverRodFTFOSMC.m               # Main MATLAB script implementing the controller
README.md                           # Documentation for the repository
```

### Key MATLAB Functions

- **CantileverRodFTFOSMC.m**: This is the main script where the control system is implemented. It uses the Cosserat rod model and applies the fractional order sliding mode control to the soft continuum manipulator.

- **Simulation Files (`.mat`)**: These files contain the results of the simulations under different conditions (e.g., with model errors, boundary disturbances, or fault conditions). These data files can be loaded into MATLAB for further analysis or visualization.

### Data Description

- **CosseratFTFOSMC_modelerror.mat**: Contains the simulation results under model errors.
- **CosseratFTFOSMC_F001.mat** and **CosseratFTFOSMC_F005.mat**: Contain the simulation results under specific fault conditions (`F001`, `F005`).

## Installation

To run the simulations and use the controller, make sure that you have MATLAB installed with the following toolboxes:

- **MATLAB R2021 or later** 
- **Control Systems Toolbox**
- **Optimization Toolbox**

## Usage

1. Open `CantileverRodFTFOSMC.m` in MATLAB.
2. Modify the parameters as needed for your specific use case.
3. Run the script to simulate the behavior of the soft continuum manipulator under different fault conditions and disturbances.
4. Use the saved `.mat` files to load and visualize the results.

## Citation

If you use this code in your work, please cite the following paper:

**Wang, Shuopeng, et al.**
Robust Fault-Tolerant Fractional Order Sliding Mode Boundary Control for a Soft Continuum Manipulator.

## Acknowledgments

This work was supported by the National Natural Science Foundation of China under grant number 62073063, 62461160260, and the Fundamental Research Funds for the Central Universities under grant number N2403020.

## License

This code is licensed under the MIT License. 
