# Modeling and Simulation of Dynamic Systems

This repository contains coursework for "Modeling and Simulation of Dynamic Systems" (ECE AUTH 2024-2025), 
featuring two lab assignments and a final project. 
All implementations are developed in MATLAB and explore various modeling and simulation techniques for dynamic systems analysis.

## Repository Structure

The work is organized into three main assingments, each focusing on different aspects of dynamic systems modeling:

### Lab 01: Pendulum Dynamics and Parameter Estimation

Explores fundamental concepts in dynamic system modeling using a pendulum system:

- **ODE-based Simulation**: Implements numerical solutions for pendulum equations of motion using the `ode45`
- **Parameter Estimation**: Applies the Least Squares Method to identify system parameters from observed data
- **Noise Analysis**: Investigates the effects of measurement noise on both state estimation and parameter identification accuracy

### Lab 02: Online Parameter Estimation

Focuses on real-time (online) parameter estimation techniques:

- **Online Parameter Estimation Methods**: Implements and compares multiple approaches including:
  - Gradient Method for iterative parameter updates
  - Lyapunov method in parallel configuration
  - Lyapunov method in series-parallel configuration
- **Noise Impact Assessment**: Analyzes how different estimation methods perform under various noise conditions

### Final Project: Constrained Estimation and Black Box System

- **Constrained Parameters Estimation**: Utilizes the series-parallel configuration Lyapunov method with projection techniques to ensure solutions remain within feasible parameters' bounds
- **Polar Error Integration**: Implements a gradient-based method to estimate parameters in a system affected by polar error ω(t), using the continuous version of the σ-modification technique
- **Black Box System Identification**: Addresses model structure selection and online parameter estimation for an unknown system

Each section of the repository is self-contained, with source code, utility functions, and output plots. 
