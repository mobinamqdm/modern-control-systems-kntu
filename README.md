# Modern Control Project — Quadruple Tank System

This repository contains the main project for the **Modern Control** course at **K. N. Toosi University of Technology**, instructed by **Dr. Bijan Moaveni**.  
The project focuses on the **analysis, linearization, and control design** of a nonlinear **Quadruple Tank Process (QTP)** system using **MATLAB** and **Simulink**.

---

## 📘 Project Overview
The project models, analyzes, and controls a **nonlinear multi-input multi-output (MIMO)** system — the **Quadruple Tank Process**.  
It aims to:
- Derive the nonlinear model based on physical laws.  
- Linearize the system around equilibrium points.  
- Evaluate stability, controllability, and observability.  
- Design and simulate state feedback controllers, observers, and compensators.  
- Compare the performance of **linear vs nonlinear** models and control methods.

This implementation demonstrates the transition from **theoretical modern control concepts** to **practical simulation-based validation** using MATLAB and Simulink.

---

## 🎯 Objectives
- Develop nonlinear and linearized models of the Quadruple Tank Process.  
- Analyze eigenvalues, stability, controllability, and observability.  
- Simulate open-loop and closed-loop responses under different inputs (step, ramp, impulse).  
- Design:
  - State feedback controllers  
  - Static pre-compensator  
  - Integral controller  
  - Luenberger observer and observer-based controller  
- Compare results between **MATLAB** and **Simulink**, validating accuracy across both platforms.

---

## 🧩 Repository Structure

modern_project/
├─ MATLAB-code/
│ ├─ control-modern-project.m # Main script: nonlinear, linear, and control design
│ ├─ quadruple-tank-script.m # Nonlinear & linearized model simulation
│
├─ simulink/
│ ├─ quadruple-tank-model.slx # Nonlinear system model
│ ├─ sym-project.slx # Linearized or symbolic-based control model
│
├─ report/
│ └─ modern-control-project-report.pdf
└─paper/
    └─
