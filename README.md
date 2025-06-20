# Black-Scholes Equation Solver using Theta Scheme

This MATLAB project implements a numerical solution for the **Black-Scholes partial differential equation** using the **finite difference method** with the θ-scheme.

##  Problem Description

We solve the PDE:
∂V/∂t = (σ² / 2) ∂²V/∂x² + (r - σ²/2) ∂V/∂x - rV
with given initial and boundary conditions, where:
- `V(x,t)` is the option price,
- `σ` is the volatility,
- `r` is the risk-free interest rate,
- `x = log(S/E)` is the transformed space variable (`S` is the asset price, `E` is the strike price).

---
- Uses a **uniform grid** in log-space for spatial discretization.
- Implements **θ-scheme** (with θ = 0.5, i.e. Crank-Nicolson method).
- Compares the **numerical solution** to the **exact analytical solution** of the European call option using the Black-Scholes formula.
- Visualizes results in both `(x,t)` and `(S,t)` coordinates.
- Computes and visualizes the **error surface** between numerical and exact solutions.
