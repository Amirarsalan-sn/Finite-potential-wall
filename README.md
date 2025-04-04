# Finite Potential Wall

In this project I've developed a program using python that simulates the behavior of a quantum particle in the Finite Potential Wall problem, a problem in which quantum tunneling is observed and we can see this phenomenon in practice via this simulation program.

It is easy to use the program. Hence, this README mainly focuses on the formalism of the problem and driving the equations from scratch (Schrödinger equation). To run a simulation, you just have to adjust the global variables at the top of the code to your desired value and run the main program, after that a video of the simulation will be saved on your computer.

## Table of Content

## Problem Specification
We first start by writing the time dependent and independent Schrödinger equations:

```math
i\hbar\frac{\partial\Psi(x,t)}{\partial t} = \left[\frac{-\hbar^2}{2m}\times \frac{d^2}{dx^2}  + V(x)\right]\Psi(x, t), \left(\Psi(x, t) = \psi(x)\phi(t)\right)
```
```math
\implies \left[\frac{-\hbar^2}{2m}\times \frac{d^2}{dx^2}  + V(x)\right]\psi(x) = E\psi(x), \left( \phi(t) = e^{\frac{-i}{\hbar}Et}\right)
```

Therefore, the problem is to solve the time independent Schrödinger equations for this potential:

```math
V(x) = \begin{cases}
V_0 & |x| \le a \\
0 & \text{else where}
\end{cases}
```
Which defines a finite (positive) potential barrier with finite width $$2a$$. In the next section, the solutions to the time independent Schrödinger equations will be explored for $$E < V_0$$ and $$E > V_0$$.

## Eigenfunctions of Hamiltonian



## Transmition and Reflection Coefficients for E < V0

## Transmition and Reflection Coefficients for E > V0

## Simple Demonstration

## Conclusion
