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
Now, to run a simulation, first, we need to find the eigenfunctions of hamiltonian and expand the initial state as a superposition of eigenstates. After that, the evolution of the system is given by multiplying the time component $$\phi(t)$$ to each eigenstate:

```math
\Psi(x, 0) = \int c(k) \psi_k (x) dk
```
```math
c(k') = \int_{-\infty}^{+\infty} \psi_{k'}^{*} (x) \Psi(x, 0) dx
```
```math
\Psi(x, t) = \int c(k) \psi_k (x) e^{\frac{-i}{\hbar}E_k t} dk
```
In this problem, the energy spectrum is divided into two parts $E < V_0$ and $E > V_0$, and the position has three sections ( $x < -a, -a \le x \le a, x > a$ ) and the time independent Schrödinger equation will be solved in each section and each energy spectrum seperately.

### E < V0
In this spectrum, Schrödinger equation is solved like below:

```math
\begin{cases}
\left[\frac{-\hbar^2}{2m}\times \frac{d^2}{dx^2}\right]\psi(x) = E\psi(x) & x < -a \\
\left[\frac{-\hbar^2}{2m}\times \frac{d^2}{dx^2} + V_0 \right]\psi(x) = E\psi(x) & -a \le x \le a \\
\left[\frac{-\hbar^2}{2m}\times \frac{d^2}{dx^2}\right]\psi(x) = E\psi(x) & x > -a \\
\end{cases}
```
```math
\implies
\begin{cases}
\frac{d^2}{dx^2}\psi(x) = \frac{-2m}{\hbar^2}E\psi(x) & x < -a \\
\frac{d^2}{dx^2}\psi(x) = \frac{-2m}{\hbar^2}\left(E - V_0\right)\psi(x) & -a \le x \le a \\
\frac{d^2}{dx^2}\psi(x) = \frac{-2m}{\hbar^2}E\psi(x) & x > -a \\
\end{cases}
```
The term $E - V_0$ is negative in this spectrum, we define $\kappa = \sqrt{\frac{2m}{\hbar^2}E}$ and $k = \sqrt{\frac{-2m}{\hbar^2}\left(E - V_0\right)}$, therefore, eigenstates of hamiltonian for $E < V_0$ would be:
```math
\psi_{\kappa}(x) = \begin{cases}
Ae^{i\kappa x} + Be^{-i\kappa x} & x < -a \\
Ce^{kx} + De^{-kx} & -a \le x \le a \\
Fe^{i\kappa x} + Ge^{-i\kappa x} & x > -a \\
\end{cases}
```
### E > V0
Morover, in this spectrum, Schrödinger equation is solved like below:

```math
\begin{cases}
\left[\frac{-\hbar^2}{2m}\times \frac{d^2}{dx^2}\right]\psi(x) = E\psi(x) & x < -a \\
\left[\frac{-\hbar^2}{2m}\times \frac{d^2}{dx^2} + V_0 \right]\psi(x) = E\psi(x) & -a \le x \le a \\
\left[\frac{-\hbar^2}{2m}\times \frac{d^2}{dx^2}\right]\psi(x) = E\psi(x) & x > -a \\
\end{cases}
```
```math
\implies
\begin{cases}
\frac{d^2}{dx^2}\psi(x) = \frac{-2m}{\hbar^2}E\psi(x) & x < -a \\
\frac{d^2}{dx^2}\psi(x) = \frac{-2m}{\hbar^2}\left(E - V_0\right)\psi(x) & -a \le x \le a \\
\frac{d^2}{dx^2}\psi(x) = \frac{-2m}{\hbar^2}E\psi(x) & x > -a \\
\end{cases}
```
The term $E - V_0$ is positive in this spectrum, we define $\kappa = \sqrt{\frac{2m}{\hbar^2}E}$ and $k = \sqrt{\frac{2m}{\hbar^2}\left(E - V_0\right)}$, therefore, eigenstates of hamiltonian for $E < V_0$ would be:
```math
\psi_{\kappa}(x) = \begin{cases}
Ae^{i\kappa x} + Be^{-i\kappa x} & x < -a \\
C\cos(kx) + D\sin(kx) & -a \le x \le a \\
Fe^{i\kappa x} + Ge^{-i\kappa x} & x > -a \\
\end{cases}
```

In this section, we found the eigenfunctions of hamiltonian which form an complete orthogonal basis for our hilbert space. In the next section, the transmission and reflection coefficients are going to be calculated by applying the boundary conditions to $\psi_{\kappa} (x)$.

## Transmission and Reflection Coefficients for E < V0
Before diving deep into calculations, it is worth noting that singularities and discontinuities in $\psi_{\kappa}$ is always two steps better than singularities and discontinuities in $V(x)$ and that's because of the time independent Schrödinger equation:

```math
\frac{d^2}{dx^2}\psi(x) = \frac{-2m}{\hbar^2}\left(E - V(x)\right)\psi(x)
```
As it is seen from the equation, in order to get $\psi_{\kappa}$, one should integrate $\frac{-2m}{\hbar^2}\left(E - V(x)\right)\psi(x)$ two times. Hence, if the potential has simple discontinuities (like the finite potential wall), $\psi_{\kappa}$ and $\frac{d}{dx}\psi_{\kappa}$ should be continueous, and that's what gives us the *boundary conditions*. With these conditions we are able to calculate the transmission and reflection coefficients.

Another important thing to note is that we normaly assume the eigenfunctions to scatter from left to right (or right to left) to make the calculation of coefficients easier and more interpretable. Without this assumption, it would be to hard too find any relationship between the coefficients of $\psi$ just by applying the boundary conditions. Thus, the value of $G$ is assumed to be zero throughout the calculations and simulations.

The eigenstates and boundary conditions for $E < V_0$ are:

```math
\psi_{\kappa}(x) = \begin{cases}
Ae^{i\kappa x} + Be^{-i\kappa x} & x < -a \\
Ce^{kx} + De^{-kx} & -a \le x \le a \\
Fe^{i\kappa x} & x > -a \\
\end{cases}
```
```math
\lim_{\epsilon \to 0}
\begin{cases}
\psi_{\kappa}(-a -\epsilon) = \psi_{\kappa}(-a +\epsilon) & (\text{I}) \\
\psi_{\kappa}(+a -\epsilon) = \psi_{\kappa}(+a +\epsilon) & (\text{II})\\
\frac{d}{dx}\psi_{\kappa}(-a -\epsilon) = \frac{d}{dx}\psi_{\kappa}(-a +\epsilon) & (\text{III}) \\
\frac{d}{dx}\psi_{\kappa}(+a -\epsilon) = \frac{d}{dx}\psi_{\kappa}(+a +\epsilon) & (\text{IV}) \\
\end{cases} 
```
```math
\implies \begin{cases}
Ae^{-i\kappa a} + Be^{i\kappa a} = Ce^{-ka} + De^{ka} & (\text{I}) \\
Fe^{i\kappa a} = Ce^{ka} + De^{-ka} & (\text{II}) \\ 
i\kappa\left( Ae^{-i\kappa a} - Be^{i\kappa a}\right) = k\left(Ce^{-ka} - De^{ka}\right) & (\text{III}) \\
i\kappa\left(Fe^{i\kappa a}\right) = k\left(Ce^{ka} - De^{-ka}\right) & (\text{IV}) \\
\end{cases}
```
```math
\implies \frac{\text{IV}}{\text{II}} = i\kappa = \frac{ k\left(Ce^{ka} - De^{-ka}\right)}{Ce^{ka} + De^{-ka}}
```
```math
\implies i\kappa\left(Ce^{ka} + De^{-ka}\right) = k\left(Ce^{ka} - De^{-ka}\right) 
```
```math
\implies C\left(i\kappa - k\right)e^{ka} = D\left(-i\kappa - k\right)e^{-ka},      \left(\lambda = i\kappa - k\implies i\kappa - k = \lambda^*\right)
```
```math
\implies D = C\frac{\lambda}{\lambda^*}e^{2ka}
```
Now that the relation between D and C is specified, it is time to calculate the relation between other coefficients:
```math
\begin{cases}
Ae^{-i\kappa a} + Be^{i\kappa a} = Ce^{-ka} + C\frac{\lambda}{\lambda^*}e^{3ka} \\
i\kappa\left( Ae^{-i\kappa a} - Be^{i\kappa a}\right) = k\left(Ce^{-ka} - C\frac{\lambda}{\lambda^*}e^{3ka}\right) \\
\end{cases}
```
```math
\implies \frac{i\kappa\left( Ae^{-i\kappa a} - Be^{i\kappa a}\right)}{Ae^{-i\kappa a} + Be^{i\kappa a}} = \frac{k\left(Ce^{-ka} - C\frac{\lambda}{\lambda^*}e^{3ka}\right)}{Ce^{-ka} + C\frac{\lambda}{\lambda^*}e^{3ka}}
```
```math
,      \left( \frac{\left(e^{-ka} - \frac{\lambda}{\lambda^*}e^{3ka}\right)}{e^{-ka} + \frac{\lambda}{\lambda^*}e^{3ka}} = \Lambda \right)
```
```math
\implies i\kappa\left( Ae^{-i\kappa a} - Be^{i\kappa a}\right) = k\Lambda\left(Ae^{-i\kappa a} + Be^{i\kappa a}\right)
```
```math
A\left(i\kappa - k\Lambda\right)e^{-i\kappa a} = B\left(i\kappa + k\Lambda\right)e^{i\kappa a}
```
```math
\implies \frac{B}{A} = \frac{i\kappa - k\Lambda}{i\kappa + k\Lambda}e^{-2i\kappa a}
```
```math
\begin{cases}
Fe^{i\kappa a} = Ce^{ka} + De^{-ka} = Ce^{ka} + C\frac{\lambda}{\lambda^*}e^{ka}  \\
Ae^{-i\kappa a} + \frac{i\kappa - k\Lambda}{i\kappa + k\Lambda}Ae^{-i\kappa a} = Ce^{-ka} + C\frac{\lambda}{\lambda^*}e^{3ka} \\
\end{cases}
```
```math
\implies F = C\left(\frac{\lambda^* + \lambda}{\lambda^*}\right)e^{ka - i\kappa a}
```
```math
\implies A = \frac{C\left(e^{-ka} + \frac{\lambda}{\lambda^*}e^{3ka}\right)}{\frac{2i\kappa}{i\kappa + k\Lambda}e^{-i\kappa a}}
```
```math
\implies \frac{F}{A} = \frac{\left(\frac{\lambda^* + \lambda}{\lambda^*}\right)e^{ka - i\kappa a}}{\frac{e^{-ka} + \frac{\lambda}{\lambda^*}e^{3ka}}{\frac{2i\kappa}{i\kappa + k\Lambda}e^{-i\kappa a}}}
```
```math
= \frac{\left(\frac{2i\kappa}{i\kappa + k\Lambda}\right)\left(\frac{\lambda^* + \lambda}{\lambda^*}\right)e^{ka - 2i\kappa a}}{e^{-ka} + \frac{\lambda}{\lambda^*}e^{3ka}}
```
Now, it is time to get rid of $\lambda$ and $\Lambda$ and make these relations look less complicated.
```math
\frac{B}{A} = \frac{i\kappa - k\Lambda}{i\kappa + k\Lambda}e^{-2i\kappa a} = \left(\frac{i\kappa - \frac{k\left(e^{-ka} - \frac{\lambda}{\lambda^*}e^{3ka}\right)}{e^{-ka} + \frac{\lambda}{\lambda^*}e^{3ka}}}{i\kappa + \frac{k\left(e^{-ka} - \frac{\lambda}{\lambda^*}e^{3ka}\right)}{e^{-ka} + \frac{\lambda}{\lambda^*}e^{3ka}}}\right)e^{-2i\kappa a}
```
```math
= \left[\frac{i\kappa\left(e^{-ka} + \frac{\lambda}{\lambda^*}e^{3ka}\right) - k\left(e^{-ka} - \frac{\lambda}{\lambda^*}e^{3ka}\right)}{i\kappa\left(e^{-ka} + \frac{\lambda}{\lambda^*}e^{3ka}\right) + k\left(e^{-ka} - \frac{\lambda}{\lambda^*}e^{3ka}\right)}\right]e^{-2i\kappa a}
```
```math
= \frac{i\kappa\left[\left(-i\kappa - k\right)e^{-ka} + \left(i\kappa - k\right)e^{3ka}\right] - k\left[\left(-i\kappa - k\right)e^{-ka} - \left(i\kappa - k\right)e^{3ka}\right]}{i\kappa\left[\left(-i\kappa - k\right)e^{-ka} + \left(i\kappa - k\right)e^{3ka}\right] + k\left[\left(-i\kappa - k\right)e^{-ka} - \left(i\kappa - k\right)e^{3ka}\right]}e^{-2i\kappa a}
```
```math
= \frac{\left(\kappa^2 - i\kappa k\right)e^{-ka} + \left(-\kappa^2 - i\kappa k\right)e^{3ka} - \left(-ik\kappa - k^2\right)e^{-ka} + \left(ik\kappa - k^2\right)e^{3ka}}{\left(\kappa^2 - i\kappa k\right)e^{-ka} + \left(-\kappa^2 - i\kappa k\right)e^{3ka} + \left(-ik\kappa - k^2\right)e^{-ka} - \left(ik\kappa - k^2\right)e^{3ka}}e^{-2i\kappa a}
```
```math
\implies \frac{B}{A} = \left[\frac{e^{-ka}\left(\kappa^2 + k^2\right) + e^{3ka}\left(-\kappa^2 - k^2\right)}{e^{-ka}\left(\kappa^2 - k^2 - 2i\kappa k\right) + e^{3ka}\left(k^2 - \kappa^2 - 2i\kappa k\right)}\right]e^{-2i\kappa a}
```
For $\frac{F}{A}$ we have:
```math
 \frac{F}{A} = \frac{\left(\frac{2i\kappa}{i\kappa + k\Lambda}\right)\left(\frac{\lambda^* + \lambda}{\lambda^*}\right)e^{ka}}{e^{-ka} + \frac{\lambda}{\lambda^*}e^{3ka}}e^{- 2i\kappa a}
```
```math
= \frac{\left(\frac{2i\kappa}{i\kappa + k\Lambda}\right)\left(\frac{\lambda^* + \lambda}{\lambda^*}\right)e^{ka}}{\lambda^* e^{-ka} + \lambda e^{3ka}}e^{- 2i\kappa a}
```
```math
= \frac{\left[\frac{2i\kappa}{i\kappa + \frac{k\left(\lambda^* e^{-ka} - \lambda e^{3ka}\right)}{\lambda^* e^{-ka} + \lambda e^{3ka}}}\right]\left(\lambda^* + \lambda\right)e^{ka}}{\lambda^* e^{-ka} + \lambda e^{3ka}}e^{- 2i\kappa a}
```
```math
= \frac{2i\kappa\left(-2k\right)e^{ka}}{i\kappa\left(\lambda^* e^{-ka} + \lambda e^{3ka}\right) + k\left(\lambda^* e^{-ka} - \lambda e^{3ka}\right)}e^{-2i\kappa a}
```
```math
\implies \frac{F}{A} = \left[\frac{2i\kappa\left(-2k\right)e^{ka}}{e^{-ka}\left(\kappa^2 - k^2 - 2i\kappa k\right) + e^{3ka}\left(k^2 - \kappa^2 - 2i\kappa k\right)}\right]e^{-2i\kappa a}
```
The interesting fact is that if you calculate $\left\Vert\frac{F}{A}\right\Vert^2 + \left\Vert\frac{B}{A}\right\Vert^2$, the result will always be 1. That's why we call $\left\Vert\frac{F}{A}\right\Vert^2$ the *Transmission Rate* ($T$) and $\left\Vert\frac{B}{A}\right\Vert^2$ the *Reflection Rate* ($R$). One can interpret this as a free particle that is scattering from left to right and at some point in the space, it hits the wall and after hitting the wall, some of it reflects, and some of it transmits and passes the wall.

Another interesting fact is that despite the fact that $E < V_0$ in this spectrum, yet some part of the quantum particle manages to pass the potential wall!!

In this section we explored the relationship between the coefficients inside the eigenstates of hamiltonian by applying the boundary conditions on them. We calculated the relation between $A$ and $B$, $A$ and $F$, $C$ and $D$ and $C$ and $F$. This means by knowing only one of these coefficients (lets say $A$), all the other coefficients can be calculated using the relations discovered earlier. In the next section, the exact same thing will be calculated but for eigenfunctions which their energies exceed $V_0$.

## Transmission and Reflection Coefficients for E > V0
In this spectrum, the boundary conditions are the same, but the shape of eigenstates is different:

```math
\psi_{\kappa}(x) = \begin{cases}
Ae^{i\kappa x} + Be^{-i\kappa x} & x < -a \\
C\cos(kx) + D\sin(kx) & -a \le x \le a \\
Fe^{i\kappa x} & x > -a \\
\end{cases}
```
```math
\lim_{\epsilon \to 0}
\begin{cases}
\psi_{\kappa}(-a -\epsilon) = \psi_{\kappa}(-a +\epsilon) & (\text{I}) \\
\psi_{\kappa}(+a -\epsilon) = \psi_{\kappa}(+a +\epsilon) & (\text{II})\\
\frac{d}{dx}\psi_{\kappa}(-a -\epsilon) = \frac{d}{dx}\psi_{\kappa}(-a +\epsilon) & (\text{III}) \\
\frac{d}{dx}\psi_{\kappa}(+a -\epsilon) = \frac{d}{dx}\psi_{\kappa}(+a +\epsilon) & (\text{IV}) \\
\end{cases} 
```
```math
\implies \begin{cases}
Ae^{-i\kappa a} + Be^{i\kappa a} = C\cos(-ka) + D\sin(-ka) & (\text{I}) \\
Fe^{i\kappa a} = C\cos(ka) + D\sin(ka) & (\text{II}) \\ 
i\kappa\left( Ae^{-i\kappa a} - Be^{i\kappa a}\right) = k\left(-C\sin(-ka) + D\cos(-ka)\right) & (\text{III}) \\
i\kappa\left(Fe^{i\kappa a}\right) = k\left(-C\sin(ka) + D\cos(ka)\right) & (\text{IV}) \\
\end{cases}
```
```math
\implies \begin{cases}
Ae^{-i\kappa a} + Be^{i\kappa a} = C\cos(ka) - D\sin(ka) & (\text{I}) \\
Fe^{i\kappa a} = C\cos(ka) + D\sin(ka) & (\text{II}) \\ 
i\kappa\left( Ae^{-i\kappa a} - Be^{i\kappa a}\right) = k\left(C\sin(ka) + D\cos(ka)\right) & (\text{III}) \\
i\kappa\left(Fe^{i\kappa a}\right) = k\left(-C\sin(ka) + D\cos(ka)\right) & (\text{IV}) \\
\end{cases}
```
```math
\implies \frac{\text{IV}}{\text{II}} = i\kappa = \frac{k\left(-C\sin(ka) + D\cos(ka)\right)}{C\cos(ka) + Dsin(ka)} \implies i\kappa\left(C\cos(ka) + Dsin(ka)\right) = k\left(-C\sin(ka) + D\cos(ka)\right)
```
```math
\implies C\left(k\sin(ka) + i\kappa\cos(ka)\right) = D\left(k\cos(ka) - i\kappa\sin(ka)\right) \implies \frac{C\left(k\sin(ka) + i\kappa\cos(ka)\right)}{k\cos(ka) - i\kappa\sin(ka)} = D
```
```math
(\text{I}) \rightarrow Ae^{-i\kappa a} + Be^{i\kappa a} = C\cos(ka) - \frac{C\left(k\sin^2(ka) + i\kappa\cos(ka)\sin(ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}
```
```math
\implies Ae^{-i\kappa a} + Be^{i\kappa a} = \frac{C\left(k\cos^2(ka) - i\kappa\sin(ka)cos(ka)\right) - C\left(ksin^(ka) + i\kappa\sin(ka)\cos(ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}
```
```math
= \frac{C\left(k\cos(2ka) - i\kappa\sin(2ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}
```
```math
(\text{III}) \rightarrow i\kappa\left( Ae^{-i\kappa a} - Be^{i\kappa a}\right) = k\left[C\sin(ka) + \frac{C\left(k\sin(ka)\cos(ka) + i\kappa\cos^2(ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}\right]
```
```math
= k\left[\frac{C\left(k\sin(ka)\cos(ka) - i\kappa\sin^2(ka)\right) + C\left(k\sin(ka)\cos(ka) + i\kappa\cos^2(ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}\right] = k\left[\frac{C\left(k\sin(2ka) + i\kappa\cos(2ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}\right]
```
```math
\begin{cases}
i\kappa\left( Ae^{-i\kappa a} - Be^{i\kappa a}\right) = k\left[\frac{C\left(k\sin(2ka) + i\kappa\cos(2ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}\right] \\
i\kappa\left( Ae^{-i\kappa a} + Be^{i\kappa a}\right) = i\kappa\left[\frac{C\left(k\cos(2ka) - i\kappa\sin(2ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}\right] \\
\end{cases}
```
```math
\implies \begin{cases}
2i\kappa\left(Ae^{-i\kappa a}\right) = \frac{C\left(\left(k^2 + \kappa^2 \right)\sin(2ka) + 2i\kappa k\cos(2ka)\right)}{k\cos(ka) - i\kappa\sin(ka)} \\
2i\kappa\left(Be^{i\kappa a}\right) = \frac{C\left(\left(\kappa^2 - k^2 \right)\sin(2ka)\right)}{k\cos(ka) - i\kappa\sin(ka)} \\
\end{cases}
```
```math
\implies \frac{B}{A} = \left[\frac{\left(\kappa^2 - k^2 \right)\sin(2ka)}{\left(k^2 + \kappa^2 \right)\sin(2ka) + 2i\kappa k\cos(2ka)}\right]e^{-2i\kappa a}
```
Now, it is time to calculate $\frac{F}{A}$:
```math
(\text{IV}) \rightarrow i\kappa\left(Fe^{i\kappa a}\right) = k\left(-C\sin(ka) + \frac{C\left(k\sin(ka)\cos(ka) + i\kappa\cos^2(ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}\right)
```
```math
= k\left[\frac{C\left(-k\sin(ka)\cos(ka) + i\kappa\sin^2(ka)\right) + C\left(k\sin(ka)\cos(ka) + i\kappa\cos^2(ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}\right] 
```
```math
\implies i\kappa\left(Fe^{i\kappa a}\right) = \frac{Ci\kappa k}{k\cos(ka) - i\kappa\sin(ka)}
```
```math
\implies \frac{F}{A} = \frac{2i\kappa k}{\left(k^2 + \kappa^2 \right)\sin(2ka) + 2i\kappa k\cos(2ka)}e^{-2i\kappa a}
```

Again, if you compute $\left\Vert\frac{F}{A}\right\Vert^2 + \left\Vert\frac{B}{A}\right\Vert^2$, the result will always be 1. The interesting fact here, is that despite the fact that $E > V_0$ in this spectrum, yet some part of the quantum particle cannot pass the potential wall and gets reflected!!

In this section the *Transmission* and *Reflection* rates were calculated and in the next section, the Normalization factor for eigenstates of hamiltonian will be examined.
## Normalization of Eigenstates
As you may know, for us to be able to work with eigenvectors of hamiltonian and create wave-packets, they should be normalized, meaning:
```math
\langle \psi_{\kappa} | \psi_{\kappa'}\rangle = \delta(\kappa - \kappa')
```
As an example, consider a free particle ($\psi(x) = e^{i\kappa x}$), in this case the normalization factor is calculated like below:
```math
\langle \psi_{\kappa} | \psi_{\kappa'}\rangle = \int_{-\infty}^{+\infty} Ne^{-i\kappa x}Ne^{i\kappa x} dx = N^2\int_{-\infty}^{+\infty} e^{-i\kappa x}e^{i\kappa x} dx = 2\pi\delta(\kappa - \kappa')
```
```math
\implies N = \frac{1}{\sqrt{2\pi}}
```
Therefore, $\frac{1}{\sqrt{2\pi}}e^{i\kappa x}$ is a normalized eigenket. However, the problem here is that if we want to calculate the normalization factor for eigenkets of hamiltonian for this problem, we would encounter so many integrals that may not be evaluated easily. This is because the wavefunction has three parts and each of them should be integrated seperately. Even the first and the third part with exponential terms cannot be evaluated like free particle because the boundaries of integration is not $-\infty \to +\infty$. So, what is the normalization factor in this specific problem?
```math
\psi_\kappa(x) = N \begin{cases}
e^{i\kappa x} + re^{-i\kappa x} & x < -a \\
\text{(some expression at the middle region)} & -a \le x \le a \\
te^{i\kappa x} & x > a
\end{cases}
```
```math
\implies N = ?
```
I'll give the answer straightforward and then give two proofs stating why the answer is correct. First, the answer is $N = \frac{1}{\sqrt{2\pi}}$, but why?

Before diving into the proofs, it is very important to review an interesting fact in quantum mechanics. And that is the probablity current conservation law which comes from the Schrödinger equation itself. 
```math
\partial_t\rho + \nabla \cdot \vec{j} = 0,      (\rho = \Psi^*\Psi,        \vec{j} = \frac{i\hbar}{2m}\nabla\cdot\left(\nabla\Psi^*\Psi - \Psi^*\nabla\Psi\right)
```
```math
\implies \frac{d}{dt}\int \rho dv = -\oint\vec{j}\cdot\hat{n}dA
```
This means if a wavefunction is an answer to the Schrödinger equation, its norm ($\langle\psi|\psi\rangle$) doesn't change over time as long as no probability is destroyed or created. This was the prerequisite to the proofs that are going to be explained below.

### Proof 1
As stated earlier, one can interpret the eigenstates for the finite potential wall problem, as free particles that were moving freely in the space until, at some point, they hit the wall and their shape changes to what we obtained earlier ($\psi_\kappa$). Now, throughout process, no probability is created or destroyed. Therefore, the norm of $\psi_\kappa$ is as same as the norm of $\frac{1}{\sqrt{2\pi}}e^{i\kappa x}$. Hence, the normalization factor that we are looking for is $\frac{1}{\sqrt{2\pi}}$.

One can object this proof because a free particle's variance in position is infinite ($e^{i\kappa x}$). The particle exists everywhere all the time and no one can specify a start and an end for its movement. Hence, we can't assume that a free particle starts moving from an approximate point b and hits the wall at a certain point in time and space. It is already everywhere and obviously not free anymore since we have a potential wall. That's why I came up with another explanation.
### Proof 2
Suppose there is a free particle (absolutely free, meaning there is no potential wall anywhere). It exsits everywhere, but, at some time (lets say t), a potential barrier starts to raise up from the ground at the origin of the space with a width of $2a$. Therefore, the free particle changes its shape to $\psi_kappa$ and, again, since no probability is created or destroyed, the normalization factor remains the same.

Now, another question that may arise is, how we know the changed shape of $\frac{1}{\sqrt{2\pi}}e^{i\kappa x}$ is exactly $\psi_\kappa$? To verify that $\lim_{V_0 \to 0} \psi_\kappa$ should be calculated and checked whether it is equal to $\frac{1}{\sqrt{2\pi}}e^{i\kappa x}$ or not. Suppose there is an eigenket with energy $E$, if $V_0$ is moved towards zero, no matter how small the energy is, at some point $V_0$ will become lower than $E$. Thus, the form of the wavefunction for $\lim_{V_0 \to 0} \psi_\kappa$, is as same the form of eigenstates with $E > V_0$:

```math
\lim_{V_0 \to 0} \frac{1}{\sqrt{2\pi}}\begin{cases}
e^{i\kappa x} + re^{-i\kappa x} & x < -a \\
C\cos(kx) + D\sin(kx) & -a \le x \le a \\
te^{i\kappa x} & x > -a \\
\end{cases}
```
```math
, \kappa = \sqrt{\frac{2m}{\hbar^2}E} , k = \sqrt{\frac{2m}{\hbar^2}\left(E - V_0\right)}
```
```math
\lim_{V_0 \to 0} \sqrt{\frac{2m}{\hbar^2}\left(E - V_0\right)} = \sqrt{\frac{2m}{\hbar^2}E} \implies k = \kappa
```
```math
\lim_{V_0 \to 0} r = \lim_{V_0 \to 0} \left[\frac{\left(\kappa^2 - k^2 \right)\sin(2ka)}{\left(k^2 + \kappa^2 \right)\sin(2ka) + 2i\kappa k\cos(2ka)}\right]e^{-2i\kappa a} = \left[\frac{\left(0\right)\sin(2ka)}{\left(k^2 + \kappa^2 \right)\sin(2ka) + 2i\kappa k\cos(2ka)}\right]e^{-2i\kappa a} = 0
```
```math
2i\kappa\left(Ae^{-i\kappa a}\right) = \frac{C\left(\left(k^2 + \kappa^2\right)\sin(2ka) + 2i\kappa k\cos(2ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}, (A = 1,  V_0 \to 0 \implies k \to \kappa)
```
```math
\implies 2i\kappa\left(e^{-i\kappa a}\right) = \frac{C\left(2\kappa^2\sin(2\kappa a) + 2i\kappa^2\cos(2\kappa a)\right)}{\kappa\cos(\kappa a) - i\kappa\sin(\kappa a)}
```
```math
\implies 2i\kappa\left(e^{-i\kappa a}\right) = \frac{C2i\kappa^2\left(e^{-2i\kappa a}\right)}{\kappa\left(e^{-i\kappa a}\right)} \implies 2i\kappa\left(e^{-i\kappa a}\right) = 2iC\kappa\left(e^{-i\kappa a}\right) \implies C = 1
```
```math
\frac{C\left(k\sin(ka) + i\kappa\cos(ka)\right)}{k\cos(ka) - i\kappa\sin(ka)} = D, (C = 1, V_0 \to 0 \implies k \to \kappa)
```
```math
\frac{i\kappa\left(-i\sin(ka) + cos(ka)\right)}{\kappa\left(\cos(ka) - i\sin(ka)\right)} = D \implies D = i
```
```math
\frac{F}{A} = \frac{2i\kappa k}{\left(k^2 + \kappa^2 \right)\sin(2ka) + 2i\kappa k\cos(2ka)}e^{-2i\kappa a}, (A = 1,  V_0 \to 0 \implies k \to \kappa) 
```
```math
\implies \frac{F}{A} = \frac{2i\kappa^2}{2i\kappa^2\left(-i\sin(2\kappa a) + \cos(2\kappa a)\right)}e^{-2i\kappa a} \implies \frac{2i\kappa^2}{2i\kappa^2e^{-2i\kappa a}}e^{-2i\kappa a} \implies F = 1
```
```math
\lim_{V_0 \to 0} \frac{1}{\sqrt{2\pi}}\begin{cases}
e^{i\kappa x} + re^{-i\kappa x} & x < -a \\
C\cos(kx) + D\sin(kx) & -a \le x \le a \\
te^{i\kappa x} & x > -a \\
\end{cases}
= \frac{1}{\sqrt{2\pi}}\begin{cases}
e^{i\kappa x} + 0e^{-i\kappa x} & x < -a \\
1\cos(\kappa x) + i\sin(\kappa x) & -a \le x \le a \\
1e^{i\kappa x} & x > -a \\
\end{cases} = \frac{1}{\sqrt{2\pi}}e^{i\kappa x}
```
As it is seen, $\psi_\kappa$ converges to $\frac{1}{\sqrt{2\pi}}e^{i\kappa x}$ as $V_0$ goes to zero.

In this section, the normalization factor for eigenkets of hamiltonian was specified, next, a minimum error wave packet is created and the coefficient of each eigenstate ($c(\kappa)$) is calculated. This is necessary because in this simulation all the initial states are minimum error wave packets.
## Minimum Error Wave Packet
A minimum error wave packet is a wave packet whose variance in position and momentum is minimum. According to the Heisenberg uncertainty principle $\sigma^2_x\sigma^2_p \ge \frac{\hbar^2}{4}$ and a minimum error wave packet is a wave packet with $\sigma^2_x = \sigma^2_p = \frac{\hbar}{2}$ and its formula is as follows:

```math
\psi_{MEWP}(x) = \left(\frac{1}{2\pi\sigma^2_x}\right)^{\frac{1}{4}}e^{\frac{-\left(x-\overline{x}\right)^2}{4\sigma^2_x}}e^{\frac{i}{\hbar}\overline{p}x}
```

Having $\Psi(x, 0)$ as $\psi_{MEWP}(x)$, we can calculate $c(\kappa) = \langle \psi_\kappa | \psi_{MEWP} \rangle$ explicitly and use it in the simulation. Moreover, if the initial wave packet starts its movement far enough from the wall (for example $\overline{x} = -a - 12\sigma_x$), the effect of the second and the third part of the eigenfunction is negligible in the integration, meaning:


```math
\langle \psi_\kappa | \psi_{MEWP} \rangle \approx \int_{-\infty}^{+\infty} \left(\frac{1}{2\pi\sigma^2_x}\right)^{\frac{1}{4}}e^{\frac{-\left(x-\overline{x}\right)^2}{4\sigma^2_x}}e^{\frac{i}{\hbar}\overline{p}x} \left(\frac{1}{\sqrt{2\pi}}e^{-i\kappa x} + \frac{r}{\sqrt{2\pi}}e^{i\kappa x}\right) dx
```
```math
= \int_{-\infty}^{+\infty} \left(\frac{1}{2\pi\sigma^2_x}\right)^{\frac{1}{4}}e^{\frac{-\left(x-\overline{x}\right)^2}{4\sigma^2_x}}e^{\frac{i}{\hbar}\overline{p}x}\frac{1}{\sqrt{2\pi}}e^{-i\kappa x} dx + \int_{-\infty}^{+\infty} \left(\frac{1}{2\pi\sigma^2_x}\right)^{\frac{1}{4}}e^{\frac{-\left(x-\overline{x}\right)^2}{4\sigma^2_x}}e^{\frac{i}{\hbar}\overline{p}x}\frac{r}{\sqrt{2\pi}}e^{i\kappa x} dx 
```
```math
= \frac{1}{\sqrt{2\pi}}\left(\frac{1}{2\pi\sigma^2_x}\right)^{\frac{1}{4}}\int_{-\infty}^{+\infty}e^{\left[\frac{-\left(x-\overline{x}\right)^2}{4\sigma^2_x} + i\left(\frac{\overline{p}}{\hbar} - \kappa\right)x\right]} dx + \frac{r}{\sqrt{2\pi}}\left(\frac{1}{2\pi\sigma^2_x}\right)^{\frac{1}{4}}\int_{-\infty}^{+\infty}e^{\left[\frac{-\left(x-\overline{x}\right)^2}{4\sigma^2_x} + i\left(\frac{\overline{p}}{\hbar} + \kappa\right)x\right]} dx
```

Recall:


```math
\int_{-\infty}^{+\infty} e^{-ax^2 + bx} dx = \int_{-\infty}^{+\infty} e^{-a\left(x - \frac{b}{2a}\right)^2 + \frac{b^2}{4a}} dx = e^{\frac{b^2}{4a}}\sqrt{\frac{\pi}{a}}
```

Therefore:

```math
= \frac{1}{\sqrt{2\pi}}\left(\frac{1}{2\pi\sigma^2_x}\right)^{\frac{1}{4}}\int_{-\infty}^{+\infty}e^{\left[\frac{-\left(x-\overline{x}\right)^2}{4\sigma^2_x} + i\left(\frac{\overline{p}}{\hbar} - \kappa\right)\left(x - \overline{x}\right) + i\left(\frac{\overline{p}}{\hbar} - \kappa\right)\overline{x}\right]} dx + \frac{r}{\sqrt{2\pi}}\left(\frac{1}{2\pi\sigma^2_x}\right)^{\frac{1}{4}}\int_{-\infty}^{+\infty}e^{\left[\frac{-\left(x-\overline{x}\right)^2}{4\sigma^2_x} + i\left(\frac{\overline{p}}{\hbar} + \kappa\right)\left(x - \overline{x}\right) + i\left(\frac{\overline{p}}{\hbar} + \kappa\right)\overline{x}\right]} dx
```

```math
\implies \langle \psi_\kappa | \psi_{MEWP} \rangle = \frac{e^{i\left(\frac{\overline{p}}{\hbar} - \kappa\right)\overline{x}}}{\sqrt{2\pi}}\left(\frac{1}{2\pi\sigma^2_x}\right)^{\frac{1}{4}}e^{\left[i\left(\frac{\overline{p}}{\hbar} - \kappa\right)\right]^2\sigma^2_x}\sqrt{4\pi\sigma^2_x} + \frac{re^{i\left(\frac{\overline{p}}{\hbar} + \kappa\right)\overline{x}}}{\sqrt{2\pi}}\left(\frac{1}{2\pi\sigma^2_x}\right)^{\frac{1}{4}}e^{\left[i\left(\frac{\overline{p}}{\hbar} + \kappa\right)\right]^2\sigma^2_x}\sqrt{4\pi\sigma^2_x}
```

```math
c(\kappa) = \langle \psi_\kappa | \psi_{MEWP} \rangle = \left[\frac{e^{i\left(\frac{\overline{p}}{\hbar} - \kappa\right)\overline{x}}}{\sqrt[4]{2\pi\left(\frac{1}{4\sigma^2_x}\right)}}\right]e^{\frac{-\left(\frac{\overline{p}}{\hbar} - \kappa\right)^2}{4\left(\frac{1}{4\sigma^2_x}\right)}} + \left[\frac{re^{i\left(\frac{\overline{p}}{\hbar} + \kappa\right)\overline{x}}}{\sqrt[4]{2\pi\left(\frac{1}{4\sigma^2_x}\right)}}\right]e^{\frac{-\left(\frac{\overline{p}}{\hbar} + \kappa\right)^2}{4\left(\frac{1}{4\sigma^2_x}\right)}}
```


Using this formula, the coefficient of each eigenstate in the superposition can be calculated.


## Conclusion
In this article, the Finite Potential Wall problem was studied and analysed. The eigenstates of hamiltonian were calculated along with the relationship between the coefficients inside each. The normalization factor was specified and two proofs were given for that. After that, a minimum error wave packet was created as the initial state, then, the coefficients of each eigenfunction inside the superposition was explicitly evaluated. Knowing all of these equations and formulas, enables us to simulate the movement of a minimum error wave packet (with average monemtum of $\overline{p}$) along with its interaction with the potential wall demonstrating the quantum tunneling phenomenon.


## Simple Demonstration
For veiw a run sample of the simulation code you can see my [linkedin post](). Feel free to share your opinion on this project with me. Thank you for reading this article, hope you've enjoyed this simulation project!