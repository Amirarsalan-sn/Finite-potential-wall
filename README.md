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
C\cos(kx) + D\sin(-kx) & -a \le x \le a \\
Fe^{i\kappa x} + Ge^{-i\kappa x} & x > -a \\
\end{cases}
```

In this section, we found the eigenfunctions of hamiltonian which form an orthonormal basis for our hilbert space. In the next section, the transmission and reflection coefficients are going to be calculated by applying the boundry conditions to $\psi_{\kappa} (x)$.

## Transmission and Reflection Coefficients for E < V0
Before diving deep into calculations, it is worth noting that singularities and discontinuities in $\psi_{\kappa}$ is always two steps better than singularities and discontinuities in $V(x)$ and that's because of the time independent Schrödinger equation:

```math
\frac{d^2}{dx^2}\psi(x) = \frac{-2m}{\hbar^2}\left(E - V(x)\right)\psi(x)
```
As it is seen from the equation, in order to get $\psi_{\kappa}$, one should integrate $\frac{-2m}{\hbar^2}\left(E - V(x)\right)\psi(x)$ two times. Hence, if the potential has simple discontinuities (like the finite potential wall), $\psi_{\kappa}$ and $\frac{d}{dx}\psi_{\kappa}$ should be continueous, and that's what gives us the *boundry conditions*. With these conditions we are able to calculate the transmission and reflection coefficients.

Another important thing to note is that we normaly assume the eigenfunctions to scatter from left to right (or right to left) to make the calculation of coefficients easier and more interpretable. Without this assumption, it would be to hard too find any relationship between the coefficients of $\psi$ just by applying the boundry conditions. Thus, the value of $G$ is assumed to be zero throughout the calculations and simulations.

The eigenstates and boundry conditions for $E < V_0$ are:

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

In this section we explored the relationship between the coefficients inside the eigenstates of hamiltonian by applying the boundry conditions on them. We calculated the relation between $A$ and $B$, $A$ and $F$, $C$ and $D$ and $C$ and $F$. This means by knowing only one of these coefficients (lets say $A$), all the other coefficients can be calculated using the relations discovered earlier. In the next section, the exact same thing will be calculated but for eigenfunctions which their energies exceed $V_0$.

## Transmission and Reflection Coefficients for E > V0
In this spectrum, the boundry conditions are the same, but the shape of eigenstates is different:

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
\implies Ae^{-i\kappa a} + Be^{i\kappa a} = \frac{C\left(k\cos^2(ka) - i\kappa\sin(ka)cos(ka)\right) - C\left(ksin^(ka) + i\kappa\sin(ka)\cos(ka)\right)}{k\cos(ka) - i\kappa\sin(ka)} = \frac{C\left(k\cos(2ka) - i\kappa\sin(2ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}
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
= k\left[\frac{C\left(-k\sin(ka)\cos(ka) + i\kappa\sin^2(ka)\right) + C\left(k\sin(ka)\cos(ka) + i\kappa\cos^2(ka)\right)}{k\cos(ka) - i\kappa\sin(ka)}\right] \implies i\kappa\left(Fe^{i\kappa a}\right) = \frac{Ci\kappa k}{k\cos(ka) - i\kappa\sin(ka)}
```
```math
\implies \frac{F}{A} = \frac{2i\kappa k}{\left(k^2 + \kappa^2 \right)\sin(2ka) + 2i\kappa k\cos(2ka)}e^{-2i\kappa a}
```

Again, if you compute $\left\Vert\frac{F}{A}\right\Vert^2 + \left\Vert\frac{B}{A}\right\Vert^2$, the result will always be 1. The interesting fact here, is that despite the fact that $E > V_0$ in this spectrum, yet some part of the quantum particle cannot pass the potential wall and gets reflected!!

In this section the *Transmission* and *Reflection* rates were calculated and in the next section, the Normalization factor for eigenstates of hamiltonian will be examined.
## Normalization of Eigenstates

## Minimum Error Wave Packet

## Simple Demonstration

## Conclusion
