# Quantum information and computing
## Exercise 4: Time Independent Schr&ouml;dinger equation
### Evolve numerically a system given the Hamiltonian:

$$\hat{H}=\hat{p}^2+\omega^2\hat{x}^2$$

1. $\hat{p}\to -i\hbar\frac{\partial}{\partial x}, \quad \hat{x}\to x$

2. From $\hat{H}\psi = E\psi$ you get: $$\left(-\frac{\hbar^2}{2}\frac{\partial^2}{\partial x^2} + \frac{1}{2}\omega^2 x^2\right)\psi(x) = E_n \psi(x)$$

$\psi(x)\equiv\psi_x$

3. The second derivative must be discretized applying finite difference method:
$$\begin{cases}
\psi_{k+1}=\psi_{k}+\psi^\prime_{k}dx+\frac{1}{2}\psi^{\prime\prime}_{k}dx^2+\frac{1}{6}\psi^{\prime\prime\prime}_k+ O(dx^4)\\
\psi_{k-1}=\psi_{k}-\psi^\prime_{k}dx+\frac{1}{2}\psi^{\prime\prime}_{k}dx^2-\frac{1}{6}\psi^{\prime\prime\prime}_k+ O(dx^4)
\end{cases}$$
$$\psi_k^{\prime\prime}=\frac{\psi_{k+1}-2\psi_{k}+\psi_{k-1}}{dx^2}$$
So we get:
$$-\frac{\hbar^2}{2}\left[\frac{\psi_{k+1}-2\psi_k + \psi_{k-1}}{dx^2}\right] + \frac{1}{2}\omega^2 x_k\psi_k = E \psi_k$$

4. Considering $H_{ij} =\,<x_i|H|x_j>$ the solution of the equation is equivalent to consider the following matrix in a eigenvalue problem:

$$H = \frac{\hbar^2}{2}\left(\begin{matrix}
\frac{2}{dx^2} + \frac{1}{2}\omega^2x_1^2 & \frac{1}{dx^2} & 0 & \dotsm & 0 \\
&&&&\\
\frac{1}{dx^2} & \frac{2}{dx^2} + \frac{1}{2}\omega^2x_2^2 & \frac{1}{dx^2} & \dotsm & 0 \\
&&&&\\
0 & \frac{1}{dx^2} &  \frac{2}{dx^2} + \frac{1}{2}\omega^2x_3^2 & \dotsm & 0 \\
\vdots & \vdots & \vdots &\ddots  & \vdots\\
0 & 0 & 0 & \dotsm  &  \frac{2}{dx^2} + \frac{1}{2}\omega^2x_N^2
\end{matrix}\right)$$


5. To get the eigenvalues this matrix must be then diagonalized.

