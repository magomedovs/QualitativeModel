# Slow deflagration solution

The following nonlinear eigenvalue boundary value problem is solved
$$\begin{align}
    \Theta_{0 \rho \rho} + c_1 \Theta_{0 \rho} &\left(1 + \dfrac{\epsilon}{2\left(1 + \epsilon U_0(\Theta_{0}) \right)} \right) = - q \omega_0, \\
    &\dfrac{1}{\mathrm{Le}} \Lambda_{0 \rho \rho} + c_1 \Lambda_{0 \rho} = -\omega_0.
 \end{align}$$
where $\omega_0 = k (1 - \Lambda_0) \exp(-\theta / \Theta_0) \, \mathrm{H}(\Theta_0 - T_{\mathrm{ign}})$, $U_0(\Theta_0) = \left(-1 + \sqrt{1 - \epsilon^2 \Theta_0 + 2 \epsilon u_r + \epsilon^2 \left( u_r^2 + T_{0 r} \right)} \right) / \epsilon$, and $c_1$ is an unknown eigenvalue.

Dirichlet boundary conditions are imposed ($T_{0r}$ is computed, $T_l$ is chosen, $\lambda(-\infty) = 1$, $\lambda(-\infty) = 0$) with an additional constraint $\Theta_0(0) = T_{\mathrm{ign}}$.