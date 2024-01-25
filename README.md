# A qualitative model for detonations and flames

Solver for computation of the traveling-wave solutions of the following system [1]

$$\begin{align}
  u_t + u_x + \epsilon u u_x &= -\frac{\epsilon}{2} T_x + \frac{4 \delta}{3 \epsilon} \Pr u_{xx}, \\
  T_t - u_t &= q \omega + \delta T_{xx}, \\
  \lambda_t &= \omega + \frac{\delta}{\mathrm{Le}} \lambda_{xx}.
 \end{align}$$

Here, $u$ is a pressure, $T$ – temperature, $\lambda$ – the reaction-progress variable, varying smoothly from $0$ in the fresh mixture to $1$ in the products, and $\omega = k (1 - \lambda) \exp(-\theta / T) \, \mathrm{H}(T - T_{\mathrm{ign}})$ is a reaction rate.

The stationary solution is computed using an open source finite element library [deal.II](https://www.dealii.org).

[1] [Faria, L. M., Lau-Chapdelaine, S., Kasimov, A. R., & Rosales, R. R. (2017). A toy model for detonations and flames. 26th International Colloquium on the Dynamics of Explosions and Reactive Systems (ICDERS).](http://www.icders.org/ICDERS2017/abstracts/ICDERS2017-1122.pdf)