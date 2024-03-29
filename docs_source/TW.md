* Intoduction
  Here, we formulate the artificial viscosity as proposed by Tscharnuter &
  Winkler 1979 (https://ui.adsabs.harvard.edu/abs/1979CoPhC..18..171T/abstract)
  for polar coordinate system. Additionally, we implement the proposed changes
  from Stone & Norman 1992 (Appendix B)
  (https://ui.adsabs.harvard.edu/abs/1992ApJS...80..753S/abstract) that read:

  1:
  the artificial viscosity shock smoothing length must be the same in all
  directions ( $l = c \cdot \mathrm{max}(\Delta r, r \Delta \varphi)$) To keep
  the artificial viscosity tensor isotropic and trace free.

  2:
  The off-diagonal (shear) components are set to zero to prevent artificial
  angular momentum transport and to not broaden large gradients of shear.

* Initialize spherical

  coordinates : $r, \theta, z\;\;$ $r \in [0, inf)$, $\theta \in [0, 2\pi)$, $z \in [0, inf)$
  
  $u^l_{;k}$ is the covariant derivative of $u_l$.
  
  $u_{(k;l)} = \frac{1}{2} (u_{k;l} + u_{l;k})$

  Metric:
  $g_{ik} = [1, r^2, 1]$
  
  $g^{ik} = [1, 1/r^2, 1]$

  $u_k = (u, rv, w)$
  
  $u^k = (u, v/r, w)$

  Christoffel symbols:
  
  $\Gamma^r_{\phi\phi} = -r$
  
  $\Gamma^\theta_{r\theta} = \Gamma^\theta_{\theta r} = \frac{1}{r}$

* symmetrized gradient of the velocity field

  $\epsilon_{\;l}^{m} = g^{mk}u_{(k;l)}$
  
  $\epsilon_{\;l}^{m} = g^{mk} \frac{1}{2}[u_{k;l} + u_{l;k}]$

  $u_{k;l} = \frac{\partial u_k}{\partial x^l} - \Gamma^i_{kl}u_i$
 
** $\epsilon_{\;\theta}^{\theta}$

$$
\begin{align}
  \epsilon^{\theta}_{\;\theta} &= g^{\theta\theta} \frac{1}{2}\left[\frac{\partial u_\theta}{\partial \theta} - \Gamma^i_{\theta\;\theta} u_i + \frac{\partial u_\theta}{\partial \theta} - \Gamma^i_{\theta\theta} u_i\right] \\
  &= g^{\theta\theta} \frac{1}{2}[\frac{\partial u_\theta}{\partial \theta} - \Gamma^r_{\theta\;\theta} u_r + \frac{\partial u_\theta}{\partial \theta} - \Gamma^r_{\theta\theta} u_r] \\
  &= \frac{1}{r^2} \frac{1}{2}[\frac{\partial u_\theta}{\partial \theta} - (-r) u_r + \frac{\partial u_\theta}{\partial \theta} - (-r) u_r] \\
  &= \frac{1}{r^2} \frac{1}{2}[\frac{\partial rv}{\partial \theta} + r u + \frac{\partial vr}{\partial \theta} + r u] \\
  &= \frac{1}{r} [\frac{\partial v}{\partial \theta} + u]
\end{align}
$$

** $\epsilon_{\;\theta}^{r}$
$$
   \begin{align}
  \epsilon_{\;\theta}^{r} &= g^{rr} \frac{1}{2}[\frac{\partial u_\theta}{\partial r} - \Gamma^i_{\theta r} u_i + \frac{\partial u_r}{\partial \theta} - \Gamma^i_{r\theta} u_i] \\
  &= \frac{1}{2}[\frac{\partial rv}{\partial r} - \frac{1}{r} rv + \frac{\partial u_r}{\partial \theta} - \frac{1}{r} vr] \\
  &= \frac{1}{2}[v\frac{\partial r}{\partial r} + r\frac{\partial v}{\partial r} - \frac{1}{r} rv + \frac{\partial u_r}{\partial \theta} - \frac{1}{r} vr] \\
  &= \frac{r}{2}[\frac{\partial v}{\partial r} + \frac{1}{r}\frac{\partial u}{\partial \theta} - \frac{v}{r}]
\end{align}
$$

** $\epsilon_{\;r}^{\theta}$
$$
\begin{align}
  \epsilon_{\;r}^{\theta} &= g^{\theta\theta} \frac{1}{2}[\frac{\partial u_r}{\partial \theta} - \Gamma^i_{r\theta} u_i + \frac{\partial u_\theta}{\partial r} - \Gamma^i_{\theta\;r} u_i] \\
  &= g^{\theta\theta} \frac{1}{2}[\frac{\partial u_r}{\partial \theta} - \Gamma^\theta_{r\theta} u_\theta + \frac{\partial u_\theta}{\partial r} - \Gamma^\theta_{\theta\;r} u_\theta] \\
  &= \frac{1}{r^2} \frac{1}{2}[\frac{\partial u}{\partial \theta} - \frac{1}{r} rv + \frac{\partial rv}{\partial r} - \frac{1}{r} rv] \\
  &= \frac{1}{r^2} \frac{1}{2}[\frac{\partial u}{\partial \theta} - \frac{1}{r} rv + \frac{r\partial v}{\partial r} + \frac{v\partial r}{\partial r} - \frac{1}{r} rv] \\
  &= \frac{1}{2r}[\frac{1}{r}\frac{\partial u}{\partial \theta} - \frac{v}{r} + \frac{\partial v}{\partial r}]
\end{align}
$$

  
  
** $\epsilon_{\;r}^{r}$
$$ 
\begin{align}
  \epsilon^{r}_{\;r} &= g^{rr} \frac{1}{2}[\frac{\partial u_r}{\partial r} - \Gamma^i_{r\;r} u_i + \frac{\partial u_r}{\partial r} - \Gamma^i_{rr} u_i] \\
  &= 1 \frac{1}{2}[\frac{\partial u}{\partial r} - 0 u + \frac{\partial u}{\partial r} - 0 u] \\
  &= \frac{\partial u}{\partial r}
\end{align}
$$

** Divergence
   $\nabla \vec{\mathrm{v}} = \mathrm{div}(u) = (\frac{\partial}{\partial r} + \frac{1}{r}) u + \frac{1}{r}\frac{\partial}{\partial \theta} v$

* Viscous Pressure Tensor
$Q^{m}_{\;l} = l^2 \rho \mathrm{div}(u)[\epsilon^{\;m}_l - \frac{1}{3} \mathrm{div}(u) \delta^{\;m}_l]$ if $\mathrm{div}(u) < 0$


$Q_{\;r}^{r} = l^2 \rho \mathrm{div}(u) [\frac{\partial u}{\partial r} - \frac{1}{3}\mathrm{div}(u)]$

$Q_{\;\theta}^{\theta} = l^2 \rho \mathrm{div}(u) [\frac{1}{r} [\frac{\partial v}{\partial \theta} + u] - \frac{1}{3}\mathrm{div}(u)]$

$Q_{\;z}^{z} = l^2 \rho \mathrm{div}(u) [- \frac{1}{3}\mathrm{div}(u)]$

$Q^{r}_{\;\theta} = l^2 \rho \mathrm{div}(u)[\frac{r}{2}[\frac{\partial v}{\partial r} + \frac{1}{r}\frac{\partial u}{\partial \theta} - \frac{v}{r}]]$

$Q^{\theta}_{\;r} = l^2 \rho \mathrm{div}(u)[\frac{1}{2r}[\frac{1}{r}\frac{\partial u}{\partial \theta} - \frac{v}{r} + \frac{\partial v}{\partial r}]]$


We ignore off diagonal terms to prevent artificial angular momentum transfer:

$Q_{\;r}^{r} = l^2 \rho \mathrm{div}(u) [\frac{\partial u}{\partial r} - \frac{1}{3}\mathrm{div}(u)]$

$Q_{\;\theta}^{\theta} = l^2 \rho \mathrm{div}(u) [\frac{1}{r} [\frac{\partial v}{\partial \theta} + u] - \frac{1}{3}\mathrm{div}(u)]$

$Q_{\;z}^{z} = l^2 \rho \mathrm{div}(u) [- \frac{1}{3}\mathrm{div}(u)]$

$Q^{r}_{\;\theta} = 0$

$Q^{\theta}_{\;r} = 0$


$Q_{\;r}^{r} = l^2 \rho \mathrm{div}(u) [\epsilon^{r}_{\;r} - \frac{1}{3}\mathrm{div}(u)]$

$Q_{\;\theta}^{\theta} = l^2 \rho \mathrm{div}(u) [\epsilon^{\theta}_{\;\theta} - \frac{1}{3}\mathrm{div}(u)]$


$\epsilon^r_{\; r} + \epsilon^\theta_{\; \theta} =  \frac{\partial u}{\partial r} + \frac{1}{r} \frac{\partial v}{\partial \theta} + \frac{u}{r} = \mathrm{div}(u)$


$Q_{\;r}^{r} = l^2 \rho \mathrm{div}(u) [\epsilon^{r}_{\;r} - \frac{1}{3}(\epsilon^r_{\; r} + \epsilon^\theta_{\; \theta})]$

$Q_{\;\theta}^{\theta} = l^2 \rho \mathrm{div}(u) [\epsilon^{\theta}_{\;\theta} - \frac{1}{3}(\epsilon^r_{\; r} + \epsilon^\theta_{\; \theta})]$

* Azimuthal Viscous Forces
$\mathrm{div}(u) = (\frac{\partial}{\partial r} + \frac{1}{r}) u + \frac{1}{r}\frac{\partial}{\partial \theta} v$


$\Gamma^r_{\theta\theta} = -r$

$\Gamma^\theta_{r\theta} = \Gamma^\theta_{\theta r} = \frac{1}{r}$


$Q^k_{\; i;k} = Q^k_{\;i, k} + \Gamma^k_{\mu k} Q^\mu_{\;i} - \Gamma^\mu_{i k} Q^k_{\; \mu}$

$Q^k_{\;i, k} = \frac{\partial Q^k_{\; i}}{\partial x_k}$

$$
\begin{align}
v_Q &= Q^k_{\;\theta;k}\\
 &= \frac{\partial Q^k_{\; \theta}}{\partial x_k} + \Gamma^k_{\mu k} Q^\mu_{\;\theta} - \Gamma^\mu_{\theta k} Q^k_{\; \mu} \\
 &= \frac{\partial Q^k_{\; \theta}}{\partial x_k} + \Gamma^k_{\mu k} Q^\mu_{\;\theta} - \Gamma^\mu_{\theta k} Q^k_{\; \mu} \\
 &= \frac{\partial Q^\theta_{\; \theta}}{\partial x_\theta} + \Gamma^k_{\theta k} Q^\theta_{\;\theta} - \Gamma^r_{\theta \theta} Q^\theta_{\; r} - \Gamma^\theta_{\theta r} Q^k_{\; \theta} \\
 &= \frac{\partial Q^\theta_{\; \theta}}{\partial \theta}\\
\end{align}
$$

** Azimuthal Velocity update
$\frac{\partial w}{\partial t} = \frac{1}{\rho} \frac{1}{r} u_Q$

$\frac{\partial w}{\partial t} = \frac{1}{\rho} \frac{1}{r} \frac{\partial Q^\theta_{\; \theta}}{\partial \theta}$


*** Finite difference conservative formulation
$\int^{R_s^{i}}_{R_a^{i}}\int^{\Delta \theta / 2}_{-\Delta \theta / 2}\frac{\partial w}{\partial t} r \mathrm{d}r \mathrm{d}\theta = \frac{1}{2}((R_s^i)^2 - (R_a^i)^2) \Delta \theta \frac{\partial w}{\partial t}$

$\int^{R_s^{i}}_{R_a^{i}}\int^{\Delta \theta / 2}_{-\Delta \theta / 2} \frac{1}{\rho} \frac{1}{r} \frac{\partial Q^\theta_{\; \theta}}{\partial \theta} r \mathrm{d}r \mathrm{d}\theta$

$\int^{\Delta \theta / 2}_{-\Delta \theta / 2} \frac{1}{\rho} \frac{\partial Q^\theta_{\; \theta}}{\partial \theta} \mathrm{d}\theta \Delta r$

$\frac{1}{2}((R_s^i)^2 - (R_a^i)^2) \Delta \theta \frac{\partial w}{\partial t} = \frac{1}{\rho} \Delta Q^\theta_{\; \theta} \Delta r$

$\frac{1}{2}(R_s^i - R_a^i) (R_s^i + R_a^i) \Delta \theta \frac{\partial w}{\partial t} = \frac{1}{\rho} \Delta Q^\theta_{\; \theta} \Delta r$

$\frac{\partial w}{\partial t} = \frac{2}{R_s^i + R_a^i} \frac{1}{\rho} \frac{\Delta Q^\theta_{\; \theta}}{\Delta \theta}$

* Radial Viscous Forces
$$
\begin{align}
u_Q &= Q^k_{\;r;k}\\
 &= \frac{\partial Q^k_{\; r}}{\partial x_k} + \Gamma^k_{\mu k} Q^\mu_{\;r} - \Gamma^\mu_{r k} Q^k_{\; \mu} \\
 &=  \frac{\partial Q^r_{\; r}}{\partial x_r} + \frac{\partial Q^\theta_{\; r}}{\partial x_\theta} + \Gamma^k_{r k} Q^r_{\;r} + \Gamma^k_{\theta k} Q^\theta_{\;r} - \Gamma^\theta_{r \theta} Q^\theta_{\; \theta}\\
 &=  \frac{\partial Q^r_{\; r}}{\partial x_r} + \frac{\partial Q^\theta_{\; r}}{\partial x_\theta} + \Gamma^\theta_{r \theta} Q^r_{\;r} + 0 Q^\theta_{\;r} - \Gamma^\theta_{r \theta} Q^\theta_{\; \theta} \\
 &=  \frac{\partial Q^r_{\; r}}{\partial x_r} + \frac{\partial 0}{\partial x_\theta} + \Gamma^\theta_{r \theta} Q^r_{\;r} + 0 0 - \Gamma^\theta_{r \theta} Q^\theta_{\; \theta} \\
 &=  \frac{\partial Q^r_{\; r}}{\partial x_r} + \frac{1}{r} Q^r_{\;r} - \frac{1}{r} Q^\theta_{\; \theta} \\
%# %&=  \frac{\partial Q^r_{\; r}}{\partial x_r} + \frac{1}{r} l^2 \rho \mathrm{div}(u) [\frac{\partial u}{\partial r} - \frac{1}{3}\mathrm{div}(u)] - \frac{1}{r} l^2 \rho \mathrm{div}(u) [\frac{1}{r} [\frac{\partial v}{\partial \theta} + u] - \frac{1}{3}\mathrm{div}(u)] \\
%# %&=  \frac{Q_{\;r}^{r}}{\partial r} - \frac{1}{r} l^2 \rho \mathrm{div}(u) [\frac{1}{r} \frac{\partial v}{\partial \theta} + \frac{u}{r} - \frac{\partial u}{\partial r}] \\
%# %&=  \frac{Q_{\;r}^{r}}{\partial r} - \frac{1}{r} l^2 \rho \mathrm{div}(u) [\frac{1}{r} \frac{\partial v}{\partial \theta} + \frac{u}{r} + \frac{\partial u}{\partial r} - 2\frac{\partial u}{\partial r}] \\
%# %&=  \frac{Q_{\;r}^{r}}{\partial r} - \frac{1}{r} l^2 \rho \mathrm{div}(u) [\mathrm{div}(u) - 2\frac{\partial u}{\partial r}] \\
\end{align}
$$

** Radial Velocity update
$\frac{\partial u}{\partial t} = \frac{1}{\rho} u_Q$

$\frac{\partial u}{\partial t} = \frac{1}{\rho} [\frac{\partial Q^r_{\; r}}{\partial x_r} + \frac{1}{r} Q^r_{\;r} - \frac{1}{r} Q^\theta_{\; \theta}]$


$\int^{R_b^{i}}_{R_b^{i-1}}\int^{\Delta \theta / 2}_{-\Delta \theta / 2}\frac{\partial u}{\partial t} r \mathrm{d}r \mathrm{d}\theta = \int^{R_b^{i}}_{R_b^{i-1}}\int^{\Delta \theta / 2}_{-\Delta \theta / 2}\frac{1}{\rho} [\frac{\partial Q^r_{\; r}}{\partial x_r} + \frac{1}{r} Q^r_{\;r} - \frac{1}{r} Q^\theta_{\; \theta}] r \mathrm{d}r \mathrm{d}\theta$

$\int^{R_b^{i}}_{R_b^{i-1}}\int^{\Delta \theta / 2}_{-\Delta \theta / 2}\frac{\partial u}{\partial t} r \mathrm{d}r \mathrm{d}\theta = \frac{\partial u}{\partial t} \frac{1}{2}((R_b^i)^2 - (R_b^{i-1})^2) \Delta \theta$


$\int^{R_b^{i}}_{R_b^{i-1}} \frac{\partial Q^r_{\; r}}{\partial x_r} r \mathrm{d}r$

$\frac{\partial Q r}{\partial r} = \frac{\partial Q}{\partial r} r + Q \frac{\partial r}{\partial r}$


$\int\frac{\partial Q r}{\partial r} \mathrm{d}r = \int\frac{\partial Q}{\partial r} r \mathrm{d}r + \int Q \frac{\partial r}{\partial r} \mathrm{d}r$


$\Delta (Q r) = \int\frac{\partial Q}{\partial r} r \mathrm{d}r + Q \Delta r$

$\int\frac{\partial Q}{\partial r} r \mathrm{d}r = \Delta (Q r) - Q \Delta r$

$\int\frac{Q}{r} r \mathrm{d}r = \int Q \mathrm{d}r = Q \Delta r$

$\int^{R_b^{i}}_{R_b^{i-1}}\frac{1}{\rho} [\frac{\partial Q^r_{\; r}}{\partial x_r} + \frac{1}{r} Q^r_{\;r} - \frac{1}{r} Q^\theta_{\; \theta}] r \mathrm{d}r \Delta \theta = \frac{\Delta \theta}{\rho}[\Delta (Q^{r}_{\; r} r) - Q^{r}_{\; r} \Delta r + Q^{r}_{\; r} \Delta r - Q^{\theta}_{\; \theta} \Delta r]$

$=\frac{\Delta \theta}{\rho}[\Delta (Q^{r}_{\; r} r) - Q^{\theta}_{\; \theta} \Delta r]$


$\frac{\partial u}{\partial t} \frac{1}{2}((R_b^i)^2 - (R_b^{i-1})^2) \Delta \theta = \frac{\Delta \theta}{\rho}[\Delta (Q^{r}_{\; r} r) - Q^{\theta}_{\; \theta} \Delta r]$


$\frac{\partial u}{\partial t}  = \frac{2}{(R_b^i)^2 - (R_b^{i-1})^2} \frac{1}{\rho}[\Delta (Q^{r}_{\; r} r) - Q^{\theta}_{\; \theta} \Delta r]$


* Artificial Energy Dissipation

  $\epsilon_q = - \frac{1}{\rho} Q^i_{\; k} \epsilon^k_{\; i}$
  
  $\epsilon_q = - \frac{1}{\rho} [Q^r_{\; r} \epsilon^r_{\; r} + Q^\theta_{\; \theta} \epsilon^\theta_{\; \theta} + Q^z_{\; z} \epsilon^z_{\; z}]$
  
  $Q^r_{\; r} \epsilon^r_{\; r} = l^2 \rho \mathrm{div}(u) [(\epsilon^r_{\; r})^2 - \epsilon^r_{\; r}\frac{1}{3}\mathrm{div}(u)]$
  
  $Q^\theta_{\; \theta} \epsilon^\theta_{\; \theta} = l^2 \rho \mathrm{div}(u) [(\epsilon^\theta_{\; \theta})^2 - \epsilon^\theta_{\; \theta}\frac{1}{3}\mathrm{div}(u)]$

  
  $Q^z_{\; z} \epsilon^z_{\; z} = 0$

  $l^2 = C_2^2 \delta x^2$
  
  $C_2 = l / \delta x$
  

  $\epsilon_q = - \frac{1}{\rho} [l^2 \rho \mathrm{div}(u) [(\epsilon^r_{\; r})^2 - \epsilon^r_{\; r}\frac{1}{3}\mathrm{div}(u)] + l^2 \rho \mathrm{div}(u) [(\epsilon^\theta_{\; \theta})^2 - \epsilon^\theta_{\; \theta}\frac{1}{3}\mathrm{div}(u)]]$
  
  $\epsilon_q = - \frac{1}{\rho} l^2 \rho \mathrm{div}(u)[(\epsilon^r_{\; r})^2 - \epsilon^r_{\; r}\frac{1}{3}\mathrm{div}(u) + (\epsilon^\theta_{\; \theta})^2 - \epsilon^\theta_{\; \theta}\frac{1}{3}\mathrm{div}(u)]$
  
  $\epsilon_q = - l^2 \mathrm{div}(u)[(\epsilon^r_{\; r})^2 + (\epsilon^\theta_{\; \theta})^2 -\frac{1}{3}\mathrm{div}(u) (\epsilon^r_{\; r} + \epsilon^\theta_{\; \theta})]$
  
  Last part are off diagonal terms, which we want to ignore:

  The formulation shown in TW produces:
  $\epsilon_q = - l^2 \mathrm{div}(u)\{\frac{1}{3}[3(\epsilon^r_{\; r})^2 + 3(\epsilon^\theta_{\; \theta})^2 - \mathrm{div}(u)^2]\}$

  But they say it is important to write it as 
  $\epsilon_q = - l^2 \mathrm{div}(u)\frac{1}{3}\{(\epsilon^r_{\; r} - \epsilon^\theta_{\; \theta})^2 + (\epsilon^r_{\; r})^2 + (\epsilon^\theta_{\; \theta})^2\}$
  which is numerically guaranteed to have the correct sign.
  Here, we show that these two formulations are equal.
  
  $\epsilon^r_{\; r} + \epsilon^\theta_{\; \theta} =  \frac{\partial u}{\partial r} + \frac{1}{r} \frac{\partial v}{\partial \theta} + \frac{u}{r} = \mathrm{div}(u)$
  
  $\epsilon_q = - l^2 \mathrm{div}(u)\frac{1}{3}\{(\epsilon^r_{\; r} - \epsilon^\theta_{\; \theta})^2 + (\epsilon^r_{\; r})^2 + (\epsilon^\theta_{\; \theta})^2\}$
  
  $\epsilon_q = - l^2 \mathrm{div}(u)\{\frac{1}{3}[2(\epsilon^r_{\; r})^2 + 2(\epsilon^\theta_{\; \theta})^2 - 2\epsilon^r_{\; r} \epsilon^\theta_{\; \theta}]\}$
  
  $\epsilon_q = - l^2 \mathrm{div}(u)\{\frac{1}{3}[3(\epsilon^r_{\; r})^2 + 3(\epsilon^\theta_{\; \theta})^2 - ((\epsilon^r_{\; r})^2 + (\epsilon^\theta_{\; \theta})^2 + 2\epsilon^r_{\; r} \epsilon^\theta_{\; \theta})]\}$
  
  $\epsilon_q = - l^2 \mathrm{div}(u)\{\frac{1}{3}[3(\epsilon^r_{\; r})^2 + 3(\epsilon^\theta_{\; \theta})^2 - (\epsilon^r_{\; r} + \epsilon^\theta_{\; \theta})^2]\}$
  
  $\epsilon_q = - l^2 \mathrm{div}(u)\{\frac{1}{3}[3(\epsilon^r_{\; r})^2 + 3(\epsilon^\theta_{\; \theta})^2 - \mathrm{div}(u)^2]\}$
  
  $\epsilon_q = - l^2 \mathrm{div}(u)\{(\epsilon^r_{\; r})^2 + (\epsilon^\theta_{\; \theta})^2 - \frac{1}{3}\mathrm{div}(u)^2\}$
