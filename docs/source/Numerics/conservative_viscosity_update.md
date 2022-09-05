# Conservative form

  Baruteau 2008 writes the viscosity updates as:

## Radial velocity
  
  $\frac{\partial \Sigma v_r}{\partial \mathrm{d}t} = \frac{1}{r}[\frac{\partial (r \tau_{rr})}{\partial r} + \frac{\partial (\tau_{r\varphi})}{\partial \varphi} - \tau_{\varphi \varphi}]$

Looking at the $r$ part, the non-conservative form is discretized as

$\frac{\partial \Sigma v_r}{dt} = \frac{1}{r}\frac{\partial (r \tau_{rr})}{\partial r} = \frac{1}{R_a^i}\frac{R_b^i \tau_{rr}^i - R_b^{i-1} \tau_{rr}^{i-1}}{R_b^i - R_b^{i-1}}$

The conservative form can is derived as follows:

$\int \frac{\partial \Sigma v_r}{\partial t} \mathrm{d}V = \int \frac{1}{r}\frac{\partial r \tau_{rr}}{\partial r} \mathrm{d}V$

with $\mathrm{d}V = r \mathrm{d}r \matrm{d} \varphi$, we have: where we
integrate $\varphi$ over $\varphi^{j-1/2}$ to $\varphi^{j+1/2}$
and integrate $r$ from $R_b^{i-1}$ to $R_b^i$ such that the middle of the
integration is at $R_a^i$, $\varphi^j$,  where $v_r$ is located.

$\int \frac{\partial \Sigma v_r}{\partial t} r \mathrm{d}r \matrm{d} \varphi = \int \frac{1}{r}\frac{\partial r \tau_{rr}}{\partial r} r \mathrm{d}r \matrm{d} \varphi$

$\Sigma \frac{1}{2}((R_b^i)^2 - (R_b^{i-1})^2) \Delta \varphi \frac{\partial v_r}{\partial t} = \int^{R_b^i}_{R_b^{i-1}} \partial (r \tau_{rr}) \Delta\varphi$

$\Sigma \frac{1}{2}((R_b^i)^2 - (R_b^{i-1})^2) \Delta \varphi \frac{\partial v_r}{\partial t} = (R_b^i \tau_{rr}^i - R_b^{i-1} \tau_{rr}^{i-1}) \Delta \varphi$

$\frac{\partial v_r}{dt} = \frac{2}{\Sigma}\frac{R_b^i \tau_{rr}^i - R_b^{i-1} \tau_{rr}^{i-1}}{(R_b^i)^2 - (R_b^{i-1})^2}$

$\frac{\partial v_r}{dt} = \frac{1}{\Sigma}\frac{2}{R_b^i + R_b^{i-1}} \frac{R_b^i \tau_{rr}^i - R_b^{i-1} \tau_{rr}^{i-1}}{R_b^i - R_b^{i-1}}$


$\int \frac{\partial \Sigma v_r}{\partial t} \mathrm{d}V = \int \frac{1}{r}\frac{\partial \tau_{r\varphi}}{\partial \varphi} \mathrm{d}V$

$\Sigma \frac{1}{2}((R_b^i)^2 - (R_b^{i-1})^2) \Delta \varphi \frac{\partial v_r}{\partial t} = \int \frac{1}{r}\frac{\partial \tau_{r\varphi}}{\partial \varphi}  r \mathrm{d}r \matrm{d} \varphi$

$\Sigma \frac{1}{2}((R_b^i)^2 - (R_b^{i-1})^2) \Delta \varphi \frac{\partial v_r}{\partial t} = \int \partial^{\varphi} \tau_{r\varphi}  \mathrm{d}r$


$\Sigma \frac{1}{2}((R_b^i)^2 - (R_b^{i-1})^2) \Delta \varphi \frac{\partial v_r}{\partial t} = (R_b^{i} - R_b^{i-1})(\tau_{r\varphi}^{j+1} - \tau_{r\varphi}^j)$

$\Sigma \Delta \varphi \frac{\partial v_r}{\partial t} = \frac{2(R_b^{i} - R_b^{i-1})}{((R_b^i)^2 - (R_b^{i-1})^2)}(\tau_{r\varphi}^{j+1} - \tau_{r\varphi}^j)$

$\frac{\partial v_r}{\partial t} = \frac{1}{\Sigma} \frac{2}{R_b^i + R_b^{i-1}} \frac{\tau_{r\varphi}^{j+1} - \tau_{r\varphi}^j} {\Delta \varphi}$



$\int \frac{\partial \Sigma v_r}{\partial t} \mathrm{d}V = \int \frac{1}{r} \tau_{\varphi\varphi} \mathrm{d}V$

$\Sigma \frac{1}{2}((R_b^i)^2 - (R_b^{i-1})^2) \Delta \varphi \frac{\partial v_r}{\partial t} = \int \frac{1}{r}\tau_{\varphi \varphi}  r \mathrm{d}r \matrm{d} \varphi$

$\Sigma \frac{1}{2}((R_b^i)^2 - (R_b^{i-1})^2) \Delta \varphi \frac{\partial v_r}{\partial t} = \tau_{\varphi \varphi}  (R_b^i - R_b^{i-1}) \Delta \varphi$

$\frac{\partial v_r}{\partial t} = \frac{1}{\Sigma} \frac{2(R_b^i - R_b^{i-1})}{(R_b^i)^2 - (R_b^{i-1})^2} \tau_{\varphi \varphi}$

$\frac{\partial v_r}{\partial t} = \frac{1}{\Sigma} \frac{2}{R_b^i + R_b^{i-1}} \tau_{\varphi \varphi}$

The total $v_r$ update is just the sum of the 3 parts we just covered:


$\frac{\partial v_r}{\partial t} = \frac{2}{\Sigma^i + \Sigma^{i-1}} \frac{2}{R_b^i + R_b^{i-1}} (\frac{R_b^i \tau_{rr}^i - R_b^{i-1} \tau_{rr}^{i-1}}{R_b^i - R_b^{i-1}} + \tau_{\varphi \varphi} + \frac{\tau_{r\varphi}^{j+1} - \tau_{r\varphi}^j} {\Delta \varphi})$


## Azimuthal velocity

   For the azimuthal velocity, it is important to do the update for the angular
   momentum such that it is conserved (see D'Angelo et al. 2002  Eq. 4).
   
The term $\frac{(\partial r \tau_{r\varphi})}{\partial r} + \tau_{r \varphi}$ can be rewritten to $\frac{1}{r}\frac{(\partial r^2 \tau_{r\varphi})}{\partial r}$ and we can write for the update:

Considering both these things, the velocity update in Masset 2002:

  $\frac{\partial \Sigma v_\varphi}{\partial t} = \frac{1}{r}[\frac{\partial (r \tau_{r\varphi})}{\partial r} + \frac{\partial (\tau_{\varphi \varphi})}{\partial \varphi} + \tau_{r \varphi}]$

  can be rewritten to the angular momentum update seen in D'Angelo et al. 2002

  $\frac{\partial \Sigma l}{\partial t} = \frac{1}{r}\frac{\partial (r^2 \tau_{r\varphi})}{\partial r} + \frac{\partial (\tau_{\varphi \varphi})}{\partial \varphi}$
  
Again we integrate over $\mathrm{d}V = r \mathrm{d}r \matrm{d} \varphi$, we
have: we have $\varphi$ ranging over $\varphi^{j+1}$ to $\varphi^{j}$ and
integrate $r$ from $R_a^{i+1}$ to $R_a^i$ such that the middle of the
integration is at $R_b^i$, $\varphi^{j-1/2}$, where $v_\varphi$ is located.


  
  $\int \frac{\partial \Sigma l}{\partial t} \mathrm{d}V = \int \frac{1}{r}\frac{\partial (r^2 \tau_{r\varphi})}{\partial r} \mathrm{d}V$

  $\int \frac{\partial \Sigma l}{\partial t} r \mathrm{d}r \matrm{d} \varphi = \int \partial (r^2 \tau_{r\varphi})\matrm{d} \varphi$

  
  $\frac{1}{2}((R_a^{i+1})^2 - (R_a^i)^2) \Delta \varphi \frac{\partial \Sigma l}{\partial t} = ((R_a^{i+1})^2 \tau_{r\varphi}^{i+1} - (R_a^{i})^2 \tau_{r\varphi}^{i}) \Delta \varphi$

  $\frac{\partial \Sigma l}{\partial t} = \frac{2}{((R_a^{i+1})^2 - (R_a^i)^2)} ((R_a^{i+1})^2 \tau_{r\varphi}^{i+1} - (R_a^{i})^2 \tau_{r\varphi}^{i})$

  $\Sigma R_b^i \frac{\partial v_\varphi}{\partial t} = \frac{2}{((R_a^{i+1})^2 - (R_a^i)^2)} ((R_a^{i+1})^2 \tau_{r\varphi}^{i+1} - (R_a^{i})^2 \tau_{r\varphi}^{i})$

  
  $\frac{\partial v_\varphi}{\partial t} = \frac{1}{\Sigma R_b^i} \frac{2}{((R_a^{i+1})^2 - (R_a^i)^2)} ((R_a^{i+1})^2 \tau_{r\varphi}^{i+1} - (R_a^{i})^2 \tau_{r\varphi}^{i})$


  
  $\int \frac{\partial \Sigma l}{\partial t} \mathrm{d}V = \int \frac{\partial (\tau_{\varphi\varphi})} {\partial \varphi}  \mathrm{d}V$
  
  $\frac{1}{2}((R_a^{i+1})^2 - (R_a^i)^2) \Delta \varphi \frac{\partial \Sigma l}{\partial t} = \int \frac{\partial (\tau_{\varphi\varphi})} {\partial \varphi} r \mathrm{d}r \mathrm{d}\varphi$

  
  $\frac{1}{2}((R_a^{i+1})^2 - (R_a^i)^2) \Delta \varphi \frac{\partial \Sigma l}{\partial t} = \frac{1}{2}((R_a^{i+1})^2 - (R_a^i)^2) (\tau_{\varphi\varphi}^{j} - \tau_{\varphi\varphi}^{j-1})$
  
  $\frac{\partial \Sigma l}{\partial t} = \frac{(\tau_{\varphi\varphi}^{j} - \tau_{\varphi\varphi}^{j-1})}{\Delta \varphi}$

  $\frac{\partial v_\varphi}{\partial t} = \frac{1}{\Sigma R_b^i} \frac{(\tau_{\varphi\varphi}^{j} - \tau_{\varphi\varphi}^{j-1})}{\Delta \varphi}$

  The total $v_\varphi$ update is then just the sum of the 2 components:

  
  $\frac{\partial v_\varphi}{\partial t} = \frac{2}{\Sigma^{j} + \Sigma^{j-1}}\frac{1}{R_b^i} (\frac{2}{((R_a^{i+1})^2 - (R_a^i)^2)} ((R_a^{i+1})^2 \tau_{r\varphi}^{i+1} - (R_a^{i})^2 \tau_{r\varphi}^{i}) + \frac{(\tau_{\varphi\varphi}^{j} - \tau_{\varphi\varphi}^{j-1})}{\Delta \varphi})$
