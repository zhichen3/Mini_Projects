************************
Nonlinear 1-D advection
************************

Finite-volume Discretization
=============================

Rewrite advection equation in conservation form:

.. math::

   a_t + \nabla[f(a)] = 0 

:math:`f(a) = ua` is the flux of the quantity a. In conservation form, the time derivative of a quantity is related to the divergence of flux.
       
For finite-volume discretization:

.. math::

   \frac{\partial a}{\partial t} = -\frac{\partial F(a)}{\partial x} 

.. math::
   
   \frac{1}{\Delta x} \int^{x_{i+1/2}}_{x_{i-1/2}} a_t dx = -\frac{1}{\Delta x}\int^{x_{i+1/2}}_{x_{i-1/2}} \frac{\partial}{\partial x} f(a) dx

.. math::
   
   \frac{\partial}{\partial t} \langle a\rangle_i = -\frac{1}{\Delta x} \left\{[f(a)]_{i+1/2} - [f(a)]_{i-1/2} \right\}

We have two approaches: 

* Predictor-corrector method or characteristic tracing method: Discretize :math:`\partial a_i / \partial t` as :math:`(a_i^{n+1} -a^n_i)/\Delta t`. Then acheive second-order accuracy by evaluating the righthand of Equation above at midpoint in time :math:`(n+1/2)`.

* Method of lines: recognize PDE becomes ODE, then use standard ODE methods like Runge-Kutta to integrate the system in time.


Predictor-Corrector Scheme:
============================

Take :math:`\partial a_i/\partial t` as :math:`(a^{n+1}_i - a^n_i)/\Delta t` (Euler's method). Then we have equation:

.. math::

   \frac{a^{n+1}_i - a^n_i}{\Delta t} = -\frac{[f(a)]^{n+1/2}_{i+1/2} - [f(a)]^{n+1/2}_{i-1/2}}{\Delta x}
   
To evaluate the fluxes at half-time:

.. math::

   [f(a)]^{n+1/2}_{i+1/2} = f(a^{n+1/2}_{i+1/2})
   
Construct a second-order accurate approximation to :math:`a^{n+1/2}_{i+1/2}` by Taylor expanding. 
Notice there are two interface states that can be constructed, either using data from left, e.g. L, or right, e.g. R.

Left: 

.. math::

   a^{n+1/2}_{i+1/2,L} &= a^n_i + \left .\frac{\Delta x}{2}\frac{\partial a}{\partial x}\right|_i + \left. \frac{\Delta t}{2} \frac{\partial a}{\partial t}\right|_i + \mathcal{O}(\Delta x^2) + \mathcal{O}(\Delta t^2) 

   &= a^n_i + \left. \frac{\Delta x}{2}\frac{\partial a}{\partial x} \right|_i + \left. \frac{\Delta t}{2}\left(-u\frac{\partial a}{\partial x} \right|_i\right)+...

   &= a^n_i + \left. \frac{\Delta x}{2}\left(1-\frac{\Delta t}{\Delta x}u\right) \frac{\partial a}{\partial x}\right|_i + ... 

Right:

.. math::

   a^{n+1/2}_{i+1/2,R} = a^n_{i+1} - \left. \frac{\Delta x}{2} \left( 1 + \frac{\Delta t}{\Delta x}u\right) \frac{\partial a}{\partial x}\right|_{i+1}

Method of Lines
================

Treat it as ODE, and by setting :math:`f(a) = ua`, and solve using Runge-Kutta:

.. math::

   \frac{\partial a_i}{\partial t} = -\frac{ua_{i+1/2} - ua_{i-1/2}}{\Delta x} = g(a_i) = \dot{a_i}

And we have the same Riemann probalem and define:

.. math::

   a_{i+1/2,L} &= a(x_i-\Delta x/2) = a_i + \frac{1}{2}\frac{\partial a_i}{\partial x} \Delta x

   a_{i+1/2,R} &= a(x_i + \Delta x/2) = a_{i+1} - \frac{1}{2}\frac{\partial a_{i+1}}{\partial x}\Delta x

And we use 4th order Runge Kutta to solve this ODE:

.. math::

   \frac{dy}{dt} = f(t,y)  

.. math::
   
   y_{n+1} = y_n + \frac{1}{6}dt(k_1+2k_2+2k_3+k_4) 

.. math::

   t_{n+1} = t_n + dt

with:

.. math::

   k_1 = f(t_n,y_n) 

.. math::

   k_2 = f(t_n+\frac{dt}{2}, y_n + dt\frac{k_1}{2})

.. math::
   
   k_3 = f(t_n+\frac{dt}{2}, y_n + dt\frac{k_2}{2}) 

.. math::
   
   k_4 = f(t_n + dt, y_n + dtk_3)

.. tip::

   Note at least 2nd order Runge Kutta method must be used for this method.


Slope choices
===============

From the two methods above, we left out the choices of the slope :math:`\frac{\partial a_i}{\partial x}`.
We could set :math:`\frac{\partial a_i}{\partial x} = 0`, which implies a *picewise constant* slope. This method does not account for overshoot or undershoots at discontinuity states.

.. tip::

   If one uses piecewise constant slope for predictor corrector method, it gets back the same upwinding method.

Using slope limiters will modify the piecewise linear slopes near extrema to prevent overshoots and undershoots. Generally check whether slope of the left and slope on the right have the same sign, by invoking :math:`a \cdot b > 0`. If they are not the same sign, we are an extrema. Therefore limit the slope to force it be 0.

Mindmod Limiter
----------------

.. math::

   \left. \frac{\partial a}{\partial x} \right|_i = \textrm{minmod}\left( \frac{a_i - a_{i-1}}{\Delta x}, \frac{a_{i+1}-a_i}{\Delta x}\right)

where minmod(a,b) is defined as:

.. math::

   \textrm{minmod}(a,b) = 
   \begin{array}{lll}
   a & if & |a| < |b| & \textrm{and} & a \cdot b > 0 \\
   b & if & |b| < |a| & \textrm{and} & a \cdot b > 0 \\
   0 && \textrm{otherwise} \\
   \end{array}

Monotonized central difference Limiter (MC limiter)
-----------------------------------------------------

Define:

.. math::

   \xi = (a_{i+1} - a_i) \cdot (a_i - a_{i-1})

Then:

.. math::

   \left. \frac{\partial a}{\partial x}\right|_i = 
   \begin{array}{ll}
   \textrm{min} \left[\frac{\left| a_{i+1}-a_{i-1} \right|}{2\Delta x} , 2 \frac{\left| a_{i+1} - a_i \right|}{\Delta x}, 2 \frac{\left| a_i - a_{i-1} \right|}{\Delta x} \right] \textrm{sign} (a_{i+1}-a_{i-1}) & \xi > 0 \\
   0 & \textrm{otherwise}
   \end{array}
