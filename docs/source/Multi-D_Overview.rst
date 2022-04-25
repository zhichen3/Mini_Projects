****************************
Overview
****************************

Theory
=======
Two Dimensional advection equation:

.. math::

   a_t + ua_x + va_y = 0

where :math:`u` is the velocity in the x-direction and :math:`v` is the velocity in the y-direction. Let :math:`a_{i,j}` denote the average of :math:`a(x,y,t)` in zone :math:`i,j`. Extend the domain with a perimeter of ghost cells to set the boundary condtions.

To put everything in Conservative form:

.. math::

   a_t + (ua)_x + (va)_y = 0

Where :math:`ua` and :math:`va` are the flux of advection.

Then we define the average of :math:`a` in a zone by integrating over volume:

.. math::

   a_{i,j} = \frac{1}{\Delta x \Delta y} \int^{x_{i+1/2}}_{x_{i-1/2}} \int^{y_{j+1/2}}_{y_{j-1/2}} a(x,y,t) dxdy

Integrating over the conservative form:

.. math::

   \frac{1}{\Delta x \Delta y} \int^{x_{i+1/2}}_{x_{i-1/2}} \int^{y_{j+1/2}}_{y_{j-1/2}} (a_t+(ua)_x+(va)_y) dxdy = 0 

Or using the divergence theorem:

.. math::

   \frac{\partial a_{i,j}}{\partial t} = -\frac{1}{\Delta x \Delta y}\left \{ \int^{y_{j+1/2}}_{y_{j-1/2}} \left[ (ua)_{i+1/2,j} - (ua)_{i-1/2,j} \right] dy + \int^{x_{i+1/2}}_{x_{i-1/2}} \left[(va)_{i,j+1/2} - (va)_{i,j-1/2}   \right]dx    \right \}


Then redefine :math:`\frac{\partial a}{\partial t}` using Euler equation:

.. math::

   \frac{\partial a}{\partial t} = \frac{a^{n+1}_{i,j} - a^{n}_{i,j}}{\partial t}

Define the flux through the interface as the average over the face of that interface and time:

* x-face:

  .. math::
     \langle(ua)_{i+1/2,j}\rangle_{(t)} = \frac{1}{\Delta y \Delta t} \int^{t^{n+1}}_{t^n} \int^{y_{j+1/2}}_{y_{j-1/2}} (ua)_{i+1/2,j} dydt

* y-face:

  .. math::
     \langle (va)_{i,j+1/2} \rangle_{(t)} = \frac{1}{\Delta x \Delta t} \int^{t^{n+1}}_{t^n} \int^{x_{i+1/2}}_{x_{i-1/2}} (va)_{i+1/2,j} dydt


For second-ortder accurate in time, replace the average in time with the flux at the midpoint in time and the average over the face with the flux at the center of the face, i.e:

.. math::

   \langle(ua)_{i+1/2,j}\rangle \approx (ua)^{u+1/2}_{i+1/2,j}

Thus, which is the same as finite volume method but an additional direction term:

There are a few approaches:
1. Dimensional Splitting (finite volume): treat each dimension to be operated independently from the other
2. Unsplit (finite volume): Treat each dimension to be dependent on the other, and is more accurate and less susceptible to grid effects
3. Method of lines: Treat it as ODE and solve using ODE solvers like Runge-Kutta.

   
