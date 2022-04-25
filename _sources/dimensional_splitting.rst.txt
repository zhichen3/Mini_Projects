**********************
Dimensional Splitting
**********************

Strang Splitting
=================

We alternate the order of the dimensiona updates in each timestep.

Use :math:`\star` to represent an intermediate state. Then the process of updating a timestep consists of two parts:

.. math::

   a^\star_{i,j} = a^n_{i,j} - \Delta t \frac{ua^{n+1/2}_{i+1/2,j}-ua^{n+1/2}_{i-1/2,j}}{\Delta x}

   a^{n+1}_{i,j} = a^\star_{i,j} - \Delta t \frac{va^{\star,n+1/2}_{i,j+1/2} - va^{\star, n+1/2}_{i,j-1/2}}{\Delta y}

Here we do update in the x-direction first and stored to the intermediate state :math:`a^\star_{i,j}`. Then the state is updated again using information in the y-direction with state :math:`a^\star_{i,j}`.

.. warning::

   One must alternate the direction for each timestep, i.e. x-y then y-x.

.. Note::

   Here we doing two different updates independently, i.e. solving two 1-D advection equations:

   .. math::

      a_t + ua_x = 0

      a_t + va_y = 0


Riemann Problem
================

For :math:`u > 0`, i.e. traveling from left to right:

.. math::

   a^{n+1/2}_{i+1/2,L}  = a^n_i + \left. \frac{\Delta x}{2}\left(1-\frac{\Delta t}{\Delta x}u\right) \frac{\partial a}{\partial x}\right|_i + ...

If :math:`u < 0`, i.e. traveling from right to left:

.. math::
   a^{n+1/2}_{i+1/2,R} = a^n_{i+1} - \left. \frac{\Delta x}{2} \left( 1 + \frac{\Delta t}{\Delta x}u\right) \frac{\partial a}{\partial x}\right|_{i+1}

The same appplies to y-direction.

Timestep
===========

Define :math:`\Delta t` for multi-dimensional advection:

.. math::

   \Delta t = C \left( \sum^d_{i=1} \frac{|U \cdot e_d|}{\Delta x_d} \right)^{-1}

where U is :math:`|U \cdot e_d|` is the velocity of magnitude of a direction :math:`d`. 
