*********************
Linear 1-D Advection
*********************

1-D advection equation, where :math:`a` represents the state of advection and :math:`u` is the velocity of advection.

.. math::

   a_t + ua_x = 0

With initial condition, :math:`a(x,t = 0)` and boundary condition,
The solution to the equation has form:

.. math::

   a(x,t) = a(x-ut)

FTCS Method
=============

Consider first order approximation on the left but second order approximation on the right:

.. math::

   \frac{a^{n+1}_i - a^n_i}{\Delta t} = -u\frac{a^n_{i+1}} - a^n_{i-1}{2\Delta x}

Rewriting it:

.. math::

   a^{n+1}_i = a^n_i - \frac{C}{2} \left( a^n_{i+1} - a^n_{i-1} \right)
  
where C is the Courant Number.

.. math::

   C = \frac{u\Delta t}{\Delta x}

.. tip::

   For FTCS method, smaller C gives better answer.
   FTCS method is also unstable.

Stability Check
----------------

Stability check through Fourier Mode Discretization:

We substitute, use :math:`I` for imaginary number.

.. math::

   a^n_i = A^ne^{-Ii\theta}

Method is stable if:

.. math::

   \left| \frac{A^{n+1}}{A_n} \right| \leq 1

For FTCS:

.. math::

   a^{n+1}_i = a^n_i - \frac{C}{2}(a^n_{i+1} - a^n_{i-1})
   
   A^{n+1}e^{Ii\theta} = A^ne^{Ii\theta} - \frac{C}{2} \left( A^ne^{I(i+1)\theta} - A^ne^{I(i-1)\theta} \right) 

   A^{n+1} = A^n - \frac{C}{2} A^n \left( e^{I\theta}-e^{-I\theta} \right) 
   
   A^{n+1} = A^n\left( 1 - IC\sin{\theta} \right) 

   The magnitude is therefore:

.. math::

   \left| \frac{A^{n+1}}{A^n} \right|^2 = 1+C^2\sin^2{\theta} \geq 1

Since its always :math:`\geq 1`, there is no way to make it stable.



	   
Upwinding Method
==================
Consider first order expansion on both sides instead of second order expansion on the rhs for FTCS method.

.. math::

   (a_x)_i \approx \frac{a_i - a_{i-1}}{\Delta x} + \mathcal{O}

Or

.. math::

   (a_x)_i \approx \frac{a_{i+1} - a_{i}}{\Delta x} + \mathcal{O}

If :math:`u > 0`, suggesting traveling from left to right, we choose upwinding, i.e. :

.. math::

   \frac{a^{n+1}_i - a^n_i }{\Delta t} = -u\frac{a^n_i - a^n_{i-1}}{\Delta x}

Rewriting:

.. math::

   a^{n+1}_i = a^n_i - C(a^n_i - a^n_{i-1})

.. tip::

   Upwinidng method gives exact solution if :math:`C = 1`. Upwinding is also stable if :math:`C \leq 1`.


Ghost Cells/Ghost Points
=========================

Imagine have N points to describe a fluid from 0, 1, 2, ... N. When we do the update for :math:`a[0]`, but we need information from :math:`a[-1]`, which is outside of the domain. We need to extend the domain past the boundary, by creating ghost points: 

Boundary conditon describes what happens to the domain when you go outside of the scope. 

Ghost points are used to implement boundary conditions, if periodic: 

.. math::

   a^n_{hi+1} = a^n_{lo+1}

   a^n_{lo-1} = a^n_{hi-1} 

if outflow:

.. math::

   a^n_{hi+1} = a^n_{hi+1} 

   a^n_{lo-1} = a^n_{lo-1} 
