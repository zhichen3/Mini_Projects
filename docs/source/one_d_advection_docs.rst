********************
1-D advection in C++
********************

See :ref:`1-D advection` under API for the implementation of solving both linear and nonlinear 1-D advection method in C++.

Descriptions
==============

_1DGrid
--------
``_1DGrid.H`` contains the class for creating a 1D Grid for advection to take place. It creates a grid on the physical domain [xmin, xmax] with finite resolution of grid points nx and ghost points ng. 

solver
-------
``solver.H`` contains the 1D advection solver with linear methods like *ftcs* and *upwinding* methods, and nonlinear methods like *predictor corrector* and *method of lines* schemes. Some slope are also provides such as *piecewise constant*, *centered*, *minmod*, and *MC*.

Initial Condition
------------------
``init_cond.H`` contains three example initial conditions, including *tophat*, *sine*, and *gaussian*.

Initial condition function should take the physical domain, e.g. :math:`x` as input parameter.

Verification
=============

To verify the validity of the advection solver, states shoudl return back to the original state after full cycles. Therefore, enter ``num_periods`` as a whole number, and see the evolved state comparing with the initial state. Initial state and final state should appear in the same directory for plotting purposes.

* If ``num_periods`` is a whole number ``solver::solve()`` also prints out the relative error between the final and the initial state.

* More accurate comparison is to plot the final and initial state. There is a pre-written python script that would plot all the dat files. Simply run the following command, where [dat files] are the dat files after running the advection solver. One could use ``*.dat`` to plot all the files.

    .. prompt:: bash

       ./plot [dat files]

* One can also do a runtime visualization of the data. In order to create runtime data after each time step, set ``runtime_dat = true`` for ``advection_solver``. runtime data are stored inside folder ``runtime_data``. See the visual by the following command:

  .. prompt:: bash
     ./plot runtime_data/*.dat --animate
       
There is a pre-wrritten ``unit_test.cpp`` that tests this mechanism and how to use the code. Simply run this file and plot all the output files. Or the following command

.. prompt:: bash

   make unit_test DEBUG=TRUE
   ./unit_test
   ./plot *dat

There is also a pre-written ``runtime_test.cpp`` used for testing runtime visualization. To run the ``runtime_test.cpp`` for test, use the following command:

.. prompt:: bash

   make runtime_test.cpp DEBUG=TRUE
   ./runtime_test
   ./plot runtime_data/* --animate

  
