LCPs
***************

In the previous tutorial, we introduced the LCP problem. Here we give a brief overview of :cpp:class:`MathOpt::LCP`.
This class provides the basic tools to solve LCPs. Since the solution of LCPs are union of polyhedra, the inheritor class :cpp:class:`MathOpt::PolyLCP` provides support for the polyhedral aspect of LCPs.

The types of problems that are handled by the class :cpp:func:`MathOpt::LCP` are in the following form:

.. math::

  Ax + By \leq b

  0 \leq x \perp Mx + Ny + q \geq 0


Yet we use a different notation. Instead of using ``y`` to refer to the variables that don't have matching complementary equations, we call *all* the variables as ``x`` and we keep track of the position of variables which are not complementary to any equation.

.. note::

    - The set of indices of ``x`` which are not complementary to any equation should be a consecutive set of indices. For the sake of clarity, these components will be called as *Leader vars components* of ``x``.
    - Suppose the leader vars components of ``x`` are removed from ``x``, in the remaining components, the first component should be complementary to the first row defined by @p M, second component should be complementary to the second row defined by ``M`` and so on.



====================================
Modeling the problem
====================================

Now consider the following linear complementarity problem.

.. math::
        x_1 + x_2 + x_3 \le 12

        0\le x_1 \perp x_4 - 1 \ge 0

        0\le x_2 \le 2 

        0 \le x_3 \perp 2x_3 + x_5 \ge 0

        0 \le x_4 \perp -x_1 + x_2 + 10 \ge 0

        0 \le x_5 \perp x_2 - x_3 + 5 \ge 0


The variable :math:`x_2` has no complementarity equation. This problem can be entered into the :cpp:class:`MathOpt::LCP` class as follows.

.. code-block:: c

    // We have four complementarity equations and 5 variables.
    arma::vec q(4); M.zeros();
    arma::sp_mat M(4, 5);
    // First equation
    M(0, 3) = 1;
    q(0) = -1;
    // Second equation
    M(1, 2) = 2;
    M(1, 4)  = 1;
    q(1) = 0;
    // Third equation
    M(2, 0) = -1;
    M(2, 1) = 1;
    q(2) = 10;
    // Fourth equation
    M(3, 1) = 1 ;
    M(3, 2) = -1;
    q(3) = 5;
    // Other common constraints
    arma::sp_mat A(2, 5); arma::vec b;
    A.zeros();
    // x_2 <= 2 constraint
    A(0, 1) = 1;
    b(0) = 2;
    // x_1 + x_2 + x_3 <= 12 constraint
    A(1, 0) = 1;
    A(1, 1) = 1;
    A(1, 2) = 1;
    b(1) = 12;

Since the variable with no complementarity pair is :math:`x_2` which is in position ``1`` (counting from 0) of the vector ``x``, the arguments ``LeadStart`` and ``LeadEnd`` in the constructor, :cpp:func:`MathOpt::LCP::LCP` are ``1`` as below.

.. code-block:: c
   
   GRBEnv env;
   LCP lcp = LCP(&env, M, q, 1, 1, A, b);

====================================
Computing solutions
====================================

This problem can be solved either with a MIP, a MINLP, or with PATH (:cpp:enum:`Data::LCP::Algorithms`). You refer :cpp:func:`MathOpt::LCP::solve` for various solution method.

.. code-block:: c

 // Solve using PATH
 arma::vec x;
 arma::vec z;
 auto indModel = lcp.solve(Data::LCP::Algorithms::PATH,x,z,-1,1);


This LCP as multiple solutions. In fact the solution set can be parameterized as below.

.. math::

 x_1 &= 10 + t
 
 x_2 &= t
 
 x_3 &= 0
 
 x_4 &= 1
 
 x_5 &= 0 
 
 \text{for}\;\; t \in [0, 1]
 
====================================
Utilities
====================================

Two functions :cpp:func:`MathOpt::LCP::LCPasMILP` and :cpp:func:`MathOpt::LCP::LCPasMIQP` allows to to optimize a linear objective function or a convex quadratic
objective function over the set of solutions. Also, note that :cpp:func:`MathOpt::LCP::setMIPLinearObjective`, :cpp:func:`MathOpt::LCP::setMIPQuadraticObjective`, :cpp:func:`MathOpt::LCP::setMIPFeasibilityObjective` can change the objective function of the MIP model (if one is called for solving the LCP).
In general, we recomend to use :cpp:func:`MathOpt::LCP::solve`, which is a general methods that delegates the solution to either one available solver.

