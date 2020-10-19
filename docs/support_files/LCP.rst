LCP Example
***************

 Before reading this page, please ensure you are aware of the functionalities
 described in @link NashGame_Example Game::NashGame tutorial @endlink before
 following this page.

 Consider the Following linear complementarity problem with constraints

.. math::

  Ax + By \leq b

  0 \leq x \perp Mx + Ny + q \geq 0`


These are the types of problems that are handled by the class :py:func:`Game::LCP`, but we use a different notation. Instead of using ``y`` to refer to the variables that don't have matching complementary equations, we call *all* the variables as ``x`` and we keep track of the position of variables which are not complementary to any equation.

**Points to note:**

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


Here indeed :math:`x_2` is the leader vars component with no complementarity equation. This problem can be entered into the :py:class:`Game::LCP` class as follows.

.. code-block:: c

        arma::sp_mat M(4, 5); // We have four complementarity eqns and 5 variables. arma::vec q(4); M.zeros();
        // First eqn
        M(0, 3) = 1;
        q(0) = -1;
        // Second eqn
        M(1, 2) = 2;
        M(1, 4)  = 1;
        q(1) = 0;
        // Third eqn
        M(2, 0) = -1;
        M(2, 1) = 1;
        q(2) = 10;
        // Fourth eqn
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

Now, since the variable with no complementarity pair is :math:`x_2` which is in position ``1`` (counting from 0) of the vector ``x``, the arguments ``LeadStart`` and ``LeadEnd`` in the constructor, :py:func:`Game::LCP::LCP` are ``1`` as below.

.. code-block:: c
   
   GRBEnv env;
   LCP lcp = LCP(&env, M, q, 1, 1, A, b);

==========
Solving
==========

This problem can be solved either using big-M based disjunctive formulation with the value of the @p bigM can also be chosen.
But a more preferred means of solving is by using SOS1 constraints, or indicator constraints. Use the former option, only if you are very confident of  your choice of a small value of ``bigM``

.. code-block:: c

 // Solve using bigM constraints
 lcp.useIndicators = false;
 lcp.bigM = 1e5;
 auto bigMModel = lcp.LCPasMIP(true);

 // Solve using indicator constraints
 lcp.useIndicators = true;
 auto indModel = lcp.LCPasMIP(true);


Both ``bigMModel`` and ``indModel`` are :py:func:`std::unique_ptr`  to ``GRBModel`` objects. So all native gurobi operations can be performed on these objects.
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

However, sometimes one might want to solve an MPEC. i.e., optimize over the feasible region of the set as decribed above. For this purpose, two functions :py:func:`Game::LCP::MPECasMILP` and :py:func:`Game::LCP::MPECasMIQP` are available, depending upon whether one wants to optimize a linear objective function or a convex quadratic
objective function over the set of solutions.


