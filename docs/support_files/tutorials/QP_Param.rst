Parametrized QPs
*****************

QP_Param stands for *Parametrized Quadratic Program*, a mathematical program in the following form:

.. math::
    \min_y \frac{1}{2}y^TQy + c^Ty + (Cx)^T y \\
    \text{s.t.} \quad  Ax + By \le b\\
    \quad \quad y \ge 0

Where :math:`y` are the decision variables for the program, and :math:`x` are parameters.
You can find the API information in :cpp:class:`MathOpt::QP_Param`. This class is an inheritor of :cpp:class:`MathOpt::MP_Param`.

====================================
Modeling the problem
====================================

Consider the following quadratic program.

.. math::

 \min_{y_1, y_2, y_3} (y_1 + y_2 - 2y_3)^2 + 2 x_1y_1 + 2 x_2y_1 + 3 x_1y_3 + y_1-y_2+y_3

 \text{s.t.}\;\;\;\;\;  y_1, y_2, y_3 &\ge 0

 \;\;\;\;\;\;\;\; y_1 + y_2 + y_3 &\le 10

 \;\;\;\;\;\;\;-y_1 +y_2 -2y_3 &\le -1 + x_1 + x_2


We can model the problem above as follows:

.. code-block:: c

    #include "zero.h"
    //[...] Your other code here
    unsigned int numParams = 2, numVars = 3, numConstr = 2;

    //Create Q (from dense to sparse)
    arma::mat Qd(3, 3);
    Qd << 1 << 1 << -2 << arma::endr
    << 1 << 1 << -2 << arma::endr
    << -2 << -2 << 4 << arma::endr;
    //Convert from dense to sparse
    arma::sp_mat Q = sp_mat(2 * Qd);

    // The matrix for x and y interaction
    arma::sp_mat C(3, 2);
    C.zeros(); C(0, 0) = 2; C(0, 1) = 2; C(2, 0) = 3;
    // The vector for linear terms in y
    arma::vec c(3);
    c << 1 << arma::endr << -1 << arma::endr << 1 << arma::endr;
    // Constraint matrix for x terms
    arma::sp_mat A(2, 2);
    A.zeros(); A(1, 0) = -1; A(1, 1) = -1;
    // Constraint matrix for y terms
    arma::mat Bd(2, 3); Bd << 1 << 1 << 1 << arma::endr << -1 << 1 << -2 << arma::endr;
    arma::sp_mat B = sp_mat(Bd);
    //RHSs
    arma::vec b(2);
    b(0) = 10;
    b(1) = -1;
    // Create a Gurobi environment to handle any solving related calls.
    GRBEnv env = GRBEnv();

Now the required object can be constructed in multiple ways.

.. code-block:: c

    // Method 1: Make a call to the constructor
    MathOpt::QP_Param q1(Q, C, A, B, c, b, &env);

    // Method 2: Using QP_Param::set member function
    MathOpt::QP_Param q2(&env);
    q2.set(Q, C, A, B, c, b);

    // Method 3: Reading from a file. This requires that such an object is saved to a file at first.
    q1.save("dat/q1dat.dat"); // Saving the file so it can be retrieved.
    MathOpt::QP_Param q3(&env);
    q3.load("dat/q1dat.dat");

    // Checking they are the same
    assert(q1==q2);
    assert(q2==q3);

We can now feed some values for the parameters and compute the corresponding optimal solution.
Assume :math:`(x_1, x_2) = (-1, 0.5)`. Then, we the problem becomes a standard QP as:

.. math::

 \min_{y_1, y_2, y_3} (y_1 + y_2 - 2y_3)^2  -y_2 -2y_3

 \text{s.t.}\;\;\;\;  y_1, y_2, y_3 &\ge 0

 \;\;\;\;\;\;\;\;y_1 + y_2 + y_3 &\le 10

 \;\;\;\;\;\;\;\;-y_1 +y_2 -2y_3 &\le -1.5

Correspondingly, we have the following code:

.. code-block:: c

    // Enter the value of x in an arma::vec
    arma::vec x(2);
    x(0) = -1;
    x(1) = 0.5;

    // Uses Gurobi to solve the model, returns a unique_ptr to GRBModel.
    // With the second parameters, we require a model which has already been solved
    auto FixedModel = q2.solveFixed(x,true);

====================
Computing solutions
====================

``FixedModel`` holds the ``GRBModel`` object, and all operations native to ``GRBModel``, like accessing the value of a variable, a dual multiplier, saving the problem to an .lp file or a .mps file. In particular, we can compare the solution with a hand-calculated one as below.

.. code-block:: c

    arma::vec sol(3);
    // Hard-coding the solution as calculated outside
    sol << 0.5417 << arma::endr << 5.9861 << arma::endr << 3.4722;
    for (unsigned int i = 0; i < numVars; i++)
        assert(abs(sol(i)- FixedModel->getVar(i).get(GRB_DoubleAttr_X)) <= 1e-5);
    cout<<FixedModel->get(GRB_DoubleAttr_ObjVal<<endl; // Will print -12.757


In many cases, one might need the "KKT" conditions of a convex quadratic program.
We can then use the function :py:func:`MathOpt::QP_Param::KKT`:

The function returns ``M``, ``N`` and ``q``, where the KKT conditions can be written as :math:`0 \leq y \perp Mx + Ny + q \geq 0`.

.. code-block:: c

    arma::sp_mat M, N; arma::vec q;
    q1.KKT(M, N, q);
    M.print("M");
    N.print("N");
    q.print("q");
