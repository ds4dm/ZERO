Parametrized IPs
*****************


IP_Param stands for *Parametrized Integer Program*, a mathematical program in the following form:

.. math::
    \min_y c^Ty + (Cx)^T y + d^Tx \\
    \text{s.t.} \quad By \le b \\
    \quad \quad y \ge 0 \\
    \quad \quad y_i \in  \mathbb{Z} \quad \forall i \in I


Where :math:`y` are the decision variables for the program, :math:`x` are parameters, and :math:`I` is the set of indices of integer variables.
ZERO supports Integer Programs with linear constraints and bi-linear objectives (linear in :math:`y`).
You can find the API information in :cpp:class:`MathOpt::IP_Param`. This class is an inheritor of :cpp:class:`MathOpt::MP_Param`, and we use the same notation of the abstract.


====================================
Modeling the problem
====================================

Following the previous example on parametrized Quadratic Programs, consider the following parametrized Integer Program:

.. math::

 \max_{y_1, y_2} y_1 + 2y_2 - 2y_1x_1 -3y_2x_2

 \text{s.t.} \quad   3y_1+4y_2 &\le 5

 \quad \quad  y_1, y_2 &\in \{0,1\}

This is a Knapsack Problem in :math:`y`, with an objective function parametrized in the  :math:`x` variables.
We assume the matrix C is a matrix with:

* A number of rows corresponding to the number of :math:`y` variables
* A number of columns corresponding to the number of :math:`x` parameters

.. code-block:: c

    #include "zero.h"
    //[...] Your other code here
    //The linear objectives in y
    arma::vec                _c(2), _d(2);
    //The matrix of interaction coefficients among parameters and variables
    arma::mat                _C = arma::zeros(2, 2);
    //The constraint matrix
    arma::mat                _A = arma::zeros(1, 2);

    //The list of indexes of variable that are integer
    arma::vec _integers(_c.size());

    //All variables are integer here
    for (unsigned int i = 0; i < _c.size(); ++i)
        _integers.at(i) = i;


    //Invert the sign since we have a maximization
    _C(0, 0) = 2;
    _C(1, 1) = 3;
    _c(0) = -1;
    _c(1) = -2;
    _d(0) = 0;
    _d(1) = 0;

    //Knapsack Constraint
    _A(0, 0) = 3;
    _A(0, 1) = 4;

    arma::vec _b = arma::vec{5};

    //Convert to sparse objects
    arma::sp_mat      _C2  = arma::sp_mat{_C};
    arma::sp_mat      _A2  = arma::sp_mat{_A};

    //Explicitly add the bounds without having to put additional constraints.
    VariableBounds    bnds = {{0, 1}, {0, 1}};
    //Gurobi Environment
    GRBEnv      test;
    MathOpt::IP_Param ipParam(_C2, _A2, _b, _c, _d, _integers, bnds, &test);


====================================
Computing solutions
====================================
With :math:`(x_1, x_2) = (-1, 0.5)`, we can provide the value of :math:`f` and solve it through Gurobi

.. code-block:: c

    // Enter the value of x in an arma::vec
    arma::vec x(2);
    x(0) = -1;
    x(1) = 0.5;

    // Uses Gurobi to solve the model, returns a unique_ptr to GRBModel
    auto FixedModel = ipParam.solveFixed(x,true);
