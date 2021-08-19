IP_Param Example
*****************
Following the previous example on quadratic programs, consider the following Integer Program.

.. math::

 \max_{y_1, y_2} y_1 + 2y_2 - 2y_1x_1 -3y_2x_2

 \text{s.t.}\;\;\;\;\;  3y_1+4y_2 &\le 5
 
 \;\;\;\;\;\;\;\; y_1, y_2 &\in \{0,1\}

This is a Knapsack Problem in :math:`y` with an objective function parametrized in some :math:`x` variables.

====================================
Modeling the problem
====================================

This data can be entered as follows. We assume the matrix Q is a matrix with:
- A number of rows corresponding to the number of :math:`y` variables
- A number of columns corresponding to the number of :math:`x` parameters

.. code-block:: c

     //The linear objectives in y
	 arma::vec                _c(2);
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
	 MathOpt::IP_Param ipParam(_C2, _A2, _b, _c, _integers, bnds, &test);

With :math:`(x_1, x_2) = (-1, 0.5)`, we can provide the value of :math:`f` and solve it through Gurobi

.. code-block:: c

        vec x(2);			// Enter the value of x in an arma::vec
        x(0) = -1;
        x(1) = 0.5;

        auto FixedModel = ipParam.solveFixed(x);	// Uses Gurobi to solve the model, returns a unique_ptr to GRBModel


Analogously to `QP_Param`, this `IP_Param` class inherit most of the methods of the `MP_Param` abstract class.