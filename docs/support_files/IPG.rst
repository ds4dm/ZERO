IPG Example
***************
Consider the following Integer Programming Game: The first player is the **x** player, where its
leader's decision variables are :math:`x`. The second player is the **y** player where its variables are :math:`y`.

.. image:: IPG.png
  :width: 500
  :alt: The Integer Programming Game for this example is a Knapsack Game.

====================================
The Players
====================================
The x player's optimization problem is given below

.. math::

 \max_{x_1, x_2} x_1 + 2x_2 - 2x_1y_1 -3x_2y_2

 \text{s.t.}\;\;\;\;\;  3x_1+4x_2 &\le 5

 \;\;\;\;\;\;\;\; x_1, x_2 &\in \{0,1\}


While the y player's optimization problem is:

.. math::

 \max_{y_1, y_2} 3y_1 + 5y_2 - 5y_1x_1 -4y_2x_2

 \text{s.t.}\;\;\;\;\;  2y_1+5y_2 &\le 5

 \;\;\;\;\;\;\;\; y_1, y_2 &\in \{0,1\}


====================================
Equilibria
====================================
The problem has two pure Nash equilibria
:math:`(x_1, x_2, y_1, y_2) = (0, 1, 1, 0)`, and :math:`(x_1, x_2, y_1, y_2) = (1, 0, 0, 1)`, and a mixed Equilibrium :math:`(x_1, x_2, y_1, y_2) = (0.22, 0.78, 0.4, 0.6)`.

====================================
Modeling the problem
====================================

The first step in modeling this Integer Programming Game is to include `zero.h` and create a derived class of :cpp:class:`Game::IPG`. The minimal constructor for :cpp:class:`Game::IPG` involves passing a pointer to `GRBEnv` (Check Gurobi's C++ `reference manual <https://www.gurobi.com/documentation/8.1/refman/cpp_api_overview.html>`_
). The derived class should indeed instantiate the base class (Game::IPG) using such a constructor. The code below gives an example.


.. code-block:: c

 #include "zero.h"
 GRBEnv      test;
 Models::IPG::IPGInstance Instance;
 try {
     //First player
	 arma::vec                _c(2);
	 arma::mat                _C = arma::zeros(2, 2);
	 arma::mat                _A = arma::zeros(1, 2);

	 arma::vec _integers(_c.size());

	 for (unsigned int i = 0; i < _c.size(); ++i)
		_integers.at(i) = i;

	 _C(0, 0) = 2;
	 _C(1, 1) = 3;
	 _A(0, 0) = 3;
	 _A(0, 1) = 4;
	 _c(0) = -1;
	 _c(1) = -2;
	 arma::vec _b = arma::vec{5};
	 arma::sp_mat      _C2  = arma::sp_mat{_C};
	 arma::sp_mat      _A2  = arma::sp_mat{_A};
	 VariableBounds    bnds = {{0, 1}, {0, 1}};
	 MathOpt::IP_Param ipParam(_C2, _A2, _b, _c, _integers, bnds, &test);
	 Instance.addIPParam(ipParam, BasePath + "raw/CustomJulia-p1");

     //Second player
	 arma::sp_mat _Ct(2, 2);
	 arma::mat    _At = arma::zeros(1, 2);

	 _Ct(0, 0) = 5;
	 _Ct(1, 1) = 4;

	 _At(0, 0) = 2;
	 _At(0, 1) = 5;
	 //_At(0, 2) = 6;

	 _c(0) = -3;
	 _c(1) = -5;
	 //_c(2) = -1;

	 arma::sp_mat      _C22 = arma::sp_mat{_Ct};
	 arma::sp_mat      _A22 = arma::sp_mat{_At};
	 MathOpt::IP_Param ipParam2(_C22, _A22, _b, _c, _integers, bnds, &test);
	 Instance.addIPParam(ipParam2, BasePath + "raw/CustomJulia-p2");

	 Models::IPG::IPG Test(&test, Instance);
	 Test.setAlgorithm(Data::IPG::Algorithms::Oracle);
	 Test.setDeviationTolerance(3e-4);
	 Test.setNumThreads(4);
	 Test.setLCPAlgorithm(Data::LCP::Algorithms::MIP);
	 Test.setGameObjective(Data::IPG::Objectives::Quadratic);
	 // Test.setTimeLimit(600);
	 Test.finalize();
	 Test.findNashEq();


- With the method `setAlgorithm` of :cpp:class:`Game::IPG`, we set the algorithm that will solve the Integer Programming Game. So far, only :cpp:class:`Algorithms::IPG::Oracle` is available.
- The method `setLCPAlgorithm` specifies the algorithm used to solve the LCPs. It can be either :cpp:class:`Data::LCP::Algorithms::MIP`, :cpp:class:`Data::LCP::Algorithms::PATH`, or :cpp:class:`Data::LCP::Algorithms::MINLP`.
- The game's objective (not supported by PATH) forces an objective into the LCP problem as to increase the chances of finding a good equilibrium given the objective. Values can be :cpp:class:`Data::IPG::Objectives::Quadratic` :cpp:class:`Data::IPG::Objectives::Linear` :cpp:class:`Data::IPG::Objectives::Feasibility`.
- Other options can be found in the documentation of :cpp:class:`Game::IPG`