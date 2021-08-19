Nash LCP Example
*****************

Assume we have two players:

**Player 1**

.. math::

 \min_{q_1}: 10 q_1 + 0.1 q_1^2 - (100 - (q_1+q_2)) q_1  =  1.1 q_1^2 - 90q_1 + q_1q_2

 \text{s.t:} q_1 >= 0


**Player 2**

.. math::

 \min_{q_2}: 5 q_2 + 0.2 q_2^2 - (100 - (q_1+q_2)) q_2 	= 1.2 q_2^2 - 95 q_2 + q_2q_1

 \text{s.t:} q_2 >= 0

The above problem corresponds to a `Cournot Competition <https://en.wikipedia.org/wiki/Cournot_competition>`_
where the demand curve is given by :math:`P = a-BQ` where ``P`` is the market price and ``Q`` is the quantity in the market. The cost of production of both the producers are given by a convex quadratic function in the quantity they produce. The solution to the problem is to find a Nash Equilibrium from which neither producer is incentivized to deviate.

====================================
Modeling the problem
====================================

To handle this problem, first we create two objects of :cpp:class:`MathOpt::QP_Param` to model each player's optimization problem, as parameterized by the other.

.. code-block:: c

        arma::sp_mat Q(1, 1), A(0, 1), B(0, 1), C(1, 1);
        arma::vec b, c(1);
        b.set_size(0);

        Q(0, 0) = 2 * 1.1;
        C(0, 0) = 1;
        c(0) = -90;
        auto q1 = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, &env);

        Q(0, 0) = 2 * 1.2;
        c(0) = -95;
        auto q2 = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, &env);

        std::vector<shared_ptr<MathOpt::QP_Param>> q{q1, q2}; // Making a vector shared_ptr to the individual players' problem


Next, since we do not have any Market clearing constraints, we set empty matrices for them. Note that, if the problem does not have market clearing constraints, still the matrices have to be input with zero rows and appropriate number of columns.

.. code-block:: c

        sp_mat MC(0, 2);
        vec MCRHS;
        MCRHS.set_size(0);

Finally now, we can make the :cpp:class:`Game::NashGame` object by invoking the constructor.

.. code-block:: c

 	GRBEnv env;
        Game::NashGame Nash = Game::NashGame(&env, q, MC, MCRHS);


Using traditional means, one can write a linear complementarity problem (LCP)
to solve the above problem. The LCP is given as follows.

.. math::

 0 \le q_1 \perp 2.2 q_1 + q_2 - 90 \geq 0

 0 \le q_2 \perp q_1 + 2.4 q_2 - 95 \geq 0

To observe the LCP formulation of this NashGame, one can use :cpp:func:`Game::NashGame::FormulateLCP` member function.

.. code-block:: c

 	arma::sp_mat M;
 	arma::vec q;
 	perps Compl;		// Stores the complementarity pairs relationships.
 	Nash.FormulateLCP(M, q, Compl);	// Compute the LCP
        M.print();
        q.print(); 

Here ``M`` and ``q`` are such that the solution to the LCP :math:`0 \le x \perp Mx + q \ge 0` solves the original NashGame. These matrices can be written to a file and solved externally now.
Alternatively, one can pass it to the :cpp:class:`Game::LCP` class, and solve it natively. To achieve this, one can pass the above matrices to the constructor of the :cpp:class:`Game::LCP` class.

.. code-block:: c

        GRBEnv env = GRBEnv();
        Game::LCP lcp = Game::LCP(&env, M, q, 1, 0);

More concisely, the class :cpp:class:`Game::LCP` offers a constructor with the NashGame itself as an argument. This way, one need not explicitly compute ``M``, ``q`` etc., to create the ``Game::LCP`` object.

.. code-block:: c

        Game::LCP lcp2 = Game::LCP(&env, Nash);


Now the ``Game::LCP`` object can be solved. And indeed the solution helps obtain the Nash equilibrium of the original Nash game.

.. code-block:: c

 auto model = lcp.LCPasMIP();
 model.optimize();
 // Alternatively, auto model = lcp.LCPasMIP(true); will already optimize and solve the model.

As was the case with :cpp:func:`MathOpt::QP_Param::solveFixed`, the above function returns a
``unique_ptr`` to ``GRBModel``. And all native operations to the ``GRBModel` can be performed and the solution be obtained.

====================================
Checking the solution
====================================

The solution to this problem can be obtained as :math:`q_1=28.271028, q_2=27.803728`. To indeed check that this solution is correct, one can create a solution vector and solve each player's :cpp:class:`MathOpt::QP_Param` and check that the solution indeed matches.
 
.. code-block:: c

  arma::vec Nashsol(2);
  Nashsol(0) = model->getVarByName("x_0").get(GRB_DoubleAttr_X); // This is 28.271028 
  Nashsol(1) = model->getVarByName("x_1").get(GRB_DoubleAttr_X); // This is 27.803728

  auto nashResp1 = Nash.respond(0, Nashsol);
  auto nashResp2 = Nash.respond(1, Nashsol);

  cout<<nashResp1->getVarByName("y_0").get(GRB_DoubleAttr_X)<<endl; // Should print 28.271028
  cout<<nashResp2->getVarByName("y_0").get(GRB_DoubleAttr_X)<<endl; // Should print 27.803728


One can, thus check that the values match the solution values obtained earlier. If only does not want the individual ``GRBModel`` handles, but just want to confirm either that the problem is solved or to provide a player with profitable deviation, one can just use :cpp:func:`Game::NashGame::isSolved` function as follow.

.. code-block:: c

	unsigned int temp1 ; arma::vec temp2;
	cout<<Nash.isSolved(Nashsol, temp1, temp2); // This should be true.


If the :cpp:func:`Game::NashGame::isSolved` function returns false, then ``temp1`` and ``temp2`` respectively contain the player with profitable deviation, and the more profitable strategy of the player.

And note that, just like :cpp:class:`MathOpt::QP_Param`, :cpp:class:`Game::NashGame` can also be saved to and loaded from an external file.

.. code-block:: c

        Nash.save("dat/Nash.dat"); //Saves the object
        Game::NashGame Nash2(&env);
        Nash2.load("dat/Nash.dat"); // Loads the object into memory.