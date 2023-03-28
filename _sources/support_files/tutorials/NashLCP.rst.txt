Simultaneous Games
**********************

In this class, we use a Linear Complementarity Problem (*LCP*) to solve a simultaneous non-cooperative game among :math:`n` players (aka Nash Games in the optimization community).
Specifically, each player solves Parametrized Quadratic Program -- stored in an instance of :cpp:class:`MathOpt::QP_Param` -- where the parameters are the other players' decisions the variables are the player's decision variables.
We will use members of the class :cpp:class:`Game::NashGame` to model the game. This class will extensively invoke :cpp:class:`MathOpt::LCP` to find the Nash equilibria.

====================================
A theory primer
====================================
Assume we have a simultaneous non-cooperative game among :math:`n` players, so that each of them solves a **convex** Parametrized Quadratic Program such as:

.. math::
    \min_y \frac{1}{2}y^TQy + c^Ty + (Cx)^T y \\
    \text{s.t.} \quad  Ax + By \le b, \quad  y \ge 0

Then, one can compute a Nash equilibrium for this game by considering the KKT conditions of each player and grouping them.
For each player, the KKT conditions are equivalent to the following LCP problem.

.. math::
    q=\begin{bmatrix} c \\ b \end{bmatrix} \quad  M=\begin{bmatrix} Q & B^{\top}  \\ -B & 0 \end{bmatrix} \quad  N=\begin{bmatrix} C  \\ -A \end{bmatrix}\\
    0 \le x \perp Mx + Ny + q \ge 0

As a standard game-theory notation, let the operator :math:`(\cdot)^i` refer to an object of player :math:`i` with :math:`i \in \{ 1,2,\dots,n\}`, and :math:`(\cdot)^{-i}` be :math:`(\cdot^1,\dots, \cdot^{i-1},\cdot^{i+1},\dots,\cdot^{n})`.
The :math:`y` variables for player :math:`i` are then just the vector :math:`x^{-i}`, namely the variables of other players. By renaming each player's :math:`i` variables :math:`y` to :math:`x^i`, and parameters :math:`x` to :math:`x^{-i}`, a Nash equilibrium for the game corresponds to a solution to following LCP:

.. math::
    0 \le x \perp M^ix^i + N^ix^{-i} + q^i \ge 0 \quad \forall i=\{1,2,\dots,n\}

====================================
Modeling the problem
====================================

We now model a simple simultaneous game among two players.
Their two optimization problem are as follows.

**Player 1**

.. math::

 \min_{q_1}: 10 q_1 + 0.1 q_1^2 - (100 - (q_1+q_2)) q_1  =  1.1 q_1^2 - 90q_1 + q_1q_2

 \text{s.t:} \quad  q_1 >= 0


**Player 2**

.. math::

 \min_{q_2}: 5 q_2 + 0.2 q_2^2 - (100 - (q_1+q_2)) q_2     = 1.2 q_2^2 - 95 q_2 + q_2q_1

 \text{s.t:} \quad  q_2 >= 0

The above problem corresponds to a `Cournot Competition <https://en.wikipedia.org/wiki/Cournot_competition>`_ where the demand curve is given by :math:`P = a-BQ` where ``P`` is the market price and ``Q`` is the quantity in the market. A convex quadratic function gives the cost of production of both the producers in the quantity they produce. The solution to the problem is to find a Nash Equilibrium from which neither producer can deviate.
To handle this problem, first, we create two objects of :cpp:class:`MathOpt::QP_Param` to model each player's optimization problem, as parameterized by the other.

.. code-block:: c

    #include "zero.h"
    //[...] Your other code here
    GRBEnv env;
    arma::sp_mat Q(1, 1), A(0, 1), B(0, 1), C(1, 1);
    arma::vec b, c(1), d(1);
    b.set_size(0);

    Q(0, 0) = 2 * 1.1;
    C(0, 0) = 1;
    c(0) = -90;
    auto q1 = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, d, &env);

    Q(0, 0) = 2 * 1.2;
    c(0) = -95;
    auto q2 = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, d, &env);

    // We create a vector with the two QP_Params
    std::vector<shared_ptr<MathOpt::QP_Param>> q{q1, q2};

    / /Cast to abstract MP_Param
    std::shared_ptr<MathOpt::MP_Param>> MPCasted=std::dynamic_cast<MathOpt::MP_Param>(q);

Since we do not have any Market clearing constraints (`more mathematical details <https://faculty.wcas.northwestern.edu/~mdo738/textbook/dls_ch5.pdf>`_), we set empty matrices for them. If the problem does not have market-clearing constraints, the matrices must be input with zero rows and the appropriate number of columns.

.. code-block:: c

    arma::sp_mat MC(0, 2);
    arma::vec MCRHS;
    MCRHS.set_size(0);

Finally, we can instantiate the :cpp:class:`Game::NashGame` object by invoking the constructor.

.. code-block:: c

    Game::NashGame Nash = Game::NashGame(&env, MPCasted, MC, MCRHS);


The LCP problem to solve this nash game is then:

.. math::

 0 \le q_1 \perp 2.2 q_1 + q_2 - 90 \geq 0
 0 \le q_2 \perp q_1 + 2.4 q_2 - 95 \geq 0

The method :cpp:func:`Game::NashGame::FormulateLCP` formulates the above LCP.

.. code-block:: c

    arma::sp_mat M;
    arma::vec q;
    // Stores the complementarity pairs relationships
    perps Compl;
    // Compute the LCP conditions
    Nash.FormulateLCP(M, q, Compl);
    M.print();
    q.print();

Here ``M`` and ``q`` are such that the solution to the LCP :math:`0 \le x \perp Mx + q \ge 0` solves the Nash Game. These matrices can be written to a file and solved externally now.
Alternatively, one can pass it to the :cpp:class:`Game::LCP` class, and solve it natively. To achieve this, one can pass the above matrices to the constructor of the :cpp:class:`Game::LCP` class.

.. code-block:: c

        Game::LCP lcp = Game::LCP(&env, M, q, 1, 0);

More concisely, the class :cpp:class:`Game::LCP` offers a constructor with a NashGame as an argument. This way, one need not explicitly compute ``M``, ``q``.

.. code-block:: c

        Game::LCP lcp2 = Game::LCP(&env, Nash);

====================================
Computing solutions
====================================

We can now solve the instance of `Game::LCP`.

.. code-block:: c

    auto model = lcp.LCPasMIP();
    model.optimize();

As was the case with :cpp:func:`MathOpt::QP_Param::solveFixed`, the above function returns a
``unique_ptr`` to ``GRBModel``.

====================================
Checking solutions
====================================

The solution to this problem is  :math:`q_1=28.271028, q_2=27.803728`. In order to verify the solution, one can create a solution vector and solve each player's :cpp:class:`MathOpt::QP_Param` and check that the solution indeed matches.

.. code-block:: c

    arma::vec Nashsol(2);
    Nashsol(0) = model->getVarByName("x_0").get(GRB_DoubleAttr_X); // This is 28.271028
    Nashsol(1) = model->getVarByName("x_1").get(GRB_DoubleAttr_X); // This is 27.803728

    auto nashResp1 = Nash.respond(0, Nashsol);
    auto nashResp2 = Nash.respond(1, Nashsol);

    cout<<nashResp1->getVarByName("y_0").get(GRB_DoubleAttr_X)<<endl; // Should print 28.271028
    cout<<nashResp2->getVarByName("y_0").get(GRB_DoubleAttr_X)<<endl; // Should print 27.803728


If only one does not want the individual ``GRBModel`` handles but just wants to confirm that the problem is solved or provide a player with profitable deviation, one can just use cpp:func:`Game::NashGame::isSolved` function as follow.

.. code-block:: c

    unsigned int temp1 ; arma::vec temp2;
    // This should be true.
    cout<<Nash.isSolved(Nashsol, temp1, temp2);

If the :cpp:func:`Game::NashGame::isSolved` function returns false, then ``temp1`` and ``temp2`` respectively contain the player with profitable deviation and the more profitable strategy of the player.
Furthermore, note that just like :cpp:class:`MathOpt::QP_Param`, :cpp:class:`Game::NashGame` can also be saved and loaded from an external file.

.. code-block:: c

    Nash.save("dat/Nash.dat");
    Game::NashGame Nash2(&env);
    Nash2.load("dat/Nash.dat");
