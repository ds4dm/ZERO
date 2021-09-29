NASPs
***************
NASPs, an acronym for Nash games among Stackelberg Players, are games among bilevel leaders.
A detailed description of these games is available `in this paper <https://arxiv.org/abs/1910.06452>`_.

Each player :math:`i` solves the optimization problem in :math:`x^i`, where :math:`c^i` is a real-valued vector of dimension :math:`n_i`, and :math:`Q_o^i` is a :math:`n_o \times n_i` real-valued matrix encapsulating the interactions between any two distinct players :math:`i` and :math:`o`.
Any leader :math:`i` has :math:`m_i` followers, each of which solves a convex continuous optimization problem.
For a given leader :math:`i`, and its respective follower :math:`j \in \{1,\dots,m_i\}`, :math:`f^{i,j}`, :math:`e^{i,j}` and :math:`D^{i,j}`, :math:`E^{i,j}`, :math:`F^{i,j}`, :math:`G^{i,j}`, :math:`H^{i,j}` are respectively vectors and matrices of conformable dimensions.

.. math::

    \min_{x^i} (c^i)^\top x^i + \sum_{o=1, o\neq i}^{n} (x^o)^\top Q^i_o x^i \\
    \text{s.t.}  \quad A^ix^i =b^i \\
    \quad \quad y^{i,j} \in \text{arg}\min_{y^{i,j}} \{ (0.5\cdot{}(y^{i,j})^\top D^{i,j} y^{i,j} + (f^{i,j}+ \sum_{k=1,k\neq j}^{m_i}(y^{i,k})^\top E^{i,j})y^{i,j} + (Fw^i)^\top y^{i,j} : \\
    \quad \qquad \qquad \qquad \quad G^{i,j}w^i + H^{i,j}y^{i,j} \le e^{i,j}, \; y^{i,j} \ge 0 \} \quad \forall j \in \{1,\dots,m_i\} \\
    \quad \quad x^i \ge 0, \; x^i=(w^i,y^i)


In the above model, we assume the variables :math:`x^i` are partitioned into leader's variables :math:`w^i` and followers' ones :math:`y^i`.
In plain English, :math:`n` leaders play a simultaneous non-cooperative game. Each player is a leader, with :math:`m_i` followers. While leaders interact, each follower interacts only with the other followers of its leader (and not those of other leaders).


====================================
A quick example
====================================

The first player is the **u-v** player, where the leader's decision variables are :math:`u` and the follower's decision variables
are :math:`v`. The second player is the **x-y** player where the leader's and the
follower's variables are :math:`x` and :math:`y` respectively.

The u-v player solves the following problem:

.. math::

  min_{u,v} v_1 -u + y_1v_2&\qquad

  \text{s.t.} \quad u \quad&\ge\quad 0

  \quad \quad v_1+v_2+u \quad&\leq\quad 5

  \quad \quad v \quad&\in\quad \arg \min _v \left \{ -v_1+v_2 : v \ge 0; 2v_1+v_2 \leq u; v_1 -2v_2 \leq -u \right \}


While the x-y player solves:

.. math::

   \min_{x,y}  y_1 - x + uy_2&\qquad

   \text{s.t.} \quad  x \quad&\ge\quad 0

    \quad \quad y_1 + y_2 + x \quad&\le\quad 7

    \quad \quad -y_1 + y_2 \quad&\le\quad 0

    \quad \quad y\quad&\in\quad\arg\min_y \left\{ y_1 - y_2: y \ge 0; -y_1 + y_2 \le 5-x; -y_1 + y_2 \le x-3 \right\}


The problem has a pure-strategy Nash equilibrium given by :math:`(u, v_1, v_2) = (2.78, 0.56, 1.67)`, and :math:`(x, y_1, y_2) = (1.67, 1.33, 0)`.

====================================
Modeling the problem
====================================

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Creating an inheritor class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step in modeling a game between Stackelberg leaders is to include `zero.h` and create a derived class of :cpp:class:`Game::EPEC`. The minimal constructor for :cpp:class:`Game::EPEC` involves passing a pointer to `GRBEnv` (Check Gurobi's C++ `reference manual <https://www.gurobi.com/documentation/8.1/refman/cpp_api_overview.html>`_). The derived class should indeed instantiate the base class (Game::EPEC) using such a constructor.

.. code-block:: c

    #include "zero.h"
    class my_Prob : public Game::EPEC
    {
    public:
        my_Prob(GRBEnv *e) : Game::EPEC(e) {}
    };


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The followers' problems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We define the lower level of each leader (u-v leader as well as the x-y leader) as a :cpp:class:`Game::NashGame` object. For convenience, we write the following two functions that return a ``std::shared_ptr<Game::NashGame>``.

.. note::
    #. The referred object contains the follower's game along with any constraint in the leader level.
    #. The referred object does not contain the follower's objective (which could depend on other leaders' variables).
    #. We create the object, **without** assuming the presence of other leaders.

The following code returns the ``std::shared_ptr<>`` as required. Refer to the previous tutorial on simultaneous games to learn how to create a :cpp:class:`Game::NashGame` object.

.. code-block:: c

    std::shared_ptr<Game::NashGame> uv_leader(GRBEnv *env) {
      // 2 variable and 2 constraints
      arma::sp_mat Q(2, 2), C(2, 1), A(2, 1), B(2, 2);
      arma::vec    c(2, arma::fill::zeros);
      arma::vec    b(2, arma::fill::zeros);
      // Q remains as 0
      // C remains as 0
      // c
      c(0) = -1;
      c(1) = 1;
      // A
      A(0, 0) = -1;
      A(1, 0) = 1;
      // B
      B(0, 0)   = 2;
      B(0, 1)   = 1;
      B(1, 0)   = 1;
      B(1, 1)   = -2;
      auto foll = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, env);

      // Lower level Market clearing constraints - empty
      arma::sp_mat MC(0, 3);
      arma::vec    MCRHS(0, arma::fill::zeros);

      arma::sp_mat LeadCons(1, 3);
      arma::vec    LeadRHS(1);
      LeadCons(0, 0) = 1;
      LeadCons(0, 1) = 1;
      LeadCons(0, 2) = 1;
      LeadRHS(0)     = 5;

      std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
      MPCasted.push_back(std::dynamic_pointer_cast<MathOpt::MP_Param>(foll));

      auto N = std::make_shared<Game::NashGame>(env, MPCasted, MC,
                MCRHS, 1, LeadCons, LeadRHS);
      return N;
 }


Next, we have a similar procedure for the x-y leader.

.. code-block:: c

    std::shared_ptr<Game::NashGame> xy_leader(GRBEnv *env) {
      // 2 variable and 2 constraints
      arma::sp_mat Q(2, 2), C(2, 1), A(2, 1), B(2, 2);
      arma::vec    c(2, arma::fill::zeros);
      arma::vec    b(2, arma::fill::zeros);
      // Q remains as 0
      // C remains as 0
      // c
      c(0) = 1;
      c(1) = -1;
      // A
      A(0, 0) = 1;
      A(1, 0) = -1;
      // B
      B(0, 0) = -1;
      B(0, 1) = 1;
      B(1, 0) = -1;
      B(1, 1) = 1;
      // b
      b(0)      = 5;
      b(1)      = -3;
      auto foll = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, env);

      // Lower level Market clearing constraints - empty
      arma::sp_mat MC(0, 3);
      arma::vec    MCRHS(0, arma::fill::zeros);

      arma::sp_mat LeadCons(2, 3);
      arma::vec    LeadRHS(2);
      LeadCons(0, 0) = 1;
      LeadCons(0, 1) = 1;
      LeadCons(0, 2) = 1;
      LeadRHS(0)     = 7;
      // Comment the following four lines for another example ;)
      LeadCons(1, 0) = -1;
      LeadCons(1, 1) = 1;
      LeadCons(1, 2) = 0;
      LeadRHS(1)     = 0;

      std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
      MPCasted.push_back(std::dynamic_pointer_cast<MathOpt::MP_Param>(foll));

      auto N = std::make_shared<Game::NashGame>(env, MPCasted, MC, MCRHS, 1, LeadCons, LeadRHS);
      return N;
 }


We introduce a member function to add the leaders to the class.

.. code-block:: c

  void addLeader(std::shared_ptr<Game::NashGame> N, const unsigned int i) {
    this->PlayersLowerLevels.push_back(N);
    ends[i] = N->getNprimals() + N->getNumLeaderVars();
    this->LocEnds.push_back(&ends[i]);
  }


.. note::

    The above code performs the following operations, which should always be present:
    * The lower-level `Game::NashGame` is pushed to ``Game::EPEC::PlayersLowerLevels``
    * Variables that track the number of variables in the current leader (``ends[i]``) is set and is tracked by ``Game::EPEC::LocEnds`` at the appropriate position.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Re-implementing methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:cpp:class:`Game::EPEC` is a pure virtual (abstract) class and it is mandatory to define two functions by every derived class that it has. First, we define :cpp:func:`Game::EPEC::makeObjectivePlayer`. This function has the following signature in its definition in :cpp:class:`Game::EPEC`.

.. code-block:: c

    virtual void makeObjectivePlayer(const unsigned int i, Game::QP_objective &QP_obj) = 0;


The parameter ``i``takes the leader's position and `QP_obj` is an out-parameter, which should be filled with an object of ``MathOpt::QP_objective``, which has the i-th leader's objective. Note that this should assume the form of :math:`c^T x + (Cx)^T x^{oth}`, where :math:`x` is the current player's set of variables and :math:`x^{oth}` is the remaining set of variables.

.. code-block:: c

    void my_Prob::makeObjectivePlayer(const unsigned int i, Game::QP_objective &QP_obj) override
    {
        QP_obj.Q.zeros(3, 3);
        QP_obj.C.zeros(3, 3);
        QP_obj.c.zeros(3);
        switch (i)
        {
            case 0: // uv_leader's objective
                QP_obj.C(1, 0) = 1;
                QP_obj.c(0) = 1;
                QP_obj.c(2) = -1;
            break;
            case 1: // xy_leader's objective
                QP_obj.C(1, 2) = 1;
                QP_obj.c(0) = 1;
                QP_obj.c(2) = 1;
            break;
             default: // Not strictly required, but for safety
                throw std::string("Invalid makeObjectivePlayer");
        }
    }


Finally,  `Game::EPEC::updateLocations` needs to be implemented.
For small toy examples, this function can only update the location of the last variable as the total number of variables defined by the user plus any convex hull variables. But, for more complicated examples, we refer the user to check :cpp:func:`Models::EPEC::updateLocations`.
Also, :cpp:func:`Game::EPEC::preFinalize` and :cpp:func:`Game::EPEC::postFinalize` are required in the derived class. These methods are called before and after :cpp:func:`Game::EPEC::finalize`.

.. code-block:: c

    void My_EPEC_Prob::updateLocations() override {
        ends[0] = this->convexHullVariables.at(0) + 3;
        ends[1] = this->convexHullVariables.at(1) + 3;
    }
    void My_EPEC_Prob::postFinalize() override { std::cout << "Pre finalized!\n"; }
    void My_EPEC_Prob::preFinalize() override { std::cout << "Post finalized!\n"; };


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Computing solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Now that the derived class is ready, the EPEC can be solved using an instantiation of the class.

To start, with set up a Gurobi environment like we did for :cpp:class:`MathOpt::QP_Param` and :cpp:class:`Game::NashGame`.

.. code-block:: c

    GRBEnv env;

We can then specify the log level via `loguru`.

.. code-block:: c

    //0 is info. The greater, the more verbose
    loguru::g_stderr_verbosity = 0;

Next, we create an object for the class and add the lower level :cpp:class:`Game::NashGame` using functions defined earlier.

.. code-block:: c

    // Create the class object
    My_EPEC_Prob epec(&env);
    // Adding uv_leader
    auto uv_lead = uv_leader(&env);
    epec.addLeader(uv_lead, 0);
    // Adding xy_leader
    auto xy_lead = xy_leader(&env);
    epec.addLeader(xy_lead, 1);


Once all the leaders' lower levels are in, we tell the program that we are adding no more players, and the code can do certain pre-processing and space allocation using :cpp:func:`Game::EPEC::finalize`. We can also optionally tell the program to do other operations before/after finalizing, by defining an override for :cpp:func:`Game::EPEC::preFinalize` and :cpp:func:`Game::EPEC::postFinalize` in the derived class.

.. code-block:: c

    // Finalize
    epec.finalize();

One can optionally choose the algorithm to solve the problem. Not setting this, chooses the default algorithm cpp:class:`Algorithms::EPEC::FullEnumeration`

.. code-block:: c

    epec.setAlgorithm(Data::EPEC::Algorithms::InnerApproximation);


Finally, we can solve the problem.

.. code-block:: c

      // Solve
      try {
         epec.findNashEq();
      } catch (ZEROException &e) {
         std::cerr << e.what() << " -- " << std::to_string(e.which()) << std::endl;
      }


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fetching solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, one can write the GRBModel (Gurobi model) solved in the last iteration or acquire a copy of the model. For the model writing, any extension allowed by Gurobi will work in the solver.

.. code-block:: c

    // Writes the model to a file. The model can then be loaded externally, resolved and analyzed.
    // Writes to an LP file, in a human readable format
    epec.writeLCPModel("my_model.lp");
    // Writes to an MPS file, in a machine readable format
    epec.writeLCPModel("my_model.sol");
    // Writes the solution to the same model.

    // Human and machine readable.
    epec.writeLCPModel("my_model.sol");


Alternatively, without saving the model, one can directly print the solution to the model.
Note that an EPEC does not necessarily have a pure-strategy Nash equilibrium or a mixed-strategy Nash equilibrium.
However, should it have one, we print the multiple pure strategies and the associated probability for that strategy. One can perform such queries with:

- :cpp:func:`Algorithms::EPEC::PolyBase::getValProbab`
- :cpp:func:`Algorithms::EPEC::PolyBase::getValLeadLeadPoly`
- :cpp:func:`Algorithms::EPEC::PolyBase::getValLeadFollPoly`

.. code-block:: c

    // Get the set of pure strategies that the leaders will play
    auto uv_strats = epec.mixedStrategyPoly(0);
    // Now print the probability of each such pure strategy and the actual strategy too.
    std::for_each(
    std::begin(uv_strats), std::end(uv_strats), [&epec](const unsigned int i) {
        // epec.getValProbab (a, b) gives the probability used to play b-th pure strategy by the player at position a.
        std::cout << "With probability  " << epec.getValProbab(0, i) << '\n';
        // epec.getValLeadLeadPoly(a, b, c) gives the bth variable of a-th leader in c-th poly.
        std::cout << "(" << epec.getValLeadLeadPoly(0, 0, i) << ", "
        // epec.getValLeadFollPoly(a, b, c) gives the bth follower variable of a-th leader in c-th poly.
            << epec.getValLeadFollPoly(0, 0, i) << ", "
            << epec.getValLeadFollPoly(0, 1, i) << ")\n";
    });

Similarly, for the x-y leader:

.. code-block:: c

    auto xy_strats = epec.mixedStrategyPoly(1);
    std::for_each(
    std::begin(xy_strats), std::end(xy_strats), [&epec](const unsigned int i) {
        std::cout << "With probability  " << epec.getValProbab(1, i) << '\n';
        std::cout << "(" << epec.getValLeadLeadPoly(1, 0, i) << ", "
            << epec.getValLeadFollPoly(1, 0, i) << ", "
            << epec.getValLeadFollPoly(1, 1, i) << ")\n";
    });

The entire example source code is as follows:

.. code-block:: c

    #include "zero.h"

    class My_EPEC_Prob : public Game::EPEC {
    public:
      My_EPEC_Prob(GRBEnv *e) : EPEC(e) {}
      void addLeader(std::shared_ptr<Game::NashGame> N, const unsigned int i) {
         this->PlayersLowerLevels.push_back(N);
         ends[i] = N->getNprimals() + N->getNumLeaderVars();
         this->LocEnds.push_back(&ends[i]);
      }
      void postFinalize() override { std::cout << "Pre finalized!\n"; }
      void preFinalize() override { std::cout << "Post finalized!\n"; };

    private:
      unsigned int ends[2];
      void         updateLocations() override {
        ends[0] = this->ConvexHullVariables.at(0) + 3;
        ends[1] = this->ConvexHullVariables.at(1) + 3;
      }
      void makeObjectivePlayer(const unsigned int i, MathOpt::QP_objective &QP_obj) override {
         QP_obj.Q.zeros(3, 3);
         QP_obj.C.zeros(3, 3);
         QP_obj.c.zeros(3);
         switch (i) {
         case 0: // uv_leader's objective
            QP_obj.C(1, 0) = 1;
            QP_obj.c(0)    = 1;
            QP_obj.c(2)    = -1;
            break;
         case 1: // xy_leader's objective
            QP_obj.C(1, 2) = 1;
            QP_obj.c(0)    = 1;
            QP_obj.c(2)    = 1;
            break;
         default:
            throw std::string("Invalid makeObjectivePlayer");
         }
      }
    };

    std::shared_ptr<Game::NashGame> uv_leader(GRBEnv *env) {
      // 2 variable and 2 constraints
      arma::sp_mat Q(2, 2), C(2, 1), A(2, 1), B(2, 2);
      arma::vec    c(2, arma::fill::zeros);
      arma::vec    b(2, arma::fill::zeros);
      // Q remains as 0
      // C remains as 0
      // c
      c(0) = -1;
      c(1) = 1;
      // A
      A(0, 0) = -1;
      A(1, 0) = 1;
      // B
      B(0, 0)   = 2;
      B(0, 1)   = 1;
      B(1, 0)   = 1;
      B(1, 1)   = -2;
      auto foll = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, env);

      // Lower level Market clearing constraints - empty
      arma::sp_mat MC(0, 3);
      arma::vec    MCRHS(0, arma::fill::zeros);

      arma::sp_mat LeadCons(1, 3);
      arma::vec    LeadRHS(1);
      LeadCons(0, 0) = 1;
      LeadCons(0, 1) = 1;
      LeadCons(0, 2) = 1;
      LeadRHS(0)     = 5;

      std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
      MPCasted.push_back(std::dynamic_pointer_cast<MathOpt::MP_Param>(foll));

      auto N = std::make_shared<Game::NashGame>(env, MPCasted, MC, MCRHS, 1, LeadCons, LeadRHS);
      return N;
    }

    std::shared_ptr<Game::NashGame> xy_leader(GRBEnv *env) {
      // 2 variable and 2 constraints
      arma::sp_mat Q(2, 2), C(2, 1), A(2, 1), B(2, 2);
      arma::vec    c(2, arma::fill::zeros);
      arma::vec    b(2, arma::fill::zeros);
      // Q remains as 0
      // C remains as 0
      // c
      c(0) = 1;
      c(1) = -1;
      // A
      A(0, 0) = 1;
      A(1, 0) = -1;
      // B
      B(0, 0) = -1;
      B(0, 1) = 1;
      B(1, 0) = -1;
      B(1, 1) = 1;
      // b
      b(0)      = 5;
      b(1)      = -3;
      auto foll = std::make_shared<MathOpt::QP_Param>(Q, C, A, B, c, b, env);

      // Lower level Market clearing constraints - empty
      arma::sp_mat MC(0, 3);
      arma::vec    MCRHS(0, arma::fill::zeros);

      arma::sp_mat LeadCons(2, 3);
      arma::vec    LeadRHS(2);
      LeadCons(0, 0) = 1;
      LeadCons(0, 1) = 1;
      LeadCons(0, 2) = 1;
      LeadRHS(0)     = 7;
      // Comment the following four lines for another example ;)
      LeadCons(1, 0) = -1;
      LeadCons(1, 1) = 1;
      LeadCons(1, 2) = 0;
      LeadRHS(1)     = 0;

      std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
      MPCasted.push_back(std::dynamic_pointer_cast<MathOpt::MP_Param>(foll));

      auto N = std::make_shared<Game::NashGame>(env, MPCasted, MC,
                MCRHS, 1, LeadCons, LeadRHS);
      return N;
    }

    int main() {
      GRBEnv env;
      loguru::g_stderr_verbosity = 0;
      My_EPEC_Prob epec(&env);
      // Adding uv_leader
      auto uv_lead = uv_leader(&env);
      epec.addLeader(uv_lead, 0);
      // Adding xy_leader
      auto xy_lead = xy_leader(&env);
      epec.addLeader(xy_lead, 1);
      // Finalize
      epec.finalize();
      epec.setAlgorithm(Data::EPEC::Algorithms::InnerApproximation);
      // Solve
      try {
         epec.findNashEq();
      } catch (ZEROException &e) {
         std::cerr << e.what() << " -- " << std::to_string(e.which()) << std::endl;
      }

      std::cout << "\nUV LEADER\n";
      std::cout << "u: " << epec.getValLeadLead(0, 0) << '\n';
      std::cout << "v_1: " << epec.getValLeadFoll(0, 0) << '\n';
      std::cout << "v_2: " << epec.getValLeadFoll(0, 1) << '\n';
      auto uv_strats = epec.mixedStrategyPoly(0);
      std::for_each(std::begin(uv_strats), std::end(uv_strats), [&epec](const unsigned int i) {
         std::cout << "With probability  " << epec.getValProbab(0, i) << '\n';
         std::cout << "(" << epec.getValLeadLeadPoly(0, 0, i) << ", " << epec.getValLeadFollPoly(0, 0, i)
                      << ", " << epec.getValLeadFollPoly(0, 1, i) << ")\n";
      });
      std::cout << '\n';
      std::cout << "\nXY LEADER\n";
      std::cout << "x: " << epec.getValLeadLead(1, 0) << '\n';
      std::cout << "y_1: " << epec.getValLeadFoll(1, 0) << '\n';
      std::cout << "y_2: " << epec.getValLeadFoll(1, 1) << '\n';
      auto xy_strats = epec.mixedStrategyPoly(1);
      std::for_each(std::begin(xy_strats), std::end(xy_strats), [&epec](const unsigned int i) {
         std::cout << "With probability  " << epec.getValProbab(1, i) << '\n';
         std::cout << "(" << epec.getValLeadLeadPoly(1, 0, i) << ", " << epec.getValLeadFollPoly(1, 0, i)
                      << ", " << epec.getValLeadFollPoly(1, 1, i) << ")\n";
      });
      std::cout << '\n';
      return 0;
    }
