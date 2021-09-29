##################
Overview
##################

ZERO's components -- or *modules* -- are abstract objects defined inside a suitable *namespace*,  which groups modules with a similar function or goal.
 
************************************
MathOpt
************************************

The `MathOpt <../api/namespace_MathOpt.html#namespace-MathOpt>`_ namespace contains the utilities and class managing the mathematical optimization objects that ZERO supports. The main class are:


MP_Program(s)
===========================

-  :cpp:class:`MathOpt::MP_Param` is an abstract class managing the parametrized convex quadratic problems (quadratic objective function with linear constraints).
-  :cpp:class:`MathOpt::QP_Param` an inheritor class managing parametrized quadratic programs.
-  :cpp:class:`MathOpt::IP_Param` an inheritor class managing parametrized bi-linear integer programs.



LCPs
===========================

The class :cpp:class:`MathOpt::LCP` manages `Linear Complementarity Problems <https://en.wikipedia.org/wiki/Linear_complementarity_problem>`_ (*LCPs*). Such problems are fundamental components for the computation of Nash Equilibria in games among quadratic linear programs. _LCPs_ can be seen as "feasibility" mathematical programs whose feasible region is given by a finite union of polyhedra.

-  :cpp:class:`MathOpt::LCP` is the base class for LCP problems.
-  :cpp:class:`MathOpt::PolyLCP` extends the LCP class to handle its polyhedral aspects. For instance, outer and inner approximations of the feasible region.

Additionally, the namespace includes some utilities and shared components (e.g., :cpp:func:`MathOpt::convexHull`, :cpp:func:`MathOpt::getDualMembershipLP`).



************************************
Games
************************************

The `Game <../api/namespace_Game.html#namespace-Game>`_  namespace contains the definitions of various types of games.

EPECs
===========================
Equilibrium Problems with Equilibrium Constraints (*EPECs*) are powerful modeling paradigms for modeling two-level hierarchical games involving equilibria both in the upper and lower level.
ZERO handles EPECs where:

-	Each ("player") is the leader of a Bilevel program
-   Each leader has a set of followers playing a simultaneous game among each other. Each follower solves a convex quadratic continuous problem
- 	While leaders can interact among themselves, followers can only interact with other followers from the same leader.

For more information, take a look at the material `here <https://dragotto.net/when-nash-meets-stackelberg>`_.
:cpp:class:`Game::EPEC` has a generic implementation for EPEC problems. You may want to look at :cpp:class:`Models::EPEC` for the associated modeling paradigm.

ZERO currently supports the following algorithms for EPECs:

-  :cpp:class:`Algorithms::EPEC::FullEnumeration` a full-enumeration algorithm for EPECs (more here_). This is an exact method capable of finding an equilibrium maximizing a given function in the players' variables.
-  :cpp:class:`Algorithms::EPEC::InnerApproximation` an inner approximation algorithm for EPECs (more here_). This algorithm inner-approximate each player feasible region (namely, it approximates an LCP and hence a finite union of polyhedra) with increasingly bigger representations.
-  :cpp:class:`Algorithms::EPEC::CombinatorialPNE` a sort of inner approximation algorithm for EPECs (more here_). This algorithm only computes pure equilibria by inner approximating each player's finite union of polyhedra with a single polyhedron. All the combinations are tested for a PNE.
-  :cpp:class:`Algorithms::EPEC::CutAndPlay` outer approximates each player's feasible region with an increasingly tight polyhedral approximation. It uses the `CutAndPlay <https://dragotto.net/the-equilibrium-oracle>`_ algorithm scheme.


IPG
===========================
Integer Programming Games (*IPGs*) are a class of simultaneous non-cooperative games where each player solves an integer program fully parametrized in the player's variable. In other words, the constraints of each player solely depend on the variables of that player. ZERO currently supports IPGs where the objective functions are at most bilinear (w.r.t. the other players) and quadratic in each player's variables.
:cpp:class:`Game::IPG` implements such games, and there is only one algorithm available:

-  :cpp:class:`Algorithms::IPG::CutAndPlay` an outer approximation algorithm (similar to :cpp:class:`Algorithms::EPEC::CutAndPlay`) again using the `CutAndPlay <https://dragotto.net/the-equilibrium-oracle>`_ algorithm scheme.. Each player's integer program is approximated with its linear relaxation. A sequence of cutting planes and branching decisions guide the algorithm towards the computation of an equilibrium.

************************************
Algorithms
************************************
The namespace `Algorithms <../api/namespace_Algorithms.html#namespace-Algorithms>`_ works as a container for any algorithm of ZERO.
The children namespace are named after the respective game (e.g, `IPG <../api/namespace_Algorithms__IPG.html#namespace-algorithms-ipg>`_)


************************************
Models
************************************
The namespace `Models <../api/namespace_Models.html#namespace-Models>`_  implements some high-level APIs to access the :cpp:class:`Game` modules.

-  :cpp:class:`Models::EPEC::EPEC` high-level access to :cpp:class:`Game::EPEC`.
-  :cpp:class:`Models::IPG::IPG` high-level access to :cpp:class:`Game::IPG`.


************************************
Utils
************************************
The namespace `Utils <../api/namespace_Utils.html>`_  is an utility namespace with some common methods (e.g., numerical tolerances, comparisons, etc).


************************************
Solvers
************************************
The namespace `Solvers <../api/namespace_Solvers.html#namespace-Solvers>`_  implements some external solvers through a customized interface (e.g., PATH).
In future releases, we hope to move also Gurobi and other solvers such as SCIP in this namespace and abstract their interface with ZERO.


************************************
Shells
************************************
There are currently two command-line interfaces for IPGs and EPECs. You can find them in the directory *shells*:

- The `EPEC Shell <../api/file_src_shells_epec.cpp.html>`_
- The `IPG Shell <../api/file_src_shells_ipg.cpp.html>`_


************************************
Other stuff
************************************
- :cpp:class:`ZEROAlgorithmData` is an abstract class to manage the data of the various algorithms, for instance, :cpp:class:`Data::EPEC::DataObject` and :cpp:class:`Data::IPG::DataObject`.
- :cpp:class:`ZEROException` is custom exception class. :cpp:func:`ZEROAssert` is a custom assertion global method.
- :cpp:enum:`ZEROStatus` is an enum of status codes for successful computations/executions, and :cpp:enum:`ZEROErrorCodes` is an enum of error codes.
- `version <../api/file_include_support_version.h.html#file-include-support-version-h>`_ defines some macros related to the versioning 
