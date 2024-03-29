What is an RBG?
=====================================
A **Reciprocally Bilinear Game** is a simultaneous non-cooperative game among :math:`n` players with each player :math:`i=1,2,\dots,n` solving the optimization problem

.. math::
    \min{x^i} (c^i)^\top x^i + (x^{-i})^\top C^ix^i + d^Tx^{-i} \\
    \text{s.t.} \quad  x^i \in \mathcal{X}^i

where :math:`\mathcal{X}^i` is a set (not necessarily closed), :math:`C` and :math:`c` are a matrix and a vector of appropriate dimensions. An RBG is *polyhedrally representable* if :math:`\text{cl conv}(\mathcal{X}^i)` is a polyhedron for every :math:`i`, and one can optimize an arbitrary linear funciton on :math:`\mathcal{X}^i`.
As a standard game-theory notation, the operator :math:`(\cdot)^i` refers to an object of player :math:`i` with :math:`i \in \{ 1,2,\dots,n\}`, and :math:`(\cdot)^{-i}` be :math:`(\cdot^1,\dots, \cdot^{i-1},\cdot^{i+1},\dots,\cdot^{n})`.

From the definition, the following properties hold for each player :math:`i`:

    * its objective function is *reciprocally bilinear*, namely, it is linear in its variables :math:`x^i`, and bilinear in the other players' ones :math:`x^{-i}`
    * its constraint set :math:`\mathcal{X}^i` is not parametrized in :math:`x^{-i}`, i.e., the interaction takes place at the objective level thus the problem is not a generalized Nash equilibrium problem.

The set :math:`\text{cl conv}(\mathcal{X}^i)` as it represents the set of all **mixed strategies** that the player can adopt.
We refer to [CNP]_ for a detailed review on the mathematics.

.. [CNP] Margarida Carvalho, Gabriele Dragotto, Andrea Lodi, Sriram Sankaranarayanan. The Cut and Play Algorithm: Computing Nash Equilibria via Outer Approximations.
