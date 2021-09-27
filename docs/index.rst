Welcome
=====================================

ZERO is a modular C++ library interfacing Mathematical Programming and Game Theory.
It provides a comprehensive toolkit of modeling interfaces, helper tools, and algorithms to design games and find their Nash equilibria.
In particular, the software supports Reciprocally Bilinear Games (RBGs), simultaneous non-cooperative games where each player solves a mathematical program with a linear objective in the player's variable and bilinear in its opponents' variables.
This class elementary generalizes to a multi-agent setting the tasks of classical Operations Research problems. ZERO also supports integer non-convexities, linear bilevel problems, and linear equilibrium problems with equilibrium constraints.


License
**********************
This code is distributed under the `CC BY-NC-SA 4.0 License <https://github.com/ds4dm/ZERO/blob/master/LICENSE>`_.


Citations
**********************
Cite ZERO through [ZERO]_, and any of the algorithms via their respective paper (e.g., [CNP]_ for the Cut-And-Play algorithm).

.. [ZERO] Sriram Sankaranarayanan, Gabriele Dragotto, Margarida Carvalho, Andrea Lodi. ZERO: Playing Mathematical Programming Games.
.. [WNMS] Margarida Carvalho, Gabriele Dragotto, Felipe Feijoo, Andrea Lodi, Sriram Sankaranarayanan. When Nash Meets Stackelberg.
.. [CNP] Margarida Carvalho, Gabriele Dragotto, Andrea Lodi, Sriram Sankaranarayanan. The Cut and Play Algorithm: Computing Nash Equilibria via Outer Approximations.


.. toctree::
   :caption: Quick Start
   :name: quickstart
   :hidden:

   self
   support_files/overview
   support_files/compiling

.. toctree::
   :caption: Tutorials
   :name: Tutorials
   :hidden:

   support_files/tutorials/QP_Param
   support_files/tutorials/IP_Param
   support_files/tutorials/NashLCP
   support_files/tutorials/LCP
   support_files/tutorials/NASP
   support_files/tutorials/IPG

.. toctree::
   :caption: Reference
   :name: Reference
   :hidden:

   support_files/reference/MathOpt
   support_files/reference/Game
   support_files/reference/Algorithms
   support_files/reference/Models
   support_files/reference/Utils
   support_files/reference/Codes




.. toctree::
   :caption: Extended Reference
   :maxdepth: 5
   :hidden:


   api/library_root

