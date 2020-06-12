Compiling
=========

We assume you have downloaded the project from the github repository.

CMake
-----

-  Open ``CMakeList.txt``.
-  You will find two configurations. ``ROSSOBIANCO`` represents the
   remote server one, while the other is the local one.
-  Set ``BOOST_ROOT``, ``ARMA_ROOT``, and ``ARMA_LIB`` to your
   customized folders.
-  Set ``GUROBI_PATH_%OS`` to your local gurobi path (and version). Note
   that %OS corresponds to your local (or remote) operatin system.