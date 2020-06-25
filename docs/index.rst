Welcome
=====================================

ZERO is a multi-purpose game solver written in C++. Currently, it manages NASPs (EPECs) and Integer Programming Games.

This project is a private fork of `EPECSolve <https://github.com/ssriram1992/EPECsolve/>`.

This is a private fork of EPECSolve, that handles an outer approximation
for NASPs and Integer Programming games. Code to compute
mixed-equilibrium in linear EPECs. - `Base
Code <https://github.com/ssriram1992/EPECsolve/>`__ -
`Instances <https://github.com/ds4dm/EPECInstances>`__ with detailed
`mathematical
description <https://github.com/ds4dm/EPECInstances/blob/master/Description.pdf>`__
- `arXiv <https://arxiv.org/abs/1910.06452>`__ pre-print for the
original code


**Mandatory packages** to make ZERO work are:

-  `Armadillo <http://arma.sourceforge.net/>`__ (Version 9.6 or later
   recommended. Minimum 8.5 required.)

   -  BLAS
   -  ARPACK
   -  LAPACK

-  `Gurobi <https://www.gurobi.com/registration/download-reg>`__
   (Version 8.1 or later)
-  `gcc/g++ <https://gcc.gnu.org/>`__ (Tested on version 4.8. Must
   support C++14 and be compatible with your version of Gurobi)
   ``sudo apt install gcc`` will install gcc/g++ on an Ubuntu machine.
-  `CMake <https://cmake.org/>`__
   ``sudo apt install cmake`` will do the job for Debian-based distro.
-  `Boost <https://www.boost.org/>`__ Required for logging, commandline
   interface to solve files etc. Can produce a boost-free version if
   there is significant interest.

Recommended but not mandatory for the algorithm. (Some examples might have these dependancies)

-  `Rapid JSON <http://rapidjson.org/>`__ To export results and save
   example problem instances.
-  `DOxygen <http://www.doxygen.nl>`__ Only if you need documentation.
   ``sudo apt install doxygen`` will install DOxygen on an Ubuntu
   machine.
   
**LICENSE**: This code is distributed under the CC0 1.0 Universal License.


.. toctree::
   :caption: Quick Start
   :hidden:
 
   self
   support_files/compiling
   support_files/command_line

.. toctree::
   :caption: Tutorials
   :hidden:

   support_files/EPEC
   support_files/LCP
   support_files/NashLCP
   support_files/QP_Param
   
.. toctree::
   :caption: Quick References
   :maxdepth: 5
   :hidden:

   api/namespace_Algorithms
   api/namespace_Game
   api/namespace_Models
   api/namespace_Utils
