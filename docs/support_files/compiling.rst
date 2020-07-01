Compiling
=====================================


Dependencies
*************

ZERO requires a few external libraries to work. Most of these are managed via the `Conan <https://conan.io/>`__ package manager for C++.
A mandatory requirement to use the software as-it-is is `CMake <https://cmake.org/>`__, which will manage the integration with Conan by itself.

The required packages for ZERO are:

-  `CMake <https://cmake.org/>`__ (Version 3.14 or later)
-  `Conan <https://conan.io/>`__  (Version 1.26 or later)
-  `Gurobi <https://www.gurobi.com/registration/download-reg>`__ (Version 8.1 or later) as a MIP solver.

The following packages will be installed by Conan

-  `Armadillo <http://arma.sourceforge.net/docs.html>`__  (Version 9.8 or later) for the linear algebra.

-  `openBLAS <https://www.openblas.net>`__  (Version 0.3.10 or later) as a requirement for Armadillo.

-  `Boost <https://www.boost.org/>`__ (Version 1.72 or later) required for logging, shell interfaces and utilities.

-  `RapidJSON <http://rapidjson.org>`__ (Version 1.1.0 or later) suggested but not mandatory for easy management of JSON output files.

Recommended but not mandatory for the documentation

-  `DOxygen <http://www.doxygen.nl>`__ (Version 1.8 or later) to generate the C++ documentation
-  `Sphinx <http://www.sphinx-doc.org>`__ (Version 1.8 or later) to render the HTML into a ReadTheDocs template
-  `Exhale <https://exhale.readthedocs.io/en/latest/>`__ (Version 0.2.3 or later) and `Breathe <https://breathe.readthedocs.io/en/latest/>`__ (Version 4.10 or later)  to integrate Doxygen, Sphinx and the ReadTheDocs template



Downloading and Compiling
*********************************
First, clone the repository of ZERO from GitHub with:

.. code-block:: python

   git clone https://github.com/ds4dm/ZERO

To build the targets with CMake, then:

.. code-block:: python

   cd ZERO/
   cmake .. && make