
# EPECsolve for Outer Games
This is a private fork of EPECSolve, that handles an outer approximation for NASPs and Integer Programming games.
Code to compute mixed-equilibrium in linear EPECs. 
- [Base Code](https://github.com/ssriram1992/EPECsolve/)
- [Instances](https://github.com/ds4dm/EPECInstances) with detailed [mathematical description](https://github.com/ds4dm/EPECInstances/blob/master/Description.pdf)
- [arXiv](https://arxiv.org/abs/1910.06452) pre-print for the original code

## LICENSE
This code is distributed under the CC0 1.0 Universal license.

# Prerequisites

## Mandatory packages
- [Armadillo](http://arma.sourceforge.net/) (Version 9.6 or later recommended. Minimum 8.5 required.)
	* BLAS
	* ARPACK
	* LAPACK
- [Gurobi](https://www.gurobi.com/registration/download-reg) (Version 8.1 or later)
- [gcc/g++](https://gcc.gnu.org/) (Tested on version 4.8. Must support C++11 and be compatible with your version of Gurobi) `sudo apt install gcc ` will install gcc/g++ on an Ubuntu machine.
- [GNU make](https://www.gnu.org/software/make/) `sudo apt install make` will install GNU make on an Ubuntu machine.
- [Boost](https://www.boost.org/) Required for logging, commandline interface to solve files etc. Can produce a boost-free version if there is significant interest.

## Recommended but not mandatory for the algorithm. (Some examples might have these dependancies)
- [Rapid JSON](http://rapidjson.org/) To export results and save example problem instances.
- [DOxygen](http://www.doxygen.nl) Only if you need documentation. `sudo apt install doxygen ` will install DOxygen on an Ubuntu machine.

# Getting the documentation
One can generate two versions of documentation for this project..
- Use the simple version of documentation, if you are only interested in using this as a predefined library which you don't intend to edit. This version of the documentation gives a sufficiently detailed explanation of every class or function you might every have to use. To avail this version, run
```bash
make docSimple
```
- Use the complete documentation if you are interested in every implementation detail of the code.This gives a complete description of every private member and fields in every class, all of which might be useful if you want to edit the code in here. To avail this version, run
```bash
make docDetailed
```
- You can alternatively use the [online documentation](https://ssriram1992.github.io/EPECsolve/html/index.html)

# Compiling
- Download the [project](https://github.com/ssriram1992/EPECsolve/). If you do not have access, please email [sriram.sankaranarayanan@polymtl.ca](mailto:sriram.sankaranarayanan@polymtl.ca).

## Makefile
- Open `Makefile`. 
- Enter the path where you downladed in `EPEC_HOME`. This folder should contain folders like `docs/`, `src/`, `test/` etc.
- Enter the path to Boost in `BOOST_HOME`. Ensure that the path to corresponding include files and boost libraries are correct.
- Enter the path of your armadillo-installation in the line defining the variable `ARMA`. Typically, the location would be like `/opt/armadillo-code`.
- Enter the path of your Gurobi-installation in the line defining the variable `GUR`. Typically, the location would be like `/opt/gurobi/gurobi811/<Your OS>`.
- Run `make install` to make the necessary folders to store temporary files, binary files etc.
- Run `make compileEPEC` to compile. 
- Optionally run `make EPECtest` to run the set of unit tests. It should all succeed without a problem.
- Run `make` to create the binary.
- Run `./bin/EPEC -h` to get a list of command line options with which you can run the executable.

## CMake
- Open `CMakeList.txt`. 
- You will find two configurations. `ROSSOBIANCO` represents the remote server one, while the other is the local one.
- Set `BOOST_ROOT`, `ARMA_ROOT`, and `ARMA_LIB` to your customized folders. 
- Set `GUROBI_PATH_%OS` to your local gurobi path (and version). Note that %OS corresponds to your local (or remote) operatin system.

# Command-line interface
EPECSolve has a command line interface built on top of the standard modeling paradigm (see documentation for src/Models.cpp).
Once the executable is compiled via `make compileEPEC`, the user can run EPECSolve with: `./EPEC $options`. 

The following I/O options are available:

* `-i (--input) PathToInstance`: This **required** option specifies the path to the JSON instance file. Extension .json is automatically appended to the argument.
* `-s (--solution) PathAndFile`: Sets the output (path) filename of the JSON solution file, if any (see option `-w. (writelevel)`. If a JSON output is required, extension .json is automatically appended to the argument.
* `-l (--log) PathAndFile`: Sets the output path/filename of the log file.
* `-w (--writelevel) intValue`: Sets the verbosity parameter for the output solution file. *0*: output only a JSON solution file (if any). *1* Only a human readable file. *2* Both.
* `-m (--message) intValue`: Sets the verbosity parameter for the program itself. *1*: only info messages, *2* info and debug, *3* info, debug and trace messages.


Algorithmic parameters:

* `-a (--algorithm) intValue`: Specifies the type of algorithm employed to solve the instance. *0* is **fullEnumeration**, *1* is **innerApproximation**, *2* is **combinatorialPNE**. Note that the latter only works with pure equilibria.
* `-p (--pure) intValue`: Specifies whether the algorithm should seek for a pure equilibrium or not (note that **combinatorialPNE** will always seek for a pure solution) . *0* no requirement (either pure or mixed), *1* only pure.
* `-bigM intValue`: Specifies the whether LCP problems should be formulated with a LCP. *0* (default) for indicator constraints, *1* for bigM formulation.
* `-b (--bound) intValue`: Specifies if the final LCP model should be bounded by a *bigM* constant (see option `-boundBigM`) in their primal variables. *0* no bounding, *1* bounding with *bigM*.
* `-boundBigM intValue`: Specifies the *bigM* constant for the above param (see option `-b (--bound)`.
* `-t (--threads) intValue`: Specifies the number of threads Gurobi uses.
* `-tl (--timelimit) intValue`: Specifies the timeLimit (in seconds) for running the whole algorithmic procedure. 

Algorithmic parameters tailored on **innerApproximation**:

* `-ad (--add) intValue`: Specifies the *EPECAddPolyMethod* strategy employed to expand the inner approximation when no information about deviations is available. *0* add polyhedra in lexographic order, *1* reverse lexographic, *2* random.
* `-ag (--aggr) intValue`: Specifies the *EPECAddPolyMethod* aggressivity (see option `-ad (--add)`). *int* specifies the number of lower-level polyhedra to be added to each Stackelberg leader at each iteration in which *EPECAddPolyMethod* is triggered.
* `-r (--recover) intValue`: Specifies the recover stragegy for the **innerApproximation** method. When the algorithm finds a mixed equilibrium and the parameter `-p (--pure)` is set to *1*, the algorithm will search for a pure equilibrium either by randomly expanding the inner approximation (*EPECAddPolyMethod*) or by triggering a **combinatorialPNE** procedure warmstarted with the current information. *0* sets the recover strategy to *incrementalEnumeration*, *1* sets it to *combinatorial*.

# Maintenance
[@gdragotto](https://github.com/gdragotto) - Contact: [gabriele.dragotto@polymtl.ca](mailto:gabriele.dragotto@polymtl.ca)
[@ssriram1992](https://github.com/ssriram1992/) - Contact: [sriram.sankaranarayanan@polymtl.ca](mailto:sriram.sankaranarayanan@polymtl.ca)

