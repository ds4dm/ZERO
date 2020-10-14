Command-line interface
=====================================

EPECSolve has a command line interface built on top of the standard
modeling paradigm (see documentation for src/Models.cpp). Once the
executable is compiled via ``make compileEPEC``, the user can run
EPECSolve with: ``./EPEC $options``.

The following I/O options are available:

-  ``-i (--input) PathToInstance``: This **required** option specifies
   the path to the JSON instance file. Extension .json is automatically
   appended to the argument.
-  ``-s (--solution) PathAndFile``: Sets the output (path) filename of
   the JSON solution file, if any (see option ``-w. (writelevel)``. If a
   JSON output is required, extension .json is automatically appended to
   the argument.
-  ``-l (--log) PathAndFile``: Sets the output path/filename of the log
   file.
-  ``-w (--writelevel) intValue``: Sets the verbosity parameter for the
   output solution file. *0*: output only a JSON solution file (if any).
   *1* Only a human readable file. *2* Both.
-  ``-m (--message) intValue``: Sets the verbosity parameter for the
   program itself. *1*: only info messages, *2* info and debug, *3*
   info, debug and trace messages.

Algorithmic parameters:

-  ``-a (--algorithm) intValue``: Specifies the type of algorithm
   employed to solve the instance. *0* is **fullEnumeration**, *1* is
   **innerApproximation**, *2* is **combinatorialPNE**. Note that the
   latter only works with pure equilibria.
-  ``-p (--pure) intValue``: Specifies whether the algorithm should seek
   for a pure equilibrium or not (note that **combinatorialPNE** will
   always seek for a pure solution) . *0* no requirement (either pure or
   mixed), *1* only pure.
-  ``-bigM intValue``: Specifies the whether LCP problems should be
   formulated with a LCP. *0* (default) for indicator constraints, *1*
   for bigM formulation.
-  ``-b (--bound) intValue``: Specifies if the final LCP model should be
   bounded by a *bigM* constant (see option ``-boundBigM``) in their
   primal variables. *0* no bounding, *1* bounding with *bigM*.
-  ``-boundBigM intValue``: Specifies the *bigM* constant for the above
   param (see option ``-b (--bound)``.
-  ``-t (--threads) intValue``: Specifies the number of threads Gurobi
   uses.
-  ``-tl (--timelimit) intValue``: Specifies the timeLimit (in seconds)
   for running the whole algorithmic procedure.

Algorithmic parameters tailored on **innerApproximation**:

-  ``-ad (--add) intValue``: Specifies the *EPECAddPolyMethod* strategy
   employed to expand the inner approximation when no information about
   deviations is available. *0* add polyhedra in lexographic order, *1*
   reverse lexographic, *2* random.
-  ``-ag (--aggr) intValue``: Specifies the *EPECAddPolyMethod*
   aggressivity (see option ``-ad (--add)``). *int* specifies the number
   of lower-level polyhedra to be added to each Stackelberg leader at
   each iteration in which *EPECAddPolyMethod* is triggered.
-  ``-r (--recover) intValue``: Specifies the recover stragegy for the
   **innerApproximation** method. When the algorithm finds a mixed
   equilibrium and the parameter ``-p (--pure)`` is set to *1*, the
   algorithm will search for a pure equilibrium either by randomly
   expanding the inner approximation (*EPECAddPolyMethod*) or by
   triggering a **combinatorialPNE** procedure warmstarted with the
   current information. *0* sets the recover strategy to
   *incrementalEnumeration*, *1* sets it to *combinatorial*.