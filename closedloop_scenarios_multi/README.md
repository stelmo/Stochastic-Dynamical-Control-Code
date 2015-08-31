# Closed Loop Scenarios

**It is necessary to run the script `init_sys.jl` in the root directory to load all the external and internal libraries.**

This folder contains stand alone scripts for illustration purposes. They are derived from the scripts in the folder `closedloop`. Since there are many scripts in this folder it is easiest to explain the naming convenction. The obvious abbreviations are what one would expect, e.g. LQG = Linear Quadratic Gaussian etc. The less obvious ones are:

1) `_broken` => the plant controller does not work as expected

2) `_E3` => the observer uses 3 models for inference

3) `_VAR_conf*` => the controller uses chance constraints to a confidence level indicated by the script name

4_ `_mc` => some monte carlo simulation

5) `_break` => the controller does not work as expected (the same as `_broken`)
