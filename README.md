# ADS_B_refractivity

Adjoint-based inversion of refracted ADS-B signal observations to retrieve atmospheric refractivity structure.

Run "make" in the command line, then use ./optim [command] to generate data.  
For the python files in "tools" and "data", run ./X.py -h to see a list of optional arguments.

The C++ code makes use of the spline library: https://kluge.in-chemnitz.de/opensource/spline/
