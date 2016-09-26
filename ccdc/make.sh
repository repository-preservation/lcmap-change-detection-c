#!/bin/bash

# This parent bash script to make allows for dynamic env configuration.
# Ror example, gsl is now a dependency.

# for example: clean; install; -n

make_arg=$1

BIN=/alt/local/run/contrib/ccdc
GSL_HOME=/alt/local/gsl/1.16.ccdc

GSL_SCI_INC=$GSL_HOME/include GSL_SCI_LIB=$GSL_HOME/lib make $make_arg >& make.out

more make.out

exit


