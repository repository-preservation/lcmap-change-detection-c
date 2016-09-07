#!/bin/bash

# for example: clean; install; -n
make_arg=$1
BIN=/alt/local/run/contrib/ccdc

GSL_SCI_INC=/alt/local/gsl/1.16.ccdc/include \
GSL_SCI_LIB=/alt/local/gsl/1.16.ccdc/lib \
MATIO_INC=/alt/local/run/contrib/ccdc/matio/include \
MATIO_LIB=/alt/local/run/contrib/ccdc/matio/lib \
HDF5INC=/alt/local/include \
HDF5LIB=/alt/local/hdf5/1.8.15p1/lib \
make $make_arg >& make.out

more make.out

exit
