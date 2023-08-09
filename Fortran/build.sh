#!/bin/bash

set -xe

FFLAGS="-J ./build/"

mkdir -p ./build/
gfortran $FFLAGS main.f90 -o ./build/fft
