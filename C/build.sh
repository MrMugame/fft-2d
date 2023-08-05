#!/bin/bash

set -xe

CFLAGS="-O3 -lm -lc -DNDEBUG"

mkdir -p ./build/
gcc $CFLAGS -o ./build/fft main.c
