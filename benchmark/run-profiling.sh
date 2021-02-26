#!/bin/sh

rm -f *.0 *.rs sys.msg
nsys profile -t nvtx,cuda ../src/flac test-big.inp
