#!/bin/bash -x
mpirun -np 8 ./exampleWindmill 2>&1 | tee log.txt



