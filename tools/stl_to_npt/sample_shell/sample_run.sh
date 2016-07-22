#!/bin/bash -x

# sample shell
mpirun -np 4 ../stl_to_npt sphere.stl 2>&1 | tee log.txt


#degree=40
#mpirun -np 4 ../stl_to_npt sphere.stl ${degree} 2>&1 | tee log.txt

