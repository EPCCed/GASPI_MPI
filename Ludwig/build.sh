#!/bin/bash 

if [ $# -gt 0 ]
then
    if [ $1 == "GPU" ]
    then
	echo; echo "Building GPU version"
	export TARGETDP=GPU
    fi
fi

echo; echo "Building dummy MPI library:"
cd mpi_s; make clean; make; cd ..

echo; echo "Building targetDP:"
cd targetDP; make clean; make; cd ..

echo; echo "Building Ludwig:"
cd src_GASPI; make clean; make mpi; cd ..






