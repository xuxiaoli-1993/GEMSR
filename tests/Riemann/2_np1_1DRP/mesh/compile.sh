#!/bin/bash

ifort pgrid.f90 -o pgrid.exe -g -traceback
ifort make_1Dnonuniform_mesh.f90 -o mksh.exe -g -trace
./mksh.exe
./pgrid.exe
cp 1D_mesh.dat.1 ../
