#!/bin/bash

ifort pgrid.f90 -o pgrid.exe -g -traceback
ifort make_2Dnonuniform_mesh.f90 -o makesh.exe -g -trace
./makesh.exe
./pgrid.exe
cp 2D_mesh.dat.1 ../
