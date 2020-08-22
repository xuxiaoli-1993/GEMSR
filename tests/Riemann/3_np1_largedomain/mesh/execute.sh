#!/bin/bash

./mksh.exe
./pgrid.exe
cp ./2D_mesh.dat.* ../
cp ./gems.par ../
preplot ./tecplot_mesh.dat ./tecplot_mesh.plt
