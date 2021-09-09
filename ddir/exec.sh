#!/bin/bash 
time mpiexec -np 4 ./wabbit PARAMS_sphere.ini --memory=16.0GB
