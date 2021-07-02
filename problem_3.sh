#!/bin/bash
if [ ! -d "bin" ]; then
	mkdir bin
fi
g++ src/problem_3.cpp -o bin/problem_3.exe -Wall -fopenmp -O3 -std=c++11
if [ $? -eq 0 ]; then
	echo "[Success] Binary build"
else
	echo "[Failure] Binary not build"
	exit -1
fi 
if [[ $# -ne 1 ]]; then
	exit 0
fi
bin/problem_3.exe $1
if [ $? -eq 0 ]; then
	echo "[Success] Result files are in ./result/"
else
	echo "[Failure] Program returned abnormally"
fi

