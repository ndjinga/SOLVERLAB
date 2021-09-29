#!/bin/bash

for filepath in ./*.dat
do
	filename=$(basename $filepath | cut -d "." -f1)
	gnuplot -c ./epsplotsoln.gnu $filename
done
