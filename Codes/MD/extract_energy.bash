#!/bin/bash

for f in $1/*.log
do
  INPUT="$(basename $f .log)"
  echo $INPUT | awk -F"_" '{printf "%s",$2 >> "energy.txt";}'
  awk -F" " '$1=="500000"{printf " %s %s\n",$(NF-1),$NF >> "energy.txt";}' $f 	
done
