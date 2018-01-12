#!/bin/bash

for f in $1/*.log
do
  #awk -F" " '$2=="TW"{print substr($NF,1,3) >> "Lr.txt";}' $f
  awk -F" " '$2=="TW"{if (substr($NF,1,1)=="4" || substr($NF,1,1)=="5") print substr($NF,1,6) >> "Lr.txt";}' $f
done
