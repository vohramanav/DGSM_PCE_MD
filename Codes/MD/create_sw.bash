#!/bin/bash

filename=$1
i=1
si="Si"
eps=2.1683
sig=2.0951
tol=0.0
cost0=-0.333333333333

while read -r line
do
  new_file="Si$i.sw"
#  if [ $i -lt 10  ]
#  then
#    new_file="Si0$i.sw"
#  fi
  printf "%s %s %s %s %s",$si,$si,$si,$eps,$sig >> "Si.sw"  
  echo $line | awk -F" " '{printf " %s",$5 >> "Si.sw"}'
  echo $line | awk -F" " '{printf " %s",$6 >> "Si.sw"}'
  echo $line | awk -F" " '{printf " %s",$7 >> "Si.sw"}'
  printf " %s",$cost0 >> "Si.sw"
  echo $line | awk -F" " '{printf " %s",$1 >> "Si.sw"}'
  echo $line | awk -F" " '{printf " %s",$2 >> "Si.sw"}'
  echo $line | awk -F" " '{printf " %s",$3 >> "Si.sw"}'
  echo $line | awk -F" " '{printf " %s",$4 >> "Si.sw"}'
  printf " %s",$tol >> "Si.sw"
  sed -i 's/,/ /g' "Si.sw"  
  sed -i 's/     //g' "Si.sw"  
  mv Si.sw $new_file
  let "i += 1"
done < "$filename"

mkdir sw_files
mv *.sw sw_files/.

