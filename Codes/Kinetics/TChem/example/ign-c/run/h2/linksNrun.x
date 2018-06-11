#!/bin/bash

RUNDIR=`pwd`
SRCDIR="$RUNDIR/../../src"
DATDIR="$RUNDIR/../../../../data"
PYDIR="$RUNDIR/../../../python"

echo "--------------------------------------------------------"
echo "run directory   : $RUNDIR"
echo "src directory   : $SRCDIR"
echo "python directory: $PYDIR"
echo "--------------------------------------------------------"
echo

cd $SRCDIR; make clean; make output=yes gettig=yes; cd $RUNDIR
ln -fs $DATDIR/periodictable.dat .

if [ -f $DATDIR/chem_h2.inp ]; then
  ln -fs $DATDIR/chem_h2.inp chem.inp
else
  echo " Missing $DATDIR/chem_h2.inp"
  exit
fi
if [ -f $DATDIR/therm_h2.dat ]; then
  ln -fs $DATDIR/therm_h2.dat therm.dat
else
  echo " Missing $DATDIR/therm_h2.inp"
  exit
fi

ln -fs $SRCDIR/ign .
/bin/cp -f input.setup input.dat
./ign 

