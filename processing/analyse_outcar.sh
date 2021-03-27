#!/bin/bash

vasp_outcar="OUTCAR"
if [ "$#" -gt 0 ]; then
    vasp_outcar=$1
fi

LAST=$(awk '{if(/Iteration/) print (NR)}' $vasp_outcar | awk 'END{print $1}')

awk -v bool=0 -v last=$LAST '{if (/ spin component [1-2]/) bool=1; if (/ *k-point *[0-9*]* *: +[0-9.-]* +[0-9.-]* +[0-9.-]* *$/) bool = 1; if (/^$/) bool=0; if (bool && NR>last) print}' $vasp_outcar  > occupation.ground

BANDS=$(awk 'END{print $1}' occupation.ground)

SPIN1=$(awk '{if(/ spin component 1/) print (NR)}' $vasp_outcar | awk 'END{print $1}')

awk -v start=$((SPIN1+3)) -v stop=$((SPIN1+3+BANDS)) '{if (NR > start && NR <= stop) print}' $vasp_outcar  > gamma.1.ground

SPIN2=$(awk '{if(/ spin component 2/) print (NR)}' $vasp_outcar | awk 'END{print $1}')

awk -v start=$((SPIN2+3)) -v stop=$((SPIN2+3+BANDS)) '{if (NR > start && NR <= stop) print}' $vasp_outcar  > gamma.2.ground
