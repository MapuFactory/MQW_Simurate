#!/bin/bash
YEAR=$(date +%Y)
MONTH=$(date +%m)
DAY=$(date +%d)
NUM=$1
./a.out $NUM $2 $3> $NUM.csv
echo calc_fin
POTENTIAL_FILE="potential_${YEAR}_${MONTH}_${DAY}_${NUM}.csv"
WAVE_FILE="wave_${YEAR}_${MONTH}_${DAY}_${NUM}.csv"
OUTPUT_FILE="result_${YEAR}_${MONTH}_${DAY}_${NUM}.pdf"
gnuplot <<EOF
set terminal pdf
set output "$OUTPUT_FILE"
set title "Potential vs Wave"
set xlabel "z [nm]"
set ylabel "Potential [eV] / Wave [a.u.]"
set datafile separator ','
set nokey
plot "$POTENTIAL_FILE" using 1:2 with lines , for [i=2:$3] "$WAVE_FILE" using 1:i with lines 
EOF