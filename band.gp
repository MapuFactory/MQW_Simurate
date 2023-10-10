set terminal pdf
set output out
set title "Potential vs Wave"
set xlabel "n"
set ylabel "Potential [eV] / Wave [a.u.]"
set datafile separator ','
plot pot using 1:2 with lines title 'Potential', for [i=2:16] wav using 1:i with lines title sprintf('Wave %d', i-1)