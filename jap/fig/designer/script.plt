reset
set terminal windows

set grid
set xr[0:5]
set yr[-100:0]
set size squar
unset key

plot "plot2.txt" using 1:2 lt -1 w l,\
"plot2.txt" using 1:3 lt -1 w l,\
"plot.txt" using 1:2 lt -1 w l,\
"plot.txt" using 1:3 lt -1 w l

