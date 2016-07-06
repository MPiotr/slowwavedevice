set datafile separator ','
set key autotitle columnhead
set xlabel "Current, A"
set ylabel "Power, W"
plot 'result_2.dat' u 3:9 with lines smooth unique lw 2 title "theory: ouput power[W]", \
      'cst_res.dat' u 1:2 with lines lw 2