set datafile separator ','
#set key autotitle columnhead
plot 'cst_result1.txt' u 1:2 with lines lw 2 lc 1 , '' notitle lc 1, \
	 'result_5.dat' u 3:10 with lines,\
	 '../Ka_bwo_axial_1/result_3.dat' u 3:10 with lines,\