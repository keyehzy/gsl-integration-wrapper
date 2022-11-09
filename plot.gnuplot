set terminal svg enhanced size 500,500
set encoding iso_8859_1
set output 'quantum_metric_2d.svg'

set dgrid3d 100,100
set pm3d explicit
set surface
set view map
set contour
set key outside
# set cntrparam cubicspline
# set cntrparam levels discrete 3.197,3.552
unset colorbox
# set cbrange [0:7000]  # Set the color range of contour values.
# set palette model RGB defined ( 0 'white', 1 'black' )
# set style line 1 lc rgb '#4169E1' pt 7 ps 2
splot 'results/quantum_metric_2d.dat' using 1:2:3 with pm3d notitle
