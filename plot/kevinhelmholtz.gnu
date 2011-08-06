#########################################################
file    = 'out/out.xy.slice'
epsfile = 'eps/kevinhelmholtz.eps'
#########################################################
set terminal postscript enhanced size 8cm,8cm eps color

set output epsfile

set key top right

set style line 1 lt 1 lc rgb 'red' lw 2
set style line 2 lt 1 lc rgb 'black' lw 1

set palette defined ( 0 "blue", 1 "white", 2 "red" )

#set lmargin screen 0.1
#set rmargin screen 0.85
#set bmargin screen 0.15
#set tmargin screen 0.9
#set xlabel "x"
#set ylabel "y"
#set xtics format "%g"
#set ytics format "%g"

set lmargin screen 0
set rmargin screen 1
set bmargin screen 0
set tmargin screen 1
unset tics
unset colorbox
set tics format ""

set xrange [0:1]
set yrange [0:1]
set cbrange [0.8:2.2]

plot file u 1:2:3 w image
#########################################################