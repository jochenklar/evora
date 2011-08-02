########################################################################
#file    = "out/out.xy.slice"
#epsfile = "eps/out.xy.eps"
########################################################################
set terminal postscript enhanced size 12cm,6cm eps color
unset key
###
set output epsfile
set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
set xtics format ""
set ytics format ""
set style line 1 lt 1 lc rgb "white" lw 1
###
set multiplot
###
v  = 0.8
v0 = 0.1
vg = 0.1
h  = 0.4
h0 = 0.025
hg = 0.1
###
set tmargin screen v0
set bmargin screen v0+v
set lmargin screen h0
set rmargin screen h0+h
###
set palette defined ( 0 "blue", 1 "white", 2 "red" )
#set logscale cb
set cbrange [0:20]
#set title "log({/Symbol d}(x,y))"
set title "{/Symbol d}(x,y)"
###
plot file u 1:2:3 w image
###
set tmargin screen v0
set bmargin screen v0+v
set lmargin screen h0+h+hg
set rmargin screen h0+2*h+hg
###
set palette defined \
( 0 "white", 25 "blue" , 50 0 0.5 1 , 75 "red", 100 "yellow")
set logscale cb
set cbrange [1:1e2]
set cbtics format "%T"
set ytics format ""
###
set title "log(T(x,y))"
plot file u 1:2:4 w image
###
unset multiplot