########################################################################
file    = "out/out.xy.slice"
epsfile = "eps/out.xy.eps"
########################################################################
set terminal postscript enhanced size 10cm,15cm eps color

unset key

set output epsfile

set xrange [-0.5:0.5]
set yrange [-0.5:0.5]

set multiplot

set style line 1 lt 1 lc rgb "white" lw 1

### density plots ####################################################

set bmargin screen 0.68
set tmargin screen 0.96

set palette defined ( 0 "blue", 1 "white", 2 "red" )
set logscale cb
set cbrange [1e-2:1e4]
unset colorbox

set xtics format ""
set ytics format ""

set lmargin screen 0.03
set rmargin screen 0.45

set title "log({/Symbol d}(x,y))"

plot 'out/out.xy.slice' u 1:2:3 w image

unset ylabel
set ytics format ""
set colorbox
set cblabel "log({/Symbol d})" offset -1,0
set cbtics format "%T"

set lmargin screen 0.47
set rmargin screen 0.89

set title "log({/Symbol d}(x,z))"

plot 'out/out.xz.slice' u 1:2:3 w image

### temperature plots ###############################################

set bmargin screen 0.35
set tmargin screen 0.63

set palette defined \
( 0 "white", 25 "blue" , 50 0 0.5 1 , 75 "red", 100 "yellow")
set logscale cb
set cbrange [1:1e8]
unset colorbox

set ytics format ""

set lmargin screen 0.03
set rmargin screen 0.45

set title "log(T(x,y))"

plot 'out/out.xy.slice' u 1:2:14 w image

unset ylabel
set ytics format ""
set colorbox
set cblabel "log(T)" offset 0,0
set cbtics format "%T"

set lmargin screen 0.47
set rmargin screen 0.89

set title "log(T(x,z))"

plot 'out/out.xz.slice' u 1:2:14 w image

### pressure and velocity ###########################################

vscale = 3e8
e = 8

set bmargin screen 0.02
set tmargin screen 0.3

set palette defined \
( 0 "white", 25 "blue" , 50 0 0.5 1 , 75 "red", 100 "yellow")
set logscale cb
set cbrange [1e-5:1e5]
unset colorbox

set ytics format ""

set lmargin screen 0.03
set rmargin screen 0.45

set title "log(p(x,y)) and (v_x,v_y)(x,y)"

plot 'out/out.xy.slice' u 1:2:13 w image,\
     'out/out.xy.slice' u 1:2:($10/vscale):($11/vscale) every e:e w vectors ls 1

unset ylabel
set ytics format ""
set colorbox
set cblabel "log(p)" offset -1,0
set cbtics format "%T"

set lmargin screen 0.47
set rmargin screen 0.89

set title "log(p(x,z)) and (v_x,v_z)(x,z)"

plot 'out/out.xz.slice' u 1:2:13 w image,\
     'out/out.xz.slice' u 1:2:($10/vscale):($12/vscale) every e:e w vectors ls 1

#####################################################################

unset multiplot