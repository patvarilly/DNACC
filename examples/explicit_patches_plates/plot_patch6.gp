reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
#set term epslatex standalone color size 16cm,12cm solid 12
set term epslatex standalone color size 8.5cm,10cm solid 9

unset xlabel
unset ylabel
unset xtics
unset ytics
unset colorbox

set output 'Patch6.tex'
set key off

#set palette rgb 33,13,10
set palette rgbformulae 22,13,-31
#set palette rgbformulae 33,13,10
set cbrange [-10:10]

set contour
set cntrparam levels incr -10,1,0
unset clabel

set multiplot


# Frame
set origin 0, 0.0
set size 1.0, 1.0
set tmargin at screen 0.975
set bmargin at screen 0.150
set colorbox
set xlabel 'Parallel displacement (nm)' offset 1
set ylabel 'Plate separation (nm)' offset -3
set xtics (-200,-100,0,100,200)
set ytics 0,15,0
set cbtics -10,2,10

set view map
splot [-240:240][10:40] 10 with pm3d


# Rest
set ytics 0,10,39
set mytics 2
unset colorbox
unset xlabel
unset ylabel
unset xtics

set origin 0, 0.67
set size 1.0, 0.33
set tmargin at screen 0.975
set bmargin at screen 0.700

set view map
splot [-240:240][10:40] \
      'patch_results_L20nm_S1_dg-6kT_upperR100nm_lowerR100nm.dat' with pm3d lc rgbcolor "#000000"


set key off

set origin 0, 0.33
set size 1.0, 0.33
set tmargin at screen 0.700
set bmargin at screen 0.425

splot [-240:240][10:40] \
      'patch_results_L20nm_S1_dg-7kT_upperR100nm_lowerR100nm.dat' with pm3d lc rgbcolor "#000000"

set origin 0, 0.00
set size 1.0, 0.33
set tmargin at screen 0.425
set bmargin at screen 0.15

splot [-240:240][10:40] \
      'patch_results_L20nm_S1_dg-8kT_upperR100nm_lowerR100nm.dat' with pm3d lc rgbcolor "#000000"

unset multiplot

unset output
set terminal X

!epstopdf Patch6-inc.eps
!pdflatex Patch6.tex

!rm Patch6.aux
!rm Patch6.log
!rm Patch6.tex
!rm Patch6-inc.eps
!rm Patch6-inc.pdf
!xpdf Patch6.pdf&
