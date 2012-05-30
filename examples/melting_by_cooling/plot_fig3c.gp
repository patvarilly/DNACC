reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 18cm,9cm 12

set output 'Fig3c.tex'
set xlabel '$\beta \Delta G_{0\alpha}$'
set ylabel 'Minimum of $\beta \Delta F/A$ ($1/L_\alpha^2$)' offset 1

#set key left Left reverse spacing 1.3
set key off

set multiplot

set origin 0,0
set size 0.5,1

set title '$X$-$X''$ pair'

plot [-30:2][-25:1] \
     'fig3c.txt' index 0 using 1:2 with lines lw 4 lc rgbcolor "red", \
     'fig3c.txt' index 1 using 1:2 with lines lw 4 lc rgbcolor "orange", \
     'fig3c.txt' index 2 using 1:2 with lines lw 4 lc rgbcolor "yellow", \
     'fig3c.txt' index 3 using 1:2 with lines lw 4 lc rgbcolor "blue", \
     0 lc rgbcolor "black"

set origin 0.5,0
set size 0.5,1

set title '$X$-$X$ pair'
unset ylabel
set key bottom right spacing 1.3

plot [-30:2][-25:1] \
     'fig3c.txt' index 0 using 1:3 with lines lw 4 lc rgbcolor "red"\
       title '$\beta\delta = 1$', \
     'fig3c.txt' index 1 using 1:3 with lines lw 4 lc rgbcolor "orange"\
       title '$\beta\delta = 3$', \
     'fig3c.txt' index 2 using 1:3 with lines lw 4 lc rgbcolor "yellow"\
       title '$\beta\delta = 5$', \
     'fig3c.txt' index 3 using 1:3 with lines lw 4 lc rgbcolor "blue"\
       title '$\beta\delta = 7$', \
     0 lc rgbcolor "black" notitle

set nomultiplot

unset output
set terminal X

!epstopdf Fig3c-inc.eps
!pdflatex Fig3c.tex

!rm Fig3c.aux
!rm Fig3c.log
!rm Fig3c.tex
!rm Fig3c-inc.eps
!rm Fig3c-inc.pdf
!xpdf Fig3c.pdf&
