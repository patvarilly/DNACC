reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 18cm,9cm 12

set output 'Fig3a.tex'
set xlabel '$\beta \Delta G_{0\alpha}$'
set ylabel 'Minimum of $\beta \Delta F/A$ ($1/L_\alpha^2$)' offset 1

#set key left Left reverse spacing 1.3
set key off

set multiplot

set origin 0,0
set size 0.5,1

set title '$X$-$X''$ pair'

plot [-30:2][-15:1] \
     'fig3a.txt' index 0 using 1:2 with lines lw 4 lc rgbcolor "red", \
     'fig3a.txt' index 1 using 1:2 with lines lw 4 lc rgbcolor "orange", \
     'fig3a.txt' index 2 using 1:2 with lines lw 4 lc rgbcolor "yellow", \
     'fig3a.txt' index 3 using 1:2 with lines lw 4 lc rgbcolor "blue", \
     0 lc rgbcolor "black"

set origin 0.5,0
set size 0.5,1

set title '$X$-$X$ pair'
unset ylabel
set key bottom right spacing 1.3

plot [-30:2][-15:1] \
     'fig3a.txt' index 0 using 1:3 with lines lw 4 lc rgbcolor "red"\
       title '$S_\alpha / L_\alpha = 1.5$', \
     'fig3a.txt' index 1 using 1:3 with lines lw 4 lc rgbcolor "orange"\
       title '$S_\alpha / L_\alpha = 1.0$', \
     'fig3a.txt' index 2 using 1:3 with lines lw 4 lc rgbcolor "yellow"\
       title '$S_\alpha / L_\alpha = 0.7$', \
     'fig3a.txt' index 3 using 1:3 with lines lw 4 lc rgbcolor "blue"\
       title '$S_\alpha / L_\alpha = 0.5$', \
     0 lc rgbcolor "black" notitle

set nomultiplot

unset output
set terminal X

!epstopdf Fig3a-inc.eps
!pdflatex Fig3a.tex

!rm Fig3a.aux
!rm Fig3a.log
!rm Fig3a.tex
!rm Fig3a-inc.eps
!rm Fig3a-inc.pdf
!xpdf Fig3a.pdf&
