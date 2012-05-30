reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 18cm,9cm 12

set output 'Fig3b.tex'
set xlabel '$\beta \Delta G_{0\alpha}$'
set ylabel 'Minimum of $\beta \Delta F/A$ ($1/L_\alpha^2$)' offset 1

#set key left Left reverse spacing 1.3
set key off

set multiplot

set origin 0,0
set size 0.5,1

set title '$X$-$X''$ pair'

plot [-30:2][-4.5:1] \
     'fig3b.txt' index 0 using 1:2 with lines lw 4 lc rgbcolor "red", \
     'fig3b.txt' index 1 using 1:2 with lines lw 4 lc rgbcolor "orange", \
     'fig3b.txt' index 2 using 1:2 with lines lw 4 lc rgbcolor "yellow", \
     'fig3b.txt' index 3 using 1:2 with lines lw 4 lc rgbcolor "blue", \
     0 lc rgbcolor "black"

set origin 0.5,0
set size 0.5,1

set title '$X$-$X$ pair'
unset ylabel
set key bottom right spacing 1.3

plot [-30:2][-4.5:1] \
     'fig3b.txt' index 0 using 1:3 with lines lw 4 lc rgbcolor "red"\
       title '$L_{\mathrm{incrt}} / L_\alpha = 1.0$', \
     'fig3b.txt' index 1 using 1:3 with lines lw 4 lc rgbcolor "orange"\
       title '$L_{\mathrm{incrt}} / L_\alpha = 1.2$', \
     'fig3b.txt' index 2 using 1:3 with lines lw 4 lc rgbcolor "yellow"\
       title '$L_{\mathrm{incrt}} / L_\alpha = 1.5$', \
     'fig3b.txt' index 3 using 1:3 with lines lw 4 lc rgbcolor "blue"\
       title '$L_{\mathrm{incrt}} / L_\alpha = 1.8$', \
     0 lc rgbcolor "black" notitle

set nomultiplot

unset output
set terminal X

!epstopdf Fig3b-inc.eps
!pdflatex Fig3b.tex

!rm Fig3b.aux
!rm Fig3b.log
!rm Fig3b.tex
!rm Fig3b-inc.eps
!rm Fig3b-inc.pdf
!xpdf Fig3b.pdf&
