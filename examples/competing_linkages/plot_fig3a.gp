reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,16cm 12

set output 'Fig3a.tex'
set xlabel '$\Delta G_{0\alpha} / k_{\mathrm{B}}T$'

#set key left Left reverse spacing 1.3
#set key off
set key bottom left Left reverse spacing 1.3

set multiplot

set origin 0,0.52
set size 1,0.48
set ylabel '$n_\alpha/N$ ($h=1$)' offset 1

plot [-30:0][0:1] \
     'fig3a.txt' index 0 using 1:2 with linespoints lw 4 lc rgbcolor "red" title '$S_\alpha = 1.06$', \
     'fig3a.txt' index 2 using 1:2 with linespoints lw 4 lc rgbcolor "blue" title '$S_\alpha = 0.75$', \
     'fig3a.txt' index 4 using 1:2 with linespoints lw 4 lc rgbcolor "black" title '$S_\alpha = 0.53$'

set origin 0,0.0
set size 1,0.48
set ylabel '$n_\alpha/N$ ($h=1.5$)' offset 1

plot [-30:0][0:1] \
     'fig3a.txt' index 1 using 1:3 with linespoints lw 4 lc rgbcolor "red" title '$S_\alpha = 1.06$', \
     'fig3a.txt' index 3 using 1:3 with linespoints lw 4 lc rgbcolor "blue" title '$S_\alpha = 0.75$', \
     'fig3a.txt' index 5 using 1:3 with linespoints lw 4 lc rgbcolor "black" title '$S_\alpha = 0.53$'

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
