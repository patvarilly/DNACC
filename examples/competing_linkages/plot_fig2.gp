reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,8cm 12

set output 'Fig2.tex'
set xlabel '$\Delta G_{0\alpha} / k_{\mathrm{B}}T$'
set ylabel '$n_\alpha/N \qquad\qquad n_\beta/N$' offset 1

#set key left Left reverse spacing 1.3
set key off

set label '$\beta\delta\Delta G = 3$'  at -12,0.5
set label '$\beta\delta\Delta G = 5$'  at -13,0.6
set label '$\beta\delta\Delta G = 8$'  at -14,0.7
set label '$\beta\delta\Delta G = 11$' at -15,0.8
set label '$\beta\delta\Delta G = 14$' at -16,0.9

plot [-50:0][0:1] \
     'fig2.txt' index 0 using 1:3 with linespoints lw 4 lc rgbcolor "red", \
     'fig2.txt' index 1 using 1:3 with linespoints lw 4 lc rgbcolor "blue", \
     'fig2.txt' index 2 using 1:3 with linespoints lw 4 lc rgbcolor "grey", \
     'fig2.txt' index 3 using 1:3 with linespoints lw 4 lc rgbcolor "grey", \
     'fig2.txt' index 4 using 1:3 with linespoints lw 4 lc rgbcolor "grey", \
     \
     'fig2.txt' index 0 using 1:2 with linespoints lw 4 lc rgbcolor "black", \
     'fig2.txt' index 1 using 1:2 with linespoints lw 4 lc rgbcolor "black", \
     'fig2.txt' index 2 using 1:2 with linespoints lw 4 lc rgbcolor "black", \
     'fig2.txt' index 3 using 1:2 with linespoints lw 4 lc rgbcolor "black", \
     'fig2.txt' index 4 using 1:2 with linespoints lw 4 lc rgbcolor "black"

unset output
set terminal X

!epstopdf Fig2-inc.eps
!pdflatex Fig2.tex

!rm Fig2.aux
!rm Fig2.log
!rm Fig2.tex
!rm Fig2-inc.eps
!rm Fig2-inc.pdf
!xpdf Fig2.pdf&
