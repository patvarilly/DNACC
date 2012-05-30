reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,8cm 12

set output 'Fig5a.tex'
set xlabel '$\Delta G_{0\alpha} / k_{\mathrm{B}}T$'
set ylabel '$n_\alpha/N \qquad\qquad n_\beta/N$' offset 1

#set key left Left reverse spacing 1.3
set key off

set label '$\beta\delta\Delta G = 3$'  at -12,0.5
set label '$\beta\delta\Delta G = 5$'  at -13,0.6
set label '$\beta\delta\Delta G = 8$'  at -14,0.7
set label '$\beta\delta\Delta G = 11$' at -15,0.8
set label '$\beta\delta\Delta G = 14$' at -16,0.9

plot [-40:0][0:1] \
     'fig5a.txt' index 0 using 1:3 with lines lw 4 lc rgbcolor "red", \
     'fig5a.txt' index 1 using 1:3 with lines lw 4 lc rgbcolor "blue", \
     'fig5a.txt' index 2 using 1:3 with lines lw 4 lc rgbcolor "brown", \
     'fig5a.txt' index 3 using 1:3 with lines lw 4 lc rgbcolor "purple", \
     'fig5a.txt' index 4 using 1:3 with lines lw 4 lc rgbcolor "magenta", \
     \
     'fig5a.txt' index 0 using 1:2 with lines lw 4 lc rgbcolor "black", \
     'fig5a.txt' index 1 using 1:2 with lines lw 4 lc rgbcolor "black", \
     'fig5a.txt' index 2 using 1:2 with lines lw 4 lc rgbcolor "black", \
     'fig5a.txt' index 3 using 1:2 with lines lw 4 lc rgbcolor "black", \
     'fig5a.txt' index 4 using 1:2 with lines lw 4 lc rgbcolor "black"

unset output
set terminal X

!epstopdf Fig5a-inc.eps
!pdflatex Fig5a.tex

!rm Fig5a.aux
!rm Fig5a.log
!rm Fig5a.tex
!rm Fig5a-inc.eps
!rm Fig5a-inc.pdf
!xpdf Fig5a.pdf&
