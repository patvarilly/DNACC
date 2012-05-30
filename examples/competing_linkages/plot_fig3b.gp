reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,8cm 12

set output 'Fig3b.tex'
set xlabel '$\Delta G_{0\alpha} / k_{\mathrm{B}}T$'
set ylabel '$n_\alpha/N \qquad\qquad n_\beta/N$' offset 1

#set key left Left reverse spacing 1.3
set key off

plot [-40:0][0:1] \
     'fig3b.txt' using 1:2 with linespoints lw 4 lc rgbcolor "black", \
     'fig3b.txt' using 1:3 with linespoints lw 4 lc rgbcolor "black"

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
