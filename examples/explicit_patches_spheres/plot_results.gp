reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,8cm 12

set output 'Results.tex'
set xlabel 'Sphere separation (nm)'
set ylabel '$F_{\mathrm{sphere}}$ ($k_{\mathrm{B}} T$)' offset 1

set key off

plot \
     'results.txt' using 1:2 with linespoints lw 4 lc rgbcolor "black", \
     0 lw 1 lc rgb "black" notitle

unset output
set terminal X

!epstopdf Results-inc.eps
!pdflatex Results.tex

!rm Results.aux
!rm Results.log
!rm Results.tex
!rm Results-inc.eps
!rm Results-inc.pdf
!xpdf Results.pdf&
