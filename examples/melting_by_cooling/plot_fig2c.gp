reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 14cm,9cm 12 header '\usepackage{color}'

set output 'Fig2c.tex'
set xlabel '$h / L_\alpha$'
set ylabel '$\beta \Delta F/A$ ($1/L_\alpha^2$)' offset 1

#set key left Left reverse spacing 1.3
set key off

plot \
     'fig2c.txt' using 1:2 with linespoints lw 4 lc rgbcolor "black", \
     0 lc rgbcolor "black"

unset output
set terminal X

!epstopdf Fig2c-inc.eps
!pdflatex Fig2c.tex

!rm Fig2c.aux
!rm Fig2c.log
!rm Fig2c.tex
!rm Fig2c-inc.eps
!rm Fig2c-inc.pdf
!xpdf Fig2c.pdf&
