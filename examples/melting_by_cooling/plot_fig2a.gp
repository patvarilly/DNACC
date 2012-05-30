reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 14cm,9cm 12 header '\usepackage{color}'

set output 'Fig2a.tex'
set xlabel '$\beta G_{0\alpha}$'
set ylabel '$N_{\mathrm{bond}}/N_{\mathrm{max}}$' offset 1
set y2label 'Minimum of $\beta \Delta F/A$ ($1/L_\alpha^2$)' offset 1

#set key left Left reverse spacing 1.3
set key off

set label '{Strong bonds}' at -14,1.1
set label '{Weak bonds}' at -35,1.25
set label '{$\Delta F_{\mathrm{min}}$}' at -22,-3

plot [-41:7][-3.6:1.5] \
     'fig2a.txt' using 1:2 with linespoints lw 4 lc rgbcolor "red", \
     'fig2a.txt' using 1:3 with linespoints lw 4 lc rgbcolor "blue", \
     'fig2a.txt' using 1:4 with linespoints lw 4 lc rgbcolor "black", \
     0 lc rgbcolor "black"

unset output
set terminal X

!epstopdf Fig2a-inc.eps
!pdflatex Fig2a.tex

!rm Fig2a.aux
!rm Fig2a.log
!rm Fig2a.tex
!rm Fig2a-inc.eps
!rm Fig2a-inc.pdf
!xpdf Fig2a.pdf&
