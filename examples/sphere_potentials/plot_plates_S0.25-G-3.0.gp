reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,8cm 12

set output 'FigPlatesS0p25G-3p0.tex'
set xlabel '$h / L$'
set ylabel '$\beta F_{\mathrm{plate}} / A$ ($1/L^2$)' offset 2

set key top right spacing 1.3

set label '$S = 0.25 L$' at 1.2,1
set label '$\beta\Delta G_0 = -3.0$' at 1.2,0.8

plot [][:2] \
     'MirjamResults/results/plates-S0.25-G-3.0.dat' using 1:4 with linespoints lw 4 lc rgbcolor "black" title 'Monte Carlo (Mirjam)', \
     'plates-S0.25-G-3.0.dat' using 1:4 with linespoints lw 4 lc rgbcolor "red" title 'Self-Consistent MF', \
     0 lw 1 lc rgb "black" notitle

unset output
set terminal X

!epstopdf FigPlatesS0p25G-3p0-inc.eps
!pdflatex FigPlatesS0p25G-3p0.tex

!rm FigPlatesS0p25G-3p0.aux
!rm FigPlatesS0p25G-3p0.log
!rm FigPlatesS0p25G-3p0.tex
!rm FigPlatesS0p25G-3p0-inc.eps
!rm FigPlatesS0p25G-3p0-inc.pdf
!xpdf FigPlatesS0p25G-3p0.pdf&
