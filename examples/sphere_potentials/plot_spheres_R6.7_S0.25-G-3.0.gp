reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,8cm 12

set output 'FigSpheresR6p7S0p25G-3p0.tex'
set xlabel '$h / L$'
set ylabel '$\beta F_{\mathrm{sphere}}$' offset 1

set key top right spacing 1.3

set label '$R = 6.7 L$' at 0.1,-3.0
set label '$S = 0.25 L$' at 0.1,-4.0
set label '$\beta\Delta G_0 = -3.0$' at 0.1,-5.0

plot [][:2] \
     'MirjamResults/results/spheres-R6.7-S0.25-G-3.0.dat' using 1:4 with linespoints lw 4 lc rgbcolor "black" title 'Monte Carlo (Mirjam)', \
     'spheres-R6.7-S0.25-G-3.0.dat' using 1:4 with linespoints lw 4 lc rgbcolor "red" title 'Self-Consistent MF + Derjaguin', \
     0 lw 1 lc rgb "black" notitle

unset output
set terminal X

!epstopdf FigSpheresR6p7S0p25G-3p0-inc.eps
!pdflatex FigSpheresR6p7S0p25G-3p0.tex

!rm FigSpheresR6p7S0p25G-3p0.aux
!rm FigSpheresR6p7S0p25G-3p0.log
!rm FigSpheresR6p7S0p25G-3p0.tex
!rm FigSpheresR6p7S0p25G-3p0-inc.eps
!rm FigSpheresR6p7S0p25G-3p0-inc.pdf
!xpdf FigSpheresR6p7S0p25G-3p0.pdf&
