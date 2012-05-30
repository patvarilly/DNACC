reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,8cm 12

set output 'FigSpheresR25p0S0p75G-5p0.tex'
set xlabel '$h / L$'
set ylabel '$\beta F_{\mathrm{sphere}}$' offset 1

set key top right spacing 1.3

set label '$R = 25.0 L$' at 0.1,-1.0
set label '$S = 0.75 L$' at 0.1,-1.5
set label '$\beta\Delta G_0 = -5.0$' at 0.1,-2.0

plot [][:2] \
     'MirjamResults/results/spheres-R25.0-S0.75-G-5.0.dat' using 1:4 with linespoints lw 4 lc rgbcolor "black" title 'Monte Carlo (Mirjam)', \
     'spheres-R25.0-S0.75-G-5.0.dat' using 1:4 with linespoints lw 4 lc rgbcolor "red" title 'Self-Consistent MF + Derjaguin', \
     0 lw 1 lc rgb "black" notitle

unset output
set terminal X

!epstopdf FigSpheresR25p0S0p75G-5p0-inc.eps
!pdflatex FigSpheresR25p0S0p75G-5p0.tex

!rm FigSpheresR25p0S0p75G-5p0.aux
!rm FigSpheresR25p0S0p75G-5p0.log
!rm FigSpheresR25p0S0p75G-5p0.tex
!rm FigSpheresR25p0S0p75G-5p0-inc.eps
!rm FigSpheresR25p0S0p75G-5p0-inc.pdf
!xpdf FigSpheresR25p0S0p75G-5p0.pdf&
