reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,8cm 12

set output 'Bonds.tex'
set xlabel 'Plate separation (nm)'
set ylabel 'Fraction of bonds formed'

set key off

plot [][0:1] \
     0 with lines lw 1 lc rgbcolor "black" notitle, \
     'bonds.dat' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines

unset output
set terminal X

!epstopdf Bonds-inc.eps
!pdflatex Bonds.tex

!rm Bonds.aux
!rm Bonds.log
!rm Bonds.tex
!rm Bonds-inc.eps
!rm Bonds-inc.pdf
!gs -dBATCH -dNOPAUSE -r100 -sDEVICE=png16m -sOutputFile=Bonds.png Bonds.pdf
!xpdf Bonds.pdf&
