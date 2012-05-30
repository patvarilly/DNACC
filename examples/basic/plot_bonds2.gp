reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,8cm 12

set output 'Bonds2.tex'
set xlabel 'Solution hybridisation free energy ($k_B T$)'
set ylabel 'Fraction of bonds formed'

set key off

plot [][0:1] \
     0 with lines lw 1 lc rgbcolor "black" notitle, \
     'bonds2.dat' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines

unset output
set terminal X

!epstopdf Bonds2-inc.eps
!pdflatex Bonds2.tex

!rm Bonds2.aux
!rm Bonds2.log
!rm Bonds2.tex
!rm Bonds2-inc.eps
!rm Bonds2-inc.pdf
!gs -dBATCH -dNOPAUSE -r100 -sDEVICE=png16m -sOutputFile=Bonds2.png Bonds2.pdf
!xpdf Bonds2.pdf&
