reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,8cm 12

set output 'Spheres.tex'
set xlabel 'Sphere separation (nm)'
set ylabel 'Potential of mean force ($k_B T$)'

set key off

plot \
     0 with lines lw 1 lc rgbcolor "black" notitle, \
     'spheres.dat' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines

unset output
set terminal X

!epstopdf Spheres-inc.eps
!pdflatex Spheres.tex

!rm Spheres.aux
!rm Spheres.log
!rm Spheres.tex
!rm Spheres-inc.eps
!rm Spheres-inc.pdf
!gs -dBATCH -dNOPAUSE -r100 -sDEVICE=png16m -sOutputFile=Spheres.png Spheres.pdf
!xpdf Spheres.pdf&
