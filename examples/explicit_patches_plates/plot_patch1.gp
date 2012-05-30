reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 16cm,12cm solid 12
set output 'Patch1.tex'
set xlabel 'Parallel displacement (nm)'
set ylabel 'Plate separation (nm)'
set key off

#set palette rgb 33,13,10
set palette rgbformulae 22,13,-31
set cbrange [:100]

set contour
set cntrparam levels incr -100,10,0
unset clabel

set view map
splot \
      'patch_results_L20nm_S1_dg-10kT_upperR100nm_lowerR100nm.dat' with pm3d lc rgbcolor "#000000"

unset output
set terminal X

!epstopdf Patch1-inc.eps
!pdflatex Patch1.tex

!rm Patch1.aux
!rm Patch1.log
!rm Patch1.tex
!rm Patch1-inc.eps
!rm Patch1-inc.pdf
!xpdf Patch1.pdf&
