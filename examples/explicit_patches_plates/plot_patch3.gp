reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 16cm,12cm solid 12
set output 'Patch3.tex'
set xlabel 'Parallel displacement (nm)'
set ylabel 'Plate separation (nm)'
set key off

#set palette rgb 33,13,10
set palette rgbformulae 22,13,-31
set cbrange [:5]

set contour
set cntrparam levels incr -100,1,0
unset clabel

set view map
splot \
      'patch_results_L10nm_S1_dg-10kT_upperR10nm_lowerR10nm.dat' with pm3d lc rgbcolor "#000000"

unset output
set terminal X

!epstopdf Patch3-inc.eps
!pdflatex Patch3.tex

!rm Patch3.aux
!rm Patch3.log
!rm Patch3.tex
!rm Patch3-inc.eps
!rm Patch3-inc.pdf
!xpdf Patch3.pdf&
