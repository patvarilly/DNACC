reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 16cm,12cm solid 12
set output 'Patch4.tex'
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
      'patch_results_L10nm_S0.8_dg-10kT_upperR10nm_lowerR10nm.dat' with pm3d lc rgbcolor "#000000"

unset output
set terminal X

!epstopdf Patch4-inc.eps
!pdflatex Patch4.tex

!rm Patch4.aux
!rm Patch4.log
!rm Patch4.tex
!rm Patch4-inc.eps
!rm Patch4-inc.pdf
!xpdf Patch4.pdf&
