reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,8cm 12

set output 'ABPlatePotentials.tex'
set xlabel 'Plate separation (nm)'
set ylabel 'Potential ($10^{-3} \times k_{\mathrm{B}} T / $nm$^2$)' offset 1

set key right spacing 1.3
#set key off

plot [][:1]\
     '../results/A_B_plate_potentials.txt' using 2:($3*1e3) lw 4 lt 1 lc rgbcolor "black" with lines title 'Self-consistent, Explicit Tethers', \
     '../results/A_B_plate_potentials.txt' using 2:($4*1e3) lw 4 lt 2 lc rgbcolor "black" with lines title 'Self-consistent, Mean Field'


unset output
set terminal X

!epstopdf ABPlatePotentials-inc.eps
!pdflatex ABPlatePotentials.tex

!rm ABPlatePotentials.aux
!rm ABPlatePotentials.log
!rm ABPlatePotentials.tex
!rm ABPlatePotentials-inc.eps
!rm ABPlatePotentials-inc.pdf
!mkdir -p ../figures
!mv ABPlatePotentials.pdf ../figures
!xpdf ../figures/ABPlatePotentials.pdf&
