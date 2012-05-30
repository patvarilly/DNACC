reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
#set term epslatex standalone color size 12cm,8cm 12
set term epslatex standalone color colortext size 8.5cm,5.67cm 8
set lmargin 7
set rmargin 1.3
set tmargin 0.6
set bmargin 3

set output 'ABABSpherePotentials.tex'
set xlabel 'Sphere separation (nm)'
set ylabel 'Potential ($k_{\mathrm{B}} T$)' offset 1

set key right spacing 1.3
#set key off

plot [:35][-15:30]\
     '../results/AB_AB_sphere_potentials.txt' every :::0::0 using 2:3 lw 4 lt 1 lc rgbcolor "black" with lines title 'SCT Plates + Derjaguin', \
     '../results/AB_AB_sphere_potentials.txt' every :::2::2 using 2:3 lw 4 lt 1 lc rgbcolor "black" with lines notitle, \
     '../results/AB_AB_sphere_potentials.txt' every :::4::4 using 2:3 lw 4 lt 1 lc rgbcolor "black" with lines notitle, \
     \
     '../results/AB_AB_sphere_potentials.txt' every :::0::0 using 2:4 lw 4 lt 2 lc rgbcolor "black" with lines title 'SCTMF Plates + Derjaguin', \
     '../results/AB_AB_sphere_potentials.txt' every :::2::2 using 2:4 lw 4 lt 2 lc rgbcolor "black" with lines notitle, \
     '../results/AB_AB_sphere_potentials.txt' every :::4::4 using 2:4 lw 4 lt 2 lc rgbcolor "black" with lines notitle, \
     \
     '../mc_results_bortolo/AB_AB_results/spheres_mc_derjaguin_POT/R32' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines title 'MC Plates + Derjaguin', \
     '../mc_results_bortolo/AB_AB_results/spheres_mc_derjaguin_POT/R27' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/AB_AB_results/spheres_mc_derjaguin_POT/R24' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines notitle


unset output
set terminal X

!epstopdf ABABSpherePotentials-inc.eps
!pdflatex ABABSpherePotentials.tex

!rm ABABSpherePotentials.aux
!rm ABABSpherePotentials.log
!rm ABABSpherePotentials.tex
!rm ABABSpherePotentials-inc.eps
!rm ABABSpherePotentials-inc.pdf
!mkdir -p ../figures
!mv ABABSpherePotentials.pdf ../figures
!xpdf ../figures/ABABSpherePotentials.pdf&
