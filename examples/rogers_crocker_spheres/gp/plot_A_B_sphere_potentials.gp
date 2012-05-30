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

set output 'ABSpherePotentials.tex'
set xlabel 'Sphere separation (nm)'
set ylabel 'Potential ($k_{\mathrm{B}} T$)' offset 1

set key right spacing 1.3
#set key off

plot [:35][-15:30]\
     '../results/A_B_sphere_potentials.txt' every :::0::0 using 2:3 lw 4 lt 1 lc rgbcolor "black" with lines title 'SCT Plates + Derjaguin', \
     '../results/A_B_sphere_potentials.txt' every :::2::2 using 2:3 lw 4 lt 1 lc rgbcolor "black" with lines notitle, \
     '../results/A_B_sphere_potentials.txt' every :::4::4 using 2:3 lw 4 lt 1 lc rgbcolor "black" with lines notitle, \
     \
     '../results/A_B_sphere_potentials.txt' every :::0::0 using 2:4 lw 4 lt 2 lc rgbcolor "black" with lines title 'SCTMF Plates + Derjaguin', \
     '../results/A_B_sphere_potentials.txt' every :::2::2 using 2:4 lw 4 lt 2 lc rgbcolor "black" with lines notitle, \
     '../results/A_B_sphere_potentials.txt' every :::4::4 using 2:4 lw 4 lt 2 lc rgbcolor "black" with lines notitle, \
     \
     '../mc_results_bortolo/A_B_results/spheres_mc_derjaguin_POT/R36.0' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines title 'MC Plates + Derjaguin', \
     '../mc_results_bortolo/A_B_results/spheres_mc_derjaguin_POT/R30.5' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_derjaguin_POT/R33.0' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_derjaguin_POT/R30.5' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     \
     '../results/rogers_crocker_A_B_avg_results_by_h.txt' every ::0::0 using 2:3 pt 6 ps 1 lw 4 lt 1 lc rgbcolor "blue" with points title 'SCT Spheres', \
     '../results/rogers_crocker_A_B_avg_results_by_h.txt' every ::2::2 using 2:3 pt 6 ps 1 lw 4 lt 1 lc rgbcolor "blue" with points notitle, \
     '../results/rogers_crocker_A_B_avg_results_by_h.txt' every ::4::4 using 2:3 pt 6 ps 1 lw 4 lt 1 lc rgbcolor "blue" with points notitle, \
     \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_POT/V_36.0' using 1:2 pt 6 ps 1 lw 4 lt 1 lc rgbcolor "black" with points title 'MC Spheres', \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_POT/V_33.0' using 1:2 pt 6 ps 1 lw 4 lt 1 lc rgbcolor "black" with points notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_POT/V_30.5' using 1:2 pt 6 ps 1 lw 4 lt 1 lc rgbcolor "black" with points notitle


unset output
set terminal X

!epstopdf ABSpherePotentials-inc.eps
!pdflatex ABSpherePotentials.tex

!rm ABSpherePotentials.aux
!rm ABSpherePotentials.log
!rm ABSpherePotentials.tex
!rm ABSpherePotentials-inc.eps
!rm ABSpherePotentials-inc.pdf
!mkdir -p ../figures
!mv ABSpherePotentials.pdf ../figures
!xpdf ../figures/ABSpherePotentials.pdf&
