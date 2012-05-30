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

set output 'AvgNumBonds.tex'
set xlabel 'Solution hybridisation free energy $\Delta G^0$'
set ylabel 'Number of bonds' offset 2

#set key right spacing 1.3
#set key off
set key bottom left Left reverse spacing 1.3

set multiplot
 
plot [][] \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/Nav.1' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points title 'MC', \
     '../results/general_A_B_avg_results_by_h.txt' every :::0::0 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines title 'SCT', \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/Nav.2' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_avg_results_by_h.txt' every :::1::1 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/Nav.3' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_avg_results_by_h.txt' every :::2::2 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/Nav.4' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_avg_results_by_h.txt' every :::3::3 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/Nav.5' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_avg_results_by_h.txt' every :::4::4 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/Nav.6' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_avg_results_by_h.txt' every :::5::5 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle
     
set size 0.475,0.525
set origin 0.50,0.45
set key off
unset xlabel
unset ylabel
set xtics (-30, -20, -10, 0)
set ytics (0, 50, 100)

plot [][] \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n1.1' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points title 'Monte Carlo', \
     '../results/general_A_B_results_1_by_h.txt' every :::0::0 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines title 'Self-consistent Explicit Tethers', \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n2.1' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_results_1_by_h.txt' every :::1::1 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n3.1' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_results_1_by_h.txt' every :::2::2 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n4.1' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_results_1_by_h.txt' every :::3::3 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n5.1' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_results_1_by_h.txt' every :::4::4 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n6.1' using 1:2 pt 7 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_results_1_by_h.txt' every :::5::5 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle


unset multiplot
unset output
set terminal X

!epstopdf AvgNumBonds-inc.eps
!pdflatex AvgNumBonds.tex

!rm AvgNumBonds.aux
!rm AvgNumBonds.log
!rm AvgNumBonds.tex
!rm AvgNumBonds-inc.eps
!rm AvgNumBonds-inc.pdf
!mkdir -p ../figures
!mv AvgNumBonds.pdf ../figures
!xpdf ../figures/AvgNumBonds.pdf&
