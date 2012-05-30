reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,8cm 12

set output 'AvgNumBonds1.tex'
set xlabel '$\beta \Delta G_0$'
set ylabel 'Number of bonds' offset 1

set key right spacing 1.3
#set key off

plot [][] \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n1.1' using 1:2 pt 6 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points title 'Monte Carlo', \
     '../results/general_A_B_results_1_by_h.txt' every :::0::0 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines title 'Self-consistent Explicit Tethers', \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n2.1' using 1:2 pt 6 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_results_1_by_h.txt' every :::1::1 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n3.1' using 1:2 pt 6 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_results_1_by_h.txt' every :::2::2 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n4.1' using 1:2 pt 6 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_results_1_by_h.txt' every :::3::3 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n5.1' using 1:2 pt 6 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_results_1_by_h.txt' every :::4::4 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../mc_results_bortolo/A_B_results/spheres_mc_explicit_NBONDS/n6.1' using 1:2 pt 6 ps 0.5 lw 1 lt 1 lc rgbcolor "black" with points notitle, \
     '../results/general_A_B_results_1_by_h.txt' every :::5::5 using 1:6 lw 4 lt 1 lc rgbcolor "red" with lines notitle
     


unset output
set terminal X

!epstopdf AvgNumBonds1-inc.eps
!pdflatex AvgNumBonds1.tex

!rm AvgNumBonds1.aux
!rm AvgNumBonds1.log
!rm AvgNumBonds1.tex
!rm AvgNumBonds1-inc.eps
!rm AvgNumBonds1-inc.pdf
!mkdir -p ../figures
!mv AvgNumBonds1.pdf ../figures
!xpdf ../figures/AvgNumBonds1.pdf&
