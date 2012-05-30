reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
# set term epslatex standalone color size 12cm,8cm 12
set term epslatex standalone color colortext size 8.5cm,5.67cm 8
set lmargin 7
set rmargin 1.3
set tmargin 0.6
set bmargin 3

set output 'LCEComparison816.tex'
set xlabel '$\beta \Delta G_0$'
set ylabel 'Number of bonds' offset 1

set key right spacing 1.3
#set key off

plot [-14:-4][:550] \
     '../results/results_8_16.txt' using 2:3 lw 4 lt 1 lc rgbcolor "black" with lines title 'Self-Consistent Theory', \
     '../mc_results_bortolo/sym_A_B_plates/LCE_wrong/lce3' using 1:2 lw 4 lt 2 lc rgbcolor "black" with lines title 'Local Chemical Equilibrium', \
     '../mc_results_bortolo/sym_A_B_plates/LCE_wrong/lce5' using 1:2 lw 4 lt 2 lc rgbcolor "black" with lines notitle, \
     '../mc_results_bortolo/sym_A_B_plates/Run_500_500_n16/N_Dg/n_Dg_h3.0' using 1:2 pt 7 ps 1 lw 4 lt 1 lc rgbcolor "black" with points title 'Monte Carlo', \
     '../mc_results_bortolo/sym_A_B_plates/Run_500_500_n16/N_Dg/n_Dg_h5.0' using 1:2 pt 7 ps 1 lw 4 lt 1 lc rgbcolor "black" with points notitle

#     '../results/results_8_16.txt' using 2:4 lw 4 lt 1 lc rgbcolor "black" with lines title 'Self-consistent Mean Field', \


unset output
set terminal X

!epstopdf LCEComparison816-inc.eps
!pdflatex LCEComparison816.tex

!rm LCEComparison816.aux
!rm LCEComparison816.log
!rm LCEComparison816.tex
!rm LCEComparison816-inc.eps
!rm LCEComparison816-inc.pdf
!mkdir -p ../figures
!mv LCEComparison816.pdf ../figures
!xpdf ../figures/LCEComparison816.pdf&
