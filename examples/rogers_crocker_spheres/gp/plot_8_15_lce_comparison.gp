reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,8cm 12

set output 'LCEComparison815.tex'
set xlabel '$\beta \Delta G_0$'
set ylabel 'Number of bonds' offset 1

set key right spacing 1.3
#set key off

plot [-14:-4][:550] \
     '../results/results_8_15.txt' using 2:3 lw 4 lt 1 lc rgbcolor "black" with lines title 'Self-consistent with Explicit Tethers', \
     '../mc_results_bortolo/sym_A_B_plates/LCE_wrong/lce3' using 1:2 lw 4 lt 2 lc rgbcolor "black" with lines title 'Local Chemical Equilibrium', \
     '../mc_results_bortolo/sym_A_B_plates/LCE_wrong/lce5' using 1:2 lw 4 lt 2 lc rgbcolor "black" with lines notitle, \
     '../mc_results_bortolo/sym_A_B_plates/Run_500_500_n15/N_Dg/n_Dg_h3.0' using 1:2 pt 7 ps 1 lw 4 lt 1 lc rgbcolor "black" with points title 'Monte Carlo', \
     '../mc_results_bortolo/sym_A_B_plates/Run_500_500_n15/N_Dg/n_Dg_h5.0' using 1:2 pt 7 ps 1 lw 4 lt 1 lc rgbcolor "black" with points notitle

#     '../results/results_8_15.txt' using 2:4 lw 4 lt 1 lc rgbcolor "black" with lines title 'Self-consistent Mean Field', \


unset output
set terminal X

!epstopdf LCEComparison815-inc.eps
!pdflatex LCEComparison815.tex

!rm LCEComparison815.aux
!rm LCEComparison815.log
!rm LCEComparison815.tex
!rm LCEComparison815-inc.eps
!rm LCEComparison815-inc.pdf
!mkdir -p ../figures
!mv LCEComparison815.pdf ../figures
!xpdf ../figures/LCEComparison815.pdf&
