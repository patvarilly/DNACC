reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,10cm 12

set output 'CompetingSimple.tex'
set xlabel '$\beta \Delta G^{\,0}_{\mathrm{strong}}$'
set ylabel 'Number of bonds' offset 1

set border 3
set xtics nomirror
set ytics nomirror

set key right spacing 1.3
#set key off

#set label '$\Delta G^{\,0}_{\alpha\alpha\prime} = -5\,k_{\mathrm{B}}T$' at 1.2, 1
#set label '$\Delta G^{\,0}_{\alpha\alpha\prime} = -10\,k_{\mathrm{B}}T$' at 1.0, -2.5
#set label '$\Delta G^{\,0}_{\alpha\alpha\prime} = -15\,k_{\mathrm{B}}T$' at 0.7, -4.5

set label '$\beta\Delta G^{\,0}_{\mathrm{weak}} = \beta\Delta G^{\,0}_{\mathrm{weak}} + 8$' at -22,600

plot [-45:0][0:1300] \
     'CompetingData/na_8' using 1:($2*1000) pt 6 ps 1 lw 4 lt 1 lc rgbcolor "black" with points title 'Monte Carlo', \
     'competing_deltaDelta8.txt' using 1:2 lw 4 lt 1 lc rgbcolor "red" with lines title 'Self-consistent with Explicit Tethers', \
     '../competing_linkages/fig2.txt' index 2 using 1:($2*1000) lw 4 lt 1 lc rgbcolor "blue" with lines title 'Self-consistent Mean Field', \
     'CompetingData/nb_8' using 1:($2*1000) pt 6 ps 1 lw 4 lt 1 lc rgbcolor "black" with points notitle, \
     'competing_deltaDelta8.txt' using 1:3 lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     '../competing_linkages/fig2.txt' index 2 using 1:($3*1000) lw 4 lt 1 lc rgbcolor "blue" with lines notitle


unset output
set terminal X

!epstopdf CompetingSimple-inc.eps
!pdflatex CompetingSimple.tex

!rm CompetingSimple.aux
!rm CompetingSimple.log
!rm CompetingSimple.tex
!rm CompetingSimple-inc.eps
!rm CompetingSimple-inc.pdf
!xpdf CompetingSimple.pdf&
