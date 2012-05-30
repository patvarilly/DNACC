reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color colortext size 8.5cm,5.67cm 8
set lmargin 7
set rmargin 1.3
set tmargin 0.2
set bmargin 3

area=(14.14*20)**2
set output 'Output.tex'
set key top left Left reverse spacing 1.3

#Lines for explicit tethers
set style line 1 lt 1 lw 4 lc rgbcolor 'red' 
set style line 2 lt 1 lw 4 lc rgbcolor 'blue' 
set style line 3 lt 1 lw 4 lc rgbcolor 'black' 

#Lines for mean-field
set style line 4 lt 2 lw 4 lc rgbcolor 'red' 
set style line 5 lt 2 lw 4 lc rgbcolor 'blue' 
set style line 6 lt 2 lw 4 lc rgbcolor 'black' 

set multiplot
set size 1,1
set origin 0,0
set xtics nomirror 
set ytics nomirror (0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
set xlabel 'Solution hybridisation free energy $\Delta G^0_a$ ($k_{\mathrm{B}} T$)'
set ylabel 'Bond density (a.u.)' offset 1
set label '[bl]{Strong}' at -7,0.3 textcolor rgbcolor "blue"
set label '[bl]{Weak}' at -19,1.2 textcolor rgbcolor "red"
plot [-25:10][0:2.19] \
     'BONDS_mean_h=1.00-weak' using 1:($2*2*area/200) with lines ls 4 notitle, \
     'BONDS_mean_h=1.00-strong' using 1:($2*area/200) with lines ls 5 notitle, \
     'BONDS_explicit_h=1.00-weak' using 1:($2*2*area/200) with lines ls 1 notitle, \
     'BONDS_explicit_h=1.00-strong' using 1:($2*area/200) with lines ls 2 notitle, \
     './monte_carlo/filebonds.10' using 2:($5/200) w po ps 1.2 pt 7 lt 1 lw 4 lc rgbcolor 'red' notitle, \
     './monte_carlo/filebonds.10' using 2:($4/200) w po ps 1.2 pt 7 lt 1 lw 4 lc rgbcolor 'blue' notitle, \
     0 ls 3 notitle, \
     -1 w po ps 1.2 pt 7 lt 1 lw 4 lc rgbcolor 'black' title 'MC', \
     -1 with lines ls 1 lc rgbcolor "black" title 'SCT', \
     -1 with lines ls 4 lc rgbcolor "black" title 'SCTMF'

set size 0.575,0.6
set origin 0.40,0.375
set xtics nomirror (0.0,0.5,1.0,1.5,2.00)
set ytics nomirror (-6,-4,-2,0,2)
set xlabel  '$h / L_{\mathrm{max}}$' offset 1.5 
set ylabel  '$F / A$ ($k_{\mathrm{B}}T / L^2$)' offset 1
set label  '[bl]{$-5\,k_{\mathrm{B}}T$}' at 1,0.7 textcolor rgbcolor "red"
set label  '[bl]{$-25\,k_{\mathrm{B}}T$}' at 0.8,-1.8 textcolor rgbcolor "black"
set label  '[bl]{$-15\,k_{\mathrm{B}}T$}' at 0.05,-3.5 textcolor rgbcolor "blue" 

set key off
plot [0:2.0][-4:4] \
     'OUTPUT_explicit_-5kt' using 1:4 with lines ls 1, \
     'OUTPUT_explicit_-15kt' using 1:4 with lines ls 2, \
     'OUTPUT_explicit_-25kt' using 1:4 with lines ls 3, \
     'OUTPUT_mean_-5kt' using 1:4 with lines ls 4, \
     'OUTPUT_mean_-15kt' using 1:4 with lines ls 5, \
     'OUTPUT_mean_-25kt' using 1:4 with lines ls 6, \
     './monte_carlo/potential.-5' using 1:2 w points ps 1 pt 7 lt 1 lw 4 lc rgbcolor 'red',\
     './monte_carlo/potential.-15' using 1:2 w po ps 1 pt 7 lt 1 lw 4 lc rgbcolor 'blue',\
     './monte_carlo/potential.-25' using 1:2 w po  ps 1 pt 7 lt 1 lw 4 lc rgbcolor 'black',\
     0 ls 3 notitle

unset multiplot
unset output
set terminal X

!epstopdf Output-inc.eps
!pdflatex Output.tex

!rm Output.aux
!rm Output.log
!rm Output.tex
!rm Output-inc.eps
!rm Output-inc.pdf
!xpdf Output.pdf&



