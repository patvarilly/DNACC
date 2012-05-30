reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
#set term epslatex standalone color colortext size 12cm,8cm 10
set term epslatex standalone color colortext size 8.5cm,5.67cm 8
set lmargin 6
set rmargin 1.3
set tmargin 0.2
set bmargin 3

set output 'Output.tex'
set xlabel 'Solution hybridisation free energy $\Delta G^0$ ($k_{\mathrm{B}}T$)'
set ylabel 'Bond density (a.u.)' offset 2

set key bottom left Left reverse spacing 1.3

area=200*200
#Lines for explicit tethers
set style line 1 lt 1 lw 4 lc rgbcolor 'red' 
set style line 2 lt 1 lw 4 lc rgbcolor 'blue' 
set style line 3 lt 1 lw 4 lc rgbcolor 'black' 

#Lines for mean-field
set style line 4 lt 2 lw 4 lc rgbcolor 'red' 
set style line 5 lt 2 lw 4 lc rgbcolor 'blue' 
set style line 6 lt 2 lw 4 lc rgbcolor 'black' 

set multiplot
 
set ytics nomirror (0.0,0.2,0.4,0.6,0.8,1.0)
set xtics nomirror

plot [-25:15][0:1.2] \
     'monte_carlo/bonds_h_20nm.dat' using 1:($2/100) with points ps 1.2 pt 7 lw 4 lc rgbcolor "black" title 'MC',\
     'BONDS-h-20nm-explicit.dat' using 1:($2*area/100) with lines ls 3 title 'SCT', \
     'BONDS-h-20nm-meanfield.dat' using 1:($2*area/100) with lines ls 6 title 'SCTMF'

set size 0.525,0.675
set origin 0.45,0.30
set xtics nomirror (0.0,0.5,1.0,1.5,2.00)
set ytics auto
set mytics 5
set xlabel  '$h/L_{\mathrm{max}}$' offset 1.5 
set ylabel  '$F / A$ ($k_{\mathrm{B}}T / L^2$)' offset 1.5
set label  '[bl]{$-4\,k_{\mathrm{B}}T$}' at 0.7,1.3 textcolor rgbcolor "red"
set label  '[tl]{$-14\,k_{\mathrm{B}}T$}' at 0.45,-5.5 textcolor rgbcolor "blue"
set label  '[bl]{$-24\,k_{\mathrm{B}}T$}' at 0.25,-12.5 textcolor rgbcolor "black"
set key off
plot [0:2][-16:4] \
     'POTENTIAL_dg-4kT-explicit.dat' using 1:4 with lines ls 1, \
     'POTENTIAL_dg-14kT-explicit.dat' using 1:4 with lines ls 2, \
     'POTENTIAL_dg-24kT-explicit.dat' using 1:4 with lines ls 3, \
     'POTENTIAL_dg-4kT-meanfield.dat' using 1:4 with lines ls 4, \
     'POTENTIAL_dg-14kT-meanfield.dat' using 1:4 with lines ls 5, \
     'POTENTIAL_dg-24kT-meanfield.dat' using 1:4 with lines ls 6, \
     './monte_carlo/potential.-4' using 1:2 w po ps 1 pt 7 lt 1 lw 4 lc rgbcolor 'red',\
     './monte_carlo/potential.-14' using 1:2 w po ps 1 pt 7 lt 1 lw 4 lc rgbcolor 'blue',\
     './monte_carlo/potential.-24' using 1:2 w po ps 1 pt 7 lt 1 lw 4 lc rgbcolor 'black',\
     0 ls 3

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


