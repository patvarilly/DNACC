reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
#set term epslatex standalone color colortext size 12cm,8cm 10
set term epslatex standalone color colortext size 8.5cm,5.67cm 8
set lmargin 7
set rmargin 1.3
set tmargin 0.6
set bmargin 3

# Monte Carlo potentials are given in units of kT/L1**2, L1=20nm
# 'Theory' potentials are given in units of kT/L2**2, L1=30nm
# Need rescaling for correct comparison!
scale=20.0**2/30.0**2
max_num_strong_bonds=144

area=(10*20)**2
set output 'Output.tex'
set key top left Left reverse spacing 1.3

#set key at -15,2.1 

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
set xlabel 'Solution hybridisation energy $\Delta G^0_\alpha$'
set ylabel 'Bond density (a.u.)' offset 1
set label '[bl]{Strong}' at -7,0.3 textcolor rgbcolor "blue"
set label '[bl]{Weak}' at -19,0.9 textcolor rgbcolor "red"
set label "[bl]{$c$--$c'$}" at -22.5,0.03 textcolor rgbcolor "black"

plot [-25:10][0.0:1.6] \
     'BONDS_mean_h=1.00-weak' using 1:($2*area/max_num_strong_bonds) with lines ls 4 notitle,\
     'BONDS_mean_h=1.00-strong' using 1:($2*area/max_num_strong_bonds) with lines ls 5 notitle,\
     'BONDS_mean_h=1.00-ggp' using 1:($2*area/max_num_strong_bonds) with lines ls 6 notitle,\
     'BONDS_explicit_h=1.00-weak' using 1:($2*area/max_num_strong_bonds) with lines ls 1 notitle,\
     'BONDS_explicit_h=1.00-strong' using 1:($2*area/max_num_strong_bonds) with lines ls 2 notitle,\
     'BONDS_explicit_h=1.00-ggp' using 1:($2*area/max_num_strong_bonds) with lines ls 3 notitle,\
     './monte_carlo/filebonds.5' using 2:($5/max_num_strong_bonds) w po ps 1.2 pt 7 lt 1 lw 2 lc rgbcolor 'red' notitle, \
     './monte_carlo/filebonds.5' using 2:($4/max_num_strong_bonds) w po ps 1.2 pt 7 lt 1 lw 2 lc rgbcolor 'blue' notitle, \
     './monte_carlo/ggpbonds.5' using ($4-5):($3/max_num_strong_bonds) w po ps 1.2 pt 7 lt 1 lw 2 lc rgbcolor 'black' notitle, \
     -1 w po ps 1.2 pt 7 lt 1 lw 4 lc rgbcolor 'black' title 'MC', \
     -1 with lines ls 1 lc rgbcolor "black" title 'SCT', \
     -1 with lines ls 4 lc rgbcolor "black" title 'SCTMF'

set size 0.575,0.6
set origin 0.40,0.375
set xtics nomirror (0.0,0.5,1.0,1.5,2.00)
set ytics nomirror (-6,-4,-2,0,2)
set xlabel  '$h / L_{\mathrm{max}}$' offset 1.5 
set ylabel  '$F / A (k_{\mathrm{B}}T / L^2)$' offset 1.5
set label  '$-5\,k_{\mathrm{B}}T$' at 1.2,1.5 textcolor rgbcolor "red"
set label  '$-15\,k_{\mathrm{B}}T$' at 1.2,-2.3 textcolor rgbcolor "blue"
set label  '$-25\,k_{\mathrm{B}}T$' at 1.1,-1.2 textcolor rgbcolor "black"
set key off
plot [0:2.2][-3.5:3.5] \
     'OUTPUT_explicit_-5kt' using ($1*1.5):($4*scale) with lines ls 1, \
     'OUTPUT_explicit_-15kt' using ($1*1.5):($4*scale) with lines ls 2, \
     'OUTPUT_explicit_-25kt' using ($1*1.5):($4*scale) with lines ls 3, \
     'OUTPUT_mean_-5kt' using ($1*1.5):($4*scale) with lines ls 4, \
     'OUTPUT_mean_-15kt' using ($1*1.5):($4*scale) with lines ls 5, \
     'OUTPUT_mean_-25kt' using ($1*1.5):($4*scale) with lines ls 6, \
     './monte_carlo/potential.-5' using 1:2 w points ps 1 pt 7 lt 1 lw 4 lc rgbcolor 'red',\
     './monte_carlo/potential.-15' using 1:2 w points ps 1 pt 7 lt 1 lw 4 lc rgbcolor 'blue',\
     './monte_carlo/potential.-25' using 1:2 w points ps 1 pt 7 lt 1 lw 4 lc rgbcolor 'black'

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
