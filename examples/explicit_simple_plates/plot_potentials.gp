reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,8cm 12

set output 'Potentials.tex'
set xlabel '$h / L$'
set ylabel '$\beta F / L^2$'

set key right spacing 1.3
#set key off

set label '$\Delta G^{\,0}_{\alpha\alpha\prime} = -5\,k_{\mathrm{B}}T$' at 1.2, 1
set label '$\Delta G^{\,0}_{\alpha\alpha\prime} = -10\,k_{\mathrm{B}}T$' at 1.0, -2.5
set label '$\Delta G^{\,0}_{\alpha\alpha\prime} = -15\,k_{\mathrm{B}}T$' at 0.7, -4.5

plot [][:6] \
     0 with lines lw 1 lc rgbcolor "black" notitle, \
     \
     'mc_pot_DG-5.00.txt' using 1:($5*20**2) pt 7 ps 1.5 lw 4 lt 1 lc rgbcolor "black" with points title 'Monte Carlo', \
     'pot_DG-5.00.txt' using 1:($3*20**2) lw 4 lt 1 lc rgbcolor "red" with lines title 'Self-consistent with Explicit Tethers', \
     'mf_pot_DG-5.00.txt' using 1:($3*20**2) lw 4 lt 1 lc rgbcolor "blue" with lines title 'Self-consistent Mean Field', \
     \
     'mc_pot_DG-10.00.txt' using 1:($5*20**2) pt 7 ps 1.5 lw 4 lt 1 lc rgbcolor "black" with points notitle, \
     'pot_DG-10.00.txt' using 1:($3*20**2) lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     'mf_pot_DG-10.00.txt' using 1:($3*20**2) lw 4 lt 1 lc rgbcolor "blue" with lines notitle, \
     \
     'mc_pot_DG-15.00.txt' using 1:($5*20**2) pt 7 ps 1.5 lw 4 lt 1 lc rgbcolor "black" with points notitle, \
     'pot_DG-15.00.txt' using 1:($3*20**2) lw 4 lt 1 lc rgbcolor "red" with lines notitle, \
     'mf_pot_DG-15.00.txt' using 1:($3*20**2) lw 4 lt 1 lc rgbcolor "blue" with lines notitle

unset output
set terminal X

!epstopdf Potentials-inc.eps
!pdflatex Potentials.tex

!rm Potentials.aux
!rm Potentials.log
!rm Potentials.tex
!rm Potentials-inc.eps
!rm Potentials-inc.pdf
!xpdf Potentials.pdf&
