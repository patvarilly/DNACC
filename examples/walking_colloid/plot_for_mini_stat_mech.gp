reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,8cm 12 header '\usepackage{amsmath}'

set output 'WalkForStatMech.tex'
unset xlabel
set ylabel 'Free energy (1000s $k_{\mathrm{B}}T$)' offset 1

set border 3
set xtics nomirror
set ytics nomirror

#set key left Left reverse spacing 1.5
set key off
#set key top right spacing 1.5

set label 'Mean $\beta G^{\,0} = -6$' at 0.05,2
set arrow from 0.15,3 to 0.2,5
set label 'Mean $\beta G^{\,0} = -12$' at 0.1,18
set arrow from 0.17,17 to 0.15,15

plot \
     'walk-S0.25-G0Mid-12.0-delta8.0.dat' using 1:($2/1000) lw 4 lt 1 lc rgbcolor 'black' with lines title '$(G^{\,0}_{\alpha\beta}+G^{\,0}_{\alpha\gamma})/2 = -12\,k_{\text{B}}T$', \
     'walk-S0.25-G0Mid-10.0-delta8.0.dat' using 1:($2/1000) lw 4 lt 1 lc rgbcolor 'blue' with lines title '$(G^{\,0}_{\alpha\beta}+G^{\,0}_{\alpha\gamma})/2 = -10\,k_{\text{B}}T$', \
     'walk-S0.25-G0Mid-8.0-delta8.0.dat' using 1:($2/1000) lw 4 lt 1 lc rgbcolor 'red' with lines title '$(G^{\,0}_{\alpha\beta}+G^{\,0}_{\alpha\gamma})/2 = -8\,k_{\text{B}}T$', \
     'walk-S0.25-G0Mid-6.0-delta8.0.dat' using 1:($2/1000) lw 4 lt 1 lc rgbcolor '#004000' with lines title '$(G^{\,0}_{\alpha\beta}+G^{\,0}_{\alpha\gamma})/2 = -6\,k_{\text{B}}T$'

unset output
set terminal X

!epstopdf WalkForStatMech-inc.eps
!pdflatex WalkForStatMech.tex

!rm WalkForStatMech.aux
!rm WalkForStatMech.log
!rm WalkForStatMech.tex
!rm WalkForStatMech-inc.eps
!rm WalkForStatMech-inc.pdf
!xpdf WalkForStatMech.pdf&
