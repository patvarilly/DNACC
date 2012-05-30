reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,8cm 12 header '\usepackage{amsmath}'

set output 'sigma0p25_delta5.tex'
set xlabel '$x$'
set ylabel 'Min $\beta F(x)$ - Min $\beta F(x=1)$'

#set key left Left reverse spacing 1.3
#set key off
set key top right spacing 1.3

plot \
     'walk-S0.25-G0Mid0.0-delta5.0.dat' using 1:2 lw 4 lt 1 lc rgbcolor 'red' with lines title '$\Delta = 0\,k_{\text{B}}T$', \
     'walk-S0.25-G0Mid-2.0-delta5.0.dat' using 1:2 lw 4 lt 1 lc rgbcolor 'blue' with lines title '$\Delta = -2\,k_{\text{B}}T$', \
     'walk-S0.25-G0Mid-4.0-delta5.0.dat' using 1:2 lw 4 lt 1 lc rgbcolor 'black' with lines title '$\Delta = -4\,k_{\text{B}}T$'

unset output
set terminal X

!epstopdf sigma0p25_delta5-inc.eps
!pdflatex sigma0p25_delta5.tex

!rm sigma0p25_delta5.aux
!rm sigma0p25_delta5.log
!rm sigma0p25_delta5.tex
!rm sigma0p25_delta5-inc.eps
!rm sigma0p25_delta5-inc.pdf
!xpdf sigma0p25_delta5.pdf&
