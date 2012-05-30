reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone color size 12cm,8cm 12 header '\usepackage{amsmath}'

set output 'next_sigma0p25_delta5.tex'
set xlabel '$x$'
set ylabel 'Min $\beta F(x)$ - Min $\beta F(x=1)$'

#set key left Left reverse spacing 1.3
#set key off
set key top right spacing 1.3

plot \
     'walk-S0.25-G0Mid-10.0-delta5.0.dat' using 1:2 lw 4 lt 1 lc rgbcolor 'red' with lines title '$\Delta = -10\,k_{\text{B}}T$', \
     'walk-S0.25-G0Mid-15.0-delta5.0.dat' using 1:2 lw 4 lt 1 lc rgbcolor 'blue' with lines title '$\Delta = -15\,k_{\text{B}}T$', \
     'walk-S0.25-G0Mid-20.0-delta5.0.dat' using 1:2 lw 4 lt 1 lc rgbcolor 'black' with lines title '$\Delta = -20\,k_{\text{B}}T$'

unset output
set terminal X

!epstopdf next_sigma0p25_delta5-inc.eps
!pdflatex next_sigma0p25_delta5.tex

!rm next_sigma0p25_delta5.aux
!rm next_sigma0p25_delta5.log
!rm next_sigma0p25_delta5.tex
!rm next_sigma0p25_delta5-inc.eps
!rm next_sigma0p25_delta5-inc.pdf
!xpdf next_sigma0p25_delta5.pdf&
