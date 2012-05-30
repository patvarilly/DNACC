reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,8cm 12

set output 'Fig4a.tex'
set xlabel '$\Delta G_{0\beta} / k_{\mathrm{B}}T$'
set ylabel '$F_{\mathrm{min}}$ [$k_{\mathrm{B}}T/A$]' offset 1

#set key left Left reverse spacing 1.3
#set key off
set key bottom right spacing 1.3

plot [-20:8][-15:0] \
     'fig4a.txt' index 0 using 1:2 with lines lw 4 lc rgbcolor "black" title 'only $n_{\beta''}$, $\Delta G_{0\alpha} = \infty$', \
     'fig4a.txt' index 1 using 1:2 with lines lw 4 lc rgbcolor "red" title '$\delta\Delta G = 8\,k_{\mathrm{B}}T$', \
     'fig4a.txt' index 2 using 1:2 with lines lw 4 lc rgbcolor "green" title '$\delta\Delta G = 5\,k_{\mathrm{B}}T$', \
     'fig4a.txt' index 3 using 1:2 with lines lw 4 lc rgbcolor "blue" title '$\delta\Delta G = 3\,k_{\mathrm{B}}T$'

unset output
set terminal X

!epstopdf Fig4a-inc.eps
!pdflatex Fig4a.tex

!rm Fig4a.aux
!rm Fig4a.log
!rm Fig4a.tex
!rm Fig4a-inc.eps
!rm Fig4a-inc.pdf
!xpdf Fig4a.pdf&
