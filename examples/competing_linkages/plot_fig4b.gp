reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,8cm 12

set output 'Fig4b.tex'
set xlabel '$f$'
set ylabel '$F_{\mathrm{min}}$ [$k_{\mathrm{B}}T/A$]' offset 1

#set key left Left reverse spacing 1.3
#set key off
set key bottom left Left reverse spacing 1.2

plot [0:1][-9:0] \
     'fig4b.txt' index 0 using 1:2 with lines lw 4 lc rgbcolor "grey" title 'only $n_{\alpha''}$, $\Delta G_{0\beta} = \infty$', \
     'fig4b.txt' index 1 using 1:2 with lines lw 4 lc rgbcolor "black" title 'only $n_{\beta''}$, $\Delta G_{0\alpha} = \infty$', \
     'fig4b.txt' index 4 using 1:2 with lines lw 4 lc rgbcolor "blue" title '$\delta\Delta G = 3\,k_{\mathrm{B}}T$', \
     'fig4b.txt' index 3 using 1:2 with lines lw 4 lc rgbcolor "green" title '$\delta\Delta G = 5\,k_{\mathrm{B}}T$', \
     'fig4b.txt' index 2 using 1:2 with lines lw 4 lc rgbcolor "red" title '$\delta\Delta G = 8\,k_{\mathrm{B}}T$'

unset output
set terminal X

!epstopdf Fig4b-inc.eps
!pdflatex Fig4b.tex

!rm Fig4b.aux
!rm Fig4b.log
!rm Fig4b.tex
!rm Fig4b-inc.eps
!rm Fig4b-inc.pdf
!xpdf Fig4b.pdf&
